"""Graphical front end: a web page served by this computer.

`pgc-risk gui` starts a small HTTP server on 127.0.0.1, opens the page in the browser and
scores whatever genome file is dropped on it.  The file goes from the browser to this
process, is written to a private temporary folder, scored, and deleted; the page loads no
outside resources.  If the program is killed while a file is being scored, the copy is
removed the next time it starts.

`pgc-risk gui --lan` also answers other devices on the local network, so a file can be
uploaded from a phone or another computer.  Because the file then crosses the network:

* other devices must use HTTPS (a self-signed certificate made on first use; its fingerprint
  is shown so it can be compared with what the browser reports).  Plain HTTP is answered
  only for this computer itself; anyone else is redirected to HTTPS;
* other devices must enter the access code shown on this computer before they get the page.
  Five wrong codes lock the form for a minute.

In both modes a request must name this server as its Host (which stops other web sites from
reaching it through DNS tricks) and the scoring API needs the random token that only the
served page knows.
"""
from __future__ import annotations

import ipaddress
import json
import os
import re
import secrets
import shutil
import socket
import ssl
import subprocess
import tempfile
import threading
import time
import webbrowser
from http.cookies import SimpleCookie
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote

from . import __version__
from .catalog import RISK, SOURCES
from .genome import GenomeFormatError

PAGE = Path(__file__).with_name("gui.html")
LOGIN_PAGE = Path(__file__).with_name("gui_login.html")
MAX_UPLOAD = 64 * 1024 ** 3
DEFAULT_PORT = 8765
TMP_PREFIX = "pgc-risk-gui-"
MAX_WRONG_CODES = 5
LOCK_SECONDS = 60
SESSION_SECONDS = 12 * 3600
CSP = ("default-src 'none'; style-src 'unsafe-inline'; script-src 'unsafe-inline'; connect-src 'self'; "
       "img-src data:; base-uri 'none'; form-action 'none'; frame-ancestors 'none'")


def sweep_stale_folders() -> int:
    """Delete upload folders left behind by an instance that was killed before it could tidy up.

    Each instance keeps a `lock` file in its folder: held open on Windows (so it cannot be
    deleted while the instance lives), holding the process id elsewhere.
    """
    removed = 0
    for folder in Path(tempfile.gettempdir()).glob(TMP_PREFIX + "*"):
        lock = folder / "lock"
        try:
            if os.name == "nt":
                lock.unlink()                       # PermissionError while its owner is running
            else:
                os.kill(int(lock.read_text()), 0)   # raises ProcessLookupError if the owner is gone
                continue
        except (PermissionError, ValueError):
            continue
        except OSError:
            pass                                    # no lock file, or its owner is gone: stale
        shutil.rmtree(folder, ignore_errors=True)
        removed += 1
    return removed


def local_addresses() -> list[str]:
    """This computer's IPv4 addresses on its networks (not loopback, not self-assigned)."""
    found: list[str] = []
    try:
        candidates = socket.gethostbyname_ex(socket.gethostname())[2]
    except OSError:
        candidates = []
    try:                                            # the address used for outgoing traffic
        with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as s:
            s.connect(("10.255.255.255", 1))
            candidates.append(s.getsockname()[0])
    except OSError:
        pass
    for a in candidates:
        ip = ipaddress.ip_address(a)
        if ip.version == 4 and not ip.is_loopback and not ip.is_link_local and a not in found:
            found.append(a)
    return found


_VIRTUAL_ADAPTER = re.compile(r"vethernet|virtual|vpn|wsl|hyper-v|loopback|bluetooth|docker|vmware|vbox|"
                              r"tailscale|wireguard|zerotier", re.I)


def adapter_names() -> dict[str, str]:
    """On Windows, the network adapter each IPv4 address belongs to ({} if that cannot be asked)."""
    if os.name != "nt":
        return {}
    try:
        out = subprocess.run(
            ["powershell", "-NoProfile", "-NonInteractive", "-Command",
             "Get-NetIPAddress -AddressFamily IPv4 | ForEach-Object { $_.IPAddress + '|' + $_.InterfaceAlias }"],
            capture_output=True, text=True, timeout=20, creationflags=subprocess.CREATE_NO_WINDOW).stdout
    except (OSError, subprocess.SubprocessError):
        return {}
    return dict(line.strip().split("|", 1) for line in out.splitlines() if "|" in line)


class App:
    """State shared by the request handlers: access secrets, the jobs and the temp folder."""

    def __init__(self, weights_dir: Path, panel_dir: Path | None = None, lan: bool = False,
                 trust_loopback: bool = True):
        self.weights_dir = Path(weights_dir)
        self.panel_dir = Path(panel_dir) if panel_dir else None
        self.lan = lan
        self.trust_loopback = trust_loopback         # False only in tests, to exercise the remote path
        self.token = secrets.token_urlsafe(24)
        self.code = f"{secrets.randbelow(10 ** 6):06d}" if lan else ""
        self.sessions: dict[str, float] = {}
        self.wrong_codes = 0
        self.locked_until = 0.0
        self.addresses = local_addresses() if lan else []
        self.adapters = adapter_names() if lan else {}
        self.host_name = socket.gethostname()
        self.tls = False
        self.fingerprint = ""
        self.port = 0
        sweep_stale_folders()
        self.tmp = Path(tempfile.mkdtemp(prefix=TMP_PREFIX))
        self._lock_file = open(self.tmp / "lock", "w")
        self._lock_file.write(str(os.getpid()))
        self._lock_file.flush()
        self.jobs: dict[str, dict] = {}
        self.lock = threading.Lock()
        self.one_at_a_time = threading.Lock()      # scoring a sequenced genome needs a few GB

    def close(self) -> None:
        self._lock_file.close()
        shutil.rmtree(self.tmp, ignore_errors=True)

    def uploads(self) -> list[Path]:
        """Uploaded files currently on disk (none, unless a file is being scored right now)."""
        return [p for p in self.tmp.iterdir() if p.name != "lock"]

    # ------------------------------------------------------------------ access
    def host_names(self) -> set[str]:
        names = {"127.0.0.1", "localhost"}
        if self.lan:
            names.update(self.addresses)
            names.update({self.host_name.lower(), self.host_name.lower() + ".local"})
        return names

    def check_code(self, code: str) -> tuple[bool, str]:
        """Compare an access code; returns (accepted, message for the person typing it)."""
        with self.lock:
            now = time.time()
            if now < self.locked_until:
                return False, f"Too many wrong codes. Try again in {int(self.locked_until - now) + 1} seconds."
            if self.code and secrets.compare_digest(code.strip(), self.code):
                self.wrong_codes = 0
                return True, ""
            self.wrong_codes += 1
            if self.wrong_codes >= MAX_WRONG_CODES:
                self.wrong_codes = 0
                self.locked_until = now + LOCK_SECONDS
                return False, f"Too many wrong codes. Try again in {LOCK_SECONDS} seconds."
            return False, "That is not the access code shown on the computer running pgc-risk."

    def new_session(self) -> str:
        sid = secrets.token_urlsafe(24)
        with self.lock:
            now = time.time()
            for old in [k for k, t in self.sessions.items() if t < now]:
                del self.sessions[old]
            self.sessions[sid] = now + SESSION_SECONDS
        return sid

    def session_ok(self, sid: str) -> bool:
        with self.lock:
            return bool(sid) and self.sessions.get(sid, 0) > time.time()

    def reachable_addresses(self) -> list[str]:
        """Addresses on real network adapters (Wi-Fi, Ethernet): the ones another device can use.

        VPN tunnels and virtual switches are left out when the adapters are known; the server
        still listens on them, but nothing on the local network can reach them.
        """
        real = [a for a in self.addresses
                if a in self.adapters and not _VIRTUAL_ADAPTER.search(self.adapters[a])]
        return real or self.addresses

    def other_adapters(self) -> list[str]:
        shown = set(self.reachable_addresses())
        return [f"{a} ({self.adapters.get(a, 'unknown adapter')})" for a in self.addresses if a not in shown]

    def vpn_adapters(self) -> list[str]:
        return sorted({n for n in self.adapters.values() if re.search("vpn|wireguard", n, re.I)})

    def urls(self) -> list[str]:
        scheme = "https" if self.tls else "http"
        return [f"{scheme}://{a}:{self.port}/" for a in self.reachable_addresses()]

    # ------------------------------------------------------------------ jobs
    def new_job(self, name: str) -> tuple[str, Path]:
        job_id = secrets.token_urlsafe(12)
        folder = self.tmp / job_id
        folder.mkdir()
        with self.lock:
            # results are kept only until the page has had time to collect them
            cutoff = time.time() - 3600
            for old in [k for k, j in self.jobs.items() if j["created"] < cutoff]:
                del self.jobs[old]
            self.jobs[job_id] = {"state": "uploading", "done": 0, "total": 0, "label": "",
                                 "created": time.time(), "name": name}
        return job_id, folder / name

    def update(self, job_id: str, **fields) -> None:
        with self.lock:
            if job_id in self.jobs:
                self.jobs[job_id].update(fields)

    def snapshot(self, job_id: str) -> dict | None:
        with self.lock:
            job = self.jobs.get(job_id)
            return dict(job) if job else None

    def score(self, job_id: str, path: Path) -> None:
        """Runs in a worker thread: score the uploaded file and delete it, then publish the outcome."""
        from .report import render_text, run

        name = path.name
        try:
            self.update(job_id, state="queued", label="Waiting for the previous file to finish")
            with self.one_at_a_time:
                rep = run(path, self.weights_dir, panel_dir=self.panel_dir,
                          progress=lambda stage, done, total, label: self.update(
                              job_id, state=stage, done=done, total=total, label=label))
            result = rep.to_dict()
            result["genome"]["file"] = name
            result["text"] = render_text(rep)
            result["sources"] = {
                s.key: {"citation": s.citation, "ref": s.ref, "notes": s.notes,
                        "prevalence_source": RISK[s.key].prevalence_source,
                        "accuracy_source": RISK[s.key].extra.get("source", "")}
                for s in SOURCES if s.key in RISK}
            outcome = {"state": "done", "result": result}
        except (GenomeFormatError, FileNotFoundError) as e:
            outcome = {"state": "error", "error": str(e)}
        except Exception as e:                      # a bug should reach the page, not vanish in a thread
            outcome = {"state": "error", "error": f"Scoring failed: {type(e).__name__}: {e}"}
        # the copy is gone before the page is told the result
        if not self.remove(path.parent) and outcome["state"] == "done":
            outcome["result"]["warnings"].append(
                f"The temporary copy of the file could not be deleted ({path.parent}). It will be "
                "removed the next time pgc-risk starts; you can also delete that folder yourself.")
        self.update(job_id, **outcome)

    @staticmethod
    def remove(folder: Path) -> bool:
        for _ in range(10):
            shutil.rmtree(folder, ignore_errors=True)
            if not folder.exists():
                return True
            time.sleep(0.2)                         # e.g. a virus scanner still holding the file
        return False


def safe_name(raw: str) -> str:
    """The uploaded file's own name, reduced to something safe to create on disk."""
    name = unquote(raw or "").replace("\\", "/").rsplit("/", 1)[-1]
    name = re.sub(r"[^A-Za-z0-9._ ()+-]", "_", name).strip(" .")
    return name[:120] or "genome"


def _escape(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


class Server(ThreadingHTTPServer):
    daemon_threads = True
    tls_context: ssl.SSLContext | None = None

    def handle_error(self, request, client_address):
        # a browser refusing the self-signed certificate, a port scan, a dropped phone: not our errors
        pass


def make_handler(app: App):
    class Handler(BaseHTTPRequestHandler):
        server_version = f"pgc-risk/{__version__}"
        protocol_version = "HTTP/1.1"
        timeout = 300                               # a stalled connection does not hold a thread forever

        def log_message(self, fmt, *args):         # keep the console quiet; file names are private
            pass

        def setup(self):
            # With TLS enabled the same port answers both protocols: a TLS hello starts with 0x16.
            self.encrypted = False
            ctx = self.server.tls_context
            if ctx is not None:
                self.request.settimeout(15)
                first = self.request.recv(1, socket.MSG_PEEK)
                if first == b"\x16":
                    self.request = ctx.wrap_socket(self.request, server_side=True)
                    self.encrypted = True
            super().setup()

        # ------------------------------------------------------------------ helpers
        def _send(self, status: int, body: bytes, ctype: str, extra: dict | None = None) -> None:
            self.send_response(status)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.send_header("Referrer-Policy", "no-referrer")
            self.send_header("Content-Security-Policy", CSP)   # the page may talk to this server only
            for k, v in (extra or {}).items():
                self.send_header(k, v)
            self.end_headers()
            self.wfile.write(body)

        def _json(self, status: int, obj: dict, extra: dict | None = None) -> None:
            self._send(status, json.dumps(obj).encode("utf-8"), "application/json; charset=utf-8", extra)

        def _drain(self) -> None:
            """Before refusing a request, take a small body off the wire.

            Closing with unread data makes the system reset the connection, and the client may
            lose the answer.  Large bodies are not read: the connection is simply closed.
            """
            self.close_connection = True
            try:
                left = int(self.headers.get("Content-Length") or 0)
                if 0 < left <= 1 << 20:
                    self.connection.settimeout(3)
                    while left > 0:
                        chunk = self.rfile.read(min(65536, left))
                        if not chunk:
                            break
                        left -= len(chunk)
            except (ValueError, OSError):
                pass

        def _refuse(self, status: int, message: str) -> None:
            self._drain()
            self._json(status, {"error": message})

        def _host(self) -> str | None:
            """The Host header's name if it is one of ours on our port, else None."""
            name, sep, port = self.headers.get("Host", "").rpartition(":")
            if not sep or port != str(self.server.server_address[1]):
                return None
            return name.lower() if name.lower() in app.host_names() else None

        def _this_computer(self) -> bool:
            return app.trust_loopback and ipaddress.ip_address(self.client_address[0]).is_loopback

        def _signed_in(self) -> bool:
            if self._this_computer():
                return True
            cookie = SimpleCookie(self.headers.get("Cookie", ""))
            return app.lan and "pgc_session" in cookie and app.session_ok(cookie["pgc_session"].value)

        def _gate(self) -> bool:
            """Checks every request must pass.  Sends the refusal itself; True means carry on."""
            host = self._host()
            if host is None:
                self._refuse(403, "this server does not answer to that host name")
                return False
            if not self._this_computer():
                if not app.lan:
                    self._refuse(403, "this server only answers the computer it runs on")
                    return False
                if self.server.tls_context is not None and not self.encrypted:
                    # another device on plain HTTP: send it to HTTPS before anything is exchanged
                    self._drain()
                    self._send(301, b"", "text/plain", {"Location": f"https://{self.headers['Host']}{self.path}"})
                    return False
            return True

        def _page(self) -> bytes:
            html = PAGE.read_text(encoding="utf-8")
            if self._this_computer():
                privacy = ("This page is served by the program running on this computer (address 127.0.0.1). "
                           "Your file is handed to that program, scored, and removed from its temporary folder "
                           "as soon as scoring ends. Nothing is sent over the internet, and this page loads "
                           "nothing from it.")
            else:
                how = ("over an encrypted connection" if self.encrypted else
                       "without encryption, so other devices on the same network could read it on the way")
                privacy = (f"This page comes from pgc-risk running on the computer {app.host_name} on your "
                           f"network. Your file travels from this device to that computer {how}. It is scored "
                           "there and its temporary copy is removed as soon as scoring ends. Nothing is sent "
                           "over the internet, and this page loads nothing from it.")
            share = ""
            if app.lan and self._this_computer():
                links = "".join(f"<li><code>{_escape(u)}</code></li>" for u in app.urls()) or \
                    "<li>No network address found - is this computer connected to a network?</li>"
                cert = (f"<p>The first visit shows a certificate warning, because the certificate was made by "
                        f"this computer rather than bought from an authority. To be sure the device is talking "
                        f"to this computer, compare the SHA-256 fingerprint the browser shows with "
                        f"<code>{_escape(app.fingerprint)}</code>, then continue.</p>" if app.tls else
                        "<p><b>The connection is not encrypted.</b> Other devices on the same network could "
                        "read the file on the way.</p>")
                vpn = app.vpn_adapters()
                blockers = (f"{_escape(' / '.join(vpn))} is active on this computer: if the other device cannot "
                            f"connect, switch on its \"allow LAN connections\" setting or pause it. " if vpn else "")
                share = (f'<div class="card share"><h2>Upload from another device</h2>'
                         f"<p>On a phone or computer on the same network, open:</p><ul>{links}</ul>"
                         f'<p>and enter the access code <b class="code">{_escape(app.code)}</b></p>{cert}'
                         f"<p class=\"sub\">{blockers}Windows Firewall must allow Python to accept connections, "
                         f"and some Wi-Fi networks (guest and public ones) keep devices from reaching each "
                         f"other at all.</p></div>")
            return (html.replace("__TOKEN__", app.token).replace("__VERSION__", __version__)
                        .replace("__PRIVACY__", _escape(privacy)).replace("__SHARE__", share)).encode("utf-8")

        # ------------------------------------------------------------------ routes
        def do_GET(self):
            if not self._gate():
                return
            path = self.path.split("?", 1)[0]
            if path == "/":
                if not self._signed_in():
                    login = LOGIN_PAGE.read_text(encoding="utf-8").replace("__HOST__", _escape(app.host_name))
                    return self._send(200, login.encode("utf-8"), "text/html; charset=utf-8")
                return self._send(200, self._page(), "text/html; charset=utf-8")
            if path == "/api/session":            # lets the page ask "am I still signed in?" before an upload
                if not self._signed_in():
                    return self._json(401, {"error": "sign in again"})
                if not secrets.compare_digest(self.headers.get("X-Token", ""), app.token):
                    return self._json(403, {"error": "missing or wrong token"})
                return self._json(200, {"ok": True})
            if path.startswith("/api/job/"):
                if not self._signed_in():
                    return self._json(401, {"error": "sign in again"})
                if not secrets.compare_digest(self.headers.get("X-Token", ""), app.token):
                    return self._json(403, {"error": "missing or wrong token"})
                job = app.snapshot(path.rsplit("/", 1)[-1])
                if job is None:
                    return self._json(404, {"error": "no such job"})
                job.pop("created", None)
                return self._json(200, job)
            return self._json(404, {"error": "not found"})

        def do_POST(self):
            if not self._gate():
                return
            path = self.path.split("?", 1)[0]
            try:
                length = int(self.headers.get("Content-Length", ""))
            except ValueError:
                return self._refuse(411, "Content-Length required")
            if path == "/api/login":
                if not app.lan or length > 200:
                    return self._refuse(403, "not allowed")
                try:
                    code = str(json.loads(self.rfile.read(length)).get("code", ""))
                except (ValueError, AttributeError):
                    code = ""
                ok, message = app.check_code(code)
                if not ok:
                    return self._json(403, {"error": message})
                cookie = (f"pgc_session={app.new_session()}; Path=/; HttpOnly; SameSite=Strict; "
                          f"Max-Age={SESSION_SECONDS}" + ("; Secure" if self.encrypted else ""))
                return self._json(200, {"ok": True}, {"Set-Cookie": cookie})
            if path != "/api/score":
                return self._refuse(404, "not found")
            if not self._signed_in():
                return self._refuse(401, "sign in again")
            if not secrets.compare_digest(self.headers.get("X-Token", ""), app.token):
                return self._refuse(403, "missing or wrong token")
            if length <= 0 or length > MAX_UPLOAD:
                return self._refuse(413, "empty file, or larger than 64 GB")
            job_id, dest = app.new_job(safe_name(self.headers.get("X-Filename", "")))
            try:
                remaining = length
                with open(dest, "wb") as out:
                    while remaining:
                        chunk = self.rfile.read(min(1 << 20, remaining))
                        if not chunk:
                            raise ConnectionError("upload interrupted")
                        out.write(chunk)
                        remaining -= len(chunk)
            except OSError as e:
                shutil.rmtree(dest.parent, ignore_errors=True)
                app.update(job_id, state="error", error=f"Upload failed: {e}")
                return self._refuse(400, "upload failed")
            threading.Thread(target=app.score, args=(job_id, dest), daemon=True).start()
            return self._json(200, {"job": job_id})

    return Handler


def ensure_certificate(directory: Path, names: list[str], addresses: list[str]) -> tuple[Path, Path, str]:
    """A self-signed certificate for this computer's names and addresses; made once and reused.

    Returns (certificate file, key file, SHA-256 fingerprint).  Needs the `cryptography` package.
    """
    import datetime

    from cryptography import x509
    from cryptography.hazmat.primitives import hashes, serialization
    from cryptography.hazmat.primitives.asymmetric import ec
    from cryptography.x509.oid import NameOID

    directory.mkdir(parents=True, exist_ok=True)
    cert_file, key_file = directory / "cert.pem", directory / "key.pem"
    now = datetime.datetime.now(datetime.timezone.utc)
    wanted_ips = {ipaddress.ip_address(a) for a in addresses} | {ipaddress.ip_address("127.0.0.1")}
    wanted_names = set(names) | {"localhost"}
    if cert_file.exists() and key_file.exists():
        try:
            cert = x509.load_pem_x509_certificate(cert_file.read_bytes())
            san = cert.extensions.get_extension_for_class(x509.SubjectAlternativeName).value
            covered = (wanted_ips <= set(san.get_values_for_type(x509.IPAddress))
                       and wanted_names <= set(san.get_values_for_type(x509.DNSName)))
            if covered and cert.not_valid_after_utc > now + datetime.timedelta(days=14):
                return cert_file, key_file, _fingerprint(cert.fingerprint(hashes.SHA256()))
        except Exception:
            pass                                    # unreadable or incomplete: make a new one
    key = ec.generate_private_key(ec.SECP256R1())
    subject = x509.Name([x509.NameAttribute(NameOID.COMMON_NAME, f"pgc-risk on {names[0]}")])
    san = x509.SubjectAlternativeName([x509.DNSName(n) for n in sorted(wanted_names)]
                                      + [x509.IPAddress(a) for a in sorted(wanted_ips, key=str)])
    cert = (x509.CertificateBuilder().subject_name(subject).issuer_name(subject).public_key(key.public_key())
            .serial_number(x509.random_serial_number())
            .not_valid_before(now - datetime.timedelta(minutes=5))
            .not_valid_after(now + datetime.timedelta(days=365))
            .add_extension(san, critical=False)
            .add_extension(x509.BasicConstraints(ca=False, path_length=None), critical=True)
            .add_extension(x509.ExtendedKeyUsage([x509.oid.ExtendedKeyUsageOID.SERVER_AUTH]), critical=False)
            .sign(key, hashes.SHA256()))
    key_file.write_bytes(key.private_bytes(serialization.Encoding.PEM, serialization.PrivateFormat.PKCS8,
                                           serialization.NoEncryption()))
    cert_file.write_bytes(cert.public_bytes(serialization.Encoding.PEM))
    return cert_file, key_file, _fingerprint(cert.fingerprint(hashes.SHA256()))


def _fingerprint(digest: bytes) -> str:
    return ":".join(f"{b:02X}" for b in digest)


def make_server(weights_dir: Path, panel_dir: Path | None = None, port: int = DEFAULT_PORT,
                lan: bool = False, cert_dir: Path | None = None, trust_loopback: bool = True):
    """Bind to this computer only, or to every interface with `lan`.  With `cert_dir`, serve HTTPS too.

    If the port is taken, any free one is used.
    """
    app = App(weights_dir, panel_dir, lan=lan, trust_loopback=trust_loopback)
    bind = "0.0.0.0" if lan else "127.0.0.1"
    try:
        server = Server((bind, port), make_handler(app))
    except OSError:
        server = Server((bind, 0), make_handler(app))
    app.port = server.server_address[1]
    if lan and cert_dir is not None:
        names = [app.host_name.lower(), app.host_name.lower() + ".local"]
        cert, key, app.fingerprint = ensure_certificate(Path(cert_dir), names, app.addresses)
        ctx = ssl.SSLContext(ssl.PROTOCOL_TLS_SERVER)
        ctx.minimum_version = ssl.TLSVersion.TLSv1_2
        ctx.load_cert_chain(cert, key)
        server.tls_context = ctx
        app.tls = True
    return server, app


def serve(weights_dir: Path, port: int = DEFAULT_PORT, open_browser: bool = True, lan: bool = False,
          encrypt: bool = True) -> int:
    from .panel import Panel

    weights_dir = Path(weights_dir)
    if not any(weights_dir.glob("*.npz")) or not Panel.exists(weights_dir.parent / "panel"):
        print(f"error: no weight files or LD panel under {weights_dir.parent} - run pgc-risk fetch "
              "and pgc-risk build first")
        return 2
    cert_dir = None
    if lan and encrypt:
        try:
            import cryptography  # noqa: F401
            cert_dir = weights_dir.parent / "gui-cert"
        except ImportError:
            print("error: network mode encrypts the connection and needs the 'cryptography' package:\n"
                  "         uv pip install --python .venv\\Scripts\\python.exe cryptography\n"
                  "       (or add --no-encryption to serve plain HTTP on the network)")
            return 2
    server, app = make_server(weights_dir, port=port, lan=lan, cert_dir=cert_dir)
    local = f"http://127.0.0.1:{app.port}/"
    say = lambda text="": print(text, flush=True)
    say(f"pgc-risk {__version__} is running.")
    say(f"  On this computer:    {local}")
    if lan:
        for url in app.urls() or ["(no network address found)"]:
            adapter = app.adapters.get(url.split("//")[1].split(":")[0], "")
            say(f"  From other devices:  {url}" + (f"   ({adapter})" if adapter else ""))
        say(f"  Access code:         {app.code}")
        if app.tls:
            say(f"  Certificate SHA-256: {app.fingerprint}")
            say("  Other devices see a certificate warning on the first visit; compare the fingerprint, then continue.")
        else:
            say("  WARNING: not encrypted - other devices on this network could read an uploaded file on the way.")
        if app.other_adapters():
            say(f"  (also listening on {', '.join(app.other_adapters())})")
        say("  If another device cannot connect: Windows Firewall must allow Python, a VPN on this computer")
        say("  must allow LAN connections, and the Wi-Fi must let its devices reach each other.")
    else:
        say("  It answers this computer only (use --lan to upload from other devices).")
    say("Close this window (or press Ctrl+C) to stop it.")
    if open_browser:
        threading.Timer(0.4, webbrowser.open, args=(local,)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
        app.close()
    return 0
