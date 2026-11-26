"""Network mode of the web interface: HTTPS, the access code, sessions and host checks.

The test client connects over loopback, which the server normally treats as "this computer".
`trust_loopback=False` switches that off so the client is handled like another device.
"""
import http.client
import json
import ssl
import threading
import time

import numpy as np
import pytest

from pgc_risk import gui, report
from test_panel import N, _make_project, _write_23andme

pytest.importorskip("cryptography")


def start(tmp_path, monkeypatch, **kw):
    rng = np.random.default_rng(21)
    wdir, panel_dir, g, chrom, pos = _make_project(tmp_path, rng)
    monkeypatch.setattr(report, "FIT_GOOD", 1.25)        # thresholds assume ~500,000 SNPs
    monkeypatch.setattr(report, "FIT_POOR", 1.3)
    genome = tmp_path / "genome.txt"
    _write_23andme(genome, chrom, pos, g[:, N // 3])
    server, app = gui.make_server(wdir, panel_dir, port=0, **kw)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return server, app, genome


@pytest.fixture
def remote(tmp_path, monkeypatch):
    """Network mode with HTTPS, seen from 'another device'."""
    server, app, genome = start(tmp_path, monkeypatch, lan=True, cert_dir=tmp_path / "cert",
                                trust_loopback=False)
    yield app, genome
    server.shutdown()
    server.server_close()
    app.close()


def request(app, method, path, body=None, headers=None, tls=True, host=None):
    """One request, no redirect following, certificate not verified (it is self-signed)."""
    if tls:
        conn = http.client.HTTPSConnection("127.0.0.1", app.port, timeout=30,
                                           context=ssl._create_unverified_context())
    else:
        conn = http.client.HTTPConnection("127.0.0.1", app.port, timeout=30)
    h = dict(headers or {})
    if host:
        h["Host"] = host
    conn.request(method, path, body=body, headers=h)
    r = conn.getresponse()
    out = (r.status, r.read(), {k.lower(): v for k, v in r.getheaders()})
    conn.close()
    return out


def sign_in(app):
    status, _, headers = request(app, "POST", "/api/login", json.dumps({"code": app.code}),
                                 {"Content-Type": "application/json"})
    assert status == 200
    return headers["set-cookie"].split(";")[0]


def test_other_devices_are_sent_to_https(remote):
    app, _ = remote
    status, body, headers = request(app, "GET", "/", tls=False)
    assert status == 301 and headers["location"] == f"https://127.0.0.1:{app.port}/"
    assert app.token.encode() not in body
    # an upload attempted over plain HTTP is turned away before its body is read
    status, _, headers = request(app, "POST", "/api/score", b"ACGT" * 10, {"X-Token": app.token}, tls=False)
    assert status == 301 and app.uploads() == []


def test_page_needs_the_access_code(remote):
    app, _ = remote
    status, body, _ = request(app, "GET", "/")
    page = body.decode()
    assert status == 200 and "Access code" in page
    assert app.token not in page and app.code not in page
    assert request(app, "GET", "/api/job/x", headers={"X-Token": app.token})[0] == 401
    assert request(app, "POST", "/api/score", b"x", {"X-Token": app.token})[0] == 401
    assert app.uploads() == []


def test_wrong_codes_lock_the_form(remote):
    app, _ = remote
    post = lambda code: request(app, "POST", "/api/login", json.dumps({"code": code}),
                                {"Content-Type": "application/json"})
    for _ in range(gui.MAX_WRONG_CODES - 1):
        status, body, headers = post("000000" if app.code != "000000" else "111111")
        assert status == 403 and "set-cookie" not in headers and b"not the access code" in body
    status, body, _ = post("xxxxxx")
    assert status == 403 and b"Too many wrong codes" in body
    status, body, headers = post(app.code)                # the right code does not help while locked
    assert status == 403 and "set-cookie" not in headers
    app.locked_until = time.time() - 1                    # the minute passes
    status, _, headers = post(app.code)
    assert status == 200
    cookie = headers["set-cookie"]
    assert "HttpOnly" in cookie and "SameSite=Strict" in cookie and "Secure" in cookie
    assert request(app, "POST", "/api/login", b"x" * 500)[0] == 403        # oversized body refused


def test_signed_in_device_gets_the_page_and_can_score(remote):
    app, genome = remote
    cookie = sign_in(app)
    status, body, _ = request(app, "GET", "/", headers={"Cookie": cookie})
    page = body.decode()
    assert status == 200 and app.token in page
    assert "over an encrypted connection" in page and "__PRIVACY__" not in page and "__SHARE__" not in page
    assert app.code not in page                           # the code is shown on the host computer only
    # the cookie alone is not enough for the API, nor is the token alone
    assert request(app, "POST", "/api/score", b"x", {"Cookie": cookie})[0] == 403
    assert request(app, "POST", "/api/score", b"x", {"X-Token": app.token})[0] == 401
    auth = {"Cookie": cookie, "X-Token": app.token}
    status, body, _ = request(app, "POST", "/api/score", genome.read_bytes(), dict(auth, **{"X-Filename": "me.txt"}))
    assert status == 200
    job = json.loads(body)["job"]
    for _ in range(200):
        state = json.loads(request(app, "GET", f"/api/job/{job}", headers=auth)[1])
        if state["state"] in ("done", "error"):
            break
        time.sleep(0.05)
    assert state["state"] == "done", state
    direct = report.run(genome, app.weights_dir, panel_dir=app.panel_dir).results[0]
    assert state["result"]["disorders"][0]["z"] == pytest.approx(direct.z)
    assert app.uploads() == []
    assert request(app, "GET", f"/api/job/{job}", headers={"X-Token": app.token})[0] == 401   # no cookie


def test_host_names(remote):
    app, _ = remote
    assert request(app, "GET", "/", host=f"evil.example:{app.port}")[0] == 403
    assert request(app, "GET", "/", host="127.0.0.1:1")[0] == 403           # right name, wrong port
    assert request(app, "GET", "/", host=f"{app.host_name}:{app.port}")[0] == 200
    for address in app.addresses:
        assert request(app, "GET", "/", host=f"{address}:{app.port}")[0] == 200


def test_certificate_is_reused_and_extended(tmp_path):
    from cryptography import x509
    cert, key, fp = gui.ensure_certificate(tmp_path, ["host", "host.local"], ["192.168.1.20"])
    san = x509.load_pem_x509_certificate(cert.read_bytes()).extensions.get_extension_for_class(
        x509.SubjectAlternativeName).value
    assert {str(a) for a in san.get_values_for_type(x509.IPAddress)} == {"127.0.0.1", "192.168.1.20"}
    assert set(san.get_values_for_type(x509.DNSName)) == {"host", "host.local", "localhost"}
    assert len(fp.split(":")) == 32 and "PRIVATE KEY" in key.read_text()
    assert gui.ensure_certificate(tmp_path, ["host", "host.local"], ["192.168.1.20"])[2] == fp     # reused
    assert gui.ensure_certificate(tmp_path, ["host", "host.local"], ["10.0.0.7"])[2] != fp         # new address


def test_host_computer_sees_how_to_connect_and_needs_no_code(tmp_path, monkeypatch):
    server, app, genome = start(tmp_path, monkeypatch, lan=True, cert_dir=tmp_path / "cert")
    try:
        status, body, _ = request(app, "GET", "/", tls=False)             # plain HTTP is fine on the host itself
        page = body.decode()
        assert status == 200 and app.token in page
        assert "Upload from another device" in page and app.code in page and app.fingerprint in page
        for url in app.urls():
            assert url in page and url.startswith("https://")
        assert "served by the program running on this computer" in page
    finally:
        server.shutdown()
        server.server_close()
        app.close()


def test_local_only_mode_refuses_other_devices(tmp_path, monkeypatch):
    server, app, _ = start(tmp_path, monkeypatch, trust_loopback=False)    # no lan=True
    try:
        assert server.server_address[0] == "127.0.0.1"                      # not even bound to the network
        status, body, _ = request(app, "GET", "/", tls=False)
        assert status == 403 and app.token.encode() not in body
        assert request(app, "POST", "/api/login", json.dumps({"code": ""}), tls=False)[0] == 403
    finally:
        server.shutdown()
        server.server_close()
        app.close()
