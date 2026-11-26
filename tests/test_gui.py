"""The local web interface: upload, progress, result, clean-up and the access checks."""
import json
import threading
import time
import urllib.error
import urllib.request

import numpy as np
import pytest

from pgc_risk import gui, report
from test_panel import N, _make_project, _write_23andme


@pytest.fixture
def running(tmp_path, monkeypatch):
    """A server on a free port over a small synthetic weight set, plus one person's genome file."""
    rng = np.random.default_rng(11)
    wdir, panel_dir, g, chrom, pos = _make_project(tmp_path, rng)
    # the real thresholds assume ~500,000 SNPs; with 3,200 the ancestry statistic is far noisier
    monkeypatch.setattr(report, "FIT_GOOD", 1.25)
    monkeypatch.setattr(report, "FIT_POOR", 1.3)
    genome = tmp_path / "my genome.txt"
    _write_23andme(genome, chrom, pos, g[:, N // 2])
    server, app = gui.make_server(wdir, panel_dir, port=0)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{server.server_address[1]}"
    yield base, app, genome
    server.shutdown()
    server.server_close()
    app.close()


def call(url, data=None, headers=None):
    req = urllib.request.Request(url, data=data, headers=headers or {})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return r.status, r.read(), dict(r.headers)
    except urllib.error.HTTPError as e:
        return e.code, e.read(), dict(e.headers)


def wait_for(base, app, job):
    for _ in range(200):
        _, body, _ = call(f"{base}/api/job/{job}", headers={"X-Token": app.token})
        state = json.loads(body)
        if state["state"] in ("done", "error"):
            return state
        time.sleep(0.05)
    raise AssertionError("job did not finish")


def test_page_is_self_contained(running):
    base, app, _ = running
    status, body, headers = call(base + "/")
    page = body.decode()
    assert status == 200 and app.token in page and "__TOKEN__" not in page
    for outside in ("http://", "https://", "//cdn", "@import", "<link"):
        assert outside not in page, outside            # nothing is fetched from anywhere else
    assert "default-src 'none'" in headers["Content-Security-Policy"]
    assert "connect-src 'self'" in headers["Content-Security-Policy"]


def test_upload_scores_and_deletes_the_file(running):
    base, app, genome = running
    status, body, _ = call(base + "/api/score", data=genome.read_bytes(),
                           headers={"X-Token": app.token, "X-Filename": "my%20genome.txt"})
    assert status == 200
    state = wait_for(base, app, json.loads(body)["job"])
    assert state["state"] == "done", state
    res = state["result"]
    assert res["valid"] and res["genome"]["file"] == "my genome.txt" and res["genome"]["format"] == "23andme"
    d = res["disorders"][0]
    assert d["key"] == "scz" and 0 < d["probability"] < 1 and d["prevalence"] == 0.01
    assert "Schizophrenia" in res["text"] and "scz" in res["sources"]
    # same answer as scoring the file directly
    direct = report.run(genome, app.weights_dir, panel_dir=app.panel_dir).results[0]
    assert d["z"] == pytest.approx(direct.z)
    assert app.uploads() == []                          # the uploaded copy is gone


def test_unreadable_file_reports_an_error_and_is_deleted(running):
    base, app, _ = running
    _, body, _ = call(base + "/api/score", data=b"this is not a genome\nat all\n",
                      headers={"X-Token": app.token, "X-Filename": "..%2F..%2Fnotes.txt"})
    state = wait_for(base, app, json.loads(body)["job"])
    assert state["state"] == "error" and state["error"].startswith("notes.txt: no autosomal SNP genotypes found")
    assert app.uploads() == []


def test_requests_without_the_token_or_from_another_host_are_refused(running):
    base, app, genome = running
    assert call(base + "/api/score", data=b"x", headers={"X-Filename": "a.txt"})[0] == 403
    assert call(base + "/api/score", data=b"x", headers={"X-Token": "wrong"})[0] == 403
    assert call(base + "/api/job/anything")[0] == 403
    assert call(base + "/api/job/anything", headers={"X-Token": app.token})[0] == 404
    # a web page that reached this port under another host name (DNS rebinding) gets nothing
    assert call(base + "/", headers={"Host": "evil.example"})[0] == 403
    assert call(base + "/api/score", data=b"x", headers={"X-Token": app.token, "Host": "evil.example"})[0] == 403
    assert app.uploads() == []


def test_safe_name():
    assert gui.safe_name("genome_Jane_Doe_v5_Full_20240101.zip") == "genome_Jane_Doe_v5_Full_20240101.zip"
    assert gui.safe_name("..%2F..%2Fetc%2Fpasswd") == "passwd"
    assert gui.safe_name("C%3A%5CUsers%5Cme%5Cmy%20file.vcf.gz") == "my file.vcf.gz"
    assert gui.safe_name("") == "genome" and gui.safe_name("...") == "genome"
    assert gui.safe_name("a<b>:c|d?.txt") == "a_b__c_d_.txt"


def test_leftovers_of_a_killed_instance_are_swept_but_a_live_one_is_left_alone(tmp_path, monkeypatch):
    import tempfile
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))
    live = gui.App(tmp_path / "weights")                 # holds its lock
    (live.tmp / "job").mkdir()
    (live.tmp / "job" / "in progress.txt").write_text("genome being scored")
    dead = tmp_path / (gui.TMP_PREFIX + "killed")
    (dead / "job").mkdir(parents=True)
    (dead / "job" / "genome.txt").write_text("left behind")
    (dead / "lock").write_text("999999999")              # nobody holds this one
    unrelated = tmp_path / "pgc-risk-something-else"
    unrelated.mkdir()
    assert gui.sweep_stale_folders() == 1
    assert not dead.exists()
    assert (live.tmp / "job" / "in progress.txt").exists() and unrelated.exists()
    live.close()
    assert not live.tmp.exists()


def test_zip_upload_is_scored_and_its_copy_is_deleted(running, tmp_path):
    import zipfile
    base, app, genome = running
    archive = tmp_path / "genome_export.zip"
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.write(genome, "genome_Jane_Doe_v5_Full.txt")
    _, body, _ = call(base + "/api/score", data=archive.read_bytes(),
                      headers={"X-Token": app.token, "X-Filename": "genome_export.zip"})
    state = wait_for(base, app, json.loads(body)["job"])
    assert state["state"] == "done" and state["result"]["warnings"] == []
    assert app.uploads() == []            # an archive left open would keep Windows from deleting it


def test_session_endpoint(running):
    base, app, _ = running
    assert call(base + "/api/session", headers={"X-Token": app.token})[0] == 200
    assert call(base + "/api/session")[0] == 403
