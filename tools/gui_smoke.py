"""Drive the running GUI server the way the page does: upload a file, follow progress, print the result.

Also writes a stand-alone copy of the page with the result already rendered (for a screenshot).

usage: python tools/gui_smoke.py <genome file> [port] [out.html]
"""
from __future__ import annotations

import json
import re
import sys
import time
import urllib.request
from pathlib import Path
from urllib.parse import quote


def main() -> None:
    path = Path(sys.argv[1])
    port = int(sys.argv[2]) if len(sys.argv) > 2 else 8765
    base = f"http://127.0.0.1:{port}"
    page = urllib.request.urlopen(base + "/").read().decode()
    token = re.search(r'const TOKEN = "([^"]+)"', page).group(1)
    t0 = time.time()
    req = urllib.request.Request(base + "/api/score", data=path.read_bytes(),
                                 headers={"X-Token": token, "X-Filename": quote(path.name)})
    job = json.loads(urllib.request.urlopen(req).read())["job"]
    seen = []
    while True:
        req = urllib.request.Request(f"{base}/api/job/{job}", headers={"X-Token": token})
        state = json.loads(urllib.request.urlopen(req).read())
        tag = (state["state"], state.get("label"))
        if tag not in seen:
            seen.append(tag)
        if state["state"] in ("done", "error"):
            break
        time.sleep(0.2)
    print(f"{path.name}: {state['state']} in {time.time() - t0:.1f}s; stages seen: {len(seen)}")
    if state["state"] == "error":
        print("error:", state["error"])
        return
    res = state["result"]
    print("valid:", res["valid"], "| fit:", round(res["ancestry_fit"]["value"], 3), "| warnings:", len(res["warnings"]))
    for d in res["disorders"]:
        print(f"  {d['key']:<5} z={d['z']:+.4f}  p={d['probability']:.4f}")
    if len(sys.argv) > 3:
        static = page.replace("</body>", "<script>render(" + json.dumps(res) + ");</script></body>")
        Path(sys.argv[3]).write_text(static, encoding="utf-8")
        print("wrote", sys.argv[3])


if __name__ == "__main__":
    main()
