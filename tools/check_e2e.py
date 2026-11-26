"""End-to-end check: score the test genome files and compare with the direct matrix computation.

tools/calibrate.py computes every 1000 Genomes sample's z straight from the genotype matrix.
Here the same people go through the real file parsers and matching code, in three formats;
the standardised scores must agree.

usage: python tools/check_e2e.py   (after tools/calibrate.py and tools/make_test_genomes.py)
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pgc_risk import panel as pn  # noqa: E402
from pgc_risk.report import run  # noqa: E402

GENOMES = ROOT / "data" / "validation" / "genomes"


def main() -> None:
    ref = np.load(ROOT / "data" / "validation" / "calibration_scores.npz")
    names = ref["samples"].tolist()
    info = pn.samples()
    worst = 0.0
    print(f"{'sample':<9}{'pop':<5}{'file':<14}{'secs':>6}{'fit':>7}{'ref fit':>9}{'max |dz|':>10}   notes")
    for path in sorted(GENOMES.glob("*")):
        sid, kind = path.name.split(".", 1)
        cfg = "array" if kind.startswith("23andme") else "all"
        t0 = time.time()
        rep = run(path, pn.WEIGHTS)
        secs = time.time() - t0
        i = names.index(sid)
        dz = max(abs(r.z - float(ref[f"{r.key}_{cfg}"][i])) for r in rep.results)
        worst = max(worst, dz)
        flags = "; ".join(w.split(":")[0] for w in rep.warnings) or "-"
        print(f"{sid:<9}{info[sid][0]:<5}{kind:<14}{secs:>6.1f}{rep.fit:>7.3f}{float(ref['fit'][i]):>9.3f}"
              f"{dz:>10.4f}   {flags}")
    print(f"\nlargest difference in z between file-based and matrix-based scoring: {worst:.4f}")


if __name__ == "__main__":
    main()
