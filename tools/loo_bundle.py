"""Hold-out accuracy of the bundle's selection rule (best available SNP per pre-computed clump).

Same designs and measurement as tools/loo.py, but the SNPs are chosen the way a bundle chooses
them: the training GWAS is clumped once over all its candidates, and a genome that has only
array SNPs is represented, clump by clump, by its most significant available SNP (among the
48 most significant of the clump).  Printed next to the full method's result on the same
array SNPs.

usage: python tools/loo_bundle.py [MDD|OCD|SCZ]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import loo  # noqa: E402

from pgc_risk import build  # noqa: E402
from pgc_risk.panel import Panel  # noqa: E402
from pgc_risk.score import ld_clump, ld_inflation, load_weights  # noqa: E402

CAP = 48


def measure(w, tb, ts, sel, panel, k):
    infl = ld_inflation(w, sel, panel)
    b, t, s, f = w.beta[sel], tb[sel], ts[sel], w.freq[sel]
    a = (b * t / s ** 2).sum() / (b * b / s ** 2).sum()
    r2_obs = a * a * (2 * f * (1 - f) * b * b).sum() / 4 / infl
    return len(sel), loo.liability_r2(r2_obs, k) * (1 if a > 0 else -1)


def main(only) -> None:
    con = build._connect()
    build._ensure_freq(con, lambda *a: None)
    panel = Panel()
    array_rs = np.load(loo.ROOT / "data" / "validation" / "array_like_rsids.npy")
    print(f"{'design':<40}{'full method':>22}{'bundle rule':>22}")
    for label, train, target, k in loo.DESIGNS:
        if only and not any(label.upper().startswith(o.upper()) for o in only):
            continue
        w = load_weights(train.key, loo.OUT)
        stats = loo.target_stats(con, train, target)
        tb, ts = np.full(w.n, np.nan), np.full(w.n, np.nan)
        hit = np.array([r in stats for r in w.rsid.tolist()])
        vals = np.array([stats[r] for r in w.rsid[hit].tolist()])
        tb[hit], ts[hit] = vals[:, 0], vals[:, 1]
        has = np.isfinite(tb) & np.isfinite(ts) & (ts > 0) & (w.prow >= 0)
        avail = has & np.isin(w.rsid, array_rs)
        full = measure(w, tb, ts, ld_clump(w, np.flatnonzero(avail), panel), panel, k)
        # bundle rule: clump every candidate once, then take each clump's best available SNP
        cand = np.flatnonzero(has)
        owners = np.full(w.n, -1, dtype=np.int64)
        ld_clump(w, cand, panel, owners=owners)
        own = owners[cand]
        order = np.lexsort((w.p[cand], cand != own, own))
        rows, own = cand[order], own[order]
        start = np.flatnonzero(np.r_[True, own[1:] != own[:-1]])
        rank = np.arange(len(rows)) - np.repeat(start, np.diff(np.r_[start, len(rows)]))
        ok = avail[rows] & (rank < CAP)
        pos = np.where(ok, np.arange(len(rows)), len(rows))
        firsts = np.minimum.reduceat(pos, start)
        sel = rows[firsts[firsts < len(rows)]]
        bundle = measure(w, tb, ts, sel, panel, k)
        print(f"{label:<40}{full[0]:>10,} SNPs {100 * full[1]:>5.2f}%{bundle[0]:>10,} SNPs {100 * bundle[1]:>5.2f}%",
              flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
