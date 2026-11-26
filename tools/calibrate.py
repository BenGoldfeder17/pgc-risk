"""Genome-wide check of the score standardisation and the ancestry check on 1000 Genomes.

For every disorder the production LD clumping runs on all 22 autosomes, once with every SNP
(as for a sequenced or imputed genome) and once with only SNPs found on genotyping arrays,
and the standardised score z is computed for all 2,504 people of 1000 Genomes phase 3.

The European LD panel is those same 503 Europeans, so to keep the test honest the LD
(clumping and the variance factor) is estimated on one half of them and z is evaluated on the
other half: they should show mean ~ 0 and sd ~ 1.  Other ancestries are expected to drift,
which is what the ancestry check is for; its values per population are printed as well.

usage: python tools/calibrate.py [--jobs N] [--chroms 1-22]
"""
from __future__ import annotations

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pgc_risk import panel as pn  # noqa: E402
from pgc_risk.catalog import SOURCES  # noqa: E402
from pgc_risk.score import (DEFAULT_P, DEFAULT_R2, DEFAULT_WINDOW, ld_clump, ld_variance_terms,  # noqa: E402
                            load_weights)

SUPER = ["EUR", "EAS", "SAS", "AMR", "AFR"]
CONFIGS = ["all", "array"]


class HalfPanel:
    """A view of the panel restricted to some of its samples."""

    def __init__(self, panel, cols):
        self.panel, self.cols, self.n = panel, cols, len(cols)

    def dosage(self, c, rows):
        return self.panel.dosage(c, rows)[:, self.cols]


def one_chrom(c: int):
    panel = pn.Panel()
    half_a = HalfPanel(panel, np.arange(0, panel.n, 2))
    cache = np.load(pn._cache(c))
    cpos, geno = cache["pos"], cache["geno"]
    array_rs = np.load(ROOT / "data" / "validation" / "array_like_rsids.npy")
    out = {}
    fit_rows: dict[int, tuple[bool, float, bool]] = {}
    for src in SOURCES:
        if not (pn.WEIGHTS / f"{src.key}.npz").exists():
            continue
        w = load_weights(src.key, pn.WEIGHTS)
        on_c = (w.chrom == c) & (w.prow >= 0)
        on_array = np.isin(w.rsid, array_rs)
        for cfg in CONFIGS:
            cand = np.flatnonzero(on_c & (w.p <= DEFAULT_P) & (on_array if cfg == "array" else True))
            for label, ld in (("half", half_a), ("full", panel)):
                sel = ld_clump(w, cand, ld, DEFAULT_WINDOW, DEFAULT_R2)
                total, indep = ld_variance_terms(w, sel, ld)
                rows = np.searchsorted(cpos, w.pos[sel])
                g = geno[rows].astype(np.float64)
                g = np.where(w.palt[sel][:, None], g, 2 - g)
                b, f = w.beta[sel], w.freq[sel]
                out[(src.key, cfg, label)] = (b @ g, float((2 * f * b).sum()),
                                              float((2 * f * (1 - f) * b * b).sum()), total, indep, len(sel))
        idx = np.flatnonzero(on_c & (w.freq >= 0.05) & (w.freq <= 0.95))
        for i, r in zip(idx.tolist(), np.searchsorted(cpos, w.pos[idx]).tolist()):
            fit_rows.setdefault(r, (bool(w.palt[i]), float(w.freq[i]), bool(on_array[i])))
    rows = np.array(sorted(fit_rows), dtype=np.int64)
    alt = np.array([fit_rows[r][0] for r in rows.tolist()])
    f = np.array([fit_rows[r][1] for r in rows.tolist()])
    arr = np.array([fit_rows[r][2] for r in rows.tolist()])
    fit = np.zeros((2, geno.shape[1]))                         # [all SNPs, array SNPs]
    for s in range(0, len(rows), 20000):                       # in slices, to bound memory
        g = geno[rows[s:s + 20000]].astype(np.float64)
        g = np.where(alt[s:s + 20000, None], g, 2 - g)
        ff = f[s:s + 20000, None]
        dev = (g - 2 * ff) ** 2 / (2 * ff * (1 - ff))
        fit[0] += dev.sum(axis=0)
        fit[1] += dev[arr[s:s + 20000]].sum(axis=0)
    return c, out, fit, (len(rows), int(arr.sum())), cache["samples"].tolist()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--chroms", default="1-22")
    args = ap.parse_args()
    lo, _, hi = args.chroms.partition("-")
    chroms = list(range(int(lo), int(hi or lo) + 1))
    with ProcessPoolExecutor(args.jobs) as pool:
        results = list(pool.map(one_chrom, chroms))
    names = results[0][4]
    info = pn.samples()
    sup = np.array([info[s][1] for s in names])
    pop = np.array([info[s][0] for s in names])
    eur_idx = np.flatnonzero(sup == "EUR")
    half_b = eur_idx[1::2]                 # the Europeans NOT used for LD in the "half" runs
    keys = [s.key for s in SOURCES if (pn.WEIGHTS / f"{s.key}.npz").exists()]

    def z_of(key, cfg, label):
        parts = [r[1][(key, cfg, label)] for r in results]
        raw = sum(p[0] for p in parts)
        mean, var = sum(p[1] for p in parts), sum(p[2] for p in parts)
        infl = sum(p[3] for p in parts) / sum(p[4] for p in parts)
        return (raw - mean) / np.sqrt(var * infl), infl, sum(p[5] for p in parts)

    print(f"window {DEFAULT_WINDOW:,}, r2 {DEFAULT_R2}, p <= {DEFAULT_P}, chromosomes {args.chroms}")
    zs = {}
    for cfg in CONFIGS:
        print(f"\n[{cfg} SNPs]  LD from half the Europeans, tested on the other half ({len(half_b)} people)")
        print(f"{'':<6}{'SNPs':>8}{'LD var':>8}{'mean z':>9}{'sd z':>7}   | full panel: "
              + " ".join(f"{s:>12}" for s in SUPER))
        sds, means = [], []
        for key in keys:
            zh, infl, n_snps = z_of(key, cfg, "half")
            zf, infl_f, _ = z_of(key, cfg, "full")
            zs[(key, cfg)] = zf
            sds.append(zh[half_b].std(ddof=1))
            means.append(zh[half_b].mean())
            print(f"{key:<6}{n_snps:>8,}{infl:>8.2f}{means[-1]:>+9.2f}{sds[-1]:>7.2f}   |             "
                  + " ".join(f"{zf[sup == s].mean():>+6.2f}/{zf[sup == s].std():>5.2f}" for s in SUPER))
        print(f"held-out Europeans: mean z {np.mean(means):+.3f} (range {min(means):+.2f}..{max(means):+.2f}), "
              f"sd z {np.mean(sds):.3f} (range {min(sds):.2f}..{max(sds):.2f}); "
              f"sampling error of each sd is about {1 / np.sqrt(2 * len(half_b)):.2f}")
    print("\nEuropean sub-population mean z (all SNPs, full panel)")
    eur_pops = sorted(set(pop[sup == "EUR"]))
    print(f"{'':<6}" + "".join(f"{p:>8}" for p in eur_pops))
    for key in keys:
        z = zs[(key, "all")]
        print(f"{key:<6}" + "".join(f"{z[pop == p].mean():>+8.2f}" for p in eur_pops))
    fits = {}
    for j, cfg in enumerate(CONFIGS):
        n_fit = sum(r[3][j] for r in results)
        fit = sum(r[2][j] for r in results) / n_fit
        fits[cfg] = fit
        e = fit[sup == "EUR"]
        print(f"\nancestry fit, {cfg} SNPs ({n_fit:,}); 1.00 = matches the European reference")
        print(f"Europeans: median {np.median(e):.3f}, 95% below {np.percentile(e, 95):.3f}, "
              f"99% below {np.percentile(e, 99):.3f}, max {e.max():.3f}")
        print(f"{'group':<10}{'n':>5}{'min':>8}{'median':>8}{'max':>8}")
        for s_ in SUPER:
            m = sup == s_
            print(f"{s_:<10}{m.sum():>5}{fit[m].min():>8.3f}{np.median(fit[m]):>8.3f}{fit[m].max():>8.3f}")
        for p_ in sorted(set(pop), key=lambda q: (sup[pop == q][0], q)):
            m = pop == p_
            print(f"  {sup[m][0]:<4}{p_:<4}{m.sum():>5}{fit[m].min():>8.3f}{np.median(fit[m]):>8.3f}{fit[m].max():>8.3f}")
    fit = fits["all"]
    np.savez(ROOT / "data" / "validation" / "calibration_scores.npz", samples=np.array(names), fit=fit,
             fit_array=fits["array"],
             **{f"{k}_{cfg}": zs[(k, cfg)] for k in keys for cfg in CONFIGS})


if __name__ == "__main__":
    main()
