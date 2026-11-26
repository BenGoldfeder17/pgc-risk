"""Does scoring from a bundle agree with the full method?  Checked on 1000 Genomes.

For all 2,504 people the bundle's scoring rule (best available SNP per pre-computed clump,
per-SNP LD factor) is applied to the genotype matrix, once with every SNP available and once
with only array SNPs, and compared with what tools/calibrate.py computed with the full method
(LD clumping on exactly the SNPs available, LD factor measured on the chosen SNPs).

Reported per disorder: correlation and mean absolute difference of z among Europeans, the
spread of z among Europeans, and the bundle's LD factor next to the one measured directly on
the SNPs it selected.  Also the ancestry check's values per population.

usage: python tools/check_bundle.py [full|chips] [--jobs 8]      (run tools/calibrate.py first)
"""
from __future__ import annotations

import argparse
import glob
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pgc_risk import panel as pn  # noqa: E402
from pgc_risk.bundle import read_bundle  # noqa: E402
from pgc_risk.score import LD_BAND  # noqa: E402

BUNDLES = ROOT / "data" / "bundle"
CONFIGS = ["all", "chips", "array", "random30", "random10"]      # which SNPs the genome file has
REFERENCE = ("all", "array")                                    # configurations tools/calibrate.py also ran
SUPER = ["EUR", "EAS", "SAS", "AMR", "AFR"]


def snp_positions(bundle) -> tuple[np.ndarray, np.ndarray]:
    """Chromosome and position of every bundle SNP, looked up in the weight files."""
    rs, ch, po = [], [], []
    for f in sorted(glob.glob(str(pn.WEIGHTS / "*.npz"))):
        z = np.load(f)
        rs.append(z["rsid"]); ch.append(z["chrom"]); po.append(z["pos"])
    rs, ch, po = np.concatenate(rs), np.concatenate(ch), np.concatenate(po)
    u, first = np.unique(rs, return_index=True)
    j = np.searchsorted(u, bundle.rsid.astype(np.int64))
    assert np.array_equal(u[j], bundle.rsid.astype(np.int64))
    return ch[first][j], po[first][j]


def one_chrom(args):
    c, variant = args
    bundle = read_bundle(BUNDLES / f"pgc-risk-reference-{variant}-v1.pgcr")
    chrom, pos = snp_positions(bundle)
    snps = np.flatnonzero(chrom == c)
    cache = np.load(pn._cache(c))
    names = cache["samples"].tolist()
    info = pn.samples()
    eur = np.array([info[s][1] == "EUR" for s in names])
    j = np.searchsorted(cache["pos"], pos[snps])
    j[j == len(cache["pos"])] = 0
    ref, alt = cache["ref"][j].astype(np.int16), cache["alt"][j].astype(np.int16)
    other = (ref != 0) & (alt != 0)                           # no A on this strand: read the other one
    ref, alt = np.where(other, 3 - ref, ref), np.where(other, 3 - alt, alt)
    b = bundle.b_allele[snps].astype(np.int16)
    a_is_alt = (alt == 0) & (ref == b)
    a_is_ref = (ref == 0) & (alt == b)
    present = (cache["pos"][j] == pos[snps]) & (a_is_alt | a_is_ref)
    counts = cache["geno"][j].astype(np.int8)
    counts = np.where(a_is_alt[:, None], counts, 2 - counts)   # A-allele counts, SNPs x people
    array_rs = np.load(ROOT / "data" / "validation" / "array_like_rsids.npy")
    chips_rs = np.load(ROOT / "data" / "validation" / "chip_rsids.npy")
    rs = bundle.rsid[snps].astype(np.int64)
    pick = ((rs * 2654435761) % 2 ** 32) / 2 ** 32                # a fixed pseudo-random number per SNP
    avail = {"all": present, "chips": present & np.isin(rs, chips_rs), "array": present & np.isin(rs, array_rs),
             "random30": present & (pick < 0.30), "random10": present & (pick < 0.10)}
    local = np.full(len(bundle.rsid), -1, dtype=np.int64)
    local[snps] = np.arange(len(snps))
    freq = bundle.freq[snps]
    out = {"scores": {}, "fit": {}}
    for cfg in CONFIGS:
        m = np.flatnonzero(avail[cfg] & (freq >= 0.05) & (freq <= 0.95))
        total = np.zeros(counts.shape[1])
        for s0 in range(0, len(m), 10000):                     # in slices, to bound memory
            f = freq[m[s0:s0 + 10000]][:, None]
            total += (((counts[m[s0:s0 + 10000]] - 2 * f) ** 2) / (2 * f * (1 - f))).sum(axis=0)
        out["fit"][cfg] = (total, len(m))
    big = np.iinfo(np.int64).max
    for d in bundle.disorders:
        loc = local[d.idx]
        on_c = np.flatnonzero(loc >= 0)                        # this chromosome's rows are contiguous
        if not len(on_c):
            continue
        lo, hi = on_c[0], on_c[-1] + 1
        starts = d.clump_start[(d.clump_start >= lo) & (d.clump_start < hi)] - lo
        loc, beta, ld = loc[lo:hi], d.beta[lo:hi], d.ld[lo:hi]
        for cfg in CONFIGS:
            rows = np.where(avail[cfg][loc], np.arange(hi - lo), big)
            first = np.minimum.reduceat(rows, starts)
            rep = first[first < big]
            s = loc[rep]
            bb, ff = beta[rep], freq[s]
            v = 2 * ff * (1 - ff) * bb * bb
            # the LD factor measured directly on the chosen SNPs (pairs within the band), on the panel
            order = np.argsort(pos[snps][s], kind="stable")
            y = (counts[s][order][:, eur] - counts[s][order][:, eur].mean(axis=1, keepdims=True)) * bb[order][:, None]
            p_sorted = pos[snps][s][order].astype(np.int64)
            csum = np.vstack([np.zeros((1, y.shape[1])), np.cumsum(y, axis=0)])
            w_lo = np.searchsorted(p_sorted, p_sorted - LD_BAND)
            w_hi = np.searchsorted(p_sorted, p_sorted + LD_BAND, side="right")
            out["scores"][(d.meta["key"], cfg)] = (
                bb @ counts[s].astype(np.float64), float((2 * ff * bb).sum()), float((v * (ld[rep] - 1)).sum()),
                float(v.sum()), float((y * (csum[w_hi] - csum[w_lo])).sum()), float((y * y).sum()), len(rep))
    return c, out, names


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("variant", nargs="?", default="full")
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args()
    with ProcessPoolExecutor(args.jobs) as pool:
        results = list(pool.map(one_chrom, [(c, args.variant) for c in range(1, 23)]))
    names = results[0][2]
    info = pn.samples()
    sup = np.array([info[s][1] for s in names])
    pop = np.array([info[s][0] for s in names])
    eur = sup == "EUR"
    ref = np.load(ROOT / "data" / "validation" / "calibration_scores.npz")
    assert ref["samples"].tolist() == names
    bundle = read_bundle(BUNDLES / f"pgc-risk-reference-{args.variant}-v1.pgcr")
    print(f"bundle: {args.variant} ({bundle.header['snps']:,} SNPs)")
    from pgc_risk.bundle import LD_SHARE_POWER
    points = []
    for cfg in CONFIGS:
        print(f"\n[{cfg}]  Europeans; bundle scoring, and where available the full method on the same SNPs")
        print(f"{'':<6}{'SNPs':>8}{'share':>7}{'LD factor':>11}{'measured':>10}{'mean z':>8}{'sd z':>7}"
              f"{'corr':>8}{'mean |dz|':>11}")
        sds, errs = [], []
        for d in bundle.disorders:
            k = d.meta["key"]
            parts = [r[1]["scores"][(k, cfg)] for r in results if (k, cfg) in r[1]["scores"]]
            raw = sum(p[0] for p in parts)
            mean, cross, indep = (sum(p[i] for p in parts) for i in (1, 2, 3))
            measured = sum(p[4] for p in parts) / sum(p[5] for p in parts)
            share = min(1.0, indep / d.meta["full_indep"])
            factor = 1 + share ** LD_SHARE_POWER * cross / indep
            points.append((share, cross / indep, measured))
            z = (raw - mean) / np.sqrt(indep * factor)
            sds.append(z[eur].std(ddof=1)); errs.append(factor / measured - 1)
            extra = ""
            if cfg in REFERENCE:
                z_ref = ref[f"{k}_{cfg}"]
                extra = (f"{np.corrcoef(z[eur], z_ref[eur])[0, 1]:>8.3f}"
                         f"{np.abs(z[eur] - z_ref[eur]).mean():>11.3f}")
            print(f"{k:<6}{sum(p[6] for p in parts):>8,}{share:>7.2f}{factor:>11.2f}{measured:>10.2f}"
                  f"{z[eur].mean():>+8.2f}{z[eur].std(ddof=1):>7.2f}{extra}")
        print(f"LD factor vs measured: {100 * np.mean(errs):+.1f}% on average ({100 * min(errs):+.1f}% .. "
              f"{100 * max(errs):+.1f}%); sd of z {np.mean(sds):.3f} ({min(sds):.2f}..{max(sds):.2f})")
        n_fit = sum(r[1]["fit"][cfg][1] for r in results)
        fit = sum(r[1]["fit"][cfg][0] for r in results) / n_fit
        e = fit[eur]
        print(f"ancestry fit over {n_fit:,} SNPs - Europeans: median {np.median(e):.3f}, 95% below "
              f"{np.percentile(e, 95):.3f}, 99% below {np.percentile(e, 99):.3f}, max {e.max():.3f}; "
              + "; ".join(f"{s_} min {fit[sup == s_].min():.3f}" for s_ in SUPER[1:]))
    # which power of the represented share best turns the stored cross terms into the measured ones
    sh, cr, me = (np.array(x) for x in zip(*points))
    part = sh < 0.999
    grid = np.linspace(0.3, 1.5, 121)
    loss = [np.mean((np.log(1 + sh[part] ** a * cr[part]) - np.log(me[part])) ** 2) for a in grid]
    best = grid[int(np.argmin(loss))]
    now = np.sqrt(np.mean((np.log(1 + sh[part] ** LD_SHARE_POWER * cr[part]) - np.log(me[part])) ** 2))
    print(f"\nbest power of the represented share: {best:.2f} (rms error {100 * np.sqrt(min(loss)):.1f}%); "
          f"in use: {LD_SHARE_POWER} (rms error {100 * now:.1f}%)")


if __name__ == "__main__":
    main()
