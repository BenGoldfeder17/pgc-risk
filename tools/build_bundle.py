"""Build the portable reference bundle (see pgc_risk/bundle.py) from the weights and the LD panel.

Two variants:
  full   every candidate SNP of every disorder (works with any genotyping chip; larger)
  chips  only SNPs that sit on genotyping chips: genotyped directly, rather than imputed, in
         at least three cohorts of the schizophrenia or the bipolar GWAS (smaller; meant for
         phones and slow connections)

For each disorder the candidates are LD-clumped once with the production code; every SNP is
filed under its clump's index SNP, most significant first, and given its LD factor.  Within a
clump only the `--cap` most significant SNPs are kept: a genome that has none of those would
have to rely on a much weaker stand-in anyway.

usage: python tools/build_bundle.py [full] [chips] [--cap 48] [--jobs 6]
"""
from __future__ import annotations

import argparse
import datetime
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from statistics import NormalDist

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pgc_risk import __version__  # noqa: E402
from pgc_risk import panel as pn  # noqa: E402
from pgc_risk.bundle import FREQ_SCALE, LD_SHARE_POWER, LD_STEP, write_bundle  # noqa: E402
from pgc_risk.catalog import BY_KEY, LOW_POWER_R2, RISK, SOURCES  # noqa: E402
from pgc_risk.fetch import local_files  # noqa: E402
from pgc_risk.score import LD_BAND, ld_clump, load_weights  # noqa: E402

OUT = ROOT / "data" / "bundle"
FIT_GOOD, FIT_POOR, MIN_FIT_SNPS = 1.04, 1.075, 2000


def chip_rsids(min_cohorts: int = 3) -> np.ndarray:
    """rs numbers genotyped directly in at least `min_cohorts` cohorts of the SCZ or BIP GWAS."""
    import duckdb

    con = duckdb.connect()
    found = []
    for key in ("scz", "bip"):
        flist = "[" + ",".join("'" + str(f).replace(chr(92), "/") + "'" for f in local_files(BY_KEY[key])) + "]"
        found.append(np.asarray(con.execute(f"""
            SELECT DISTINCT rs FROM (
                SELECT TRY_CAST(regexp_extract(SNP, '^rs([0-9]+)', 1) AS BIGINT) AS rs
                FROM read_parquet({flist}, union_by_name=true)
                WHERE TRY_CAST(ngt AS INTEGER) >= {min_cohorts} AND CHR BETWEEN 1 AND 22)
            WHERE rs IS NOT NULL""").fetchnumpy()["rs"], dtype=np.int64))
    return np.unique(np.concatenate(found))


def effect_dosage(w, panel, c: int, rows: np.ndarray) -> np.ndarray:
    """Centred effect-allele counts of the panel, times the SNP's effect: one row per SNP."""
    g = panel.dosage(c, w.prow[rows]).astype(np.float64)
    g = np.where(w.palt[rows][:, None], g, 2 - g)
    return (g - g.mean(axis=1, keepdims=True)) * w.beta[rows][:, None]


def clump_disorder(args) -> dict:
    """LD-clump one disorder's candidates and lay them out clump by clump, with LD factors."""
    key, universe, cap = args
    t0 = time.time()
    w = load_weights(key, pn.WEIGHTS)
    panel = pn.Panel()
    mask = w.prow >= 0
    if universe is not None:
        mask &= np.isin(w.rsid, universe)
    cand = np.flatnonzero(mask)
    owners = np.full(w.n, -1, dtype=np.int64)
    index = np.sort(ld_clump(w, cand, panel, owners=owners))
    own = owners[cand]
    order = np.lexsort((w.p[cand], cand != own, own))            # by clump; index SNP first; then by p
    rows, own = cand[order], own[order]
    start = np.flatnonzero(np.r_[True, own[1:] != own[:-1]])
    rank = np.arange(len(rows)) - np.repeat(start, np.diff(np.r_[start, len(rows)]))
    sizes_all = np.diff(np.r_[start, len(rows)])
    keep = rank < cap
    rows, own = rows[keep], own[keep]
    start = np.flatnonzero(np.r_[True, own[1:] != own[:-1]])
    clump_len = np.diff(np.r_[start, len(rows)])
    assert len(clump_len) == len(index) and np.array_equal(rows[start], index)

    # LD factor: 1 + cov(this SNP, the index SNPs of the other clumps within the band) / var(this SNP)
    ld = np.ones(len(rows))
    for c in np.unique(w.chrom[rows]).tolist():
        ii = index[w.chrom[index] == c]
        yi = effect_dosage(w, panel, c, ii)
        csum = np.vstack([np.zeros((1, yi.shape[1])), np.cumsum(yi, axis=0)])
        ipos = w.pos[ii].astype(np.int64)
        on_c = np.flatnonzero(w.chrom[rows] == c)
        for s in range(0, len(on_c), 4000):
            k = on_c[s:s + 4000]
            ym = effect_dosage(w, panel, c, rows[k])
            pos = w.pos[rows[k]].astype(np.int64)
            lo = np.searchsorted(ipos, pos - LD_BAND)
            hi = np.searchsorted(ipos, pos + LD_BAND, side="right")
            others = csum[hi] - csum[lo]
            mine = np.searchsorted(ii, own[k])
            inside = (mine >= lo) & (mine < hi)
            others -= np.where(inside[:, None], yi[mine], 0.0)   # its own clump's index SNP does not count
            self_var = np.einsum("ij,ij->i", ym, ym)
            cross = np.einsum("ij,ij->i", ym, others)
            ld[k] = 1 + np.divide(cross, self_var, out=np.zeros_like(cross), where=self_var > 0)
    return {"key": key, "rsid": w.rsid[rows], "chrom": w.chrom[rows], "pos": w.pos[rows], "a1": w.a1[rows],
            "a2": w.a2[rows], "beta": w.beta[rows], "freq": w.freq[rows], "ld": ld, "clump_len": clump_len,
            "candidates": int(len(cand)), "largest_clump": int(sizes_all.max()),
            "dropped_by_cap": int(len(cand) - len(rows)), "seconds": round(time.time() - t0)}


def assemble(name: str, parts: list[dict], description: str, cap: int) -> tuple[dict, dict]:
    """Merge the disorders onto one strand-free SNP table and quantise; returns (header, arrays)."""
    canon = []
    for p in parts:
        a1, a2 = p["a1"].astype(np.int16), p["a2"].astype(np.int16)
        other_strand = (a1 != 0) & (a2 != 0)                    # neither allele is A: read the other strand
        a1 = np.where(other_strand, 3 - a1, a1)
        a2 = np.where(other_strand, 3 - a2, a2)
        effect_is_a = a1 == 0
        b = np.where(effect_is_a, a2, a1)
        assert ((a1 == 0) ^ (a2 == 0)).all() and np.isin(b, (1, 2)).all()
        canon.append((b.astype(np.uint8), np.where(effect_is_a, p["freq"], 1 - p["freq"]),
                      np.where(effect_is_a, p["beta"], -p["beta"])))
    rs = np.concatenate([p["rsid"] for p in parts])
    key = (np.concatenate([p["chrom"] for p in parts]).astype(np.int64) << 32) | np.concatenate([p["pos"] for p in parts])
    b_all = np.concatenate([c[0] for c in canon])
    f_all = np.concatenate([c[1] for c in canon])
    uniq, first = np.unique(rs, return_index=True)
    order = np.argsort(key[first], kind="stable")               # SNP table in genome order
    snp_rs = uniq[order]
    snp_b, snp_f = b_all[first][order], f_all[first][order]
    where = np.argsort(snp_rs, kind="stable")
    lookup = lambda r: where[np.searchsorted(snp_rs[where], r)]
    pos_of = lookup(rs)
    assert np.array_equal(snp_b[pos_of], b_all), "a SNP has different alleles in two GWAS"
    assert np.abs(snp_f[pos_of] - f_all).max() < 1e-4, "a SNP has different frequencies in two GWAS"
    freq_q = np.clip(np.rint(snp_f * FREQ_SCALE), 1, FREQ_SCALE - 1).astype(np.uint16)
    freq = freq_q / FREQ_SCALE
    arrays = {"snp.rsid": snp_rs.astype(np.uint32), "snp.b": snp_b, "snp.freq": freq_q}
    disorders = []
    for p, (_, _, beta_a) in zip(parts, canon):
        src, risk = BY_KEY[p["key"]], RISK[p["key"]]
        idx = lookup(p["rsid"])
        scale = float(np.abs(beta_a).max() / 32767)
        beta_q = np.rint(beta_a / scale).astype(np.int16)
        ld_q = np.clip(np.rint(p["ld"] / LD_STEP), 0, 255).astype(np.uint8)
        k = p["key"]
        arrays[f"{k}.idx_delta"] = np.diff(idx, prepend=0).astype(np.int32)
        arrays[f"{k}.beta"] = beta_q
        arrays[f"{k}.ld"] = ld_q
        arrays[f"{k}.clump_len"] = p["clump_len"].astype(np.uint32)
        start = np.concatenate([[0], np.cumsum(p["clump_len"])[:-1]])
        b, f = beta_q[start] * scale, freq[idx[start]]          # as a reader will see the index SNPs
        full_indep = float((2 * f * (1 - f) * b * b).sum())    # variance if every clump has its index SNP
        disorders.append({
            "key": k, "name": src.name, "prevalence": risk.prevalence, "prevalence_source": risk.prevalence_source,
            "liability_threshold": NormalDist().inv_cdf(1 - risk.prevalence),
            "r2": risk.r2_liability, "r2_basis": risk.extra.get("basis", ""), "r2_source": risk.extra.get("source", ""),
            "citation": src.citation, "ref": src.ref, "notes": src.notes,
            "rows": int(len(idx)), "clumps": int(len(start)), "beta_scale": scale, "full_indep": full_indep})
    header = {
        "name": name, "version": 1, "tool": f"pgc-risk {__version__}",
        "created": datetime.date.today().isoformat(), "description": description,
        "snps": int(len(snp_rs)), "ld_step": LD_STEP, "freq_scale": FREQ_SCALE, "ld_share_power": LD_SHARE_POWER,
        "fit_good": FIT_GOOD, "fit_poor": FIT_POOR, "min_fit_snps": MIN_FIT_SNPS, "low_power_r2": LOW_POWER_R2,
        "method": {"p_threshold": 0.05, "clump_r2": 0.1, "clump_window_bp": 250000, "ld_band_bp": LD_BAND,
                   "max_snps_per_clump": cap, "ld_panel": pn.Panel().meta["source"]},
        "attribution": "GWAS summary statistics: Psychiatric Genomics Consortium (https://pgc.unc.edu/), CC BY 4.0; "
                       "cite the paper listed for each disorder. LD: 1000 Genomes Project phase 3.",
        "disorders": disorders}
    return header, arrays


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("variants", nargs="*", default=["full", "chips"])
    ap.add_argument("--cap", type=int, default=48)
    ap.add_argument("--jobs", type=int, default=6)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    keys = [s.key for s in SOURCES if (pn.WEIGHTS / f"{s.key}.npz").exists()]
    for variant in args.variants:
        if variant == "chips":
            universe = chip_rsids()
            np.save(ROOT / "data" / "validation" / "chip_rsids.npy", universe)
            description = ("SNPs that sit on genotyping chips (genotyped directly in at least three cohorts of "
                           "the PGC schizophrenia or bipolar GWAS)")
        else:
            universe = None
            description = "every candidate SNP (works with any genotyping chip)"
        t0 = time.time()
        with ProcessPoolExecutor(args.jobs) as pool:
            parts = list(pool.map(clump_disorder, [(k, universe, args.cap) for k in keys]))
        header, arrays = assemble(f"pgc-risk reference data, {variant}", parts, description, args.cap)
        path = OUT / f"pgc-risk-reference-{variant}-v1.pgcr"
        size = write_bundle(path, header, arrays)
        raw = sum(a.nbytes for a in arrays.values())
        print(f"\n[{variant}] {header['snps']:,} SNPs, {sum(d['rows'] for d in header['disorders']):,} rows, "
              f"{raw / 1e6:.1f} MB raw -> {size / 1e6:.1f} MB at {path}  ({time.time() - t0:.0f}s)")
        print(f"{'':<6}{'candidates':>11}{'rows':>9}{'clumps':>8}{'largest':>9}{'cut by cap':>11}")
        for p, d in zip(parts, header["disorders"]):
            print(f"{p['key']:<6}{p['candidates']:>11,}{d['rows']:>9,}{d['clumps']:>8,}{p['largest_clump']:>9,}"
                  f"{p['dropped_by_cap']:>11,}")


if __name__ == "__main__":
    main()
