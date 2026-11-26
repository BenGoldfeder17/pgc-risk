"""Out-of-sample accuracy of the scorer, measured from summary statistics alone.

PGC publishes some GWAS twice: with and without one cohort (UK Biobank), or split by sex.
That gives a training GWAS and an independent target:

* with/without UK Biobank: the UK Biobank-only association is recovered by undoing the
  fixed-effects meta-analysis,  w_ukb = w_full - w_train,
  b_ukb = (w_full*b_full - w_train*b_train) / w_ukb,   w = 1/se^2.  This needs the held-out
  cohort to carry a few percent of the total weight, otherwise the rounding of the published
  standard errors swamps the difference (which rules it out for bipolar disorder, where
  UK Biobank is about 2% of the meta-analysis);
* male/female: train on one sex, test on the other.

A score built from the training GWAS with the production matching + clumping code is then
regressed on the target summary statistics (the `gtx` approach for independent SNPs):
    a = sum(b_i t_i / s_i^2) / sum(b_i^2 / s_i^2),   se(a) = 1 / sqrt(sum(b_i^2 / s_i^2))
with b = training log-odds, t and s = target log-odds and standard error.  The variance the
score explains in a balanced case-control sample is a^2 var(score) / 4, which the formula of
Lee et al. (2012) moves to the liability scale.

usage: python tools/loo.py fetch | build | run [MDD|OCD|SCZ]
(after build, run: pgc-risk panel --also data/validation/loo)
"""
from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path
from statistics import NormalDist

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pgc_risk import build  # noqa: E402
from pgc_risk.catalog import BY_KEY  # noqa: E402
from pgc_risk.fetch import fetch_sources, local_files  # noqa: E402
from pgc_risk.panel import Panel  # noqa: E402
from pgc_risk.score import ld_clump, ld_inflation, load_weights  # noqa: E402

OUT = ROOT / "data" / "validation" / "loo"
ND = NormalDist()

MDD_NOUKB = replace(BY_KEY["mdd"], key="mdd_noukb", series=714,
                    source_file="pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11_formatted.tsv.gz")
OCD_NOUKB = replace(BY_KEY["ocd"], key="ocd_noukb", series=676, study_freq="FRQ_U_988884",
                    source_file="daner_OCDmeta_wo23andMe_LOOUKBB_080425.gz")
SCZ_MALE = replace(BY_KEY["scz"], key="scz_male", series=759, study_freq="FRQ_U_41159",
                   source_file="daner_PGC_SCZ_w3_81_0618a_all_male.gz")
SCZ_FEMALE = replace(BY_KEY["scz"], key="scz_female", series=755, study_freq="FRQ_U_45626",
                     source_file="daner_PGC_SCZ_w3_81_0618a_all_female.gz")
EXTRA = [MDD_NOUKB, OCD_NOUKB, SCZ_MALE, SCZ_FEMALE]

# (label, training source, target: ("minus", full source) or ("direct", source), prevalence K)
DESIGNS = [
    ("MDD: all but UK Biobank -> UK Biobank", MDD_NOUKB, ("minus", BY_KEY["mdd"]), 0.15),
    ("OCD: all but UK Biobank -> UK Biobank", OCD_NOUKB, ("minus", BY_KEY["ocd"]), 0.01),
    ("SCZ: men -> women", SCZ_MALE, ("direct", SCZ_FEMALE), 0.01),
    ("SCZ: women -> men", SCZ_FEMALE, ("direct", SCZ_MALE), 0.01),
]


def liability_r2(r2_obs: float, k: float, p: float = 0.5) -> float:
    """Lee et al. 2012: observed-scale R2 in a case-control sample -> liability scale."""
    t = ND.inv_cdf(1 - k)
    z = ND.pdf(t)
    m = z / k
    c = (k * (1 - k)) ** 2 / (z * z * p * (1 - p))
    theta = m * (p - k) / (1 - k) * (m * (p - k) / (1 - k) - t)
    return r2_obs * c / (1 + r2_obs * theta * c)


def target_stats(con, train, target) -> dict[int, tuple[float, float]]:
    """rsid -> (log-odds for the training file's a1, standard error) in the independent target."""
    kind, other = target
    build.load_raw(con, train, "t0")
    build.load_raw(con, other, "f0")
    rows = con.execute("""
        SELECT t.rsid, t.beta AS bt, t.se AS st, f.se AS sf,
               CASE WHEN f.a1 = t.a1 AND f.a2 = t.a2 THEN f.beta
                    WHEN f.a1 = t.a2 AND f.a2 = t.a1 THEN -f.beta END AS bf
        FROM t0 t JOIN f0 f USING (rsid) WHERE t.se > 0 AND f.se > 0""").fetchnumpy()
    con.execute("DROP TABLE t0; DROP TABLE f0")
    rs = np.asarray(rows["rsid"])
    bt, st = np.ma.filled(rows["bt"].astype(float), np.nan), np.ma.filled(rows["st"].astype(float), np.nan)
    bf, sf = np.ma.filled(rows["bf"].astype(float), np.nan), np.ma.filled(rows["sf"].astype(float), np.nan)
    ok = np.isfinite(bf)
    if kind == "direct":
        b, s = bf, sf
    else:
        wf, wt = 1 / sf ** 2, 1 / st ** 2
        wd = wf - wt
        ok &= wd > 0.02 * wf                 # the held-out cohort must actually contribute
        with np.errstate(divide="ignore", invalid="ignore"):
            b = (wf * bf - wt * bt) / wd
            s = np.sqrt(1 / wd)
    return dict(zip(rs[ok].tolist(), zip(b[ok].tolist(), s[ok].tolist())))


def evaluate(w, tb, ts, avail, panel, window, r2, p_thr, k):
    cand = np.flatnonzero(avail & (w.p <= p_thr) & (w.prow >= 0))
    sel = ld_clump(w, cand, panel, window, r2)
    if len(sel) < 5:
        return None
    infl = ld_inflation(w, sel, panel)
    b, t, s, f = w.beta[sel], tb[sel], ts[sel], w.freq[sel]
    den = (b * b / s ** 2).sum()
    a = (b * t / s ** 2).sum() / den
    # cov(score, y) = sum b_i t_i var(g_i) holds whatever the LD; var(score) needs the LD factor
    r2_obs = a * a * (2 * f * (1 - f) * b * b).sum() / 4 / infl
    neff = float(np.median(4 / (s ** 2 * 2 * f * (1 - f))))
    return {"snps": len(sel), "slope": a, "z": a * den ** 0.5 / infl ** 0.5, "r2_obs": r2_obs,
            "inflation": infl, "r2_liab": liability_r2(r2_obs, k) * (1 if a > 0 else -1),
            "target_neff": neff}


def build_training(only=()) -> None:
    con = build._connect()
    build._ensure_freq(con, lambda *a: None)
    for src in EXTRA:
        if only and src.key not in only:
            continue
        if local_files(src):
            build.build_source(con, src, log=print, out_dir=OUT)


def run(only=(), windows=(250_000, 500_000), r2s=(0.1, 0.2), p_thrs=(5e-8, 1e-4, 1e-2, 0.05)) -> None:
    con = build._connect()
    build._ensure_freq(con, lambda *a: None)
    panel = Panel()
    array_rs = np.load(ROOT / "data" / "validation" / "array_like_rsids.npy")
    for label, train, target, k in DESIGNS:
        if only and not any(label.upper().startswith(o.upper()) for o in only):
            continue
        if not (OUT / f"{train.key}.npz").exists():
            print(f"\n{label}: training weights not built")
            continue
        w = load_weights(train.key, OUT)
        stats = target_stats(con, train, target)
        tb = np.full(w.n, np.nan)
        ts = np.full(w.n, np.nan)
        hit = np.array([r in stats for r in w.rsid.tolist()])
        vals = np.array([stats[r] for r in w.rsid[hit].tolist()])
        tb[hit], ts[hit] = vals[:, 0], vals[:, 1]
        has = np.isfinite(tb) & np.isfinite(ts) & (ts > 0)
        on_array = np.isin(w.rsid, array_rs)
        print(f"\n{label}   (K = {k:g}; {int(has.sum()):,} of {w.n:,} training SNPs have a target "
              f"estimate, {int((w.prow >= 0).sum()):,} are in the LD panel)")
        print(f"  {'SNP set':<12}{'window':>8}{'r2':>5}{'p <=':>8}{'SNPs':>9}{'slope':>7}{'LD var':>8}"
              f"{'z':>7}{'R2 liab':>9}{'target Neff':>13}")
        for name, avail in (("all", has), ("array-like", has & on_array)):
            for window in windows:
                for r2 in r2s:
                    for p_thr in p_thrs:
                        r = evaluate(w, tb, ts, avail, panel, window, r2, p_thr, k)
                        if r is None:
                            continue
                        print(f"  {name:<12}{window // 1000:>7}k{r2:>5}{p_thr:>8.0e}{r['snps']:>9,}"
                              f"{r['slope']:>7.3f}{r['inflation']:>8.2f}{r['z']:>7.1f}"
                              f"{100 * r['r2_liab']:>8.2f}%{r['target_neff']:>13,.0f}", flush=True)


if __name__ == "__main__":
    if sys.argv[1:2] == ["fetch"]:
        fetch_sources(EXTRA)
    elif sys.argv[1:2] == ["build"]:
        build_training(sys.argv[2:])
    else:
        run(only=[a for a in sys.argv[1:] if a != "run"])
