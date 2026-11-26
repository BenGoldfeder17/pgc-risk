"""Two independent checks that every GWAS's effect allele is read the right way round.

1. Within a study that reports case and control allele frequencies, a risk-increasing allele
   must be more common in cases: sign(effect) should equal sign(freq_cases - freq_controls).
   Only strong effects are used: in a large meta-analysis the pooled case and control
   frequencies also differ because cohorts contribute cases and controls in different
   proportions, which swamps the tiny shift a weak association causes.
2. Across studies, psychiatric disorders are positively genetically correlated, so on a panel
   of evenly spaced SNPs the association Z-scores of two disorders (aligned to the same
   allele) correlate positively.  A study read with flipped alleles shows up as a row of
   negative correlations.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pgc_risk import build  # noqa: E402
from pgc_risk.catalog import SOURCES  # noqa: E402
from pgc_risk.fetch import local_files  # noqa: E402


def case_control_check(con, src) -> str:
    files = local_files(src)
    flist = "[" + ",".join("'" + str(f).replace("\\", "/") + "'" for f in files) + "]"
    scan = f"read_parquet({flist}, union_by_name=true)"
    cols = [r[0] for r in con.execute(f"DESCRIBE SELECT * FROM {scan}").fetchall()]
    fa = next((c for c in cols if re.match(r"FRQ_A_\d+$", c)), None)
    fu = next((c for c in cols if re.match(r"FRQ_U_\d+$", c)), None)
    if src.key == "mdd":
        fa, fu = "effect_allele_frequency_cases", "effect_allele_frequency"
    if not fa or not fu:
        return "no case/control frequencies"
    eff = f'ln("{src.effect}")' if src.effect_kind == "or" else f'"{src.effect}"'
    n, agree = con.execute(f"""
        SELECT count(*), avg((sign({eff}) = sign("{fa}" - "{fu}"))::INT)
        FROM {scan} WHERE TRY_CAST("{src.p}" AS DOUBLE) < 1e-6 AND "{fa}" <> "{fu}"
              AND abs({eff}) >= 0.05""").fetchone()
    if not n or n < 20:
        return "too few strong associations to tell"
    return f"{agree:.3f} of {n:,} SNPs (p<1e-6, odds ratio beyond 1.05) have the effect allele enriched in cases"


def main() -> None:
    con = build._connect()
    build._ensure_freq(con, print)
    # panel: one common SNP per 100 kb
    con.execute("""CREATE TABLE panel AS SELECT rsid, a1, a2 FROM eurfreq
                   WHERE least(freq, 1 - freq) >= 0.1
                   QUALIFY row_number() OVER (PARTITION BY chr, pos // 100000 ORDER BY hash(rsid)) = 1""")
    n_panel = con.execute("SELECT count(*) FROM panel").fetchone()[0]
    print(f"panel: {n_panel:,} SNPs\n\n1. case/control frequency direction")
    zs, keys = {}, []
    for src in SOURCES:
        if not local_files(src):
            continue
        print(f"   {src.key:<5} {case_control_check(con, src)}")
        build.load_raw(con, src, "s0")
        # Z aligned to the panel's a1 (a study on the other strand is complemented first)
        comp = build._COMP
        rows = con.execute(f"""
            SELECT f.rsid, CASE
                WHEN s.a1 = f.a1 AND s.a2 = f.a2 THEN 1
                WHEN s.a1 = f.a2 AND s.a2 = f.a1 THEN -1
                WHEN {comp.format(x='s.a1')} = f.a1 AND {comp.format(x='s.a2')} = f.a2 THEN 1
                WHEN {comp.format(x='s.a1')} = f.a2 AND {comp.format(x='s.a2')} = f.a1 THEN -1 END AS sgn,
                s.beta, s.se, s.p
            FROM s0 s JOIN panel f USING (rsid)""").fetchnumpy()
        from statistics import NormalDist
        sgn = np.asarray(rows["sgn"], dtype=float)
        beta = np.asarray(rows["beta"], dtype=float)
        p = np.clip(np.asarray(rows["p"], dtype=float), 1e-300, 1.0)
        nd = NormalDist()
        absz = np.array([-nd.inv_cdf(x / 2) for x in p])        # |Z| from the p-value: works for every release
        zs[src.key] = dict(zip(np.asarray(rows["rsid"]).tolist(), (sgn * np.sign(beta) * absz).tolist()))
        keys.append(src.key)
        con.execute("DROP TABLE s0")
    print("\n2. correlation of Z-scores across the panel (x100); expected positive")
    print("      " + "".join(f"{k:>6}" for k in keys))
    worst = {}
    for a in keys:
        cells = []
        for b in keys:
            if a == b:
                cells.append(f"{'':>6}")
                continue
            common = zs[a].keys() & zs[b].keys()
            x = np.array([zs[a][r] for r in common])
            y = np.array([zs[b][r] for r in common])
            ok = np.isfinite(x) & np.isfinite(y)
            r = float(np.corrcoef(x[ok], y[ok])[0, 1])
            worst.setdefault(a, []).append(r)
            cells.append(f"{100 * r:>6.1f}")
        print(f"{a:<6}" + "".join(cells))
    print(f"\n(standard error of each entry is about {100 / np.sqrt(n_panel):.1f})")
    for a in keys:
        m = float(np.mean(worst[a]))
        print(f"   {a:<5} mean correlation with the others {100 * m:+.1f}" + ("   <-- CHECK" if m < 0 else ""))


if __name__ == "__main__":
    main()
