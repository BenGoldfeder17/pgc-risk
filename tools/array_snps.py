"""List of SNPs that sit on genotyping arrays, for testing how array files score.

The schizophrenia release records, per SNP, in how many of its cohorts the SNP was genotyped
directly rather than imputed (column ngt).  SNPs typed in at least 20 cohorts are on the
common Illumina and Affymetrix chips; that gives about 826,000 autosomal SNPs, comparable to
a consumer array.  Written to data/validation/array_like_rsids.npy.

usage: python tools/array_snps.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pgc_risk.catalog import BY_KEY  # noqa: E402
from pgc_risk.fetch import local_files  # noqa: E402


def main() -> None:
    import duckdb

    files = local_files(BY_KEY["scz"])
    flist = "[" + ",".join("'" + str(f).replace(chr(92), "/") + "'" for f in files) + "]"
    rs = duckdb.connect().execute(f"""
        SELECT DISTINCT rs FROM (
            SELECT TRY_CAST(regexp_extract(SNP, '^rs([0-9]+)', 1) AS BIGINT) AS rs
            FROM read_parquet({flist}, union_by_name=true)
            WHERE TRY_CAST(ngt AS INTEGER) >= 20 AND CHR BETWEEN 1 AND 22)
        WHERE rs IS NOT NULL""").fetchnumpy()["rs"]
    out = ROOT / "data" / "validation" / "array_like_rsids.npy"
    out.parent.mkdir(parents=True, exist_ok=True)
    np.save(out, np.sort(np.asarray(rs, dtype=np.int64)))
    print(f"{len(rs):,} array-typed autosomal SNPs -> {out}")


if __name__ == "__main__":
    main()
