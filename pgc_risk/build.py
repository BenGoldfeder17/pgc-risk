"""Turn the downloaded PGC shards into one compact weight file per disorder.

For every source the build
  1. checks every row really comes from the intended PGC file,
  2. keeps autosomal, bi-allelic, strand-unambiguous SNPs with an rs number, good imputation
     quality and a near-complete sample size,
  3. attaches one shared European allele frequency (and drops SNPs whose alleles or
     frequency disagree with it),
  4. keeps the SNPs with p <= P_STORE and adds the hg19 reference base and the hg38 position.

Output: data/weights/<key>.npz (+ .json with the QC trail).
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

from .catalog import BY_KEY, FREQ_SOURCES, SOURCES, Source
from .fetch import DATA, local_files

WEIGHTS = DATA / "weights"
BUILD = DATA / "build"
REF_DIR = DATA / "reference" / "ucsc"

P_STORE = 0.05          # SNPs above this p-value are never used, so they are not stored
MIN_INFO = 0.8          # imputation quality
MIN_MAF = 0.01
MIN_N_FRACTION = 0.6    # per-SNP effective sample size relative to the study's 90th percentile
MAX_FREQ_DIFF = 0.2     # study frequency vs shared European reference

CODE = {"A": 0, "C": 1, "G": 2, "T": 3}


def _q(name: str) -> str:
    return '"' + name.replace('"', '""') + '"'


def _connect():
    import duckdb

    con = duckdb.connect()
    con.execute("SET preserve_insertion_order=false")
    return con


def _se_column(src: Source, columns: list[str]) -> str | None:
    for cand in ("SE", "se", "standard_error", "StdErr"):
        if cand in columns:
            return cand
    return None


def load_raw(con, src: Source, table: str) -> dict:
    """Read one source into a DuckDB table with uniform column names; return row counts."""
    files = local_files(src)
    if not files:
        raise FileNotFoundError(f"[{src.key}] no downloaded shards - run `pgc-risk fetch` first")
    flist = "[" + ",".join("'" + str(f).replace("\\", "/").replace("'", "''") + "'" for f in files) + "]"
    scan = f"read_parquet({flist}, union_by_name=true)"
    columns = [r[0] for r in con.execute(f"DESCRIBE SELECT * FROM {scan}").fetchall()]
    need = [src.snp, src.chr, src.pos, src.a1, src.a2, src.effect, src.p]
    missing = [c for c in need if c not in columns]
    if missing:
        raise RuntimeError(f"[{src.key}] columns {missing} not found; shards have {columns}")
    stats = {"files": len(files)}
    if "_source_file" in columns:
        rows = con.execute(f"SELECT _source_file, count(*) FROM {scan} GROUP BY 1").fetchall()
        if [r[0] for r in rows] != [src.source_file]:
            raise RuntimeError(f"[{src.key}] expected rows from {src.source_file!r}, found {rows}")
    se = _se_column(src, columns)
    dbl = lambda c: f"TRY_CAST({_q(c)} AS DOUBLE)"
    if src.effect_kind == "or":
        beta = f"CASE WHEN {dbl(src.effect)} > 0 THEN ln({dbl(src.effect)}) END"
    else:
        beta = dbl(src.effect)        # log-odds, or a Z score that is rescaled once frequencies are known
    sel = f"""
        SELECT TRY_CAST(regexp_extract(CAST({_q(src.snp)} AS VARCHAR), '^rs([0-9]+)', 1) AS BIGINT) AS rsid,
               TRY_CAST(regexp_replace(CAST({_q(src.chr)} AS VARCHAR), '\\.0$', '') AS INTEGER) AS chr,
               TRY_CAST({dbl(src.pos)} AS BIGINT) AS pos,
               upper(CAST({_q(src.a1)} AS VARCHAR)) AS a1,
               upper(CAST({_q(src.a2)} AS VARCHAR)) AS a2,
               {beta} AS beta,
               {dbl(se) if se and src.effect_kind != 'z' else 'NULL::DOUBLE'} AS se,
               {dbl(src.p)} AS p,
               {dbl(src.info) if src.info else 'NULL::DOUBLE'} AS info,
               {dbl(src.study_freq) if src.study_freq else 'NULL::DOUBLE'} AS sfreq,
               {dbl(src.n_cases) if src.n_cases else 'NULL::DOUBLE'} AS ncas,
               {dbl(src.n_controls) if src.n_controls else 'NULL::DOUBLE'} AS ncon
        FROM {scan}"""
    con.execute(f"CREATE OR REPLACE TABLE {table}_all AS {sel}")
    stats["rows"] = con.execute(f"SELECT count(*) FROM {table}_all").fetchone()[0]
    # Basic per-SNP filters that need nothing but the row itself.
    con.execute(f"""
        CREATE OR REPLACE TABLE {table}_snv AS
        SELECT * FROM {table}_all
        WHERE rsid IS NOT NULL AND chr BETWEEN 1 AND 22 AND pos > 0
          AND a1 IN ('A','C','G','T') AND a2 IN ('A','C','G','T') AND a1 <> a2
          AND NOT ((a1='A' AND a2='T') OR (a1='T' AND a2='A') OR (a1='C' AND a2='G') OR (a1='G' AND a2='C'))
          AND beta IS NOT NULL AND isfinite(beta) AND p IS NOT NULL AND p >= 0 AND p <= 1""")
    stats["autosomal_unambiguous_snvs"] = con.execute(f"SELECT count(*) FROM {table}_snv").fetchone()[0]
    con.execute(f"""
        CREATE OR REPLACE TABLE {table} AS
        SELECT * FROM {table}_snv
        WHERE (info IS NULL OR info >= {MIN_INFO})
          AND rsid NOT IN (SELECT rsid FROM {table}_snv GROUP BY rsid HAVING count(*) > 1)""")
    stats["after_info_and_duplicates"] = con.execute(f"SELECT count(*) FROM {table}").fetchone()[0]
    con.execute(f"DROP TABLE {table}_all; DROP TABLE {table}_snv")
    return stats


# SQL fragment: frequency of the row's a1 given a reference row (r1, r2, rf = frequency of r1).
# Strand flips are resolved by complementing; anything else is an allele mismatch (NULL).
_COMP = "CASE {x} WHEN 'A' THEN 'T' WHEN 'T' THEN 'A' WHEN 'C' THEN 'G' ELSE 'C' END"
_ALIGN = f"""CASE
    WHEN s.a1 = f.a1 AND s.a2 = f.a2 THEN f.freq
    WHEN s.a1 = f.a2 AND s.a2 = f.a1 THEN 1 - f.freq
    WHEN {_COMP.format(x='s.a1')} = f.a1 AND {_COMP.format(x='s.a2')} = f.a2 THEN f.freq
    WHEN {_COMP.format(x='s.a1')} = f.a2 AND {_COMP.format(x='s.a2')} = f.a1 THEN 1 - f.freq
    END"""


def build_frequency_table(con, log=print) -> dict:
    """One European allele-frequency per SNP, merged from the releases that publish one."""
    BUILD.mkdir(parents=True, exist_ok=True)
    out = BUILD / "eur_freq.parquet"
    info = {"sources": []}
    con.execute("CREATE OR REPLACE TABLE eurfreq (rsid BIGINT, chr INTEGER, pos BIGINT, a1 VARCHAR, "
                "a2 VARCHAR, freq DOUBLE, src VARCHAR)")
    for key in FREQ_SOURCES:
        src = BY_KEY[key]
        t0 = time.time()
        load_raw(con, src, "fsrc")
        before = con.execute("SELECT count(*) FROM eurfreq").fetchone()[0]
        if before:
            # agreement with what is already in the table, as a check on allele coding
            r = con.execute(f"""
                SELECT count(*), corr(s.sfreq, {_ALIGN}), avg(abs(s.sfreq - ({_ALIGN}))),
                       sum(CASE WHEN ({_ALIGN}) IS NULL THEN 1 ELSE 0 END)
                FROM fsrc s JOIN eurfreq f USING (rsid) WHERE s.sfreq IS NOT NULL""").fetchone()
            agree = {"shared": r[0], "freq_correlation": r[1], "mean_abs_diff": r[2], "allele_mismatch": r[3]}
        else:
            agree = {}
        con.execute(f"""
            INSERT INTO eurfreq
            SELECT rsid, chr, pos, a1, a2, sfreq, '{key}' FROM fsrc s
            WHERE sfreq > 0 AND sfreq < 1 AND rsid NOT IN (SELECT rsid FROM eurfreq)""")
        after = con.execute("SELECT count(*) FROM eurfreq").fetchone()[0]
        info["sources"].append({"key": key, "column": src.study_freq, "added": after - before, **agree})
        log(f"[freq] {key}: +{after - before:,} SNPs ({time.time() - t0:.0f}s) {agree}")
        con.execute("DROP TABLE fsrc")
    info["snps"] = con.execute("SELECT count(*) FROM eurfreq").fetchone()[0]
    con.execute(f"COPY eurfreq TO '{str(out).replace(chr(92), '/')}' (FORMAT parquet)")
    (BUILD / "eur_freq.json").write_text(json.dumps(info, indent=1))
    return info


def _ensure_freq(con, log=print) -> None:
    f = BUILD / "eur_freq.parquet"
    if f.exists():
        con.execute(f"CREATE OR REPLACE TABLE eurfreq AS SELECT * FROM read_parquet('{str(f).replace(chr(92), '/')}')")
    else:
        build_frequency_table(con, log)


def build_source(con, src: Source, ref=None, lift=None, log=print, p_store: float = P_STORE,
                 out_dir: Path = WEIGHTS, files_override=None) -> dict:
    t0 = time.time()
    stats = load_raw(con, src, "s0")
    # shared European frequency; drops allele mismatches and position disagreements
    con.execute(f"""
        CREATE OR REPLACE TABLE s1 AS
        SELECT s.*, {_ALIGN} AS freq, (s.chr = f.chr AND s.pos = f.pos) AS same_pos
        FROM s0 s JOIN eurfreq f USING (rsid)""")
    stats["with_reference_frequency"] = con.execute("SELECT count(*) FROM s1").fetchone()[0]
    r = con.execute("""SELECT sum(CASE WHEN freq IS NULL THEN 1 ELSE 0 END),
                              sum(CASE WHEN NOT same_pos THEN 1 ELSE 0 END),
                              corr(sfreq, freq), avg(abs(sfreq - freq)) FROM s1""").fetchone()
    stats["allele_mismatch"] = int(r[0] or 0)
    stats["position_mismatch"] = int(r[1] or 0)
    stats["study_vs_reference_freq_corr"] = r[2]
    stats["study_vs_reference_freq_mean_abs_diff"] = r[3]
    if src.effect_kind == "z":
        # log-odds per allele from a Z score: beta = z * se, se = 2 / sqrt(2pq * Neff)
        con.execute("""
            CREATE OR REPLACE TABLE s2 AS
            SELECT * REPLACE (2.0 / sqrt(2*freq*(1-freq) * neff) AS se,
                              beta * 2.0 / sqrt(2*freq*(1-freq) * neff) AS beta)
            FROM (SELECT *, 4.0 / (1.0/ncas + 1.0/ncon) AS neff FROM s1
                  WHERE freq > 0 AND freq < 1 AND ncas > 0 AND ncon > 0)""")
    else:
        con.execute("""
            CREATE OR REPLACE TABLE s2 AS
            SELECT *, CASE WHEN se > 0 THEN 4.0 / (se*se * 2*freq*(1-freq)) END AS neff
            FROM s1 WHERE freq > 0 AND freq < 1""")
    n90 = con.execute("SELECT quantile_cont(neff, 0.9) FROM s2 WHERE neff IS NOT NULL").fetchone()[0]
    stats["effective_n_p90"] = n90
    con.execute(f"""
        CREATE OR REPLACE TABLE s3 AS
        SELECT * FROM s2
        WHERE same_pos AND least(freq, 1 - freq) >= {MIN_MAF}
          AND (sfreq IS NULL OR abs(sfreq - freq) <= {MAX_FREQ_DIFF})
          AND neff IS NOT NULL AND neff >= {MIN_N_FRACTION} * {n90}
          AND abs(beta) < 1""")
    stats["after_qc"] = con.execute("SELECT count(*) FROM s3").fetchone()[0]
    # genomic inflation of the test statistics, for the record
    med = con.execute("SELECT median(p) FROM s3").fetchone()[0]
    stats["median_p"] = med
    stats["genome_wide_significant"] = con.execute("SELECT count(*) FROM s3 WHERE p < 5e-8").fetchone()[0]
    rows = con.execute(f"""
        SELECT rsid, chr, pos, a1, a2, beta, se, p, freq FROM s3 WHERE p <= {p_store}
        ORDER BY chr, pos""").fetchnumpy()
    stats["stored"] = int(len(rows["rsid"]))
    stats["p_store"] = p_store
    n = stats["stored"]
    chrom = rows["chr"].astype(np.int8)
    pos = rows["pos"].astype(np.int32)
    lut = np.full(256, -1, dtype=np.int8)
    for b, c in CODE.items():
        lut[ord(b)] = c
    a1 = lut[np.frombuffer("".join(rows["a1"]).encode("ascii"), dtype=np.uint8)]
    a2 = lut[np.frombuffer("".join(rows["a2"]).encode("ascii"), dtype=np.uint8)]
    refbase = np.full(n, -1, dtype=np.int8)
    pos38 = np.full(n, -1, dtype=np.int32)
    if ref is not None or lift is not None:
        for c in range(1, 23):
            m = chrom == c
            if not m.any():
                continue
            if ref is not None:
                refbase[m] = ref.bases(c, pos[m])
            if lift is not None:
                pos38[m] = lift.lift(c, pos[m])
        if ref is not None:
            on_ref = (refbase == a1) | (refbase == a2)
            flipped = (refbase == 3 - a1) | (refbase == 3 - a2)
            stats["reference_base_matches_an_allele"] = int(on_ref.sum())
            stats["reference_base_matches_other_strand"] = int((~on_ref & flipped).sum())
        if lift is not None:
            stats["lifted_to_hg38"] = int((pos38 > 0).sum())
    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez(out_dir / f"{src.key}.npz",
             rsid=rows["rsid"].astype(np.int64), chrom=chrom, pos=pos, pos38=pos38,
             a1=a1, a2=a2, ref=refbase,
             beta=rows["beta"].astype(np.float32), se=rows["se"].astype(np.float32),
             p=rows["p"].astype(np.float64), freq=rows["freq"].astype(np.float32))
    stats["seconds"] = round(time.time() - t0, 1)
    meta = {"key": src.key, "name": src.name, "repo": src.repo, "revision": src.revision,
            "config": src.config or src.file, "source_file": src.source_file,
            "citation": src.citation, "ref": src.ref, "cases": src.cases, "controls": src.controls,
            "notes": src.notes, "qc": stats}
    (out_dir / f"{src.key}.json").write_text(json.dumps(meta, indent=1))
    log(f"[{src.key}] rows={stats['rows']:,} qc={stats['after_qc']:,} stored={n:,} "
        f"gws={stats['genome_wide_significant']:,} allele_mismatch={stats['allele_mismatch']:,} "
        f"freq_corr={stats['study_vs_reference_freq_corr']} ({stats['seconds']}s)")
    for t in ("s0", "s1", "s2", "s3"):
        con.execute(f"DROP TABLE IF EXISTS {t}")
    return meta


def load_reference(log=print):
    """hg19 bases and the hg19->hg38 chain (downloaded from UCSC on first use)."""
    from .reference import LiftOver, TwoBit, ensure

    try:
        two, chain = ensure(REF_DIR, log)
    except OSError as e:
        log(f"[build] could not download the UCSC reference files ({e}): sequencing VCFs will not be "
            "scorable; array files are unaffected")
        return None, None
    return TwoBit(two), LiftOver(chain)


def build(keys: list[str] | None = None, log=print) -> None:
    con = _connect()
    _ensure_freq(con, log)
    ref, lift = load_reference(log)
    for src in SOURCES:
        if keys and src.key not in keys:
            continue
        build_source(con, src, ref, lift, log)
    log("build complete")


if __name__ == "__main__":
    build(sys.argv[1:] or None)
