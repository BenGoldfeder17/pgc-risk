"""Which PGC GWAS feeds each disorder, where it lives on Hugging Face, and how to read it.

The Psychiatric Genomics Consortium (PGC) summary statistics were published on Hugging Face
as `OpenMed/pgc-*`.  Those repos are no longer publicly reachable, so this catalog points at
the byte-identical public mirrors (`introvoyz041/pgc-*`; `tunahanf/pgc-*` carries the same
sha256 for every file) pinned to a commit, plus one harmonized re-release for the only
eating-disorder GWAS that survives anywhere on the Hub.

Every entry is a European-ancestry (or, for schizophrenia, ~80% European) case-control GWAS
with genome-wide coverage.  Releases whose Hugging Face conversion is broken (the PGC
"sumstats VCF" files: bip2021-all, ptsd2024, panic2019, hoarding2022, scz2022-european) are
not usable and are skipped in favour of the newest release that is.
"""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Source:
    key: str                      # short id used on the command line and for file names
    name: str                     # display name
    repo: str                     # Hugging Face dataset repo
    revision: str                 # pinned commit
    # Shard series: data/<config>/train-NNNNN-of-<series>.parquet, optionally a sub-range.
    config: str = ""
    series: int = 0
    first: int | None = None
    last: int | None = None
    file: str | None = None       # or a single parquet file instead of a shard series
    source_file: str = ""         # PGC's original file name (checked against `_source_file`)
    # Column names in the parquet shards.
    snp: str = "SNP"
    chr: str = "CHR"
    pos: str = "BP"
    a1: str = "A1"                # effect allele
    a2: str = "A2"
    effect: str = "OR"
    effect_kind: str = "or"       # "or" | "beta" (log-odds) | "z" (needs n_cases/n_controls)
    p: str = "P"
    info: str | None = "INFO"
    study_freq: str | None = None  # allele frequency of a1 in the study (QC cross-check only)
    n_cases: str | None = None     # only needed for effect_kind == "z"
    n_controls: str | None = None
    # What the release is.
    citation: str = ""
    ref: str = ""                 # DOI or PubMed id of the paper
    cases: int = 0                # where the release states them
    controls: int = 0
    notes: str = ""


_M = "introvoyz041"
_REV = {
    "pgc-adhd": "7a58dc9aad1401b038ba557466efdaccf1f88743",
    "pgc-anxiety": "3d6582a073d795bd952cd0ce5ef9a4d736fcb7ec",
    "pgc-autism": "72daadaa31a6dd1e31ce071d8db22a0b29ad04ca",
    "pgc-bipolar": "69e42efb44c1119897f31442ce9b627e23d85558",
    "pgc-mdd": "85ec296f372e65b40a74cf4b86515ada63ac31b6",
    "pgc-ocd-tourette": "fd68b66ee86734b17b357898da82e83a15493734",
    "pgc-other": "18a5025a2deb022abf7db94f232cbf394d241a03",
    "pgc-ptsd": "351cfa7cca2303a8ba0db21219029ae9874fa731",
    "pgc-schizophrenia": "1de77741028afb6e6bed53b77c9daf97c6369cd7",
    "pgc-substance-use": "fe551dbe6193017d7702ea93211f620c64bacf07",
}


def _mirror(repo: str, **kw) -> dict:
    return dict(repo=f"{_M}/{repo}", revision=_REV[repo], **kw)


SOURCES: list[Source] = [
    Source(
        key="scz", name="Schizophrenia",
        **_mirror("pgc-schizophrenia", config="scz2022", series=764),
        source_file="daner_PGC_SCZ_w3_90_0418b_ukbbdedupe.gz",
        study_freq="FRQ_U_93456",
        citation="Trubetskoy et al. 2022, Nature 604:502 (PGC3 core GWAS)",
        ref="doi:10.1038/s41586-022-04434-5", cases=67390, controls=94015,
        notes="Primary PGC3 meta-analysis: about 80% European and 20% East Asian samples.",
    ),
    Source(
        key="bip", name="Bipolar disorder",
        **_mirror("pgc-bipolar", config="bip2024", series=694),
        source_file="bip2024_eur_no23andMe.gz",
        study_freq="HRC_FRQ_A1",
        citation="O'Connell et al. 2025, Nature 639:968 (European ancestry, without 23andMe)",
        ref="doi:10.1038/s41586-024-08468-9",
    ),
    Source(
        key="mdd", name="Major depression",
        **_mirror("pgc-mdd", config="mdd2025", series=737),
        source_file="pgc-mdd2025_no23andMe_eur_v3-49-24-11_formatted.tsv.gz",
        snp="rsid", chr="chromosome", pos="base_pair_location",
        a1="effect_allele", a2="other_allele", effect="beta", effect_kind="beta",
        p="p_value", info="info", study_freq="effect_allele_frequency",
        citation="Adams et al. 2025, Cell (European ancestry, without 23andMe)",
        ref="doi:10.1016/j.cell.2024.12.002",
        notes="One of the 737 shards is missing from the Hugging Face copy (about 0.1% of SNPs).",
    ),
    Source(
        key="anx", name="Anxiety disorders",
        **_mirror("pgc-anxiety", config="anx2026", series=724),
        source_file="ANX_2026_daner_fullANX_v12_woUTAH_11022026.gz",
        study_freq="FRQ_U_729602",
        citation="Strom et al. 2026, Nature Genetics 58:275 (any anxiety disorder)",
        ref="PMID 41634414", cases=122083, controls=729602,
    ),
    Source(
        key="adhd", name="ADHD",
        **_mirror("pgc-adhd", config="adhd2022", series=678),
        source_file="ADHD2022_iPSYCH_deCODE_PGC.meta.gz",
        study_freq="FRQ_U_186843",
        citation="Demontis et al. 2023, Nature Genetics", ref="PMID 36702997",
        cases=38691, controls=186843,
    ),
    Source(
        key="asd", name="Autism spectrum disorder",
        **_mirror("pgc-autism", config="asd2019", series=912),
        source_file="iPSYCH-PGC_ASD_Nov2017.gz",
        citation="Grove et al. 2019, Nature Genetics", ref="PMID 30804558",
        cases=18381, controls=27969,
    ),
    Source(
        key="ptsd", name="PTSD",
        **_mirror("pgc-ptsd", config="ptsd2019", series=926, first=539, last=734),
        source_file="pts_eur_freeze2_overall.results.gz",
        study_freq="FRQ_U_151447",
        citation="Nievergelt et al. 2019, Nature Communications (freeze 2, European ancestry)",
        ref="doi:10.1038/s41467-019-12576-w", cases=23212, controls=151447,
        notes="The 2024 freeze-3 release is unreadable in the Hugging Face copy.",
    ),
    Source(
        key="ocd", name="Obsessive-compulsive disorder",
        **_mirror("pgc-ocd-tourette", config="ocd2025", series=680),
        source_file="daner_OCD_full_wo23andMe_190522.gz",
        study_freq="FRQ_U_1114613",
        citation="Strom et al. 2025, Nature Genetics (without 23andMe)", ref="PMID 40360802",
        cases=23493, controls=1114613,
    ),
    Source(
        key="ts", name="Tourette syndrome",
        **_mirror("pgc-ocd-tourette", config="ts2019", series=827),
        source_file="TS_Oct2018.gz",
        citation="Yu et al. 2019, American Journal of Psychiatry", ref="PMID 30818990",
        cases=4819, controls=9488,
    ),
    Source(
        key="an", name="Anorexia nervosa",
        repo="lighteternal/pgc-psychiatric-gwas-harmonized",
        revision="6c3241136562df51a4e942d84a466a4202866c99",
        file="data/eating_disorders.parquet",
        source_file="OpenMed/pgc-eating-disorders:an2017",
        snp="snp_id", chr="chr", pos="pos", a1="effect_allele", a2="other_allele",
        effect="effect_size", effect_kind="beta", p="pval", info="info",
        citation="Duncan et al. 2017, American Journal of Psychiatry", ref="PMID 28494655",
        cases=3495, controls=10982,
        notes="Only the small 2017 GWAS survives on the Hub (via a harmonized re-release).",
    ),
    Source(
        key="bpd", name="Borderline personality disorder",
        **_mirror("pgc-other", config="bpd2025", series=605),
        source_file="daner_bpd15_2025_all_24b.gz",
        study_freq="FRQ_U_1041717",
        citation="Streit et al. 2026, Nature Genetics", ref="PMID 42477112",
        cases=12339, controls=1041717,
    ),
    Source(
        key="alc", name="Alcohol dependence",
        **_mirror("pgc-substance-use", config="sud2018-alc", series=10893, first=7597, last=8511),
        source_file="pgc_alcdep.eur_unrel_genotyped.aug2018_release.txt.gz",
        citation="Walters et al. 2018, Nature Neuroscience (European ancestry, unrelated genotyped subset)",
        ref="PMID 30482948",
        notes="72 of this release's 915 shards are missing from the Hugging Face copy (about 8% of SNPs).",
    ),
    Source(
        key="cud", name="Cannabis use disorder",
        **_mirror("pgc-substance-use", config="sud2020-cud", series=3461, first=2565, last=3460),
        source_file="CUD_EUR_full_public_11.14.2020.gz",
        effect="Z", effect_kind="z", info=None, n_cases="N_CAS", n_controls="N_CON",
        citation="Johnson et al. 2020, Lancet Psychiatry (European ancestry)", ref="PMID 33096046",
        notes="Released as Z-scores; effect sizes are reconstructed from Z, N and allele frequency.",
    ),
]

BY_KEY = {s.key: s for s in SOURCES}

# European allele frequencies come from these releases, in priority order.  One shared table
# keeps every disorder on the same population reference, including the releases that ship no
# frequencies at all (autism, Tourette, anorexia, alcohol, cannabis) and schizophrenia, whose
# own frequency column mixes European and East Asian controls.
FREQ_SOURCES = ["bip", "mdd", "anx", "adhd"]


@dataclass(frozen=True)
class RiskModel:
    """Numbers that turn a standardized polygenic score into an absolute risk."""
    prevalence: float             # lifetime risk in the general population (K)
    r2_liability: float           # share of liability variance a polygenic score captures
    prevalence_source: str = ""
    r2_source: str = ""
    low_power: bool = False       # GWAS too small for the score to carry much information
    extra: dict = field(default_factory=dict)


_PAPER = "lifetime risk assumed in the GWAS paper"

# Lifetime prevalence K per disorder.  The score accuracy (r2_liability) is filled in from
# pgc_risk/accuracy.json, which tools/pseudo_r2.py measures on the data itself.
PREVALENCE: dict[str, tuple[float, str]] = {
    "scz": (0.01, _PAPER),
    "bip": (0.02, _PAPER),
    "mdd": (0.15, _PAPER),
    "anx": (0.20, _PAPER),
    "adhd": (0.05, _PAPER),
    "asd": (0.0122, _PAPER),
    "ptsd": (0.068, "US lifetime prevalence, National Comorbidity Survey Replication (Kessler et al. 2005)"),
    "ocd": (0.01, _PAPER),
    "ts": (0.008, _PAPER),
    "an": (0.01, "approximate lifetime prevalence quoted in the GWAS paper"),
    "bpd": (0.015, _PAPER),
    "alc": (0.125, "US lifetime prevalence of DSM-IV alcohol dependence, NESARC (Hasin et al. 2007)"),
    "cud": (0.063, "US lifetime prevalence of DSM-5 cannabis use disorder, NESARC-III (Hasin et al. 2016)"),
}

# A score explaining less than this share of liability variance is flagged as uninformative.
LOW_POWER_R2 = 0.01


def _load_risk() -> dict[str, RiskModel]:
    import json
    from pathlib import Path

    path = Path(__file__).with_name("accuracy.json")
    measured = json.loads(path.read_text()) if path.exists() else {}
    risk = {}
    for key, (k, source) in PREVALENCE.items():
        m = measured.get(key)
        r2 = float(m["r2_liability"]) if m else 0.0
        risk[key] = RiskModel(k, max(r2, 0.0), source, m.get("source", "") if m else "not measured",
                              low_power=r2 < LOW_POWER_R2, extra=m or {})
    return risk


RISK: dict[str, RiskModel] = _load_risk()
