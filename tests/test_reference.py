"""Checks of the reference helpers against the real UCSC files (skipped if not downloaded)."""
from pathlib import Path

import numpy as np
import pytest

from pgc_risk.reference import LiftOver, TwoBit

REF = Path(__file__).resolve().parent.parent / "data" / "reference" / "ucsc"
needs_2bit = pytest.mark.skipif(not (REF / "hg19.2bit").exists(), reason="hg19.2bit not downloaded")
needs_chain = pytest.mark.skipif(not (REF / "hg19ToHg38.over.chain.gz").exists(), reason="chain not downloaded")

# Well-known SNPs: (rs, chromosome, hg19 position, hg19 reference base, hg38 position)
KNOWN = [
    ("rs429358", 19, 45411941, "T", 44908684),    # APOE e4
    ("rs7412", 19, 45412079, "C", 44908822),      # APOE e2
    ("rs1801133", 1, 11856378, "G", 11796321),    # MTHFR C677T
    ("rs4680", 22, 19951271, "G", 19963748),      # COMT Val158Met
    ("rs6265", 11, 27679916, "C", 27658369),      # BDNF Val66Met
    ("rs9939609", 16, 53820527, "T", 53786615),   # FTO
]
CODE = {"A": 0, "C": 1, "G": 2, "T": 3}


@needs_2bit
def test_twobit_reference_bases():
    tb = TwoBit(REF / "hg19.2bit")
    for rs, c, pos, base, _ in KNOWN:
        assert tb.bases(c, np.array([pos]))[0] == CODE[base], rs
    # start of chr1 is an N block; position 0 and past-the-end are out of range
    assert tb.bases(1, np.array([1, 5000, 0, 400_000_000])).tolist() == [-1, -1, -1, -1]


@needs_chain
def test_liftover_known_positions():
    lo = LiftOver(REF / "hg19ToHg38.over.chain.gz")
    for rs, c, pos, _, pos38 in KNOWN:
        assert lo.lift(c, np.array([pos]))[0] == pos38, rs
