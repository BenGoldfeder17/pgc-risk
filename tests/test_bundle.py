"""The portable reference bundle: file format, strand-free matching and the scoring rule."""
from statistics import NormalDist

import numpy as np
import pytest

from pgc_risk import bundle as B

A, C, G, T = 0, 1, 2, 3


def make(tmp_path, power=0.72):
    """Two clumps on three SNPs.  SNP 0 (A/C) and SNP 1 (A/G) form clump one; SNP 2 (A/G) is clump two."""
    header = {
        "name": "test", "created": "2026-01-01", "description": "test", "snps": 3,
        "ld_step": B.LD_STEP, "freq_scale": B.FREQ_SCALE, "ld_share_power": power,
        "fit_good": 1.04, "fit_poor": 1.075, "min_fit_snps": 0, "low_power_r2": 0.01,
        "disorders": [{"key": "x", "name": "Disorder X", "prevalence": 0.01, "r2": 0.07, "r2_basis": "published",
                       "liability_threshold": NormalDist().inv_cdf(0.99), "beta_scale": 0.001,
                       "full_indep": 2 * 0.3 * 0.7 * 0.1 ** 2 + 2 * 0.5 * 0.5 * 0.05 ** 2}],
    }
    arrays = {
        "snp.rsid": np.array([11, 22, 33], dtype=np.uint32),
        "snp.b": np.array([C, G, G], dtype=np.uint8),
        "snp.freq": np.rint(np.array([0.3, 0.4, 0.5]) * B.FREQ_SCALE).astype(np.uint16),
        "x.idx_delta": np.array([0, 1, 1], dtype=np.int32),           # rows: SNP 0, SNP 1, SNP 2
        "x.beta": np.array([100, 80, -50], dtype=np.int16),            # 0.1, 0.08, -0.05 per A allele
        "x.ld": np.rint(np.array([1.5, 1.2, 1.0]) / B.LD_STEP).astype(np.uint8),
        "x.clump_len": np.array([2, 1], dtype=np.uint32),
    }
    path = tmp_path / "test.pgcr"
    B.write_bundle(path, header, arrays)
    return B.read_bundle(path)


def test_roundtrip(tmp_path):
    b = make(tmp_path)
    assert b.rsid.tolist() == [11, 22, 33] and b.b_allele.tolist() == [C, G, G]
    assert b.freq == pytest.approx([0.3, 0.4, 0.5], abs=1e-4)
    d = b.disorders[0]
    assert d.idx.tolist() == [0, 1, 2] and d.clump_start.tolist() == [0, 2]
    assert d.beta == pytest.approx([0.1, 0.08, -0.05]) and d.ld == pytest.approx([1.5, 1.2, 1.0])
    (tmp_path / "junk.pgcr").write_bytes(b"not a bundle at all")
    with pytest.raises(B.BundleError):
        B.read_bundle(tmp_path / "junk.pgcr")


def test_counts_are_strand_free_and_duplicates_are_dropped(tmp_path):
    b = make(tmp_path)
    rs = np.array([11, 22, 33, 33, 44, -1], dtype=np.int64)
    # SNP 11 is A/C: "TG" is the same call on the other strand (A and C) -> one A
    # SNP 22 is A/G: "AA" -> two A.  SNP 33 is listed twice -> not used.  44 is not in the bundle.
    a1 = np.array([T, A, A, G, A, A], dtype=np.int8)
    a2 = np.array([G, A, G, G, A, A], dtype=np.int8)
    assert B.a_allele_counts(b, rs, a1, a2).tolist() == [1, 2, -1]
    # a pair that fits neither strand (A and T at an A/C SNP), and a missing call
    assert B.a_allele_counts(b, np.array([11, 22]), np.array([A, -1], np.int8), np.array([T, -1], np.int8)).tolist() == [-1, -1, -1]
    # an A/G SNP typed on the other strand as T/C: "CC" is G/G -> no A
    assert B.a_allele_counts(b, np.array([22]), np.array([C], np.int8), np.array([C], np.int8)).tolist() == [-1, 0, -1]


def test_scoring_uses_the_best_available_snp_of_each_clump(tmp_path):
    b = make(tmp_path)
    nd = NormalDist()
    # everything present: clump one is represented by its first SNP (row 0), clump two by SNP 2
    res = B.score_counts(b, np.array([2, 0, 1], dtype=np.int8))["disorders"][0]
    v0, v2 = 2 * 0.3 * 0.7 * 0.1 ** 2, 2 * 0.5 * 0.5 * 0.05 ** 2
    raw, mean = 0.1 * 2 - 0.05 * 1, 2 * 0.3 * 0.1 - 2 * 0.5 * 0.05
    var = (v0 + v2) + 1.0 * (v0 * 0.5 + v2 * 0.0)                     # every clump present: full LD factor
    assert res["snps_used"] == 2 and res["coverage"] == pytest.approx(1.0, abs=1e-3)
    assert res["z"] == pytest.approx((raw - mean) / var ** 0.5, rel=1e-3)
    t = nd.inv_cdf(0.99)
    assert res["probability"] == pytest.approx(1 - nd.cdf((t - 0.07 ** 0.5 * res["z"]) / (1 - 0.07) ** 0.5))
    # the first SNP is missing: clump one falls back to its second SNP; clump two is missing altogether
    res = B.score_counts(b, np.array([-1, 1, -1], dtype=np.int8))["disorders"][0]
    v1 = 2 * 0.4 * 0.6 * 0.08 ** 2
    share = v1 / (v0 + v2)
    var = v1 + share ** 0.72 * (v1 * 0.2)                              # neighbours only partly present
    assert res["snps_used"] == 1 and res["coverage"] == pytest.approx(share, rel=1e-3)
    assert res["z"] == pytest.approx((0.08 * 1 - 2 * 0.4 * 0.08) / var ** 0.5, rel=1e-3)


def test_ancestry_gate(tmp_path):
    b = make(tmp_path)
    b.header["min_fit_snps"] = 2
    ok = B.score_counts(b, np.array([1, 1, 1], dtype=np.int8))
    assert ok["valid"] and ok["ancestry_fit"]["snps"] == 3
    few = B.score_counts(b, np.array([1, -1, -1], dtype=np.int8))
    assert not few["valid"] and "too few" in few["warnings"][0]
    far = B.score_counts(b, np.array([2, 2, 2], dtype=np.int8))        # homozygous everywhere: a poor fit
    assert not far["valid"] and far["ancestry_fit"]["value"] > 1.075
