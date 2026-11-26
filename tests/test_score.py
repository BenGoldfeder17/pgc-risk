"""Matching, clumping, standardisation and the risk model, on synthetic weights."""
from statistics import NormalDist

import numpy as np
import pytest

from pgc_risk.genome import Genome, read_genome
from pgc_risk.score import (GenomeIndex, Weights, ancestry_fit, clump, liability_risk,
                            score_disorder)

A, C, G, T = 0, 1, 2, 3
BASES = "ACGT"


def make_weights(n=4000, seed=1, spacing=400_000):
    """n independent SNPs spread over the autosomes, unambiguous alleles, with hg38 positions."""
    rng = np.random.default_rng(seed)
    per = n // 22
    chrom = np.repeat(np.arange(1, 23), per).astype(np.int8)
    pos = np.tile(np.arange(1, per + 1) * spacing, 22).astype(np.int32)
    n = len(chrom)
    pairs = np.array([(A, C), (A, G), (C, A), (C, T), (G, A), (G, T), (T, C), (T, G)], dtype=np.int8)
    pick = pairs[rng.integers(0, len(pairs), n)]
    a1, a2 = pick[:, 0].copy(), pick[:, 1].copy()
    ref = np.where(rng.random(n) < 0.5, a1, a2).astype(np.int8)
    return Weights("syn", {"name": "Synthetic"}, rsid=np.arange(1, n + 1, dtype=np.int64) * 7,
                   chrom=chrom, pos=pos, pos38=(pos + 12345).astype(np.int32), a1=a1, a2=a2, ref=ref,
                   beta=rng.normal(0, 0.03, n), p=rng.uniform(1e-9, 0.05, n),
                   freq=rng.uniform(0.05, 0.95, n))


def sample_dosage(w, rng):
    return rng.binomial(2, w.freq).astype(np.int8)


def genome_from_dosage(w, dosage, build=37, with_rsid=True, strand_flip=None):
    """An array-style genome (every SNP listed) carrying the given effect-allele counts."""
    u1 = np.where(dosage >= 1, w.a1, w.a2).astype(np.int8)
    u2 = np.where(dosage == 2, w.a1, w.a2).astype(np.int8)
    if strand_flip is not None:
        u1 = np.where(strand_flip, 3 - u1, u1).astype(np.int8)
        u2 = np.where(strand_flip, 3 - u2, u2).astype(np.int8)
    pos = w.pos if build == 37 else w.pos38
    rsid = w.rsid.copy() if with_rsid else np.full(w.n, -1, dtype=np.int64)
    return Genome("mem", "23andme", build, rsid, w.chrom.copy(), pos.copy(), u1, u2)


def test_dosage_roundtrip_by_rsid_position_and_strand():
    w = make_weights()
    rng = np.random.default_rng(2)
    d = sample_dosage(w, rng)
    for kwargs in ({}, {"with_rsid": False}, {"with_rsid": False, "build": 38},
                   {"strand_flip": rng.random(w.n) < 0.3}):
        avail, dosage, info = GenomeIndex(genome_from_dosage(w, d, **kwargs)).match(w)
        assert avail.all(), kwargs
        assert np.array_equal(dosage, d), kwargs
    assert info["strand_flipped"] > 0


def test_no_calls_and_allele_mismatches_are_unavailable():
    w = make_weights()
    d = sample_dosage(w, np.random.default_rng(3))
    g = genome_from_dosage(w, d)
    g.a1[:10] = -1
    g.a2[:10] = -1
    for i in range(10, 20):                     # a pair that fits the GWAS alleles on neither strand
        g.a1[i] = w.a1[i]
        g.a2[i] = 3 - w.a1[i]
    avail, dosage, info = GenomeIndex(g).match(w)
    assert not avail[:20].any() and avail[20:].all()
    assert info["no_call"] == 10 and info["allele_mismatch"] == 10


def test_clump_keeps_best_snp_and_enforces_spacing():
    rng = np.random.default_rng(4)
    n = 20_000
    w = make_weights(n=22)
    w.chrom = np.ones(n, dtype=np.int8)
    w.pos = np.sort(rng.integers(100_000_000, 150_000_000, n)).astype(np.int32)   # no special region
    w.p = rng.uniform(1e-8, 0.05, n)
    sel = clump(w, np.arange(n), window=250_000)
    kept = np.sort(w.pos[sel])
    assert np.diff(kept).min() >= 250_000 - 1000
    assert int(np.argmin(w.p)) in sel.tolist()
    # every dropped SNP has a more significant kept SNP within the window
    dropped = np.setdiff1d(np.arange(n), sel)
    for i in rng.choice(dropped, 300, replace=False):
        near = sel[np.abs(w.pos[sel].astype(np.int64) - int(w.pos[i])) <= 251_000]
        assert len(near) and w.p[near].min() <= w.p[i]


def test_clump_takes_one_snp_from_the_mhc():
    n = 500
    w = make_weights(n=22)
    w.chrom = np.full(n, 6, dtype=np.int8)
    w.pos = np.linspace(20_000_000, 40_000_000, n).astype(np.int32)
    w.p = np.random.default_rng(5).uniform(1e-8, 0.05, n)
    sel = clump(w, np.arange(n), window=250_000)
    in_mhc = (w.pos[sel] >= 25_000_000) & (w.pos[sel] <= 34_000_000)
    assert in_mhc.sum() == 1
    assert (~in_mhc).sum() > 20


def test_liability_risk_properties():
    nd = NormalDist()
    for k, r2 in ((0.01, 0.07), (0.15, 0.04), (0.02, 0.0)):
        risks = [liability_risk(z, k, r2) for z in (-3, -1, 0, 1, 3)]
        assert all(b >= a for a, b in zip(risks, risks[1:]))
        # averaging the risk over the population's scores must give back the prevalence
        zs = np.linspace(-8, 8, 16001)
        dens = np.array([nd.pdf(z) for z in zs])
        avg = np.trapezoid(np.array([liability_risk(z, k, r2) for z in zs]) * dens, zs)
        assert avg == pytest.approx(k, rel=1e-4)
    assert liability_risk(2.5, 0.02, 0.0) == 0.02
    # worked example: K = 1%, R2 = 7.3%, z = +2  ->  1 - Phi((2.3263 - 0.5404) / 0.9628)
    assert liability_risk(2.0, 0.01, 0.073) == pytest.approx(1 - nd.cdf((2.326348 - 0.540370) / 0.962808),
                                                             rel=1e-5)


def test_scores_of_simulated_people_are_standard_normal():
    w = make_weights(n=6600)
    rng = np.random.default_rng(6)
    zs = []
    for _ in range(300):
        res, _ = score_disorder(GenomeIndex(genome_from_dosage(w, sample_dosage(w, rng))), w)
        zs.append(res.z)
    zs = np.array(zs)
    assert abs(zs.mean()) < 0.2 and 0.88 < zs.std() < 1.12
    assert 0.9 * w.n < res.snps_used < w.n   # only the long-range-LD regions are thinned


def test_score_matches_hand_computation():
    w = make_weights(n=220)
    d = sample_dosage(w, np.random.default_rng(7))
    res, _ = score_disorder(GenomeIndex(genome_from_dosage(w, d)), w)
    raw = float((w.beta * d).sum())
    mean = float((2 * w.freq * w.beta).sum())
    sd = float(np.sqrt((2 * w.freq * (1 - w.freq) * w.beta ** 2).sum()))
    assert res.raw_score == pytest.approx(raw) and res.z == pytest.approx((raw - mean) / sd)


def _write_vcf(path, w, dosage, build, with_rsid, contig_header=True, extra=""):
    pos = w.pos if build == 37 else w.pos38
    length = 249250621 if build == 37 else 248956422
    lines = ["##fileformat=VCFv4.2"]
    if contig_header:
        lines.append(f"##contig=<ID=chr1,length={length}>")
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tME")
    order = np.lexsort((pos, w.chrom))
    for i in order:
        ref = int(w.ref[i])
        alt = int(w.a2[i] if ref == w.a1[i] else w.a1[i])
        n_alt = int(dosage[i]) if alt == w.a1[i] else 2 - int(dosage[i])
        if n_alt == 0:
            continue                               # variants-only: reference sites are not listed
        gt = "0/1" if n_alt == 1 else "1/1"
        rid = f"rs{w.rsid[i]}" if with_rsid else "."
        lines.append(f"chr{w.chrom[i]}\t{pos[i]}\t{rid}\t{BASES[ref]}\t{BASES[alt]}\t60\tPASS\t.\tGT\t{gt}")
    path.write_text("\n".join(lines) + "\n" + extra)


@pytest.mark.parametrize("build,with_rsid,contig", [(37, True, True), (37, False, True),
                                                    (38, False, True), (38, False, False),
                                                    (37, False, False)])
def test_variants_only_vcf_scores_like_the_full_genotype(tmp_path, build, with_rsid, contig):
    w = make_weights(n=2200)
    d = sample_dosage(w, np.random.default_rng(8))
    want, _ = score_disorder(GenomeIndex(genome_from_dosage(w, d)), w)
    p = tmp_path / "me.vcf"
    _write_vcf(p, w, d, build, with_rsid, contig_header=contig)
    g = read_genome(p)
    assert g.variants_only and g.build == (build if contig else None)
    idx = GenomeIndex(g)
    got, detail = score_disorder(idx, w)
    assert idx.pos_build == build
    assert got.snps_used == want.snps_used
    assert got.z == pytest.approx(want.z)
    assert got.match["assumed_reference"] == int((np.where(w.ref == w.a1, 2, 0) == d).sum())


def test_variants_only_vcf_does_not_assume_reference_under_another_record(tmp_path):
    w = make_weights(n=220)
    d = np.where(w.ref == w.a1, 2, 0).astype(np.int8)       # homozygous reference everywhere
    p = tmp_path / "me.vcf"
    i = 5                                                    # a filtered call sits on SNP i
    alt = BASES[int(w.a2[i] if w.ref[i] == w.a1[i] else w.a1[i])]
    extra = f"chr{w.chrom[i]}\t{w.pos[i]}\t.\t{BASES[int(w.ref[i])]}\t{alt}\t3\tLowQual\t.\tGT\t0/1\n"
    extra += "chr1\t999\t.\tA\tG\t60\tPASS\t.\tGT\t0/1\n"   # so the file is not empty
    _write_vcf(p, w, d, 37, False, extra=extra)
    avail, dosage, info = GenomeIndex(read_genome(p)).match(w)
    assert not avail[i] and avail.sum() == w.n - 1
    assert np.array_equal(dosage[avail], d[avail])


def test_ancestry_fit_is_one_for_matching_population_and_high_otherwise():
    w = make_weights(n=6600)
    rng = np.random.default_rng(9)
    same = sample_dosage(w, rng)
    fit, n = ancestry_fit([(w.rsid, same, w.freq), (w.rsid, same, w.freq)])   # duplicates collapse
    assert n == w.n and 0.93 < fit < 1.07
    # a population whose frequencies drifted away from the reference
    drift = np.clip(w.freq + rng.normal(0, 0.2, w.n), 0.02, 0.98)
    other = rng.binomial(2, drift).astype(np.int8)
    fit2, _ = ancestry_fit([(w.rsid, other, w.freq)])
    assert fit2 > 1.25
