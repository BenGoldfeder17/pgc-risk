"""The portable reference bundle: one file holding everything needed to score an array genome.

The full pipeline (`pgc-risk fetch` + `build`) downloads about 20 GB and scores with an LD
panel of 4 million SNPs.  The bundle is the same information boiled down so that it can be
scored anywhere, including in a web browser, without that panel:

* LD clumping is done once, when the bundle is built.  Every candidate SNP is filed under the
  index SNP of its clump, most significant first.  At scoring time each clump is represented
  by its most significant SNP that the genome file actually has.
* Each SNP carries an LD factor: how much its contribution to the score's variance grows
  because of its correlation with the neighbouring clumps (measured on the panel, with the
  neighbours' index SNPs standing in for whatever represents them).

Layout: gzip of  b"PGCRISK1" | uint32 header length | JSON header | arrays (8-byte aligned).
The header names every array with its type, offset and length, and carries the disorders,
their prevalence and accuracy figures, the ancestry thresholds and the data attribution, so a
reader needs no built-in knowledge of what the bundle contains.

SNPs are keyed by rs number and stored strand-free: every non-palindromic SNP is A/C or A/G
on one of the two strands, so one flag (C or G) fixes its alleles, and effects and
frequencies refer to the A allele.
"""
from __future__ import annotations

import gzip
import json
import struct
from dataclasses import dataclass
from pathlib import Path
from statistics import NormalDist

import numpy as np

MAGIC = b"PGCRISK1"
FORMAT = 1
LD_STEP = 0.02                 # the LD factor is stored in one byte, in steps of 0.02 (0 .. 5.1)
# The stored LD factor assumes every neighbouring clump is represented.  When a genome file
# covers only a share of them, the part of the factor that comes from the neighbours shrinks;
# on 1000 Genomes it follows share ** 0.72 closely for chips of very different density
# (tools/check_bundle.py).
LD_SHARE_POWER = 0.72
FREQ_SCALE = 65535
_DTYPES = {"u1": np.uint8, "u2": np.uint16, "u4": np.uint32, "i2": np.int16, "i4": np.int32}
_N = NormalDist()


class BundleError(ValueError):
    pass


def write_bundle(path: Path, header: dict, arrays: dict[str, np.ndarray]) -> int:
    """Write a bundle; returns the compressed size in bytes."""
    codes = {np.dtype(v).str.lstrip("<|"): k for k, v in _DTYPES.items()}
    sections, blobs, offset = {}, [], 0
    for name, a in arrays.items():
        a = np.ascontiguousarray(a)
        code = codes[a.dtype.str.lstrip("<|")]
        raw = a.astype(a.dtype.newbyteorder("<"), copy=False).tobytes()
        sections[name] = {"type": code, "offset": offset, "length": int(a.size)}
        pad = (-len(raw)) % 8
        blobs.append(raw + b"\0" * pad)
        offset += len(raw) + pad
    head = json.dumps(dict(header, format=FORMAT, sections=sections), separators=(",", ":")).encode("utf-8")
    head += b" " * ((-(len(MAGIC) + 4 + len(head))) % 8)        # arrays start on an 8-byte boundary
    with gzip.open(path, "wb", compresslevel=9) as fh:
        fh.write(MAGIC + struct.pack("<I", len(head)) + head)
        for b in blobs:
            fh.write(b)
    return Path(path).stat().st_size


@dataclass
class BundleDisorder:
    meta: dict
    idx: np.ndarray            # row -> SNP index, clump by clump, most significant first
    beta: np.ndarray           # log-odds per A allele
    ld: np.ndarray             # LD factor
    clump_start: np.ndarray    # first row of each clump


@dataclass
class Bundle:
    header: dict
    rsid: np.ndarray           # uint32, in genome order
    b_allele: np.ndarray       # 1 = the SNP is A/C, 2 = A/G
    freq: np.ndarray           # frequency of the A allele among Europeans
    disorders: list[BundleDisorder]


def read_bundle(path: str | Path) -> Bundle:
    raw = Path(path).read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    if raw[:8] != MAGIC:
        raise BundleError(f"{Path(path).name}: not a pgc-risk reference bundle")
    (hlen,) = struct.unpack("<I", raw[8:12])
    header = json.loads(raw[12:12 + hlen])
    if header.get("format") != FORMAT:
        raise BundleError(f"{Path(path).name}: bundle format {header.get('format')} is not supported "
                          f"by this version (it reads format {FORMAT})")
    base = 12 + hlen

    def arr(name):
        s = header["sections"][name]
        dt = np.dtype(_DTYPES[s["type"]]).newbyteorder("<")
        return np.frombuffer(raw, dtype=dt, count=s["length"], offset=base + s["offset"])

    disorders = []
    for d in header["disorders"]:
        k = d["key"]
        lens = arr(f"{k}.clump_len").astype(np.int64)
        disorders.append(BundleDisorder(
            d, np.cumsum(arr(f"{k}.idx_delta").astype(np.int64)),
            arr(f"{k}.beta").astype(np.float64) * d["beta_scale"],
            arr(f"{k}.ld").astype(np.float64) * LD_STEP,
            np.concatenate([[0], np.cumsum(lens)[:-1]])))
    return Bundle(header, arr("snp.rsid"), arr("snp.b"), arr("snp.freq").astype(np.float64) / FREQ_SCALE,
                  disorders)


def liability_risk(z: float, prevalence: float, r2: float) -> float:
    if r2 <= 0:
        return prevalence
    t = _N.inv_cdf(1 - prevalence)
    return 1 - _N.cdf((t - (r2 ** 0.5) * z) / (1 - r2) ** 0.5)


def a_allele_counts(bundle: Bundle, rsid: np.ndarray, a1: np.ndarray, a2: np.ndarray) -> np.ndarray:
    """Count of the A allele for every bundle SNP (-1 where the genome has no usable call).

    `rsid`, `a1`, `a2` are the genome's records (alleles coded A=0, C=1, G=2, T=3, -1 = no
    call).  An rs number listed more than once in the genome is not used.
    """
    counts = np.full(len(bundle.rsid), -1, dtype=np.int8)
    ok = (rsid >= 0) & (rsid < 2 ** 32)
    rs, u1, u2 = rsid[ok], a1[ok].astype(np.int16), a2[ok].astype(np.int16)
    uniq, first, n = np.unique(rs, return_index=True, return_counts=True)
    once = first[n == 1]
    rs, u1, u2 = rs[once], u1[once], u2[once]
    order = np.argsort(bundle.rsid, kind="stable")
    sorted_rs = bundle.rsid[order].astype(np.int64)
    j = np.searchsorted(sorted_rs, rs)
    j[j == len(sorted_rs)] = 0
    hit = sorted_rs[j] == rs
    snp, u1, u2 = order[j[hit]], u1[hit], u2[hit]
    b = bundle.b_allele[snp].astype(np.int16)
    called = (u1 >= 0) & (u2 >= 0)
    direct = called & ((u1 == 0) | (u1 == b)) & ((u2 == 0) | (u2 == b))
    f1, f2 = 3 - u1, 3 - u2                                    # the same call read on the other strand
    flipped = called & ~direct & ((f1 == 0) | (f1 == b)) & ((f2 == 0) | (f2 == b))
    counts[snp[direct]] = ((u1 == 0).astype(np.int8) + (u2 == 0))[direct]
    counts[snp[flipped]] = ((f1 == 0).astype(np.int8) + (f2 == 0))[flipped]
    return counts


def score_counts(bundle: Bundle, counts: np.ndarray, r2_scale: float = 1.0) -> dict:
    """Score every disorder from A-allele counts.  The browser version computes exactly this."""
    h = bundle.header
    have = counts >= 0
    common = have & (bundle.freq >= 0.05) & (bundle.freq <= 0.95)
    f, d = bundle.freq[common], counts[common].astype(np.float64)
    fit = float(np.mean((d - 2 * f) ** 2 / (2 * f * (1 - f)))) if common.any() else float("nan")
    fit_snps = int(common.sum())
    valid, warnings = True, []
    if fit_snps < h["min_fit_snps"]:
        valid = False
        warnings.append(f"Only {fit_snps:,} of this file's SNPs are in the reference data - too few to score.")
    elif fit > h["fit_poor"]:
        valid = False
        warnings.append(
            f"Ancestry fit {fit:.3f}: this genome does not match the European-ancestry reference the "
            "scores are calibrated on (non-European or strongly admixed ancestry). In such genomes the "
            "scores sit several standard deviations off for reasons unrelated to risk, so no "
            "percentiles or probabilities are given.")
    elif fit > h["fit_good"]:
        warnings.append(
            f"Ancestry fit {fit:.3f}: this genome sits some distance from the European-ancestry "
            "reference (about 1 European in 20 does; so do partly non-European genomes). "
            "Treat the numbers as rough.")
    out = []
    big = np.iinfo(np.int64).max
    for dis in bundle.disorders:
        m = dis.meta
        rows = np.where(have[dis.idx], np.arange(len(dis.idx)), big)
        first = np.minimum.reduceat(rows, dis.clump_start)         # each clump's best available SNP
        rep = first[first < big]
        snp = dis.idx[rep]
        b, fr, dose = dis.beta[rep], bundle.freq[snp], counts[snp].astype(np.float64)
        v = 2 * fr * (1 - fr) * b * b
        raw, mean, indep = float((b * dose).sum()), float((2 * fr * b).sum()), float(v.sum())
        cross = float((v * (dis.ld[rep] - 1)).sum())
        coverage = min(1.0, indep / m["full_indep"]) if m["full_indep"] > 0 else 0.0
        var = indep + coverage ** h["ld_share_power"] * cross
        z = (raw - mean) / var ** 0.5 if var > 0 else 0.0
        r2 = m["r2"] * r2_scale
        prob = liability_risk(z, m["prevalence"], r2)
        res = {"key": m["key"], "name": m["name"], "snps_used": int(len(rep)), "clumps": int(len(dis.clump_start)),
               "coverage": coverage, "z": z, "percentile": _N.cdf(z), "prevalence": m["prevalence"],
               "r2_liability": r2, "r2_basis": m["r2_basis"], "probability": prob,
               "relative_risk": prob / m["prevalence"], "low_power": m["r2"] < h["low_power_r2"],
               "warnings": []}
        if coverage < 0.2:
            res["warnings"].append("this file has few of the SNPs the score is built from; the estimate is "
                                   "less accurate than for a typical genotyping chip")
        out.append(res)
    return {"valid": valid, "warnings": warnings, "ancestry_fit": {"value": fit, "snps": fit_snps,
            "good_below": h["fit_good"], "invalid_above": h["fit_poor"]}, "disorders": out}


def score_genome(bundle: Bundle, genome, r2_scale: float = 1.0) -> dict:
    """Score a parsed genome file (pgc_risk.genome.Genome) against a bundle."""
    if genome.fmt == "vcf" and (genome.variants_only or genome.ref_blocks):
        raise BundleError(
            "A sequencing VCF lists only the sites that differ from the reference, and the bundle has no "
            "reference genome to fill in the rest. Use the full build (pgc-risk fetch, pgc-risk build) for it.")
    res = score_counts(bundle, a_allele_counts(bundle, genome.rsid, genome.a1, genome.a2), r2_scale)
    res["genome"] = {"file": Path(genome.path).name, "format": genome.fmt, "build": genome.build,
                     "called": int(genome.called.sum()), "notes": genome.notes}
    return res


def render_text(res: dict, bundle: Bundle) -> str:
    """The result of score_genome as the table the command line prints."""
    def pct(x: float) -> str:
        v = 100 * x
        return f"{v:.2f}%" if v < 0.095 else (f"{v:.1f}%" if v < 9.95 else f"{int(v + 0.5)}%")

    def ordinal(p: float) -> str:
        v = min(99, max(1, round(100 * p)))
        return f"{v}{'th' if 10 <= v % 100 <= 20 else {1: 'st', 2: 'nd', 3: 'rd'}.get(v % 10, 'th')}"

    g, fit = res["genome"], res["ancestry_fit"]
    out = [f"pgc-risk, scored from {bundle.header['name']} (built {bundle.header['created']})", "",
           f"Genome file  : {g['file']}", f"Genotypes    : {g['called']:,} autosomal SNPs with a call",
           f"Ancestry fit : {fit['value']:.3f} over {fit['snps']:,} SNPs (1.00 = typical, above "
           f"{fit['invalid_above']:g} = not scorable)"]
    for w in res["warnings"]:
        out += ["", f"WARNING: {w}"]
    head = f"{'Disorder':<34}{'Population':>11}{'Your':>10}{'Relative':>10}{'Score':>12}{'SNPs':>8}"
    out += ["", head, f"{'':<34}{'baseline':>11}{'estimate':>10}{'to avg.':>10}{'percentile':>12}{'used':>8}",
            "-" * len(head)]
    for d in res["disorders"]:
        name = d["name"] + (" *" if d["low_power"] else "")
        if res["valid"]:
            out.append(f"{name:<34}{pct(d['prevalence']):>11}{pct(d['probability']):>10}{d['relative_risk']:>9.2f}x"
                       f"{ordinal(d['percentile']):>12}{d['snps_used']:>8,}")
        else:
            out.append(f"{name:<34}{pct(d['prevalence']):>11}{'-':>10}{'-':>10}{'-':>12}{d['snps_used']:>8,}")
        out += [f"    ! {w}" for w in d["warnings"]]
    out += ["", "* the GWAS behind this score is small: the estimate stays close to the baseline.",
            "A statistical expectation for people with a similar score, not a diagnosis. European ancestry only."]
    return "\n".join(out)
