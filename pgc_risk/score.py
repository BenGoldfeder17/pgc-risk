"""Match a genome to the GWAS weights, build a polygenic score and convert it to a risk.

Method: clumping + thresholding (the PGC's own polygenic-score recipe).

* Candidate SNPs are the GWAS SNPs with p <= p_threshold that the genome file actually
  has a usable genotype for.
* LD clumping against the European reference panel: walking the candidates from the most to
  the least significant, a SNP is kept and every candidate within `window` base pairs that is
  correlated with it (r^2 >= r2) is dropped.  Three regions with megabase-scale linkage
  disequilibrium contribute a single SNP each, and the other long-range-LD regions use a
  four times wider window.
* score = sum(log-odds * effect-allele count).  It is standardised against the score
  distribution expected in a European-ancestry population: mean = sum 2pb from the GWAS
  allele frequencies, variance = sum 2p(1-p)b^2 times an LD factor measured on the panel
  (clumped SNPs are still weakly correlated, which widens the real distribution).
* The standardised score z becomes an absolute risk through the liability-threshold model:
  with lifetime prevalence K and a score that explains a share R2 of liability variance,
  risk(z) = 1 - Phi((T - sqrt(R2) z) / sqrt(1 - R2)),  T = Phi^-1(1 - K).

Without a panel the code falls back to clumping by distance alone (used by the unit tests).
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from statistics import NormalDist

import numpy as np

from .genome import Genome

_N = NormalDist()

DEFAULT_P = 0.05
DEFAULT_WINDOW = 250_000      # clumping reach on each side of an index SNP
DEFAULT_R2 = 0.1
LD_BAND = 2_000_000           # SNP pairs further apart than this are treated as uncorrelated
GRID = 1000

# hg19 regions where linkage disequilibrium runs over megabases.
ONE_SNP_REGIONS = [            # (chromosome, start, end): one SNP for the whole region
    (6, 25_000_000, 34_000_000),       # MHC
    (8, 7_000_000, 13_000_000),        # 8p23.1 inversion
    (17, 43_000_000, 45_500_000),      # 17q21.31 inversion
]
WIDE_WINDOW_REGIONS = [        # long-range LD regions (Price et al. 2008): 4x window
    (1, 48_000_000, 52_000_000), (2, 86_000_000, 100_500_000), (2, 134_500_000, 138_000_000),
    (2, 183_000_000, 190_000_000), (3, 47_500_000, 50_000_000), (3, 83_500_000, 87_000_000),
    (3, 89_000_000, 97_500_000), (5, 44_500_000, 50_500_000), (5, 98_000_000, 100_500_000),
    (5, 129_000_000, 132_000_000), (5, 135_500_000, 138_500_000), (6, 57_000_000, 64_000_000),
    (6, 140_000_000, 142_500_000), (7, 55_000_000, 66_000_000), (8, 43_000_000, 50_000_000),
    (8, 112_000_000, 115_000_000), (10, 37_000_000, 43_000_000), (11, 46_000_000, 57_000_000),
    (11, 87_500_000, 90_500_000), (12, 33_000_000, 40_000_000), (12, 109_500_000, 112_000_000),
    (20, 32_000_000, 34_500_000),
]


@dataclass
class Weights:
    key: str
    meta: dict
    rsid: np.ndarray
    chrom: np.ndarray
    pos: np.ndarray
    pos38: np.ndarray
    a1: np.ndarray
    a2: np.ndarray
    ref: np.ndarray
    beta: np.ndarray
    p: np.ndarray
    freq: np.ndarray
    prow: np.ndarray | None = None      # row in the LD panel (-1 = not in the panel)
    palt: np.ndarray | None = None      # True where the effect allele is the panel's ALT allele

    @property
    def n(self) -> int:
        return len(self.rsid)


def load_weights(key: str, directory: Path) -> Weights:
    z = np.load(directory / f"{key}.npz")
    meta = json.loads((directory / f"{key}.json").read_text())
    return Weights(key, meta, z["rsid"], z["chrom"], z["pos"], z["pos38"], z["a1"], z["a2"], z["ref"],
                   z["beta"].astype(np.float64), z["p"], z["freq"].astype(np.float64),
                   z["prow"] if "prow" in z.files else None, z["palt"] if "palt" in z.files else None)


def _key(chrom: np.ndarray, pos: np.ndarray) -> np.ndarray:
    return (chrom.astype(np.int64) << 32) | pos.astype(np.int64)


class GenomeIndex:
    """Lookup tables over a genome file, shared by every disorder."""

    def __init__(self, g: Genome):
        self.g = g
        has_rs = np.flatnonzero(g.rsid >= 0)
        order = has_rs[np.argsort(g.rsid[has_rs], kind="stable")]
        rs = g.rsid[order]
        # an rs number listed twice with different calls is ambiguous: drop it
        dup = np.zeros(len(rs), dtype=bool)
        if len(rs) > 1:
            same = rs[1:] == rs[:-1]
            dup[1:] |= same
            dup[:-1] |= same
        self.rs_sorted, self.rs_index = rs[~dup], order[~dup]
        k = _key(g.chrom, g.pos)
        order = np.argsort(k, kind="stable")
        ks = k[order]
        dup = np.zeros(len(ks), dtype=bool)
        if len(ks) > 1:
            same = ks[1:] == ks[:-1]
            dup[1:] |= same
            dup[:-1] |= same
        self.key_sorted, self.key_index = ks[~dup], order[~dup]
        self.key_all = ks                       # with duplicates: "is anything at this position?"
        self.pos_build = g.build if g.build in (37, 38) else None
        self._build_checked = g.fmt != "vcf" or g.build in (37, 38)

    @staticmethod
    def _lookup(sorted_keys: np.ndarray, index: np.ndarray, query: np.ndarray) -> np.ndarray:
        if not len(sorted_keys):
            return np.full(len(query), -1, dtype=np.int64)
        j = np.searchsorted(sorted_keys, query)
        j[j == len(sorted_keys)] = 0
        return np.where(sorted_keys[j] == query, index[j], -1)

    def _infer_build(self, w: Weights) -> None:
        """A VCF that does not name its assembly: see which coordinates its sites agree with."""
        self._build_checked = True
        h37 = np.isin(_key(w.chrom, w.pos), self.key_all).mean()
        ok38 = w.pos38 > 0
        h38 = np.isin(_key(w.chrom[ok38], w.pos38[ok38]), self.key_all).mean() if ok38.any() else 0.0
        best, other = (37, h38) if h37 >= h38 else (38, h37)
        if max(h37, h38) > 0.02 and max(h37, h38) > 5 * other:
            self.pos_build = best
            self.g.notes.append(f"Genome build not stated in the file; positions match GRCh{best}.")

    def match(self, w: Weights) -> tuple[np.ndarray, np.ndarray, dict]:
        """For every weight SNP: is a genotype available, and the effect-allele count."""
        g = self.g
        if not self._build_checked:
            self._infer_build(w)
        gi = self._lookup(self.rs_sorted, self.rs_index, w.rsid)
        by_rs = gi >= 0
        wpos = None
        if self.pos_build == 37:
            wpos = w.pos
        elif self.pos_build == 38:
            wpos = w.pos38
        by_pos = np.zeros(w.n, dtype=bool)
        if wpos is not None:
            todo = ~by_rs & (wpos > 0)
            hit = self._lookup(self.key_sorted, self.key_index, _key(w.chrom[todo], wpos[todo]))
            idx = np.flatnonzero(todo)
            gi[idx] = hit
            by_pos[idx] = hit >= 0
        found = gi >= 0
        u1 = np.where(found, g.a1[np.clip(gi, 0, None)], -1).astype(np.int16)
        u2 = np.where(found, g.a2[np.clip(gi, 0, None)], -1).astype(np.int16)
        called = (u1 >= 0) & (u2 >= 0)
        a1, a2 = w.a1.astype(np.int16), w.a2.astype(np.int16)
        direct = called & ((u1 == a1) | (u1 == a2)) & ((u2 == a1) | (u2 == a2))
        f1, f2 = 3 - u1, 3 - u2
        flipped = called & ~direct & ((f1 == a1) | (f1 == a2)) & ((f2 == a1) | (f2 == a2))
        dosage = np.zeros(w.n, dtype=np.int8)
        dosage[direct] = ((u1 == a1).astype(np.int8) + (u2 == a1))[direct]
        dosage[flipped] = ((f1 == a1).astype(np.int8) + (f2 == a1))[flipped]
        avail = direct | flipped
        info = {"by_rsid": int((by_rs & avail).sum()), "by_position": int((by_pos & avail).sum()),
                "no_call": int((found & ~called).sum()), "allele_mismatch": int((called & ~avail).sum()),
                "strand_flipped": int(flipped.sum()), "assumed_reference": 0}
        # Sites a sequencing file does not list are homozygous for the reference base.
        if (g.variants_only or g.ref_blocks) and wpos is not None:
            absent = ~found & (wpos > 0) & (w.ref >= 0)
            if g.variants_only:
                # ...unless the file has some other record there (an indel, a filtered call)
                absent &= ~np.isin(_key(w.chrom, wpos), self.key_all)
            else:
                covered = np.zeros(w.n, dtype=bool)
                for c, (starts, ends) in g.ref_blocks.items():
                    m = np.flatnonzero(absent & (w.chrom == c))
                    if len(m):
                        k = np.searchsorted(starts, wpos[m], side="right") - 1
                        covered[m] = (k >= 0) & (wpos[m] <= ends[np.clip(k, 0, None)])
                absent &= covered
            ref = w.ref.astype(np.int16)
            fwd = absent & ((ref == a1) | (ref == a2))
            rev = absent & ~fwd & ((3 - ref == a1) | (3 - ref == a2))
            dosage[fwd] = 2 * (ref == a1)[fwd]
            dosage[rev] = 2 * (3 - ref == a1)[rev]
            avail |= fwd | rev
            info["assumed_reference"] = int((fwd | rev).sum())
        info["available"] = int(avail.sum())
        return avail, dosage, info


def clump(w: Weights, candidates: np.ndarray, window: int = DEFAULT_WINDOW) -> np.ndarray:
    """Greedy distance clumping; returns the indices of the SNPs kept (most significant first)."""
    if not len(candidates):
        return candidates
    order = candidates[np.argsort(w.p[candidates], kind="stable")]
    chroms = w.chrom[order].tolist()
    cells = (w.pos[order] // GRID).tolist()
    half = max(1, window // GRID)
    blocked: dict[int, bytearray] = {}
    special: dict[int, list[tuple[int, int, int]]] = {}
    for c, s, e in ONE_SNP_REGIONS:
        special.setdefault(c, []).append((s // GRID, e // GRID, 0))
    for c, s, e in WIDE_WINDOW_REGIONS:
        special.setdefault(c, []).append((s // GRID, e // GRID, 4 * half))
    keep = []
    for i, c, cell in zip(order.tolist(), chroms, cells):
        b = blocked.get(c)
        if b is None:
            b = blocked[c] = bytearray(260_000)
        if b[cell]:
            continue
        keep.append(i)
        lo, hi = cell - half, cell + half + 1
        for s, e, wide in special.get(c, ()):
            if s <= cell <= e:
                if wide:
                    lo, hi = cell - wide, cell + wide + 1
                else:
                    lo, hi = min(lo, s), max(hi, e + 1)
                break
        lo = max(lo, 0)
        hi = min(hi, len(b))
        b[lo:hi] = b"\x01" * (hi - lo)
    return np.array(keep, dtype=np.int64)


def _special_span(c: int, pos: int, window: int) -> tuple[int, int, bool]:
    """Base-pair span an index SNP clears, and whether it stands for a whole one-SNP region."""
    for cc, s, e in ONE_SNP_REGIONS:
        if cc == c and s <= pos <= e:
            return min(pos - window, s), max(pos + window, e), True
    for cc, s, e in WIDE_WINDOW_REGIONS:
        if cc == c and s <= pos <= e:
            return pos - 4 * window, pos + 4 * window, False
    return pos - window, pos + window, False


def ld_clump(w: Weights, candidates: np.ndarray, panel, window: int = DEFAULT_WINDOW,
             r2: float = DEFAULT_R2, owners: np.ndarray | None = None) -> np.ndarray:
    """LD clumping against the reference panel; returns the indices of the index SNPs.

    Every candidate must be in the panel (w.prow >= 0).  Weights are sorted by chromosome and
    position, so each chromosome's candidates are already in position order.  If `owners`
    (an int array as long as the weights) is given, it receives for every candidate the index
    SNP whose clump it ended up in (an index SNP owns itself).
    """
    keep_all = []
    special_chroms = {c for c, _, _ in ONE_SNP_REGIONS} | {c for c, _, _ in WIDE_WINDOW_REGIONS}
    cand_chrom = w.chrom[candidates]
    for c in np.unique(cand_chrom).tolist():
        idx = np.sort(candidates[cand_chrom == c])
        pos = w.pos[idx].astype(np.int64)
        g = panel.dosage(c, w.prow[idx]).astype(np.float32)
        g -= g.mean(axis=1, keepdims=True)
        norm = np.sqrt((g * g).sum(axis=1))
        g /= np.where(norm > 0, norm, 1)[:, None]
        lo = np.searchsorted(pos, pos - window).tolist()
        hi = np.searchsorted(pos, pos + window, side="right").tolist()
        alive = np.ones(len(idx), dtype=bool)
        order = np.argsort(w.p[idx], kind="stable").tolist()
        special = c in special_chroms
        pos_l = pos.tolist()
        kept = []
        for i in order:
            if not alive[i]:
                continue
            kept.append(i)
            l, h = lo[i], hi[i]
            if special:
                s_bp, e_bp, whole = _special_span(c, pos_l[i], window)
                if whole or e_bp - s_bp != 2 * window:
                    l = int(np.searchsorted(pos, s_bp))
                    h = int(np.searchsorted(pos, e_bp, side="right"))
                if whole:
                    if owners is not None:
                        owners[idx[l:h][alive[l:h]]] = idx[i]
                    alive[l:h] = False             # the region is represented by this one SNP
                    continue
            r = g[l:h] @ g[i]
            if owners is not None:
                owners[idx[l:h][alive[l:h] & ((r * r) >= r2)]] = idx[i]
                owners[idx[i]] = idx[i]
            alive[l:h] &= (r * r) < r2
            alive[i] = False
        keep_all.append(idx[np.array(kept, dtype=np.int64)])
    return np.concatenate(keep_all) if keep_all else candidates[:0]


def ld_variance_terms(w: Weights, selected: np.ndarray, panel, band: int = LD_BAND) -> tuple[float, float]:
    """Panel variance of the score counting SNP pairs within `band` bp, and counting none."""
    total = indep = 0.0
    sel_chrom = w.chrom[selected]
    for c in np.unique(sel_chrom).tolist():
        idx = np.sort(selected[sel_chrom == c])
        pos = w.pos[idx].astype(np.int64)
        g = panel.dosage(c, w.prow[idx]).astype(np.float64)
        g = np.where(w.palt[idx][:, None], g, 2 - g)
        y = (g - g.mean(axis=1, keepdims=True)) * w.beta[idx][:, None]
        csum = np.vstack([np.zeros((1, y.shape[1])), np.cumsum(y, axis=0)])
        lo = np.searchsorted(pos, pos - band)
        hi = np.searchsorted(pos, pos + band, side="right")
        total += float((y * (csum[hi] - csum[lo])).sum())
        indep += float((y * y).sum())
    return total, indep


def ld_inflation(w: Weights, selected: np.ndarray, panel, band: int = LD_BAND) -> float:
    """Variance of the score in the panel relative to the variance if SNPs were independent.

    Uses only SNP pairs within `band` base pairs, so the estimate is not swamped by the chance
    correlations between unlinked SNPs that any 503-person panel contains.
    """
    total, indep = ld_variance_terms(w, selected, panel, band)
    if indep <= 0:
        return 1.0
    return max(total / indep, 0.25)


def liability_risk(z: float, prevalence: float, r2: float) -> float:
    """Lifetime risk for a standardised polygenic score under the liability-threshold model."""
    if r2 <= 0:
        return prevalence
    t = _N.inv_cdf(1 - prevalence)
    return 1 - _N.cdf((t - (r2 ** 0.5) * z) / (1 - r2) ** 0.5)


@dataclass
class DisorderResult:
    key: str
    name: str
    snps_used: int
    snps_available: int
    snps_in_gwas: int
    raw_score: float
    expected_mean: float
    expected_sd: float
    ld_inflation: float
    z: float
    percentile: float
    prevalence: float | None = None
    r2_liability: float | None = None
    r2_basis: str = ""
    probability: float | None = None
    relative_risk: float | None = None
    low_power: bool = False
    match: dict = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


def score_disorder(index: GenomeIndex, w: Weights, panel=None, p_threshold: float = DEFAULT_P,
                   window: int = DEFAULT_WINDOW, r2: float = DEFAULT_R2,
                   matched: tuple | None = None) -> tuple[DisorderResult, dict]:
    avail, dosage, info = matched if matched is not None else index.match(w)
    usable = avail & (w.p <= p_threshold)
    if panel is not None:
        if w.prow is None:
            raise ValueError(f"{w.key}: weight file has no LD-panel annotation - run pgc-risk panel")
        sel = ld_clump(w, np.flatnonzero(usable & (w.prow >= 0)), panel, window, r2)
        inflation = ld_inflation(w, sel, panel) if len(sel) else 1.0
    else:
        sel = clump(w, np.flatnonzero(usable), window)
        inflation = 1.0
    b, f, d = w.beta[sel], w.freq[sel], dosage[sel].astype(np.float64)
    raw = float((b * d).sum())
    mean = float((2 * f * b).sum())
    var = float((2 * f * (1 - f) * b * b).sum()) * inflation
    sd = var ** 0.5
    z = (raw - mean) / sd if sd > 0 else 0.0
    res = DisorderResult(w.key, w.meta.get("name", w.key), int(len(sel)), int(avail.sum()), w.n,
                         raw, mean, sd, inflation, z, _N.cdf(z), match=info)
    if len(sel) < 200:
        res.warnings.append(f"only {len(sel)} usable SNPs - too few for a meaningful score")
    detail = {"avail": avail, "dosage": dosage, "selected": sel}
    return res, detail


def ancestry_fit(parts: list[tuple[np.ndarray, np.ndarray, np.ndarray]]) -> tuple[float, int]:
    """How well the genotypes fit the European reference frequencies.

    Returns mean((count - 2p)^2 / (2p(1-p))) over distinct common SNPs.  It is 1.0 for a person
    drawn from the reference population and grows with genetic distance from it.
    """
    if not parts:
        return float("nan"), 0
    rs = np.concatenate([p[0] for p in parts])
    d = np.concatenate([p[1] for p in parts]).astype(np.float64)
    f = np.concatenate([p[2] for p in parts])
    _, first = np.unique(rs, return_index=True)
    d, f = d[first], f[first]
    common = (f >= 0.05) & (f <= 0.95)
    d, f = d[common], f[common]
    if not len(d):
        return float("nan"), 0
    return float(np.mean((d - 2 * f) ** 2 / (2 * f * (1 - f)))), int(len(d))
