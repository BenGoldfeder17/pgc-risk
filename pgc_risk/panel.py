"""European linkage-disequilibrium reference panel, built from 1000 Genomes phase 3.

The scorer needs to know which SNPs are correlated: to pick one SNP per association signal
(LD clumping) and to get the spread of a score right when neighbouring SNPs are not fully
independent.  This module downloads the 1000 Genomes phase-3 genotypes, keeps the SNPs that
occur in any weight file, and stores the 503 European-ancestry samples in a compact form
(two bit-planes per SNP, about 126 bytes).

    download()      VCFs + sample list from the 1000 Genomes public bucket (about 15 GB)
    extract_all()   genotypes of all 2,504 samples at the weight SNPs (cached per chromosome)
    write_panel()   the compact European panel in data/panel/
    annotate()      adds each weight SNP's panel row to the weight files
"""
from __future__ import annotations

import gzip
import json
import shutil
import time
import urllib.request
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from pathlib import Path

import numpy as np

from .fetch import DATA

KG_DIR = DATA / "reference" / "1000g"
PANEL_DIR = DATA / "panel"
WEIGHTS = DATA / "weights"
_S3 = "https://1000genomes.s3.amazonaws.com/release/20130502/"
VCF_URL = _S3 + "ALL.chr{c}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
SAMPLES_URL = _S3 + "integrated_call_samples_v3.20130502.ALL.panel"
SAMPLES_FILE = "integrated_call_samples_v3.20130502.ALL.panel"
_CODE = {b"A": 0, b"C": 1, b"G": 2, b"T": 3}


def vcf_path(chrom: int) -> Path:
    return KG_DIR / f"chr{chrom}.phase3_v5a.vcf.gz"


def _download_one(url: str, dest: Path) -> str:
    if dest.exists() and dest.stat().st_size > 0:
        return "cached"
    part = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url) as resp, open(part, "wb") as out:
        shutil.copyfileobj(resp, out, 1 << 20)
    part.replace(dest)
    return "downloaded"


def download(chroms=range(1, 23), workers: int = 6, log=print) -> None:
    KG_DIR.mkdir(parents=True, exist_ok=True)
    jobs = [(SAMPLES_URL, KG_DIR / SAMPLES_FILE)] + [(VCF_URL.format(c=c), vcf_path(c)) for c in chroms]
    t0 = time.time()
    with ThreadPoolExecutor(workers) as pool:
        states = list(pool.map(lambda j: _download_one(*j), jobs))
    log(f"[panel] 1000 Genomes: {states.count('downloaded')} file(s) downloaded, "
        f"{states.count('cached')} already present ({time.time() - t0:.0f}s)")


def samples() -> dict[str, tuple[str, str]]:
    """Sample -> (population, super-population)."""
    info = {}
    with open(KG_DIR / SAMPLES_FILE) as fh:
        next(fh)
        for line in fh:
            f = line.split()
            if len(f) >= 3:
                info[f[0]] = (f[1], f[2])
    return info


def wanted_positions(chrom: int, weight_dirs: list[Path]) -> np.ndarray:
    want = []
    for d in weight_dirs:
        for f in sorted(Path(d).glob("*.npz")):
            z = np.load(f)
            want.append(z["pos"][z["chrom"] == chrom])
    return np.unique(np.concatenate(want)) if want else np.empty(0, dtype=np.int32)


def _cache(chrom: int) -> Path:
    return KG_DIR / f"chr{chrom}.genotypes.npz"


def extract(chrom: int, weight_dirs: list[Path] | None = None, force: bool = False) -> dict:
    """ALT-allele counts of every 1000 Genomes sample at the weight SNPs of one chromosome."""
    weight_dirs = weight_dirs or [WEIGHTS]
    want_arr = wanted_positions(chrom, weight_dirs)
    cache = _cache(chrom)
    if cache.exists() and not force:
        z = np.load(cache)
        if np.isin(want_arr, z["wanted"]).all():   # the cache covers at least what is asked for
            return {k: z[k] for k in z.files}
    want = set(want_arr.tolist())
    pos, ref, alt, rows = [], [], [], []
    names: list[str] = []
    n = 0
    with gzip.open(vcf_path(chrom), "rb") as fh:
        for line in fh:
            if line[:1] == b"#":
                if line.startswith(b"#CHROM"):
                    names = line.decode().rstrip("\r\n").split("\t")[9:]
                    n = len(names)
                continue
            t1 = line.find(b"\t")
            t2 = line.find(b"\t", t1 + 1)
            p = int(line[t1 + 1:t2])
            if p not in want:
                continue
            f = line.split(b"\t", 9)
            r, a = _CODE.get(f[3]), _CODE.get(f[4])
            if r is None or a is None:
                continue                          # indel or multi-allelic record
            buf = np.frombuffer(f[9], dtype=np.uint8)
            if len(buf) != 4 * n:
                continue                          # not plain phased diploid calls
            g = (buf[0::4] - 48) + (buf[2::4] - 48)
            if g.max() > 2:
                continue
            pos.append(p)
            ref.append(r)
            alt.append(a)
            rows.append(g.astype(np.int8))
    pos_a = np.array(pos, dtype=np.int32)
    uniq, counts = np.unique(pos_a, return_counts=True)
    keep = np.isin(pos_a, uniq[counts == 1])      # a position listed twice is multi-allelic
    geno = np.array(rows, dtype=np.int8)[keep] if rows else np.empty((0, n), dtype=np.int8)
    out = {"pos": pos_a[keep], "ref": np.array(ref, dtype=np.int8)[keep],
           "alt": np.array(alt, dtype=np.int8)[keep], "geno": geno,
           "samples": np.array(names), "wanted": want_arr}
    order = np.argsort(out["pos"], kind="stable")
    for k in ("pos", "ref", "alt", "geno"):
        out[k] = out[k][order]
    np.savez(cache, **out)
    return out


def _extract_job(args) -> tuple[int, int, int, float]:
    chrom, dirs = args
    t0 = time.time()
    d = extract(chrom, dirs)
    return chrom, len(d["pos"]), len(d["wanted"]), time.time() - t0


def extract_all(weight_dirs: list[Path] | None = None, jobs: int = 8, log=print) -> None:
    dirs = [Path(d) for d in (weight_dirs or [WEIGHTS])]
    with ProcessPoolExecutor(jobs) as pool:
        for chrom, found, wanted, secs in pool.map(_extract_job, [(c, dirs) for c in range(1, 23)]):
            log(f"[panel] chr{chrom}: {found:,} of {wanted:,} weight SNPs in 1000 Genomes ({secs:.0f}s)")


def write_panel(log=print) -> dict:
    """Compact European panel: per chromosome pos / ref / alt and two packed bit-planes."""
    PANEL_DIR.mkdir(parents=True, exist_ok=True)
    info = samples()
    total = 0
    n_eur = 0
    for c in range(1, 23):
        z = np.load(_cache(c))
        eur = np.array([info[s][1] == "EUR" for s in z["samples"].tolist()])
        g = z["geno"][:, eur]
        n_eur = int(eur.sum())
        poly = g.std(axis=1) > 0                  # monomorphic in Europeans: no LD information
        g = g[poly]
        np.savez(PANEL_DIR / f"chr{c}.npz", pos=z["pos"][poly], ref=z["ref"][poly], alt=z["alt"][poly],
                 p1=np.packbits(g >= 1, axis=1), p2=np.packbits(g == 2, axis=1))
        total += int(poly.sum())
    meta = {"source": "1000 Genomes phase 3 (release 20130502, v5a), European super-population",
            "samples": n_eur, "snps": total}
    (PANEL_DIR / "panel.json").write_text(json.dumps(meta, indent=1))
    log(f"[panel] wrote {total:,} SNPs x {n_eur} European samples to {PANEL_DIR}")
    return meta


class Panel:
    def __init__(self, directory: Path = PANEL_DIR):
        self.dir = Path(directory)
        self.meta = json.loads((self.dir / "panel.json").read_text())
        self.n = int(self.meta["samples"])
        self._chrom: dict[int, dict] = {}

    @staticmethod
    def exists(directory: Path = PANEL_DIR) -> bool:
        return (Path(directory) / "panel.json").exists()

    def chrom(self, c: int) -> dict:
        if c not in self._chrom:
            z = np.load(self.dir / f"chr{c}.npz")
            self._chrom[c] = {k: z[k] for k in z.files}
        return self._chrom[c]

    def lookup(self, c: int, pos: np.ndarray, a1: np.ndarray, a2: np.ndarray):
        """Panel row of each SNP (-1 if absent or alleles disagree) and whether a1 is the ALT allele."""
        d = self.chrom(c)
        row = np.full(len(pos), -1, dtype=np.int32)
        is_alt = np.zeros(len(pos), dtype=bool)
        if not len(d["pos"]) or not len(pos):
            return row, is_alt
        j = np.searchsorted(d["pos"], pos)
        j[j == len(d["pos"])] = 0
        hit = d["pos"][j] == pos
        ref, alt = d["ref"][j], d["alt"][j]
        alt_like = ((a1 == alt) & (a2 == ref)) | ((3 - a1 == alt) & (3 - a2 == ref))
        ref_like = ((a1 == ref) & (a2 == alt)) | ((3 - a1 == ref) & (3 - a2 == alt))
        ok = hit & (alt_like | ref_like)
        row[ok] = j[ok]
        is_alt[ok] = alt_like[ok]
        return row, is_alt

    def dosage(self, c: int, rows: np.ndarray) -> np.ndarray:
        """ALT-allele counts, shape (len(rows), n samples), int8."""
        d = self.chrom(c)
        g = np.unpackbits(d["p1"][rows], axis=1, count=self.n)
        g += np.unpackbits(d["p2"][rows], axis=1, count=self.n)
        return g.view(np.int8)


def annotate(weight_dir: Path = WEIGHTS, panel: Panel | None = None, log=print) -> None:
    """Add `prow` (panel row) and `palt` (effect allele is the panel's ALT) to each weight file."""
    panel = panel or Panel()
    for f in sorted(Path(weight_dir).glob("*.npz")):
        z = np.load(f)
        arrays = {k: z[k] for k in z.files}
        n = len(arrays["rsid"])
        prow = np.full(n, -1, dtype=np.int32)
        palt = np.zeros(n, dtype=bool)
        for c in range(1, 23):
            m = np.flatnonzero(arrays["chrom"] == c)
            if len(m):
                prow[m], palt[m] = panel.lookup(c, arrays["pos"][m], arrays["a1"][m], arrays["a2"][m])
        arrays["prow"], arrays["palt"] = prow, palt
        np.savez(f, **arrays)
        log(f"[panel] {f.stem}: {int((prow >= 0).sum()):,} of {n:,} SNPs are in the panel")


def build_panel(weight_dirs: list[Path] | None = None, jobs: int = 8, log=print) -> None:
    dirs = [Path(d) for d in (weight_dirs or [WEIGHTS])]
    download(log=log)
    extract_all(dirs, jobs, log)
    write_panel(log)
    panel = Panel()
    for d in dirs:
        annotate(d, panel, log)
