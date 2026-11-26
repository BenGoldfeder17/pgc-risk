"""Reference-genome helpers: read bases from a UCSC .2bit file and lift hg19 positions to hg38.

Both are needed only for sequencing VCFs.  A variants-only VCF omits every site where the
person matches the reference, so to score such a site we must know which of the two GWAS
alleles is the reference base (hg19.2bit) and, for hg38 files, where the SNP sits in hg38
(hg19ToHg38 chain).
"""
from __future__ import annotations

import gzip
import struct
from pathlib import Path

import numpy as np

UCSC_2BIT = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit"
UCSC_CHAIN = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"

def ensure(directory: Path, log=print) -> tuple[Path, Path]:
    """Download hg19.2bit (816 MB) and the hg19->hg38 chain from UCSC if they are not there yet."""
    import shutil
    import urllib.request

    directory.mkdir(parents=True, exist_ok=True)
    paths = []
    for url in (UCSC_2BIT, UCSC_CHAIN):
        dest = directory / url.rsplit("/", 1)[1]
        if not dest.exists() or dest.stat().st_size == 0:
            log(f"[build] downloading {url}")
            part = dest.with_suffix(dest.suffix + ".part")
            with urllib.request.urlopen(url) as resp, open(part, "wb") as out:
                shutil.copyfileobj(resp, out, 1 << 20)
            part.replace(dest)
        paths.append(dest)
    return paths[0], paths[1]


# .2bit packs T=0, C=1, A=2, G=3; this maps to our A=0, C=1, G=2, T=3.
_TWOBIT_TO_CODE = np.array([3, 1, 0, 2], dtype=np.int8)


class TwoBit:
    def __init__(self, path: str | Path):
        self.mm = np.memmap(path, dtype=np.uint8, mode="r")
        sig, version, count, _ = struct.unpack("<IIII", self.mm[:16].tobytes())
        if sig != 0x1A412743 or version != 0:
            raise ValueError(f"{path}: not a little-endian .2bit file")
        self.index: dict[str, int] = {}
        off = 16
        for _ in range(count):
            n = int(self.mm[off])
            name = self.mm[off + 1: off + 1 + n].tobytes().decode("ascii")
            (self.index[name],) = struct.unpack("<I", self.mm[off + 1 + n: off + 5 + n].tobytes())
            off += 5 + n
        self._seq: dict[str, tuple[int, np.ndarray, np.ndarray, int]] = {}

    def _u32(self, off: int, count: int = 1) -> np.ndarray:
        return np.frombuffer(self.mm[off: off + 4 * count].tobytes(), dtype="<u4")

    def _header(self, name: str):
        if name not in self._seq:
            o = self.index[name]
            size, n_blocks = (int(x) for x in self._u32(o, 2))
            starts = self._u32(o + 8, n_blocks).astype(np.int64)
            sizes = self._u32(o + 8 + 4 * n_blocks, n_blocks).astype(np.int64)
            o2 = o + 8 + 8 * n_blocks
            n_mask = int(self._u32(o2)[0])
            packed = o2 + 4 + 8 * n_mask + 4
            self._seq[name] = (size, starts, starts + sizes, packed)
        return self._seq[name]

    def bases(self, chrom: int, pos1: np.ndarray) -> np.ndarray:
        """Allele codes at 1-based positions of one chromosome (-1 for N or out of range)."""
        size, n_start, n_end, packed = self._header(f"chr{chrom}")
        p = np.asarray(pos1, dtype=np.int64) - 1
        ok = (p >= 0) & (p < size)
        q = np.where(ok, p, 0)
        byte = np.asarray(self.mm[packed + q // 4]).astype(np.uint8)
        two = (byte >> (6 - 2 * (q % 4)).astype(np.uint8)) & 3
        out = _TWOBIT_TO_CODE[two]
        if len(n_start):
            k = np.searchsorted(n_start, q, side="right") - 1
            in_n = (k >= 0) & (q < n_end[np.clip(k, 0, None)])
            out = np.where(in_n, -1, out)
        return np.where(ok, out, -1).astype(np.int8)


class LiftOver:
    """Point liftover along a UCSC chain file (same chromosome, forward strand only).

    Positions that fall in a gap, move to another chromosome or land on the reverse strand
    come back as -1; that is a fraction of a percent of common SNPs between hg19 and hg38.
    """

    def __init__(self, chain_path: str | Path):
        # per source chromosome: list of (score, block_start, block_end, target_start) arrays
        self.chains: dict[str, list[tuple[int, np.ndarray, np.ndarray, np.ndarray]]] = {}
        opener = gzip.open if str(chain_path).endswith(".gz") else open
        with opener(chain_path, "rt") as fh:
            cur = None
            for line in fh:
                if line.startswith("chain"):
                    f = line.split()
                    score, t_name, t_strand, t_start = int(f[1]), f[2], f[4], int(f[5])
                    q_name, q_strand, q_start = f[7], f[9], int(f[10])
                    keep = t_name == q_name and t_strand == "+" and q_strand == "+"
                    cur = [score, t_name, t_start, q_start, [], [], [], keep]
                elif cur is not None:
                    f = line.split()
                    if not f:
                        self._close(cur)
                        cur = None
                        continue
                    size = int(f[0])
                    cur[4].append(cur[2])
                    cur[5].append(cur[2] + size)
                    cur[6].append(cur[3])
                    if len(f) == 3:
                        cur[2] += size + int(f[1])
                        cur[3] += size + int(f[2])
                    else:
                        self._close(cur)
                        cur = None
            if cur is not None:
                self._close(cur)
        for lst in self.chains.values():
            lst.sort(key=lambda c: -c[0])

    def _close(self, cur) -> None:
        score, t_name, _, _, starts, ends, targets, keep = cur
        arrs = (np.array(starts, dtype=np.int64), np.array(ends, dtype=np.int64),
                np.array(targets, dtype=np.int64))
        # a chain we do not follow still claims its region, so a lower-scoring chain cannot
        self.chains.setdefault(t_name, []).append((score, *arrs) if keep else (score, arrs[0], arrs[1], None))

    def lift(self, chrom: int, pos1: np.ndarray) -> np.ndarray:
        p = np.asarray(pos1, dtype=np.int64) - 1
        out = np.full(len(p), -1, dtype=np.int64)
        done = np.zeros(len(p), dtype=bool)
        for _, starts, ends, targets in self.chains.get(f"chr{chrom}", []):
            todo = np.flatnonzero(~done)
            if not len(todo):
                break
            x = p[todo]
            k = np.searchsorted(starts, x, side="right") - 1
            hit = (k >= 0) & (x < ends[np.clip(k, 0, None)])
            idx = todo[hit]
            if targets is not None:
                kk = k[hit]
                out[idx] = targets[kk] + (x[hit] - starts[kk]) + 1
            done[idx] = True
        return out.astype(np.int32)
