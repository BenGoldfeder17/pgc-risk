"""Read a personal genome file into arrays of autosomal SNP genotypes.

Supported: 23andMe, AncestryDNA, MyHeritage, FamilyTreeDNA and LivingDNA raw-data exports
(plain, .gz or .zip) and VCF / gVCF (plain or .gz, single- or multi-sample).

Alleles are coded A=0, C=1, G=2, T=3, so the complementary base is ``3 - code`` and a
strand-ambiguous pair (A/T, C/G) sums to 3.  A missing call is -1.
"""
from __future__ import annotations

import gzip
import io
import re
import zipfile
from array import array
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

BASE = {"A": 0, "C": 1, "G": 2, "T": 3}
AUTOSOMES = {str(i): i for i in range(1, 23)}
AUTOSOMES.update({f"chr{i}": i for i in range(1, 23)})

# chr1 length identifies the assembly in a VCF header.
_CHR1_LEN = {249250621: 37, 248956422: 38, 247249719: 36}


class GenomeFormatError(ValueError):
    pass


@dataclass
class Genome:
    path: str
    fmt: str                      # "23andme" | "ancestry" | "csv" | "vcf"
    build: int | None             # 36 / 37 / 38, or None when the file does not say
    rsid: np.ndarray              # int64, -1 when the record has no rs number
    chrom: np.ndarray             # int8, 1..22
    pos: np.ndarray               # int32
    a1: np.ndarray                # int8 allele codes; -1 = no call
    a2: np.ndarray
    sample: str = ""
    # VCF only.  A variants-only VCF lists just the sites that differ from the reference, so
    # a SNP that is absent is homozygous reference; in an all-sites file it is simply missing.
    variants_only: bool = False
    ref_blocks: dict[int, tuple[np.ndarray, np.ndarray]] = field(default_factory=dict)  # gVCF
    notes: list[str] = field(default_factory=list)

    @property
    def n(self) -> int:
        return len(self.rsid)

    @property
    def called(self) -> np.ndarray:
        return (self.a1 >= 0) & (self.a2 >= 0)


def _open_text(path: Path) -> io.TextIOBase:
    with open(path, "rb") as fh:
        magic = fh.read(4)
    if magic[:2] == b"\x1f\x8b":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    if magic[:2] == b"PK":
        zf = zipfile.ZipFile(path)
        members = [m for m in zf.infolist() if not m.is_dir()]
        if not members:
            raise GenomeFormatError(f"{path.name}: the zip archive is empty")
        biggest = max(members, key=lambda m: m.file_size)
        return io.TextIOWrapper(zf.open(biggest), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def _rs_number(text: str) -> int:
    if len(text) > 2 and text[0] in "rR" and text[1] in "sS":
        digits = text[2:]
        if digits.isdigit():
            return int(digits)
        m = re.match(r"\d+", digits)          # "rs123;rs456" or "rs123:1000:A:G"
        if m:
            return int(m.group())
    return -1


def _build_from_text(text: str) -> int | None:
    t = text.lower()
    for build, marks in ((38, ("grch38", "hg38", "build 38", "b38", "hs38", "_v38")),
                         (37, ("grch37", "hg19", "build 37", "b37", "hs37", "_v37")),
                         (36, ("ncbi36", "hg18", "build 36", "b36"))):
        if any(m in t for m in marks):
            return build
    return None


def read_genome(path: str | Path, sample: str | None = None) -> Genome:
    path = Path(path)
    if not path.is_file():
        raise GenomeFormatError(f"{path}: no such file")
    fh = _open_text(path)
    try:
        head = []
        for line in fh:
            head.append(line)
            if len(head) >= 400 or (line.strip() and not line.startswith("#")):
                break
        fh.close()
        fh = _open_text(path)
        first = head[0] if head else ""
        if first.startswith("##fileformat=VCF") or any(l.startswith("#CHROM") for l in head):
            return _read_vcf(path, fh, sample)
        return _read_array(path, fh)
    finally:
        fh.close()


# --------------------------------------------------------------------------- array exports
def _read_array(path: Path, fh) -> Genome:
    rsid, chrom, pos = array("q"), array("b"), array("i")
    a1, a2 = array("b"), array("b")
    comments: list[str] = []
    fmt = None
    delim = "\t"
    n_cols = 0
    get = BASE.get
    autos = AUTOSOMES
    for line in fh:
        if not line or line[0] == "#":
            if len(comments) < 60:
                comments.append(line)
            continue
        if fmt is None:
            stripped = line.strip()
            if not stripped:
                continue
            delim = "\t" if "\t" in stripped else ("," if "," in stripped else None)
            cols = [c.strip().strip('"') for c in (stripped.split(delim) if delim else stripped.split())]
            n_cols = len(cols)
            if n_cols not in (4, 5):
                raise GenomeFormatError(
                    f"{path.name}: not a recognised genome file "
                    f"(expected a 23andMe/AncestryDNA/MyHeritage/FTDNA export or a VCF)")
            fmt = "ancestry" if n_cols == 5 else ("csv" if delim == "," else "23andme")
            if cols[0].lower() in ("rsid", "rs", "snp", "name", "# rsid"):
                continue                       # header row
        parts = line.rstrip("\r\n").split(delim) if delim else line.split()
        if len(parts) != n_cols:
            continue
        if delim == ",":
            parts = [p.strip().strip('"') for p in parts]
        c = autos.get(parts[1])
        if c is None:
            continue
        try:
            p = int(parts[2])
        except ValueError:
            continue
        if n_cols == 5:
            g1, g2 = get(parts[3].strip().upper(), -1), get(parts[4].strip().upper(), -1)
        else:
            g = parts[3].strip().upper()
            if len(g) == 2:
                g1, g2 = get(g[0], -1), get(g[1], -1)
            else:
                g1 = g2 = -1
        if g1 < 0 or g2 < 0:
            g1 = g2 = -1
        rsid.append(_rs_number(parts[0]))
        chrom.append(c)
        pos.append(p)
        a1.append(g1)
        a2.append(g2)
    if fmt is None or not len(rsid):
        raise GenomeFormatError(
            f"{path.name}: no autosomal SNP genotypes found (expected a 23andMe/AncestryDNA/MyHeritage/"
            f"FTDNA export or a VCF)")
    build = _build_from_text(" ".join(comments)) or 37
    g = Genome(str(path), fmt, build,
               np.frombuffer(rsid, dtype=np.int64).copy(), np.frombuffer(chrom, dtype=np.int8).copy(),
               np.frombuffer(pos, dtype=np.int32).copy(), np.frombuffer(a1, dtype=np.int8).copy(),
               np.frombuffer(a2, dtype=np.int8).copy())
    if build != 37:
        g.notes.append(f"File header says genome build {build}; SNPs are matched by rs number only.")
    return g


# --------------------------------------------------------------------------- VCF / gVCF
def _read_vcf(path: Path, fh, sample: str | None) -> Genome:
    rsid, chrom, pos = array("q"), array("b"), array("i")
    a1, a2 = array("b"), array("b")
    blocks: dict[int, tuple[array, array]] = {}
    header_text: list[str] = []
    build = None
    samples: list[str] = []
    col = 9
    n_homref = n_filtered = n_records = 0
    is_gvcf = False
    get = BASE.get
    autos = AUTOSOMES
    for line in fh:
        if line[0] == "#":
            if line.startswith("##"):
                if line.startswith("##contig") and build is None:
                    m = re.search(r"ID=(?:chr)?1,.*?length=(\d+)", line)
                    if m:
                        build = _CHR1_LEN.get(int(m.group(1)))
                elif len(header_text) < 200 and line.startswith(("##reference", "##assembly")):
                    header_text.append(line)
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\r\n").split("\t")[9:]
                if not samples:
                    raise GenomeFormatError(f"{path.name}: the VCF has no sample (genotype) column")
                if sample is not None:
                    if sample not in samples:
                        raise GenomeFormatError(
                            f"{path.name}: sample {sample!r} not in the VCF "
                            f"(it has: {', '.join(samples[:8])}{' ...' if len(samples) > 8 else ''})")
                    col = 9 + samples.index(sample)
                else:
                    col = 9
                continue
            continue
        parts = line.split("\t", col + 1)
        if len(parts) <= col:
            continue
        c = autos.get(parts[0])
        if c is None:
            continue
        ref = parts[3]
        if len(ref) != 1:
            continue
        r = get(ref.upper(), -1)
        if r < 0:
            continue
        alt = parts[4]
        fmt_keys = parts[8]
        gt_field = parts[col].rstrip("\r\n")
        if fmt_keys != "GT" and not fmt_keys.startswith("GT:"):
            keys = fmt_keys.split(":")
            if "GT" not in keys:
                continue
            gt = gt_field.split(":")[keys.index("GT")]
        else:
            gt = gt_field.split(":", 1)[0]
        n_records += 1
        p = int(parts[1])
        flt = parts[6]
        ok = flt == "PASS" or flt == "."
        if alt[0] == "<" or alt == ".":
            # gVCF reference block or a monomorphic all-sites record
            sep = "|" if "|" in gt else "/"
            homref = gt in ("0/0", "0|0")
            if alt[0] == "<" and "END=" in parts[7]:
                is_gvcf = True
                if homref and ok:
                    m = re.search(r"(?:^|;)END=(\d+)", parts[7])
                    if m:
                        b = blocks.setdefault(c, (array("i"), array("i")))
                        b[0].append(p)
                        b[1].append(int(m.group(1)))
                continue
            if "," in gt or sep not in gt:
                continue
            g1 = g2 = r if (homref and ok) else -1
            n_homref += homref
        else:
            if "|" in gt:
                i1, _, i2 = gt.partition("|")
            else:
                i1, _, i2 = gt.partition("/")
            if not ok:
                n_filtered += 1
                g1 = g2 = -1
            elif not i1.isdigit() or not i2.isdigit():
                g1 = g2 = -1                       # "./." or haploid
            else:
                if i1 == "0" and i2 == "0":
                    n_homref += 1
                alts = alt.split(",") if "," in alt else (alt,)
                k1, k2 = int(i1), int(i2)
                try:
                    s1 = ref if k1 == 0 else alts[k1 - 1]
                    s2 = ref if k2 == 0 else alts[k2 - 1]
                except IndexError:
                    s1 = s2 = ""
                g1 = get(s1.upper(), -1) if len(s1) == 1 else -1
                g2 = get(s2.upper(), -1) if len(s2) == 1 else -1
                if g1 < 0 or g2 < 0:
                    g1 = g2 = -1                   # indel, symbolic allele or '*'
        rsid.append(_rs_number(parts[2]))
        chrom.append(c)
        pos.append(p)
        a1.append(g1)
        a2.append(g2)
    if not len(rsid) and not blocks:
        raise GenomeFormatError(f"{path.name}: no autosomal SNP records with genotypes found")
    if build is None:
        build = _build_from_text(" ".join(header_text))
    g = Genome(str(path), "vcf", build,
               np.frombuffer(rsid, dtype=np.int64).copy(), np.frombuffer(chrom, dtype=np.int8).copy(),
               np.frombuffer(pos, dtype=np.int32).copy(), np.frombuffer(a1, dtype=np.int8).copy(),
               np.frombuffer(a2, dtype=np.int8).copy(),
               sample=samples[col - 9] if samples else "")
    for c, (s, e) in blocks.items():
        s_arr = np.frombuffer(s, dtype=np.int32).copy()
        e_arr = np.frombuffer(e, dtype=np.int32).copy()
        order = np.argsort(s_arr, kind="stable")
        g.ref_blocks[c] = (s_arr[order], e_arr[order])
    # Sequencing VCFs list variant sites only; array-derived / imputed VCFs list every site.
    homref_share = n_homref / max(1, n_records)
    g.variants_only = (not is_gvcf) and homref_share < 0.05
    if is_gvcf:
        g.notes.append("gVCF: sites inside reference blocks are read as homozygous reference.")
    elif g.variants_only:
        g.notes.append("Variants-only VCF: SNPs absent from the file are read as homozygous reference.")
    if n_filtered:
        g.notes.append(f"{n_filtered:,} calls failing the VCF FILTER column were treated as missing.")
    if len(samples) > 1:
        g.notes.append(f"Multi-sample VCF: scored sample {g.sample!r} "
                       f"(choose another with --sample; {len(samples)} available).")
    return g
