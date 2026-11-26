"""Write real people's genotypes (1000 Genomes) as genome files in the formats the scorer reads.

For each requested sample three files are written to data/validation/genomes/:
  <id>.23andme.txt   array-style export: only SNPs found on genotyping arrays, rs numbers, GRCh37
  <id>.hg19.vcf      sequencing-style VCF: variant sites only, no rs numbers, GRCh37
  <id>.hg38.vcf.gz   the same lifted to GRCh38, without a contig header (build must be inferred)

Only the SNPs that occur in a weight file are written, which is all the scorer looks at.

usage: python tools/make_test_genomes.py NA12878 HG00096 ...
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pgc_risk import panel as pn  # noqa: E402

OUT = ROOT / "data" / "validation" / "genomes"
BASES = "ACGT"


def snp_table(c: int):
    """pos -> (rsid, pos38) for every weight SNP on a chromosome."""
    rs, p38 = {}, {}
    for f in sorted(pn.WEIGHTS.glob("*.npz")):
        z = np.load(f)
        m = z["chrom"] == c
        for p, r, q in zip(z["pos"][m].tolist(), z["rsid"][m].tolist(), z["pos38"][m].tolist()):
            rs[p] = r
            p38[p] = q
    return rs, p38


def main(ids: list[str]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    array_rs = set(np.load(ROOT / "data" / "validation" / "array_like_rsids.npy").tolist())
    files = {}
    for sid in ids:
        a = open(OUT / f"{sid}.23andme.txt", "w", newline="\n")
        a.write("# This data file was generated from 1000 Genomes phase 3 for testing.\n"
                "# We are using reference human assembly build 37 (also known as Annotation Release 104).\n"
                "# rsid\tchromosome\tposition\tgenotype\n")
        v19 = open(OUT / f"{sid}.hg19.vcf", "w", newline="\n")
        v19.write("##fileformat=VCFv4.2\n##contig=<ID=1,length=249250621>\n"
                  f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sid}\n")
        v38 = gzip.open(OUT / f"{sid}.hg38.vcf.gz", "wt", newline="\n")
        v38.write(f"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sid}\n")
        files[sid] = (a, v19, v38)
    for c in range(1, 23):
        z = np.load(pn._cache(c))
        names = z["samples"].tolist()
        rs, p38 = snp_table(c)
        pos, ref, alt = z["pos"].tolist(), z["ref"].tolist(), z["alt"].tolist()
        for sid in ids:
            g = z["geno"][:, names.index(sid)].tolist()
            a, v19, v38 = files[sid]
            lines38 = []
            for p, r, al, d in zip(pos, ref, alt, g):
                rsid = rs.get(p)
                if rsid is None:
                    continue
                rb, ab = BASES[r], BASES[al]
                if rsid in array_rs:
                    a.write(f"rs{rsid}\t{c}\t{p}\t{(rb + rb, rb + ab, ab + ab)[d]}\n")
                if d:
                    gt = "0/1" if d == 1 else "1/1"
                    v19.write(f"{c}\t{p}\t.\t{rb}\t{ab}\t50\tPASS\t.\tGT\t{gt}\n")
                    q = p38.get(p, -1)
                    if q > 0:
                        lines38.append((q, f"chr{c}\t{q}\t.\t{rb}\t{ab}\t50\tPASS\t.\tGT\t{gt}\n"))
            lines38.sort()
            v38.writelines(line for _, line in lines38)
        print(f"chr{c} done", flush=True)
    for fs in files.values():
        for fh in fs:
            fh.close()


if __name__ == "__main__":
    main(sys.argv[1:] or ["NA12878"])
