"""The browser scorer (web/pgcrisk-core.js, run under Node) against the Python bundle scorer.

Both read the same reference bundle and the same genome files; every number must agree.
Also checks the JavaScript normal distribution function against Python's, the handling of
zipped and gzipped exports, and the messages for files the browser version does not take.

usage: python tools/check_web.py [full|chips]     (needs Node 18+, and tools/make_test_genomes.py output)
"""
from __future__ import annotations

import gzip
import json
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path
from statistics import NormalDist

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pgc_risk.bundle import read_bundle, score_genome  # noqa: E402
from pgc_risk.genome import read_genome  # noqa: E402

GENOMES = ROOT / "data" / "validation" / "genomes"


def main() -> None:
    variant = sys.argv[1] if len(sys.argv) > 1 else "full"
    bundle_path = ROOT / "data" / "bundle" / f"pgc-risk-reference-{variant}-v1.pgcr"
    bundle = read_bundle(bundle_path)
    files = sorted(GENOMES.glob("*.23andme.txt"))
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        # the same person as AncestryDNA-style, MyHeritage-style, zipped and gzipped exports
        src = files[0]
        rows = [l.rstrip("\n").split("\t") for l in src.read_text().splitlines() if l and not l.startswith("#")]
        anc = tmp / "AncestryDNA.txt"
        anc.write_text("#AncestryDNA raw data download\n#human reference build 37\nrsid\tchromosome\tposition\tallele1\tallele2\n"
                       + "".join(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3][0]}\t{r[3][1]}\n" for r in rows))
        csv = tmp / "MyHeritage.csv"
        csv.write_text('# MyHeritage DNA raw data.\nRSID,CHROMOSOME,POSITION,RESULT\n'
                       + "".join(f'"{r[0]}","{r[1]}","{r[2]}","{r[3]}"\n' for r in rows))
        gz = tmp / "genome.txt.gz"
        gz.write_bytes(gzip.compress(src.read_bytes()))
        zp = tmp / "genome_export.zip"
        with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("readme.txt", "hello")
            zf.write(src, "genome_Full.txt")
        dup = tmp / "duplicates.txt"                              # an rs number listed twice, a no-call, an indel
        text = src.read_text().splitlines()
        dup.write_text("\n".join(text[:5000] + [text[4000], text[4500].rsplit("\t", 1)[0] + "\t--",
                                                 "rs999999999\t1\t5\tDI"] + text[5000:]) + "\n")
        refused = [GENOMES / files[0].name.replace("23andme.txt", "hg19.vcf"), ROOT / "README.md", bundle_path]
        todo = files + [anc, csv, gz, zp, dup]
        out = subprocess.run(["node", str(ROOT / "web" / "test" / "score.mjs"), str(bundle_path)]
                             + [str(f) for f in todo + refused], capture_output=True, text=True, check=True)
        lines = [json.loads(l) for l in out.stdout.splitlines()]
        worst = 0.0
        print(f"bundle: {variant}  ({bundle.header['snps']:,} SNPs)")
        print(f"{'file':<26}{'js ms':>7}{'fit':>8}{'max |dz|':>11}{'max |dprob|':>13}{'same SNP counts':>17}")
        for f, js in zip(todo, lines):
            py = score_genome(bundle, read_genome(f))
            dz = max(abs(a["z"] - b["z"]) for a, b in zip(js["disorders"], py["disorders"]))
            dp = max(abs(a["probability"] - b["probability"]) for a, b in zip(js["disorders"], py["disorders"]))
            same = all(a["snps_used"] == b["snps_used"] for a, b in zip(js["disorders"], py["disorders"]))
            dfit = abs(js["ancestry_fit"]["value"] - py["ancestry_fit"]["value"])
            ok = (js["valid"] == py["valid"] and js["ancestry_fit"]["snps"] == py["ancestry_fit"]["snps"]
                  and js["genome"]["called"] == py["genome"]["called"] and js["genome"]["format"] == py["genome"]["format"])
            worst = max(worst, dz, dp, dfit, 0 if (same and ok) else 1)
            print(f"{f.name:<26}{js['ms']:>7}{js['ancestry_fit']['value']:>8.3f}{dz:>11.2e}{dp:>13.2e}"
                  f"{str(same and ok):>17}")
        # the four re-encodings of the first person must give that person's result
        base = lines[0]
        for f, js in zip(todo[len(files):len(files) + 4], lines[len(files):len(files) + 4]):
            assert all(abs(a["z"] - b["z"]) < 1e-12 for a, b in zip(js["disorders"], base["disorders"])), f.name
        print("\nfiles the browser version does not take:")
        for f, js in zip(refused, lines[len(todo):]):
            print(f"  {f.name:<34} -> {js.get('error', 'SCORED (unexpected)')[:92]}")
            worst = max(worst, 0 if "error" in js else 1)
        nd = NormalDist()
        grid = lines[-1]["normalCdf"]
        err = max(abs(v - nd.cdf(-9 + 0.125 * i)) / max(nd.cdf(-9 + 0.125 * i), 1e-300) for i, v in enumerate(grid[:73]))
        err = max(err, max(abs(v - nd.cdf(-9 + 0.125 * i)) for i, v in enumerate(grid)))
        print(f"\nnormal distribution function, JavaScript vs Python on [-9, 9]: largest error {err:.1e}")
        worst = max(worst, err)
    print("PASS" if worst < 1e-9 else f"FAIL (largest discrepancy {worst:.2e})")
    raise SystemExit(0 if worst < 1e-9 else 1)


if __name__ == "__main__":
    main()
