"""Command line: pgc-risk fetch | build | gui | score <genome file>."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .catalog import SOURCES


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    # `pgc-risk my_genome.txt` is shorthand for `pgc-risk score my_genome.txt`
    if argv and argv[0] not in ("fetch", "build", "panel", "score", "gui", "list", "-h", "--help", "--version"):
        argv.insert(0, "score")
    keys_help = "disorder keys (default: all): " + ", ".join(s.key for s in SOURCES)
    ap = argparse.ArgumentParser(prog="pgc-risk", description=__doc__)
    ap.add_argument("--version", action="version", version=f"pgc-risk {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fetch", help="download the PGC summary statistics from Hugging Face")
    f.add_argument("keys", nargs="*", help=keys_help)
    f.add_argument("--workers", type=int, default=6)

    b = sub.add_parser("build", help="turn the downloaded shards into weight files, then build the LD panel")
    b.add_argument("keys", nargs="*", help=keys_help)
    b.add_argument("--no-panel", action="store_true", help="skip the 1000 Genomes LD panel step")
    b.add_argument("--jobs", type=int, default=8)

    pn = sub.add_parser("panel", help="(re)build the European LD panel from 1000 Genomes (about 15 GB download)")
    pn.add_argument("--also", action="append", default=[], metavar="DIR", help="extra weight directory to cover")
    pn.add_argument("--jobs", type=int, default=8)

    sub.add_parser("list", help="show the disorders and the GWAS behind each")

    gu = sub.add_parser("gui", help="open the graphical interface in the browser")
    gu.add_argument("--port", type=int, default=8765)
    gu.add_argument("--no-browser", action="store_true", help="start the server without opening a browser")
    gu.add_argument("--lan", action="store_true",
                    help="also answer other devices on the local network (HTTPS, with an access code)")
    gu.add_argument("--no-encryption", action="store_true",
                    help="with --lan: serve plain HTTP on the network instead of HTTPS (not recommended)")
    gu.add_argument("--weights", default=None, help="directory holding the built weight files")

    s = sub.add_parser("score", help="estimate disorder probabilities for a genome file")
    s.add_argument("genome", help="23andMe / AncestryDNA / MyHeritage / FTDNA export, or a VCF")
    s.add_argument("-d", "--disorders", help="comma-separated keys (default: all)")
    s.add_argument("--sample", help="sample name to score in a multi-sample VCF")
    s.add_argument("--json", metavar="FILE", help="also write the full result as JSON ('-' = stdout)")
    s.add_argument("--r2-scale", type=float, default=1.0,
                   help="multiply every score's assumed accuracy (R2) by this factor, to see how "
                        "much the estimates depend on it")
    s.add_argument("--weights", default=None, help="directory holding the built weight files")
    s.add_argument("--bundle", metavar="FILE",
                   help="score from a downloaded reference bundle (.pgcr) instead of the full local build; "
                        "genotyping-chip files only")
    args = ap.parse_args(argv)

    if args.cmd == "fetch":
        from .fetch import fetch
        fetch(args.keys or None, workers=args.workers)
        return 0
    if args.cmd == "build":
        from .build import build
        build(args.keys or None)
        if not args.no_panel:
            from .panel import build_panel
            build_panel(jobs=args.jobs)
        return 0
    if args.cmd == "panel":
        from .panel import WEIGHTS, build_panel
        build_panel([WEIGHTS] + [Path(d) for d in args.also], jobs=args.jobs)
        return 0
    if args.cmd == "gui":
        from .fetch import DATA
        from .gui import serve
        return serve(Path(args.weights) if args.weights else DATA / "weights", args.port, not args.no_browser,
                     lan=args.lan, encrypt=not args.no_encryption)
    if args.cmd == "list":
        from .catalog import RISK
        for src in SOURCES:
            m = RISK.get(src.key)
            extra = f"  baseline {100 * m.prevalence:g}%, score R2 {100 * m.r2_liability:g}%" if m else ""
            print(f"{src.key:<5} {src.name:<32} {src.citation}{extra}")
        return 0

    from .fetch import DATA
    from .genome import GenomeFormatError
    from .report import render_json, render_text, run

    if args.bundle:
        from .bundle import BundleError, read_bundle, render_text as render_bundle, score_genome
        from .genome import read_genome
        try:
            bundle = read_bundle(args.bundle)
            res = score_genome(bundle, read_genome(args.genome, sample=args.sample), args.r2_scale)
        except (GenomeFormatError, BundleError, FileNotFoundError) as e:
            print(f"error: {e}", file=sys.stderr)
            return 2
        import json
        if args.json == "-":
            print(json.dumps(res, indent=2))
        else:
            print(render_bundle(res, bundle))
            if args.json:
                Path(args.json).write_text(json.dumps(res, indent=2), encoding="utf-8")
                print(f"\nJSON written to {args.json}")
        return 0
    weights = Path(args.weights) if args.weights else DATA / "weights"
    keys = [k.strip() for k in args.disorders.split(",")] if args.disorders else None
    unknown = [k for k in keys or [] if k not in {s.key for s in SOURCES}]
    if unknown:
        ap.error(f"unknown disorder key(s): {', '.join(unknown)}")
    try:
        rep = run(args.genome, weights, keys, sample=args.sample, r2_scale=args.r2_scale)
    except (GenomeFormatError, FileNotFoundError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    if args.json == "-":
        print(render_json(rep))
    else:
        print(render_text(rep))
        if args.json:
            Path(args.json).write_text(render_json(rep), encoding="utf-8")
            print(f"\nJSON written to {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
