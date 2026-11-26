"""Run every disorder for one genome and render the result as text or JSON."""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

from . import __version__
from .catalog import RISK, SOURCES
from .genome import Genome, read_genome
from .panel import Panel
from .score import (DEFAULT_P, DEFAULT_R2, DEFAULT_WINDOW, DisorderResult, GenomeIndex, ancestry_fit,
                    liability_risk, load_weights, score_disorder)

# ancestry_fit() is 1.0 for someone drawn from the European reference population.  Measured on
# the 2,504 people of 1000 Genomes with array SNPs: 95% of Europeans are below 1.037 and all
# are below 1.070; every South Asian is above 1.11, every East Asian above 1.36.
FIT_GOOD = 1.04
FIT_POOR = 1.075

FORMAT_NAMES = {"23andme": "23andMe-style export", "ancestry": "AncestryDNA export",
                "csv": "MyHeritage / FamilyTreeDNA export", "vcf": "VCF"}


@dataclass
class Report:
    genome: Genome
    results: list[DisorderResult]
    fit: float
    fit_snps: int
    p_threshold: float
    window: int
    valid: bool = True                 # False: the standardisation does not apply to this genome
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        g = self.genome
        return {
            "tool": f"pgc-risk {__version__}",
            "genome": {"file": g.path, "format": g.fmt, "build": g.build, "sample": g.sample,
                       "autosomal_snps": g.n, "called": int(g.called.sum()),
                       "variants_only": g.variants_only, "notes": g.notes},
            "ancestry_fit": {"value": self.fit, "snps": self.fit_snps,
                             "good_below": FIT_GOOD, "invalid_above": FIT_POOR},
            "settings": {"p_threshold": self.p_threshold, "clump_window_bp": self.window,
                         "clump_r2": DEFAULT_R2},
            "valid": self.valid,
            "warnings": self.warnings,
            "disorders": [r.to_dict() for r in self.results],
        }


def run(genome_path: str | Path, weights_dir: Path, keys: list[str] | None = None,
        sample: str | None = None, p_threshold: float = DEFAULT_P,
        window: int = DEFAULT_WINDOW, r2_scale: float = 1.0, panel_dir: Path | None = None,
        progress=None) -> Report:
    """Score one genome file.  `progress(stage, done, total, label)` is called as work advances."""
    panel_dir = panel_dir or weights_dir.parent / "panel"
    if not Panel.exists(panel_dir):
        raise FileNotFoundError(f"no LD panel in {panel_dir} - run pgc-risk build (or pgc-risk panel) first")
    panel = Panel(panel_dir)
    todo = [s for s in SOURCES if (weights_dir / f"{s.key}.npz").exists()]
    if progress:
        progress("reading", 0, len(todo), "Reading the genome file")
    genome = read_genome(genome_path, sample=sample)
    index = GenomeIndex(genome)
    results: list[DisorderResult] = []
    parts = []
    for i, src in enumerate(todo):
        if progress:
            progress("scoring", i, len(todo), src.name)
        w = load_weights(src.key, weights_dir)
        matched = index.match(w)
        # the ancestry check always uses every disorder's SNPs, so it does not depend on -d
        parts.append((w.rsid[matched[0]], matched[1][matched[0]], w.freq[matched[0]]))
        if keys and src.key not in keys:
            continue
        res, detail = score_disorder(index, w, panel, p_threshold, window, DEFAULT_R2, matched)
        res.name = src.name
        model = RISK.get(src.key)
        if model is not None:
            res.prevalence = model.prevalence
            res.r2_liability = model.r2_liability * r2_scale
            res.r2_basis = model.extra.get("basis", "")
            res.low_power = model.low_power
            res.probability = liability_risk(res.z, res.prevalence, res.r2_liability)
            res.relative_risk = res.probability / res.prevalence
        results.append(res)
    if not results:
        raise FileNotFoundError(
            f"no weight files in {weights_dir} - run pgc-risk fetch and pgc-risk build first")
    fit, fit_n = ancestry_fit(parts)
    rep = Report(genome, results, fit, fit_n, p_threshold, window)
    if fit_n < 2000:
        rep.valid = False
        rep.warnings.append(
            f"Only {fit_n:,} of this file's SNPs overlap the GWAS data - too few to score.")
    elif fit > FIT_POOR:
        rep.valid = False
        rep.warnings.append(
            f"Ancestry fit {fit:.3f}: this genome does not match the European-ancestry reference the "
            "scores are calibrated on (non-European or strongly admixed ancestry). In such genomes the "
            "scores sit several standard deviations off for reasons unrelated to risk, so no "
            "percentiles or probabilities are given. The PGC data used here cannot score this genome.")
    elif fit > FIT_GOOD:
        rep.warnings.append(
            f"Ancestry fit {fit:.3f}: this genome sits some distance from the European-ancestry "
            "reference (about 1 European in 20 does; so do partly non-European genomes). "
            "Treat the numbers below as rough.")
    if genome.fmt == "vcf" and genome.variants_only and index.pos_build is None:
        rep.valid = False
        rep.warnings.append(
            "The VCF lists variant sites only and its genome build could not be determined, so sites "
            "that match the reference could not be filled in. The scores would be biased and are not shown.")
    return rep


def _pct(x: float) -> str:
    v = 100 * x
    if v < 0.095:
        return f"{v:.2f}%"
    if v < 9.95:
        return f"{v:.1f}%"
    return f"{int(v + 0.5)}%"          # halves round up, as in the graphical interface


def _ordinal_percentile(p: float) -> str:
    v = min(99, max(1, round(100 * p)))
    suf = "th" if 10 <= v % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(v % 10, "th")
    return f"{v}{suf}"


def render_text(rep: Report) -> str:
    g = rep.genome
    out = [f"pgc-risk {__version__} - polygenic risk estimates from PGC GWAS summary statistics", ""]
    build = f"GRCh{g.build}" if g.build else "build not stated"
    sample = f", sample {g.sample}" if g.sample else ""
    out.append(f"Genome file  : {Path(g.path).name}  ({FORMAT_NAMES.get(g.fmt, g.fmt)}, {build}{sample})")
    out.append(f"Genotypes    : {int(g.called.sum()):,} autosomal SNPs with a call")
    verdict = ("matches" if rep.fit <= FIT_GOOD else
               "borderline" if rep.fit <= FIT_POOR else "DOES NOT MATCH")
    out.append(f"Ancestry fit : {rep.fit:.3f} over {rep.fit_snps:,} SNPs - {verdict} the European reference "
               f"(1.00 = typical, above {FIT_POOR:g} = not scorable)")
    for n in g.notes:
        out.append(f"Note         : {n}")
    for w in rep.warnings:
        out.append("")
        out.append(f"WARNING: {w}")
    out.append("")
    head = (f"{'Disorder':<34}{'Population':>11}{'Your':>10}{'Relative':>10}{'Score':>12}{'SNPs':>8}")
    sub = (f"{'':<34}{'baseline':>11}{'estimate':>10}{'to avg.':>10}{'percentile':>12}{'used':>8}")
    out += [head, sub, "-" * len(head)]
    for r in rep.results:
        flag = " *" if r.low_power else ""
        base = _pct(r.prevalence) if r.prevalence is not None else "-"
        if not rep.valid or r.probability is None:
            out.append(f"{r.name + flag:<34}{base:>11}{'-':>10}{'-':>10}{'-':>12}{r.snps_used:>8,}")
        else:
            out.append(f"{r.name + flag:<34}{base:>11}{_pct(r.probability):>10}"
                       f"{r.relative_risk:>9.2f}x{_ordinal_percentile(r.percentile):>12}{r.snps_used:>8,}")
        for w in r.warnings:
            out.append(f"    ! {w}")
    out.append("")
    if any(r.low_power for r in rep.results):
        out.append("* the GWAS behind this score is small: the score carries very little information,")
        out.append("  so the estimate stays close to the population baseline whatever the genotype.")
        out.append("")
    out += [
        "How to read this",
        "  Population baseline  lifetime risk for an average person.",
        "  Your estimate        lifetime risk for people with this polygenic score - a statistical",
        "                       expectation, not a diagnosis or a prediction for one person.",
        "  Score percentile     where the score falls among people of European ancestry.",
        "",
        "  Common-variant scores explain only a few percent of who develops these conditions",
        "  (about 7% of the liability for schizophrenia and bipolar disorder, 4% for depression,",
        "  2% or less for the rest).  Family history, environment and rare variants are not in",
        "  this estimate.  It is calibrated for European ancestry only, and within Europe",
        "  regional ancestry can shift a score by a few tenths of a standard deviation.",
        "  This is a research tool, not a medical test.",
    ]
    return "\n".join(out)


def render_json(rep: Report) -> str:
    return json.dumps(rep.to_dict(), indent=2)
