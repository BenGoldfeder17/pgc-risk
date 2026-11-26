"""Where each score's accuracy figure comes from; writes pgc_risk/accuracy.json.

The risk model needs, per disorder, the share of liability variance the score explains (R2).
In order of preference the figure is

  published   the GWAS paper's own hold-out result for a clumping + thresholding score,
  scaled      the paper's hold-out result for a stronger method (PRS-CS, MegaPRS), times 0.8,
  hold-out    this tool's own test on a held-out cohort (tools/loo.py),
  expected    what a GWAS of this size and heritability should deliver, from the standard
              relation R2 = h2 x / (1 + x), x = Neff h2_obs / M (Daetwyler 2008, Wray 2013).

M, the effective number of independent SNPs for this kind of score, is calibrated on the three
cases where a clumping + thresholding result is known: the published depression result and
this tool's two schizophrenia hold-outs.  They agree (about 85,000).

usage: python tools/accuracy.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from statistics import NormalDist

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from pgc_risk.catalog import PREVALENCE, SOURCES  # noqa: E402

ND = NormalDist()
WEIGHTS = ROOT / "data" / "weights"
LOO = ROOT / "data" / "validation" / "loo"
OUT = ROOT / "pgc_risk" / "accuracy.json"
SIMPLE_VS_BAYESIAN = 0.8      # a clumping score relative to PRS-CS / MegaPRS on the same data


def obs_per_liab(k: float) -> float:
    """h2 on the observed scale of a balanced case-control sample per unit liability-scale h2."""
    z = ND.pdf(ND.inv_cdf(1 - k))
    return z * z * 0.25 / (k * (1 - k)) ** 2


def neff(key: str, directory: Path = WEIGHTS) -> float:
    return float(json.loads((directory / f"{key}.json").read_text())["qc"]["effective_n_p90"])


def implied_m(r2: float, h2_liab: float, k: float, n: float) -> float:
    x = (r2 / h2_liab) / (1 - r2 / h2_liab)
    return n * h2_liab * obs_per_liab(k) / x


def expected_r2(h2_liab_paper: float, k_paper: float, k_use: float, n: float, m: float) -> float:
    h2_obs = h2_liab_paper * obs_per_liab(k_paper)
    h2_liab = h2_obs / obs_per_liab(k_use)
    x = n * h2_obs / m
    return h2_liab * x / (1 + x)


# SNP heritability on the liability scale and the prevalence it was computed for, as published.
H2 = {
    "scz": (0.24, 0.01),      # Trubetskoy 2022 (European samples)
    "mdd": (0.084, 0.15),     # Adams 2025
    "adhd": (0.14, 0.05),     # Demontis 2023
    "asd": (0.118, 0.0122),   # Grove 2019
    "ts": (0.21, 0.008),      # Yu 2019
    "an": (0.20, 0.01),       # Duncan 2017 (prevalence taken as the ~1% quoted in the paper)
    "alc": (0.090, 0.159),    # Walters 2018 (prevalence among alcohol-exposed people)
    "cud": (0.121, 0.085),    # Johnson 2020 (upper end of the range the paper reports)
}


def main() -> None:
    anchors = [
        ("depression, published (Adams 2025)", implied_m(0.039, *H2["mdd"], neff("mdd"))),
        ("schizophrenia, men -> women hold-out", implied_m(0.0665, *H2["scz"], neff("scz_male", LOO))),
        ("schizophrenia, women -> men hold-out", implied_m(0.0526, *H2["scz"], neff("scz_female", LOO))),
    ]
    for label, m in anchors:
        print(f"M implied by {label:<40} {m:>10,.0f}")
    m_eff = round(sum(m for _, m in anchors) / len(anchors), -3)
    print(f"M used: {m_eff:,.0f}\n")

    def expected(key):
        h2, k_paper = H2[key]
        return expected_r2(h2, k_paper, PREVALENCE[key][0], neff(key), m_eff)

    acc = {
        "scz": dict(r2_liability=0.073, basis="published",
                    source="Trubetskoy 2022: clumping + thresholding at p<0.05 explained a median 7.3% of "
                           "liability variance in left-out cohorts. This tool's own hold-outs (trained on one "
                           "sex, tested on the other, about half the sample): 6.7% and 5.3%."),
        "mdd": dict(r2_liability=0.039, basis="published",
                    source="Adams 2025: 3.9% for clumping + thresholding at p<=0.05 from this same public file "
                           "(mean of 42 left-out cohorts). This tool's own hold-out on UK Biobank alone: 2.2%."),
        "bip": dict(r2_liability=round(0.084 * SIMPLE_VS_BAYESIAN, 4), basis="scaled",
                    source="O'Connell 2025: 8.4% with PRS-CS-auto in clinical cohorts, from the European "
                           "GWAS without self-reported cases; times 0.8 for a clumping score."),
        "anx": dict(r2_liability=round(0.0227 * SIMPLE_VS_BAYESIAN, 4), basis="scaled",
                    source="Strom 2026: 2.27% with MegaPRS in UK Biobank; times 0.8 for a clumping score."),
        "bpd": dict(r2_liability=round((0.031 + 0.026) / 2 * SIMPLE_VS_BAYESIAN, 4), basis="scaled",
                    source="Streit 2026: 3.1% and 2.6% with PRS-CS in two independent replication cohorts "
                           "(4.6% in the discovery cohorts); mean times 0.8 for a clumping score."),
        "ocd": dict(r2_liability=0.008, basis="hold-out",
                    source="This tool's own hold-out: trained without UK Biobank, tested on UK Biobank "
                           "(776 cases): 0.7% with all SNPs, 1.0% with array SNPs. No published figure."),
        "ptsd": dict(r2_liability=0.0015, basis="published",
                     source="Nievergelt 2019: 0.15% between independent PTSD data sets (clumping score)."),
    }
    notes = {
        "adhd": "Demontis 2023 gives no hold-out figure on the liability scale.",
        "asd": "Grove 2019 reports 2.45% on the Nagelkerke scale, which corresponds to roughly 1% of liability.",
        "ts": "Yu 2019 gives no hold-out figure.",
        "an": "Duncan 2017 gives no hold-out figure; the far larger 2019 GWAS (not on the Hub) reached 1.7%.",
        "alc": "Walters 2018 gives no liability-scale hold-out figure.",
        "cud": "Johnson 2020 gives no hold-out figure for the disorder itself.",
    }
    for key, note in notes.items():
        r2 = expected(key)
        acc[key] = dict(r2_liability=round(r2, 4), basis="expected",
                        source=f"Expected from GWAS size (Neff {neff(key):,.0f}) and SNP heritability "
                               f"({H2[key][0]:g} at prevalence {H2[key][1]:g}) with M = {m_eff:,.0f}. {note}")
    # the same formula for the disorders that do have a figure, as a cross-check
    print(f"{'':<6}{'K':>7}{'R2 used':>9}  basis      (expected from GWAS size, where h2 is known)")
    out = {}
    for src in SOURCES:
        a = acc[src.key]
        check = f"{100 * expected(src.key):.2f}%" if src.key in H2 else "-"
        print(f"{src.key:<6}{PREVALENCE[src.key][0]:>7.3g}{100 * a['r2_liability']:>8.2f}%  {a['basis']:<10} {check}")
        out[src.key] = a
    out["_method"] = {"m_effective": m_eff, "m_anchors": {k: round(v) for k, v in anchors},
                      "simple_vs_bayesian": SIMPLE_VS_BAYESIAN}
    OUT.write_text(json.dumps(out, indent=1))
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
