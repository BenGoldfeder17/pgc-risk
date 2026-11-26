"""LD panel storage, LD clumping, the LD variance factor and the full report on synthetic data."""
import json

import numpy as np
import pytest

from pgc_risk import report
from pgc_risk.panel import Panel, annotate
from pgc_risk.score import Weights, ld_clump, ld_inflation

A, C, G, T = 0, 1, 2, 3
N = 300            # panel samples


def write_panel(directory, chrom_data):
    """chrom_data: {chrom: (pos, ref, alt, genotype matrix)} -> panel files as panel.write_panel makes them."""
    directory.mkdir(parents=True, exist_ok=True)
    total = 0
    for c, (pos, ref, alt, g) in chrom_data.items():
        np.savez(directory / f"chr{c}.npz", pos=np.asarray(pos, np.int32), ref=np.asarray(ref, np.int8),
                 alt=np.asarray(alt, np.int8), p1=np.packbits(g >= 1, axis=1), p2=np.packbits(g == 2, axis=1))
        total += len(pos)
    (directory / "panel.json").write_text(json.dumps({"samples": N, "snps": total, "source": "synthetic"}))
    return Panel(directory)


def blocks(rng, n_blocks, per_block, noise):
    """Genotypes in LD blocks: every SNP of a block is a noisy copy of the block's lead SNP."""
    rows = []
    for _ in range(n_blocks):
        lead = rng.binomial(2, rng.uniform(0.2, 0.8), N)
        for _ in range(per_block):
            g = lead.copy()
            flip = rng.random(N) < noise
            g[flip] = rng.binomial(2, 0.5, int(flip.sum()))
            rows.append(g)
    return np.array(rows, dtype=np.int8)


def weights_for(chrom, pos, ref, alt, beta, p, freq):
    n = len(pos)
    return Weights("syn", {"name": "Synthetic"}, rsid=np.arange(1, n + 1, dtype=np.int64),
                   chrom=np.full(n, chrom, np.int8), pos=np.asarray(pos, np.int32),
                   pos38=np.asarray(pos, np.int32) + 1000, a1=np.asarray(alt, np.int8),
                   a2=np.asarray(ref, np.int8), ref=np.asarray(ref, np.int8),
                   beta=np.asarray(beta, float), p=np.asarray(p, float), freq=np.asarray(freq, float),
                   prow=np.arange(n, dtype=np.int32), palt=np.ones(n, bool))


def test_panel_roundtrip_and_lookup(tmp_path):
    rng = np.random.default_rng(1)
    g = rng.integers(0, 3, (50, N)).astype(np.int8)
    pos = np.arange(1, 51) * 1000
    ref, alt = np.full(50, A), np.full(50, G)
    panel = write_panel(tmp_path / "panel", {3: (pos, ref, alt, g)})
    assert np.array_equal(panel.dosage(3, np.arange(50)), g)
    assert np.array_equal(panel.dosage(3, np.array([7, 2])), g[[7, 2]])
    q = np.array([5000, 6000, 7000, 8000, 999])
    row, is_alt = panel.lookup(3, q, np.array([G, A, C, C, G]), np.array([A, G, T, A, A]))
    #            effect=ALT  effect=REF  other strand (C/T = G/A)  wrong alleles  not in panel
    assert row.tolist() == [4, 5, 6, -1, -1]
    assert is_alt.tolist() == [True, False, True, False, False]


def test_ld_clump_keeps_one_snp_per_ld_block(tmp_path):
    rng = np.random.default_rng(2)
    g = blocks(rng, n_blocks=40, per_block=5, noise=0.05)          # r2 within a block is about 0.8
    n = len(g)
    pos = 100_000_000 + np.arange(n) * 2000                         # everything within one window
    panel = write_panel(tmp_path / "panel", {1: (pos, np.full(n, A), np.full(n, G), g)})
    p = rng.uniform(1e-6, 0.05, n)
    w = weights_for(1, pos, np.full(n, A), np.full(n, G), rng.normal(0, 0.05, n), p, np.full(n, 0.5))
    sel = ld_clump(w, np.arange(n), panel, window=250_000, r2=0.1)
    assert len(sel) == 40
    assert sorted(sel // 5) == list(range(40))                      # one per block ...
    best = p.reshape(40, 5).argmin(axis=1) + np.arange(40) * 5
    assert sorted(sel.tolist()) == sorted(best.tolist())            # ... and it is the most significant


def test_ld_clump_ignores_ld_beyond_the_window(tmp_path):
    rng = np.random.default_rng(3)
    lead = rng.binomial(2, 0.5, N).astype(np.int8)
    g = np.array([lead, lead, lead], dtype=np.int8)                 # three identical SNPs
    pos = np.array([50_000_000, 50_100_000, 51_000_000])            # the third is 1 Mb away
    panel = write_panel(tmp_path / "panel", {2: (pos, np.full(3, C), np.full(3, T), g)})
    w = weights_for(2, pos, np.full(3, C), np.full(3, T), [0.1, 0.1, 0.1], [1e-9, 1e-3, 1e-2], [0.5] * 3)
    assert sorted(ld_clump(w, np.arange(3), panel, window=250_000).tolist()) == [0, 2]


def test_ld_inflation(tmp_path):
    rng = np.random.default_rng(4)
    indep = rng.binomial(2, 0.4, (400, N)).astype(np.int8)
    pos = 10_000_000 + np.arange(400) * 300_000
    panel = write_panel(tmp_path / "a", {4: (pos % 180_000_000 + np.arange(400), np.full(400, A), np.full(400, G), indep)})
    w = weights_for(4, np.sort(pos % 180_000_000 + np.arange(400)), np.full(400, A), np.full(400, G),
                    rng.normal(0, 0.05, 400), np.full(400, 0.01), np.full(400, 0.4))
    assert ld_inflation(w, np.arange(400), panel) == pytest.approx(1.0, abs=0.25)   # unlinked SNPs
    lead = rng.binomial(2, 0.5, N).astype(np.int8)
    twins = np.array([lead, lead], dtype=np.int8)
    panel2 = write_panel(tmp_path / "b", {5: ([1_000_000, 1_050_000], [A, A], [G, G], twins)})
    same = weights_for(5, [1_000_000, 1_050_000], [A, A], [G, G], [0.1, 0.1], [0.01, 0.01], [0.5, 0.5])
    assert ld_inflation(same, np.arange(2), panel2) == pytest.approx(2.0)           # var(2x) = 2 * (var x + var x)
    opposite = weights_for(5, [1_000_000, 1_050_000], [A, A], [G, G], [0.1, -0.1], [0.01, 0.01], [0.5, 0.5])
    assert ld_inflation(opposite, np.arange(2), panel2) == 0.25                     # cancels; floored


def _make_project(tmp_path, rng, n_blocks=1600):
    """A weight file + panel where SNPs come in LD pairs, and a way to draw people from the panel's LD."""
    g = blocks(rng, n_blocks=n_blocks, per_block=2, noise=0.1)
    n = len(g)
    per_chrom = n // 4
    chrom = np.repeat(np.arange(1, 5), per_chrom).astype(np.int8)
    pos = np.tile(1_000_000 + np.arange(per_chrom) * 60_000, 4).astype(np.int32)
    ref, alt = np.full(n, A, np.int8), np.full(n, G, np.int8)
    panel_dir = tmp_path / "panel"
    write_panel(panel_dir, {c: (pos[chrom == c], ref[chrom == c], alt[chrom == c], g[chrom == c]) for c in range(1, 5)})
    wdir = tmp_path / "weights"
    wdir.mkdir()
    beta = np.repeat(rng.normal(0, 0.05, n // 2), 2)                # LD partners share their effect
    np.savez(wdir / "scz.npz", rsid=np.arange(1, n + 1, dtype=np.int64), chrom=chrom, pos=pos, pos38=pos + 5,
             a1=alt, a2=ref, ref=ref, beta=beta.astype(np.float32), se=np.full(n, 0.02, np.float32),
             p=rng.uniform(1e-6, 0.05, n), freq=(g.mean(axis=1) / 2).astype(np.float32))
    (wdir / "scz.json").write_text(json.dumps({"name": "Schizophrenia"}))
    annotate(wdir, Panel(panel_dir), log=lambda *a: None)
    return wdir, panel_dir, g, chrom, pos


def _write_23andme(path, chrom, pos, dosage):
    lines = ["# rsid\tchromosome\tposition\tgenotype"]
    for i, (c, p, d) in enumerate(zip(chrom.tolist(), pos.tolist(), dosage.tolist())):
        lines.append(f"rs{i + 1}\t{c}\t{p}\t{('AA', 'AG', 'GG')[d]}")
    path.write_text("\n".join(lines) + "\n")


def test_report_end_to_end_and_ancestry_gate(tmp_path, monkeypatch):
    rng = np.random.default_rng(5)
    wdir, panel_dir, g, chrom, pos = _make_project(tmp_path, rng)
    zs, fits = [], []
    f = tmp_path / "g.txt"
    with monkeypatch.context() as m:
        # the real thresholds assume ~500,000 SNPs; with 3,200 the statistic is far noisier
        m.setattr(report, "FIT_GOOD", 1.25)
        m.setattr(report, "FIT_POOR", 1.3)
        for person in range(0, N, 3):                               # panel members are valid "Europeans"
            _write_23andme(f, chrom, pos, g[:, person])
            rep = report.run(f, wdir, panel_dir=panel_dir)
            assert rep.valid and len(rep.results) == 1
            r = rep.results[0]
            assert r.snps_used == len(g) // 2                       # one SNP per LD pair survives clumping
            assert 0 < r.probability < 1 and r.prevalence == 0.01
            zs.append(r.z)
            fits.append(rep.fit)
    assert 0.97 < np.mean(fits) < 1.03
    zs = np.array(zs)
    assert abs(zs.mean()) < 0.35 and 0.75 < zs.std() < 1.25
    text = report.render_text(rep)
    assert "Schizophrenia" in text and "the European reference" in text
    assert json.loads(report.render_json(rep))["disorders"][0]["key"] == "scz"
    # a genome whose allele frequencies are far from the reference is not scored
    far = rng.binomial(2, np.clip(g.mean(axis=1) / 2 + rng.choice([-0.3, 0.3], len(g)), 0.02, 0.98))
    _write_23andme(f, chrom, pos, far)
    rep = report.run(f, wdir, panel_dir=panel_dir)
    assert not rep.valid and rep.fit > report.FIT_POOR
    line = next(l for l in report.render_text(rep).splitlines() if l.startswith("Schizophrenia"))
    assert "%" in line and line.split()[2] == "-"                   # baseline shown, estimate withheld


def test_run_without_panel_explains_what_to_do(tmp_path):
    (tmp_path / "weights").mkdir()
    with pytest.raises(FileNotFoundError, match="pgc-risk build"):
        report.run(tmp_path / "x.txt", tmp_path / "weights")
