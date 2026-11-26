# pgc-risk

Give it a genome file, get back an estimated lifetime probability for 13 psychiatric
disorders, computed from the Psychiatric Genomics Consortium (PGC) GWAS summary statistics
published on Hugging Face.

```
> pgc-risk.bat my_genome.txt

Genome file  : my_genome.txt  (23andMe-style export, GRCh37)
Genotypes    : 503,109 autosomal SNPs with a call
Ancestry fit : 1.014 over 471,968 SNPs - matches the European reference (1.00 = typical, above 1.075 = not scorable)

Disorder                           Population      Your  Relative       Score    SNPs
                                     baseline  estimate   to avg.  percentile    used
-------------------------------------------------------------------------------------
Schizophrenia                            1.0%      3.4%     3.41x        98th  23,374
Bipolar disorder                         2.0%      0.8%     0.38x        13th  20,127
Major depression                          15%       14%     0.93x        45th  24,442
...
```

(The example is NA12878, a public 1000 Genomes sample.)

**Read the limits first.** These are statistical expectations for groups of people with a
similar score, not predictions for one person and not a diagnosis. Common-variant scores
capture about 7% of the liability to schizophrenia and bipolar disorder, 4% for depression
and 2% or less for the others, so most estimates stay close to the population baseline.
Family history, environment and rare variants are not included. The scores are valid for
**European ancestry only**; the tool checks this and refuses to score a genome that does
not fit. It is a research tool, not a medical test.

## Using it

**In a web browser, with nothing installed:** <https://bengoldfeder.com/tools/polygenic-risk>
(the page is `web/` in this repository). You download a reference file once from this
repository's releases (35 MB, or a 12 MB compact one) and load it together with your genome
file. Both are read by the page on your own device. Its content security policy forbids its
code from making any network request (`connect-src 'none'`), so a file cannot be sent
anywhere, and the page keeps working with the network switched off. The reference file is
kept in the browser's own storage so it is loaded only once; a genome file is never stored.
It reads genotyping-chip exports (23andMe, AncestryDNA, MyHeritage, FamilyTreeDNA,
LivingDNA). Sequencing VCFs need the full build below.

**From a downloaded reference file, on the command line:**
`pgc-risk <genome file> --bundle pgc-risk-reference-full-v1.pgcr` gives the same numbers as
the browser without the 20 GB build described under "Rebuilding from scratch".

The rest of this section is the full build, which clumps on exactly the SNPs each file has
and also reads sequencing files.

**With the graphical interface:** double-click `pgc-risk-gui.bat`. A console window opens
and your browser shows a page where you drop the genome file; results appear as a chart
and a table, and can be saved as JSON or printed to PDF. Close the console window to stop.

The page is served from your own computer (`http://127.0.0.1:8765/`) and nothing leaves
it: the browser hands the file to the local program, which scores it and deletes its
temporary copy as soon as scoring ends (a copy left by a program killed mid-run is removed
the next time it starts). The page loads no outside resources, the server listens on
127.0.0.1 only, refuses requests addressed to any other host name, and accepts uploads only
from the page it served. A checked run in a real browser made no request to anywhere else
(`tools/gui_browser_check.mjs`).

**From another device on the same network:** double-click `pgc-risk-gui-lan.bat` instead
(`pgc-risk gui --lan`). The console window and the page on this computer show the address to
open on the phone or other computer, for example `https://192.168.1.20:8765/`, and a six-digit
access code. Because the genome file now crosses the network:

* other devices connect over HTTPS. The certificate is made by this computer on first use
  (`data/gui-cert`), so the browser shows a warning once; the certificate's SHA-256
  fingerprint is displayed so it can be compared before continuing. Plain HTTP is answered
  only on this computer itself; any other device is redirected to HTTPS before anything is
  exchanged;
* other devices must enter the access code, which changes every time the program starts.
  Five wrong codes lock the form for a minute;
* the server listens on every network adapter but offers only the Wi-Fi or Ethernet
  address; it is still reachable only inside the local network, not from the internet.

Whether another device can actually connect also depends on things outside the program:
Windows Firewall must allow Python to accept connections, a VPN on this computer must
allow LAN connections, and the Wi-Fi must let its devices reach each other (guest and
public networks often do not). The plain `pgc-risk-gui.bat` never listens on the network.

**From the command line:**

```
pgc-risk.bat <genome file>                     drag-and-drop onto the .bat also works
                                               (in PowerShell: .\pgc-risk.bat <genome file>)
.venv\Scripts\pgc-risk <genome file> [options]

  -d scz,bip,mdd        only these disorders (keys: pgc-risk list)
  --json out.json       also write every number as JSON ("-" prints JSON only)
  --sample NAME         which sample of a multi-sample VCF
  --r2-scale 0.5        rerun with every score's assumed accuracy scaled, as a sensitivity check
```

Accepted genome files (plain, `.gz` or `.zip`):

| Format | Matched by | Notes |
| --- | --- | --- |
| 23andMe, LivingDNA | rs number, then GRCh37 position | about 4 s |
| AncestryDNA | same | |
| MyHeritage, FamilyTreeDNA (CSV) | same | |
| VCF, all sites (imputed or array-derived) | rs number, then position | GRCh37 or GRCh38 |
| VCF, variants only (whole-genome sequencing) | position | sites absent from the file are read as homozygous reference; about 30 s |
| gVCF | position | sites inside reference blocks are read as homozygous reference |

A VCF that does not name its genome build is matched against both GRCh37 and GRCh38
coordinates and the one that fits is used. Only autosomal SNPs are used; indels, the sex
chromosomes and strand-ambiguous (A/T, C/G) SNPs are ignored.

## What is inside

| Key | Disorder | PGC release used (European ancestry) | Cases / controls |
| --- | --- | --- | --- |
| scz | Schizophrenia | Trubetskoy 2022, PGC3 core GWAS (about 80% European, 20% East Asian) | 67,390 / 94,015 |
| bip | Bipolar disorder | O'Connell 2025, without 23andMe | |
| mdd | Major depression | Adams 2025, without 23andMe | |
| anx | Anxiety disorders | Strom 2026 | 122,083 / 729,602 |
| adhd | ADHD | Demontis 2023 | 38,691 / 186,843 |
| asd | Autism spectrum disorder | Grove 2019 | 18,381 / 27,969 |
| ptsd | PTSD | Nievergelt 2019, freeze 2 | 23,212 / 151,447 |
| ocd | Obsessive-compulsive disorder | Strom 2025, without 23andMe | 23,493 / 1,114,613 |
| ts | Tourette syndrome | Yu 2019 | 4,819 / 9,488 |
| an | Anorexia nervosa | Duncan 2017 | 3,495 / 10,982 |
| bpd | Borderline personality disorder | Streit 2026 | 12,339 / 1,041,717 |
| alc | Alcohol dependence | Walters 2018, unrelated genotyped subset | |
| cud | Cannabis use disorder | Johnson 2020 | |

Where the data comes from. The PGC summary statistics were put on Hugging Face as
`OpenMed/pgc-*`. Those repositories are no longer publicly reachable, so the tool downloads
from `introvoyz041/pgc-*`, a public mirror, pinned to a commit. A second mirror
(`tunahanf/pgc-*`) carries the same sha256 for all 92,801 files, which is the evidence that
the mirror is an unmodified copy. Anorexia nervosa is not in the mirrors; it comes from
`lighteternal/pgc-psychiatric-gwas-harmonized`. `pgc_risk/catalog.py` lists every repo,
commit, file range and column mapping.

Not everything on the Hub is usable. Releases that PGC distributes in its "sumstats VCF"
layout were converted to empty tables (bipolar 2021 full, PTSD 2024, panic disorder,
hoarding, the European-only schizophrenia file), so the newest readable release is used
instead. One of 737 depression shards and 72 of 915 alcohol-dependence shards are missing
upstream (about 0.1% and 8% of those SNPs).

## How a probability is computed

1. **Weights** (`pgc-risk build`). Per GWAS: keep autosomal, bi-allelic, unambiguous SNPs
   with an rs number, imputation quality at least 0.8, minor-allele frequency at least 1%
   and a near-complete sample size; attach one shared European allele frequency; keep the
   SNPs with p <= 0.05. Each SNP also gets its hg19 reference base and its GRCh38 position.
2. **LD panel**. The 503 European samples of 1000 Genomes phase 3, at the 4.0 million SNPs
   the weights use (`data/panel`, about 500 MB).
3. **Score** (`pgc-risk <file>`). Clumping and thresholding, which is the PGC's own recipe:
   among the SNPs the file actually has, walk from the most to the least significant, keep
   a SNP and drop every SNP within 250 kb that is correlated with it (r2 >= 0.1 in the
   panel). The MHC and two large inversions contribute one SNP each, and the other known
   long-range-LD regions use a 1 Mb reach. The score is the sum of log-odds times
   effect-allele count.
4. **Standardise**. z = (score - mean) / sd, where mean = sum 2pb from the allele
   frequencies and sd is the independent-SNP value sqrt(sum 2p(1-p)b^2) widened by an LD
   factor measured on the panel. Clumped SNPs are still weakly correlated, which makes the
   real spread 1.1 to 1.4 times wider than the textbook formula.
5. **Risk**. Liability-threshold model: with lifetime prevalence K and a score that
   explains a share R2 of liability variance,
   `risk(z) = 1 - Phi((T - sqrt(R2) z) / sqrt(1 - R2))`, `T = Phi^-1(1 - K)`.
6. **Ancestry check**. The mean of `(count - 2p)^2 / (2p(1-p))` over the file's SNPs is
   1.00 for a person from the reference population and grows with genetic distance. Above
   1.075 no estimate is given.

### The reference bundle (browser and `--bundle`)

The full method needs the 500 MB LD panel at scoring time. The bundle
(`pgc_risk/bundle.py`, built by `tools/build_bundle.py`) carries the same information in a
form that can be scored anywhere:

* Clumping is done once, at build time, over every candidate SNP. Each SNP is filed under
  the index SNP of its clump, most significant first (at most 48 per clump). A genome file
  is then represented, clump by clump, by the most significant SNP it actually has.
* Each SNP carries an LD factor: how much its share of the score's variance grows through
  correlation with the neighbouring clumps. When a file covers only part of the clumps the
  neighbour part shrinks with the covered share to the power 0.72, a value fitted on 1000
  Genomes across five very different SNP densities (rms error 1%).
* SNPs are stored strand-free and keyed by rs number, so no positions or reference genome
  are needed; that is also why the bundle cannot score sequencing files.

How close it is to the full method, on the 503 Europeans of 1000 Genomes: with every SNP
available the two agree exactly (correlation 1.000). With array SNPs only, the correlation
of the standardised scores is 0.96 to 0.97 (mean difference 0.2 sd), the spread of z is
1.02 on average, and the LD factor is within 2.5% of the one measured directly on the
selected SNPs. Hold-out accuracy with array SNPs, bundle rule against full method:
schizophrenia 6.75% against 6.91% and 5.35% against 5.43%, depression 2.24% against 2.31%,
OCD 0.94% against 1.00%.

Two bundles are published with each release: `full` (3.1 million SNPs, 35 MB, works with any
chip) and `chips` (1.0 million SNPs that sit on genotyping chips, 12 MB).

### The two numbers the risk depends on

| Key | Baseline K | R2 used | Where R2 comes from |
| --- | --- | --- | --- |
| scz | 1.0% | 7.3% | published, same method (Trubetskoy 2022); own hold-outs 6.7% and 5.3% at half the sample |
| bip | 2.0% | 6.7% | published 8.4% with PRS-CS (O'Connell 2025), times 0.8 |
| mdd | 15% | 3.9% | published, same method and same file (Adams 2025); own hold-out on UK Biobank alone 2.2% |
| anx | 20% | 1.8% | published 2.27% with MegaPRS (Strom 2026), times 0.8 |
| bpd | 1.5% | 2.3% | published 2.6% and 3.1% with PRS-CS in replication cohorts (Streit 2026), times 0.8 |
| ocd | 1.0% | 0.8% | own hold-out on UK Biobank (0.7% to 1.0%); no published figure |
| ptsd | 6.8% | 0.15% | published (Nievergelt 2019) |
| adhd | 5.0% | 2.1% | expected from GWAS size and heritability |
| asd | 1.2% | 1.0% | expected (the paper's 2.45% Nagelkerke R2 corresponds to about 1% of liability) |
| ts | 0.8% | 1.1% | expected |
| an | 1.0% | 0.7% | expected |
| cud | 6.3% | 0.9% | expected |
| alc | 12.5% | 0.1% | expected |

K is the lifetime prevalence each GWAS paper assumed, except PTSD, alcohol dependence and
cannabis use disorder, whose papers condition on exposure; those use US population surveys
(NCS-R, NESARC, NESARC-III). "Times 0.8" is an assumption: a clumping score typically
reaches about four fifths of what PRS-CS-type methods reach on the same data. "Expected"
uses `R2 = h2 x / (1 + x)`, `x = Neff h2_obs / M`, with M = 92,000 calibrated on the three
cases where a clumping result is known (they give 87,000, 89,000 and 99,000). The six
"expected" figures are the least certain numbers here; `--r2-scale` shows how much an
estimate moves if they are off. `tools/accuracy.py` regenerates the table and
`pgc_risk/accuracy.json` holds it with sources.

## What was checked

All of this is reproducible with the scripts in `tools/`; outputs are in `data/validation/`.

* **Effect direction** (`check_directions.py`). For schizophrenia, 98.5% of 21,773 strong
  associations have the effect allele more common in cases. Every other GWAS correlates
  positively with schizophrenia and with the rest across 26,546 evenly spaced SNPs, as the
  known genetic correlations require (depression-anxiety 0.57, schizophrenia-bipolar 0.39,
  ADHD-autism 0.31); a GWAS read with flipped alleles would show a row of negatives.
* **Does the score predict in people it was not trained on** (`loo.py`). Training on one
  PGC subset and testing on another, liability R2:

  | Trained on | Tested on | All SNPs | Array SNPs |
  | --- | --- | --- | --- |
  | Schizophrenia, men | women | 6.7% | 6.9% |
  | Schizophrenia, women | men | 5.3% | 5.4% |
  | Depression without UK Biobank | UK Biobank | 2.2% | 2.3% |
  | OCD without UK Biobank | UK Biobank | 0.7% | 1.0% |

  Array files do as well as sequenced genomes. An earlier version that pruned SNPs by
  distance alone reached about half of these values, which is why the LD panel exists.
* **Is z really standard normal in Europeans** (`calibrate.py`). LD was estimated on half
  of the 1000 Genomes Europeans and z evaluated on the other 251: sd 1.04 to 1.05 averaged
  over disorders (0.94 to 1.25 individually, sampling error 0.04), mean +0.10 to +0.14.
  Without the LD factor the sd would be 1.1 to 1.4.
* **File handling** (`make_test_genomes.py`, `check_e2e.py`). Eight 1000 Genomes people
  written as a 23andMe-style file, a GRCh37 variants-only VCF and a GRCh38 VCF without
  build information: the first two reproduce the direct matrix computation exactly, the
  third to within 0.08 sd (a few SNPs cannot be lifted to GRCh38).
* **Ancestry check** on all 2,504 people of 1000 Genomes, array SNPs: Europeans median
  1.011, 95% below 1.037, maximum 1.070; South Asians 1.11 or more, East Asians 1.37 or
  more, Africans 1.15 or more. Their scores sit up to 12 sd from the European centre, which
  is why nothing is reported for them.
* **Graphical interface** (`gui_browser_check.mjs`). A headless browser opened the served
  page, chose a file through the page's own file input and waited for the result, for an
  array file, a 54 MB sequencing VCF, a non-European genome and a file that is not a
  genome: same numbers as the command line, the right message in the last two cases, no
  request to any address but the local server, and no browser error or blocked resource.
  Network mode was driven the same way through this computer's own Wi-Fi address (HTTPS,
  wrong then right access code, array file and sequencing VCF). It has not been exercised
  from a second physical device.
* **Bundle and browser version.** `check_bundle.py` and `loo_bundle.py` produce the bundle
  figures quoted above. `check_web.py` runs the browser's scorer under Node next to the
  Python one on 13 files (plain, AncestryDNA-style, CSV, zipped, gzipped, with duplicated
  and missing calls): every number agrees to 12 decimal places. `web_browser_check.mjs`
  opens the page in a headless browser, takes it offline, loads the reference file and a
  genome, and compares the page's table with the scorer run directly; it fails on any
  network request after the files are chosen, any request to another host, any browser
  error or blocked resource.
* 50 unit tests (`pytest`), including reference-base and liftover lookups for six
  well-known SNPs, the bundle format and scoring rule, and the interface's upload,
  clean-up, HTTPS, access-code and host checks.

## Known limits

* **Within Europe the centre moves.** Average z by 1000 Genomes population and disorder
  ranges from -1.0 (Finns, bipolar) to +0.9 (Utah residents, PTSD); about four in five are
  within 0.4. The tool does not correct for this, so a percentile can be off by that much
  for regional reasons.
* **Mostly-European admixed genomes can pass the ancestry check** (some Puerto Rican and
  Colombian 1000 Genomes samples do) while their scores are already shifted by up to about
  1 sd. Above a fit of 1.04 the report says to treat the numbers as rough.
* **Low-powered GWAS** (marked `*` in the output: PTSD, OCD, anorexia, alcohol, cannabis)
  barely move the estimate.
* The p <= 0.05 cut-off is fixed at build time. It is the best of the thresholds tested for
  schizophrenia and depression; weaker GWAS might do slightly better with a looser one.
* The browser version and `--bundle` approximate the full method (see "The reference
  bundle") and have not been run on a real consumer export either: the test files are 1000
  Genomes samples written in those formats.
* Sequencing VCFs: a site missing from a variants-only file is assumed to be homozygous
  reference. That is right for a 30x genome and wrong for low-coverage or exome data.
* Prevalences are single population-wide figures; there is no adjustment for sex or age.

## Rebuilding from scratch

```
uv venv .venv                    any Python 3.11 or newer
uv pip install --python .venv\Scripts\python.exe -e .[dev]
.venv\Scripts\pgc-risk fetch     about 9,000 parquet shards, 4.9 GB, from Hugging Face
.venv\Scripts\pgc-risk build     weights, then the LD panel (downloads 1000 Genomes, 15 GB)
```

`python tools/build_bundle.py` then writes the two reference bundles to `data/bundle`.

`build` also fetches `hg19.2bit` (816 MB) and the hg19-to-hg38 chain from UCSC, which
sequencing VCFs need. With everything already downloaded a rebuild takes about a minute.
After a build only `data/weights` (280 MB) and `data/panel` (500 MB) are needed to score;
`data/raw` (7 GB), `data/reference` (26 GB) and `data/validation` can be deleted.

```
pgc_risk/   catalog.py   which GWAS, where, how to read it, prevalences
            fetch.py     Hugging Face download
            build.py     QC and weight files
            panel.py     1000 Genomes LD panel
            reference.py hg19 bases and hg19 -> hg38 liftover
            genome.py    genome-file parsers
            score.py     matching, clumping, standardising, risk model
            report.py    text and JSON output
            bundle.py    the portable reference bundle: format, reader, scoring rule
            gui.py       local web server for the graphical interface
            gui.html     the page it serves (self-contained)
            gui_login.html  the access-code page other devices see first
web/        the browser version (pgcrisk-core.js is the scorer; tool.css, app.js, worker.js the page)
tools/      validation scripts described above, and build_bundle.py
tests/      unit tests
```

## Credit and terms

Summary statistics: Psychiatric Genomics Consortium, https://pgc.unc.edu/, CC BY 4.0. Using
them means citing the original paper for each GWAS (listed by `pgc-risk list` and in
`pgc_risk/catalog.py`) and not attempting to re-identify participants. LD panel: 1000
Genomes Project phase 3. Reference sequence and liftover chain: UCSC Genome Browser.
