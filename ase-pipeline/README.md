# `/ase-pipeline` — Allele-Specific Expression Skill (Stages 1 to 3)

An interactive Claude Code skill that walks you through a complete allele-specific expression (ASE) analysis of bulk RNA-seq on an HPC cluster running SLURM and Singularity: reference and genotype preparation, allele-aware alignment (STAR), allele counting (GATK ASEReadCounter), and downstream statistics as R Markdown reports plus a standalone summary page. Two designs share one wizard: **F1 crosses** (for example B6 x AJ mice) and **outbred / human** samples with per-individual genotypes.

Stages 1 to 3 cover per-sample allelic imbalance, the reciprocal F1 analysis (strain versus parent-of-origin effects), differential ASE between conditions and, in outbred mode, phASER read-backed phasing with gene-level haplotype counts (reported next to the unphased results, never instead of them).

**What Stage 2 adds**

- **Reciprocal F1 (Rmd 03)**: strain and parent-of-origin effects from crosses in both directions.
- **Differential ASE (Rmd 04)**: changes in allelic imbalance between conditions, for F1 crosses (GLM per gene and per SNP) and for outbred individuals (paired, direction-free test).
- **Thinned Rmd 02 F1 gene test**: the gene-level likelihood-ratio test of Stage 1 now uses each gene's SNPs thinned to one SNP per `THIN_BP` window, so that no read pair is counted twice.

**New requirements of Stage 2**

- `samples.csv` column `cross_direction` must be exactly `AxB` or `BxA` for the reciprocal analysis (A = the strain carrying the reference allele, the mother in `AxB`); each direction needs at least 2 samples.
- Replicates mean at least 2 samples in a condition (F1) or at least 2 individuals sampled in both conditions of a contrast (outbred); a contrast without them is not tested and is listed as such.
- Differential ASE asks for a reference condition (the baseline that every other condition is compared with). In F1 mode the condition must not be confounded with the cross direction (Rmd 04 stops when it is).

**What Stage 3 adds**

- **phASER (Rmd 05)**: a per-sample array runs phASER (repository secastel/phaser, pinned to commit `aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301`) on the WASP-filtered, duplicate-marked BAM and the same heterozygous genotype VCF that STAR and ASEReadCounter use. phASER links heterozygous SNPs into haplotype blocks from read pairs that span them, and its gene step (`phaser_gene_ae`, over the gene spans in `reference/genes_span.bed`) gives gene-level counts for haplotype A and B. Rmd 05 tests these counts per sample and gene (exact beta-binomial test against 0.5, symmetric in A and B) and compares them gene by gene with the unphased results of Rmd 02. It is reported next to Rmd 02, never instead of it: the unphased per-SNP tests and the ACAT gene test stay the main per-sample result.
- **One-time installation job** (`setup_phaser_env.sh`, a compute job of 10 to 20 minutes; it needs network access on a compute node for `git clone` and conda): a conda environment with pinned package versions under `{PHASER_HOME}`, built from the pinned commit and never the bioconda package `phaser` (which is Python 2.7). Every script exports `PYTHONNOUSERSITE=1` and puts the environment's `bin` first on `PATH`; no conda environment is activated by the scripts. The chain submits the installation only when the installed commit or package list differs from the pinned ones.
- **Wizard questions**: menu option 4 (analyses word `phaser`), whether the genotype VCFs are phased (with the unphased caveat told to the user), and where to install phASER.
- **Chain and summary page**: the phASER array runs after the alignment array, Rmd 05 after the phASER array and Rmd 02 (`afterok` on both), and the standalone summary page gets a phASER section (table, direction statement, limits, links).

**New requirements of Stage 3**

- Outbred mode and paired-end data only: phASER runs with `--paired_end 1`, the only setting tested. It is never offered in F1 mode.
- Contig names without `_` (see the limitations), and the same contig names in the GTF, the FASTA and the VCFs.
- About 650 MB for the installation directory (about 2.6 GB free while it is built), reused by later projects.

---

## Installation

```bash
cp ase-pipeline.md ~/.claude/commands/
```

Then invoke it in Claude Code:

```
/ase-pipeline
```

---

## Prerequisites

| Requirement | Notes |
|-------------|-------|
| SLURM scheduler | Every job is submitted with `sbatch -p bcc`; nothing heavy runs on the login node |
| Singularity >= 3.10 | Loaded via `module add singularity/3.10.4` for every tool; the only other module is `module add miniconda3/v4`, in the one-time phASER installation job |
| Cached biocontainers | In `$NXF_SINGULARITY_CACHEDIR` (or `~/.singularity/cache`): STAR 2.7.10b, GATK 4.4.0.0, bcftools 1.20, samtools 1.21, Picard 3.1.1 (`depot.galaxyproject.org-singularity-*.img`); a missing image is downloaded by the prep job (it fetches all five, Picard included, so the per-sample array job only checks for them), never in the foreground on the login node |
| `bulkrnaseq` image | `/net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif`, the only place R is available; needs `openxlsx`, `tidyverse`, `GenomicRanges`, `rtracklayer` (checked by a tiny `sbatch` job in Step 2) |
| Genotypes | F1: a parental-difference VCF (plain `.vcf` or bgzipped; sites-only is fine, strain names default to `C57BL_6NJ` / `A_J` when the header has none), or the built-in Mouse Genomes Project helper (needs internet on a compute node). Outbred: one VCF per individual |
| Internet access | Only for downloading missing containers, references and the Mouse Genomes Project data |
| phASER (optional, outbred) | Network access on a compute node for `git clone` and conda during the one-time installation; about 650 MB under `{PHASER_HOME}` |

---

## What the skill does

| Step | Topic |
|------|-------|
| 0 | Working directory, `{WD_NAME}`, date tag, results directory (auto) |
| 1 | Email address for SLURM notifications |
| 2 | Environment check: Singularity, cached containers, R packages in the `bulkrnaseq` image (no R on the login node) |
| 3 | Mode: F1 cross or outbred / human |
| 4 | Reference: standard Ensembl-style folder or custom FASTA and GTF |
| 5 | Sample sheet (`samples.csv`), name sanitisation, read length and STAR index per read length and mode |
| 6 | Genotype source: parental-difference VCF or mouse helper (F1); per-individual VCFs (outbred) |
| 7 | Analysis menu: per-sample (always), reciprocal F1, differential ASE (with reference condition), phASER (outbred, paired-end; phased-genotypes and install-location questions) |
| 8 | Constants: `MIN_DEPTH`, `FDR_SIG`, `ABS_DEV_SIG`, `BIAS_TOL`, `RHO_MIN` |
| 9 | Resources for the array and prep jobs |
| 10 | Reference and genotype preparation scripts (masked reference or per-individual het VCFs, STAR index, optional mouse helper) |
| 11 | Per-sample array scripts: STAR, duplicate marking, ASEReadCounter, empty-table guard |
| 12 | Rmd 01: import, site filters, reference-bias QC, WASP removal report (outbred), checkpoint |
| 13 | Rmd 02: per-SNP and per-gene tests, figures, xlsx, `summary_numbers.tsv` |
| 14 | Statistics functions (shipped code, unit-tested) |
| 15 | Submission order with dependencies, the standalone summary page, notes for the assistant |
| 16 | Rmd 03: reciprocal F1 |
| 17 | Rmd 04: differential ASE |
| 18 | phASER: installation and per-sample haplotype counts |
| 19 | Rmd 05: phASER gene-level haplotype imbalance |

---

## Output files

Everything is written under `{CWD}/results/{YYYY-MM-DD}_{WD_NAME}/` (for example `results/2026-09-29_proj/`) (raw FASTQ, BAM and VCF inputs are never modified). The STAR index goes to `{GENOME_DIR}/index/` so it can be reused across projects: the F1 masked index is named after the strain pair, the read length and a checksum of the masked sites (`star_ase_masked_{STRAIN_A}_{STRAIN_B}_sjdb{N}_{MASK_KEY}`), so a different strain pair or a corrected parental VCF never reuses a stale mask:

| Output | Description |
|--------|-------------|
| `scripts/` | Generated `sbatch` scripts: prep (`prep_f1_reference.sh` or `prep_genotypes.sh`), optional mouse helper, per-sample array (`align_count_f1.sh` or `align_wasp_count.sh`), `run_01_import_qc.sh`, `run_02_imbalance.sh`, and when selected `run_03_reciprocal.sh`, `run_04_differential.sh` |
| `bam/` | Coordinate-sorted BAMs with read groups, duplicates marked (not removed), plus indexes |
| `ase_counts/` | ASEReadCounter table per sample; outbred mode also `{sample}.wasp_stats.tsv` and `{sample}.unfiltered.table` (the same count before the WASP filter, for the before/after REF fraction in Rmd 01, reported for information, not a gate; Rmd 01 reads every table by sample name, never by glob) |
| `ase_checkpoint.rds`, `ase_imbalance_checkpoint.rds`, `ase_reciprocal_checkpoint.rds`, `ase_differential_checkpoint.rds` | Checkpoints of Rmd 01, 02, 03 and 04 (the last two only when run) |
| `{date}_{project}_01_import_qc.html`, `..._02_imbalance.html`, `..._03_reciprocal.html`, `..._04_differential.html` | The stage reports (03 and 04 only when selected) |
| `{date}_{project}_ASE_imbalance.xlsx` | Sheets `SNP`, `Gene`, `Summary` |
| `{date}_{project}_ASE_reciprocal.xlsx` | Sheets `Gene`, `SNPs_used`, `Design`, `Excluded`, `Summary` |
| `{date}_{project}_ASE_differential.xlsx` | Sheets `SNP`, `Gene`, `Design`, `Excluded`, `Summary` |
| `..._ASE_sites_vs_depth.pdf`, `..._ASE_genes.pdf`, `..._ASE_reciprocal.pdf`, `..._ASE_differential.pdf` | Figures |
| `summary_numbers.tsv` | One row per sample with the headline numbers (written by Rmd 02) |
| `summary_numbers_reciprocal.tsv`, `summary_numbers_differential.tsv` | Headline numbers of Rmd 03 and Rmd 04 |
| `reference/genes_span.bed` | phASER only: one row per gene (0-based start and end over all its exons, introns included; written at the end of `prep_genotypes.sh`) |
| `scripts/setup_phaser_env.sh`, `phaser_count.sh`, `run_05_phaser.sh` | phASER only: the installation job, the per-sample phASER array and the Rmd 05 run script |
| `phaser/{sample}.haplotypic_counts.txt` and `.gene_ae.txt` | phASER only: haplotype blocks and the gene-level counts (`aCount`, `bCount`, `totalCount`, `n_variants`, `gw_phased`); the other phASER files of the sample are kept beside them (not used by the Rmds) |
| `logs/phaser_{sample}.log` | phASER only: phASER and `phaser_gene_ae` output of one sample |
| `{date}_{project}_05_phaser.html` | The Rmd 05 report (only when phASER was selected) |
| `{date}_{project}_ASE_phaser.xlsx` | Sheets `Gene`, `Comparison`, `Samples`, `Summary` |
| `{date}_{project}_ASE_phaser.pdf` | Figures of Rmd 05 |
| `{date}_{project}_ASE_phaser_gene.tsv.gz`, `..._ASE_phaser_comparison.tsv.gz` | The full `Gene` and `Comparison` tables, tab-separated and gzip-compressed. The xlsx sheets hold the full tables only while they fit in a worksheet (1,048,575 data rows, Excel's limit); beyond that they hold only the significant rows and Rmd 05 prints a note that names these files. The `SNP` and `Gene` sheets of Rmd 02 and Rmd 04 have the same limit; their checkpoints hold the full tables |
| `ase_phaser_checkpoint.rds`, `summary_numbers_phaser.tsv` | Checkpoint and headline numbers of Rmd 05 |
| `{WD_NAME}_summary_report.html` | Standalone page (inline CSS, no R, relative links to the reports and the xlsx, reference-bias flags shown prominently) |

The jobs run in one dependency chain: prep -> per-sample array (`afterok`) -> Rmd 01 (`afterok`) -> Rmd 02 (`afterok`), and, when selected, Rmd 03 / Rmd 04 (`afterok` on Rmd 01, beside Rmd 02), preceded by the mouse helper chain when used. A wait script follows every submitted job and stops when only jobs that can never start remain. With phASER, the chain also holds the installation job (only when the installed commit or package list differs from the pinned ones) and the phASER array after the alignment array, then Rmd 05 (`afterok` on the phASER array and Rmd 02).

---

## Key design points

- **F1 mapping bias: third-allele masked reference.** At every parental SNP the genome base is replaced with a third allele (neither strain's), so neither strain is favoured in mapping; STAR runs on the masked genome and ASEReadCounter counts against the het-sites VCF.
- **Outbred mapping bias: WASP.** STAR's WASP re-mapping filter tags every alignment (`vW`); alignments with `vW` 2-7 are removed before counting, and Rmd 01 reports how many were removed and from which allele. WASP does not promise less bias, and none is claimed.
- **Per-sample beta-binomial tests.** Each sample gets its own overdispersion. F1 uses a free-mean-per-gene estimate with a bias correction (a plain H0-based estimate is inflated by real imbalance), floored at `RHO_MIN`, for both SNP-level exact tests and a gene-level likelihood-ratio test on the SNP-level counts (counts are not summed); the gene test and its own rho estimate use each gene's SNPs thinned to one SNP per `THIN_BP` window (exon coordinates), so no read pair is counted twice. Outbred (unphased, so no per-gene mean can be fitted) uses a trimmed estimate: the H0 beta-binomial fitted on the central sites only (each site's central 90 percent region), with a truncated likelihood that corrects for the trimming, floored at `RHO_MIN`. The naive all-sites H0 estimate is shown only for comparison (`rho_h0_naive`): real imbalance inflates it, and on the synthetic acceptance data it detected 0 of 16 planted outbred genes.
- **Reciprocal F1 (Rmd 03).** Per gene, `logit(p) = b0 + b1 * d`, where p is the strain-A fraction and d = +1 for `AxB`, -1 for `BxA`. b0 is the strain effect and b1 the parent-of-origin effect. Rows are the samples, with the counts summed over the gene's thinned SNPs. The dispersion is estimated between replicate animals, and each gene uses the larger of the pooled value, its own estimate (with at least 2 residual df) and `RHO_MIN`. Both effects get likelihood-ratio tests. X, Y and MT are excluded. b1 is not affected by a constant mapping bias toward one strain (it cancels between the two directions); b0 is, so the reference-bias flags are repeated next to the strain results.
- **Differential ASE (Rmd 04).** F1: a GLM per gene and per SNP with the condition term and, when both directions are present, the cross-direction term; the Rmd stops when condition and direction are confounded. Outbred: a paired, direction-free per-SNP test summed over the individuals sampled in both conditions (the allele carrying a regulatory variant differs between individuals, so no direction is reported), genes combined with ACAT. No mixed model (lme4) is used: a random individual effect is not estimable with 2-6 individuals, and a shared condition slope would cancel phase-dependent changes.
- **Base-R models.** The GLM, the pair dispersion and the tests are base-R functions shipped in Step 14 and extracted by the unit tests, so the tested code is the shipped code; Rmd 03 and 04 need no package beyond those of Rmd 01 and 02.
- **Unphased outbred data.** Without phasing, SNPs of a gene cannot be pooled: per-SNP tests are combined per gene with ACAT and labelled "unphased, no direction".
- **phASER gene test (Rmd 05).** phASER counts each read pair once per haplotype block, so a gene's counts pool its heterozygous SNPs without counting a fragment twice, and genes whose SNPs each have fewer than `MIN_DEPTH` reads can still be tested. The test is an exact beta-binomial test of `aCount` out of `totalCount` against 0.5 per sample and gene, BH within each sample. The overdispersion is the larger of the sample's own trimmed central-region fit, the cohort value (the median of the samples' own dispersions, so adding or removing samples can change the p-values of the others) and `RHO_MIN`. It is symmetric in A and B and gives no direction with unphased genotypes (the haplotype labels are arbitrary); with phased genotypes the direction is given only for genome-wide phased genes. The gene spans include introns, while Rmd 02 uses exonic SNPs only. There is no differential ASE on phASER counts.
- **Honest reference-bias diagnostic.** For every sample the mean reference-allele fraction over the filtered sites is compared with 0.5, and the sample is flagged when the deviation exceeds `max(BIAS_TOL, 3 x SE)`. It is a screen, not a test; flags are shown next to every ratio table and figure and at the top of the summary page.
- **Duplicates are marked, not removed.** ASEReadCounter skips duplicate-flagged reads itself.
- **Read group and indexed het VCF are mandatory.** Without a read group, or with a VCF that is not bgzipped, indexed and genotyped, ASEReadCounter exits 0 with an empty table; both requirements were demonstrated by the verification gate, and the scripts fail on an empty table.
- **Containers only, no R on the login node.** R and all tools run through `sbatch -p bcc` with Singularity; logs and temporary files go to the shared filesystem, never `/tmp`.

---

## Validation status

What has been exercised:

- **Verification gate on synthetic data** (STAR, GATK, samtools, bcftools in the cached containers): read group, indexed VCF, heterozygous genotype column and reference requirements of ASEReadCounter, STAR WASP behaviour, the Mouse Genomes Project access pattern. Results: `tests/fixtures/verification.md`.
- **Script smoke tests**: one F1 sample and one outbred sample run through the generated scripts.
- **Unit-tested statistics**: null size of the F1 gene-level test 0.047-0.056 pooled (worst seed 0.064) and power 0.99, simulated with independent SNPs; with SNPs that share read pairs (fragment-level simulation, `tests/r/test_ase_stats_stage2.R`) the unthinned F1 gene test had a size of 0.1450 and the thinned test that Rmd 02 now uses 0.0497 (see the thinning limitation below); outbred path (trimmed rho, per-SNP test, ACAT per gene) on unphased simulated data: null gene size 0.003-0.053 pooled across six scenarios (worst seed 0.083), power 1.00 / 0.88 for planted 0.85 / 0.70-0.30 genes on binomial data where the naive H0 rho gave 0; plus rho recovery, edge cases and ACAT (`tests/r/test_ase_stats.R`; needs a job with at least 8 CPUs).
- **Rmd rendering tests** on perturbed copies of one synthetic sample (F1 and outbred, Rmd 01 and Rmd 02), including the `summary_numbers.tsv` output.
- **End-to-end synthetic acceptance run** (2026-09-29, DONE) in both modes: the scripts and Rmds were generated from the installed skill's templates by extracting its code blocks and substituting the placeholders, then submitted following its steps; the interactive dialogue of Steps 0-9 (FASTQ scan and pairing, name sanitisation, sample-sheet scaffolding, menus) was not exercised as a user would drive it. The whole Step 15 chain ran (prep -> per-sample array -> Rmd 01 -> Rmd 02 with `afterok`); every job exited 0.
  - F1: 24 of 24 planted gene x sample tests significant with the correct direction, 0 of 42 null false positives, with the Stage 1 gene test. With the thinned F1 gene test (one SNP per `THIN_BP` window), the same count tables give 21 of 24 planted gene x sample tests and 0 of 42 null false positives, with the per-SNP results unchanged (65 significant SNPs); the three lost tests are genes thinned to one SNP.
  - Outbred: the WASP-filtered null-gene REF fraction was within `max(0.03, 3 x SE)` of 0.5 in 4 of 4 samples; planted-gene detection 10 of 16 with 0 of 28 null false positives (unphased, so direction is not applicable).
  - The unfiltered (pre-WASP) outbred null REF fraction is reported for information, not as a gate: 2 of 4 samples exceed the tolerance before WASP (0.547, 0.546), and WASP removes the excess.
  - Step 15 summary pages were written for both projects; the href-existence and no-`http` checks passed.
  - History: the first acceptance run detected 0 of 16 planted outbred genes with the naive H0 estimator, which led to the robust estimator (`bb_estimate_rho_trim`).

- **Input guards added after the acceptance run** (2026-09-29), run from scripts generated from the skill's templates on the synthetic data, each through `sbatch`:
  - F1 parental VCF with two strain columns: only `0/0` + `1/1` sites kept (49 of 55; same-genotype, heterozygous and missing sites dropped and counted).
  - Stops with a non-zero exit on: sites where strain A is `1/1` and strain B `0/0`; a VCF where both strains carry ALT (0 sites left); unknown strain names; `1` vs `chr1` contigs (F1 and outbred, before any alignment).
  - The masked STAR index is keyed on the mask: the same sites reused the index (also from a two-sample VCF), and a corrected VCF built a new one.
  - Mouse helper: one 200 kb region of chromosome 19 ran through the helper and concat jobs while the genome folder did not exist yet; the prep job then compared the contigs and stopped on the mismatch with the synthetic `chr1` genome. A `chr19` chromosome line stops the helper.
  - Rmd 01 and 02 render on count tables that mix contigs `1`, `X` and `MT` (the previous Rmd 01 stopped there).
  - The masking step now uses a genotyped mask VCF (`bcftools consensus -s MASK`): with the earlier sites-only VCF and `-H A`, bcftools 1.20 applied the sites in only 4 of 10 identical runs, which the masked-genome check caught as a failed prep job.

Limits of the outbred result: per-sample power is 10 of 16. The `RHO_MIN` floor caps detection at about 12 of 16 on these data, and two samples have only 27-29 sites. With real overdispersion and many imbalanced sites the robust dispersion estimate stays too high, and power for moderate imbalance drops.

### Stage 2 (synthetic data only)

- **Unit tests** (`tests/r/test_ase_stats_stage2.R` together with `tests/r/test_ase_stats.R`; the job needs at least 8 CPUs): final run 173 checks ok, 0 failed, exit 0. Simulation sizes are over 10 seeds; true dispersion 0.02 unless stated.
  - Reciprocal F1 (3+3 animals, 2000 null genes per seed): null sizes 0.0370 (strain) / 0.0367 (parent of origin); 0.0362 / 0.0346 (2+2), 0.0356 / 0.0341 (3+2), 0.0380 / 0.0365 (4+2). Power 1.000 for a strain effect of 0.7 and 1.000 for a maternal effect of 0.8 (3+3); signs correct in every call. The tests are conservative (about 0.036 at a nominal 0.05).
  - Differential F1 (null genes): 0.0369 (3 vs 3), 0.0363 (2 vs 2), 0.0365 (3 vs 2), 0.0372 (unbalanced directions with the cross-direction term). Power for a change from 0.5 to 0.7: 0.7813 (3 vs 3), 0.6853 (3 vs 2), 0.5983 (2 vs 2).
  - Paired outbred differential test (3000 SNPs per seed). Sizes are over the SNPs that could be tested (at least 2 informative individuals; 60 percent of individuals heterozygous per SNP): tested fraction 0.4244 / 0.6961 / 0.8434 / 0.9646 for 2 / 3 / 4 / 6 individuals. SNP / gene null sizes at dispersion 0.02: 0.0485 / 0.0490 (2 individuals), 0.0511 / 0.0505 (3), 0.0515 / 0.0526 (4), 0.0490 / 0.0514 (6). Power (4 individuals, 10 percent of SNPs changed from 0.5 to 0.8) is 0.8653 (phase-heterogeneous) and 0.9174 (consistent). The gate is restated from an absolute 0.8 to at least 0.70 and at least 90 percent of the oracle power with the dispersion fixed at its true value (0.9022 and 0.9424); the original 0.8 gate was not met (0.626 / 0.665 with the plain moment dispersion; 0.756 / 0.798 even with the true dispersion, in the first version of that simulation).
  - Fragment-level simulation (read pairs shared between SNPs): the unthinned Stage 1 F1 gene test had a size of 0.1450 (0.2084 in SNP-dense genes), the thinned test that Rmd 02 now uses 0.0497; thinned reciprocal tests 0.0314 (strain) / 0.0288 (parent of origin), 0.0331 in SNP-dense genes.
  - Runtime check: 5000 units x 12 rows, 2 tests, one core, 32 s; 60,000 units project to 0.11 h (projected, not run at genome scale).
- **Rmd 03 on synthetic counts** (60 genes), compared with the simulation truth by `tests/synthetic/evaluate_stage2.R`: 5 of 5 planted strain genes and 5 of 5 planted parent-of-origin genes found with the correct sign, 0 of 47 null genes called for either effect, exit 0.
- **Rmd 04 on synthetic counts**: F1, 4 of 4 planted genes found with the correct sign, 0 of 9 strain or parent-of-origin genes called as condition effects, 0 of 47 null genes called. Outbred (paired), 5 of 5 planted genes found, 0 of 3 imbalanced-but-unchanged genes and 0 of 52 null genes called. On these data every outbred pair dispersion sits at the `RHO_MIN` floor, so a dispersion above the floor is exercised only by the unit tests.
- **Negative tests**: a confounded F1 design (condition = cross direction) stops with a message; an unpaired outbred individual is listed and ignored; a contrast with nothing testable in a three-condition design is reported with `tested` = 0; a gene moved to chromosome X is excluded from Rmd 03 and listed in its `Excluded` sheet.
- **Rmd 02 with the thinned gene test**: see the Stage 1 acceptance data above (21 of 24 planted, 0 of 42 null, 65 significant SNPs unchanged).
- **Submission chain**: the run scripts, `submit_chain.sh` and `wait_chain.sh` were run against a stub scheduler in 10 scenarios, and the rendered Rmds on the synthetic data completed with exit 0.

Stage 2 synthetic acceptance run: DONE (2026-09-30)

Measured on the synthetic Stage 2 data, with the scripts and Rmds generated from the skill text and no hand edits to any Rmd:

- **F1 reciprocal**: 5 of 5 planted strain genes and 5 of 5 planted parent-of-origin genes found with the correct sign; 0 of 4 parent-of-origin calls among strain-only genes and 0 of 4 strain calls among parent-of-origin-only genes; 0 of 47 null genes called in each test (0/47 false positives per test).
- **F1 differential**: 4 of 4 planted genes found with the correct sign, 0 of 9 calls among strain / parent-of-origin genes, 0 of 47 null genes.
- **Outbred differential**: 5 of 5 planted genes found, 0 of 3 imbalanced-but-unchanged genes, 0 of 52 null genes; the direction labels were right ("REF lower in treat (all tested individuals)", "mixed (phase differs)" on the phase-heterogeneous SNPs).
- **Build and outputs**: all Rmds built without hand edits; 28 of 28 count tables non-empty; both summary pages passed the href and no-`http` checks.
- **Chain**: a manually adapted chain (the previous text of `submit_chain.sh` with hand renames and deleted lines, same dependencies) and `wait_chain.sh` ran prep, array, Rmd 01, Rmd 02 and Rmd 03 / 04 to COMPLETED (66 s in F1 mode, 77 s in outbred mode from submission to the last Rmd). The shipped mode-neutral `submit_chain.sh` (mode variables, conditional Rmd 03 / 04, word guard, refusal to overwrite `chain_job_ids.tsv`) was verified with a stub scheduler, not live (the Stage 3 acceptance run later ran it live). The times are short only because the synthetic genome is tiny; they say nothing about real data. `sacct` was unavailable, so `wait_chain.sh` inferred success from the result files.
- **Rmd 02 per sample on the aligned F1 data**: planted strain genes detected in 50 of 60 gene x sample tests, all in the expected direction; other planted imbalances 70 of 72; null false positives 1 of 564.
- **Aligned versus direct counts**: the aligned counts were about 95 percent (F1) and 94 percent (outbred, after WASP) of the direct reads, and the differential gene calls were the same. The F1 strain calls differ by one gene (SYN000047, a planted differential gene that is also allowed a strain call), significant only in the aligned data; the cause was not isolated.

### Stage 3 (synthetic data, and a real-data smoke test of the phASER part)

- **Unit tests** (`tests/r/test_ase_stats_stage3.R` together with the Stage 1 and Stage 2 tests; the job needs at least 8 CPUs): final run 213 checks ok, 0 failed, exit 0 (62 + 114 + 37). `hap_gene_test` was simulated with cohorts of 6 samples x 2000 genes (null size is the pooled value over the seeds; the oracle uses the true dispersion):
  - True dispersion 0.02: null size 0.0454 with no imbalanced genes, 0.0357 with 20 percent and 0.0188 with 40 percent of the genes imbalanced; power 0.8560 against an oracle of 0.8751 with 20 percent imbalanced. The test is conservative.
  - **Power is lost when many genes are imbalanced**, because real imbalance inflates the cohort dispersion (conservative, but it costs power; the estimator is the same as in Rmd 02): at a true dispersion of 0.02 and 40 percent imbalanced genes the power is 0.7843 against an oracle of 0.8748; at 0.05 and 20 percent it is 0.4882 against an oracle of 0.5957; at 0.05 and 40 percent the cohort dispersion is 0.1324 (true 0.05), the null size 0.0008 and the power 0.1460 against an oracle of 0.5961. No false positives are added, but real effects are missed.
  - **Power is lower with few genes per sample**: 0.6330 / 0.6757 / 0.7789 with 20 / 30 / 60 genes per sample (oracle about 0.877; null size 0.0233 / 0.0210 / 0.0265), because a small sample's dispersion estimate is noisy and the larger value is used.
  - Fragment-level simulation (read pairs shared between the SNPs of a gene; 6 samples x 1500 genes, 30 percent imbalanced; the number of fragments per gene is Poisson with mean 300 or 40): at a mean of 300 the haplotype test has a size of 0.0184 against 0.0068 for the unphased ACAT test, and in genes with at least 4 SNPs a power of 0.7877 against 0.5873. At a mean of 40 the size is 0.0042 and the power 0.1495 against 0.0357; 8-SNP genes are tested in 0.9962 of the cases against 0.4169 for the unphased path. **This gain is optimistic for unphased genotypes**: the simulation puts all SNPs of a gene into one block, which is what the haplotype test gets only with phased genotypes (see the limitations).
  - Runtime check: one sample of 20,000 genes, one core, 55 s and 2027 MB; 100 such samples project to 1.54 h on one core (projected, not run).
- **Synthetic data** (seed 20261001; 6 samples of 80 genes: 4 `hap_strong`, 4 `hap_moderate`, 3 `hap_lowdepth`, 2 `two_block`, 6 `null_linked` and 61 null; truth tables for the phase, the fragments and every sample-gene cell): the checks pass, 25 injected defects each fire the check named for them (swapped phase, swapped read counts, a dropped read pair, duplicated rows, a wrong contig and others), the read-level check confirms the simulated reads follow the truth haplotypes, and a regeneration gives the same file checksums (the sample sheet, which holds absolute paths, is excluded).
- **phASER installation** (compute job, `sbatch -p bcc`): the pinned solve finished in 631 s (an unpinned solve of the upstream environment file took 1252 s) and the reinstall after the package-list hash was added in 669 s; python 3.14.7, numpy 2.5.3, scipy 1.18.1, pandas 3.0.6, pysam 0.24.1, intervaltree 3.2.1, samtools and bcftools and htslib 1.24, bedtools 2.31.1; a second run of the script ends with "already installed" in under 1 s, and a non-empty directory without the ownership marker is refused. Without `PYTHONNOUSERSITE=1` importing numpy failed here (a per-user site-packages directory shadowed the environment: `GLIBC_2.27 not found`); with it, the imports and a smoke run without any conda activation worked.
- **Rmd 05 on emulated phASER tables** (exact gene tables derived from the truth, rendered after Rmd 01 and Rmd 02; both phased and unphased genotypes; every gate of `tests/synthetic/evaluate_stage3.R`): hap_strong 24 of 24 gene x sample cells significant; hap_moderate 16 of 24 (unphased ACAT 12); hap_lowdepth 15 of 18, all labelled "phASER only (no SNP tested unphased)" (ACAT 0); two_block 11 of 12 in both runs; null false positives 2 of 389 (phased) and 1 of 389 (unphased), unphased ACAT 0 of 402; phased runs: the direction of all 66 significant planted cells follows haplotype 1 (66 of 66), unphased runs: no direction on any row. Swapping the A and B labels of 230 of 470 rows left every p-value, padj and call unchanged (orientation invariance). The two_block gate was **restated**: the first version required all 12 cells, which is not a sound target at about 100 reads per cell (the missed cell had 67 against 34 of 101 reads, a realised fraction of 0.66 for a planted 0.8, and an exact test with the true dispersion misses it too). It now requires at least 90 percent of the calls of an oracle (an exact beta-binomial test with the true dispersion 0.005 on the same input counts): phased 11 against an oracle of 11 (bound 9.9), unphased 11 against 12 (bound 10.8, a margin of one cell). 21 tampered copies of the results each fail with the gate named for them.
- **phASER on aligned synthetic data** (scripts and Rmds cut from the skill text by an assembler script and submitted in chain order, both genotype modes, every job exit 0; the shipped `submit_chain.sh` was not used in this earlier run, only in the stub dry run below; it was used live in the acceptance run after it): 115 to 125 haplotype blocks per sample, 47 to 61 of them with several SNPs, and every multi-SNP block phased as in the truth (for example 52 of 52 and 57 of 57); 80 of 80 gene spans in the prep log; 77 to 80 genes with counts per sample; in the phased run `aCount` is haplotype 1 in all 12 to 13 planted genes of every sample, including genes whose haplotype 1 carries the reference allele and trans genes; the unphased run reports a single block per gene. One phASER array task took 16 to 21 s and about 100 MB on this tiny genome (single-threaded; this says nothing about real data). The phased two_block gate was restated as well: the second SNP cluster is only partly observable, because the first SNP of each two_block gene sits at transcript offset 2 and no read covers it after the WASP filter, so the maximum is `n_variants` 5,5 (not 6). The gate now requires `n_variants` to equal the number of the gene's heterozygous SNPs with at least one read in the sample's ASEReadCounter table (5,5 in all 6 samples), SNPs from both clusters and `gw_phased` 1; the unphased gate (at most 3 variants) is unchanged. 13 tampered copies each fail with the gate named for them.
- **Negative tests**: a contig-name mismatch (stale outputs are removed), a false phased-genotype claim, a contig name containing `_`, an install marker with only its first line, a sample name with `/` and `..`, `PHASED_GT` not 0 or 1 and a single-end BAM each stop `phaser_count.sh` with their message and leave no output (committed in `tests/synthetic/test_phaser_count_guards.sh`); Rmd 05 stops on a missing sample table, all-zero counts, a header-only table, empty count fields and a wrong column set (`tests/synthetic/test_rmd05_stops.R`), and attributes a warning of the dispersion fit to its sample.
- **Chain**: `submit_chain.sh` and `wait_chain.sh` were run against a stub scheduler in the dry run `tests/chain/dry_run_chain.sh` (41 checks ok, including the install-or-skip decisions, the F1 refusal and the recipes for re-submitting after a failure); this dry run is with a stub scheduler, not live; the acceptance run below then ran both scripts live.

Stage 3 synthetic acceptance run: DONE (2026-10-01)

Results of the live run (two outbred projects of the synthetic data, one with unphased and one with phased genotypes; the shipped `submit_chain.sh` and `wait_chain.sh` ran unmodified; all numbers are from synthetic data):

- **Projects and chain**: the unphased-genotype project used the installed phASER and skipped the installation (no installation job); the phased-genotype project used a fresh installation that ran live in 645 s. Every job exited 0. The dependencies were as designed: the phASER array after the alignment array (and, in the phased project, after the installation), and Rmd 05 after the phASER array and Rmd 02. Rmd 01, 02 and 05 were built without hand edits in both projects, and all 12 count tables per project were non-empty.
- **Evaluator**: 102 evaluator checks PASS and 0 FAIL across the six evaluator runs (phaser mode and Rmd mode, both projects; exit 0 in all).
- **phaser mode**: 47 to 61 multi-SNP blocks per sample, phase agreement 1.0000 with the truth, trans blocks present; two_block genes have 3 variants with unphased genotypes and 5 variants with phased genotypes (equal to the ASEReadCounter-covered SNPs); orientation was right in 12 to 13 of the 12 to 13 planted genes per sample, with genes of both phases.
- **Rmd 05 (both projects unless stated)**: input_complete 470 of 470; hap_strong 24 of 24; the restated two_block gate observed / oracle / bound 11 / 12 / at least 10.8 in both projects; null cells significant 1 of 388 (unphased) and 2 of 388 (phased) with a limit of 6; phased hap_lowdepth 15 of 18, all labelled "phASER only"; phased orientation and labels 66 of 66; unphased: no significant row carries a direction.
- **Reported, not gated**:
  - hap_moderate (24 gene x sample cells): phASER detects 16 with phased genotypes but only 13 with unphased genotypes, while the unphased ACAT test of Rmd 02 detects 14. With unphased genotypes phASER was therefore not better than ACAT on moderate genes (consistent with the single-block limit). An oracle with the true fraction detects 20 (phased) and 19 (unphased).
  - hap_lowdepth with unphased genotypes: 15 of 18 found only by phASER (no SNP tested unphased in those cells), unphased ACAT 0 of 18.
  - Null false positives: phASER 2 of 388 (phased) and 1 of 388 (unphased), unphased ACAT 0 of 402.
  - Aligned versus emulated gene counts (totalCount over the simulated fragments of both haplotypes, median per class; never above 1.000): phased aligned 0.932 (hap_lowdepth), 0.901 (hap_moderate), 0.890 (hap_strong), 0.933 (null), 0.894 (null_linked), 0.930 (two_block), sum ratio 0.908, against 1.000 in every emulated class; unphased aligned 0.913, 0.901, 0.890, null 0.765, null_linked 0.894, two_block 0.560 (sum 0.760), against emulated 1.000, 1.000, 1.000, null 0.861, null_linked 1.000, two_block 0.621 (sum 0.843); the unphased ratios below 1 for null and two_block come from the single-block rule. Aligned counts are about 0.89 to 0.93 of the simulated fragments.
  - Wall times: installation 645 s (measured once); each phASER array task 3 to 17 s. The synthetic genome is tiny (80 genes), so these timings and the resource request say nothing about real data.
- **Limits of what this run exercised**: `wait_chain.sh` infers the status of the prep, alignment, installation and phASER jobs from log text because `sacct` is unreachable, so a job killed without an error line would read "finished"; its failure, cancel and timeout branches were not exercised live (only in the stub dry run). The gates have thin margins (unphased two_block 11 against a bound of 10.8; phased hap_lowdepth 15 against a gate of 15). Everything in this run is synthetic data; the real-data smoke test follows.

Stage 3 real-data smoke test: DONE (2026-10-01)

Results of the real-data smoke test. Only the phASER part of the skill ran: the gene-span block of `prep_genotypes.sh` and `phaser_count.sh`, both cut from the skill text, on chr22 of two GEUVADIS RNA-seq samples (NA06986 and NA20808) with their phased 1000 Genomes phase 3 genotypes (SHAPEIT2) and the Ensembl release 75 GRCh37 GTF. Every heavy step ran as a SLURM job.

- **Inputs and adaptations** (each one limits what the test shows): the BAMs are the public GEUVADIS alignments, made with the GEM aligner on GRCh37, not STAR with WASP filtering as in the skill. Only chr22 and only two samples were used. The BAM contigs were renamed from `chr22` to `22` to match the VCF and the GTF. The BAMs had no duplicates marked, so the test marked them with `samtools markdup`. The only flag change: `--mapq 91` was chosen for GEM's MAPQ scale, while the skill's `--mapq 255` is for STAR. GEM gives its repeat alignments a MAPQ of at most 90 and its unique ones at least 98, and `--mapq 255` kept only 27% of GEM's unique reads. The skill's `--mapq 255` is right for STAR, which the skill uses, and wrong for BAMs from other aligners.
- **Gene spans** (whole GTF): 63,677 genes, 1,263 of them on chr22, and no gene left out for lying on more than one contig or strand. The GTF has 204 contigs whose names contain `_` (patch and haplotype contigs) that the BAMs do not have. The block then warned that phASER could not be used and printed an empty example list. This was a false alarm (phASER ran normally) and is now fixed: the span file keeps only genes on contigs of the FASTA, and the warning looks only at the FASTA's contig names, as the wizard's check in Step 7 does.
- **`phaser_count.sh`**: all four runs (two samples, phased and unphased genotypes) exited 0. There were 28,902 and 29,658 heterozygous sites, of which 2,360 and 2,406 were covered, and 284 and 331 multi-SNP haplotype blocks. 382 and 391 genes had counts; in the phased runs 377 and 385 of them were genome-wide phased.
- **Accuracy**: 1,099 of 1,140 directly connected SNP pairs (96.4%) of phASER's read-backed phase (unphased run) agree with the 1000 Genomes statistical phase. By read support, 96.5% / 97.1% / 98.9% agree for pairs linked by at least 2 / 5 / 10 reads. Whole multi-SNP blocks agree up to a flip in 94.6% to 95.8%. Agreement falls with distance: of the 13 pairs more than 10 kb apart, 9 agree. The cause of these discordant pairs is unexplained (it fits GEM's very long spliced gaps, but this was not proven).
- **Unphased genotypes on real data**: the single-best-block rule of `phaser_gene_ae` lost 76-77% of the covered SNPs of the covered genes (1,692 of 2,218 and 1,914 of 2,478), against 25 to 35 percent on the synthetic data. Gene totals were a median 0.74 of the phased totals. With these 2x75 reads only 2.4% to 2.7% of the heterozygous sites were phased with another site. The skill therefore recommends phASER mainly with phased genotypes (Step 7).
- **Time and memory, measured on chr22**: `phaser.py` took 2 to 3.6 s and about 117 MB per sample on one core, and `phaser_gene_ae` about 1.5 s. A separate run with 2.57 million heterozygous sites (the chr22 sites copied onto the other chromosomes, with the same reads) took 32 s and 1.59 GB, so memory grows mainly with the number of heterozygous sites. Without thread limits, numpy used 71 s of CPU time; with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` it used 6.6 s, with no loss of wall time. `phaser_count.sh` now sets both variables.
- **Whole genome, projected and not measured**: a whole-genome sample like these needs about 3 to 5 min and 2 to 3 GB. The projection uses an **assumed** genome-wide count of heterozygous sites (chr22's share of the genome length) and scales the read part linearly from chr22. It agrees with the upstream benchmark's 88 s for NA06986, but no whole-genome run was made. The large-genome request is now `-n 1 --mem=16G -t 2:00:00`, which leaves headroom over the projection.
- **Not run on real data**: STAR, WASP, ASEReadCounter, Rmd 01, 02 and 05, `submit_chain.sh` and `wait_chain.sh`, GRCh38, and a whole-genome sample.

What was **not exercised** in Stage 3: real data beyond the phASER smoke test above (STAR, WASP, ASEReadCounter and the Rmds have not run on real data), a whole-genome phASER run (the large-genome request `-n 1 --mem=16G -t 2:00:00` rests on chr22 measurements and a projection, so it is **not verified on real data** at whole-genome scale), single-end data, phased genotypes other than the 1000 Genomes statistical phase on chr22, contig names containing `_` on contigs of the FASTA (guarded only by a stop and a warning), population phasing (`phaser_pop`), and a Python 3.13 or older environment. Genes with alternative isoforms and overlapping genes were exercised only on chr22 of the real-data test.

What was **not exercised** in Stages 1 and 2: real biological data, designs with more than two conditions beyond the relabelled three-condition render test above, unbalanced real designs, multi-transcript genes (thinning uses the exon union), genome-scale runtime (projected from the unit-test benchmark only), the Mouse Genomes Project helper scripts end to end (`extract_mgp_parental_vcf.sh` and `concat_mgp_parental_vcf.sh`; the acceptance run used a user parental VCF), the standard Ensembl-style reference folder (the run used a custom FASTA and GTF), and SLURM mail notification options.

**Nothing on real biological data has been run for Stages 1 and 2; for Stage 3 only the phASER smoke test above (chr22, two samples) has.**

---

## Known limitations

- The overdispersion correction is a first-order approximation; on real data the gene-level test size is only approximate.
- **F1 gene test thinned to one SNP per `THIN_BP` window (approximate).** ASEReadCounter counts a fragment at every SNP it overlaps. Rmd 02's F1 gene likelihood-ratio test therefore uses each gene's SNPs thinned to one SNP per `THIN_BP` window in exon coordinates (deepest first), so no read pair is counted twice; the gene table reports `n_snps` and `n_snps_used`. Its rho (`rho_gene`) is the free-mean estimate on the thinned SNPs, with a fallback to the unthinned estimate when fewer than 5 genes keep 2 or more thinned SNPs (the fallback with the unthinned dispersion has not been verified when the true overdispersion is above RHO_MIN; because SNPs that share reads push that estimate low, the gene test may be anti-conservative in that case); the per-SNP tests keep the Stage 1 rho unchanged. In a fragment-level simulation the unthinned gene test had a size of 0.1450 at a nominal 0.05 (0.2084 in SNP-dense genes), the thinned test 0.0497. A correct gene test costs power in genes whose SNPs are close together: they keep one or two SNPs. On the Stage 1 synthetic acceptance F1 data (12 genes, SNPs 16-375 exonic bases apart) most genes keep a single SNP, and the thinned test finds 21 of 24 planted gene x sample tests (0 of 42 null false positives); the Stage 1 figure of 24 of 24 relied on counting shared read pairs several times and no longer applies. Other limits: the exon coordinates ignore alternative splicing (distant exons can share a fragment), and the kept SNP is chosen on the depth pooled over all samples of the project, so adding samples or conditions can change which SNP is kept in a dense gene. The outbred ACAT gene test is valid under dependence and is not thinned.
- One overdispersion (rho) per sample cannot capture gene-to-gene variation.
- Outbred rho is a trimmed fit on the central sites: with real overdispersion and many imbalanced sites it stays above the truth (lower power for moderate imbalance, no extra false positives), and with only a few dozen sites it sometimes cannot separate the imbalanced ones and stays near the naive value.
- The requirement of at least 5 usable genes (2 or more SNPs) to estimate rho at gene level is weak.
- Unphased outbred data cannot give gene-level direction.
- F1 mode requires strain A to carry the reference assembly's allele: sites where a substrain differs from the assembly are dropped (for example C57BL_6NJ against the B6J assembly GRCm39). A cross of two non-reference strains is not supported. With a parental VCF that has the two strains' genotype columns, only sites with strain A `0/0` and strain B `1/1` are used and the prep job stops on sites where strain A is `1/1` and strain B `0/0` (unless the user confirms a substrain and sets `A_ALT_SITES=drop`); a sites-only parental VCF cannot be checked and is taken on trust.
- Contig names of the parental or genotype VCFs must match the FASTA (`1` vs `chr1`); the prep job stops before any alignment if they do not, and the skill does not rename contigs itself.
- Paths containing spaces or commas are not supported; the skill stops and asks for another path.
- Indels and multi-allelic sites are excluded.
- RNA-derived genotypes (for example the nf-core/rnavar VCF) are circular and biased toward balance; external genotypes are preferred.
- **Stage 3 (phASER) and its limits.**
  - **Third-party tool.** phASER is maintained as a "Fast Beta" revival (PEJ Lab, repository secastel/phaser) with no releases and no tags, so the skill pins the commit `aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301` and records the package versions. Upstream has open issues that can matter (for example `phaser_gene_ae` index errors or a mis-parsed column on some inputs, #78 and #85, and zero counts when contig names differ between the BED and the BAM, #65 and #86). Bumping the pin needs the procedure of Step 18 and a new test run.
  - **Single-threaded.** With the pinned Python 3.14 environment `phaser.py` with `--threads` above 1 fails with `NameError: name 'args' is not defined` (Python 3.14's default process start method no longer passes phASER's settings to the worker processes), so every sample runs single-threaded (`--threads 1`, `-n 1`, samples side by side as array tasks); a Python 3.13 or older pin was never tested, so whether it would allow more threads is unknown.
  - **Unphased genotypes** (what the prep jobs write and what an RNA-derived VCF gives): `phaser_gene_ae` keeps only the most-covered haplotype block of each gene and drops the other blocks, which costs coverage (on the synthetic test, 31 of the 50 genes with at least 2 covered SNPs lost SNPs in one sample, and the gene step used 113 of its 153 heterozygous SNPs; phased genotypes gave about 25 to 35 percent more coverage per gene). On real data the loss is much larger: on chr22 of two GEUVADIS samples 76-77% of the covered SNPs of the covered genes were lost and the gene totals were a median 0.74 of the phased totals, so phASER is recommended mainly with phased genotypes. The haplotype-versus-ACAT power gain measured by the fragment-level simulation is optimistic there, because the simulation puts all SNPs of a gene into one block. There is no direction.
  - **Direction and comparisons.** A direction is given only for genome-wide phased rows (phased genotypes, `gw_phased` TRUE); otherwise the haplotype labels are arbitrary. Haplotype A of one individual is unrelated to haplotype A of another, so across individuals only the size of the imbalance can be compared, and there is no differential ASE on phASER counts (use Rmd 04). The synthetic data do not exercise the rule that genes of a phased run that are not genome-wide phased get no direction.
  - **Contig names containing `_` are not supported** on contigs of the FASTA: `phaser_gene_ae` splits variant IDs on `_`. The prep job warns, the wizard does not offer phASER for such a reference, and `phaser_count.sh` stops as a last guard; no site is filtered silently. Contigs with `_` that exist only in the GTF (the patch contigs of an Ensembl GTF) are harmless: their genes are left out of the gene-span file.
  - **Gene spans and dispersion.** The gene spans include introns (Rmd 02 uses exonic SNPs only) and overlapping genes share reads. Rmd 05's cohort dispersion is the median of the samples' own dispersions, so a sample's p-values depend on the other samples of the run (adding or removing a sample can change them; do not split the samples into several runs), and the dispersion stays too high when many genes are imbalanced or when a sample has few genes (power numbers above).
  - **Resources measured on chr22, projected for a whole genome.** The phASER request is `-n 1 --mem=4G -t 0:30:00` for small genomes and `-n 1 --mem=16G -t 2:00:00` otherwise. On real data it was measured on chr22 of two samples only (2 to 3.6 s and about 117 MB for `phaser.py`; 32 s and 1.59 GB with 2.57 million heterozygous sites). The whole-genome figures (about 3 to 5 min and 2 to 3 GB per sample) are projected from an assumed count of heterozygous sites and were not verified on real data by a whole-genome run; the request keeps about 5 times the projected memory and about 24 times the projected time. The Rmd 05 runtime and memory per million haplotype reads are projected from one measurement.
  - **Installation.** It needs network access on a compute node (git clone and conda; 10 to 20 minutes) and uses a cluster-specific conda module (`module add miniconda3/v4` and `/home/software/conda/miniconda3/bin/condainit`); on another cluster this line has to be adapted.
  - **Large outputs.** Excel holds at most 1,048,575 data rows per sheet: beyond that the Rmd 05 xlsx sheets keep only the significant rows and the full tables are in the `.tsv.gz` files. The `SNP` and `Gene` sheets of Rmd 02 and Rmd 04 have the same limit (their checkpoints hold the full tables).
  - **Real data only for the phASER part, on chr22**: the real-data smoke test (see the validation section) ran the gene-span block and `phaser_count.sh` on chr22 of two GEUVADIS samples aligned with GEM on GRCh37; STAR, WASP, ASEReadCounter, Rmd 05 and the chain have run on synthetic data only.
- **Stage 2 tests and their limits.**
  - Paired outbred differential test: power is limited with few individuals and depends on the fraction of SNPs that change. With 4 individuals and 10 percent of SNPs changed the power is 0.8653 / 0.9174 (phase-heterogeneous / consistent) against an oracle of 0.9022 / 0.9424; with 30 percent changed (phase-heterogeneous) it is 0.6698 against an oracle of 0.9163, because the pair dispersion is then inflated (0.0454 for a true 0.02). This is conservative (no false positives added). The pair dispersion includes truly changed SNPs, and a SNP is tested only with at least 2 informative individuals (tested fraction 0.4244 / 0.6961 / 0.8434 / 0.9646 for 2 / 3 / 4 / 6 individuals at 60 percent heterozygosity); gene calls have no direction.
  - The F1 GLM tests are conservative: null sizes are about 0.036 at a true dispersion of 0.02 (nominal 0.05), because of the maximum rule and the `RHO_MIN` floor.
  - One dispersion per analysis plus each gene's own estimate when larger (with at least 2 residual df); with few replicates the per-gene estimate is noisy.
  - Thinning uses exon-union coordinates, so alternative splicing can bring two distant exons into one fragment; `THIN_BP` should match the longest fragment. The kept SNP is chosen on the depth pooled over all samples of the project, so adding samples or conditions can change which SNP is kept in a dense gene.
  - Rmd 01 and Rmd 02 do not set X aside in F1 mode: male F1 animals carry one X, so X genes can look imbalanced and X sites can trip the reference-bias flag; read X rows with care. Rmd 03 and 04 exclude X/Y/MT.
  - X, Y and MT are not analysed in Rmd 03 and 04. Alt, random and unplaced contigs are treated as autosomes.
  - A reciprocal design needs both directions with replicates; condition-by-direction interactions are not modelled, and b1 is one effect common to all conditions.
- The unit-test job needs at least 8 CPUs.
- `sacct` is unavailable on the test cluster, so `wait_chain.sh` infers that an Rmd job succeeded from its result files and takes failures from the logs; where `sacct` works, the same script has not been compared with it.
- The mouse helper was only partly exercised (see the validation section).
