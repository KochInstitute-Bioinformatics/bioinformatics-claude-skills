# `/ase-pipeline` — Allele-Specific Expression Skill (Stages 1 and 2)

An interactive Claude Code skill that walks you through a complete allele-specific expression (ASE) analysis of bulk RNA-seq on an HPC cluster running SLURM and Singularity: reference and genotype preparation, allele-aware alignment (STAR), allele counting (GATK ASEReadCounter), and downstream statistics as R Markdown reports plus a standalone summary page. Two designs share one wizard: **F1 crosses** (for example B6 x AJ mice) and **outbred / human** samples with per-individual genotypes.

Stages 1 and 2 cover per-sample allelic imbalance, the reciprocal F1 analysis (strain versus parent-of-origin effects) and differential ASE between conditions. phASER phasing (Stage 3) is not included.

**What Stage 2 adds**

- **Reciprocal F1 (Rmd 03)**: strain and parent-of-origin effects from crosses in both directions.
- **Differential ASE (Rmd 04)**: changes in allelic imbalance between conditions, for F1 crosses (GLM per gene and per SNP) and for outbred individuals (paired, direction-free test).
- **Thinned Rmd 02 F1 gene test**: the gene-level likelihood-ratio test of Stage 1 now uses each gene's SNPs thinned to one SNP per `THIN_BP` window, so that no read pair is counted twice.

**New requirements of Stage 2**

- `samples.csv` column `cross_direction` must be exactly `AxB` or `BxA` for the reciprocal analysis (A = the strain carrying the reference allele, the mother in `AxB`); each direction needs at least 2 samples.
- Replicates mean at least 2 samples in a condition (F1) or at least 2 individuals sampled in both conditions of a contrast (outbred); a contrast without them is not tested and is listed as such.
- Differential ASE asks for a reference condition (the baseline that every other condition is compared with). In F1 mode the condition must not be confounded with the cross direction (Rmd 04 stops when it is).

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
| Singularity >= 3.10 | Loaded via `module add singularity/3.10.4`; the only module the skill ever loads |
| Cached biocontainers | In `$NXF_SINGULARITY_CACHEDIR` (or `~/.singularity/cache`): STAR 2.7.10b, GATK 4.4.0.0, bcftools 1.20, samtools 1.21, Picard 3.1.1 (`depot.galaxyproject.org-singularity-*.img`); a missing image is downloaded by the prep job (it fetches all five, Picard included, so the per-sample array job only checks for them), never in the foreground on the login node |
| `bulkrnaseq` image | `/net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif`, the only place R is available; needs `openxlsx`, `tidyverse`, `GenomicRanges`, `rtracklayer` (checked by a tiny `sbatch` job in Step 2) |
| Genotypes | F1: a parental-difference VCF (plain `.vcf` or bgzipped; sites-only is fine, strain names default to `C57BL_6NJ` / `A_J` when the header has none), or the built-in Mouse Genomes Project helper (needs internet on a compute node). Outbred: one VCF per individual |
| Internet access | Only for downloading missing containers, references and the Mouse Genomes Project data |

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
| 7 | Analysis menu: per-sample (always), reciprocal F1, differential ASE (with reference condition); phASER later |
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
| `{WD_NAME}_summary_report.html` | Standalone page (inline CSS, no R, relative links to the reports and the xlsx, reference-bias flags shown prominently) |

The jobs run in one dependency chain: prep -> per-sample array (`afterok`) -> Rmd 01 (`afterok`) -> Rmd 02 (`afterok`), and, when selected, Rmd 03 / Rmd 04 (`afterok` on Rmd 01, beside Rmd 02), preceded by the mouse helper chain when used. A wait script follows every submitted job and stops when only jobs that can never start remain.

---

## Key design points

- **F1 mapping bias: third-allele masked reference.** At every parental SNP the genome base is replaced with a third allele (neither strain's), so neither strain is favoured in mapping; STAR runs on the masked genome and ASEReadCounter counts against the het-sites VCF.
- **Outbred mapping bias: WASP.** STAR's WASP re-mapping filter tags every alignment (`vW`); alignments with `vW` 2-7 are removed before counting, and Rmd 01 reports how many were removed and from which allele. WASP does not promise less bias, and none is claimed.
- **Per-sample beta-binomial tests.** Each sample gets its own overdispersion. F1 uses a free-mean-per-gene estimate with a bias correction (a plain H0-based estimate is inflated by real imbalance), floored at `RHO_MIN`, for both SNP-level exact tests and a gene-level likelihood-ratio test on the SNP-level counts (counts are not summed); the gene test and its own rho estimate use each gene's SNPs thinned to one SNP per `THIN_BP` window (exon coordinates), so no read pair is counted twice. Outbred (unphased, so no per-gene mean can be fitted) uses a trimmed estimate: the H0 beta-binomial fitted on the central sites only (each site's central 90 percent region), with a truncated likelihood that corrects for the trimming, floored at `RHO_MIN`. The naive all-sites H0 estimate is shown only for comparison (`rho_h0_naive`): real imbalance inflates it, and on the synthetic acceptance data it detected 0 of 16 planted outbred genes.
- **Reciprocal F1 (Rmd 03).** Per gene, `logit(p) = b0 + b1 * d`, where p is the strain-A fraction and d = +1 for `AxB`, -1 for `BxA`. b0 is the strain effect and b1 the parent-of-origin effect. Rows are the samples, with the counts summed over the gene's thinned SNPs. The dispersion is estimated between replicate animals, and each gene uses the larger of the pooled value, its own estimate (with at least 2 residual df) and `RHO_MIN`. Both effects get likelihood-ratio tests. X, Y and MT are excluded. b1 is not affected by a constant mapping bias toward one strain (it cancels between the two directions); b0 is, so the reference-bias flags are repeated next to the strain results.
- **Differential ASE (Rmd 04).** F1: a GLM per gene and per SNP with the condition term and, when both directions are present, the cross-direction term; the Rmd stops when condition and direction are confounded. Outbred: a paired, direction-free per-SNP test summed over the individuals sampled in both conditions (the allele carrying a regulatory variant differs between individuals, so no direction is reported), genes combined with ACAT. No mixed model (lme4) is used: a random individual effect is not estimable with 2-6 individuals, and a shared condition slope would cancel phase-dependent changes.
- **Base-R models.** The GLM, the pair dispersion and the tests are base-R functions shipped in Step 14 and extracted by the unit tests, so the tested code is the shipped code; Rmd 03 and 04 need no package beyond those of Rmd 01 and 02.
- **Unphased outbred data.** Without phasing, SNPs of a gene cannot be pooled: per-SNP tests are combined per gene with ACAT and labelled "unphased, no direction".
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

Stage 2 synthetic acceptance run: PENDING

What was **not exercised**: real biological data, phASER (Stage 3), designs with more than two conditions beyond the relabelled three-condition render test above, unbalanced real designs, multi-transcript genes (thinning uses the exon union), genome-scale runtime (projected from the unit-test benchmark only), the Mouse Genomes Project helper scripts end to end (`extract_mgp_parental_vcf.sh` and `concat_mgp_parental_vcf.sh`; the acceptance run used a user parental VCF), the standard Ensembl-style reference folder (the run used a custom FASTA and GTF), and SLURM mail notification options.

**Nothing on real biological data has been run.**

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
- phASER (Stage 3) is not implemented; the menu lists it as available in a later stage.
- **Stage 2 tests and their limits.**
  - Paired outbred differential test: power is limited with few individuals and depends on the fraction of SNPs that change. With 4 individuals and 10 percent of SNPs changed the power is 0.8653 / 0.9174 (phase-heterogeneous / consistent) against an oracle of 0.9022 / 0.9424; with 30 percent changed (phase-heterogeneous) it is 0.6698 against an oracle of 0.9163, because the pair dispersion is then inflated (0.0454 for a true 0.02). This is conservative (no false positives added). The pair dispersion includes truly changed SNPs, and a SNP is tested only with at least 2 informative individuals (tested fraction 0.4244 / 0.6961 / 0.8434 / 0.9646 for 2 / 3 / 4 / 6 individuals at 60 percent heterozygosity); gene calls have no direction.
  - The F1 GLM tests are conservative: null sizes are about 0.036 at a true dispersion of 0.02 (nominal 0.05), because of the maximum rule and the `RHO_MIN` floor.
  - One dispersion per analysis plus each gene's own estimate when larger (with at least 2 residual df); with few replicates the per-gene estimate is noisy.
  - Thinning uses exon-union coordinates, so alternative splicing can bring two distant exons into one fragment; `THIN_BP` should match the longest fragment. The kept SNP is chosen on the depth pooled over all samples of the project, so adding samples or conditions can change which SNP is kept in a dense gene.
  - X, Y and MT are not analysed in Rmd 03 and 04. Alt, random and unplaced contigs are treated as autosomes.
  - A reciprocal design needs both directions with replicates; condition-by-direction interactions are not modelled, and b1 is one effect common to all conditions.
- The unit-test job needs at least 8 CPUs.
- `sacct` is unavailable on the test cluster, so `wait_chain.sh` infers that an Rmd job succeeded from its result files and takes failures from the logs; where `sacct` works, the same script has not been compared with it.
- The mouse helper was only partly exercised (see the validation section).
