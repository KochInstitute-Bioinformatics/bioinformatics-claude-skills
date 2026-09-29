# `/ase-pipeline` — Allele-Specific Expression Skill (Stage 1)

An interactive Claude Code skill that walks you through a complete allele-specific expression (ASE) analysis of bulk RNA-seq on an HPC cluster running SLURM and Singularity: reference and genotype preparation, allele-aware alignment (STAR), allele counting (GATK ASEReadCounter), and downstream statistics as two R Markdown reports plus a standalone summary page. Two designs share one wizard: **F1 crosses** (for example B6 x AJ mice) and **outbred / human** samples with per-individual genotypes.

Stage 1 covers **per-sample allelic imbalance only**. Reciprocal F1 analysis, differential ASE between conditions and phASER phasing are later stages.

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
| Cached biocontainers | In `$NXF_SINGULARITY_CACHEDIR` (or `~/.singularity/cache`): STAR 2.7.10b, GATK 4.4.0.0, bcftools 1.20, samtools 1.21, Picard 3.1.1 (`depot.galaxyproject.org-singularity-*.img`); a missing image is downloaded by the generated helper script, never in the foreground on the login node |
| `bulkrnaseq` image | `/net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif`, the only place R is available; needs `aod`, `lme4`, `openxlsx`, `tidyverse`, `GenomicRanges`, `rtracklayer` (checked by a tiny `sbatch` job in Step 2) |
| Genotypes | F1: a parental-difference VCF, or the built-in Mouse Genomes Project helper (needs internet on a compute node). Outbred: one VCF per individual |
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
| 7 | Analysis menu (only the per-sample analysis is implemented in Stage 1) |
| 8 | Constants: `MIN_DEPTH`, `FDR_SIG`, `ABS_DEV_SIG`, `BIAS_TOL`, `RHO_MIN` |
| 9 | Resources for the array and prep jobs |
| 10 | Reference and genotype preparation scripts (masked reference or per-individual het VCFs, STAR index, optional mouse helper) |
| 11 | Per-sample array scripts: STAR, duplicate marking, ASEReadCounter, empty-table guard |
| 12 | Rmd 01: import, site filters, reference-bias QC, WASP removal report (outbred), checkpoint |
| 13 | Rmd 02: per-SNP and per-gene tests, figures, xlsx, `summary_numbers.tsv` |
| 14 | Statistics functions (shipped code, unit-tested) |
| 15 | Submission order with dependencies, the standalone summary page, notes for the assistant |

---

## Output files

Everything is written under `{CWD}/results/{YYMMDD}_{WD_NAME}/` (raw FASTQ, BAM and VCF inputs are never modified):

| Output | Description |
|--------|-------------|
| `scripts/` | Generated `sbatch` scripts: prep (`prep_f1_reference.sh` or `prep_genotypes.sh`), optional mouse helper, per-sample array (`align_count_f1.sh` or `align_wasp_count.sh`), `run_01_import_qc.sh`, `run_02_imbalance.sh` |
| `bam/` | Coordinate-sorted BAMs with read groups, duplicates marked (not removed), plus indexes |
| `ase_counts/` | ASEReadCounter table per sample (and `wasp_stats.tsv` per sample in outbred mode) |
| `ase_checkpoint.rds`, `ase_imbalance_checkpoint.rds` | Checkpoints of Rmd 01 and Rmd 02 |
| `{date}_{project}_01_import_qc.html`, `..._02_imbalance.html` | The two stage reports |
| `{date}_{project}_ASE_imbalance.xlsx` | Sheets `SNP`, `Gene`, `Summary` |
| `..._ASE_sites_vs_depth.pdf`, `..._ASE_genes.pdf` | Figures |
| `summary_numbers.tsv` | One row per sample with the headline numbers (written by Rmd 02) |
| `{WD_NAME}_summary_report.html` | Standalone page (inline CSS, no R, relative links to the reports and the xlsx, reference-bias flags shown prominently) |

The jobs run in one dependency chain: prep -> per-sample array (`afterok`) -> Rmd 01 (`afterok`) -> Rmd 02 (`afterok`), preceded by the mouse helper chain when used.

---

## Key design points

- **F1 mapping bias: third-allele masked reference.** At every parental SNP the genome base is replaced with a third allele (neither strain's), so neither strain is favoured in mapping; STAR runs on the masked genome and ASEReadCounter counts against the het-sites VCF.
- **Outbred mapping bias: WASP.** STAR's WASP re-mapping filter tags every alignment (`vW`); alignments with `vW` 2-7 are removed before counting, and Rmd 01 reports how many were removed and from which allele. WASP does not promise less bias, and none is claimed.
- **Per-sample beta-binomial tests.** Each sample gets its own overdispersion. F1 uses a free-mean-per-gene estimate with a bias correction (a plain H0-based estimate is inflated by real imbalance), floored at `RHO_MIN`, for both SNP-level exact tests and a gene-level likelihood-ratio test on the SNP-level counts (counts are not summed). Outbred uses the H0-based estimate.
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
- **Unit-tested statistics**: null size of the F1 gene-level test 0.047-0.056 pooled (worst seed 0.064) and power 0.99, plus rho recovery, edge cases and ACAT (`tests/r/test_ase_stats.R`; needs a job with at least 8 CPUs).
- **Rmd rendering tests** on perturbed copies of one synthetic sample (F1 and outbred, Rmd 01 and Rmd 02), including the `summary_numbers.tsv` output.

What is **PENDING**: the full end-to-end synthetic acceptance run (all samples, all scripts chained with the dependencies of Step 15, recovery of the planted genes with the correct direction) has not been completed yet; its result will be recorded here when it is completed.

What was **never run**: the Step 15 summary HTML page (it has never been generated; only its input `summary_numbers.tsv` was checked against the Summary sheet of the xlsx); the run scripts `run_01_import_qc.sh` and `run_02_imbalance.sh` (never submitted; the Rmd renders above used direct `rmarkdown::render` calls); the full sbatch dependency chain of Step 15 (never submitted end to end); and the href-existence check of the summary page.

**Nothing on real biological data has been run.**

---

## Known limitations

- The overdispersion correction is a first-order approximation; on real data the gene-level test size is only approximate.
- One overdispersion (rho) per sample cannot capture gene-to-gene variation.
- The requirement of at least 5 usable genes (2 or more SNPs) to estimate rho at gene level is weak.
- Unphased outbred data cannot give gene-level direction.
- F1 mode requires strain A to carry the reference assembly's allele: sites where a substrain differs from the assembly are dropped (for example C57BL_6NJ against the B6J assembly GRCm39).
- Indels and multi-allelic sites are excluded.
- RNA-derived genotypes (for example the nf-core/rnavar VCF) are circular and biased toward balance; external genotypes are preferred.
- Stage 1 covers only per-sample imbalance; reciprocal F1, differential ASE and phASER are later stages.
- The unit-test job needs at least 8 CPUs.
