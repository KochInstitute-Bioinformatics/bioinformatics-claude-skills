# ase-pipeline Stage 3 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add Stage 3 to the shipped `ase-pipeline` skill: optional phASER read-backed phasing in outbred mode. It has four parts: a pinned, one-time installation job; a per-sample phASER array job with gene-level haplotype counts; a new Rmd 05 that tests the haplotype counts per sample and gene and reports the result next to the unphased ACAT result of Rmd 02 (never instead of it); and the wizard, chain and summary integration. As in Stages 1 and 2, the statistics come first and are shown correct by simulation, and a live synthetic acceptance run comes last.

**Architecture:** phASER runs outside the Singularity containers, from a conda environment that a compute job installs under a user-chosen directory `{PHASER_HOME}`. The repository is pinned to commit `aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301`, and the package versions are pinned to the resolution the feasibility gate tested. The run-time scripts use the environment's `bin` directory directly, with `PYTHONNOUSERSITE=1`. The per-sample array (`phaser_count.sh`) reads the same WASP-filtered, duplicate-marked BAM and heterozygous VCF that ASEReadCounter uses, and writes `phaser/{sample}.gene_ae.txt`. The one new statistics function, `hap_gene_test`, lives in the Step 14 block and is unit-tested by extraction. Its test is symmetric in the two haplotypes, because without phased genotypes the labels A and B are arbitrary. Rmd 05 reads the Rmd 01 and Rmd 02 checkpoints plus the gene tables, and runs `afterok` on the phASER array and Rmd 02. Rmds 01-04 do not change.

**Tech Stack:** Markdown skill file; bash/awk; R in `/net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif` (base R, `parallel`, `Biostrings` for the synthetic generator, `openxlsx`, `tidyverse`); phASER at commit aa1f8ec (Python 3.14.7 conda environment: numpy 2.5.3, scipy 1.18.1, pandas 3.0.6, pysam 0.24.1, intervaltree 3.2.1, samtools/bcftools/htslib 1.24, bedtools 2.31.1); conda 4.14 from `module add miniconda3/v4`, used by the installation job only; SLURM (`sbatch -p bcc`); the Stage 1 upstream (STAR 2.7.10b WASP, GATK 4.4.0.0, bcftools 1.20, samtools 1.21, Picard 3.1.1 containers).

**Spec:** `docs/superpowers/specs/2026-09-29-ase-pipeline-design.md` (binding; sections "Rmd 02" (outbred with phASER), "Stage 3 — phASER", "Staging", "Testing and acceptance"). Task 1 Step 9 amends the spec where the planner decisions below depart from it. Read it together with the feasibility gate report `.superpowers/sdd/stage3-phaser-gate-report.md` (the source of every proven phASER fact below), the Stage 2 plan `docs/superpowers/plans/2026-09-30-ase-pipeline-stage2.md`, and the ledgers `.superpowers/sdd/2026-09-29-ase-pipeline-stage1/progress.md` and `.superpowers/sdd/2026-09-30-ase-pipeline-stage2/progress.md`.

## Planner decisions (need user review)

Each decision can be reversed independently before execution.

1. **Separate Rmd 05 (new Step 19), not an extension of Rmd 02.** phASER results are reported next to the unphased per-SNP/ACAT results of Rmd 02, never instead of them. Reason: the Stage 1/2 outputs stay byte-for-byte the same, and a phASER failure cannot break Rmd 02.
2. **Test:** per sample and gene, an exact two-sided beta-binomial test of `aCount` out of `totalCount` against 0.5, symmetric in A and B; BH within each sample; `sig` = padj < `FDR_SIG` and major-haplotype fraction − 0.5 ≥ `ABS_DEV_SIG`. Reason: without phased genotypes the labels are arbitrary, so only the folded size can be tested; per sample matches Rmd 02.
3. **Dispersion:** `rho_own` = the trimmed central-region estimator `bb_estimate_rho_trim` over the sample's genes; `rho_cohort` = the median of `rho_own` over the samples; `rho_used = max(rho_cohort, rho_own, RHO_MIN)`. Reason: an H0 fit is inflated by real imbalance (Stage 1 D1); a free mean per gene is impossible because the orientation is arbitrary; the maximum rule (as in Stage 2) protects noisy samples and samples with few genes. Cost: conservative for clean samples in a noisy cohort. A pooled trimmed fit over all sample-gene rows was rejected: it expands every count of every row and would need tens of GB on a real cohort.
4. **Installation:** one compute job clones the repository and checks out commit `aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301`. It creates the environment with `conda env create -p {PHASER_HOME}/env` from a skill-owned environment file that pins the exact versions resolved and tested in the gate (python 3.14.7 ...). It uses neither the unpinned upstream `environment.yml` nor a `python=3.12` pin; the 3.12 pin was never tested. The job also writes a bump procedure. Reason: the installation is reproducible and uses only tested versions.
5. **The `module` exception:** `module add miniconda3/v4` (plus `source /home/software/conda/miniconda3/bin/condainit`) appears only in `setup_phaser_env.sh`; it is the single exception to the "only `singularity/3.10.4`" rule. The run-time scripts load no module and do not activate conda. Instead they put `{PHASER_HOME}/env/bin` first on `PATH` and export `PYTHONNOUSERSITE=1`. Task 3 verifies this without activation; if it fails, `conda activate` (the gate-proven way) is used and the ledger records the change.
6. **phASER input and flags:** `bam/{sample}.bam` (WASP-filtered, duplicates marked), `genotypes/{individual}.het.vcf.gz` (the heterozygous VCF that STAR and ASEReadCounter use), `--paired_end 1 --mapq 255 --baseq {MIN_BASEQ} --pass_only 0 --unique_ids 1 --gw_phase_vcf {PHASED_GT}`. phASER is offered only for paired-end data. Reason: these flags were proven in the gate; `--paired_end 0` was never tested.
7. **Phased genotypes:** a Step 7 question sets `{PHASED_GT}` (1 or 0). `phaser_count.sh` stops if `{PHASED_GT}` = 1 and fewer than 90 % of the het genotypes are phased. For unphased genotypes the single-best-block limitation is told to the user. A direction ("haplotype A higher / B higher", never REF/ALT) is reported only for genome-wide phased rows with phased genotypes. Haplotype A is the left GT allele's haplotype (read from `phaser_gene_ae.py` at aa1f8ec, and checked in Tasks 5 and 8). Across individuals, only the magnitudes can be compared.
8. **Gene-span BED:** `prep_genotypes.sh` always builds it from the GTF **exons**: the span of each `gene_id`, 0-based, one row per gene, genes on more than one contig or strand left out. It is written in every outbred project. Reason: this avoids a conditional template edit (Stage 2 lesson), works for GTFs without `gene` lines, and gives the same `gene_id`s as Rmd 02.
9. **Contigs with `_` stop `phaser_count.sh`.** `phaser_gene_ae` parses variant IDs with `int(id.split("_")[1])`. Changing `--id_separator` was never tested. Ensembl naming, the skill's default, has no such contig.
10. **New synthetic generator** `simulate_ase_stage3.R` (new file, seed 20261001), instead of reusing the Stage 2 outbred data. Reason: the Stage 2 data hold only `cis` configurations, and no low-depth or two-block genes, so they cannot show that phasing or the pooling of reads across SNPs is used. The existing generators are not touched.
11. **Rmd 05 runs `afterok` on both the phASER array and Rmd 02,** because it needs Rmd 02's gene table for the side-by-side comparison. Cost: if Rmd 02 fails, Rmd 05 does not run.
12. **Out of scope:** differential ASE on phASER counts (the orientation is arbitrary, and even with phased genotypes haplotype A of one individual is unrelated to A of another), cohort-level phASER tests, phaser_pop and population phasing.
13. **The real-data smoke test (Task 9) is optional; its dataset is a user decision.** Every other task is independent of it. Until it runs, the phASER resource request (`-n 4 --mem=32G -t 2:00:00`) is labelled "not verified on real data".

## Proven facts (do not re-derive)

**From the feasibility gate** (`.superpowers/sdd/stage3-phaser-gate-report.md`, 2026-10-01; spike workspace `/net/bmc-lab3/data/bcc/ase_scratch/phaser/`):

- **Repository:**
  - https://github.com/secastel/phaser, default branch `master`, GPL-3.0, maintained by PEJ Lab as a "Fast Beta".
  - There are no releases and no tags.
  - The pinned commit is `aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301` (2026-03-21, "Update gene file links to Google Drive"). `phaser.py` prints "phASER v1.2.0".
  - The code is Python 3 only. No Cython build is needed: the `read_variant_map.so` check is commented out, and the README "Setup" step is stale.
- **Bioconda `phaser`** is version 0.1.1 for Python 2.7. Never use it.
- **Installation on a compute node** (job 11377680):
  - Compute nodes have internet: GitHub returns 200 and conda-forge 206.
  - `module add miniconda3/v4; source /home/software/conda/miniconda3/bin/condainit` gives conda 4.14.0 with the classic solver (no mamba).
  - `conda env create -p <prefix> -f environment.yml` succeeded in 1252 s (about 21 min).
  - Versions resolved: python 3.14.7, numpy 2.5.3, scipy 1.18.1, pandas 3.0.6, pysam 0.24.1, intervaltree 3.2.1, samtools / bcftools / htslib 1.24, bedtools 2.31.1. They work with glibc 2.17.
  - Sizes: the environment is 632 MB; the package cache was 1.9 GB and was deleted afterwards.
- **`PYTHONNOUSERSITE=1` is mandatory.** Without it, `~/.local/lib/python3.14` (from `pip --user`) shadows the environment's numpy and the import fails with `GLIBC_2.27 not found`.
- **`/usr/bin/time` does not exist** on the compute nodes.
- **Default flags fail on the skill's VCFs:**
  - `--pass_only 1` keeps 0 sites, because the FILTER column is `.` ("FATAL ERROR: No heterozygous sites that passed all filters"). Fix: `--pass_only 0`.
  - The ID column is `.`, so `--unique_ids 1` is needed; `phaser_gene_ae` parses IDs as `chr_pos_ref_alt`.
- **Working command** (job 11377683; all 8 samples rc 0):

  ```
  phaser.py --vcf <het.vcf.gz> --bam <bam> --sample <vcf sample> --paired_end 1 --mapq 255 --baseq 10 --pass_only 0 --unique_ids 1 --threads 1 --temp_dir <dir> --o <prefix>
  phaser_gene_ae.py --haplotypic_counts <prefix>.haplotypic_counts.txt --features <bed> --o <out>
  ```
- **Outputs** (each prefixed with `<prefix>`):
  - `.haplotypic_counts.txt`, with columns `contig start stop variants variantCount variantsBlacklisted variantCountBlacklisted haplotypeA haplotypeB aCount bCount totalCount blockGWPhase gwStat max_haplo_maf bam aReads bReads`;
  - `.haplotypes.txt`, `.allelic_counts.txt`, `.variant_connections.txt`, `.allele_config.txt`, `.vcf.gz` and `.tbi`.
  - The gene_ae columns are `contig start stop name aCount bCount totalCount log2_aFC n_variants variants gw_phased bam`.
- **phASER on the synthetic outbred data** (Stage 2 `obd`):
  - It phased 59-90 variant connections per sample, giving 27-40 multi-SNP blocks.
  - 100 % of the phased pairs matched the truth. Every configuration was `cis`, because that generator puts every ALT allele on one haplotype.
  - `--mapq 255` and `--mapq 0` gave the same result (12590 and 12592 reads).
  - The WASP-filtered BAM is the right input: phASER does no WASP of its own.
- **Unphased genotypes:** every block has `blockGWPhase` `0/1`, so `phaser_gene_ae` keeps only the most-covered block of each gene. On ob_s1, 31 of 50 genes lost SNPs and 113 of 153 SNPs were used.
- **Phased genotypes plus `--gw_phase_vcf 1`** (job 11377687):
  - Every covered gene is `gw_phased` = 1 and every SNP is used.
  - Gene coverage was 25-35 % higher than with unphased genotypes, with the same accuracy against the truth (mean |fraction − truth| 0.025-0.031).
- **Genes without any covered het SNP** are written with `aCount` = `bCount` = 0, `log2_aFC` `inf` and `gw_phased` 1. Filter on `totalCount > 0`.
- **Contig-name mismatches** between the BAM, the VCF and the BED give silent zero counts (upstream issues #65 / #86).
- **Runtime:**
  - On the tiny synthetic data: 1.4-1.9 s and about 100 MB per sample.
  - Upstream benchmark: 88-466 s per real sample on one thread.
  - Suggested request: `-n 4 --mem=32G -t 2:00:00`. It is not verified on real data.

**From reading `phaser_gene_ae/phaser_gene_ae.py` at aa1f8ec** (lines 95-218; reading only, no run):

- **Orientation:**
  - In a block with `blockGWPhase` `0|1` (and `gwStat` ≥ `--gw_cutoff` 0.9), A is "GW haplotype 0", that is, the haplotype of the **left** allele of the phased GT.
  - For `1|0` blocks the A and B counts are swapped, so A is still the left allele's haplotype.
- **How blocks are summed:**
  - Genome-wide phased blocks overlapping a gene are summed.
  - Of the unphased blocks, only the one with the most reads is kept.
  - Whichever total is larger is written: `gw_phased` 1 (phased sum) or 0 (single block).
- **Read deduplication:** within one block, reads are deduplicated as a set of read IDs over the block's variants, so a read pair that covers two SNPs of a block counts once.
- **Variant IDs:** they are split on `--id_separator` and read with `int(fields[1])`. A contig name containing `_` therefore crashes the gene step.
- **Features file:** BED columns 1-4 (chr, 0-based start, stop, name); every row is one feature.

**Carried over from Stages 1 and 2:**

- **Dispersion estimation:**
  - A rho fitted under H0 on all sites is inflated by real imbalance (0.08-0.10 against a binomial truth; it detected 0 of 16 planted outbred genes).
  - `bb_estimate_rho_trim` fixes this for unphased data. It is symmetric in the two alleles, because it folds the counts around n/2.
  - `RHO_MIN` = 0.01 caps the power.
- **ACAT** caps p at 0.99.
- **Simulation gates:** 10 seeds, with a pooled gate (≤ 0.065) and a worst-seed gate (≤ 0.09). The unit-test job needs 8 CPUs.
- **R** runs only in the `bulkrnaseq` image, only through `sbatch -p bcc`.
- **`sacct` is unreachable;** job states come from `scontrol show job`.
- **Upstream facts:**
  - STAR WASP needs a het-only, single-sample VCF.
  - ASEReadCounter needs a read group plus a bgzipped, indexed VCF with genotypes.
  - Rmd 01 reads tables by sample name, never by glob.
- **"Wizard-by-script"** means scripts and Rmds generated from the skill's code blocks by extraction and placeholder substitution; the interactive dialogue is not driven.
- **Lesson:** in Stages 1 and 2 the acceptance run found defects that the unit tests missed (Stage 1 D1 and I1-I6; Stage 2 T1-T3), and every whole-branch review found real bugs. So the Stage 3 synthetic data exercise every new branch, and the acceptance runs both genotype modes live.

## Global Constraints

- **Skill file:** `ase-pipeline/ase-pipeline.md` stays self-contained. Every finite-choice question is a numbered list. Open questions are allowed only for the email address and for paths (`{PHASER_HOME}` is a path).
- **Paths and data safety:**
  - Results go to `results/{TODAY}_{WD_NAME}`, with `{TODAY}` = `YYYY-MM-DD`.
  - Raw FASTQ, BAM and VCF inputs are read-only.
  - Paths with spaces or commas are refused (Step 0 rule, applied to `{PHASER_HOME}` too).
- **Cluster use:**
  - Never run R, Singularity, conda, python or any heavy work on the login node. Everything goes through `sbatch -p bcc`, and debugging uses `salloc`/`srun`.
  - Requests are at most 64 G and 4 h by default.
  - Resource requests must be honest.
- **Tools:**
  - Upstream tools come from Singularity containers; their only module is `module add singularity/3.10.4`, checked with `|| exit 1`.
  - The single exception is `module add miniconda3/v4` in `setup_phaser_env.sh` (planner decision 5).
  - phASER runs only from `{PHASER_HOME}/env`, with `PYTHONNOUSERSITE=1` exported in every script that runs its Python.
  - Never use bioconda `phaser`.
  - `/tmp` is node-local, so every output and temporary file goes under `{RESULTS_DIR}` or `{PHASER_HOME}`.
- **Pinned commit:** `aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301`, set in exactly three script templates: `setup_phaser_env.sh`, `phaser_count.sh` and `submit_chain.sh`.
- **Generated scripts** start with `set -uo pipefail` and check every exit code themselves (no blind `set -e`). Chained jobs use `--parsable` and `--dependency=afterok`.
- **Flags:**
  - Every `--flag` in the skill text appears in a recorded fixture (`star_help.txt`, `gatk_ASEReadCounter_help.txt`, `phaser_help.txt`, `phaser_gene_ae_help.txt`) or in the checker allowlist.
  - Every flag on a `phaser.py` command must be in `phaser_help.txt`, and every flag on a `phaser_gene_ae.py` command in `phaser_gene_ae_help.txt`.
  - Never add a STAR, GATK or phASER flag to the allowlist.
- **Rmd conventions:**
  - Rmds are self-contained (no `source()`) and use `cache = FALSE` and `options(scipen = 9)`.
  - Bioconductor packages load before `tidyverse`; verbs are `dplyr::`-prefixed.
  - Each Rmd checks its constants against the Rmd 01 checkpoint.
  - Figure counts must equal the Summary counts.
  - Count files are read by sample name, never by glob.
- **Statistics code:**
  - It lives only between `# --- ase-stats-begin` and `# --- ase-stats-end`, and is pasted verbatim into Rmd 02-05.
  - The unit tests extract it from the skill text.
- **Simulation and acceptance gates are binding:**
  - The thresholds are fixed in this plan.
  - If a gate fails, the implementer reports BLOCKED with the measured numbers. They never weaken a gate, change a threshold or swap a method.
  - A gate may be restated only by a controller ruling recorded in the ledger, based on evidence (for example an oracle comparison, as in the Stage 2 Task 1 ruling).
- **Checker:**
  - Every new `need`/`forbid` line must be shown to FAIL on the pre-task skill (`git show HEAD:ase-pipeline/ase-pipeline.md > $SCR/pre.md`, then run the checker on it).
  - Every new guard or rule in a test tool must be shown to FAIL on an injected defect.
- **Scratch and data locations:**
  - Everything that is not committed goes under `SCR=/net/bmc-lab3/data/bcc/yannvrb/ase_stage3_scratch`.
  - The synthetic data go to `/net/bmc-lab3/data/bcc/yannvrb/ase_synthetic_stage3` (seed 20261001).
  - The phASER installation of the tests lives in `$SCR/phaser_home`.
  - Never write generated data into the repository.
- **Stage 1/2 outputs unchanged:**
  - Rmds 01-04 and the count tables are byte-for-byte the same.
  - The only upstream change is the extra `reference/genes_span.bed` from `prep_genotypes.sh`.
  - `simulate_ase_data.R` and `simulate_ase_stage2.R` are not modified.
- **Git:**
  - Branch `feat/ase-pipeline-stage3`.
  - Show `git status` and `git diff --stat` before each commit.
  - Never `git push`, and never merge without the user's approval.
  - Commit messages end with a blank line and `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>` (or the trailer the harness prescribes for the executing model).

## Review Focus

Each line names an input or condition, the behaviour a user would expect, and the test that pins it.

1. **Contig names differ between the gene BED, the het VCF and the BAM** (`chr1` against `1`). Expected: the job stops with a message that names both naming styles, and no empty gene table reaches Rmd 05. Test: Task 5 Step 7(a) (a renamed BED makes `phaser_count.sh` exit 1 with the guard message and leave no output). Also Task 5 Step 2, which adds the BED-versus-`.fai` guard in `prep_genotypes.sh`.
2. **The user says the genotypes are phased, but the VCF is not (or the reverse).** Expected: a stop when phased is claimed but the VCF is not phased, and a warning in the reverse case, so that no direction is ever reported from arbitrary labels. Tests: Task 5 Step 7(b) (`PHASED_GT=1` with an unphased VCF stops); Task 4 Step 6 (with unphased genotypes, no significant row carries a direction).
3. **Haplotype labels flip between genes and samples** (unphased genotypes). Expected: every p-value, dispersion and call is unchanged. Tests: Task 1 section 1 (a full swap and a random half flipped give identical p and dispersions); Task 4 Step 5 (Rmd 05 rendered on flipped and on unflipped gene tables gives identical results).
4. **Genes written with zero counts (`inf`, `gw_phased` 1) and genes below `MIN_DEPTH`.** Expected: neither kind is tested nor counted as genome-wide phased, and neither enters the dispersion fit. Tests: Task 1 section 1 (zero, NA and Inf rows); Task 4 Step 5 (the evaluator checks that no zero-count row is in the gene table and that every row below `MIN_DEPTH` has p `NA`).
5. **A missing or stale per-sample phASER output** (a failed task, or a re-run into the same results directory). Expected: `phaser_count.sh` deletes the sample's earlier outputs before running, and Rmd 05 stops with a message that names the missing samples, never reading another sample's file. Tests: Task 5 Step 7(d) (a stale file is gone after a failed run); Task 4 Step 6 (a deleted `gene_ae.txt` stops Rmd 05 with the sample named).

## File Structure

| File | Responsibility | Task |
|---|---|---|
| `ase-pipeline/ase-pipeline.md` | Step 14: `hap_gene_test` and its prose. Step 10: gene-span block in `prep_genotypes.sh`. New Step 18 (installation, per-sample phASER) and Step 19 (Rmd 05). Steps 2, 7, 9, 15 and Notes. | 1, 3, 4, 5, 6 |
| `ase-pipeline/tests/r/test_ase_stats_stage3.R` | Simulation unit tests of `hap_gene_test` (extracted from the skill). | 1 |
| `ase-pipeline/tests/r/run_stats_tests.sh` | Adds the `stage3` target. | 1 |
| `ase-pipeline/tests/synthetic/simulate_ase_stage3.R` | Stage 3 ground truth: 6 outbred individuals with cis/trans phase, low-depth and two-block genes, phased and unphased VCFs, direct count tables and emulated phASER gene tables. | 2 |
| `ase-pipeline/tests/synthetic/check_simulation_stage3.sh` | Exact-quota and consistency checks of the Stage 3 dataset. | 2 |
| `ase-pipeline/tests/synthetic/prove_checks_stage3.sh` | Injected-defect proofs of that check. | 2 |
| `ase-pipeline/tests/synthetic/check_phase_reads_stage3.R` | Read-level check that read pairs carry the planted cis/trans phase. | 2 |
| `ase-pipeline/tests/synthetic/cut_block.sh` | Prints the first fenced code block after a heading or anchor line of the skill. | 3 |
| `ase-pipeline/tests/fixtures/phaser_help.txt`, `phaser_gene_ae_help.txt`, `phaser_env.txt` | Recorded `--help` output and environment of the pinned installation. | 3 |
| `ase-pipeline/tests/fixtures/verification.md` | A Stage 3 section. | 3 |
| `ase-pipeline/tests/synthetic/evaluate_stage3.R` | Evaluates Rmd 05 (`rmd`) and the raw phASER outputs (`phaser`) against the truth. | 4, 5 |
| `ase-pipeline/tests/synthetic/assemble_outbred_from_skill.sh` | Generates every outbred script and Rmd of a project from the skill's code blocks, as the wizard would. | 5 |
| `ase-pipeline/tests/chain/dry_run_chain.sh` | Runs `submit_chain.sh` and `wait_chain.sh` against stub `sbatch`/`squeue`/`scancel`. | 6 |
| `ase-pipeline/tests/check_skill.sh` | Static checker; one block per task. | 1, 3-7 |
| `ase-pipeline/README.md`, root `README.md` | Stage 3 documentation and validation status. | 7, 8 |
| `docs/superpowers/specs/2026-09-29-ase-pipeline-design.md` | Spec amendments. | 1 |

**Step numbering** is fixed now. Steps 0-17 keep their numbers. phASER becomes `## Step 18 — phASER: installation and per-sample haplotype counts (outbred, optional)` and Rmd 05 becomes `## Step 19 — Rmd 05: phASER gene-level haplotype imbalance (outbred, optional)`. Both go after Step 17 and before `## Notes for the assistant`, Step 18 first. Step 15 refers forward to them.

All repo paths are relative to `/net/bmc-lab3/data/bcc/yannvrb/Github_repos/bioinformatics-claude-skills`.

**Task order and dependencies:**

- Task 1 (statistics) comes before any Rmd text.
- Task 2 (synthetic data) does not depend on the other tasks.
- Task 3 (installation) creates `$SCR/phaser_home`, which Tasks 5 and 8 use.
- Task 4 (Rmd 05) needs Tasks 1 and 2. It renders on the emulated gene tables, so it does not need phASER.
- Task 5 (`phaser_count.sh` on aligned data) needs Tasks 2, 3 and 4.
- Task 6 (wizard and chain) needs Tasks 3-5.
- Task 7 (README) needs Tasks 1-6.
- Task 8 (acceptance) is run by the controller.
- Task 9 is optional.

---

### Task 1: `hap_gene_test` with simulation unit tests (plus spec amendments)

**Files:**
- Modify `ase-pipeline/ase-pipeline.md`:
  - the Step 14 stats block: insert the function directly above the line `# --- ase-stats-end`;
  - the Step 14 intro sentence ("These base-R functions are pasted verbatim into Rmd 02, Rmd 03 and Rmd 04 ...");
  - a new Step 14 prose bullet after the last bullet ("**SNPs sharing read pairs.** ...").
- Create: `ase-pipeline/tests/r/test_ase_stats_stage3.R`
- Modify: `ase-pipeline/tests/r/run_stats_tests.sh`, `ase-pipeline/tests/check_skill.sh` (new block before the final `[ $fail -eq 0 ]` line), `docs/superpowers/specs/2026-09-29-ase-pipeline-design.md`

**Interfaces:**
- Consumes the existing block functions `bb_estimate_rho_trim(x, n, keep = 0.9, max_iter = 100)`, which returns `c(rho, central_frac)` (NA below 20 sites), and `bb_pvalue(x, n, rho, p0 = 0.5)`. Neither is changed.
- Produces (Tasks 4, 5 and 8 use these exact names):
  - `hap_gene_test(a, b, sample, rho_min, min_genes = 20)`. `a`, `b` and `sample` are vectors with one entry per sample-gene row. It returns `list(p, rho_used, samples)`:
    - `p`: numeric, one per row; `NA` for unusable rows (a or b missing, infinite or negative, or a + b = 0) and for every row when no sample has `min_genes` usable rows.
    - `rho_used`: per row, the dispersion of that row's sample.
    - `samples`: a data.frame with columns exactly `sample, n_genes, rho_own, central_frac, rho_cohort, rho_used`, one row per distinct sample (sorted).
  - `sbatch -p bcc -o <log> ase-pipeline/tests/r/run_stats_tests.sh [stage1|stage2|stage3|all]`.

**Gates** (fixed here; `gate()` means pooled ≤ 0.065 over 10 seeds and worst seed ≤ 0.09, the Stage 1/2 rule):

| Section | Measure | Gate |
|---|---|---|
| 1 | symmetry: full swap and random half flipped | p and dispersions equal within 1e-12 |
| 1 | binomial data, 10 seeds × 6 samples | `rho_used` = `RHO_MIN` in every sample |
| 2 | null size, cohort 6 × 2000 genes, rho ∈ {0.005, 0.02, 0.05} × imbalanced fraction ∈ {0, 0.2, 0.4} | `gate()` each |
| 2 | power, rho 0.02, 20 % imbalanced at 0.75/0.25 | ≥ 0.7 and ≥ 0.9 × oracle power (true rho) |
| 2 | few genes: 6 samples × 20 / 30 / 60 genes | `gate()` each |
| 2 | one sample with 15 genes (cohort rho used) | `gate()` on that sample's null genes |
| 2 | heterogeneous cohort (true rho 0.08, 0.05, 0.02, 0.02, 0.01, 0.01) | `gate()` overall and for the 0.08 sample |
| 3 | fragment level, null genes, haplotype test, lambda 300 and 40 | `gate()` each |
| 3 | fragment level, genes with ≥ 4 SNPs, h 0.75: power of the haplotype test against ACAT of per-SNP tests (not tested = not detected) | haplotype ≥ ACAT, both lambdas |
| 3 | lambda 40, 8-SNP genes: fraction tested by the haplotype test | ≥ 0.9 (the ACAT fraction is reported) |
| 4 | runtime: one sample of 20,000 genes, depth median 300 (max 50,000), one core | ≤ 120 s and ≤ 4096 MB max memory |

- [ ] **Step 1: Write the failing test script**

Create `ase-pipeline/tests/r/test_ase_stats_stage3.R` with exactly this content:

```r
# Stage 3 statistics tests (phASER gene-level haplotype test); hap_gene_test is extracted from the shipped skill text.
# Usage: Rscript test_ase_stats_stage3.R <ase-pipeline.md>     (needs >= 8 CPUs: parallel::mclapply, mc.cores = 8)
args <- commandArgs(trailingOnly = TRUE); skill <- args[1]
txt <- readLines(skill)
b <- grep("^# --- ase-stats-begin", txt); e <- grep("^# --- ase-stats-end", txt)
stopifnot(length(b) == 1, length(e) == 1, e > b)
eval(parse(text = txt[(b + 1):(e - 1)]))
fail <- 0
ok <- function(cond, msg) { if (!isTRUE(cond)) { cat("FAIL:", msg, "\n"); fail <<- 1 } else cat("ok  ", msg, "\n") }
if (!exists("hap_gene_test")) { cat("FAIL: missing from the stats block: hap_gene_test\n"); quit(status = 1) }
const <- function(name) {   # a Step 8 default, read from the skill text (never hard-coded here)
  l <- grep(paste0("^\\| `", name, "` \\|"), txt, value = TRUE); stopifnot(length(l) == 1)
  as.numeric(trimws(strsplit(l, "\\|")[[1]][3]))
}
RHO_MIN <- const("RHO_MIN"); MIN_DEPTH <- const("MIN_DEPTH")
ok(RHO_MIN == 0.01 && MIN_DEPTH == 10, sprintf("Step 8 defaults read from the skill: RHO_MIN %g, MIN_DEPTH %g", RHO_MIN, MIN_DEPTH))
MC <- 8
seeds_rbind <- function(seeds, f) {
  r <- parallel::mclapply(seeds, f, mc.cores = MC)
  stopifnot(!any(vapply(r, function(x) inherits(x, "try-error"), logical(1))))
  do.call(rbind, r)
}
gate <- function(m, col, label, pooled = 0.065, worst = 0.09)
  ok(mean(m[, col]) <= pooled && max(m[, col]) <= worst,
     sprintf("%s: pooled %.4f (<= %.3f), worst seed %.3f (<= %.2f)", label, mean(m[, col]), pooled, max(m[, col]), worst))
report <- function(m, label) cat(sprintf("     %s: %s\n", label, paste(sprintf("%s %.4f", colnames(m), colMeans(m, na.rm = TRUE)), collapse = ", ")))
depth_ln <- function(m, med = 200, sdlog = 1) round(exp(rnorm(m, log(med), sdlog))) + 5
draw <- function(n, h, rho) { pr <- if (rho > 0) rbeta(length(n), h * (1 - rho) / rho, (1 - h) * (1 - rho) / rho) else h; rbinom(length(n), n, pr) }

# ---- 1. interface, symmetry (haplotype labels are arbitrary), unusable rows, few genes, the RHO_MIN floor
set.seed(1); n <- depth_ln(400); a <- draw(n, 0.5, 0.02); smp <- rep(c("s1", "s2"), each = 200)
r <- hap_gene_test(a, n - a, smp, RHO_MIN)
ok(identical(names(r), c("p", "rho_used", "samples")) &&
   identical(names(r$samples), c("sample", "n_genes", "rho_own", "central_frac", "rho_cohort", "rho_used")),
   "hap_gene_test returns p, rho_used and samples (sample, n_genes, rho_own, central_frac, rho_cohort, rho_used)")
ok(length(r$p) == 400 && all(is.finite(r$p)) && all(r$p >= 0 & r$p <= 1), "one p per row, every p in [0, 1]")
r_sw <- hap_gene_test(n - a, a, smp, RHO_MIN)
ok(isTRUE(all.equal(r$p, r_sw$p, tolerance = 1e-12)) && isTRUE(all.equal(r$samples, r_sw$samples, tolerance = 1e-12)),
   "swapping the haplotype labels (a <-> b) in every row leaves p and the dispersions unchanged")
fl <- runif(400) < 0.5
r_fl <- hap_gene_test(ifelse(fl, n - a, a), ifelse(fl, a, n - a), smp, RHO_MIN)
ok(isTRUE(all.equal(r$p, r_fl$p, tolerance = 1e-12)) && isTRUE(all.equal(r$samples, r_fl$samples, tolerance = 1e-12)),
   "flipping the labels of a random half of the rows leaves p and the dispersions unchanged")
ae <- c(a, 0, 5, NA, 7, Inf, 30); be <- c(n - a, 0, 0, 3, NA, 2, 0); se <- c(smp, rep("s1", 6))
re <- hap_gene_test(ae, be, se, RHO_MIN)
ok(all(is.na(re$p[400 + c(1, 3, 4, 5)])), "zero-coverage, NA and Inf rows: p NA")
ok(is.finite(re$p[402]) && re$p[406] < 1e-3, "monoallelic rows are tested: 5 of 5 gives a finite p, 30 of 30 gives p < 1e-3")
keep <- c(1:400, 402, 406)
ok(isTRUE(all.equal(re$p[keep], hap_gene_test(ae[keep], be[keep], se[keep], RHO_MIN)$p, tolerance = 1e-12)),
   "NA, Inf and zero rows do not enter the dispersion fit (same p as without them)")
ok(identical(as.integer(re$samples$n_genes), c(202L, 200L)), "n_genes counts the usable rows per sample")
r19 <- hap_gene_test(a[1:19], (n - a)[1:19], rep("s1", 19), RHO_MIN)
ok(all(is.na(r19$p)) && is.na(r19$samples$rho_cohort), "no sample with 20 usable genes: nothing tested (p NA, rho_cohort NA)")
ix <- c(1:15, 201:400)
mix <- hap_gene_test(a[ix], (n - a)[ix], rep(c("s1", "s2"), c(15, 200)), RHO_MIN)
ok(is.na(mix$samples$rho_own[1]) && all(is.finite(mix$p[1:15])) &&
   isTRUE(all.equal(mix$samples$rho_used[1], max(mix$samples$rho_cohort[1], RHO_MIN))),
   "a sample with fewer than 20 genes is tested with max(rho_cohort, RHO_MIN)")
floor_run <- function(s) { set.seed(100 + s); nn <- depth_ln(12000); aa <- draw(nn, 0.5, 0)
  hap_gene_test(aa, nn - aa, rep(paste0("s", 1:6), each = 2000), RHO_MIN)$samples$rho_used }
fl0 <- unlist(parallel::mclapply(1:10, floor_run, mc.cores = MC))
ok(length(fl0) == 60 && all(fl0 == RHO_MIN), sprintf("binomial data: rho_used sits at the RHO_MIN floor in all 60 samples (max %.4f)", max(fl0)))

# ---- 2. cohorts: null size with real imbalance present, power against an oracle, few genes, heterogeneous dispersion
cohort_run <- function(seed, rho, frac_imb, S = 6, G = 2000, h_imb = 0.75, rho_s = NULL, G_s = NULL) {
  set.seed(seed)
  Gs <- if (is.null(G_s)) rep(G, S) else G_s; rs <- if (is.null(rho_s)) rep(rho, S) else rho_s; S <- length(Gs)
  d <- do.call(rbind, lapply(seq_len(S), function(s) {
    nn <- depth_ln(Gs[s]); imb <- runif(Gs[s]) < frac_imb
    h <- ifelse(imb, ifelse(runif(Gs[s]) < 0.5, h_imb, 1 - h_imb), 0.5)
    data.frame(sample = paste0("s", s), a = draw(nn, h, rs[s]), n = nn, imb = imb, rho_true = rs[s], stringsAsFactors = FALSE)
  }))
  r <- hap_gene_test(d$a, d$n - d$a, d$sample, RHO_MIN)
  po <- vapply(seq_len(nrow(d)), function(i) bb_pvalue(d$a[i], d$n[i], d$rho_true[i]), numeric(1))   # oracle: the true dispersion
  sz <- function(p, sel) mean(p[sel] < 0.05, na.rm = TRUE)
  c(size = sz(r$p, !d$imb), size_s1 = sz(r$p, !d$imb & d$sample == "s1"),
    power = if (frac_imb > 0) sz(r$p, d$imb) else NA_real_, oracle_power = if (frac_imb > 0) sz(po, d$imb) else NA_real_,
    oracle_size = sz(po, !d$imb), rho_cohort = r$samples$rho_cohort[1], rho_used_mean = mean(r$samples$rho_used))
}
k <- 0
for (rho in c(0.005, 0.02, 0.05)) for (fi in c(0, 0.2, 0.4)) {
  k <- k + 1
  m <- seeds_rbind(1:10, function(s) cohort_run(1000 * k + s, rho, fi))
  lab <- sprintf("cohort 6 x 2000 genes, true rho %.3f, %d%% imbalanced", rho, round(100 * fi))
  report(m, lab); gate(m, "size", paste("null size,", lab))
  if (rho == 0.02 && fi == 0.2)
    ok(mean(m[, "power"]) >= 0.7 && mean(m[, "power"]) >= 0.9 * mean(m[, "oracle_power"]),
       sprintf("power, h 0.75, rho 0.02, 20%% imbalanced: %.4f (>= 0.7 and >= 0.9 x oracle %.4f)", mean(m[, "power"]), mean(m[, "oracle_power"])))
}
avg_runs <- function(seed0, reps, ...) { x <- t(vapply(seq_len(reps), function(j) cohort_run(seed0 + j, ...), numeric(7))); colMeans(x, na.rm = TRUE) }
for (G in c(20, 30, 60)) {
  m <- seeds_rbind(1:10, function(s) avg_runs(20000 + 1000 * G + 100 * s, 30, rho = 0.02, frac_imb = 0.2, G = G))
  report(m, sprintf("few genes: 30 cohorts of 6 samples x %d genes, rho 0.02, 20%% imbalanced", G))
  gate(m, "size", sprintf("null size with %d genes per sample", G))
}
m <- seeds_rbind(1:10, function(s) avg_runs(300000 + 1000 * s, 100, rho = 0.02, frac_imb = 0.2, G_s = c(15, rep(300, 5))))
report(m, "one sample with 15 genes in a cohort of 6 (tested with max(rho_cohort, RHO_MIN))")
gate(m, "size_s1", "null size of the 15-gene sample")
m <- seeds_rbind(1:10, function(s) cohort_run(40000 + s, 0.02, 0.2, rho_s = c(0.08, 0.05, 0.02, 0.02, 0.01, 0.01)))
report(m, "heterogeneous cohort (true rho 0.08, 0.05, 0.02, 0.02, 0.01, 0.01), 20% imbalanced")
gate(m, "size", "heterogeneous cohort, null size over all samples"); gate(m, "size_s1", "heterogeneous cohort, null size of the rho 0.08 sample")

# ---- 3. fragment level: phASER-style haplotype counts (each read pair once) against the unphased per-SNP path + ACAT (Rmd 02)
frag_gene <- function(pos, h, nf, span, rl = 100) {   # one gene and sample: fragments, haplotype-1 indicator, SNPs each fragment covers
  L <- round(runif(nf, 150, 350)); st <- floor(runif(nf) * (span - L + 1)) + 1; A <- runif(nf) < h
  cv <- (outer(st, pos, "<=") & outer(st + rl - 1, pos, ">=")) | (outer(st + L - rl, pos, "<=") & outer(st + L - 1, pos, ">="))
  anyc <- rowSums(cv) > 0
  list(snp_a = colSums(cv & A), snp_n = colSums(cv), hap_a = sum(anyc & A), hap_n = sum(anyc))
}
frag_run <- function(seed, lambda, S = 6, G = 1500, rho_bio = 0.02, frac_imb = 0.3, span = 1500) {
  set.seed(seed)
  k <- sample(c(1, 2, 4, 8), G, TRUE); imb <- runif(G) < frac_imb
  pos <- lapply(k, function(kk) sort(sample(400:1000, kk)))             # the SNPs of a gene lie within 600 transcript bases
  out <- lapply(seq_len(S), function(s) {
    h <- ifelse(imb, ifelse(runif(G) < 0.5, 0.75, 0.25), 0.5)
    hs <- rbeta(G, h * (1 - rho_bio) / rho_bio, (1 - h) * (1 - rho_bio) / rho_bio)
    lapply(seq_len(G), function(g) frag_gene(pos[[g]], hs[g], rpois(1, lambda), span))
  })
  acat_p <- matrix(NA_real_, G, S); hap_a <- matrix(NA_real_, G, S); hap_n <- matrix(NA_real_, G, S)
  for (s in seq_len(S)) {   # unphased path of Rmd 01/02: sites with depth >= MIN_DEPTH, trimmed rho floored at RHO_MIN, bb_pvalue, ACAT
    sa <- unlist(lapply(out[[s]], `[[`, "snp_a")); sn <- unlist(lapply(out[[s]], `[[`, "snp_n")); sg <- rep(seq_len(G), k)
    use <- sn >= MIN_DEPTH
    rt <- bb_estimate_rho_trim(sa[use], sn[use])[["rho"]]; rh <- if (is.na(rt)) NA_real_ else max(rt, RHO_MIN)
    ps <- rep(NA_real_, length(sa))
    if (!is.na(rh)) ps[use] <- vapply(which(use), function(i) bb_pvalue(sa[i], sn[i], rh), numeric(1))
    acat_p[, s] <- vapply(seq_len(G), function(g) acat(ps[sg == g]), numeric(1))
    hap_a[, s] <- vapply(out[[s]], `[[`, numeric(1), "hap_a"); hap_n[, s] <- vapply(out[[s]], `[[`, numeric(1), "hap_n")
  }
  okn <- hap_n >= MIN_DEPTH
  r <- hap_gene_test(ifelse(okn, hap_a, NA), ifelse(okn, hap_n - hap_a, NA), rep(paste0("s", seq_len(S)), each = G), RHO_MIN)
  hp <- matrix(r$p, G, S); I <- matrix(imb, G, S); big <- matrix(k >= 4, G, S); K8 <- matrix(k == 8, G, S)
  det <- function(p, sel) mean(!is.na(p[sel]) & p[sel] < 0.05)          # not tested counts as not detected
  c(hap_size = mean(hp[!I] < 0.05, na.rm = TRUE), acat_size = mean(acat_p[!I] < 0.05, na.rm = TRUE),
    hap_power_k4 = det(hp, I & big), acat_power_k4 = det(acat_p, I & big),
    hap_tested_k8 = mean(!is.na(hp[K8])), acat_tested_k8 = mean(!is.na(acat_p[K8])))
}
for (lam in c(300, 40)) {
  m <- seeds_rbind(1:10, function(s) frag_run(50000 + 100 * lam + s, lam)); report(m, sprintf("fragment level, %d fragments per gene", lam))
  gate(m, "hap_size", sprintf("fragment level, lambda %d, haplotype test null size", lam))
  ok(mean(m[, "hap_power_k4"]) >= mean(m[, "acat_power_k4"]),
     sprintf("fragment level, lambda %d, genes with >= 4 SNPs: haplotype-test power %.4f >= unphased ACAT power %.4f", lam,
             mean(m[, "hap_power_k4"]), mean(m[, "acat_power_k4"])))
  if (lam == 40) ok(mean(m[, "hap_tested_k8"]) >= 0.9,
     sprintf("lambda 40, 8-SNP genes: tested by the haplotype test %.4f (>= 0.9); by the unphased path %.4f (reported)",
             mean(m[, "hap_tested_k8"]), mean(m[, "acat_tested_k8"])))
}

# ---- 4. runtime and memory: one sample of genome scale on one core (Rmd 05 runs the samples one after the other)
set.seed(900); nn <- pmin(round(exp(rnorm(20000, log(300), 1.2))) + 5, 50000); aa <- draw(nn, 0.5, 0.02)
invisible(gc(reset = TRUE))
tt <- system.time(rr <- hap_gene_test(aa, nn - aa, rep("s1", 20000), RHO_MIN))[["elapsed"]]
g <- gc(); mx <- sum(g[, ncol(g)])
ok(tt <= 120 && mx <= 4096, sprintf("runtime: one sample, 20,000 genes (depth median 300, max 50,000): %.0f s (<= 120), max memory %.0f MB (<= 4096); projected 100 samples %.2f h on one core",
                                    tt, mx, tt * 100 / 3600))
quit(status = fail)
```

- [ ] **Step 2: Add the `stage3` target to the runner**

In `ase-pipeline/tests/r/run_stats_tests.sh`, change the usage comment to `[stage1|stage2|stage3|all]`, and insert this before the `echo "EXIT $rc"` line:

```bash
if [ "$WHAT" = stage3 ] || [ "$WHAT" = all ]; then
  singularity exec --bind /net/bmc-lab3 "$SIF" Rscript ase-pipeline/tests/r/test_ase_stats_stage3.R ase-pipeline/ase-pipeline.md || rc=1
fi
```

- [ ] **Step 3: Run the Stage 3 tests to verify they fail**

Run: `mkdir -p $SCR/logs && sbatch -p bcc -o $SCR/logs/t1_red_%j.out ase-pipeline/tests/r/run_stats_tests.sh stage3`
Expected: the log shows `FAIL: missing from the stats block: hap_gene_test` and `EXIT 1`.

- [ ] **Step 4: Add `hap_gene_test` to the stats block**

In Step 14, insert this code directly above the line `# --- ase-stats-end`:

```r
hap_gene_test <- function(a, b, sample, rho_min, min_genes = 20) {
  # phASER gene-level haplotype counts, one row per sample and gene (a = aCount, b = bCount of phaser_gene_ae). Symmetric in a
  # and b: without phased genotypes the labels A and B are arbitrary per gene and sample, so only |a / (a + b) - 0.5| carries
  # information, and every quantity below is unchanged when a and b are swapped in any row.
  # Dispersion: rho_own = bb_estimate_rho_trim over the sample's usable genes (folded central-region H0 fit, truncation-
  # corrected, so genes with real imbalance do not inflate it; NA below min_genes genes); rho_cohort = the median of rho_own
  # over the samples; per sample rho_used = max(rho_cohort, rho_own, rho_min). Test per row: bb_pvalue against 0.5 (exact,
  # two-sided). Rows with a or b missing, infinite or negative, or a + b = 0, are not used and get p = NA, as does every row
  # when no sample has min_genes usable rows.
  n <- a + b
  use <- is.finite(a) & is.finite(b) & a >= 0 & b >= 0 & n > 0
  sample <- as.character(sample); smp <- sort(unique(sample))
  own <- vapply(smp, function(s) {
    i <- use & sample == s
    if (sum(i) >= min_genes) bb_estimate_rho_trim(a[i], n[i]) else c(rho = NA_real_, central_frac = NA_real_)
  }, c(rho = 0, central_frac = 0))
  coh <- if (all(is.na(own["rho", ]))) NA_real_ else median(own["rho", ], na.rm = TRUE)
  ru <- if (is.na(coh)) rep(NA_real_, length(smp)) else pmax(coh, ifelse(is.na(own["rho", ]), 0, own["rho", ]), rho_min)
  per <- data.frame(sample = smp, n_genes = unname(vapply(smp, function(s) sum(use & sample == s), integer(1))),
                    rho_own = unname(own["rho", ]), central_frac = unname(own["central_frac", ]), rho_cohort = coh,
                    rho_used = unname(ru), stringsAsFactors = FALSE)
  rr <- per$rho_used[match(sample, per$sample)]
  p <- rep(NA_real_, length(n))
  for (j in which(use & !is.na(rr))) p[j] <- bb_pvalue(a[j], n[j], rr[j])
  list(p = p, rho_used = rr, samples = per)
}
```

Change the Step 14 intro sentence "These base-R functions are pasted verbatim into Rmd 02, Rmd 03 and Rmd 04 (`{TODAY}_{WD_NAME}_02_imbalance.Rmd`, `..._03_reciprocal.Rmd`, `..._04_differential.Rmd`)" to "These base-R functions are pasted verbatim into Rmd 02, Rmd 03, Rmd 04 and Rmd 05 (`{TODAY}_{WD_NAME}_02_imbalance.Rmd`, `..._03_reciprocal.Rmd`, `..._04_differential.Rmd`, `..._05_phaser.Rmd`)".

- [ ] **Step 5: Run all unit tests and handle the gates**

Run: `sbatch -p bcc -o $SCR/logs/t1_green_%j.out ase-pipeline/tests/r/run_stats_tests.sh all`

Expected:
- The Stage 1 and Stage 2 files print only `ok` lines; nothing else in the block changed.
- The Stage 3 file prints only `ok` and `report` lines.
- The log ends with `EXIT 0`.

If the job fails:
- **Coding defect** (an R error, a wrong column, a wrong type): fix the function and run again.
- **Gate failure** (any size, power, symmetry, floor or runtime line): change nothing and report BLOCKED with the full log. The controller rules. Two ruling options are prepared but not pre-approved:
  - for a few-genes size failure, raise `min_genes` to the smallest tested G that passes, and treat smaller samples with the cohort value;
  - for a fragment-level power failure, record the measured comparison as the honest result and drop the claim from the prose.

Copy the log to `$SCR/logs/t1_final.log`.

- [ ] **Step 6: Document the measured behaviour in Step 14**

Append this bullet after the last Step 14 bullet ("**SNPs sharing read pairs.** ..."). Fill every `<...>` from `$SCR/logs/t1_final.log` by copying the printed values; nothing is estimated by hand.

"- **phASER gene test (Rmd 05): `hap_gene_test`.**
  - **Input:** one row per sample and gene with phASER's haplotype counts (`aCount`, `bCount` of `phaser_gene_ae`). phASER counts a read pair once per haplotype block, so the pooled counts are not inflated by SNPs that share reads.
  - **Test:** each row is tested against 0.5 with `bb_pvalue`. The test is symmetric in the two haplotypes: without phased genotypes the labels are arbitrary per gene and sample. Swapping the labels in any row changes no p-value and no dispersion; the unit tests check this.
  - **Dispersion:** `rho_own` is the trimmed central-region fit (`bb_estimate_rho_trim`) over the sample's genes. `rho_cohort` is the median of `rho_own` over the samples. `rho_used = max(rho_cohort, rho_own, RHO_MIN)`. A sample with fewer than 20 genes at `MIN_DEPTH` is tested with `max(rho_cohort, RHO_MIN)`.
  - **Null sizes** (measured; 10 seeds; cohorts of 6 samples × 2000 genes):
    - true rho 0.005 / 0.02 / 0.05 without imbalance: `<sizes>`;
    - with 20 % imbalanced genes: `<sizes>`;
    - with 40 %: `<sizes>`;
    - few genes (6 samples × 20 / 30 / 60): `<sizes>`;
    - one sample with 15 genes: `<size_s1>`;
    - heterogeneous cohort: `<size>` overall and `<size_s1>` for the noisiest sample.
  - **Power:** for 0.75/0.25 genes at rho 0.02 the power is `<power>`, against an oracle with the true rho of `<oracle_power>`.
  - **Comparison with the unphased path** (fragment-level simulation in which read pairs are shared between SNPs; genes with 4 or more SNPs within 600 transcript bases):
    - haplotype-test power `<hap_power_k4, lambda 300>` against `<acat_power_k4, lambda 300>` for the per-SNP tests combined with ACAT (Rmd 02) at 300 fragments per gene, and `<...>` against `<...>` at 40;
    - at 40 fragments per gene, `<hap_tested_k8>` of 8-SNP genes are tested by the haplotype test and `<acat_tested_k8>` by the unphased path. That path needs one SNP with at least `MIN_DEPTH` reads.
    - The haplotype test's null size is `<hap_size>` (300) and `<hap_size>` (40).
  - **Runtime:** `<s>` s and `<MB>` MB for one sample of 20,000 genes on one core."

- [ ] **Step 7: Checker lines, proven RED, then GREEN**

Append to `ase-pipeline/tests/check_skill.sh` before the final `[ $fail -eq 0 ]` line. Also delete the Stage 2 line `need "pasted verbatim into Rmd 02, Rmd 03 and Rmd 04"`: the intro sentence now names Rmd 05, and the new line below replaces it.

```bash
# --- Stage 3 statistics (hap_gene_test); each string is absent from the Stage 2 skill (de43444)
need "hap_gene_test <- function(a, b, sample, rho_min, min_genes = 20)"
need "per sample rho_used = max(rho_cohort, rho_own, rho_min)"
need "pasted verbatim into Rmd 02, Rmd 03, Rmd 04 and Rmd 05"
need "**phASER gene test (Rmd 05): \`hap_gene_test\`.**"
[ -s "$HERE/r/test_ase_stats_stage3.R" ] || { echo "FAIL: missing tests/r/test_ase_stats_stage3.R"; fail=1; }
grep -qF 'test_ase_stats_stage3.R' "$HERE/r/run_stats_tests.sh" || { echo "FAIL: run_stats_tests.sh does not run the Stage 3 tests"; fail=1; }
# --- end Stage 3 statistics
```

Proof:
- Run the checker on `git show HEAD:ase-pipeline/ase-pipeline.md > $SCR/pre.md`. Expected: 4 new `FAIL` lines (the two file checks pass, because the files exist by then).
- Run it on the working copy. Expected: `PASS`.

- [ ] **Step 8: Run the checker on the working copy**

Run: `bash ase-pipeline/tests/check_skill.sh ase-pipeline/ase-pipeline.md`
Expected: `PASS`.

- [ ] **Step 9: Amend the spec**

In `docs/superpowers/specs/2026-09-29-ase-pipeline-design.md`:

- Replace "**outbred with phASER** — the gene-level haplotypic counts replace pooled SNP counts." with "**outbred with phASER** (Stage 3) — the gene-level haplotype counts are tested in Rmd 05 and reported next to the unphased ACAT result, which they do not replace."
- Replace the paragraph that starts `**Stage 3 — phASER (outbred only, gated).**` with:

  "**Stage 3 — phASER (outbred only).** The feasibility gate (2026-10-01, `.superpowers/sdd/stage3-phaser-gate-report.md`) passed with constraints.
  - **Installation:** a one-time job clones https://github.com/secastel/phaser at the pinned commit `aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301`. It creates a conda environment under a user-chosen directory, from a skill-owned environment file that pins the versions resolved and tested in the gate (Python 3.14.7). The bioconda `phaser` package (Python 2.7) is never used, and every phASER script sets `PYTHONNOUSERSITE=1`.
  - **Per-sample array:** a job runs `phaser.py` on the WASP-filtered, duplicate-marked BAM with the heterozygous VCF that STAR used (`--paired_end 1 --mapq 255 --baseq MIN_BASEQ --pass_only 0 --unique_ids 1`; `--gw_phase_vcf 1` only for phased genotypes). It then runs `phaser_gene_ae` with a gene-span BED (0-based, one row per gene) built from the GTF exons.
  - **Rmd 05** tests each sample's gene-level haplotype counts against 0.5 with a beta-binomial test that is symmetric in the two haplotypes. Its dispersion is the trimmed central-region fit per sample: the larger of the sample's own value and the cohort median, floored at `RHO_MIN`. It reports the result next to the unphased ACAT result of Rmd 02, never instead of it.
  - **Direction:** without phased genotypes the haplotype labels are arbitrary, so no direction is reported, and differential ASE on phASER counts is not offered."

- [ ] **Step 10: Commit**

```bash
git status && git diff --stat
git add ase-pipeline/ase-pipeline.md ase-pipeline/tests/r/test_ase_stats_stage3.R ase-pipeline/tests/r/run_stats_tests.sh ase-pipeline/tests/check_skill.sh docs/superpowers/specs/2026-09-29-ase-pipeline-design.md
git commit -m "ase-pipeline: Stage 3 statistics (symmetric haplotype gene test, trimmed dispersion with cohort maximum rule) with simulation tests

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Stage 3 synthetic ground truth (phase, low-depth and two-block genes)

**Files:**
- Create:
  - `ase-pipeline/tests/synthetic/simulate_ase_stage3.R`;
  - `ase-pipeline/tests/synthetic/check_simulation_stage3.sh`;
  - `ase-pipeline/tests/synthetic/prove_checks_stage3.sh`;
  - `ase-pipeline/tests/synthetic/check_phase_reads_stage3.R`.
- Do not modify `simulate_ase_data.R` or `simulate_ase_stage2.R`: their data must stay byte-identical for their seeds.

**Interfaces:**
- Produces: `Rscript simulate_ase_stage3.R <outdir> <seed>`. Tasks 4, 5 and 8 depend on the names below. `O` = `<outdir>/outbred_phase`.

  | Path | Content |
  |---|---|
  | `<outdir>/seed.txt` | the seed |
  | `<outdir>/genome/genome.fa`, `genome.gtf` | chr1, 1,500,000 bp, 80 plus-strand genes `SYN000001`..`SYN000080` of 4-5 exons of 150-300 bp |
  | `O/samples.csv` | `sample,fastq_1,fastq_2,condition,cross_direction,individual`; 6 rows `obp_s1`..`obp_s6`, condition `ctrl`, cross_direction `NA`, individual `ind1`..`ind6` (sample k = individual k); absolute FASTQ paths |
  | `O/{sample}_1.fastq.gz`, `_2.fastq.gz` | paired reads, length 100 |
  | `O/{ind}.all.vcf`, `{ind}.het.vcf` | unphased GT (`0/1`, `0/0`, `1/1`); FILTER and ID `.`; the het file has exactly the `0/1` rows |
  | `O/{ind}.all.phased.vcf`, `{ind}.het.phased.vcf` | the same sites; het GT `0\|1` when the ALT allele is on haplotype 2 and `1\|0` when on haplotype 1 (left allele = haplotype 1); homozygous `0\|0` / `1\|1`; the het file has exactly the `0\|1` / `1\|0` rows |
  | `O/truth_genes.tsv` | `gene_id class h_abs n_snps lambda layout` |
  | `O/truth_snps.tsv` | `chrom pos ref alt gene_id tx_off cluster`, sorted by `pos` (cluster 1 or 2 for `two_block` genes, 0 otherwise) |
  | `O/truth_individual_genes.tsv` | `individual gene_id h1 n_het` (h1 = expression fraction of haplotype 1) |
  | `O/truth_phase.tsv` | `individual chrom pos gene_id gt alt_on` (`alt_on` 1 or 2 for `0/1`, `NA` otherwise) |
  | `O/truth_fragments.tsv` | `sample gene_id frag_h1 frag_h2 cover_h1 cover_h2 max_snp_depth` (`cover_*` = fragments of that haplotype covering at least one het SNP of the individual in the gene; `max_snp_depth` = the largest number of fragments covering one het SNP, 0 when none) |
  | `O/direct_counts/{sample}.table`, `.unfiltered.table`, `.wasp_stats.tsv` | as in Stage 2 (ASEReadCounter columns; only the individual's het SNPs with total ≥ 1; the unfiltered table is a copy; the WASP stats rows are `vW vA n` / `1 1 <REF fragments>` / `1 2 <ALT fragments>` / `none - 0`) |
  | `O/direct_gene_ae_phased/{sample}.gene_ae.txt` | phaser_gene_ae format (header `contig start stop name aCount bCount totalCount log2_aFC n_variants variants gw_phased bam`), one row per gene: `aCount = cover_h1`, `bCount = cover_h2`, `gw_phased` 1; a gene without coverage gets `0 0 0 inf 0 <empty> 1` |
  | `O/direct_gene_ae_unphased/{sample}.gene_ae.txt` | the same format, emulating unphased phaser_gene_ae: only the most-covered block of the gene, with A/B randomly swapped per sample and gene; `gw_phased` 0 (uncovered genes as above, `gw_phased` 1) |
  | `O/direct_flips.tsv` | `sample gene_id n_blocks best_block_n_snps flipped` (`flipped` TRUE/FALSE, `NA` when not covered) |

  The direct tables ignore sequencing errors. They test the Rmds, not the aligner or phASER.

- [ ] **Step 1: Write the check script first**

Create `ase-pipeline/tests/synthetic/check_simulation_stage3.sh`:

```bash
#!/bin/bash
# Usage: check_simulation_stage3.sh <outdir>   (after simulate_ase_stage3.R): exact quotas, layouts, phase truth, VCFs, direct tables
set -u; D=${1:?usage: check_simulation_stage3.sh <outdir>}; fail=0
bad() { echo "FAIL: $*"; fail=1; }
chk() { [ -s "$1" ] || bad "missing/empty $1"; }
export LC_ALL=C
O=$D/outbred_phase
for f in $D/seed.txt $D/genome/genome.fa $D/genome/genome.gtf $O/samples.csv $O/truth_genes.tsv $O/truth_snps.tsv \
         $O/truth_individual_genes.tsv $O/truth_phase.tsv $O/truth_fragments.tsv $O/direct_flips.tsv; do chk $f; done
[ "$(head -1 $O/samples.csv 2>/dev/null)" = "sample,fastq_1,fastq_2,condition,cross_direction,individual" ] || bad "samples.csv header"
[ "$(tail -n +2 $O/samples.csv 2>/dev/null | awk -F, '{print $1":"$4":"$5":"$6}' | paste -sd' ')" = \
  "obp_s1:ctrl:NA:ind1 obp_s2:ctrl:NA:ind2 obp_s3:ctrl:NA:ind3 obp_s4:ctrl:NA:ind4 obp_s5:ctrl:NA:ind5 obp_s6:ctrl:NA:ind6" ] || bad "samples.csv design"
# classes: exact quota with h_abs, n_snps, lambda and layout
keys=$(awk -F'\t' 'NR>1{print $2":"$3":"$4":"$5":"$6}' $O/truth_genes.tsv | sort | uniq -c | awk '{print $2"="$1}' | paste -sd' ')
want="hap_lowdepth:0.4:8:30:spread=3 hap_moderate:0.17:4:450:cluster=4 hap_strong:0.3:4:450:cluster=4 null:0:4:450:random=61 null_linked:0:4:450:cluster=6 two_block:0.3:6:450:two_cluster=2"
[ "$keys" = "$want" ] || bad "truth_genes classes: got '$keys'"
# SNP layouts in transcript offsets: cluster within 300 bases; spread every 120; two clusters >= 800 apart, cluster 1 within the
# first 150 bases; random >= 40 apart
awk -F'\t' 'NR==FNR{ if (FNR>1) { lay[$1]=$6; ns[$1]=$4 } next }
  FNR>1{ g=$5; k=++n[g]; o[g,k]=$6; c[g,k]=$7 }
  END{ for (g in lay) {
         if (n[g]!=ns[g]) { print "  " g ": " n[g] " SNPs, truth says " ns[g]; e=1; continue }
         span=o[g,n[g]]-o[g,1]
         if (lay[g]=="cluster" && span>=300) { print "  " g ": cluster span " span; e=1 }
         if (lay[g]=="spread") for (i=2;i<=n[g];i++) if (o[g,i]-o[g,i-1]!=120) { print "  " g ": spread spacing"; e=1 }
         if (lay[g]=="two_cluster" && !(c[g,1]==1 && c[g,3]==1 && c[g,4]==2 && c[g,6]==2 && o[g,3]<=150 && o[g,4]-o[g,3]>=800)) { print "  " g ": two clusters"; e=1 }
         if (lay[g]!="spread") for (i=2;i<=n[g];i++) if (o[g,i]-o[g,i-1]<40) { print "  " g ": SNPs closer than 40"; e=1 }
       } exit e }' $O/truth_genes.tsv $O/truth_snps.tsv || bad "SNP layouts"
# REF = genome base
awk 'NR==FNR{ if(/^>/) next; g=g $0; next } FNR>1{ if (substr(g,$2,1)!=$3) e=1 } END{exit e}' $D/genome/genome.fa $O/truth_snps.tsv || bad "truth_snps REF differs from the genome"
# phase truth: planted and null_linked genes heterozygous at every SNP in every individual; alt_on exactly for 0/1;
# per individual at least 2 cis and 2 trans genes among the cluster genes (hap_strong, hap_moderate, null_linked)
awk -F'\t' 'NR==FNR{ if (FNR>1) cl[$1]=$2; next }
  FNR>1{ if (($5=="0/1") != ($6=="1" || $6=="2")) e=1
         if (cl[$4]!="null" && $5!="0/1") e=1
         if (cl[$4]=="hap_strong" || cl[$4]=="hap_moderate" || cl[$4]=="null_linked") { k=$1 SUBSEP $4; a[k]=a[k] $6 } }
  END{ for (k in a) { split(k, p, SUBSEP); if (a[k] ~ /1/ && a[k] ~ /2/) tr[p[1]]++; else ci[p[1]]++ }
       for (i=1;i<=6;i++) { id="ind" i; if (tr[id]<2 || ci[id]<2) e=1 } exit e }' $O/truth_genes.tsv $O/truth_phase.tsv \
  || bad "truth_phase (planted genes all het, alt_on only for 0/1, >= 2 cis and >= 2 trans cluster genes per individual)"
# VCFs against the truth
vcf_rows() { grep -v '^#' "$1" | awk -F'\t' '{print $1"\t"$2"\t"$4"\t"$5}' | sort; }
truth_rows() { tail -n +2 "$1" | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4}' | sort; }
for i in 1 2 3 4 5 6; do
  I=ind$i
  for f in all het all.phased het.phased; do chk $O/$I.$f.vcf; done
  diff <(vcf_rows $O/$I.all.vcf) <(truth_rows $O/truth_snps.tsv) >/dev/null || bad "$I.all.vcf sites differ from truth_snps.tsv"
  diff <(vcf_rows $O/$I.all.phased.vcf) <(truth_rows $O/truth_snps.tsv) >/dev/null || bad "$I.all.phased.vcf sites differ from truth_snps.tsv"
  awk -F'\t' -v ind=$I 'NR==FNR{ if ($1==ind) g[$3]=$5; next } !/^#/{ if (g[$2]!=$10) e=1 } END{exit e}' $O/truth_phase.tsv $O/$I.all.vcf \
    || bad "$I.all.vcf GT differs from truth_phase"
  awk -F'\t' -v ind=$I 'NR==FNR{ if ($1==ind) { g[$3]=$5; a[$3]=$6 } next }
    !/^#/{ w = (g[$2]=="0/0") ? "0|0" : (g[$2]=="1/1") ? "1|1" : (a[$2]=="2") ? "0|1" : "1|0"; if ($10!=w) e=1 } END{exit e}' \
    $O/truth_phase.tsv $O/$I.all.phased.vcf || bad "$I.all.phased.vcf GT differs from truth_phase (left allele = haplotype 1)"
  diff <(grep -v '^#' $O/$I.all.vcf | awk -F'\t' '$10=="0/1"') <(grep -v '^#' $O/$I.het.vcf) >/dev/null || bad "$I.het.vcf is not the 0/1 rows of $I.all.vcf"
  diff <(grep -v '^#' $O/$I.all.phased.vcf | awk -F'\t' '$10=="0|1" || $10=="1|0"') <(grep -v '^#' $O/$I.het.phased.vcf) >/dev/null \
    || bad "$I.het.phased.vcf is not the het rows of $I.all.phased.vcf"
done
# reads, direct counts and direct gene tables per sample
for s in $(tail -n +2 $O/samples.csv | cut -d, -f1); do
  ind=$(awk -F, -v s=$s '$1==s{print $6}' $O/samples.csv)
  for f in ${s}_1.fastq.gz ${s}_2.fastq.gz direct_counts/$s.table direct_counts/$s.unfiltered.table direct_counts/$s.wasp_stats.tsv \
           direct_gene_ae_phased/$s.gene_ae.txt direct_gene_ae_unphased/$s.gene_ae.txt; do chk $O/$f; done
  n1=$(zcat $O/${s}_1.fastq.gz 2>/dev/null | awk 'END{print NR/4}'); n2=$(zcat $O/${s}_2.fastq.gz 2>/dev/null | awk 'END{print NR/4}')
  nf=$(awk -F'\t' -v s=$s '$1==s{t+=$3+$4} END{print t+0}' $O/truth_fragments.tsv)
  [ "$n1" = "$n2" ] && [ "$n1" = "$nf" ] || bad "$s read pairs ($n1 / $n2) differ from truth_fragments ($nf)"
  awk -F'\t' 'NR==1{ if ($1!="contig"||$6!="refCount"||$7!="altCount"||$8!="totalCount") exit 1; next } { if ($6+$7!=$8 || $11!=$8 || $8<1) exit 1 }' \
    $O/direct_counts/$s.table || bad "$s direct counts header or ref+alt/total/rawDepth"
  awk -F'\t' 'NR==FNR{ if (!/^#/ && $10=="0/1") h[$2]=1; next } FNR>1 && !($2 in h){e=1} END{exit e}' $O/$ind.het.vcf $O/direct_counts/$s.table \
    || bad "$s direct counts at positions that are not het in $ind"
  cmp -s $O/direct_counts/$s.table $O/direct_counts/$s.unfiltered.table || bad "$s unfiltered.table must equal the table"
  awk -F'\t' -v s=$s 'NR==FNR{ if ($1==s) { c1[$2]=$5; c2[$2]=$6 } next }
    FNR==1{ if ($4!="name"||$5!="aCount"||$6!="bCount"||$7!="totalCount"||$9!="n_variants"||$11!="gw_phased"||$12!="bam") e=1; next }
    { n++; if ($5!=c1[$4] || $6!=c2[$4] || $7!=$5+$6 || $11!=1 || $12!=s) e=1; if ($7==0 && $8!="inf") e=1 } END{ if (n!=80) e=1; exit e }' \
    $O/truth_fragments.tsv $O/direct_gene_ae_phased/$s.gene_ae.txt || bad "$s phased direct gene table differs from the truth cover counts"
  awk -F'\t' -v s=$s 'FILENAME ~ /truth_genes/ { if (FNR>1) cl[$1]=$2; next }
    FILENAME ~ /truth_fragments/ { if ($1==s) { c1[$2]=$5; c2[$2]=$6 } next }
    FILENAME ~ /direct_flips/ { if ($1==s) nb[$2]=$3; next }
    FNR>1 { n++; t=c1[$4]+c2[$4]; if ($7>t || $5+$6!=$7) e=1
            if (nb[$4]==1 && !(($5==c1[$4] && $6==c2[$4]) || ($5==c2[$4] && $6==c1[$4]))) e=1
            if (cl[$4]=="two_block" && $9>3) e=1
            if ($7>0 && $11!=0) e=1 } END{ if (n!=80) e=1; exit e }' \
    $O/truth_genes.tsv $O/truth_fragments.tsv $O/direct_flips.tsv $O/direct_gene_ae_unphased/$s.gene_ae.txt || bad "$s unphased direct gene table"
done
# hap_lowdepth cells: every SNP below depth 10 and at least 18 covering fragments, in every sample (6 x 3 = 18 cells);
# and the direct counts at those SNPs agree (total < 10)
awk -F'\t' 'NR==FNR{ if (FNR>1) cl[$1]=$2; next } FNR>1 && cl[$2]=="hap_lowdepth" { k++; if ($7>=10 || $5+$6<18) e=1 } END{ if (k!=18) e=1; exit e }' \
  $O/truth_genes.tsv $O/truth_fragments.tsv || bad "hap_lowdepth cells (max SNP depth < 10, >= 18 covering fragments, 18 cells)"
awk -F'\t' 'FILENAME ~ /truth_genes/ { if (FNR>1 && $2=="hap_lowdepth") lg[$1]=1; next } FILENAME ~ /truth_snps/ { if (FNR>1 && ($5 in lg)) lp[$2]=1; next }
  FNR>1 && ($2 in lp) && $8>=10 { e=1 } END{exit e}' $O/truth_genes.tsv $O/truth_snps.tsv $O/direct_counts/*.table || bad "direct counts: a hap_lowdepth SNP reaches depth 10"
# planted cells: realised haplotype-1 fractions follow h1 (z within +-5; sum of z^2 over the 78 planted cells <= 130)
awk -F'\t' 'FILENAME ~ /truth_genes/ { if (FNR>1) cl[$1]=$2; next } FILENAME ~ /truth_individual_genes/ { if (FNR>1) h[$1 SUBSEP $2]=$3; next }
  FNR>1 && cl[$2]!="null" && cl[$2]!="null_linked" { ind="ind" substr($1, 6); n=$5+$6; p=h[ind SUBSEP $2]
     if (n<1) { e=1; next } z=($5-n*p)/sqrt(n*p*(1-p)*(1+(n-1)*0.005)); k++; s2+=z*z; if (z>5 || z<-5) { print "  " $1, $2, z; e=1 } }
  END{ if (k!=78 || s2>130) { print "  cells " k ", sum z^2 " s2; e=1 } exit e }' \
  $O/truth_genes.tsv $O/truth_individual_genes.tsv $O/truth_fragments.tsv || bad "planted cells: realised haplotype-1 fractions do not follow h1"
# planted h1: |h1 - 0.5| = h_abs, both signs present per class; nulls exactly 0.5
awk -F'\t' 'NR==FNR{ if (FNR>1) { ha[$1]=$3; cl[$1]=$2 } next } FNR>1{ d=$3-0.5; a=(d<0)?-d:d; if ((a-ha[$2])^2>1e-10) e=1; if (d>0) up[cl[$2]]++; if (d<0) dn[cl[$2]]++ }
  END{ split("hap_strong hap_moderate hap_lowdepth two_block", c, " "); for (i in c) if (!up[c[i]] || !dn[c[i]]) e=1; exit e }' \
  $O/truth_genes.tsv $O/truth_individual_genes.tsv || bad "truth_individual_genes h1 (|h1 - 0.5| = h_abs, both signs in every planted class)"
# the unphased emulation flips labels both ways (>= 10 each over all covered cells)
awk -F'\t' 'NR>1 && $5=="TRUE"{t++} NR>1 && $5=="FALSE"{f++} END{exit !(t>=10 && f>=10)}' $O/direct_flips.tsv || bad "direct_flips: fewer than 10 flipped or unflipped cells"
[ $fail -eq 0 ] && echo PASS || exit 1
```

Run it on an empty directory: `mkdir -p $SCR/empty && bash ase-pipeline/tests/synthetic/check_simulation_stage3.sh $SCR/empty`. Expected: `FAIL` lines and exit 1.

- [ ] **Step 2: Write the generator**

`simulate_ase_stage3.R` (R with `Biostrings` and base R; about 330 lines). It is seeded by its second argument, which it writes to `seed.txt`, and follows this algorithm.

1. **Constants and genome:**
   - `READLEN <- 100L; ERR <- 0.002; PHI_BIO <- 0.005; NGENES <- 80L; GLEN <- 1500000L; MIN_DEPTH <- 10L`.
   - The genome and genes come from the Stage 2 block of `simulate_ase_stage2.R` (from `# ---- genome and genes` up to the line `stopifnot(all(vapply(genes, function(x) x$len, 1) >= 600))`), copied with these constants.
2. **Copied functions:** copy verbatim from `simulate_ase_stage2.R` the functions `alt_of`, `spliced`, `revcomp`, `add_err`, `write_fq`, `vcf_hdr`, `frags2` and `wtab`. Also copy `haplo`, with its argument changed to a full allele vector per gene.
3. **Classes:**
   - Shuffle the gene ids with `sample()`.
   - From the genes with spliced length ≥ 1100, take the first 2 as `two_block`. From the remaining genes with length ≥ 1000, take 3 as `hap_lowdepth`.
   - From the rest, in shuffled order: 4 `hap_strong`, 4 `hap_moderate`, 6 `null_linked`. All remaining genes are `null`.
   - Stop with an error if fewer genes than needed are long enough.

   | class | h_abs | n_snps | lambda | layout |
   |---|---|---|---|---|
   | hap_strong | 0.3 | 4 | 450 | cluster |
   | hap_moderate | 0.17 | 4 | 450 | cluster |
   | hap_lowdepth | 0.4 | 8 | 30 | spread |
   | two_block | 0.3 | 6 | 450 | two_cluster |
   | null_linked | 0 | 4 | 450 | cluster |
   | null | 0 | 4 | 450 | random |

   Write `h_abs` with `format(x)`, so the values are `0.3`, `0.17`, `0.4` and `0`.
4. **SNP offsets** (transcript offsets; plus strand):
   - `cluster`: 4 offsets within 300 consecutive bases, pairwise ≥ 40 apart. Use `repeat` with `w0 <- sample.int(len - 299, 1)`.
   - `spread`: `o0 + 120 * (0:7)` with `o0 <- sample(61:(len - 900), 1)`.
   - `two_cluster`: cluster 1 is 3 offsets in `1:150`, pairwise ≥ 40 apart; cluster 2 is 3 offsets in `(len - 149):len`, pairwise ≥ 40 apart. Write `cluster` 1 and 2 in `truth_snps.tsv`; every other gene gets 0.
   - `random`: 4 offsets in `1:len`, pairwise ≥ 40 apart.

   Map offsets to the genome as in Stage 2 (`cum`). The REF base is the genome base; the ALT base is `alt_of(REF)`. Sort `truth_snps.tsv` by `pos`.
5. **Individuals `ind1`..`ind6`:**
   - **Genotypes:** every SNP of the planted classes and of `null_linked` is `0/1`. A `null` SNP draws its genotype from `sample(c("0/1", "0/0", "1/1"), prob = c(0.6, 0.2, 0.2))`.
   - **`alt_on` for `0/1` SNPs:**
     - In `hap_strong`, `hap_moderate` and `null_linked` genes: with probability 0.5 the gene is *cis* for that individual (every `alt_on` equal to one value drawn from 1 or 2); otherwise it is *trans* (`alt_on` drawn per SNP and redrawn until both 1 and 2 occur).
     - In all other genes, `alt_on` is drawn per SNP from 1 or 2.
     - After all genes are drawn, `stopifnot` that every individual has at least 2 cis and 2 trans cluster genes.
   - **h1:** for the planted classes, `h1 <- 0.5 + h_abs * sample(c(-1, 1), 1)` per individual and gene; for `null` and `null_linked`, `h1 <- 0.5`. Write `truth_individual_genes.tsv` with `n_het` = the number of `0/1` SNPs.
   - **Haplotypes:** haplotype 1 carries the ALT base at a het SNP when `alt_on == 1` and the REF base when `alt_on == 2`; haplotype 2 is the opposite. Homozygous SNPs get their allele on both haplotypes. Write `truth_phase.tsv`.
   - **VCFs:** write the four VCFs of the interface table with `vcf_hdr(ind)`, FILTER `.` and ID `.`.
6. **Fragments and direct counts per sample** (sample k = individual k), for each gene:
   - Draw `ps <- rbeta(1, h1 * (1 - PHI_BIO) / PHI_BIO, (1 - h1) * (1 - PHI_BIO) / PHI_BIO)`, `n <- rpois(1, lambda)` and `n1 <- rbinom(1, n, ps)`.
   - Fragments: `f1 <- frags2(h1seq, n1, offs)` and `f2 <- frags2(h2seq, n - n1, offs)`.
   - **`hap_lowdepth` genes:** with `n` and `n1` fixed, redraw `f1` and `f2` (fragment positions only, never the alleles) until two conditions hold: every SNP has `colSums(f1$cover) + colSums(f2$cover) < MIN_DEPTH`, and at least 18 fragments cover one or more SNPs. Stop with an error after 1000 attempts.
   - **Direct counts** at het SNPs: `refCount` = coverage from the haplotype carrying REF (`f2` when `alt_on == 1`, `f1` when `alt_on == 2`); `altCount` = coverage from the other haplotype.
   - **`truth_fragments.tsv`:** `frag_h1 = n1`, `frag_h2 = n - n1`; `cover_h1` / `cover_h2` = fragments covering at least one het SNP; `max_snp_depth` = the maximum over het SNPs of the total coverage.
   - **Reads:** written through `add_err` and `write_fq` as in Stage 2.
   - **WASP stats and unfiltered copy:** as in Stage 2, with REF fragments = the sum of `refCount` over the sample's rows and ALT fragments = the sum of `altCount`.
7. **Emulated phASER gene tables.**
   - **Phased** (`direct_gene_ae_phased`):
     - `aCount = cover_h1`, `bCount = cover_h2`, `totalCount` = their sum;
     - `log2_aFC` = `log2(a / b)` written with `format(x, digits = 6)`, and the strings `inf` / `-inf` when b or a is 0;
     - `n_variants` = het SNPs with coverage; `variants` = `chr1_<pos>_<ref>_<alt>` joined by `,`;
     - `gw_phased` 1; `contig` `chr1`; `start` = gene start − 1; `stop` = gene end; `bam` = the sample.
     - An uncovered gene is written as `0 0 0 inf 0 <empty> 1`.
   - **Unphased** (`direct_gene_ae_unphased`): blocks of the covered het SNPs come from this function (copy it verbatim). The best block's counts are written, swapped when `flipped` (`flipped <- runif(1) < 0.5`, drawn per covered sample and gene), with `gw_phased` 0. Uncovered genes are written as in the phased table. Write `direct_flips.tsv`.

   ```r
   emulate_blocks <- function(cover1, cover2, het) {   # fragment x SNP logical matrices of haplotype 1 / 2; het: logical over the gene's SNPs
     cv <- rbind(cover1, cover2)[, het, drop = FALSE]; hs <- c(rep(1L, nrow(cover1)), rep(2L, nrow(cover2)))
     cov_snp <- which(colSums(cv) > 0); if (length(cov_snp) == 0) return(NULL)
     comp <- seq_along(cov_snp)                        # union-find: SNPs joined when one fragment covers both
     find <- function(i) { while (comp[i] != i) i <- comp[i]; i }
     sub <- cv[, cov_snp, drop = FALSE]
     for (f in which(rowSums(sub) >= 2)) { s <- which(sub[f, ]); r <- find(s[1]); for (x in s[-1]) comp[find(x)] <- r }
     root <- vapply(seq_along(cov_snp), find, integer(1))
     blocks <- lapply(split(seq_along(cov_snp), root), function(ix) {
       anyc <- rowSums(sub[, ix, drop = FALSE]) > 0
       list(snps = which(het)[cov_snp[ix]], a = sum(anyc & hs == 1L), b = sum(anyc & hs == 2L)) })
     tot <- vapply(blocks, function(x) x$a + x$b, numeric(1)); first <- vapply(blocks, function(x) min(x$snps), numeric(1))
     list(n_blocks = length(blocks), best = blocks[[order(-tot, first)[1]]])
   }
   ```
8. **Sample sheet:** written with absolute FASTQ paths (`normalizePath(outdir)`).

- [ ] **Step 3: Write the read-level phase check**

Create `ase-pipeline/tests/synthetic/check_phase_reads_stage3.R`:

```r
#!/usr/bin/env Rscript
# Usage: Rscript check_phase_reads_stage3.R <outdir>
# From the FASTQs alone: every read pair that carries alleles at two or more heterozygous SNPs of one gene (exact 41 bp allele
# windows, one window per SNP and allele) must carry alleles of ONE haplotype of truth_phase.tsv (cis or trans as planted).
# Pass: in every sample at least 98% of such pairs agree with the truth, and at least 50 pairs link SNPs in trans.
suppressPackageStartupMessages(library(Biostrings))
O <- file.path(commandArgs(trailingOnly = TRUE)[1], "outbred_phase"); W <- 20L
gseq <- as.character(readDNAStringSet(file.path(O, "../genome/genome.fa"))[[1]])
ex <- read.delim(file.path(O, "../genome/genome.gtf"), header = FALSE, quote = "", stringsAsFactors = FALSE)
ex <- ex[ex$V3 == "exon", ]; ex$gene_id <- sub('.*gene_id "([^"]+)".*', "\\1", ex$V9)
ts <- read.delim(file.path(O, "truth_snps.tsv"), stringsAsFactors = FALSE)
tp <- read.delim(file.path(O, "truth_phase.tsv"), stringsAsFactors = FALSE, colClasses = c(gt = "character", alt_on = "character"))
sm <- read.csv(file.path(O, "samples.csv"), stringsAsFactors = FALSE)
spl <- lapply(split(ex, ex$gene_id), function(e) { e <- e[order(e$V4), ]; strsplit(paste(substring(gseq, e$V4, e$V5), collapse = ""), "")[[1]] })
allok <- TRUE
for (i in seq_len(nrow(sm))) {
  ph <- tp[tp$individual == sm$individual[i] & tp$gt == "0/1", c("pos", "alt_on")]
  s <- merge(ts, ph, by = "pos")
  w <- do.call(rbind, lapply(seq_len(nrow(s)), function(k) {
    sp <- spl[[s$gene_id[k]]]; o <- s$tx_off[k]
    if (o <= W || o + W > length(sp)) return(NULL)
    hr <- sp; hr[o] <- s$ref[k]; ha <- sp; ha[o] <- s$alt[k]
    data.frame(snp = k, allele = c("ref", "alt"), seq = c(paste(hr[(o - W):(o + W)], collapse = ""), paste(ha[(o - W):(o + W)], collapse = "")),
               stringsAsFactors = FALSE)
  }))
  pd <- PDict(DNAStringSet(w$seq))
  h1 <- vwhichPDict(pd, readDNAStringSet(sm$fastq_1[i], format = "fastq"))
  h2 <- vwhichPDict(pd, reverseComplement(readDNAStringSet(sm$fastq_2[i], format = "fastq")))
  hits <- mapply(function(x, y) unique(c(x, y)), h1, h2, SIMPLIFY = FALSE)
  total <- 0; agree <- 0; trans <- 0
  for (q in which(lengths(hits) >= 2)) {
    hw <- w[hits[[q]], ]
    if (anyDuplicated(hw$snp) || length(unique(s$gene_id[hw$snp])) != 1) next   # both alleles of one SNP (error) or two genes
    hap <- ifelse(hw$allele == "alt", s$alt_on[hw$snp], ifelse(s$alt_on[hw$snp] == "1", "2", "1"))
    total <- total + 1; agree <- agree + (length(unique(hap)) == 1)
    if (length(unique(s$alt_on[hw$snp])) > 1) trans <- trans + 1
  }
  okk <- total > 0 && agree / total >= 0.98 && trans >= 50
  cat(sprintf("%s (%s): %d read pairs link two or more SNPs, %.4f agree with truth_phase, %d link SNPs in trans: %s\n",
              sm$sample[i], sm$individual[i], total, if (total > 0) agree / total else NA, trans, if (okk) "ok" else "FAIL"))
  allok <- allok && okk
}
cat(if (allok) "PHASE READS PASS\n" else "PHASE READS FAIL\n"); quit(status = if (allok) 0 else 1)
```

- [ ] **Step 4: Write the injected-defect proofs**

Create `ase-pipeline/tests/synthetic/prove_checks_stage3.sh`:

```bash
#!/bin/bash
# Usage: prove_checks_stage3.sh <good outdir> <scratch dir>: the check must PASS on the good data and FAIL on each injected defect
set -u; G=${1:?}; W=${2:?}; HERE=$(dirname "$0"); fail=0
bash $HERE/check_simulation_stage3.sh $G >/dev/null || { echo "FAIL: check does not pass on the good data"; exit 1; }
defect() {   # $1 = name, $2 = shell command run inside the copy ($C)
  rm -rf $W/copy; cp -r $G $W/copy; C=$W/copy/outbred_phase; eval "$2"
  if out=$(bash $HERE/check_simulation_stage3.sh $W/copy 2>&1); then echo "FAIL: defect '$1' not detected"; fail=1
  else echo "ok   defect '$1' detected: $(echo "$out" | grep "^FAIL" | cut -c1-140 | paste -sd"|")"; fi
}
defect "phased GT of one het site swapped"   'sed -i "0,/\t0|1\$/s//\t1|0/" $C/ind1.all.phased.vcf'
defect "class quota changed"                 'awk -F"\t" "BEGIN{OFS=FS} NR>1 && \$2==\"null\" && !x {\$2=\"null_linked\"; x=1} 1" $C/truth_genes.tsv > $C/t && mv $C/t $C/truth_genes.tsv'
defect "hap_lowdepth SNP depth 12"           'awk -F"\t" "BEGIN{OFS=FS} NR==FNR{ if (\$2==\"hap_lowdepth\") g[\$1]=1; next } FNR>1 && (\$2 in g) && !x {\$7=12; x=1} 1" $C/truth_genes.tsv $C/truth_fragments.tsv > $C/t && mv $C/t $C/truth_fragments.tsv'
defect "phased direct table a/b swapped"     'awk -F"\t" "BEGIN{OFS=FS} NR>1 && \$5!=\$6 && !x {t=\$5; \$5=\$6; \$6=t; x=1} 1" $C/direct_gene_ae_phased/obp_s1.gene_ae.txt > $C/t && mv $C/t $C/direct_gene_ae_phased/obp_s1.gene_ae.txt'
defect "two clusters moved together"         'awk -F"\t" "BEGIN{OFS=FS} NR==FNR{ if (\$2==\"two_block\" && !g) g=\$1; next } FNR>1 && \$5==g && \$7==2 && !x {\$6=300; x=1} 1" $C/truth_genes.tsv $C/truth_snps.tsv > $C/t && mv $C/t $C/truth_snps.tsv'
defect "one read pair dropped from mate 1"   'zcat $C/obp_s1_1.fastq.gz | tail -n +5 | gzip > $C/t && mv $C/t $C/obp_s1_1.fastq.gz'
defect "het.vcf misses a het row"            'awk "/^#/ || !x {if (!/^#/) {x=1; next}} 1" $C/ind1.het.vcf > $C/t && mv $C/t $C/ind1.het.vcf'
defect "seed.txt removed"                    'rm -f $W/copy/seed.txt'
defect "planted gene homozygous at one SNP"  'awk -F"\t" "BEGIN{OFS=FS} NR==FNR{ if (\$2==\"hap_strong\") g[\$1]=1; next } FNR>1 && (\$4 in g) && !x {\$5=\"0/0\"; \$6=\"NA\"; x=1} 1" $C/truth_genes.tsv $C/truth_phase.tsv > $C/t && mv $C/t $C/truth_phase.tsv'
defect "direct count ref+alt != total"       'awk -F"\t" "BEGIN{OFS=FS} NR==2{\$6=\$6+5} 1" $C/direct_counts/obp_s1.table > $C/t && mv $C/t $C/direct_counts/obp_s1.table'
defect "h1 flipped for one planted cell"     'awk -F"\t" "BEGIN{OFS=FS} NR==FNR{ if (\$2==\"hap_strong\") g[\$1]=1; next } FNR>1 && (\$2 in g) && !x {\$3=1-\$3; x=1} 1" $C/truth_genes.tsv $C/truth_individual_genes.tsv > $C/t && mv $C/t $C/truth_individual_genes.tsv'
defect "unphased one-block counts altered"   'awk -F"\t" "BEGIN{OFS=FS} NR==FNR{ if (\$1==\"obp_s1\" && \$3==1) g[\$2]=1; next } FNR>1 && (\$4 in g) && !x {\$5=\$5+3; \$7=\$7+3; x=1} 1" $C/direct_flips.tsv $C/direct_gene_ae_unphased/obp_s1.gene_ae.txt > $C/t && mv $C/t $C/direct_gene_ae_unphased/obp_s1.gene_ae.txt'
rm -rf $W/copy
[ $fail -eq 0 ] && echo "PROOFS PASS" || exit 1
```

- [ ] **Step 5: Generate, check, prove and test determinism on the cluster**

```bash
SCR=/net/bmc-lab3/data/bcc/yannvrb/ase_stage3_scratch; mkdir -p $SCR/logs
cat > $SCR/t2_gen.sh <<'EOF'
#!/bin/bash
#SBATCH -p bcc
#SBATCH -n 1 --mem=8G -t 1:00:00
set -uo pipefail
module add singularity/3.10.4 || exit 1
SIF=/net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif; T=ase-pipeline/tests/synthetic
cd /net/bmc-lab3/data/bcc/yannvrb/Github_repos/bioinformatics-claude-skills || exit 1
singularity exec --bind /net/bmc-lab3 $SIF Rscript $T/simulate_ase_stage3.R /net/bmc-lab3/data/bcc/yannvrb/ase_synthetic_stage3 20261001 || exit 1
singularity exec --bind /net/bmc-lab3 $SIF Rscript $T/simulate_ase_stage3.R /net/bmc-lab3/data/bcc/yannvrb/ase_stage3_scratch/regen 20261001 || exit 1
singularity exec --bind /net/bmc-lab3 $SIF Rscript $T/check_phase_reads_stage3.R /net/bmc-lab3/data/bcc/yannvrb/ase_synthetic_stage3 || exit 1
bash $T/check_simulation_stage3.sh /net/bmc-lab3/data/bcc/yannvrb/ase_synthetic_stage3 || exit 1
bash $T/prove_checks_stage3.sh /net/bmc-lab3/data/bcc/yannvrb/ase_synthetic_stage3 /net/bmc-lab3/data/bcc/yannvrb/ase_stage3_scratch/proofs || exit 1
EOF
sbatch -p bcc -o $SCR/logs/t2_gen_%j.out $SCR/t2_gen.sh
```

Expected log:
- `PHASE READS PASS`, with a per-sample line showing ≥ 0.98 agreement and ≥ 50 trans pairs;
- `PASS` from the check;
- 12 `ok` defect lines and `PROOFS PASS`.

Then check determinism: `diff <(cd /net/bmc-lab3/data/bcc/yannvrb/ase_synthetic_stage3 && find . -type f ! -name samples.csv -exec md5sum {} + | sort -k2) <(cd $SCR/regen && find . -type f ! -name samples.csv -exec md5sum {} + | sort -k2)` prints nothing.

Then prove the phase check: copy `ase_synthetic_stage3` to `$SCR/phase_bad`. In its `truth_phase.tsv`, swap `alt_on` 1↔2 at the first SNP of one `hap_strong` gene of `ind1` (one awk keyed on that gene and individual), then run `check_phase_reads_stage3.R` on it in a job. Expected: `PHASE READS FAIL`, because `obp_s1` falls below 0.98 agreement.

If the generator cannot satisfy a planted condition (the cis/trans quota or the `hap_lowdepth` redraw), fix the generator parameters only within the interface and the quota of Step 1. Never relax the check.

- [ ] **Step 6: Commit (scripts only, never the data)**

```bash
git status && git diff --stat
git add ase-pipeline/tests/synthetic/simulate_ase_stage3.R ase-pipeline/tests/synthetic/check_simulation_stage3.sh ase-pipeline/tests/synthetic/prove_checks_stage3.sh ase-pipeline/tests/synthetic/check_phase_reads_stage3.R
git commit -m "ase-pipeline: Stage 3 synthetic ground truth (cis/trans phase, low-depth and two-block genes, phased VCFs, emulated phASER tables) with checks and defect proofs

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 3: phASER installation job (Step 18, first part), fixtures and the phASER flag rule

**Files:**
- Modify `ase-pipeline/ase-pipeline.md`: new `## Step 18 — phASER: installation and per-sample haplotype counts (outbred, optional)`, between the end of Step 17 (its closing `---`) and `## Notes for the assistant`, followed by its own `---`. This task writes the introduction and the `setup_phaser_env.sh` section; Task 5 adds `phaser_count.sh`.
- Create:
  - `ase-pipeline/tests/synthetic/cut_block.sh`;
  - `ase-pipeline/tests/fixtures/phaser_help.txt`, `phaser_gene_ae_help.txt`, `phaser_env.txt`.
- Modify: `ase-pipeline/tests/fixtures/verification.md` (append a Stage 3 section), `ase-pipeline/tests/check_skill.sh`.

**Interfaces:**
- Produces:
  - `bash cut_block.sh <skill.md> <anchor>`: prints the body of the first fenced code block (an opening line `` ``` `` or `` ```bash `` / `` ```r ``, closing line `` ``` ``) that opens after the first line containing `<anchor>` (a fixed string). It exits 1 when there is none.
  - The skill block `setup_phaser_env.sh`. It is a complete script with placeholders `{PHASER_HOME}`, `{USER_EMAIL}` and `{RESULTS_DIR}`, and is cut by the anchor `` ### `setup_phaser_env.sh` ``.
  - The installation contract that Tasks 5 and 6 rely on:
    - after a successful run, `{PHASER_HOME}/install_ok.txt` exists and its first line is exactly `commit aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301`;
    - `{PHASER_HOME}/env/bin/python`, `{PHASER_HOME}/src/phaser/phaser/phaser.py` and `{PHASER_HOME}/src/phaser/phaser_gene_ae/phaser_gene_ae.py` exist.
  - The installation of the tests in `$SCR/phaser_home`.

- [ ] **Step 1: Failing checker lines**

Append to `check_skill.sh`:

```bash
# --- Stage 3 Task 3 (phASER installation, fixtures, flag rule); each check fails on the Task 2 skill
need "## Step 18 — phASER: installation and per-sample haplotype counts (outbred, optional)"
need "PHASER_COMMIT=aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301"
need "PHASER_REPO=https://github.com/secastel/phaser.git"
need 'conda env create -p "$PHASER_HOME/env" -f "$PHASER_HOME/environment.phaser.yml"'
need "  - python=3.14.7"
need "export PYTHONNOUSERSITE=1"
need "never the bioconda package"
need "**Bumping the pinned commit.**"
[ "$(grep -c '^module add miniconda3/v4' "$SKILL")" = 1 ] || { echo "FAIL: exactly one line may load miniconda3/v4 (setup_phaser_env.sh)"; fail=1; }
forbid "conda install -c bioconda phaser"
forbid "bioconda::phaser"
for f in phaser_help.txt phaser_gene_ae_help.txt phaser_env.txt; do [ -s "$FIX/$f" ] || { echo "FAIL: missing fixture $f"; fail=1; }; done
[ -s "$HERE/synthetic/cut_block.sh" ] || { echo "FAIL: missing tests/synthetic/cut_block.sh"; fail=1; }
# phASER flags: every --flag on a phaser.py / phaser_gene_ae.py command (continuation lines joined) is in the fixture of the pinned commit
ph_flags() { grep -oE -- '--[a-z][a-z0-9_]*' "$FIX/$1" 2>/dev/null | sed 's/^--//' | sort -u; }
ph_check() {   # $1 = script path fragment on the command line, $2 = fixture
  awk '{ if (sub(/\\$/, "")) { buf = buf $0 " "; next } print buf $0; buf = "" }' "$SKILL" | grep -F -- "$1" |
    grep -oE -- '(^|[ (])--[a-z][a-z0-9_]*' | sed -E 's/^[ (]?--//' | sort -u |
    while read -r f; do ph_flags "$2" | grep -qx -- "$f" || echo "FAIL: $1 flag not in fixtures/$2: --$f"; done
}
ph_out=$(ph_check "phaser/phaser.py" phaser_help.txt; ph_check "phaser_gene_ae/phaser_gene_ae.py" phaser_gene_ae_help.txt)
[ -z "$ph_out" ] || { echo "$ph_out"; fail=1; }
# --- end Stage 3 Task 3
```

Then extend the global flag check (section (a) at the top of the checker) so that phASER flags count as tool flags. Insert directly after the `tool_flags=$( ... )` assignment:

```bash
for f in phaser_help.txt phaser_gene_ae_help.txt; do
  [ -s "$FIX/$f" ] && tool_flags="$tool_flags
$(grep -oE -- '--[a-z][a-z0-9_]*' "$FIX/$f" | sed 's/^--//' | sort -u)"
done
```

Run the checker. Expected: FAIL. The fixtures and Step 18 do not exist yet.

- [ ] **Step 2: Write `cut_block.sh`**

Create `ase-pipeline/tests/synthetic/cut_block.sh`:

```bash
#!/bin/bash
# Usage: cut_block.sh <skill.md> <anchor>
# Prints the body of the first fenced code block (``` or ```<lang>, closed by ```) that opens after the first line containing
# <anchor> (fixed string). Exits 1 when the anchor or the block is missing. Four-backtick (````rmd) fences never match.
set -uo pipefail
SKILL=${1:?usage: cut_block.sh <skill.md> <anchor>}; ANCHOR=${2:?usage: cut_block.sh <skill.md> <anchor>}
awk -v a="$ANCHOR" '!seen && index($0, a) { seen = 1; next }
  seen && !inb && /^```[a-z]*$/ { inb = 1; next }
  inb && /^```$/ { found = 1; exit }
  inb { print }
  END { exit !found }' "$SKILL"
```

Self-test:
- `bash ase-pipeline/tests/synthetic/cut_block.sh ase-pipeline/ase-pipeline.md "### Shared block R"` prints exactly the 4 lines of block R (starting `FASTA="{FASTA_PATH}"`).
- `bash ase-pipeline/tests/synthetic/cut_block.sh ase-pipeline/ase-pipeline.md "no such anchor"; echo $?` prints `1`.

- [ ] **Step 3: Write Step 18 (introduction and `setup_phaser_env.sh`) into the skill**

Insert the following as the new Step 18, after Step 17's closing `---`:

`````markdown
## Step 18 — phASER: installation and per-sample haplotype counts (outbred, optional)

Only when phASER was selected in Step 7 (`{MODE}` = `outbred`, paired-end data).

phASER phases each sample's heterozygous SNPs from the reads themselves: a read pair that covers two SNPs shows which alleles lie on one molecule. It then counts the read pairs of each haplotype, once per read pair, over the SNPs of a haplotype block, and `phaser_gene_ae` turns the blocks into one haplotype-A and one haplotype-B count per gene. These counts are tested in Rmd 05 (Step 19), next to the unphased per-SNP and ACAT results of Rmd 02, never instead of them.

**Which phASER.** Use the repository https://github.com/secastel/phaser (maintained by PEJ Lab as a "Fast Beta"; Python 3; no releases or tags), pinned to commit `aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301` (phASER v1.2.0, 2026-03-21). Use this repository, never the bioconda package `phaser`, which is version 0.1.1 for Python 2.7.

**Environment.** Installed once per user by `setup_phaser_env.sh` into `{PHASER_HOME}` (Step 7) and reused by every project:
- `src/phaser`: the clone at the pinned commit;
- `env`: a conda environment, about 630 MB;
- `environment.phaser.yml`: the pinned package list;
- `conda_list.txt`;
- `install_ok.txt`: written last, after the checks pass.

The installation is a compute job: conda and `git` need the internet, which the compute nodes have, and the classic solver takes about 20 minutes. `module add miniconda3/v4` in this job is the only module besides `singularity/3.10.4` that the skill ever loads. The phASER scripts of a project load no module and do not activate the environment: they put `{PHASER_HOME}/env/bin` first on `PATH`, and `samtools`, `bcftools`, `tabix` and `python` come from it. Every script that runs the environment's Python exports `PYTHONNOUSERSITE=1`, because a `pip --user` site-packages directory of the same Python version (`~/.local/lib/python3.14`) would otherwise shadow the environment's packages (it fails with `GLIBC_2.27 not found` on this cluster).

### `setup_phaser_env.sh`

Written to `{RESULTS_DIR}/scripts/setup_phaser_env.sh` (substitute the placeholders) and submitted by `submit_chain.sh` (Step 15) only when `{PHASER_HOME}/install_ok.txt` does not start with the pinned commit. If it does, the chain skips this job. Two projects that install into the same `{PHASER_HOME}` at the same time are stopped by the lock directory `{PHASER_HOME}/.installing`. If a job was killed, remove that directory by hand.

```bash
#!/bin/bash
#SBATCH -J ase_phaser_setup
#SBATCH -N 1 -p bcc
#SBATCH -n 4 --mem=16G -t 2:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={USER_EMAIL}
#SBATCH -o {RESULTS_DIR}/logs/setup_phaser_env_%j.out
set -uo pipefail
# One-time installation of phASER (pinned commit) and its pinned conda environment under {PHASER_HOME}; reused by every project.
# never the bioconda package "phaser" (version 0.1.1, Python 2.7). Runs on a compute node (internet for git and conda).
PHASER_HOME="{PHASER_HOME}"
PHASER_REPO=https://github.com/secastel/phaser.git
PHASER_COMMIT=aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301   # no releases or tags upstream: pinned by SHA (Step 18, bumping)
die() { echo "ERROR: $*" >&2; exit 1; }
if [ "$(sed -n 1p "$PHASER_HOME/install_ok.txt" 2>/dev/null)" = "commit $PHASER_COMMIT" ]; then echo "phASER already installed in $PHASER_HOME"; exit 0; fi
mkdir -p "$PHASER_HOME" || die "cannot create $PHASER_HOME"
mkdir "$PHASER_HOME/.installing" 2>/dev/null || die "another installation is running in $PHASER_HOME (if not, remove $PHASER_HOME/.installing)"
trap 'rmdir "$PHASER_HOME/.installing" 2>/dev/null' EXIT
rm -f "$PHASER_HOME/install_ok.txt"
rm -rf "$PHASER_HOME/src/phaser" "$PHASER_HOME/env"
git clone -q "$PHASER_REPO" "$PHASER_HOME/src/phaser" || die "git clone of $PHASER_REPO failed"
git -C "$PHASER_HOME/src/phaser" checkout -q "$PHASER_COMMIT" || die "commit $PHASER_COMMIT not found"
[ "$(git -C "$PHASER_HOME/src/phaser" rev-parse HEAD)" = "$PHASER_COMMIT" ] || die "the checked-out commit differs from $PHASER_COMMIT"
cat > "$PHASER_HOME/environment.phaser.yml" <<'YML'
name: phaser-pinned
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.14.7
  - numpy=2.5.3
  - scipy=1.18.1
  - pandas=3.0.6
  - pysam=0.24.1
  - intervaltree=3.2.1
  - samtools=1.24
  - bcftools=1.24
  - htslib=1.24
  - bedtools=2.31.1
  - pip
YML
module add miniconda3/v4 || die "cannot load the miniconda3/v4 module"
source /home/software/conda/miniconda3/bin/condainit || die "cannot initialise conda"
export CONDA_PKGS_DIRS="$PHASER_HOME/pkgs"
conda env create -p "$PHASER_HOME/env" -f "$PHASER_HOME/environment.phaser.yml" || die "conda env create failed"
rm -rf "$PHASER_HOME/pkgs"     # the package cache (about 1.9 GB) is not needed after the installation
export PYTHONNOUSERSITE=1      # a pip --user site-packages directory of the same Python version would shadow the environment
export PATH="$PHASER_HOME/env/bin:$PATH"; PY="$PHASER_HOME/env/bin/python"
"$PY" -c "import numpy, scipy, pysam, pandas, intervaltree" || die "the environment's Python packages do not import"
"$PY" "$PHASER_HOME/src/phaser/phaser/phaser.py" --help > /dev/null || die "phaser.py --help failed"
"$PY" "$PHASER_HOME/src/phaser/phaser_gene_ae/phaser_gene_ae.py" --help > /dev/null || die "phaser_gene_ae.py --help failed"
for t in samtools bcftools bgzip tabix bedtools; do command -v "$t" > /dev/null || die "$t is missing from the environment"; done
conda list -p "$PHASER_HOME/env" > "$PHASER_HOME/conda_list.txt" || die "conda list failed"
{ echo "commit $PHASER_COMMIT"; "$PY" -c 'import sys; print("python", sys.version.split()[0])'; date -u; } > "$PHASER_HOME/install_ok.txt" \
  || die "cannot write $PHASER_HOME/install_ok.txt"
echo "phASER installed in $PHASER_HOME"
```

The versions in `environment.phaser.yml` are the ones conda resolved from phASER's own unpinned `environment.yml` when this installation was tested; they are pinned so that a later installation gets the tested set. A pin to another Python (for example 3.12) has not been tested.

**Bumping the pinned commit.** phASER has no releases or tags, so the skill pins a commit SHA, and a newer commit is adopted on purpose, never by accident. To adopt one:
1. Change `PHASER_COMMIT` in `setup_phaser_env.sh` and in `phaser_count.sh`, and the commit string in the phASER check of `submit_chain.sh` (Step 15). All three must agree.
2. Install into a new `{PHASER_HOME}` and re-record `ase-pipeline/tests/fixtures/phaser_help.txt` and `phaser_gene_ae_help.txt` from it (the checker compares every phASER flag in this file with them).
3. Re-run the synthetic acceptance (README, Validation status), and only then use it.

Change the version list of `environment.phaser.yml` the same way, after a test installation. A changed commit makes `setup_phaser_env.sh` reinstall into `{PHASER_HOME}`: the old clone and environment are removed, so a project that still needs the old commit must use another `{PHASER_HOME}`.
`````

- [ ] **Step 4: Install on the cluster from the skill text (the real script)**

```bash
SCR=/net/bmc-lab3/data/bcc/yannvrb/ase_stage3_scratch; T=$SCR/t3; mkdir -p $T/logs $T/scripts
bash ase-pipeline/tests/synthetic/cut_block.sh ase-pipeline/ase-pipeline.md '### `setup_phaser_env.sh`' \
  | sed -e "s|{PHASER_HOME}|$SCR/phaser_home|g" -e "s|{USER_EMAIL}|yannvrb@mit.edu|g" -e "s|{RESULTS_DIR}|$T|g" > $T/scripts/setup_phaser_env.sh
grep -nE '\{[A-Z][A-Z_0-9]*\}' $T/scripts/setup_phaser_env.sh && echo "LEFTOVER PLACEHOLDER" || true
J=$(sbatch -p bcc --parsable $T/scripts/setup_phaser_env.sh); echo $J
```

Expected:
- The job ends with exit 0, and the log ends with `phASER installed in .../phaser_home`. Read the state with `scontrol show job $J`.
- `head -1 $SCR/phaser_home/install_ok.txt` prints `commit aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301`.
- `grep -E '^(python|numpy|scipy|pandas|pysam|intervaltree|samtools|bcftools|htslib|bedtools) ' $SCR/phaser_home/conda_list.txt` shows exactly the pinned versions.
- Record the wall time from the log timestamps.

If the solve fails with the pinned versions, report BLOCKED with the log. The ruling option, proven in the gate, is the upstream `environment.yml` verbatim plus the recorded `conda list`.

Then submit the same script a second time. Expected: exit 0 within a minute, with `phASER already installed`.

- [ ] **Step 5: Record the fixtures and verify the run-time setup without activation**

The run-time setup is `PATH` plus `PYTHONNOUSERSITE`, with no module and no activation (planner decision 5). The phASER smoke run reads the Stage 2 acceptance BAM and VCF, which are read-only inputs; it writes only under `$T`.

```bash
cat > $T/scripts/fixtures.sh <<EOF
#!/bin/bash
#SBATCH -p bcc
#SBATCH -N 1 -n 1 --mem=4G -t 0:20:00
set -uo pipefail
H=$SCR/phaser_home; T=$T
export PATH="\$H/env/bin:\$PATH"; PY="\$H/env/bin/python"; mkdir -p \$T/fixtures \$T/run
echo "== import WITHOUT PYTHONNOUSERSITE (evidence only, not a gate) =="
"\$PY" -c "import numpy; print(numpy.__file__)" 2>&1 | tail -2
export PYTHONNOUSERSITE=1
"\$PY" -c "import numpy, scipy, pysam, pandas, intervaltree; print('imports ok', numpy.__file__)" || exit 1
"\$PY" \$H/src/phaser/phaser/phaser.py --help > \$T/fixtures/phaser_help.txt 2>&1 || exit 1
"\$PY" \$H/src/phaser/phaser_gene_ae/phaser_gene_ae.py --help > \$T/fixtures/phaser_gene_ae_help.txt 2>&1 || exit 1
{ echo "repository https://github.com/secastel/phaser"; echo "commit \$(git -C \$H/src/phaser rev-parse HEAD)"; sed -n 2p \$H/install_ok.txt
  grep -E '^(python|numpy|scipy|pandas|pysam|intervaltree|samtools|bcftools|htslib|bedtools) ' \$H/conda_list.txt; } > \$T/fixtures/phaser_env.txt
R=/net/bmc-lab3/data/bcc/yannvrb/ase_accept_stage2/obd/results/2026-09-30_obd
"\$PY" \$H/src/phaser/phaser/phaser.py --vcf \$R/genotypes/ind1.het.vcf.gz --bam \$R/bam/ob_s1.bam --sample ind1 --paired_end 1 --mapq 255 --baseq 10 \
  --pass_only 0 --unique_ids 1 --gw_phase_vcf 0 --threads 1 --temp_dir \$T/run --o \$T/run/ob_s1 > \$T/run/ob_s1.log 2>&1 || { tail \$T/run/ob_s1.log; exit 1; }
grep -m1 PHASED \$T/run/ob_s1.log; echo "blocks: \$(awk 'NR>1' \$T/run/ob_s1.haplotypic_counts.txt | wc -l)"
EOF
sbatch -p bcc -o $T/logs/fixtures_%j.out $T/scripts/fixtures.sh
```

Expected:
- `imports ok` with a numpy path under `phaser_home/env`.
- Two help files, each starting with `usage:`; `phaser_help.txt` lists `--paired_end`, `--pass_only`, `--unique_ids`, `--gw_phase_vcf`, `--temp_dir` and `--o`.
- The smoke run prints `PHASED 99 of 153` (or within a few of it: the gate's figure on the same BAM) and `blocks: 94`.
- The "without PYTHONNOUSERSITE" line is evidence only. Record whatever it printed.

If the smoke run fails without activation but the help works, report it. The ruling option is `conda activate` (proven); planner decision 5 is then reversed and recorded.

Copy `$T/fixtures/phaser_help.txt`, `phaser_gene_ae_help.txt` and `phaser_env.txt` into `ase-pipeline/tests/fixtures/`. Append to `ase-pipeline/tests/fixtures/verification.md` a section `## Stage 3: phASER (2026-10-01)` that lists:
- the gate report path;
- the installation job id, wall time and pinned versions (from `phaser_env.txt`);
- the idempotence re-run;
- the import evidence with and without `PYTHONNOUSERSITE`;
- the smoke-run result;
- a statement that the help fixtures come from commit aa1f8ec.

- [ ] **Step 6: Run the checker (expected PASS), then prove the phASER flag rule**

Run: `bash ase-pipeline/tests/check_skill.sh ase-pipeline/ase-pipeline.md`
Expected: `PASS`.

Proofs, each on a scratch copy `$SCR/mut.md` of the working skill:
- (a) `sed 's/phaser\.py" --help/phaser.py" --helpp/'` → `FAIL: phaser/phaser.py flag not in fixtures/phaser_help.txt: --helpp`.
- (b) `sed 's/phaser_gene_ae\.py" --help/phaser_gene_ae.py" --feature/'` → FAIL naming `--feature`.
- (c) delete the line `module add miniconda3/v4 || die "cannot load the miniconda3/v4 module"` → FAIL "exactly one line may load miniconda3/v4".

Then run the checker on `git show HEAD:...` as `$SCR/pre.md`. Expected: every new `need` line fails.

- [ ] **Step 7: Commit**

```bash
git status && git diff --stat
git add ase-pipeline/ase-pipeline.md ase-pipeline/tests/check_skill.sh ase-pipeline/tests/synthetic/cut_block.sh ase-pipeline/tests/fixtures/phaser_help.txt ase-pipeline/tests/fixtures/phaser_gene_ae_help.txt ase-pipeline/tests/fixtures/phaser_env.txt ase-pipeline/tests/fixtures/verification.md
git commit -m "ase-pipeline: phASER installation job (pinned commit and versions), help fixtures and phASER flag rule

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Rmd 05 — phASER gene-level haplotype imbalance (Step 19) and its evaluator

**Files:**
- Modify `ase-pipeline/ase-pipeline.md`: new `## Step 19 — Rmd 05: phASER gene-level haplotype imbalance (outbred, optional)`, after Step 18 and its closing `---`, before `## Notes for the assistant`.
- Create: `ase-pipeline/tests/synthetic/evaluate_stage3.R` (the `rmd` mode; Task 5 adds `phaser`).
- Modify: `ase-pipeline/tests/check_skill.sh`.

**Interfaces:**
- Consumes:
  - from Task 1: `hap_gene_test`;
  - from Rmd 01: `ase_checkpoint.rds` (`samples` with `sample, condition, individual`; `constants`; `bias`);
  - from Rmd 02: `ase_imbalance_checkpoint.rds` (`gene` with `sample, gene_id, n_snps, acat_p, padj, sig`; `constants`);
  - from Task 2: `direct_gene_ae_phased`, `direct_gene_ae_unphased`, `direct_flips.tsv` and the truth tables;
  - the existing `render_from_skill.sh <skill> <step> <values.tsv> <out.Rmd>` (it supports `19`).
- Produces:
  - `{CWD}/{TODAY}_{WD_NAME}_05_phaser.Rmd`, which renders to `{RESULTS_DIR}/{TODAY}_{WD_NAME}_05_phaser.html`;
  - `{RESULTS_DIR}/{TODAY}_{WD_NAME}_ASE_phaser.xlsx` (sheets `Gene`, `Comparison`, `Samples`, `Summary`) and `{TODAY}_{WD_NAME}_ASE_phaser.pdf`;
  - `{RESULTS_DIR}/ase_phaser_checkpoint.rds` = `list(gene, comparison, samples, summary, constants)`:
    - `gene` columns: `sample, gene_id, contig, aCount, bCount, totalCount, n_variants, gw_phased, variants, condition, individual, tested, p, rho_used, hap_A_frac, major_frac, oriented, padj, sig, direction`;
    - `comparison` columns: `sample, gene_id, totalCount, n_variants, gw_phased, major_frac, p, padj, sig, n_snps_unphased, acat_p, padj_unphased, sig_unphased, category`;
    - `constants` = Rmd 01's constants plus `RHO_MIN` and `PHASED_GT`;
  - `{RESULTS_DIR}/summary_numbers_phaser.tsv` with columns exactly `sample, genotypes, genes_with_counts, genes_tested, sig_genes, gw_phased_genes, sig_unphased, sig_both, sig_phaser_only, sig_phaser_only_untested_unphased, sig_unphased_only, rho_used, rho_own, rho_cohort`;
  - the `{PHASED_GT}` placeholder (0 or 1);
  - `Rscript evaluate_stage3.R rmd <RESULTS_DIR> <truth dir outbred_phase> <phased|unphased>`, which prints `PASS`/`FAIL`/`REPORT` lines and exits 0 or 1.

**Rmd acceptance gates** (on direct tables here; the same criteria apply in Tasks 5 and 8):

| Cells | Gate |
|---|---|
| `hap_strong` (24 cells) | all significant |
| `two_block` (12 cells) | all significant |
| `hap_lowdepth` (18 cells), phased genotypes | ≥ 15 significant (unphased genotypes: reported) |
| every significant `hap_lowdepth` cell | category "phASER only (no SNP tested unphased)" |
| null and null_linked cells | at most 6 significant |
| phased genotypes | every significant planted cell genome-wide phased, with the sign of `hap_A_frac − 0.5` equal to that of `h1 − 0.5` |
| unphased genotypes | no significant row carries a direction |
| every run | no gene row with `totalCount` 0; every row below `MIN_DEPTH` has p `NA`; one `summary_numbers_phaser.tsv` row per sample, with `sig_genes` equal to the gene table |

- [ ] **Step 1: Failing checker lines**

```bash
# --- Stage 3 Task 4 (Rmd 05, Step 19); each string is absent from the Task 3 skill
need "## Step 19 — Rmd 05: phASER gene-level haplotype imbalance (outbred, optional)"
need "_05_phaser.Rmd"
need "_ASE_phaser.xlsx"
need "ase_phaser_checkpoint.rds"
need "summary_numbers_phaser.tsv"
need "run_05_phaser.sh"
need 'PHASED_GT   <- {PHASED_GT}'
need 'ht <- hap_gene_test(ifelse(tst, genes$aCount, NA), ifelse(tst, genes$bCount, NA), genes$sample, RHO_MIN)'
need 'genes <- ga %>% dplyr::filter(totalCount > 0) %>%'
need '"no direction (haplotype labels arbitrary)"'
need "haplotype A is the haplotype of the first (left) allele of the phased genotype"
need "Rmd 05 is reported next to them, never instead of them"
need 'dplyr::full_join(unph, by = c("sample", "gene_id"))'
[ -s "$HERE/synthetic/evaluate_stage3.R" ] || { echo "FAIL: missing tests/synthetic/evaluate_stage3.R"; fail=1; }
# --- end Stage 3 Task 4
```

Run the checker. Expected: FAIL (14 lines).

- [ ] **Step 2: Write Step 19 into the skill**

Insert after Step 18's closing `---`:

``````markdown
## Step 19 — Rmd 05: phASER gene-level haplotype imbalance (outbred, optional)

Only when phASER was selected in Step 7. Write `{CWD}/{TODAY}_{WD_NAME}_05_phaser.Rmd` with the same conventions as Step 12 and the same `{AUTHOR}` and `{PROJECT_TITLE}`.

**Inputs:**
- `ase_checkpoint.rds` (Rmd 01): samples, constants, bias flags;
- `ase_imbalance_checkpoint.rds` (Rmd 02): the unphased gene results;
- `phaser/{sample}.gene_ae.txt` (Step 18), read for every sample of `samples.csv` by name.

It pastes the statistics functions of Step 14 into the chunk marked below.

**Outputs, written to `{RESULTS_DIR}`:**
- `{TODAY}_{WD_NAME}_ASE_phaser.xlsx` (sheets `Gene`, `Comparison`, `Samples`, `Summary`);
- `{TODAY}_{WD_NAME}_ASE_phaser.pdf`;
- `ase_phaser_checkpoint.rds`;
- `summary_numbers_phaser.tsv`.

- **What it adds, and what it does not replace.** The unphased per-SNP tests and the ACAT gene test of Rmd 02 stay the main per-sample result; Rmd 05 is reported next to them, never instead of them.
  - phASER counts each read pair once per haplotype block, so a gene's count pools all the heterozygous SNPs of the block without counting a fragment twice.
  - This gives one effect size per gene: the major haplotype fraction.
  - It can test genes whose SNPs each have fewer than `MIN_DEPTH` reads.
- **Test.** Per sample and gene, `hap_gene_test` (Step 14): an exact two-sided beta-binomial test of `aCount` out of `totalCount` against 0.5.
  - The dispersion comes from the trimmed central-region fit over the sample's genes (`rho_own`); the cohort value (`rho_cohort`) is the median of `rho_own`, and `rho_used = max(rho_cohort, rho_own, RHO_MIN)`.
  - BH within each sample.
  - `sig` requires `padj < FDR_SIG` and a major haplotype fraction `max(aCount, bCount) / totalCount` at least `0.5 + ABS_DEV_SIG`.
  - Genes with `totalCount` below `MIN_DEPTH` are not tested.
  - Genes that phASER wrote with zero counts (no covered heterozygous SNP; written with `log2_aFC` `inf` and `gw_phased` 1) are dropped.
  - If no sample has 20 genes at `MIN_DEPTH`, the Rmd stops.
- **Direction.** The test is symmetric in the two haplotypes, because the haplotype labels carry no meaning without phased genotypes.
  - **Unphased genotypes** (`{PHASED_GT}` = 0): A and B are arbitrary per gene and sample. Only the size of the imbalance is reported, with the direction "no direction (haplotype labels arbitrary)". Each gene's counts come from its most-covered haplotype block only, because `phaser_gene_ae` keeps a single block per gene when the genotypes are unphased. The other blocks are dropped, which costs coverage (on the synthetic test, one sample used 113 of 153 heterozygous SNPs).
  - **Phased genotypes** (`{PHASED_GT}` = 1, phASER run with `--gw_phase_vcf 1`): in genome-wide phased genes (`gw_phased` TRUE), haplotype A is the haplotype of the first (left) allele of the phased genotype (`0|1`: A carries the REF allele at that SNP), and the direction is "haplotype A higher" or "haplotype B higher". Genes that are not genome-wide phased keep arbitrary labels.
- **What comparisons the results support.** Across samples and individuals, only the size of the imbalance can be compared.
  - Haplotype A of one individual is unrelated to haplotype A of another, even with phased genotypes.
  - With unphased genotypes, even two samples of one individual may label a gene's haplotypes differently.
  - Differential ASE on phASER counts is therefore not offered (use Rmd 04).
- **Gene spans.** phASER counts every heterozygous SNP inside the gene span (`reference/genes_span.bed`, introns included), while Rmd 02 uses exonic SNPs only. Overlapping genes share reads.
- **Comparison table.** Every sample and gene of either analysis gets one row and one category:
  - "both";
  - "phASER only";
  - "phASER only (no SNP tested unphased)": no SNP of the gene reached `MIN_DEPTH` in Rmd 02;
  - "unphased only";
  - "unphased only (not tested by phASER)";
  - "neither".

`````rmd
---
title: "{PROJECT_TITLE} - ASE phASER haplotype counts"
author: "{AUTHOR}"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(cache = FALSE, echo = TRUE, message = FALSE, warning = FALSE, fig.width = 10, fig.height = 7)
options(scipen = 9)
library(openxlsx)   # no Bioconductor package is needed in this Rmd
library(tidyverse)
```

## Constants and checkpoints

```{r constants}
MODE        <- "{MODE}"
RESULTS_DIR <- "{RESULTS_DIR}"
DATE_TAG    <- "{TODAY}_{WD_NAME}"
MIN_DEPTH   <- {MIN_DEPTH}
FDR_SIG     <- {FDR_SIG}
ABS_DEV_SIG <- {ABS_DEV_SIG}
RHO_MIN     <- {RHO_MIN}             # floor for the dispersion used by the tests
BIAS_TOL    <- {BIAS_TOL}
PHASED_GT   <- {PHASED_GT}           # 1: phased genotype VCFs, phASER ran with --gw_phase_vcf 1 (Step 7); 0: unphased
if (MODE != "outbred") stop("Rmd 05 (phASER) needs MODE = outbred", call. = FALSE)
if (!PHASED_GT %in% c(0, 1)) stop("PHASED_GT must be 0 or 1 (Step 7)", call. = FALSE)
ck <- readRDS(file.path(RESULTS_DIR, "ase_checkpoint.rds"))
same <- c(MODE = identical(ck$constants$MODE, MODE),
          vapply(c("MIN_DEPTH", "FDR_SIG", "ABS_DEV_SIG", "BIAS_TOL"),
                 function(k) isTRUE(all.equal(ck$constants[[k]], get(k))), logical(1)))
if (!all(same)) stop("constants differ from the Rmd 01 checkpoint (", paste(names(same)[!same], collapse = ", "),
                     "); re-render Rmd 01 with the same values", call. = FALSE)
ck2 <- readRDS(file.path(RESULTS_DIR, "ase_imbalance_checkpoint.rds"))   # Rmd 02: the unphased per-SNP and ACAT gene results
if (!isTRUE(all.equal(ck2$constants, ck$constants)))
  stop("the Rmd 02 checkpoint does not come from the current Rmd 01 checkpoint; re-render Rmd 02", call. = FALSE)
bias <- ck$bias
knitr::kable(bias, digits = 4, caption = "Reference-bias diagnostic from Rmd 01 (a screen, not a test)")
if (any(bias$flagged)) cat("FLAGGED samples:", paste(bias$sample[bias$flagged], collapse = ", "), "- read their results with caution\n")
orient_note <- if (PHASED_GT == 1) {
  paste("Phased genotypes: in genome-wide phased genes, haplotype A is the haplotype of the first (left) allele of the phased genotype",
        "(0|1: A carries REF at that SNP); haplotype A of one individual is unrelated to haplotype A of another;",
        "genes that are not genome-wide phased have arbitrary labels.")
} else {
  paste("Unphased genotypes: the haplotype labels A and B are arbitrary per gene and sample, so only the size of the imbalance",
        "(major haplotype fraction) is reported, without a direction; each gene's counts come from its most-covered haplotype block.")
}
cat("NOTE:", orient_note, "\n")
```

## Statistics functions

```{r stats}
# <<< paste here the code of Step 14 between the two marker lines (the marker lines themselves are not pasted) >>>
```

## phASER gene counts

The gene tables are read by sample name from `samples.csv`, never by globbing.

```{r read}
samples <- unique(ck$samples[, c("sample", "condition", "individual")])
f <- file.path(RESULTS_DIR, "phaser", paste0(samples$sample, ".gene_ae.txt"))
miss <- samples$sample[!file.exists(f) | is.na(file.size(f)) | file.size(f) == 0]
if (length(miss) > 0)
  stop("phASER gene counts missing for sample(s) ", paste(miss, collapse = ", "), " (", file.path(RESULTS_DIR, "phaser"),
       "/<sample>.gene_ae.txt): check their phaser_count logs, re-run them, or drop the samples (Step 15)", call. = FALSE)
cls <- c(contig = "character", start = "numeric", stop = "numeric", name = "character", aCount = "numeric", bCount = "numeric",
         totalCount = "numeric", log2_aFC = "character", n_variants = "numeric", variants = "character", gw_phased = "character",
         bam = "character")
ga <- dplyr::bind_rows(lapply(seq_len(nrow(samples)), function(i) {
  d <- read.delim(f[i], colClasses = cls, quote = "", na.strings = character(0))
  if (!identical(names(d), names(cls))) stop("unexpected columns in ", f[i], call. = FALSE)
  dplyr::mutate(d, sample = samples$sample[i])
}))
n_rows <- nrow(ga)
genes <- ga %>% dplyr::filter(totalCount > 0) %>%   # genes without counts are written with 0 counts, log2_aFC inf and gw_phased 1
  dplyr::transmute(sample, gene_id = name, contig, aCount, bCount, totalCount, n_variants, gw_phased = gw_phased == "1", variants) %>%
  dplyr::left_join(samples, by = "sample")
if (nrow(genes) == 0) stop("no gene has haplotype counts in any sample (every totalCount is 0)", call. = FALSE)
if (any(genes$aCount + genes$bCount != genes$totalCount)) stop("aCount + bCount differs from totalCount in the phASER gene tables", call. = FALSE)
if (PHASED_GT == 1 && mean(genes$gw_phased) < 0.5)
  cat("WARNING: PHASED_GT = 1 but only", sum(genes$gw_phased), "of", nrow(genes), "genes with counts are genome-wide phased\n")
if (PHASED_GT == 0 && any(genes$gw_phased))
  cat("WARNING: PHASED_GT = 0 but", sum(genes$gw_phased), "genes are genome-wide phased; their labels are treated as arbitrary\n")
cat(n_rows, "sample-gene rows read;", nrow(genes), "with haplotype counts;", sum(genes$totalCount >= MIN_DEPTH), "at or above MIN_DEPTH\n")
```

## Haplotype-level gene test

Per sample and gene, `hap_gene_test` (Step 14) tests the haplotype counts against 0.5. The test is symmetric in the two haplotypes. The dispersion is the larger of the sample's own trimmed fit and the cohort median, floored at `RHO_MIN`.

```{r test}
tst <- genes$totalCount >= MIN_DEPTH
ht <- hap_gene_test(ifelse(tst, genes$aCount, NA), ifelse(tst, genes$bCount, NA), genes$sample, RHO_MIN)
if (all(is.na(ht$p)))
  stop("no sample has 20 genes with at least MIN_DEPTH haplotype reads, so no dispersion can be estimated; genes per sample at MIN_DEPTH: ",
       paste(sprintf("%s %d", ht$samples$sample, ht$samples$n_genes), collapse = ", "), call. = FALSE)
genes <- genes %>%
  dplyr::mutate(tested = tst, p = ht$p, rho_used = ht$rho_used, hap_A_frac = aCount / totalCount,
                major_frac = pmax(aCount, bCount) / totalCount, oriented = PHASED_GT == 1 & gw_phased) %>%
  dplyr::group_by(sample) %>% dplyr::mutate(padj = p.adjust(p, method = "BH")) %>% dplyr::ungroup() %>%
  dplyr::mutate(sig = !is.na(padj) & padj < FDR_SIG & major_frac - 0.5 >= ABS_DEV_SIG,
                direction = dplyr::case_when(!sig ~ "none", !oriented ~ "no direction (haplotype labels arbitrary)",
                                             hap_A_frac > 0.5 ~ "haplotype A higher", TRUE ~ "haplotype B higher")) %>%
  dplyr::arrange(sample, gene_id)
knitr::kable(ht$samples, digits = 4, caption = paste("Dispersion per sample: rho_own = trimmed central-region fit over the sample's genes",
             "(NA below 20 genes), rho_cohort = median of rho_own, rho_used = max(rho_cohort, rho_own, RHO_MIN)"))
knitr::kable(head(dplyr::filter(genes, sig), 40), digits = 4, caption = paste("Significant genes (first 40).", orient_note))
```

## Next to the unphased result of Rmd 02

```{r compare}
unph <- ck2$gene %>% dplyr::transmute(sample, gene_id, n_snps_unphased = n_snps, acat_p, padj_unphased = padj, sig_unphased = sig)
cmp <- genes %>% dplyr::select(sample, gene_id, totalCount, n_variants, gw_phased, major_frac, p, padj, sig) %>%
  dplyr::full_join(unph, by = c("sample", "gene_id")) %>%
  dplyr::mutate(sig = !is.na(sig) & sig, sig_unphased = !is.na(sig_unphased) & sig_unphased,
                category = dplyr::case_when(sig & sig_unphased ~ "both",
                                            sig & is.na(acat_p) ~ "phASER only (no SNP tested unphased)",
                                            sig ~ "phASER only",
                                            sig_unphased & is.na(p) ~ "unphased only (not tested by phASER)",
                                            sig_unphased ~ "unphased only",
                                            TRUE ~ "neither")) %>%
  dplyr::arrange(sample, gene_id)
knitr::kable(dplyr::count(cmp, category), caption = "Sample-gene pairs by category: phASER haplotype test next to the unphased ACAT test of Rmd 02")
```

## Summary

```{r summary}
z0 <- function(x) dplyr::coalesce(as.integer(x), 0L)
summary_ph <- samples %>% dplyr::select(sample) %>%
  dplyr::left_join(dplyr::select(ht$samples, sample, rho_used, rho_own, rho_cohort), by = "sample") %>%
  dplyr::left_join(genes %>% dplyr::group_by(sample) %>%
                     dplyr::summarise(genes_with_counts = dplyr::n(), genes_tested = sum(!is.na(p)), sig_genes = sum(sig),
                                      gw_phased_genes = sum(gw_phased), .groups = "drop"), by = "sample") %>%
  dplyr::left_join(cmp %>% dplyr::group_by(sample) %>%
                     dplyr::summarise(sig_unphased = sum(sig_unphased), sig_both = sum(category == "both"),
                                      sig_phaser_only = sum(startsWith(category, "phASER only")),
                                      sig_phaser_only_untested_unphased = sum(category == "phASER only (no SNP tested unphased)"),
                                      sig_unphased_only = sum(startsWith(category, "unphased only")), .groups = "drop"), by = "sample") %>%
  dplyr::mutate(dplyr::across(c(genes_with_counts, genes_tested, sig_genes, gw_phased_genes, sig_unphased, sig_both, sig_phaser_only,
                                sig_phaser_only_untested_unphased, sig_unphased_only), z0),
                genotypes = if (PHASED_GT == 1) "phased" else "unphased") %>%
  dplyr::select(sample, genotypes, genes_with_counts, genes_tested, sig_genes, gw_phased_genes, sig_unphased, sig_both, sig_phaser_only,
                sig_phaser_only_untested_unphased, sig_unphased_only, rho_used, rho_own, rho_cohort)
knitr::kable(summary_ph, digits = 4, caption = "Summary per sample (FDR_SIG, ABS_DEV_SIG and MIN_DEPTH as in the constants block)")
```

## Figure

```{r figure}
fig <- dplyr::filter(genes, !is.na(p))
fig_n <- vapply(summary_ph$sample, function(s) sum(fig$sig[fig$sample == s]), integer(1))
stopifnot(identical(unname(fig_n), summary_ph$sig_genes))
cat("Figure sig counts equal Summary counts: TRUE\n")
lab <- fig %>% dplyr::group_by(sample) %>%
  dplyr::summarise(lab = paste0(dplyr::first(sample), ": ", sum(sig), " significant of ", dplyr::n(), " genes"), .groups = "drop")
p1 <- dplyr::left_join(fig, lab, by = "sample") %>%
  ggplot(aes(totalCount, major_frac, colour = sig)) + geom_hline(yintercept = 0.5, linetype = 2) + geom_point(alpha = 0.7) +
  scale_x_log10() + scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick")) + facet_wrap(~lab) +
  labs(x = "haplotype reads per gene (log10)", y = "major haplotype fraction (folded, no direction)", colour = "sig",
       caption = paste0("sig: BH < ", FDR_SIG, " and major fraction - 0.5 >= ", ABS_DEV_SIG,
                        if (any(bias$flagged)) "; some samples flagged for reference bias" else "")) + theme_bw()
print(p1)
ggsave(file.path(RESULTS_DIR, paste0(DATE_TAG, "_ASE_phaser.pdf")), p1, width = 10, height = 7)
```

## Export

```{r export}
xlsx_file <- file.path(RESULTS_DIR, paste0(DATE_TAG, "_ASE_phaser.xlsx"))
wb <- openxlsx::createWorkbook()
for (nm in c("Gene", "Comparison", "Samples", "Summary")) openxlsx::addWorksheet(wb, nm)
openxlsx::writeData(wb, "Gene", genes); openxlsx::writeData(wb, "Comparison", cmp)
openxlsx::writeData(wb, "Samples", ht$samples); openxlsx::writeData(wb, "Summary", summary_ph)
openxlsx::saveWorkbook(wb, xlsx_file, overwrite = TRUE)
saveRDS(list(gene = genes, comparison = cmp, samples = ht$samples, summary = summary_ph,
             constants = c(ck$constants, list(RHO_MIN = RHO_MIN, PHASED_GT = PHASED_GT))),
        file.path(RESULTS_DIR, "ase_phaser_checkpoint.rds"))
write.table(summary_ph, file.path(RESULTS_DIR, "summary_numbers_phaser.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", xlsx_file, "and summary_numbers_phaser.tsv\n")
sessionInfo()
```
`````

Render it with `{RESULTS_DIR}/scripts/run_05_phaser.sh`. This is the same script as `run_01_import_qc.sh` (Step 12), with these differences:
- job name `ase_05_phaser`, log `run_05_phaser_%j.out`, the Rmd 05 file name;
- the first case of the `--bind` rule (`--bind {CWD}`): the Rmd reads no GTF;
- `-n 1 --mem=16G -t 4:00:00`. It runs single-threaded. The unit-test runtime check (Step 14) gives the time per sample of genome scale; the request covers about 100 samples.

If the job stops with "phASER gene counts missing for sample(s)", the named samples have no phASER output: look at their `phaser_count` logs, re-run them, or drop them (Step 15). If it stops with "the Rmd 02 checkpoint does not come from the current Rmd 01 checkpoint", re-render Rmd 02 first. Submission order: Step 15.
``````

- [ ] **Step 3: Write the evaluator (`rmd` mode)**

Create `ase-pipeline/tests/synthetic/evaluate_stage3.R`:

```r
# Usage: Rscript evaluate_stage3.R rmd    <RESULTS_DIR> <truth dir (outbred_phase)> <phased|unphased>
#        Rscript evaluate_stage3.R phaser <RESULTS_DIR> <truth dir (outbred_phase)> <phased|unphased>
# rmd: Rmd 05 results (ase_phaser_checkpoint.rds, summary_numbers_phaser.tsv) against the planted truth.
# phaser: the raw phASER outputs (RESULTS_DIR/phaser/<sample>.*) against the truth (added by the phaser_count task).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4 || !args[1] %in% c("rmd", "phaser") || !args[4] %in% c("phased", "unphased"))
  stop("usage: Rscript evaluate_stage3.R <rmd|phaser> <RESULTS_DIR> <truth dir> <phased|unphased>", call. = FALSE)
what <- args[1]; RES <- args[2]; TRUTH <- args[3]; GT <- args[4]
fail <- 0
crit <- function(cond, msg) { cat(if (isTRUE(cond)) "PASS" else "FAIL", msg, "\n"); if (!isTRUE(cond)) fail <<- 1 }
rd <- function(f, ...) read.delim(file.path(TRUTH, f), stringsAsFactors = FALSE, ...)
tg <- rd("truth_genes.tsv"); tig <- rd("truth_individual_genes.tsv"); tf <- rd("truth_fragments.tsv")
tp <- rd("truth_phase.tsv", colClasses = c(gt = "character", alt_on = "character")); ts <- rd("truth_snps.tsv")
sm <- read.csv(file.path(TRUTH, "samples.csv"), stringsAsFactors = FALSE)
cell <- merge(merge(tf, sm[, c("sample", "individual")], by = "sample"), tig, by = c("individual", "gene_id"))
cell <- merge(cell, tg[, c("gene_id", "class")], by = "gene_id")
if (what == "rmd") {
  ck <- readRDS(file.path(RES, "ase_phaser_checkpoint.rds"))
  crit(identical(as.numeric(ck$constants$PHASED_GT), if (GT == "phased") 1 else 0), sprintf("checkpoint PHASED_GT matches '%s'", GT))
  g <- ck$gene; MIN_D <- ck$constants$MIN_DEPTH
  crit(nrow(g) > 0 && all(g$totalCount > 0) && all(is.na(g$p[g$totalCount < MIN_D])),
       "no gene row with totalCount 0, and every row below MIN_DEPTH has p NA")
  m <- merge(cell, g[, c("sample", "gene_id", "aCount", "bCount", "totalCount", "gw_phased", "p", "padj", "sig", "hap_A_frac",
                         "major_frac", "direction")], by = c("sample", "gene_id"), all.x = TRUE)
  m$sig[is.na(m$sig)] <- FALSE
  cnt <- function(cls) c(sum(m$sig[m$class == cls]), sum(m$class == cls))
  for (cls in c("hap_strong", "two_block")) {
    x <- cnt(cls); crit(x[2] > 0 && x[1] == x[2], sprintf("%s cells significant: %d of %d (all)", cls, x[1], x[2]))
  }
  x <- cnt("hap_lowdepth")
  if (GT == "phased") {
    crit(x[2] == 18 && x[1] >= 15, sprintf("hap_lowdepth cells significant (phased genotypes): %d of %d (>= 15)", x[1], x[2]))
  } else {
    cat(sprintf("REPORT hap_lowdepth cells significant (unphased genotypes, single best block): %d of %d\n", x[1], x[2]))
  }
  cmp <- ck$comparison
  lo <- merge(m[m$class == "hap_lowdepth", c("sample", "gene_id", "sig")], cmp[, c("sample", "gene_id", "category")],
              by = c("sample", "gene_id"), all.x = TRUE)
  n_only <- sum(lo$sig & lo$category %in% "phASER only (no SNP tested unphased)")
  crit(n_only == sum(lo$sig), sprintf("significant hap_lowdepth cells labelled 'phASER only (no SNP tested unphased)': %d of %d", n_only, sum(lo$sig)))
  mod <- tg$gene_id[tg$class == "hap_moderate"]
  x <- cnt("hap_moderate")
  cat(sprintf("REPORT hap_moderate cells significant: phASER %d of %d; unphased ACAT (Rmd 02) %d\n", x[1], x[2],
              sum(cmp$sig_unphased[cmp$gene_id %in% mod])))
  nul <- m$class %in% c("null", "null_linked")
  crit(sum(nul) > 0 && sum(m$sig[nul]) <= 6, sprintf("null cells significant: %d of %d tested (limit 6)", sum(m$sig[nul]), sum(nul & !is.na(m$p))))
  pl <- m$sig & !nul
  if (GT == "phased") {
    agree <- sign(m$hap_A_frac[pl] - 0.5) == sign(m$h1[pl] - 0.5)
    crit(sum(pl) > 0 && all(m$gw_phased[pl]) && all(agree),
         sprintf("phased: significant planted cells genome-wide phased, haplotype A = haplotype 1 (left GT allele): %d of %d", sum(agree), sum(pl)))
  } else {
    crit(sum(g$sig) > 0 && all(g$direction[g$sig] == "no direction (haplotype labels arbitrary)"), "unphased: no significant row carries a direction")
  }
  s <- read.delim(file.path(RES, "summary_numbers_phaser.tsv"), stringsAsFactors = FALSE)
  crit(identical(names(s), c("sample", "genotypes", "genes_with_counts", "genes_tested", "sig_genes", "gw_phased_genes", "sig_unphased", "sig_both",
                             "sig_phaser_only", "sig_phaser_only_untested_unphased", "sig_unphased_only", "rho_used", "rho_own", "rho_cohort")) &&
       nrow(s) == nrow(sm) && all(s$sig_genes == vapply(s$sample, function(x) sum(g$sig[g$sample == x]), integer(1))),
       "summary_numbers_phaser.tsv: exact columns, one row per sample, sig_genes equal to the gene table")
}
quit(status = fail)
```

- [ ] **Step 4: Run the checker. Expected: PASS**

- [ ] **Step 5: Render on the direct tables (phased, unphased flipped, unphased unflipped) and evaluate**

```bash
SCR=/net/bmc-lab3/data/bcc/yannvrb/ase_stage3_scratch; SYN=/net/bmc-lab3/data/bcc/yannvrb/ase_synthetic_stage3; O=$SYN/outbred_phase
T=ase-pipeline/tests/synthetic
mkp() {   # $1 project, $2 PHASED_GT, $3 gene-table dir
  local P=$SCR/$1 RD=$SCR/$1/results/2026-10-01_$1; mkdir -p $RD/ase_counts $RD/phaser $RD/logs
  cp $O/direct_counts/* $RD/ase_counts/; cp $O/samples.csv $P/$1_samples.csv; cp -r $SYN/genome $P/genome
  for s in obp_s1 obp_s2 obp_s3 obp_s4 obp_s5 obp_s6; do cp $3/$s.gene_ae.txt $RD/phaser/; done
  printf '%s\t%s\n' MODE outbred STRAIN_A REF STRAIN_B ALT SAMPLES_CSV $P/$1_samples.csv RESULTS_DIR $RD GTF_PATH $P/genome/genome.gtf \
    MIN_DEPTH 10 FDR_SIG 0.05 ABS_DEV_SIG 0.1 BIAS_TOL 0.03 RHO_MIN 0.01 THIN_BP 500 TODAY 2026-10-01 WD_NAME $1 \
    PROJECT_TITLE "Stage 3 test" AUTHOR tester PHASED_GT $2 > $P/values.tsv
  for n in 12:01_import_qc 13:02_imbalance 19:05_phaser; do
    bash $T/render_from_skill.sh ase-pipeline/ase-pipeline.md ${n%%:*} $P/values.tsv $P/2026-10-01_$1_${n#*:}.Rmd || return 1
  done
}
mkp t4_ph 1 $O/direct_gene_ae_phased && mkp t4_unph 0 $O/direct_gene_ae_unphased && mkp t4_unflip 0 $O/direct_gene_ae_unphased
for s in obp_s1 obp_s2 obp_s3 obp_s4 obp_s5 obp_s6; do   # undo the emulated flips: swap a/b where direct_flips says TRUE
  f=$SCR/t4_unflip/results/2026-10-01_t4_unflip/phaser/$s.gene_ae.txt
  awk -F'\t' 'BEGIN{OFS=FS} NR==FNR{ if ($5=="TRUE") fl[$1 SUBSEP $2]=1; next } FNR==1{print; next}
    { if (($12 SUBSEP $4) in fl) { t=$5; $5=$6; $6=t } print }' $O/direct_flips.tsv $f > $f.tmp && mv $f.tmp $f
done
cat > $SCR/t4_render.sh <<EOF
#!/bin/bash
#SBATCH -p bcc
#SBATCH -n 1 --mem=16G -t 1:00:00
set -uo pipefail
module add singularity/3.10.4 || exit 1
R() { singularity exec --bind /net/bmc-lab3 /net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif "\$@"; }
for p in t4_ph t4_unph t4_unflip; do
  for r in 01_import_qc 02_imbalance 05_phaser; do
    R Rscript -e "rmarkdown::render('$SCR/\$p/2026-10-01_\${p}_\$r.Rmd', output_dir = '$SCR/\$p/results/2026-10-01_\$p', knit_root_dir = '$SCR/\$p')" || exit 1
  done
done
EV=/net/bmc-lab3/data/bcc/yannvrb/Github_repos/bioinformatics-claude-skills/$T/evaluate_stage3.R
R Rscript \$EV rmd $SCR/t4_ph/results/2026-10-01_t4_ph $O phased; e1=\$?
R Rscript \$EV rmd $SCR/t4_unph/results/2026-10-01_t4_unph $O unphased; e2=\$?
R Rscript -e 'a <- readRDS("$SCR/t4_unph/results/2026-10-01_t4_unph/ase_phaser_checkpoint.rds")\$gene; b <- readRDS("$SCR/t4_unflip/results/2026-10-01_t4_unflip/ase_phaser_checkpoint.rds")\$gene
  ok <- identical(a\$gene_id, b\$gene_id) && identical(a\$sample, b\$sample) && isTRUE(all.equal(a\$p, b\$p)) && isTRUE(all.equal(a\$padj, b\$padj)) &&
        identical(a\$sig, b\$sig) && isTRUE(all.equal(a\$major_frac, b\$major_frac)) && isTRUE(all.equal(a\$rho_used, b\$rho_used))
  cat(if (ok) "ORIENTATION INVARIANCE PASS\n" else "ORIENTATION INVARIANCE FAIL\n"); quit(status = !ok)'; e3=\$?
echo "EXIT \$e1 \$e2 \$e3"; [ \$e1 -eq 0 ] && [ \$e2 -eq 0 ] && [ \$e3 -eq 0 ]
EOF
sbatch -p bcc -o $SCR/logs/t4_render_%j.out $SCR/t4_render.sh
```

Expected:
- All 9 renders exit 0.
- Every evaluator line is `PASS`; `REPORT` lines are informational.
- `ORIENTATION INVARIANCE PASS`.
- The log ends with `EXIT 0 0 0`.

If an evaluator line fails, check the Rmd first (the totalCount filter, the `sig` rule, the comparison join), never the truth. A statistical gate that fails with a correct Rmd is reported BLOCKED with the numbers.

- [ ] **Step 6: Negative tests (Review Focus 2, 4 and 5)**

1. **Missing sample:**
   - Copy `$SCR/t4_unph` to `$SCR/t4_miss` (adapt `RESULTS_DIR` in `values.tsv`, regenerate Rmd 05 from the skill) and delete `phaser/obp_s3.gene_ae.txt`.
   - Render Rmd 05.
   - Expected: the job stops, and its log contains `phASER gene counts missing for sample(s) obp_s3`.
2. **Uncovered rows:** count the rows with `totalCount` 0 in the six `t4_unph` input tables (`awk -F'\t' 'FNR>1 && $7==0' .../phaser/*.gene_ae.txt | wc -l`) and confirm that the Rmd 05 log says the read row count minus that number "with haplotype counts".
3. **Direction labels:** in `t4_unph`, confirm that every significant row of the `Gene` sheet of `_ASE_phaser.xlsx` has the direction "no direction (haplotype labels arbitrary)". The evaluator line already checks this; record the count.

- [ ] **Step 7: Commit**

```bash
git status && git diff --stat
git add ase-pipeline/ase-pipeline.md ase-pipeline/tests/check_skill.sh ase-pipeline/tests/synthetic/evaluate_stage3.R
git commit -m "ase-pipeline: Rmd 05 phASER gene-level haplotype test next to the unphased result (Step 19) and its truth evaluator

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 5: `phaser_count.sh` and the gene spans (Steps 10 and 18), the outbred assembler, and a smoke run on aligned synthetic data

**Files:**
- Modify `ase-pipeline/ase-pipeline.md`:
  - Step 10: a new subsection `### Gene spans for phASER (end of \`prep_genotypes.sh\`, before block I)` directly after the `prep_genotypes.sh` paragraph that starts "The array job reads `{RESULTS_DIR}/genotypes/{individual}.het.vcf`", and that paragraph's text;
  - Step 18: the `phaser_count.sh` section after the installation section.
- Create: `ase-pipeline/tests/synthetic/assemble_outbred_from_skill.sh`.
- Modify: `ase-pipeline/tests/synthetic/evaluate_stage3.R` (the `phaser` mode), `ase-pipeline/tests/check_skill.sh`.

**Interfaces:**
- Consumes:
  - Task 3's installation contract and `cut_block.sh`;
  - Task 4's Step 19 and `render_from_skill.sh`;
  - Task 2's data.
- Produces:
  - `{RESULTS_DIR}/reference/genes_span.bed`;
  - the skill block `phaser_count.sh` (anchor `` ### `phaser_count.sh` ``), a complete script with placeholders `{ARRAY_N}`, `{PHASER_RESOURCES}`, `{USER_EMAIL}`, `{RESULTS_DIR}`, `{PHASER_HOME}`, `{PHASED_GT}`, `{SAMPLES_CSV}` and `{MIN_BASEQ}`;
  - per sample: `{RESULTS_DIR}/phaser/{sample}.haplotypic_counts.txt`, `.gene_ae.txt`, the other phASER files, and the log `{RESULTS_DIR}/logs/phaser_{sample}.log`;
  - `bash assemble_outbred_from_skill.sh <skill.md> <values.tsv> <genotype_map.tsv>`. It writes these into `{RESULTS_DIR}/scripts/` and `{CWD}`:
    - `prep_genotypes.sh`, `align_wasp_count.sh`, `run_01_import_qc.sh`, `run_02_imbalance.sh`, `run_05_phaser.sh`, `setup_phaser_env.sh`, `phaser_count.sh`, `submit_chain.sh`, `wait_chain.sh`;
    - the three Rmds 01, 02 and 05;
  - `Rscript evaluate_stage3.R phaser <RESULTS_DIR> <truth dir> <phased|unphased>`.

**Gates of the `phaser` mode** (here and in Task 8), per sample:
- at least 20 multi-SNP blocks, of which ≥ 99 % are phased as in the truth, and at least one block joining SNPs in trans;
- no gene count above the simulated fragments covering the gene's het SNPs;
- `two_block` genes: unphased genotypes give `n_variants` ≤ 3 (single block); phased genotypes give `n_variants` = 6 and `gw_phased` 1;
- phased genotypes: every covered gene genome-wide phased, and `aCount` = haplotype 1 (left GT allele) in every planted gene whose fraction deviates by more than 0.1.

- [ ] **Step 1: Failing checker lines**

```bash
# --- Stage 3 Task 5 (phaser_count.sh, gene spans, assembler); each string is absent from the Task 4 skill
need "### Gene spans for phASER (end of \`prep_genotypes.sh\`, before block I)"
need 'SPAN_BED="{RESULTS_DIR}/reference/genes_span.bed"'
need "### \`phaser_count.sh\`"
need '--paired_end 1 --mapq 255 --baseq {MIN_BASEQ}'
need '--pass_only 0 --unique_ids 1 --gw_phase_vcf "$PHASED_GT"'
need 'rm -f "$OUT".*'
need "phaser_gene_ae splits variant IDs on '_'"
need "PHASED_GT=1 but only"
need 'N_COV=$(awk -F'"'"'\t'"'"' '"'"'NR > 1 && $7 > 0'"'"' "$OUT.gene_ae.txt"'
[ "$(grep -c 'PHASER_COMMIT=aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301' "$SKILL")" = 2 ] || { echo "FAIL: the pinned commit must be set in setup_phaser_env.sh and phaser_count.sh"; fail=1; }
[ -s "$HERE/synthetic/assemble_outbred_from_skill.sh" ] || { echo "FAIL: missing tests/synthetic/assemble_outbred_from_skill.sh"; fail=1; }
# --- end Stage 3 Task 5
```

Run the checker. Expected: FAIL.

- [ ] **Step 2: Gene spans in `prep_genotypes.sh` (Step 10)**

1. **Prose.** Replace the sentence "The array job reads `{RESULTS_DIR}/genotypes/{individual}.het.vcf` (STAR) and `{individual}.het.vcf.gz` (ASEReadCounter; a plain VCF fails there). Block R has already run before the loop; then block I with ..." with: "The array job reads `{RESULTS_DIR}/genotypes/{individual}.het.vcf` (STAR) and `{individual}.het.vcf.gz` (ASEReadCounter and phASER; a plain VCF fails there). Block R has already run before the loop; then the gene-spans block below; then block I with ...". Keep the rest of the paragraph.
2. **New subsection.** Insert this after that paragraph:

`````markdown
### Gene spans for phASER (end of `prep_genotypes.sh`, before block I)

Written in every outbred project (a few seconds; phASER, if selected, reads it in Step 18). One row per gene: the span of the gene's exons on its contig, 0-based start, 4 columns (contig, start, stop, `gene_id`), as `phaser_gene_ae` needs. Genes whose exons lie on more than one contig or strand are left out and counted. Built from the exons, not from `gene` lines, so it works for GTFs without `gene` lines and uses the same `gene_id`s as Rmd 02.

```bash
SPAN_BED="{RESULTS_DIR}/reference/genes_span.bed"; mkdir -p "{RESULTS_DIR}/reference" || exit 1
awk -F'\t' '$3 == "exon" && match($9, /gene_id "[^"]+"/) {
    g = substr($9, RSTART + 9, RLENGTH - 10); k = $1 SUBSEP $7
    if (!(g in key)) { key[g] = k; lo[g] = $4; hi[g] = $5 } else if (key[g] != k) bad[g] = 1
    else { if ($4 < lo[g]) lo[g] = $4; if ($5 > hi[g]) hi[g] = $5 } }
  END { for (g in key) if (!(g in bad)) { split(key[g], c, SUBSEP); print c[1] "\t" lo[g] - 1 "\t" hi[g] "\t" g }
        n = 0; for (g in bad) n++; if (n > 0) print n " genes on more than one contig or strand left out" > "/dev/stderr" }' "{GTF_PATH}" \
  | sort -k1,1 -k2,2n > "$SPAN_BED" || { echo "ERROR: cannot write $SPAN_BED" >&2; exit 1; }
N_SPAN=$(wc -l < "$SPAN_BED")
N_SPAN_ON=$(awk 'NR == FNR {c[$1] = 1; next} ($1 in c)' "$FASTA.fai" "$SPAN_BED" | wc -l)
[ "$N_SPAN" -gt 0 ] && [ "$N_SPAN_ON" -gt 0 ] \
  || { echo "ERROR: no gene span of {GTF_PATH} lies on a contig of $FASTA (GTF contigs: $(cut -f1 "$SPAN_BED" | uniq | head -3 | paste -sd' '); FASTA contigs: $(cut -f1 "$FASTA.fai" | head -3 | paste -sd' '))" >&2; exit 1; }
echo "Gene spans for phASER: $N_SPAN genes, $N_SPAN_ON on contigs of the FASTA"
```
`````

- [ ] **Step 3: `phaser_count.sh` in Step 18**

Insert after the "**Bumping the pinned commit.**" paragraph of Step 18:

`````markdown
### `phaser_count.sh`

A SLURM array job, one task per data row of `{SAMPLES_CSV}`, submitted by `submit_chain.sh` (Step 15) after the per-sample array job (Step 11) and, when it runs, the installation job. Write it to `{RESULTS_DIR}/scripts/phaser_count.sh` and substitute the placeholders: `{PHASER_RESOURCES}` from Step 9, and `{PHASED_GT}` and `{PHASER_HOME}` from Step 7.

**Inputs.**
- The BAM is `bam/{sample}.bam`, which is WASP-filtered and duplicate-marked; phASER does no WASP of its own and removes the duplicates itself.
- The VCF is `genotypes/{individual}.het.vcf.gz`, the heterozygous VCF that STAR and ASEReadCounter use. phASER needs the sample name inside that VCF, so the script reads it from the file.

**Flags.** The flags are the tested set:
- `--paired_end 1` (the reason phASER is offered only for paired-end data);
- `--mapq 255` (STAR's unique alignments, the same reads as ASEReadCounter's `--min-mapping-quality` 10 keeps);
- `--baseq {MIN_BASEQ}`;
- `--pass_only 0` (the prepared VCF has FILTER `.`; with the default 1 phASER keeps no site);
- `--unique_ids 1` (its ID column is `.`);
- `--gw_phase_vcf` 1 for phased genotypes, 0 otherwise.

**Guards.** The script stops, and leaves no output of the sample, in these cases:
- an installation of another commit;
- a contig name containing `_` (`phaser_gene_ae` splits variant IDs on `_`);
- het sites on contigs missing from the BAM header, or no gene span on a BAM contig (mismatched names give silent zero counts in phASER);
- `{PHASED_GT}` = 1 while fewer than 90 % of the het genotypes are phased;
- an empty haplotype table, or a gene table without any gene with counts.

`````

followed by this code block (```` ```bash ```` ... ```` ``` ````):

```bash
#!/bin/bash
#SBATCH -J ase_phaser
#SBATCH -N 1 -p bcc
#SBATCH --array=1-{ARRAY_N}
#SBATCH {PHASER_RESOURCES}
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={USER_EMAIL}
#SBATCH -o {RESULTS_DIR}/logs/phaser_count_%A_%a.out
set -uo pipefail
# phASER for one sample (array task = data row of samples.csv): read-backed phasing on the WASP-filtered, duplicate-marked BAM
# with the heterozygous VCF that STAR used, then gene-level haplotype counts. No module and no container: the pinned environment
# of setup_phaser_env.sh is used through its bin directory (python, samtools, bcftools and tabix come from it).
PHASER_HOME="{PHASER_HOME}"; PHASER_COMMIT=aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301; PHASED_GT={PHASED_GT}
R="{RESULTS_DIR}"; BED="$R/reference/genes_span.bed"
export PYTHONNOUSERSITE=1      # a pip --user site-packages directory of the same Python version would shadow the environment
export PATH="$PHASER_HOME/env/bin:$PATH"; PY="$PHASER_HOME/env/bin/python"; SRC="$PHASER_HOME/src/phaser"
SAMPLE=""
die() { echo "ERROR: sample ${SAMPLE:-?}: $*" >&2; exit 1; }
[ "$(sed -n 1p "$PHASER_HOME/install_ok.txt" 2>/dev/null)" = "commit $PHASER_COMMIT" ] && [ -x "$PY" ] \
  || die "no phASER installation of commit $PHASER_COMMIT in $PHASER_HOME (run setup_phaser_env.sh, Step 18)"
ROW=$(awk -v n="$SLURM_ARRAY_TASK_ID" 'NR==n+1' "{SAMPLES_CSV}" | tr -d '\r')
[ -n "$ROW" ] || die "no row $SLURM_ARRAY_TASK_ID in {SAMPLES_CSV}"
IFS=, read -r SAMPLE FQ1 FQ2 CONDITION CROSS INDIVIDUAL <<< "$ROW"
[ -n "$INDIVIDUAL" ] || die "empty individual"
BAM="$R/bam/$SAMPLE.bam"; VCF="$R/genotypes/$INDIVIDUAL.het.vcf.gz"
OUT="$R/phaser/$SAMPLE"; LOG="$R/logs/phaser_$SAMPLE.log"; TMPD="$R/tmp/phaser_$SAMPLE"
mkdir -p "$R/phaser" "$TMPD" || die "cannot create output directories"
rm -f "$OUT".*                 # never keep outputs of an earlier run: a failed sample must have none
[ -s "$BAM" ] && [ -s "$BAM.bai" ] || die "missing $BAM or its index (run the per-sample array job first)"
[ -s "$VCF" ] && [ -s "$VCF.tbi" ] || die "missing $VCF or its index (run prep_genotypes.sh first)"
[ -s "$BED" ] || die "missing $BED (written by prep_genotypes.sh)"
NS=$(bcftools query -l "$VCF" | wc -l); VS=$(bcftools query -l "$VCF" | head -1)
[ "$NS" -eq 1 ] || die "$VCF has $NS sample columns; phASER needs exactly one"
# contig names: phaser_gene_ae splits variant IDs on '_' (a contig such as chrUn_xxx breaks it), and names that differ between
# the VCF, the BAM and the gene spans give silent zero counts
N_US=$(bcftools query -f '%CHROM\n' "$VCF" | awk '$1 ~ /_/' | wc -l)
[ "$N_US" -eq 0 ] || die "$N_US heterozygous sites lie on contigs whose names contain '_' ($(bcftools query -f '%CHROM\n' "$VCF" | awk '$1 ~ /_/' | uniq | head -3 | paste -sd' ')); phaser_gene_ae splits variant IDs on '_'. Supply genotype VCFs without these contigs (Step 6)"
BAM_CTG="$TMPD/bam_contigs.txt"
samtools view -H "$BAM" | awk -F'\t' '$1 == "@SQ" {sub(/^SN:/, "", $2); print $2}' > "$BAM_CTG" || die "cannot read the header of $BAM"
N_OFF=$(bcftools query -f '%CHROM\n' "$VCF" | awk 'NR == FNR {c[$1] = 1; next} !($1 in c)' "$BAM_CTG" - | wc -l)
[ "$N_OFF" -eq 0 ] || die "$N_OFF heterozygous sites of $VCF lie on contigs that are not in the BAM header"
N_BED=$(awk 'NR == FNR {c[$1] = 1; next} ($1 in c)' "$BAM_CTG" "$BED" | wc -l)
[ "$N_BED" -gt 0 ] || die "no gene of $BED lies on a contig of the BAM (gene contigs: $(cut -f1 "$BED" | uniq | head -3 | paste -sd' '); BAM contigs: $(head -3 "$BAM_CTG" | paste -sd' '))"
# phased genotypes: --gw_phase_vcf 1 only when the user said so (Step 7) and the VCF really is phased
N_HET=$(bcftools view -H "$VCF" | wc -l); N_PH=$(bcftools query -f '[%GT]\n' "$VCF" | grep -c '|')
if [ "$PHASED_GT" = 1 ]; then
  [ $((N_PH * 10)) -ge $((N_HET * 9)) ] || die "PHASED_GT=1 but only $N_PH of $N_HET heterozygous genotypes in $VCF are phased (0|1 or 1|0); answer 'not phased' in Step 7 or supply phased VCFs"
elif [ $((N_PH * 10)) -ge $((N_HET * 9)) ]; then
  echo "WARNING: $N_PH of $N_HET genotypes are phased but PHASED_GT=0: gene counts use the most-covered block of each gene only (Step 7)"
fi
"$PY" "$SRC/phaser/phaser.py" --vcf "$VCF" --bam "$BAM" --sample "$VS" --paired_end 1 --mapq 255 --baseq {MIN_BASEQ} \
     --pass_only 0 --unique_ids 1 --gw_phase_vcf "$PHASED_GT" --threads "${SLURM_NTASKS:-1}" --temp_dir "$TMPD" --o "$OUT" > "$LOG" 2>&1 \
  || { tail -5 "$LOG" >&2; rm -f "$OUT".*; die "phaser.py failed (log $LOG)"; }
N_BLK=$(awk 'NR > 1' "$OUT.haplotypic_counts.txt" 2>/dev/null | wc -l)
[ "$N_BLK" -gt 0 ] || { rm -f "$OUT".*; die "phASER wrote no haplotype block (check the BAM, the VCF sample $VS and the contig names; log $LOG)"; }
"$PY" "$SRC/phaser_gene_ae/phaser_gene_ae.py" --haplotypic_counts "$OUT.haplotypic_counts.txt" --features "$BED" --o "$OUT.gene_ae.txt" >> "$LOG" 2>&1 \
  || { tail -5 "$LOG" >&2; rm -f "$OUT".*; die "phaser_gene_ae.py failed (log $LOG)"; }
N_COV=$(awk -F'\t' 'NR > 1 && $7 > 0' "$OUT.gene_ae.txt" 2>/dev/null | wc -l)
[ "$N_COV" -gt 0 ] || { rm -f "$OUT".*; die "no gene has haplotype counts (totalCount > 0): the gene spans and the variants do not overlap (contig names? log $LOG)"; }
N_GW=$(awk -F'\t' 'NR > 1 && $7 > 0 && $11 == 1' "$OUT.gene_ae.txt" | wc -l)
echo "Sample $SAMPLE: $N_BLK haplotype blocks, $N_COV genes with counts ($N_GW genome-wide phased); $(grep -m1 'PHASED' "$LOG")"
rm -rf "$TMPD"
```

After the code block, add this outputs table to Step 18:

| Output | Content |
|---|---|
| `{RESULTS_DIR}/phaser/{sample}.haplotypic_counts.txt` | haplotype blocks with `aCount`, `bCount`, `blockGWPhase` |
| `{RESULTS_DIR}/phaser/{sample}.gene_ae.txt` | one row per gene of `genes_span.bed`: `aCount bCount totalCount n_variants gw_phased` (Rmd 05) |
| `{RESULTS_DIR}/phaser/{sample}.allele_config.txt`, `.variant_connections.txt`, `.haplotypes.txt`, `.allelic_counts.txt`, `.vcf.gz` | other phASER outputs (not used by the Rmds) |
| `{RESULTS_DIR}/logs/phaser_{sample}.log` | phASER and phaser_gene_ae output |

- [ ] **Step 4: Write the outbred assembler**

Create `ase-pipeline/tests/synthetic/assemble_outbred_from_skill.sh`:

```bash
#!/bin/bash
# Usage: assemble_outbred_from_skill.sh <skill.md> <values.tsv> <genotype_map.tsv>
# Builds every outbred-mode script and Rmd 01/02/05 of a project from the skill's code blocks (cut by headings or anchor lines,
# never by line numbers), as the wizard would: placeholders substituted from values.tsv (NAME<TAB>value; must define CWD,
# RESULTS_DIR, TODAY and WD_NAME); the genotype map (individual<TAB>vcf<TAB>sample in vcf) written into prep_genotypes.sh; the test-data
# settings FILTER_EXPR='' and MIN_SITES_WARN=20; the GTF bind dropped from the run scripts (the test genome lies inside CWD); the
# four mouse-helper lines of submit_chain.sh deleted. Stops on any unsubstituted placeholder.
set -uo pipefail
SKILL=${1:?}; VALS=${2:?}; MAP=${3:?}; HERE=$(cd "$(dirname "$0")" && pwd)
die() { echo "ERROR: $*" >&2; exit 1; }
val() { awk -F'\t' -v k="$1" '$1 == k {print $2; f = 1} END {exit !f}' "$VALS" || die "values file has no $1"; }
CWD=$(val CWD); RD=$(val RESULTS_DIR); TAG="$(val TODAY)_$(val WD_NAME)"; S=$RD/scripts
mkdir -p "$S" "$RD/logs" "$RD/tmp" || die "cannot create $S"
blk() { bash "$HERE/cut_block.sh" "$SKILL" "$1" || die "no code block after: $1"; }
sub() { awk -F'\t' 'NR == FNR { v[$1] = $2; next }
  { for (k in v) { t = "{" k "}"; while ((i = index($0, t)) > 0) $0 = substr($0, 1, i - 1) v[k] substr($0, i + length(t)) } print }' "$VALS" -; }
hdr() {   # $1 log name, $2 resources placeholder, $3 "array" or ""
  echo '#!/bin/bash'; echo '#SBATCH -N 1 -p bcc'; [ -n "$3" ] && echo '#SBATCH --array=1-{ARRAY_N}'; echo "#SBATCH $2"
  echo '#SBATCH --mail-type=END,FAIL'; echo '#SBATCH --mail-user={USER_EMAIL}'; echo "#SBATCH -o {RESULTS_DIR}/logs/$1"; echo 'set -uo pipefail'; }
# prep_genotypes.sh: header, block C, prep fetch, genotype VCF binds, block R, genotype loop, gene spans, block I
{ hdr 'prep_genotypes_%j.out' '{PREP_RESOURCES}' ''
  blk "### Shared block C"; blk "**Prep container fetch (both prep scripts).**"
  cut -f2 "$MAP" | xargs -n1 dirname | sort -u | sed 's/.*/add_bind "&"/'
  blk "### Shared block R"
  blk '### `prep_genotypes.sh` (outbred mode)' | sed -e "s/^FILTER_EXPR=.*/FILTER_EXPR=''/" -e 's/^MIN_SITES_WARN=.*/MIN_SITES_WARN=20/' |
    awk -v m="$MAP" '/^\{INDIVIDUAL\}/ { while ((getline l < m) > 0) print l; next } { print }'
  blk "### Gene spans for phASER"
  echo 'INDEX_FASTA="$FASTA"; STAR_INDEX="{STAR_INDEX}"; INDEX_KEY="unmasked"'
  blk "### Shared block I"; echo 'echo "prep_genotypes done"'; } | sub > "$S/prep_genotypes.sh"
# align_wasp_count.sh: header, block C with the check-only fetch_sif and its four calls, common start, WASP block, empty-table guard
{ hdr 'align_wasp_count_%A_%a.out' '{ARRAY_RESOURCES}' array
  blk "### Shared block C"
  echo 'fetch_sif() { [ -s "$1" ] || { echo "ERROR: missing container $1 (run the prep job first)" >&2; exit 1; }; }'
  echo 'fetch_sif "$STAR_SIF"; fetch_sif "$GATK_SIF"; fetch_sif "$SAMTOOLS_SIF"; fetch_sif "$PICARD_SIF"'
  blk "**Common start of both scripts**"; echo '[ -n "$INDIVIDUAL" ] || die "empty individual"'
  echo 'HET_PLAIN="$R/genotypes/$INDIVIDUAL.het.vcf"; HET_GZ="$R/genotypes/$INDIVIDUAL.het.vcf.gz"'
  blk '### `align_wasp_count.sh` (outbred mode)'; blk "### Empty-table guard"; } | sub > "$S/align_wasp_count.sh"
# Rmd run scripts from the run_01 template (first bind case: the test genome lies inside CWD)
R01=$(blk "The template below shows the second case:")
for r in 01_import_qc 02_imbalance 05_phaser; do
  printf '%s\n' "$R01" | sed -e "s/ase_01_import_qc/ase_$r/; s/run_01_import_qc_/run_${r}_/; s/_01_import_qc\.Rmd/_$r.Rmd/; s/Rmd 01 failed/Rmd ${r%%_*} failed/" \
    -e 's/--bind {CWD},{GTF_DIR} /--bind {CWD} /' | sub > "$S/run_$r.sh"
done
sed -i 's/-t 1:00:00/-t 4:00:00/' "$S/run_05_phaser.sh"     # Step 19: -n 1 --mem=16G -t 4:00:00
blk '### `setup_phaser_env.sh`' | sub > "$S/setup_phaser_env.sh"
blk '### `phaser_count.sh`' | sub > "$S/phaser_count.sh"
blk "### Submission order (one block" | awk '/^# Without the helper, delete the next four lines/ {print; skip = 4; next} skip > 0 {skip--; next} {print}' | sub > "$S/submit_chain.sh"
blk "**Waiting.**" | sub > "$S/wait_chain.sh"
for n in 12:01_import_qc 13:02_imbalance 19:05_phaser; do
  bash "$HERE/render_from_skill.sh" "$SKILL" "${n%%:*}" "$VALS" "$CWD/${TAG}_${n#*:}.Rmd" > /dev/null || die "render of Step ${n%%:*} failed"
done
if grep -nE '\{[A-Z][A-Z_0-9]*\}' "$S"/*.sh; then die "unsubstituted placeholders above"; fi
for f in "$S"/*.sh; do bash -n "$f" || die "syntax error in $f"; done
echo "assembled $(ls "$S" | wc -l) scripts in $S and 3 Rmds in $CWD"
```

- [ ] **Step 5: Add the `phaser` mode to the evaluator**

In `evaluate_stage3.R`, insert this before the final `quit(status = fail)`:

```r
if (what == "phaser") {
  for (i in seq_len(nrow(sm))) {
    s <- sm$sample[i]; ind <- sm$individual[i]
    hc <- read.delim(file.path(RES, "phaser", paste0(s, ".haplotypic_counts.txt")), stringsAsFactors = FALSE,
                     colClasses = c(variants = "character", haplotypeA = "character", haplotypeB = "character", blockGWPhase = "character"))
    ga <- read.delim(file.path(RES, "phaser", paste0(s, ".gene_ae.txt")), stringsAsFactors = FALSE,
                     colClasses = c(name = "character", gw_phased = "character", variants = "character", log2_aFC = "character"))
    ph <- tp[tp$individual == ind & tp$gt == "0/1", ]
    h1allele <- setNames(ifelse(ph$alt_on == "1", ts$alt[match(ph$pos, ts$pos)], ts$ref[match(ph$pos, ts$pos)]), ph$pos)
    multi <- hc[hc$variantCount >= 2, ]
    pos_of <- function(v) as.integer(vapply(strsplit(strsplit(v, ",")[[1]], "_"), `[`, "", 2))
    okb <- vapply(seq_len(nrow(multi)), function(j) {
      A <- strsplit(multi$haplotypeA[j], ",")[[1]]; t1 <- unname(h1allele[as.character(pos_of(multi$variants[j]))])
      !anyNA(t1) && (all(A == t1) || all(A != t1)) }, logical(1))
    ntr <- sum(vapply(seq_len(nrow(multi)), function(j) length(unique(ph$alt_on[match(pos_of(multi$variants[j]), ph$pos)])) > 1, logical(1)))
    crit(nrow(multi) >= 20 && mean(okb) >= 0.99,
         sprintf("%s: %d multi-SNP blocks (>= 20), %.4f phased as in the truth (>= 0.99)", s, nrow(multi), mean(okb)))
    crit(ntr >= 1, sprintf("%s: blocks joining SNPs in trans: %d (>= 1)", s, ntr))
    g <- merge(ga[ga$totalCount > 0, ], tf[tf$sample == s, ], by.x = "name", by.y = "gene_id")
    crit(nrow(g) > 0 && all(g$totalCount <= g$cover_h1 + g$cover_h2),
         sprintf("%s: no gene count exceeds the simulated fragments covering its het SNPs (%d genes; median ratio %.3f)", s, nrow(g),
                 median(g$totalCount / (g$cover_h1 + g$cover_h2))))
    tb <- g[g$name %in% tg$gene_id[tg$class == "two_block"], ]
    if (GT == "unphased") {
      crit(nrow(tb) == 2 && all(tb$n_variants <= 3), sprintf("%s: two_block genes use one block (n_variants %s <= 3)", s, paste(tb$n_variants, collapse = ",")))
    } else {
      crit(nrow(tb) == 2 && all(tb$n_variants == 6 & tb$gw_phased == "1"),
           sprintf("%s: two_block genes use both blocks (n_variants %s = 6, genome-wide phased)", s, paste(tb$n_variants, collapse = ",")))
      crit(all(g$gw_phased == "1"), sprintf("%s: every covered gene genome-wide phased (%d of %d)", s, sum(g$gw_phased == "1"), nrow(g)))
      h <- merge(g, tig[tig$individual == ind, c("gene_id", "h1")], by.x = "name", by.y = "gene_id")
      pl <- h[h$h1 != 0.5 & abs(h$aCount / h$totalCount - 0.5) > 0.1, ]
      ag <- sign(pl$aCount / pl$totalCount - 0.5) == sign(pl$h1 - 0.5)
      crit(nrow(pl) > 0 && all(ag), sprintf("%s: aCount is haplotype 1 (left GT allele) in %d of %d planted genes", s, sum(ag), nrow(pl)))
    }
    lo <- g[g$name %in% tg$gene_id[tg$class == "hap_lowdepth"], ]
    cat(sprintf("REPORT %s: hap_lowdepth genes, SNPs used: %s of 8\n", s, paste(lo$n_variants, collapse = ",")))
  }
}
```

- [ ] **Step 6: Run the checker (expected PASS), then smoke-run on aligned synthetic data, both genotype modes**

```bash
SCR=/net/bmc-lab3/data/bcc/yannvrb/ase_stage3_scratch; SYN=/net/bmc-lab3/data/bcc/yannvrb/ase_synthetic_stage3
C=/home/yannvrb/.singularity/cache/depot.galaxyproject.org-singularity
mkproj() {   # $1 project, $2 PHASED_GT, $3 VCF suffix (all.vcf | all.phased.vcf)
  local P=$SCR/$1 RD=$SCR/$1/results/2026-10-01_$1
  rm -rf $P; mkdir -p $P; cp -r $SYN/genome $P/genome; cp $SYN/outbred_phase/samples.csv $P/$1_samples.csv
  for i in 1 2 3 4 5 6; do printf 'ind%s\t%s/outbred_phase/ind%s.%s\tind%s\n' $i $SYN $i $3 $i; done > $P/genotype_map.tsv
  printf '%s\t%s\n' CWD $P WD_NAME $1 TODAY 2026-10-01 RESULTS_DIR $RD MODE outbred STRAIN_A REF STRAIN_B ALT \
    SAMPLES_CSV $P/$1_samples.csv GENOME_DIR $P/genome FASTA_PATH $P/genome/genome.fa GTF_PATH $P/genome/genome.gtf GTF_DIR $P/genome \
    STAR_INDEX $P/genome/index/star_ase_sjdb99 SJDB_OVERHANG 99 \
    STAR_SIF $C-star-2.7.10b--h9ee0642_0.img GATK_SIF $C-gatk4-4.4.0.0--py36hdfd78af_0.img BCFTOOLS_SIF $C-bcftools-1.20--h8b25389_0.img \
    SAMTOOLS_SIF $C-samtools-1.21--h50ea8bc_0.img PICARD_SIF $C-picard-3.1.1--hdfd78af_0.img \
    R_SIF /net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif USER_EMAIL yannvrb@mit.edu \
    PREP_RESOURCES "-n 4 --mem=8G -t 0:30:00" ARRAY_RESOURCES "-n 4 --mem=8G -t 1:00:00" PHASER_RESOURCES "-n 2 --mem=4G -t 0:30:00" \
    ARRAY_N 6 MIN_MAPQ 10 MIN_BASEQ 10 MIN_DEPTH 10 FDR_SIG 0.05 ABS_DEV_SIG 0.1 BIAS_TOL 0.03 RHO_MIN 0.01 THIN_BP 500 \
    PROJECT_TITLE "Stage 3 test" AUTHOR tester ANALYSES "per-sample phaser" PHASER_HOME $SCR/phaser_home PHASED_GT $2 REF_CONDITION ctrl > $P/values.tsv
  bash ase-pipeline/tests/synthetic/assemble_outbred_from_skill.sh ase-pipeline/ase-pipeline.md $P/values.tsv $P/genotype_map.tsv
}
submit() {   # $1 project: prep -> array -> (phaser, Rmd 01 -> Rmd 02) -> Rmd 05, by hand (the chain integration is Task 6)
  local S=$SCR/$1/results/2026-10-01_$1/scripts
  P=$(sbatch -p bcc --parsable $S/prep_genotypes.sh); A=$(sbatch -p bcc --parsable --dependency=afterok:$P $S/align_wasp_count.sh)
  PH=$(sbatch -p bcc --parsable --dependency=afterok:$A $S/phaser_count.sh)
  R1=$(sbatch -p bcc --parsable --dependency=afterok:$A $S/run_01_import_qc.sh); R2=$(sbatch -p bcc --parsable --dependency=afterok:$R1 $S/run_02_imbalance.sh)
  R5=$(sbatch -p bcc --parsable --dependency=afterok:$PH:$R2 $S/run_05_phaser.sh); echo "$1: prep $P array $A phaser $PH Rmd01 $R1 Rmd02 $R2 Rmd05 $R5"
}
mkproj t5_unph 0 all.vcf && mkproj t5_ph 1 all.phased.vcf && submit t5_unph && submit t5_ph
```

When every job has left the queue (`squeue -u $USER`; read the states with `scontrol show job <id>`), run an evaluation job: `evaluate_stage3.R phaser` and `evaluate_stage3.R rmd` for `t5_unph` (unphased) and for `t5_ph` (phased), in the R container, against `$SYN/outbred_phase`.

Expected:
- every job exits 0, and the prep log shows `Gene spans for phASER: 80 genes, 80 on contigs of the FASTA`;
- every `phaser_count` log has the `Sample ...: N haplotype blocks` line;
- every evaluator line is `PASS`, under the gates of this task and Task 4.

In particular, the phased project shows every covered gene genome-wide phased. That proves STAR, ASEReadCounter and phASER accept phased genotypes. If STAR or ASEReadCounter fails on a phased `0|1` VCF, report BLOCKED with the log; it would be a new fact.

Record the `REPORT` lines, and the per-sample phASER wall time from the `phaser_count` logs (job start and end from `scontrol`).

- [ ] **Step 7: Negative tests (Review Focus 1, 2 and 5)**

Work on a copy `$SCR/t5_neg` of `$SCR/t5_unph`. Regenerate `phaser_count.sh` with the assembler and a `values.tsv` whose `RESULTS_DIR` and `CWD` point at the copy. Run one array task with `sbatch -p bcc --array=1 <script>`.

- (a) **Contig mismatch:** `sed -i 's/^chr1\t/1\t/' .../reference/genes_span.bed`. Expected: exit 1; the log has `no gene of ... lies on a contig of the BAM`; no `phaser/obp_s1.*` file remains.
- (b) **False phased claim:** set `PHASED_GT=1` in the copied script, which still uses the unphased VCFs. Expected: exit 1 with `PHASED_GT=1 but only 0 of`.
- (c) **Contig with `_`:** in a copy of `genotypes/ind1.het.vcf.gz`, rename `chr1` to `chr1_x` (`bcftools view` piped through `sed`, `bgzip`, `tabix`, using the environment's tools in a job). Expected: exit 1 with `contigs whose names contain '_'`.
- (d) **Stale output:** before test (a), create a stale `phaser/obp_s1.gene_ae.txt` (`echo stale > ...`). Expected after (a): the file is gone.

- [ ] **Step 8: Prove the phASER flag rule on `phaser_count.sh`, then commit**

Proofs, each on a scratch copy of the skill; the checker must FAIL:
- `--paired_end 1` replaced by `--paired-end 1` → FAIL naming `--paired-end`;
- ` --gw_phase_vcf "$PHASED_GT"` replaced by ` --gw_phase "$PHASED_GT"` → FAIL naming `--gw_phase`;
- `--haplotypic_counts` replaced by `--haplotype_counts` → FAIL naming `--haplotype_counts`.

Then run the checker on the working copy (PASS) and on `git show HEAD:...` (every new `need` fails).

```bash
git status && git diff --stat
git add ase-pipeline/ase-pipeline.md ase-pipeline/tests/check_skill.sh ase-pipeline/tests/synthetic/assemble_outbred_from_skill.sh ase-pipeline/tests/synthetic/evaluate_stage3.R
git commit -m "ase-pipeline: phaser_count.sh per-sample phASER array with contig, phase and output guards; gene-span BED; outbred assembler

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Wizard, submission chain, wait, summary page and Notes (Steps 2, 7, 9, 15, Notes), with stub-scheduler dry runs

**Files:**
- Modify: `ase-pipeline/ase-pipeline.md` (Steps 2, 7, 9, 15 and Notes)
- Create: `ase-pipeline/tests/chain/dry_run_chain.sh`
- Modify: `ase-pipeline/tests/check_skill.sh`

**Interfaces:**
- Consumes:
  - the script names `setup_phaser_env.sh`, `phaser_count.sh` and `run_05_phaser.sh`;
  - the Rmd 05 files `{TODAY}_{WD_NAME}_05_phaser.html`, `_ASE_phaser.xlsx`, `_ASE_phaser.pdf` and `summary_numbers_phaser.tsv` (Task 4 columns);
  - the installation contract (Task 3);
  - `cut_block.sh` (Task 3).
- Produces:
  - the `{ANALYSES}` word `phaser`;
  - the placeholders `{PHASED_GT}` (0/1), `{PHASER_HOME}` (a path, or the word `none`) and `{PHASER_RESOURCES}`;
  - `bash ase-pipeline/tests/chain/dry_run_chain.sh <skill.md> <scratch dir>`, which prints `DRY RUN PASS` or exits 1.

- [ ] **Step 1: Failing checker lines, and retire the replaced ones**

In `check_skill.sh`, delete these lines (each text now changes):
- `need "available in a later stage"` (the Stage 1 Task 4 block);
- `need "offered only when \`{MODE}\` = \`outbred\`; available in a later stage"`;
- `need "pasted verbatim into Rmd 02, 03 and 04"`;
- the line `need 'for w in $ANALYSES; do case "$w" in per-sample|reciprocal|differential) ;; *) die "ANALYSES contains '"'"'$w'"'"'; the accepted words are per-sample, reciprocal and differential'`.

Then append:

```bash
# --- Stage 3 Task 6 (wizard, chain, wait, summary, Notes); each check fails on the Task 5 skill
forbid "available in a later stage"
need "offered only when \`{MODE}\` = \`outbred\` and the data are paired-end"
need 'menu number 1 is the word `per-sample`, 2 is `reciprocal`, 3 is `differential`, 4 is `phaser`'
need "Are the genotype VCFs phased?"
need "Store it as \`{PHASED_GT}\`"
need "Where should phASER be installed?"
need 'for w in $ANALYSES; do case "$w" in per-sample|reciprocal|differential|phaser) ;; *) die "ANALYSES contains '"'"'$w'"'"'; the accepted words are per-sample, reciprocal, differential and phaser'
need 'case " $ANALYSES " in *" phaser "*) [ "$MODE" = outbred ] || die "phASER is outbred only"'
pair_need '  case " $NEED " in *" setup_phaser_env.sh "*) PS=$(sbatch -p bcc --parsable $S/setup_phaser_env.sh); got setup_phaser_env.sh "$PS" ;; esac' \
          '  PH=$(sbatch -p bcc --parsable --dependency=afterok:$A${PS:+:$PS} $S/phaser_count.sh); got phaser_count.sh "$PH"'
pair_need '  PH=$(sbatch -p bcc --parsable --dependency=afterok:$A${PS:+:$PS} $S/phaser_count.sh); got phaser_count.sh "$PH"' \
          '  R5=$(sbatch -p bcc --parsable --dependency=afterok:$PH:$R2 $S/run_05_phaser.sh); got run_05_phaser.sh "$R5" ;;'
need 'run_05_*) echo summary_numbers_phaser.tsv ;;'
need "### phASER section (only when Rmd 05 ran)"
need "`{PHASER_RESOURCES}` for `phaser_count.sh`"
need "pasted verbatim into Rmd 02, 03, 04 and 05"
need "the one exception is \`module add miniconda3/v4\`"
[ -s "$HERE/chain/dry_run_chain.sh" ] || { echo "FAIL: missing tests/chain/dry_run_chain.sh"; fail=1; }
# --- end Stage 3 Task 6
```

`pair_need` is defined earlier in the checker (Stage 2 final-review block). It requires the second string on the line directly after a line that contains the first.

Run the checker. Expected: FAIL.

- [ ] **Step 2: Step 2 (environment)**

Replace "No conda environment is needed for this skill." with "No conda environment is needed for this skill, except for the optional phASER analysis (Steps 18-19), which installs its own pinned environment in a compute job; nothing of it runs on the login node."

- [ ] **Step 3: Step 7 (menu, word mapping, phASER questions)**

1. Replace menu item 4 with:

   "4. **phASER haplotype phasing** (Rmd 05, Steps 18 and 19): offered only when `{MODE}` = `outbred` and the data are paired-end (phASER runs with `--paired_end 1`, the only setting tested). It phases each sample's heterozygous SNPs from the reads and tests gene-level haplotype counts, reported next to (never instead of) the unphased results of Rmd 02."

2. In the paragraph after the list:
   - replace "menu number 1 is the word `per-sample`, 2 is `reciprocal`, 3 is `differential` (`per-sample` is always present); for example the answer "1,2,3" is stored as `per-sample reciprocal differential`" with "menu number 1 is the word `per-sample`, 2 is `reciprocal`, 3 is `differential`, 4 is `phaser` (`per-sample` is always present); for example the answer "1,2,3" is stored as `per-sample reciprocal differential` and "1,4" as `per-sample phaser`";
   - delete the sentence "Menu number 4 (phASER) runs nothing yet and is not stored in `{ANALYSES}`;" and the sentence "If the user chooses phASER, say "available in a later stage" and continue.";
   - change "stops when `{ANALYSES}` holds any other word" to "stops when `{ANALYSES}` holds any word other than these four".

3. Add after that paragraph:

```markdown
If phASER is selected, ask two more things:

- "Are the genotype VCFs phased? 1. Yes: phased genotypes (written `0|1` / `1|0`, for example WGS with statistical or trio phasing) · 2. No, or not sure". Store it as `{PHASED_GT}` (1 or 0). The rnavar VCFs of Step 6 option 2 are unphased: use 0. With 0, tell the user: "With unphased genotypes, phASER's gene step keeps only the most-covered haplotype block of each gene (on the synthetic test one sample used 113 of its 153 heterozygous SNPs), and the haplotype labels A and B are arbitrary per gene and sample, so Rmd 05 reports the size of the imbalance without a direction. The unphased results of Rmd 02 stay the main result." `phaser_count.sh` (Step 18) checks the answer against the VCFs and stops when phased genotypes were claimed but fewer than 90 percent are phased.
- "Where should phASER be installed? (a directory you can write to; about 650 MB; reused by later projects)". Store the path as `{PHASER_HOME}` (Step 0 rule: no space or comma). If `{PHASER_HOME}/install_ok.txt` exists and its first line is `commit aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301`, say "phASER is already installed there"; otherwise say "the chain first runs a one-time installation job (about 20-25 minutes; Step 18)".

When phASER is not selected, set `{PHASED_GT}` to 0 and `{PHASER_HOME}` to the word `none`.
```

- [ ] **Step 4: Step 9 (resources)**

Add after the array-job tiers:

"`{PHASER_RESOURCES}` for `phaser_count.sh` (phASER only; one array task per sample):
- small genome (under 100 Mb): `-n 2 --mem=4G -t 0:30:00`;
- otherwise: `-n 4 --mem=32G -t 2:00:00`.

The larger request is an estimate from upstream's benchmark (88-466 s per real sample on one thread). It is not verified on real data here. The installation job requests `-n 4 --mem=16G -t 2:00:00` (Step 18)."

- [ ] **Step 5: Step 15 (submission chain, wait, summary page)**

1. **Paragraph before the block:**
   - change "and Steps 16 and 17 for the selected analyses" to "and Steps 16-19 for the selected analyses";
   - change "holds any word other than `per-sample`, `reciprocal`, `differential` (for example menu numbers)" to "holds any word other than `per-sample`, `reciprocal`, `differential`, `phaser` (for example menu numbers), and when `phaser` is selected in F1 mode".
2. **`submit_chain.sh` block, edits in order:**
   - Change `MODE="{MODE}"; ANALYSES="{ANALYSES}"   # MODE: f1 or outbred; ANALYSES: the words selected in Step 7, e.g. "per-sample reciprocal differential"` to `MODE="{MODE}"; ANALYSES="{ANALYSES}"; PHASER_HOME="{PHASER_HOME}"   # Step 7: MODE f1 or outbred; ANALYSES the words, e.g. "per-sample phaser"; PHASER_HOME a path or none`.
   - Change the word-guard text `write the Step 7 words (per-sample, reciprocal, differential)` to `write the Step 7 words (per-sample, reciprocal, differential, phaser)`.
   - Replace the line `for w in $ANALYSES; do case "$w" in per-sample|reciprocal|differential) ;; *) die "ANALYSES contains '$w'; the accepted words are per-sample, reciprocal and differential, separated by spaces (not menu numbers or other spellings)" ;; esac; done` with `for w in $ANALYSES; do case "$w" in per-sample|reciprocal|differential|phaser) ;; *) die "ANALYSES contains '$w'; the accepted words are per-sample, reciprocal, differential and phaser, separated by spaces (not menu numbers or other spellings)" ;; esac; done`.
   - After the line `case " $ANALYSES " in *differential*) NEED="$NEED run_04_differential.sh" ;; esac`, insert:

     ```bash
     case " $ANALYSES " in *" phaser "*) [ "$MODE" = outbred ] || die "phASER is outbred only"; NEED="$NEED phaser_count.sh run_05_phaser.sh"
       [ "$(sed -n 1p "$PHASER_HOME/install_ok.txt" 2>/dev/null)" = "commit aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301" ] || NEED="$NEED setup_phaser_env.sh" ;; esac
     ```
   - Change `DEP=""; R3="not selected"; R4="not selected"` to `DEP=""; R3="not selected"; R4="not selected"; R5="not selected"; PS=""`.
   - After the Rmd 04 `case ... esac` lines, before `echo "submitted (script, job id):"`, insert:

     ```bash
     # 7. phASER only when selected (Step 7, outbred): one-time installation (only when PHASER_HOME has no installation of the pinned
     #    commit; it runs beside the prep job), the per-sample phASER array after the alignment array, then Rmd 05 after phASER and Rmd 02
     case " $ANALYSES " in *" phaser "*)
       case " $NEED " in *" setup_phaser_env.sh "*) PS=$(sbatch -p bcc --parsable $S/setup_phaser_env.sh); got setup_phaser_env.sh "$PS" ;; esac
       PH=$(sbatch -p bcc --parsable --dependency=afterok:$A${PS:+:$PS} $S/phaser_count.sh); got phaser_count.sh "$PH"
       R5=$(sbatch -p bcc --parsable --dependency=afterok:$PH:$R2 $S/run_05_phaser.sh); got run_05_phaser.sh "$R5" ;;
     esac
     ```
   - Change the last line to `echo "Rmd 03: $R3; Rmd 04: $R4; Rmd 05: $R5"`.
3. **The "Do not edit the block except for the mouse helper" paragraph:**
   - add `{PHASER_HOME}` (the Step 7 path, or `none`) to the substituted placeholders, and the word `phaser` to the list of words;
   - add the sentence: "With phASER, `submit_chain.sh` adds the installation job only when `{PHASER_HOME}/install_ok.txt` does not hold the pinned commit; `phaser_count.sh` then waits for both the alignment array and the installation. Rmd 05 depends on `phaser_count.sh` and on Rmd 02 (it compares with Rmd 02's gene results), so a failure of either leaves it cancelled."
4. **`wait_chain.sh`:** replace the two `result_of` lines with:

   ```bash
   result_of() { case $1 in run_01_*) echo ase_checkpoint.rds ;; run_02_*) echo summary_numbers.tsv ;;
                            run_03_*) echo summary_numbers_reciprocal.tsv ;; run_04_*) echo summary_numbers_differential.tsv ;;
                            run_05_*) echo summary_numbers_phaser.tsv ;; esac; }
   ```

   In the "**Waiting.**" prose, change "`finished` for the prep, helper and array jobs" to "`finished` for the prep, helper, array, phASER installation and phASER array jobs".
5. **Failure paragraphs:**
   - Change "A failure of Rmd 03 or Rmd 04 does not invalidate the per-sample results" to "A failure of Rmd 03, Rmd 04, Rmd 05 or of the phASER jobs does not invalidate the per-sample results".
   - In "Re-submitting part of the chain", after "then submit Rmd 01, Rmd 02 and the selected Rmd 03 / Rmd 04 again", add "(and `phaser_count.sh` and Rmd 05 when phASER was selected; `phaser_count.sh` without its dependency on the old array job once the BAMs exist)".
   - In "Dropping a failed sample", add: "Rmd 05 reads the checkpoints of Rmd 01 and Rmd 02 and the phASER table of every sample in `samples.csv`; re-run it after Rmd 02. `phaser_count.sh` need not run again for the remaining samples."
6. **Summary report:**
   - Change "If Rmd 03 or Rmd 04 was selected but failed" to "If Rmd 03, Rmd 04 or Rmd 05 was selected but failed". Add `run_05_phaser_<jobid>.out` to the log names, and add: "when `phaser_count.sh` failed (Rmd 05 then shows as CANCELLED), quote the first error line of its log `phaser_count_<jobid>_<task>.out`".
   - Add after the "### Differential ASE section" subsection:

```markdown
### phASER section (only when Rmd 05 ran)

Read `{RESULTS_DIR}/summary_numbers_phaser.tsv` with `awk -F'\t'` (columns `sample, genotypes, genes_with_counts, genes_tested, sig_genes, gw_phased_genes, sig_unphased, sig_both, sig_phaser_only, sig_phaser_only_untested_unphased, sig_unphased_only, rho_used, rho_own, rho_cohort`) and show one row per sample: genes with haplotype counts, genes tested, significant (phASER), significant (unphased, Rmd 02), both, phASER only (of which: no SNP tested unphased), unphased only, genome-wide phased genes, `rho_used`. Add, by the `genotypes` column: unphased — "Genotypes unphased: each gene's counts come from its most-covered haplotype block, and the haplotype labels are arbitrary per gene and sample, so no direction is given."; phased — "Genotypes phased: in genome-wide phased genes haplotype A is the haplotype of the first allele of the phased genotype; haplotype A of one individual is unrelated to haplotype A of another." Add: "The unphased results of Rmd 02 remain the main per-sample result; phASER is shown next to them. Differential ASE is not run on phASER counts." Add link cards to `{TODAY}_{WD_NAME}_05_phaser.html`, `{TODAY}_{WD_NAME}_ASE_phaser.xlsx` and `{TODAY}_{WD_NAME}_ASE_phaser.pdf`.
```

   - In "**Verify before finishing.**", extend the list of reported paths with "and, when it ran, the Rmd 05 HTML file, xlsx file and `summary_numbers_phaser.tsv`".

- [ ] **Step 6: Notes for the assistant**

- **Containers, not modules:** in that bullet, replace "the only module ever loaded is `singularity/3.10.4`, always with a checked `|| exit 1`" with "the only module loaded for the tools is `singularity/3.10.4`, always with a checked `|| exit 1`; the one exception is `module add miniconda3/v4` in the one-time phASER installation job (`setup_phaser_env.sh`, Step 18), because conda is not available as a container here."
- **Statistics block:** replace "is pasted verbatim into Rmd 02, 03 and 04" with "is pasted verbatim into Rmd 02, 03, 04 and 05".
- **Analyses available:** replace that bullet with: "- **Analyses available.** Per-sample imbalance (Rmd 01/02), reciprocal F1 (Rmd 03), differential ASE (Rmd 04) and, in outbred mode with paired-end data, phASER haplotype counts (Steps 18-19, Rmd 05, reported next to Rmd 02, never instead of it). Differential ASE on phASER counts is not offered. Never add a model or package that the Step 14 block does not contain. `lme4` and `aod` are not used by the Rmds (the Step 14 models are base R; Step 2 lists the mandatory packages)."
- **New bullet:** "- **phASER.** Always the repository https://github.com/secastel/phaser at the pinned commit (Step 18), never the bioconda package `phaser` (Python 2.7). Every phASER script exports `PYTHONNOUSERSITE=1` and uses `{PHASER_HOME}/env/bin` directly (no activation, no module). Change the commit only as Step 18 describes (all three places, fixtures, acceptance)."

- [ ] **Step 7: Write the dry-run test**

Create `ase-pipeline/tests/chain/dry_run_chain.sh`:

```bash
#!/bin/bash
# Usage: dry_run_chain.sh <skill.md> <scratch dir>
# Cuts submit_chain.sh and wait_chain.sh out of Step 15 (cut_block.sh), edits them as the wizard would (placeholders; the four
# mouse-helper lines deleted), and runs them against stub sbatch / squeue / scancel. Each scenario compares the submitted jobs and
# their afterok dependencies (as script names) with the expected list, or the stop message. Prints "DRY RUN PASS" or exits 1.
set -u
SKILL=${1:?}; W=${2:?}; HERE=$(cd "$(dirname "$0")" && pwd); CUT="$HERE/../synthetic/cut_block.sh"; fail=0
COMMIT=aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301
rm -rf "$W"; mkdir -p "$W/bin" || exit 1
cat > "$W/bin/sbatch" <<'EOF'
#!/bin/bash
# stub: prints the next job id and logs "<id> <dependency or -> <script>"; refuses a missing script and an empty id in a dependency
n=$(( $(cat "$STUB_DIR/n" 2>/dev/null || echo 100) + 1 )); echo "$n" > "$STUB_DIR/n"
dep=-; script=""
for a in "$@"; do case "$a" in --dependency=*) dep=${a#--dependency=} ;; -p|bcc|--parsable) ;; *) script=$a ;; esac; done
case "$dep" in *afterok:|*::*|*:) echo "stub sbatch: empty id in dependency '$dep'" >&2; exit 1 ;; esac
[ -s "$script" ] || { echo "stub sbatch: no such script $script" >&2; exit 1; }
echo "$n $dep $(basename "$script")" >> "$STUB_DIR/calls"; echo "$n"
EOF
printf '#!/bin/bash\nexit 0\n' > "$W/bin/squeue"; printf '#!/bin/bash\nexit 0\n' > "$W/bin/scancel"; chmod +x "$W/bin/"*
SUB=$(bash "$CUT" "$SKILL" "### Submission order (one block") || { echo "FAIL: no submit_chain.sh block"; exit 1; }
WAIT=$(bash "$CUT" "$SKILL" "**Waiting.**") || { echo "FAIL: no wait_chain.sh block"; exit 1; }
named() {   # calls file -> "script afterok:<names>" lines
  awk '{ id[$1] = $3; d = $2
         if (d != "-") { sub(/^afterok:/, "", d); n = split(d, x, ":"); d = "afterok:"; for (i = 1; i <= n; i++) d = d (i > 1 ? ":" : "") id[x[i]] }
         print $3, d }' "$1"; }
scenario() {   # $1 name, $2 MODE, $3 ANALYSES, $4 installed yes|no, $5 expected lines or "DIE:<text>", $6 scripts left out
  local D="$W/$1" got out rc; mkdir -p "$D/res/scripts" "$D/res/logs" "$D/home"
  for s in prep_genotypes.sh align_wasp_count.sh prep_f1_reference.sh align_count_f1.sh run_01_import_qc.sh run_02_imbalance.sh \
           run_03_reciprocal.sh run_04_differential.sh setup_phaser_env.sh phaser_count.sh run_05_phaser.sh; do
    case " $6 " in *" $s "*) ;; *) echo '#!/bin/bash' > "$D/res/scripts/$s" ;; esac
  done
  [ "$4" = yes ] && echo "commit $COMMIT" > "$D/home/install_ok.txt"
  printf '%s\n' "$SUB" | awk '/^# Without the helper, delete the next four lines/ {print; skip = 4; next} skip > 0 {skip--; next} {print}' |
    sed -e "s|{RESULTS_DIR}|$D/res|g" -e "s|{CWD}|$D|g" -e "s|{MODE}|$2|g" -e "s|{ANALYSES}|$3|g" -e "s|{PHASER_HOME}|$D/home|g" > "$D/submit_chain.sh"
  out=$(cd "$D" && STUB_DIR="$D" PATH="$W/bin:$PATH" bash "$D/submit_chain.sh" 2>&1); rc=$?
  case "$5" in
    DIE:*) if [ $rc -ne 0 ] && printf '%s' "$out" | grep -qF -- "${5#DIE:}" && [ ! -s "$D/calls" ]; then echo "ok   $1: stops, nothing submitted"
           else echo "FAIL: $1: expected a stop with '${5#DIE:}' and no submission (rc $rc): $out"; fail=1; fi ;;
    *) got=$([ -s "$D/calls" ] && named "$D/calls")
       if [ $rc -eq 0 ] && [ "$got" = "$5" ]; then echo "ok   $1: $(printf '%s\n' "$got" | wc -l) jobs with the expected dependencies"
       else echo "FAIL: $1 (rc $rc)"; printf '  expected:\n%s\n  got:\n%s\n  output: %s\n' "$5" "$got" "$out"; fail=1; fi ;;
  esac
}
BASE_OB='prep_genotypes.sh -
align_wasp_count.sh afterok:prep_genotypes.sh
run_01_import_qc.sh afterok:align_wasp_count.sh
run_02_imbalance.sh afterok:run_01_import_qc.sh'
scenario ob_phaser_installed outbred "per-sample phaser" yes "$BASE_OB
phaser_count.sh afterok:align_wasp_count.sh
run_05_phaser.sh afterok:phaser_count.sh:run_02_imbalance.sh" ""
scenario ob_phaser_install outbred "per-sample phaser" no "$BASE_OB
setup_phaser_env.sh -
phaser_count.sh afterok:align_wasp_count.sh:setup_phaser_env.sh
run_05_phaser.sh afterok:phaser_count.sh:run_02_imbalance.sh" ""
scenario ob_diff_phaser outbred "per-sample differential phaser" yes "$BASE_OB
run_04_differential.sh afterok:run_01_import_qc.sh
phaser_count.sh afterok:align_wasp_count.sh
run_05_phaser.sh afterok:phaser_count.sh:run_02_imbalance.sh" ""
scenario ob_no_phaser outbred "per-sample" no "$BASE_OB" "setup_phaser_env.sh phaser_count.sh run_05_phaser.sh"
scenario f1_all f1 "per-sample reciprocal differential" no 'prep_f1_reference.sh -
align_count_f1.sh afterok:prep_f1_reference.sh
run_01_import_qc.sh afterok:align_count_f1.sh
run_02_imbalance.sh afterok:run_01_import_qc.sh
run_03_reciprocal.sh afterok:run_01_import_qc.sh
run_04_differential.sh afterok:run_01_import_qc.sh' ""
scenario f1_phaser f1 "per-sample phaser" yes "DIE:phASER is outbred only" ""
scenario menu_numbers outbred "per-sample 4" yes "DIE:ANALYSES contains '4'" ""
scenario missing_rmd05 outbred "per-sample phaser" yes "DIE:MISSING in" "run_05_phaser.sh"
scenario missing_setup outbred "per-sample phaser" no "DIE:MISSING in" "setup_phaser_env.sh"
wait_case() {   # $1 name, $2 yes if summary_numbers_phaser.tsv is newer than the job list, $3 expected word for run_05_phaser.sh
  local D="$W/$1" out; mkdir -p "$D/res/logs"
  printf 'run_02_imbalance.sh\t201\nphaser_count.sh\t202\nrun_05_phaser.sh\t203\n' > "$D/res/logs/chain_job_ids.tsv"
  touch -d '2 minutes ago' "$D/res/logs/chain_job_ids.tsv"; touch "$D/res/summary_numbers.tsv"
  [ "$2" = yes ] && touch "$D/res/summary_numbers_phaser.tsv"
  printf '%s\n' "$WAIT" | sed -e "s|{RESULTS_DIR}|$D/res|g" > "$D/wait_chain.sh"
  out=$(PATH="$W/bin:$PATH" bash "$D/wait_chain.sh" 2>&1)
  if printf '%s\n' "$out" | grep -q "^run_05_phaser.sh 203 $3" && printf '%s\n' "$out" | grep -q "^phaser_count.sh 202 finished"
  then echo "ok   $1: run_05_phaser.sh $3"; else echo "FAIL: $1: $out"; fail=1; fi
}
wait_case wait_rmd05_ok yes ok
wait_case wait_rmd05_failed no FAILED
[ $fail -eq 0 ] && echo "DRY RUN PASS" || exit 1
```

Run: `bash ase-pipeline/tests/chain/dry_run_chain.sh ase-pipeline/ase-pipeline.md $SCR/dry`
Expected: 11 `ok` lines and `DRY RUN PASS`. This runs only `bash` and stubs, so it may run on the login node.

Mutation proofs, each on a scratch copy of the skill; `dry_run_chain.sh` must exit 1:
- (m1) `sed 's/per-sample|reciprocal|differential|phaser)/per-sample|reciprocal|differential)/'` → the phaser scenarios fail.
- (m2) `sed 's/--dependency=afterok:$PH:$R2 $S\/run_05_phaser.sh/--dependency=afterok:$PH $S\/run_05_phaser.sh/'` → `ob_phaser_installed` fails.
- (m3) delete the line `run_05_*) echo summary_numbers_phaser.tsv ;; esac; }` and restore the previous closing `esac; }` on the line before it → `wait_rmd05_ok` fails.

- [ ] **Step 8: Test the summary-page instructions on Task 5's outputs**

In `$SCR/t5_unph/results/2026-10-01_t5_unph/` and in `$SCR/t5_ph/...`, write `t5_unph_summary_report.html` and `t5_ph_summary_report.html` exactly as Step 15 instructs: bash plus awk, no R, reading `summary_numbers.tsv` and `summary_numbers_phaser.tsv`. Run the Step 15 href loop and the no-`http`/no-`/` check. Expected: both print nothing. Keep both pages for Task 8.

- [ ] **Step 9: Run the checker (expected PASS; every new line proven RED on `git show HEAD`), commit**

```bash
bash ase-pipeline/tests/check_skill.sh ase-pipeline/ase-pipeline.md
git status && git diff --stat
git add ase-pipeline/ase-pipeline.md ase-pipeline/tests/check_skill.sh ase-pipeline/tests/chain/dry_run_chain.sh
git commit -m "ase-pipeline: phASER in the wizard (menu word, phased-genotype and install questions), chain (install, phaser array, Rmd 05), wait, summary page; stub dry runs

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 7: README, root README and checker README block

**Files:**
- Modify: `ase-pipeline/README.md`, root `README.md` (the `/ase-pipeline` row), `ase-pipeline/tests/check_skill.sh`

**Interfaces:**
- Consumes:
  - the measured numbers of `$SCR/logs/t1_final.log`;
  - the Task 3 installation record (`ase-pipeline/tests/fixtures/verification.md`);
  - the Task 4 and Task 5 evaluator logs, including the `REPORT` lines;
  - the dry-run result.
- Produces: README text that Task 8 completes, with the marker line `Stage 3 synthetic acceptance run: PENDING`.

- [ ] **Step 1: Failing checker lines, and retire the replaced ones**

In `check_skill.sh`:
- Delete `grep -qF "(Stages 1 and 2)" "$RD" || { echo "FAIL: README title must say Stages 1 and 2"; fail=1; }` and `rneed "phASER (Stage 3) is not implemented"`.
- In the line `grep -qF "Stage 2 synthetic acceptance run: DONE (2026-09-30)" "$RD" && ! grep -qF "synthetic acceptance run: PENDING" "$RD" || ...`, change `! grep -qF "synthetic acceptance run: PENDING"` to `! grep -qF "Stage 2 synthetic acceptance run: PENDING"`.

Then append:

```bash
# --- Stage 3 README and registration; each check fails on the Task 6 README and root README
RD="$HERE/../README.md"
rneed "(Stages 1 to 3)"
rneed "**What Stage 3 adds**"
rneed "phASER (Rmd 05)"
rneed "aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301"
rneed "never the bioconda package"
rneed "PYTHONNOUSERSITE=1"
rneed "summary_numbers_phaser.tsv"
rneed "| 18 | phASER: installation and per-sample haplotype counts |"
rneed "| 19 | Rmd 05: phASER gene-level haplotype imbalance |"
rneed "haplotype labels are arbitrary"
rneed "most-covered haplotype block"
rneed "not verified on real data"
rneed "Fast Beta"
rneed "module add miniconda3/v4"
rforbid "phASER (Stage 3) is not implemented"
rforbid "the only module the skill ever loads"
grep -qE "Stage 3 synthetic acceptance run: (PENDING|DONE \(2026-[0-9-]+\))" "$RD" || { echo "FAIL: README must state the Stage 3 acceptance status"; fail=1; }
grep -q "/ase-pipeline.*phASER haplotype counts" "$HERE/../../README.md" || { echo "FAIL: root README row must mention phASER haplotype counts"; fail=1; }
! grep -q "/ase-pipeline.*phASER not included" "$HERE/../../README.md" || { echo "FAIL: root README row still says phASER not included"; fail=1; }
# --- end Stage 3 README and registration
```

Run the checker. Expected: FAIL.

- [ ] **Step 2: Edit `ase-pipeline/README.md`**

1. **Title and introduction:**
   - Title: "# `/ase-pipeline` — Allele-Specific Expression Skill (Stages 1 to 3)".
   - Replace the "Stages 1 and 2 cover ..." paragraph with: "Stages 1 to 3 cover per-sample allelic imbalance, the reciprocal F1 analysis (strain versus parent-of-origin effects), differential ASE between conditions and, in outbred mode, phASER read-backed phasing with gene-level haplotype counts (reported next to the unphased results, never instead of them)."
2. **New section** after the Stage 2 requirements, "**What Stage 3 adds**", with 4 bullets:
   - "**phASER (Rmd 05)**": a per-sample array runs phASER (pinned commit `aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301`, repository secastel/phaser, a "Fast Beta" with no releases) on the WASP-filtered BAM and the same heterozygous VCF; `phaser_gene_ae` gives gene-level haplotype counts; Rmd 05 tests them per sample and gene and compares them with Rmd 02.
   - "**One-time installation**": `setup_phaser_env.sh` in a compute job, about 20-25 minutes; a pinned conda environment under `{PHASER_HOME}`; never the bioconda package `phaser` (Python 2.7); `PYTHONNOUSERSITE=1`; `module add miniconda3/v4` only in that job.
   - "**Phased or unphased genotypes**": with unphased genotypes the haplotype labels are arbitrary and each gene's counts come from its most-covered haplotype block; with phased genotypes, haplotype A is the left GT allele's haplotype in genome-wide phased genes.
   - "**New requirements**": outbred mode, paired-end data, contigs without `_`, and the same contig names in the GTF, the FASTA and the VCFs.
3. **Prerequisites:**
   - Change the Singularity row note to "Loaded via `module add singularity/3.10.4` for every tool; the only other module is `module add miniconda3/v4`, in the one-time phASER installation job".
   - Add a row "phASER (optional, outbred)": internet on a compute node for `git clone` and conda; about 650 MB under `{PHASER_HOME}`.
4. **Wizard steps table:**
   - Row 7: "Analysis menu: per-sample (always), reciprocal F1, differential ASE (with reference condition), phASER (outbred, paired-end; phased-genotypes and install-location questions)".
   - Add rows `| 18 | phASER: installation and per-sample haplotype counts |` and `| 19 | Rmd 05: phASER gene-level haplotype imbalance |`.
5. **Outputs table:** add rows for:
   - `reference/genes_span.bed`;
   - `scripts/setup_phaser_env.sh`, `phaser_count.sh`, `run_05_phaser.sh`;
   - `phaser/{sample}.haplotypic_counts.txt` and `.gene_ae.txt` (and the other phASER files);
   - `logs/phaser_{sample}.log`;
   - `{TODAY}_{WD_NAME}_05_phaser.html`, `_ASE_phaser.xlsx` (sheets `Gene`, `Comparison`, `Samples`, `Summary`), `_ASE_phaser.pdf`;
   - `ase_phaser_checkpoint.rds` and `summary_numbers_phaser.tsv`.

   Extend the chain sentence with "... and, with phASER, the installation job (only when needed) and the phASER array after the alignment array, then Rmd 05 (`afterok` on the phASER array and Rmd 02)".
6. **Key design points:** add the bullet "**phASER gene test (Rmd 05)**". Summarise planner decisions 2, 3 and 7 and the Step 19 direction rules. State that it is symmetric in A and B, gives no direction with unphased genotypes, includes introns in the gene spans, and offers no differential ASE on phASER counts.
7. **Validation status:** add a "### Stage 3 (synthetic data only)" subsection containing:
   - **Unit tests:** copy the measured sizes, power, the fragment-level comparison and the runtime from `t1_final.log`; cite `tests/r/test_ase_stats_stage3.R`.
   - **Synthetic data:** the Stage 3 checks, 12 defect proofs, the read-level phase check, and determinism.
   - **Installation:** wall time, pinned versions, idempotence, the evidence that `PYTHONNOUSERSITE=1` is needed, and the smoke run without activation.
   - **Rmd 05 on emulated tables:** the Task 4 evaluator results, including the orientation-invariance test.
   - **phASER on aligned synthetic data:** the Task 5 evaluator results in both genotype modes, with the `REPORT` lines and the per-sample phASER wall time.
   - **Negative tests:** contig mismatch, false phased claim, `_` contig, stale output, missing sample.
   - **Chain:** 11 stub dry-run scenarios (`tests/chain/dry_run_chain.sh`).
   - The line `Stage 3 synthetic acceptance run: PENDING`.
   - **What was not exercised:**
     - real data;
     - the `-n 4 --mem=32G -t 2:00:00` phASER request ("not verified on real data");
     - single-end data;
     - real phased WGS genotypes;
     - genes with alternative isoforms or overlapping genes;
     - contig naming with `_`, which is guarded only by a stop;
     - population phasing (`phaser_pop`).
8. **Known limitations:** replace "phASER (Stage 3) is not implemented; the menu lists it as available in a later stage." with a "**Stage 3 (phASER) and its limits**" item. Its sub-bullets:
   - **Maintenance:** phASER is a "Fast Beta" with no releases or tags, pinned by SHA, with open upstream issues (#78 / #85 gene_ae parsing errors on some inputs).
   - **Unphased genotypes:** each gene uses only its most-covered block (on the synthetic test, 113 of 153 SNPs in one sample), and there is no direction.
   - **Labels across individuals:** haplotype A of one individual is unrelated to haplotype A of another, so across individuals only magnitudes can be compared, and there is no differential ASE on phASER counts.
   - **Gene spans:** they include introns, and overlapping genes share reads.
   - **Dispersion:** conservative (the maximum rule and the `RHO_MIN` floor).
   - **Resources:** the phASER resource request is not verified on real data.
   - **Contigs with `_`:** stop the phASER job.
   - **Installation:** it uses a cluster-specific conda module (`module add miniconda3/v4`, `/home/software/conda/miniconda3/bin/condainit`).
9. **Root `README.md`:** in the `/ase-pipeline` row, replace "; phASER not included" with "; in outbred mode, phASER haplotype counts reported next to the unphased results".

- [ ] **Step 3: Run the checker (expected PASS); commit**

Before committing, the controller checks that every number in the README Stage 3 subsection is present in the cited log or report (as in Stage 2 Task 7).

```bash
bash ase-pipeline/tests/check_skill.sh ase-pipeline/ase-pipeline.md
git status && git diff --stat
git add ase-pipeline/README.md README.md ase-pipeline/tests/check_skill.sh
git commit -m "ase-pipeline: README and registration for Stage 3 (phASER)

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 8: Synthetic acceptance on the cluster (controller-run, not a subagent)

The controller runs this task: it needs live cluster jobs and judgement. The scripts and Rmds come from the branch skill through the committed `assemble_outbred_from_skill.sh`, which cuts code blocks by heading and substitutes placeholders. Unlike Stage 2, the **shipped** `submit_chain.sh` and `wait_chain.sh` run live; no manually adapted chain is used. The installed copy `~/.claude/commands/ase-pipeline.md` is not overwritten without the user's approval.

**Files:** Modify `ase-pipeline/README.md` (the acceptance line and results) and `ase-pipeline/tests/check_skill.sh` (tighten the status line to DONE).

- [ ] **Step 1: Data and projects**

1. Re-check the data: `check_simulation_stage3.sh /net/bmc-lab3/data/bcc/yannvrb/ase_synthetic_stage3` must print PASS.
2. Create two projects under `/net/bmc-lab3/data/bcc/yannvrb/ase_accept_stage3/`, with Task 5 Step 6's `mkproj` values, adapted:
   - `obp_unph`: `PHASED_GT` 0, `ind*.all.vcf`, `PHASER_HOME` = `$SCR/phaser_home` (installed, so the chain skips the installation);
   - `obp_ph`: `PHASED_GT` 1, `ind*.all.phased.vcf`, `PHASER_HOME` = `/net/bmc-lab3/data/bcc/yannvrb/ase_accept_stage3/phaser_home` (fresh, so the chain runs the installation job live).
3. Use `ANALYSES` `per-sample phaser` and `TODAY` the day of the run.
4. Run Step 2's `check_r_packages.sh` job once.

- [ ] **Step 2: Submit the shipped chain in each project**

Run `bash {RESULTS_DIR}/scripts/submit_chain.sh` from each project's `CWD`, then `bash {RESULTS_DIR}/scripts/wait_chain.sh` in the background.

Check:
- the dependencies in `chain_job_ids.tsv` and with `scontrol show job` (`sacct` is unreachable);
- `obp_unph` has no `setup_phaser_env.sh` line; `obp_ph` has one.

Record the wall time of the installation job and of each phASER array task.

- [ ] **Step 3: Evaluate against truth**

Run one `sbatch` job with:
- `evaluate_stage3.R phaser` and `evaluate_stage3.R rmd` on `obp_unph` (unphased);
- the same two on `obp_ph` (phased).

Acceptance criteria, all stated with numbers in the report:
- **phASER, per sample:** all the Task 5 `phaser` gates (multi-SNP blocks ≥ 20, phase accuracy ≥ 0.99, trans blocks present, no gene count above the simulated fragments, two-block behaviour by genotype mode, phased-mode orientation).
- **Rmd 05:** all the Task 4 gates (`hap_strong` 24/24, `two_block` 12/12, `hap_lowdepth` ≥ 15/18 with phased genotypes and every one "phASER only (no SNP tested unphased)", null cells ≤ 6, phased orientation, unphased no direction, summary consistency).
- **Rmd 01 and Rmd 02** build without hand edits, and every `ase_counts/*.table` is non-empty.
- **Reported, not gated:**
  - `hap_moderate` detection, phASER against unphased ACAT;
  - `hap_lowdepth` SNPs used and detection with unphased genotypes (single block);
  - aligned against emulated gene counts (the ratio `totalCount` / truth cover per class);
  - the phASER and installation wall times.

- [ ] **Step 4: Summary pages**

Write both `{WD_NAME}_summary_report.html` pages as Step 15 instructs (awk from the TSVs, no R), including the phASER section. The href loop and the no-`http` check must print nothing.

- [ ] **Step 5: Defects, README, ledger**

Defects go through the process: one fix subagent per defect wave, a scoped re-review, and a re-run of the affected steps. When every criterion passes:

1. Replace `Stage 3 synthetic acceptance run: PENDING` with `Stage 3 synthetic acceptance run: DONE (<date>)` followed by the numbers. State the honest limits:
   - the synthetic genome is tiny;
   - wall times say nothing about real data;
   - the phASER resource request is not verified on real data.
2. In `check_skill.sh`, change the Stage 3 status check to `grep -qE "Stage 3 synthetic acceptance run: DONE \(2026-[0-9-]+\)" "$RD" && ! grep -qF "Stage 3 synthetic acceptance run: PENDING" "$RD"`.
3. Run the checker and the dry run.
4. Commit (`git status && git diff --stat` first) with the message `ase-pipeline: record Stage 3 synthetic acceptance results` plus the Co-Authored-By trailer.
5. Update the ledger `.superpowers/sdd/2026-10-01-ase-pipeline-stage3/progress.md`.

- [ ] **Step 6: Whole-branch review**

Run a whole-branch final review (fresh reviewer, base `de43444`). It covers the skill diff, the tests, the README and the acceptance evidence. Findings go through one fix wave, a scoped re-review, the checker, the dry run and, if code in the stats block changed, the unit tests. Install to `~/.claude/commands/` and merge only with the user's approval. Never push.

---

### Task 9 (optional): real-data smoke test, only if the user supplies or approves a dataset

The scope, the dataset and whether this task runs at all are a **user decision** taken at plan review. Tasks 1-8 do not depend on it. If it is not run, the README keeps "not verified on real data" for the phASER request and states that nothing ran on real data.

Example candidate, for the user to accept or replace: one public paired-end human RNA-seq sample with matching phased genotypes (for example a GEUVADIS lymphoblastoid sample and its 1000 Genomes phased VCF), restricted to chr22 and an Ensembl GRCh38 chr22 FASTA/GTF.

- [ ] **Step 1: Agree the dataset with the user**

Agree: the source URLs (verified with a HEAD request in the session, never typed from memory), the licence and data-use terms, the disk budget under `$SCR/real`, and that downloads run in a compute job.

- [ ] **Step 2: Prepare inputs in compute jobs**

- Download the data.
- Subset the reads to chr22 if the source is a BAM: `samtools view` by region, then `samtools fastq`, in the samtools container.
- Subset the VCF to the sample and chr22: `bcftools view -s -r` in the bcftools container. Keep the GT phasing, and rename contigs to the FASTA naming if needed (Step 6 rule).

- [ ] **Step 3: Run the assembled chain**

Assemble a one-sample outbred project with `assemble_outbred_from_skill.sh`, using `PHASED_GT` 1, `PHASER_HOME` `$SCR/phaser_home` and the large-genome-free tier, because chr22 is under 100 Mb. Then run `submit_chain.sh` and `wait_chain.sh`.

Record:
- the phASER task's wall time (`scontrol show job`; `/usr/bin/time` is absent);
- phASER's own "Global maximum memory usage" log line;
- the number of het sites, phased blocks, genes with counts and genome-wide phased genes;
- the Rmd 05 summary row.

- [ ] **Step 4: Update the resource text from the measurement**

With the controller's ruling, update the Step 9 `{PHASER_RESOURCES}` tier and the README resource sentence. Scale chr22 to the genome by its share of het sites, and state the scaling. Commit with the message `ase-pipeline: phASER real-data smoke test (chr22) and resource request` plus the trailer.

---

## Self-review notes

- **Spec coverage:**
  - "Stage 3 — phASER": the installation helper (Task 3), the per-sample run on the BAM and genotype VCF (Task 5), the gene-level step with a gene BED (Task 5), and the gate test run plus "the stage is kept only if the environment installs" (the gate passed; Tasks 3 and 8 re-prove it on the skill's own script).
  - Rmd 02 "outbred with phASER": re-specified as the separate Rmd 05 next to Rmd 02 (planner decision 1; spec amended in Task 1 Step 9).
  - Wizard Step 7 "phASER only in outbred mode" (Task 6).
  - The spec's "the phASER stage creates its own environment" in Step 2 (Task 6 Step 2).
  - Testing: the checker flag rule against recorded `--help` (Tasks 3 and 5), synthetic ground truth (Task 2), and a live acceptance run (Task 8).
  - The real-data smoke test: optional Task 9 (user decision).
  - Staging "Stage 3: phASER, gated": the gate passed (2026-10-01), and the plan carries the gate's constraints as proven facts.
- **Lessons carried from Stages 1 and 2:**
  - Statistics come before any Rmd text, with fixed gates: null size over 10 seeds (pooled and worst), power against an oracle, few genes, heterogeneous dispersion, fragment-level sharing, symmetry, floor and runtime.
  - The H0-inflated dispersion is avoided by the trimmed fit; Neyman-Scott shrinkage is not an issue, because no free mean is fitted.
  - The data exercise every new branch: cis and trans, low depth, two blocks, phased and unphased, and label flips.
  - Every guard has an injected-defect proof.
  - The acceptance run uses the shipped chain live, not a hand-adapted one (Stage 2 M8).
  - The README states honest limits.
- **Interface consistency:**
  - `hap_gene_test(a, b, sample, rho_min, min_genes = 20)` → `list(p, rho_used, samples[sample, n_genes, rho_own, central_frac, rho_cohort, rho_used])` is the same in Tasks 1, 4 and 8.
  - The Rmd 05 checkpoint columns of Task 4 are those that `evaluate_stage3.R` reads (`gene`: `aCount, bCount, totalCount, gw_phased, p, padj, sig, hap_A_frac, major_frac, direction`; `comparison`: `category, sig_unphased`).
  - The `summary_numbers_phaser.tsv` columns are identical in Task 4 (Rmd), Task 4 (evaluator) and Task 6 (summary page).
  - The script names `setup_phaser_env.sh`, `phaser_count.sh` and `run_05_phaser.sh` are identical in Tasks 3-6, the assembler and the dry run.
  - The anchors used by `cut_block.sh` are the headings written in Tasks 3, 5 and 6: `` ### `setup_phaser_env.sh` ``, `` ### `phaser_count.sh` ``, `### Gene spans for phASER`, `### Submission order (one block`, `**Waiting.**`, and the existing Step 10/11/12 headings.
  - The pinned commit string appears in `setup_phaser_env.sh`, `phaser_count.sh` and `submit_chain.sh`, as the checker and the bump procedure require.
- **Review Focus:** each of the 5 items has a test in its owning task: contigs (Task 5 Step 7a, plus the prep guard of Step 2); phased claim (Task 5 Step 7b, Task 4 Step 6); label flips (Task 1 section 1, Task 4 Step 5); zero and low rows (Task 1 section 1, Task 4 Step 5 evaluator); missing or stale outputs (Task 5 Step 7d, Task 4 Step 6).
- **Known risk left to execution:**
  - The R, bash and awk code in this plan has not been run. The binding parts are the interfaces, the tests and the gates. Coding defects are fixed in place; statistical or acceptance gate failures go to the controller as BLOCKED, never to a weaker gate.
  - Three facts are first established during execution, and each has a stated fallback or BLOCKED path:
    - the pinned-version solve (Task 3 Step 4);
    - the run without activation (Task 3 Step 5);
    - STAR and ASEReadCounter with phased GT (Task 5 Step 6).
