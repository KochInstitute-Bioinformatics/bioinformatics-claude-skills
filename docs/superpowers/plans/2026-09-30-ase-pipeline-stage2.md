# ase-pipeline Stage 2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add Stage 2 to the shipped `ase-pipeline` skill: Rmd 03 (reciprocal F1: strain effect and parent-of-origin effect per gene) and Rmd 04 (differential ASE between conditions, F1 and paired outbred). The statistics come first and are shown correct by simulation. Then come a synthetic reciprocal/multi-condition ground-truth dataset, the two Rmd templates, the wizard, run-script and summary integration, the README, and a live cluster acceptance run.

**Architecture:** Every new model lives in the single shipped statistics block of Step 14 (between `# --- ase-stats-begin` and `# --- ase-stats-end`). A new R test script extracts that block from the skill text and tests it by simulation, so the tested code is the shipped code. F1 gene-level models use one row per sample and gene. Each row sums the counts of SNPs thinned so that no read pair is counted twice. The dispersion is estimated from the spread between replicate animals, never from SNP-to-SNP spread. Outbred differential ASE is a per-SNP paired test within each individual, summed over individuals, with no direction. The two Rmds are self-contained templates, rendered by their own `sbatch` scripts and chained after Rmd 02 with `afterok`.

**Tech Stack:** Markdown skill file; bash/awk (no python on this host); R in `/net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif` (base R `optim`/`uniroot`, `parallel`, `aod` only as a test cross-check, `openxlsx`, `tidyverse`, `GenomicRanges`, `rtracklayer`, `Biostrings`); SLURM (`sbatch -p bcc`); the Stage 1 upstream (STAR 2.7.10b, GATK 4.4.0.0 ASEReadCounter, bcftools 1.20, samtools 1.21, Picard 3.1.1 containers).

**Spec:** `docs/superpowers/specs/2026-09-29-ase-pipeline-design.md` (binding; sections "Rmd 03", "Rmd 04", "Staging", "Testing and acceptance"). Task 1 Step 9 amends the spec where the planner decisions below depart from it. Read with the Stage 1 plan `docs/superpowers/plans/2026-09-29-ase-pipeline-stage1.md` and the ledger `.superpowers/sdd/2026-09-29-ase-pipeline-stage1/progress.md`.

## Planner decisions (need user review)

The spec is silent (or, for 3 and 4, too loose to implement correctly) on these points. Each can be reversed independently before execution.

1. **Gene-level F1 counts:** one row per sample and gene. A row sums only the SNPs kept after thinning to one SNP per `THIN_BP` = 500 exon-coordinate bases (the deepest SNP is kept). The dispersion is estimated from the spread between replicate animals, never from SNP-to-SNP spread. Reason: the animal is the unit of biological replication, and thinning stops a read pair from being counted at several SNPs. Stage 1 proved that summed counts combined with a per-SNP rho are wrong. The fragment-level simulation of Task 1 (section 8) must show that this gives the nominal size.
2. **Dispersion:** one pooled moment (Williams-type) estimate per analysis. Each gene is tested with `phi_used = max(phi_common, RHO_MIN, phi_gene)`, where `phi_gene` is used only when the gene has at least 2 residual df. Reason: moment estimates divide by the residual df, so they carry no Neyman-Scott bias. Taking the maximum protects genes that are noisier than average.
3. **Beta-binomial GLMs are fitted by base-R code** in the stats block (bounded L-BFGS-B, fixed dispersion, likelihood-ratio tests), not by `aod::betabin`. `aod` is used only in the unit tests, to cross-check the likelihood. Reason: aod estimates a free per-gene phi from 4-12 rows, which is biased low and so anti-conservative, and its unbounded fit misbehaves on monoallelic (imprinted) genes.
4. **`lme4::glmer` is not used.** Outbred differential ASE is tested per SNP and individual with a paired beta-binomial LRT (one REF fraction versus one per condition), using a per-individual dispersion. The statistics are summed over individuals (df = number of individuals), the gene level uses ACAT, and no direction is reported. Reason: the allele that carries a regulatory variant differs between individuals, so a shared condition slope (fixed or random) cancels real changes. A random-effect variance is also not estimable with 2-6 individuals.
5. **Sign conventions:**
   - p = the strain-A (REF) fraction, as in the spec. Note that Rmd 02 reports the strain-B fraction.
   - d = +1 for `AxB` (strain A is the mother).
   - b0 > 0 is reported as "`{STRAIN_A}` higher"; b1 > 0 as "maternal higher".
   - In Rmd 04, delta > 0 means that the strain-A fraction (F1) or the REF fraction (outbred) is higher in the tested condition than in the reference condition.
6. **`cross_direction` must be exactly `AxB` or `BxA`.** A is `{STRAIN_A}` and the mother is written first. Any other value stops Rmds 03 and 04. Reason: a free-text value such as `B6xAJ` could silently invert every parent-of-origin call.
7. **Replicate minimums:**
   - Reciprocal: at least 2 samples per direction, and per gene at least 2 samples with coverage in each direction.
   - F1 differential: at least 2 samples per condition.
   - Outbred differential: at least 2 individuals sampled in both conditions, and per SNP at least 2 informative individuals.

   Genes or SNPs below these minimums are listed as not tested.
8. **Single-SNP genes and overlapping genes:** genes with one SNP are tested at gene level, because the replication comes from samples, not SNPs. SNPs that overlap two genes are left out of the gene-level tests.
9. **Covariates:**
   - Rmd 03 adjusts for condition (sum-to-zero coding) when several conditions exist.
   - Rmd 04 (F1) adds the cross-direction term when both directions are among the contrast's samples.
   - Rmd 04 stops with a message when condition and direction are confounded.
10. **Contrasts:** each other condition is compared with one user-chosen reference condition, by a 1-df LRT per contrast, with BH within each contrast and level (SNP or gene). There is no omnibus test.
11. **Excluded chromosomes:** X, Y and MT (and their `chr`-prefixed forms) are excluded from the Rmd 03/04 tests and counted in a separate table. Reason: a hemizygous X in males and the maternally inherited MT mimic maternal imprinting.
12. **Effect-size threshold:** `ABS_DEV_SIG` is reused as the effect-size threshold: `|frac_A - 0.5|` (strain), `|maternal_frac - 0.5|` (parent of origin) and `|delta fraction|` (differential).
13. **Rmd 02 thinning:** the F1 gene test of Rmd 02 also moves to thinned SNPs (Task 5), and does so only if Task 1 measures the unthinned Stage 1 gene test as anti-conservative. This closes the Stage 1 limitation I6 but changes Stage 1 results.

## Proven facts (verified in Stage 1: do not re-derive)

- **ASEReadCounter** needs all of the following:
  - a read group in the BAM (STAR `--outSAMattrRGline`);
  - a bgzipped, tabix-indexed VCF with a heterozygous genotype column;
  - `-R` with `.fai` and `.dict`.

  Without them it writes an empty table and exits 0. Its columns are `contig position variantID refAllele altAllele refCount altCount totalCount lowMAPQDepth lowBaseQDepth rawDepth otherBases improperPairs`. It counts a fragment at every SNP the fragment overlaps, so nearby SNPs share read pairs.
- **STAR WASP** does not ignore homozygous genotypes (use a het-only VCF) and reads only the first sample column.
- **Packages in the image:** `aod` 1.3.3, `lme4` 2.0.6, `openxlsx`, `tidyverse`, `GenomicRanges`, `rtracklayer`, `Biostrings` and `parallel` are present. `VGAM`, `glmmTMB`, `MBASED` and `VariantAnnotation` are absent. R runs only from that image, and only through `sbatch -p bcc`.
- **Dispersion estimates from Stage 1:**
  - A rho estimated under H0 (p = 0.5) from all sites is inflated by real imbalance: 0.08-0.10 against a binomial truth, which detected 0 of 16 planted outbred genes.
  - A maximum-likelihood rho with a free mean per group is biased low by (k - 1)/k (Neyman-Scott: 0.0115 for a true 0.02 with 4-SNP genes). `bb_estimate_rho_gene` corrects it with `cf * rho_free + (cf - 1) / m`.
- **`RHO_MIN` = 0.01 floor:** with rho 0 on the acceptance data there were 3 of 28 null false positives. The floor also caps outbred power at about 12 of 16.
- **ACAT:** a p-value of exactly 1 makes ACAT collapse to 1 (`tan(-pi/2)`), so p is capped at 0.99. `bb_pvalue` returns exactly 1 at the modal count.
- **Unit-test design:**
  - Single-seed size tests are fragile. Use 10 seeds with a pooled gate and a worst-seed gate.
  - An "exactly 0 on binomial data" assertion must be a majority rule (ML sits at the boundary only about half the time).
  - The unit-test job needs at least 8 CPUs (`mc.cores = 8`).
- **Runtime at scale** (1 CPU): `bb_estimate_rho_gene` takes 3.5 s for 30,000 SNPs in 6,000 genes; `bb_estimate_rho_trim` takes 26 s and 0.7 GB for 50,000 sites. `-n 1 --mem=16G` is adequate for Rmd 02.
- **bcftools and Rmd 01 fixes:**
  - `bcftools consensus -H A` with a sites-only VCF applied the sites nondeterministically (4 of 10 runs), so the skill uses a genotyped `MASK` sample and `-s MASK`.
  - Rmd 01 reads contigs with `colClasses` character; otherwise `1`..`22` versus `X` fails to bind.
- **SLURM job state:** `sacct` is unreachable on this cluster ("Connection refused"), so job states come from `scontrol show job`.
- **Stage 1 synthetic data:**
  - The generator's SYNTHETIC genome mode (chr1, plus-strand genes) is the base for synthetic data, because the nf-core test genome has too few loci.
  - Stage 1 synthetic data: `/net/bmc-lab3/data/bcc/yannvrb/ase_synthetic/` (seed 20260929).
  - Stage 1 acceptance outputs: `/net/bmc-lab3/data/bcc/yannvrb/ase_accept2/{f1,outbred}/results/2026-09-29_{f1,outbred}/` (F1: 24 of 24 planted gene x sample calls, 0 of 42 null; outbred: 10 of 16, 0 of 28).
- **"Wizard-by-script"** means the scripts and Rmds are generated from the skill's code blocks by extraction plus placeholder substitution. The interactive dialogue of Steps 0-9 is not driven.
- **The Stage 1 acceptance run found defects the unit tests missed:** D1, the outbred rho that detected 0 of 16, and I1-I6, input paths never exercised. The Stage 2 synthetic data therefore exercises every new branch: both directions, both conditions, imprinted, strain and differential genes, monoallelic genes, SNP-dense genes, and phase-heterogeneous outbred changes.
- STAR prints "Could not move Log.out" for the index build when `--outFileNamePrefix` is set. It is harmless.

## Global Constraints

- **Skill file:** `ase-pipeline/ase-pipeline.md` is self-contained (no references to other skills' steps). Every finite-choice question is a numbered list; open questions are used only for email and paths.
- **Paths and inputs:**
  - Results go to `results/{TODAY}_{WD_NAME}` with `{TODAY}` = `YYYY-MM-DD` (user data-safety rule).
  - Raw FASTQ, BAM and VCF inputs are read-only.
  - Paths with spaces or commas are refused (Step 0).
- **Cluster use:**
  - Never run R, Singularity or any heavy work on the login node. Every job is `sbatch -p bcc`, and debugging uses `salloc`/`srun`.
  - Default requests are at most 64 G and 4 h.
  - Resources must be honest: `-n 1` unless the code parallelises.
- **Containers and modules:** tools come from Singularity containers; the only module is `module add singularity/3.10.4`, checked with `|| exit 1`. `BIND` lists each directory once. `/tmp` is node-local, so job outputs and temporary files go to the shared filesystem.
- **Generated scripts:** use `set -uo pipefail` with explicit exit-code checks (no blind `set -e`). Chained jobs use `--parsable` and `--dependency=afterok`.
- **Flags:** every `--flag` in the skill text appears in the recorded fixtures (`star_help.txt`, `gatk_ASEReadCounter_help.txt`) or in the checker allowlist (`bind mem mail-type mail-user array dependency parsable regions samples genotype types min-alleles max-alleles rename-chrs`). Never add a STAR or GATK flag to the allowlist.
- **Rmd conventions:**
  - Rmds are self-contained (no `source()`) and use `knitr::opts_chunk$set(cache = FALSE, ...)` and `options(scipen = 9)`.
  - Bioconductor packages load before `tidyverse`, and verbs are `dplyr::`-prefixed.
  - Shared constants are defined once per Rmd, and each Rmd checks them against the Rmd 01 checkpoint.
  - Tables and figures use only the constants and one `sig` column. Each Rmd checks that the counts drawn in its figures equal its Summary counts.
- **Statistics code:**
  - It lives only in the one block between `# --- ase-stats-begin` and `# --- ase-stats-end` (exactly one pair of markers).
  - That block is pasted verbatim into Rmds 02, 03 and 04 at the line `# <<< paste here the code of Step 14 between the two marker lines (the marker lines themselves are not pasted) >>>`.
  - The unit tests extract it from the skill text.
- **Simulation gates are binding:** the thresholds are fixed in this plan. If a gate fails, the implementer reports BLOCKED to the controller with the measured numbers and never weakens a gate, changes a threshold or swaps the method.
- **Checker lines:** every new `need`/`forbid` line must be shown to FAIL on the pre-task skill (`git show HEAD:ase-pipeline/ase-pipeline.md > $SCR/pre.md`, then run the checker on `$SCR/pre.md`).
- **Scratch:** `SCR=/net/bmc-lab3/data/bcc/yannvrb/ase_stage2_scratch` holds everything that is not committed: synthetic data, renders and logs. Never write generated data into the repository.
- **Git:**
  - Show `git status` and `git diff --stat` before each commit. Never `git push` (the user must approve).
  - Commit messages end with a blank line and `Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>`.
  - Branch: `feat/ase-pipeline-stage2`.
- **Stage 3 (phASER) is out of scope.** The wizard keeps saying "available in a later stage" for phASER.

## Review Focus

1. **Cross-direction coding.** `AxB` must give d = +1. A swapped or free-text value (`B6xAJ`) must stop the Rmd, never produce inverted parent-of-origin calls. Tests: Task 1 section 5 (sign of b0 and b1 in simulation), and Task 3 Step 6 (a free-text value stops Rmd 03; directions on synthetic truth are correct).
2. **Condition confounded with cross direction** (every ctrl sample is `AxB` and every treat sample `BxA`). This must stop the F1 differential with a clear message rather than report parent of origin as condition. Tests: Task 1 section 6 (confounded design is `not_estimable`) and Task 4 Step 6 (a confounded sheet makes Rmd 04 stop with the message).
3. **Monoallelic genes** (fully imprinted, or X-linked in males). These must give a finite p and status `ok_at_bound`, never an error. X, Y and MT are excluded and counted. Tests: Task 1 section 3 (complete separation) and Task 3 Step 6 (a gene moved to `chrX` is excluded and listed).
4. **Unpaired outbred individuals,** and SNPs heterozygous in only one individual. These are left out of the paired test and reported, never tested as if paired. Tests: Task 1 section 7 (driver edge cases) and Task 4 Step 6 (an extra ctrl-only individual is listed as unpaired, and the other results are unchanged).
5. **SNP-dense genes whose SNPs share read pairs.** Thinning keeps SNPs at least `THIN_BP` apart, and the gene test keeps its size. Tests: Task 1 section 8 (fragment-level simulation) and Task 3 Step 5 (the kept SNPs of the `dense_null` genes are at least `THIN_BP` apart in `SNPs_used`).

## File Structure

| File | Responsibility | Task |
|---|---|---|
| `ase-pipeline/ase-pipeline.md` | Step 14 stats block and prose, Step 8 `THIN_BP`; new Step 16 (Rmd 03) and Step 17 (Rmd 04); Steps 5, 7, 15 and Notes; Step 13 (Rmd 02 thinning, conditional) | 1, 3, 4, 5, 6 |
| `ase-pipeline/tests/r/test_ase_stats_stage2.R` | Simulation unit tests of the Stage 2 functions (extracted from the skill) | 1 |
| `ase-pipeline/tests/r/run_stats_tests.sh` | `sbatch` wrapper that runs the Stage 1 and Stage 2 unit tests in the container | 1 |
| `ase-pipeline/tests/synthetic/simulate_ase_stage2.R` | Stage 2 ground truth: reciprocal x two-condition F1 and paired two-condition outbred; FASTQs, VCFs, truth tables, direct count tables | 2 |
| `ase-pipeline/tests/synthetic/check_simulation_stage2.sh` | Exact-quota and consistency checks of the Stage 2 dataset | 2 |
| `ase-pipeline/tests/synthetic/check_biology_stage2.R` | FASTQ-level allele-content check against the realised truth | 2 |
| `ase-pipeline/tests/synthetic/prove_checks_stage2.sh` | Injected-defect proofs: each defect must make the check FAIL | 2 |
| `ase-pipeline/tests/synthetic/render_from_skill.sh` | Cuts a Step's `rmd` block out of the skill, splices the stats block, substitutes placeholders | 3 |
| `ase-pipeline/tests/synthetic/evaluate_stage2.R` | Evaluates Rmd 03/04 outputs against the truth tables (used by Tasks 3, 4, 8) | 3, 4 |
| `ase-pipeline/tests/check_skill.sh` | Static checker, one block per task | 1, 3-7 |
| `ase-pipeline/README.md`, root `README.md` | Stage 2 documentation, validation status | 7, 8 |
| `docs/superpowers/specs/2026-09-29-ase-pipeline-design.md` | Spec amendments for the planner decisions | 1 |

Step numbering is fixed now to avoid the Stage 1 renumbering conflict. The existing Steps 0-15 keep their numbers. Rmd 03 becomes `## Step 16` and Rmd 04 `## Step 17`, both placed after Step 15 and before `## Notes for the assistant`. Step 15 (submission and summary) refers forward to them.

All repo paths are relative to `/net/bmc-lab3/data/bcc/yannvrb/Github_repos/bioinformatics-claude-skills`.

---

### Task 1: Stage 2 statistics functions with simulation unit tests (plus spec amendments)

**Files:**
- Modify: `ase-pipeline/ase-pipeline.md`:
  - Step 14 stats block: insert the new functions directly above the line `# --- ase-stats-end` (currently :1327).
  - Step 14 intro sentence (:1230) and prose bullets (after :1337).
  - Step 8 constants table: new row after the `RHO_MIN` row (:206).
- Create: `ase-pipeline/tests/r/test_ase_stats_stage2.R`, `ase-pipeline/tests/r/run_stats_tests.sh`
- Modify: `ase-pipeline/tests/check_skill.sh` (new block before the final `[ $fail -eq 0 ]` line)
- Modify: `docs/superpowers/specs/2026-09-29-ase-pipeline-design.md` (Rmd 03, Rmd 04 and GLM paragraphs)

**Interfaces:**
- Consumes: the existing block functions `bb_ll`, `bb_pvalue`, `acat`, `bb_gene_lrt`, `bb_estimate_rho_gene` (unchanged), and the Step 8 row `| \`RHO_MIN\` | 0.01 |`.
- Produces (every later task uses these exact names):
  - `chrom_class(contig)` returns a character vector over `"autosome"`, `"X"`, `"Y"` and `"MT"`.
  - `thin_snps(pos, depth, window)` returns a logical keep vector.
  - `bb_ll_rows(x, n, p, rho)` returns the per-row log-likelihood.
  - `bb_glm_fit(y, n, X, rho, bound = 15)` returns `list(beta` (named by `colnames(X)`)`, loglik, conv, at_bound)`.
  - `bb_glm_lrt(y, n, X, drop_cols, rho, bound = 15)` returns `list(beta, stat, df, p, conv, at_bound)`.
  - `bb_moment_phi(units, max_iter = 25, tol = 1e-6, bound = 15)` returns a numeric scalar (`NA` if no unit has residual df). `units` is a list of `list(y, n, X)`.
  - `ase_glm_test(units, tests, rho_min, min_df_unit = 2, bound = 15)` returns a data.frame. `units` is a named list of `list(y, n, X)` with named columns and rows with n > 0; `tests` is a named list mapping a test name to the columns dropped under H0. The columns are exactly `unit, n_rows, df_resid, depth, phi_common, phi_unit, phi_used, status`, then `beta_<column>` for every X column name, then `stat_<test>` and `p_<test>` for every test. `status` is one of `ok`, `ok_at_bound`, `not_converged`, `fit_error`, `not_estimable`.
  - `ase_paired_test(d, rho_min, min_individuals = 2, bound = 15)` returns `list(snp, phi)`. `d` is a data.frame with `snp, individual, cond, y, n` (cond 0 = reference, 1 = tested; y = REF count). `snp` has columns exactly `snp, n_individuals, n_failed, stat, df, p, mean_delta, max_abs_delta, n_up, n_down, status` (status `ok` or `too_few_individuals`). `phi` has columns `individual, phi_pair`.
  - Step 8 constant row `| \`THIN_BP\` | 500 | ... |`.
  - `sbatch -p bcc -o <log> ase-pipeline/tests/r/run_stats_tests.sh [stage1|stage2|all]`.

- [ ] **Step 1: Probe the `aod` object layout (tiny job, records a fact the tests rely on)**

```bash
SCR=/net/bmc-lab3/data/bcc/yannvrb/ase_stage2_scratch; mkdir -p $SCR/logs
cat > $SCR/probe_aod.sh <<'EOF'
#!/bin/bash
#SBATCH -p bcc
#SBATCH -n 1 --mem=4G -t 0:10:00
module add singularity/3.10.4 || exit 1
singularity exec --bind /net/bmc-lab3 /net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif Rscript -e '
suppressPackageStartupMessages(library(aod)); set.seed(1); d <- rep(c(1, -1), 6); n <- rep(80, 12); y <- rbinom(12, n, plogis(0.4 + 0.8 * d))
f <- aod::betabin(cbind(y, n - y) ~ d, ~ 1, data = data.frame(y, n, d), fixpar = list(3, 0.02))
print(packageVersion("aod")); print(slotNames(f)); print(f@fixed.param); print(f@random.param); print(f@logL)'
EOF
sbatch -p bcc -o $SCR/logs/probe_aod_%j.out $SCR/probe_aod.sh
```
Expected: the version is 1.3.3, the slot names include `fixed.param`, `random.param` and `logL`, and `random.param` equals 0.02 (so the dispersion was fixed). If a slot is named differently, use the printed name in the two `aod` lines of test section 2 below. This is the only adaptation allowed. If `random.param` is not 0.02, `fixpar` did not fix the dispersion: report BLOCKED.

- [ ] **Step 2: Write the failing test script**

Create `ase-pipeline/tests/r/test_ase_stats_stage2.R` with exactly this content:

```r
# Stage 2 statistics tests (reciprocal F1, differential ASE); the functions are extracted from the shipped skill text.
# Usage: Rscript test_ase_stats_stage2.R <ase-pipeline.md>     (needs >= 8 CPUs: parallel::mclapply, mc.cores = 8)
args <- commandArgs(trailingOnly = TRUE); skill <- args[1]
txt <- readLines(skill)
b <- grep("^# --- ase-stats-begin", txt); e <- grep("^# --- ase-stats-end", txt)
stopifnot(length(b) == 1, length(e) == 1, e > b)
eval(parse(text = txt[(b + 1):(e - 1)]))
fail <- 0
ok <- function(cond, msg) { if (!isTRUE(cond)) { cat("FAIL:", msg, "\n"); fail <<- 1 } else cat("ok  ", msg, "\n") }
fns <- c("chrom_class", "thin_snps", "bb_ll_rows", "bb_glm_fit", "bb_glm_lrt", "bb_moment_phi", "ase_glm_test", "ase_paired_test")
miss <- fns[!vapply(fns, exists, logical(1))]
if (length(miss) > 0) { cat("FAIL: missing from the stats block:", paste(miss, collapse = ", "), "\n"); quit(status = 1) }
const <- function(name) {   # a Step 8 default, read from the skill text (never hard-coded here)
  l <- grep(paste0("^\\| `", name, "` \\|"), txt, value = TRUE); stopifnot(length(l) == 1)
  as.numeric(trimws(strsplit(l, "\\|")[[1]][3]))
}
RHO_MIN <- const("RHO_MIN"); THIN_BP <- const("THIN_BP")
ok(RHO_MIN == 0.01 && THIN_BP == 500, sprintf("Step 8 defaults read from the skill: RHO_MIN %g, THIN_BP %g", RHO_MIN, THIN_BP))
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

# ---- 1. helpers
ok(identical(chrom_class(c("1", "chr2", "X", "chrX", "Y", "MT", "chrM", "M", "x")),
             c("autosome", "autosome", "X", "X", "Y", "MT", "MT", "MT", "X")), "chrom_class maps Ensembl and chr-prefixed names")
ok(identical(thin_snps(c(100, 150, 700, 720, 1500), c(10, 50, 30, 30, 5), 500), c(FALSE, TRUE, TRUE, FALSE, TRUE)),
   "thin_snps keeps the deepest SNP per window (ties broken by position)")
ok(identical(thin_snps(42, 7, 500), TRUE), "thin_snps keeps a single SNP")
ok(identical(thin_snps(numeric(0), numeric(0), 500), logical(0)), "thin_snps with no SNP returns logical(0)")
ok(all(thin_snps(c(1, 2000, 4000), c(1, 1, 1), 500)), "SNPs further apart than the window are all kept")
set.seed(11); ps <- sort(sample(1:5000, 60)); kp <- ps[thin_snps(ps, rpois(60, 50), THIN_BP)]
ok(all(diff(kp) >= THIN_BP), "kept SNPs are at least THIN_BP apart")

# ---- 2. GLM fit: equals glm(binomial) at rho = 0 and aod::betabin with the dispersion fixed
set.seed(12); n <- rpois(12, 80) + 5; d <- rep(c(1, -1), 6); X <- cbind("(Intercept)" = 1, d = d)
y <- rbinom(12, n, plogis(0.4 + 0.8 * d))
g <- glm(cbind(y, n - y) ~ d, family = binomial); f0 <- bb_glm_fit(y, n, X, 0)
ok(max(abs(f0$beta - coef(g))) < 1e-4 && abs(f0$loglik - as.numeric(logLik(g))) < 1e-6,
   "bb_glm_fit with rho = 0 equals glm(binomial): coefficients and log-likelihood")
suppressPackageStartupMessages(library(aod))
dd <- data.frame(y = y, n = n, d = d)
a1 <- aod::betabin(cbind(y, n - y) ~ d, ~ 1, data = dd, fixpar = list(3, 0.02))
a0 <- aod::betabin(cbind(y, n - y) ~ 1, ~ 1, data = dd, fixpar = list(2, 0.02))
f2 <- bb_glm_fit(y, n, X, 0.02); f20 <- bb_glm_fit(y, n, X[, 1, drop = FALSE], 0.02)
ok(max(abs(f2$beta - a1@fixed.param)) < 1e-3, "bb_glm_fit at rho 0.02 matches aod::betabin (phi fixed at 0.02): coefficients")
ok(abs((f2$loglik - f20$loglik) - (a1@logL - a0@logL)) < 1e-3, "bb_glm_fit log-likelihood difference matches aod (same LRT)")
l2 <- bb_glm_lrt(y, n, X, "d", 0.02)
ok(abs(l2$stat - 2 * (f2$loglik - f20$loglik)) < 1e-6 && l2$df == 1, "bb_glm_lrt statistic = 2 x log-likelihood difference, df 1")

# ---- 3. separation and driver edge cases
set.seed(13)
mk <- function(y, n, d) list(y = y, n = n, X = cbind("(Intercept)" = 1, d = d))
d6 <- c(1, 1, 1, -1, -1, -1)
base <- lapply(1:200, function(i) { n <- rpois(6, 100) + 10; mk(rbinom(6, n, plogis(rnorm(1, 0, 0.5))), n, d6) })
names(base) <- paste0("g", 1:200)
edge <- list(sep = mk(c(60, 70, 50, 0, 0, 0), c(60, 70, 50, 55, 65, 45), d6),
             onedir = mk(c(30, 32, 28), rep(50, 3), c(1, 1, 1)),
             onerow = mk(30, 50, 1),
             empty = list(y = numeric(0), n = numeric(0), X = cbind("(Intercept)" = numeric(0), d = numeric(0))),
             broken = mk(c(30, 32, 28, 20, 22, 25), c(50, 50, Inf, 50, 50, 50), d6),
             onesnp = mk(c(40, 45, 38, 12, 10, 15), rep(50, 6), d6))
r <- ase_glm_test(c(base, edge), list(strain = "(Intercept)", parent_of_origin = "d"), RHO_MIN)
ok(identical(names(r), c("unit", "n_rows", "df_resid", "depth", "phi_common", "phi_unit", "phi_used", "status", "beta_(Intercept)",
                         "beta_d", "stat_strain", "stat_parent_of_origin", "p_strain", "p_parent_of_origin")), "ase_glm_test output columns")
st <- setNames(r$status, r$unit); rr <- function(u, col) r[[col]][r$unit == u]
ok(st[["sep"]] == "ok_at_bound" && rr("sep", "p_parent_of_origin") < 1e-6 && rr("sep", "beta_d") > 5,
   "monoallelic maternal gene (complete separation): ok_at_bound, p < 1e-6, b1 large and positive, no error")
ok(all(st[c("onedir", "onerow", "empty")] == "not_estimable") && all(is.na(r$p_strain[r$unit %in% c("onedir", "onerow", "empty")])),
   "one direction only / one row / no rows: not_estimable with p NA")
ok(st[["broken"]] == "fit_error" && is.na(rr("broken", "p_strain")), "a unit whose fit fails is fit_error with p NA")
ok(all(st[names(base)] %in% c("ok", "ok_at_bound")) && all(is.finite(r$p_strain[r$unit %in% names(base)])), "a failing unit does not affect the others")
ok(st[["onesnp"]] == "ok" && rr("onesnp", "beta_d") > 0.8 && rr("onesnp", "p_parent_of_origin") < 0.05,
   "a gene with one SNP (one row per sample) is tested; maternal-rich counts give b1 > 0")
ok(all(r$phi_used[r$status %in% c("ok", "ok_at_bound")] >= RHO_MIN), "phi_used is never below RHO_MIN")

# ---- 4. moment estimator: recovers phi, no Neyman-Scott shrinkage with 6 rows and 2 coefficients
sim_units <- function(U, d, phi, depth = function(m) round(runif(m, 30, 500))) {
  u <- lapply(seq_len(U), function(i) {
    p <- plogis(rnorm(1, 0, 0.5) + rnorm(1, 0, 0.5) * d); n <- depth(length(d))
    pr <- if (phi > 0) rbeta(length(d), p * (1 - phi) / phi, (1 - p) * (1 - phi) / phi) else p
    list(y = rbinom(length(d), n, pr), n = n, X = cbind("(Intercept)" = 1, d = d))
  })
  setNames(u, paste0("u", seq_len(U)))
}
for (ph in c(0.005, 0.02, 0.05)) {
  est <- unlist(parallel::mclapply(1:8, function(s) { set.seed(200 + s); bb_moment_phi(sim_units(1000, d6, ph)) }, mc.cores = MC))
  ok(abs(mean(est) / ph - 1) < 0.15, sprintf("bb_moment_phi recovers phi %.3f: mean %.4f over 8 seeds (within 15%%)", ph, mean(est)))
}
est0 <- unlist(parallel::mclapply(1:8, function(s) { set.seed(300 + s); bb_moment_phi(sim_units(1000, d6, 0)) }, mc.cores = MC))
ok(all(est0 < 2e-3), sprintf("binomial data: bb_moment_phi max %.2e (< 2e-3)", max(est0)))

# ---- 5. reciprocal F1 model (rows = samples, logit p_A = b0 + b1 d): size, leakage, power, sign
sim_design <- function(G, d, cond, b0, b1, bc, phi, phi_sd = 0, depth = function(m) round(runif(m, 60, 400))) {
  lapply(seq_len(G), function(g) {
    ph <- if (phi_sd > 0) min(phi * exp(rnorm(1, -phi_sd^2 / 2, phi_sd)), 0.3) else phi
    p <- plogis(b0[g] + b1[g] * d + bc[g] * cond); n <- depth(length(d))
    pr <- if (ph > 0) rbeta(length(d), p * (1 - ph) / ph, (1 - p) * (1 - ph) / ph) else p
    list(y = rbinom(length(d), n, pr), n = n)
  })
}
recip_run <- function(seed, nA, nB, phi, phi_sd = 0) {
  set.seed(seed); d <- c(rep(1, nA), rep(-1, nB))
  cls <- rep(c("null", "strain_only", "imprinted_only", "power_strain", "power_poo"), c(2000, 400, 400, 200, 200)); G <- length(cls)
  b0 <- ifelse(cls == "strain_only", sample(c(-1, 1), G, TRUE) * qlogis(0.8), ifelse(cls == "power_strain", qlogis(0.7), 0))
  b1 <- ifelse(cls == "imprinted_only", sample(c(-1, 1), G, TRUE) * qlogis(0.9), ifelse(cls == "power_poo", qlogis(0.8), 0))
  s <- sim_design(G, d, 0, b0, b1, rep(0, G), phi, phi_sd)
  u <- setNames(lapply(s, function(z) list(y = z$y, n = z$n, X = cbind("(Intercept)" = 1, d = d))), paste0("g", seq_len(G)))
  r <- ase_glm_test(u, list(strain = "(Intercept)", parent_of_origin = "d"), RHO_MIN)
  sz <- function(p, k) mean(p[cls == k] < 0.05, na.rm = TRUE)
  ss <- which(r$p_strain < 0.05 & cls %in% c("strain_only", "power_strain"))           # which(): p = NA never indexes
  sp <- which(r$p_parent_of_origin < 0.05 & cls %in% c("imprinted_only", "power_poo"))
  c(null_strain = sz(r$p_strain, "null"), null_poo = sz(r$p_parent_of_origin, "null"),
    leak_poo = sz(r$p_parent_of_origin, "strain_only"), leak_strain = sz(r$p_strain, "imprinted_only"),
    pow_strain = sz(r$p_strain, "power_strain"), pow_poo = sz(r$p_parent_of_origin, "power_poo"),
    sign_strain = mean(sign(r[["beta_(Intercept)"]][ss]) == sign(b0[ss])), sign_poo = mean(sign(r$beta_d[sp]) == sign(b1[sp])),
    bad_status = mean(!r$status %in% c("ok", "ok_at_bound")), phi_common = r$phi_common[1])
}
rc <- list(list("3+3, phi 0.005", 3, 3, 0.005, 0), list("3+3, phi 0.02", 3, 3, 0.02, 0), list("2+2, phi 0.02", 2, 2, 0.02, 0),
           list("3+2, phi 0.02", 3, 2, 0.02, 0), list("4+2, phi 0.02", 4, 2, 0.02, 0),
           list("3+3, phi 0.02 heterogeneous (lognormal sd 0.7)", 3, 3, 0.02, 0.7))
for (cs in rc) {
  m <- seeds_rbind(1:10, function(s) recip_run(4000 + s, cs[[2]], cs[[3]], cs[[4]], cs[[5]])); report(m, paste("reciprocal", cs[[1]]))
  het <- cs[[5]] > 0; pl <- if (het) 0.075 else 0.065; wo <- if (het) 0.10 else 0.09
  for (col in c("null_strain", "null_poo", "leak_poo", "leak_strain")) gate(m, col, paste("reciprocal", cs[[1]], col), pl, wo)
  ok(min(m[, "sign_strain"]) >= 0.99 && min(m[, "sign_poo"]) >= 0.99, sprintf("reciprocal %s: sign of b0 and b1 correct in >= 99%% of calls", cs[[1]]))
  ok(max(m[, "bad_status"]) <= 0.01, sprintf("reciprocal %s: at most 1%% of genes not ok (max %.4f)", cs[[1]], max(m[, "bad_status"])))
  if (cs[[1]] == "3+3, phi 0.005") ok(mean(m[, "pow_strain"]) >= 0.9 && mean(m[, "pow_poo"]) >= 0.9,
    sprintf("reciprocal 3+3 power: strain 0.7 %.3f, maternal 0.8 %.3f (>= 0.9 each)", mean(m[, "pow_strain"]), mean(m[, "pow_poo"])))
}

# ---- 6. F1 differential (logit p_A = b0 + bc cond [+ b1 d]): size with strain / imprinted genes present, power, confounding
diff_run <- function(seed, dir_ctrl, dir_treat, phi, with_d = TRUE) {
  set.seed(seed); d <- c(dir_ctrl, dir_treat); cond <- rep(c(0, 1), c(length(dir_ctrl), length(dir_treat)))
  cls <- rep(c("null", "imprinted", "power"), c(2000, 600, 300)); G <- length(cls)
  b0 <- rnorm(G, 0, 0.6)
  b1 <- ifelse(cls == "imprinted", sample(c(-1, 1), G, TRUE) * qlogis(0.9), 0)
  bc <- ifelse(cls == "power", sample(c(-1, 1), G, TRUE) * qlogis(0.7), 0)
  s <- sim_design(G, d, cond, b0, b1, bc, phi)
  X <- if (with_d && length(unique(d)) == 2) cbind("(Intercept)" = 1, cond = cond, d = d) else cbind("(Intercept)" = 1, cond = cond)
  u <- setNames(lapply(s, function(z) list(y = z$y, n = z$n, X = X)), paste0("g", seq_len(G)))
  r <- ase_glm_test(u, list(condition = "cond"), RHO_MIN); sg <- r$p_condition < 0.05; pw <- which(sg & cls == "power")
  c(null = mean(sg[cls == "null"], na.rm = TRUE), imprinted = mean(sg[cls == "imprinted"], na.rm = TRUE),
    power = mean(sg[cls == "power"], na.rm = TRUE), sign = mean(sign(r$beta_cond[pw]) == sign(bc[pw])),
    bad_status = mean(!r$status %in% c("ok", "ok_at_bound")))
}
dc <- list(list("3 vs 3, one direction", c(1, 1, 1), c(1, 1, 1), TRUE), list("2 vs 2, one direction", c(1, 1), c(1, 1), TRUE),
           list("3 vs 2, one direction", c(1, 1, 1), c(1, 1), TRUE),
           list("unbalanced directions, with d term", c(1, 1, -1), c(1, -1, -1), TRUE))
for (cs in dc) {
  m <- seeds_rbind(1:10, function(s) diff_run(5000 + s, cs[[2]], cs[[3]], 0.02, cs[[4]])); report(m, paste("differential F1", cs[[1]]))
  gate(m, "null", paste("differential F1", cs[[1]], "null genes")); gate(m, "imprinted", paste("differential F1", cs[[1]], "imprinted genes (no condition effect)"))
  ok(min(m[, "sign"]) >= 0.99 && max(m[, "bad_status"]) <= 0.01, sprintf("differential F1 %s: sign correct, <= 1%% not ok", cs[[1]]))
  if (cs[[1]] == "3 vs 3, one direction") ok(mean(m[, "power"]) >= 0.7, sprintf("differential F1 3 vs 3 power (0.5 -> 0.7): %.3f (>= 0.7)", mean(m[, "power"])))
}
m <- seeds_rbind(1:10, function(s) diff_run(5000 + s, c(1, 1, -1), c(1, -1, -1), 0.02, FALSE))
cat(sprintf("     EVIDENCE unbalanced directions WITHOUT the d term: imprinted-gene size %.4f (with the d term: gated above)\n", mean(m[, "imprinted"])))
set.seed(5100); s <- sim_design(50, c(1, 1, 1, -1, -1, -1), c(0, 0, 0, 1, 1, 1), rep(0, 50), rep(0, 50), rep(0, 50), 0.02)
uc <- setNames(lapply(s, function(z) list(y = z$y, n = z$n, X = cbind("(Intercept)" = 1, cond = c(0, 0, 0, 1, 1, 1), d = c(1, 1, 1, -1, -1, -1)))), paste0("c", 1:50))
rcf <- ase_glm_test(uc, list(condition = "cond"), RHO_MIN)
ok(all(rcf$status == "not_estimable") && all(is.na(rcf$p_condition)), "condition confounded with cross direction: every gene not_estimable (Rmd 04 stops before this)")

# ---- 7. outbred paired test: per individual LRT, summed over individuals (direction-free), ACAT per gene
sim_pairs <- function(S, I, phi_pair, frac_diff = 0, consistent = FALSE, shift = qlogis(0.8), depth = function(m) round(runif(m, 30, 200))) {
  g <- expand.grid(i = seq_len(I), s = seq_len(S)); g <- g[runif(nrow(g)) < 0.6, ]            # heterozygous individual x SNP cells
  diff_snp <- runif(S) < frac_diff; m <- nrow(g); ph <- sample(c(-1, 1), m, replace = TRUE)   # phase differs between individuals
  l0 <- ifelse(runif(m) < 0.5, ph * qlogis(0.7), 0)                                           # baseline imbalance on the high haplotype
  l1 <- l0 + ifelse(diff_snp[g$s], if (consistent) shift else ph * shift, 0)
  n <- matrix(depth(2 * m), m); p <- plogis(cbind(l0, l1))
  pr <- if (phi_pair > 0) matrix(rbeta(2 * m, p * (1 - phi_pair) / phi_pair, (1 - p) * (1 - phi_pair) / phi_pair), m) else p
  data.frame(snp = paste0("s", rep(g$s, 2)), individual = paste0("i", rep(g$i, 2)), cond = rep(c(0, 1), each = m),
             y = rbinom(2 * m, as.vector(n), as.vector(pr)), n = as.vector(n), diff = rep(diff_snp[g$s], 2), stringsAsFactors = FALSE)
}
pair_run <- function(seed, I, phi, frac_diff = 0, consistent = FALSE) {
  set.seed(seed); s <- sim_pairs(3000, I, phi, frac_diff, consistent)
  gmap <- rep(seq_len(3000), sample(1:4, 3000, TRUE))[1:3000]                                  # genes of 1-4 SNPs
  r <- ase_paired_test(s, RHO_MIN); x <- r$snp; idx <- as.integer(sub("^s", "", x$snp))
  truth <- tapply(s$diff, s$snp, any)[x$snp]; gene <- gmap[idx]
  gp <- tapply(x$p, gene, acat); gnull <- tapply(!truth, gene, all)[names(gp)]
  sg <- !is.na(x$p) & x$p < 0.05
  c(snp_size = mean(sg[!truth]), gene_size = mean(gp[gnull] < 0.05, na.rm = TRUE),
    power = if (frac_diff > 0) mean(sg[truth]) else NA_real_,
    sign_up = if (consistent) mean(x$mean_delta[sg & truth] > 0) else NA_real_,
    tested = mean(x$status == "ok"), phi = mean(r$phi$phi_pair))
}
for (I in c(2, 3, 4, 6)) for (ph in c(0.005, 0.02)) {
  m <- seeds_rbind(1:10, function(s) pair_run(6000 + 10 * I + s, I, ph)); report(m, sprintf("paired outbred I = %d, phi %.3f", I, ph))
  gate(m, "snp_size", sprintf("paired outbred I = %d, phi %.3f, SNP", I, ph)); gate(m, "gene_size", sprintf("paired outbred I = %d, phi %.3f, gene (ACAT)", I, ph))
}
m <- seeds_rbind(1:10, function(s) pair_run(7000 + s, 4, 0.02, 0.1)); report(m, "paired outbred I = 4, 10% SNPs changed, phase-heterogeneous")
gate(m, "snp_size", "paired outbred with 10% changed SNPs: size on the unchanged SNPs")
ok(mean(m[, "power"]) >= 0.8, sprintf("paired outbred power, phase-heterogeneous change 0.5 -> 0.8, I = 4: %.3f (>= 0.8)", mean(m[, "power"])))
m <- seeds_rbind(1:10, function(s) pair_run(7100 + s, 4, 0.02, 0.1, TRUE)); report(m, "paired outbred I = 4, consistent change")
ok(mean(m[, "power"]) >= 0.8 && min(m[, "sign_up"]) >= 0.95, sprintf("paired outbred consistent change: power %.3f (>= 0.8), mean_delta > 0 in %.3f of calls (>= 0.95)", mean(m[, "power"]), min(m[, "sign_up"])))
set.seed(7200); bg <- sim_pairs(300, 3, 0.01)
ed <- data.frame(snp = c("a", "a", "a", "a", "b", "b", "c", "c", "c"), individual = c("i1", "i1", "i2", "i2", "i1", "i1", "i3", "i3", "i1"),
                 cond = c(0, 1, 0, 1, 0, 1, 0, 0, 0), y = c(20, 40, 22, 38, 10, 12, 5, 6, 7), n = c(50, 50, 50, 50, 20, 20, 10, 10, 10), diff = FALSE)
pe <- ase_paired_test(rbind(bg, ed), RHO_MIN)
ok(identical(names(pe$snp), c("snp", "n_individuals", "n_failed", "stat", "df", "p", "mean_delta", "max_abs_delta", "n_up", "n_down", "status")) &&
   identical(names(pe$phi), c("individual", "phi_pair")), "ase_paired_test output columns")
pa <- pe$snp[pe$snp$snp == "a", ]; pb <- pe$snp[pe$snp$snp == "b", ]
ok(pa$status == "ok" && pa$n_individuals == 2 && pa$n_up == 2 && is.finite(pa$p), "SNP in 2 paired individuals: tested, both REF up")
ok(pb$status == "too_few_individuals" && is.na(pb$p), "SNP paired in 1 individual only: too_few_individuals, p NA")
ok(!"c" %in% pe$snp$snp, "a SNP without any paired individual is not reported as tested")
ok(all(pe$phi$phi_pair >= RHO_MIN), "per-individual pair dispersion never below RHO_MIN")

# ---- 8. fragment level: SNPs that share read pairs (ASEReadCounter counts a fragment at every SNP it covers)
frag_counts <- function(pos, p, nf, span, rl = 100) {
  L <- round(runif(nf, 150, 350)); st <- floor(runif(nf) * (span - L + 1)) + 1; A <- runif(nf) < p
  cv <- (outer(st, pos, "<=") & outer(st + rl - 1, pos, ">=")) | (outer(st + L - rl, pos, "<=") & outer(st + L - 1, pos, ">="))
  list(a = colSums(cv & A), n = colSums(cv))
}
frag_run <- function(seed, G = 1500, nA = 3, nB = 3, phi_bio = 0.005, span = 1500) {
  set.seed(seed); d <- c(rep(1, nA), rep(-1, nB)); S <- length(d)
  k <- sample(c(1, 2, 4, 8), G, TRUE); dense <- runif(G) < 0.3
  genes <- lapply(seq_len(G), function(g) {
    pos <- if (dense[g]) sort(sample(600:900, k[g])) else sort(sample(1:span, k[g]))       # dense: every SNP within 300 bp
    pr <- rbeta(S, 0.5 * (1 - phi_bio) / phi_bio, 0.5 * (1 - phi_bio) / phi_bio)        # null: animal-level variation only
    cnt <- lapply(seq_len(S), function(s) frag_counts(pos, pr[s], rpois(1, 150), span))
    list(pos = pos, a = sapply(cnt, `[[`, "a"), n = sapply(cnt, `[[`, "n"))               # SNP x sample matrices
  })
  mat <- function(v) if (is.matrix(v)) v else matrix(v, nrow = 1)
  unit <- function(gi, thin) {
    a <- mat(genes[[gi]]$a); n <- mat(genes[[gi]]$n)
    keep <- if (thin) thin_snps(genes[[gi]]$pos, rowSums(n), THIN_BP) else rep(TRUE, nrow(n))
    y <- colSums(a[keep, , drop = FALSE]); nn <- colSums(n[keep, , drop = FALSE]); ok_s <- nn > 0
    list(y = y[ok_s], n = nn[ok_s], X = cbind("(Intercept)" = 1, d = d)[ok_s, , drop = FALSE])
  }
  tests <- list(strain = "(Intercept)", parent_of_origin = "d")
  rt <- ase_glm_test(setNames(lapply(seq_len(G), unit, thin = TRUE), paste0("g", 1:G)), tests, RHO_MIN)
  ru <- ase_glm_test(setNames(lapply(seq_len(G), unit, thin = FALSE), paste0("g", 1:G)), tests, RHO_MIN)
  sz <- function(p, sel) mean(p[sel] < 0.05, na.rm = TRUE)
  multi <- k >= 2
  rmd02 <- function(thin) {   # the Stage 1 Rmd 02 F1 gene path, per sample: free-mean rho (corrected, floored) + gene LRT
    unlist(lapply(seq_len(S), function(s) {
      x <- integer(0); n <- integer(0); gid <- integer(0)
      for (gi in seq_len(G)) {
        a <- mat(genes[[gi]]$a); nn <- mat(genes[[gi]]$n)
        keep <- if (thin) thin_snps(genes[[gi]]$pos, rowSums(nn), THIN_BP) else rep(TRUE, nrow(nn))
        x <- c(x, a[keep, s]); n <- c(n, nn[keep, s]); gid <- c(gid, rep(gi, sum(keep)))
      }
      e <- bb_estimate_rho_gene(x, n, gid); rho <- max(if (is.na(e[["corrected"]])) 0 else e[["corrected"]], RHO_MIN)
      vapply(which(multi), function(gi) { i <- gid == gi; if (sum(i) < 2) NA_real_ else bb_gene_lrt(x[i], n[i], rho)$p }, numeric(1))
    }))
  }
  p02u <- rmd02(FALSE); p02t <- rmd02(TRUE); dm <- rep(dense[multi], S)
  c(thin_strain = sz(rt$p_strain, TRUE), thin_poo = sz(rt$p_parent_of_origin, TRUE), thin_dense_poo = sz(rt$p_parent_of_origin, dense & multi),
    unthin_strain = sz(ru$p_strain, TRUE), unthin_poo = sz(ru$p_parent_of_origin, TRUE), unthin_dense_poo = sz(ru$p_parent_of_origin, dense & multi),
    rmd02_unthinned = mean(p02u < 0.05, na.rm = TRUE), rmd02_thinned = mean(p02t < 0.05, na.rm = TRUE),
    rmd02_dense_unthinned = mean(p02u[dm] < 0.05, na.rm = TRUE), rmd02_dense_thinned = mean(p02t[dm] < 0.05, na.rm = TRUE))
}
m <- seeds_rbind(1:10, function(s) frag_run(8000 + s)); report(m, "fragment-level null (shared read pairs)")
gate(m, "thin_strain", "fragment level, thinned sums, strain test"); gate(m, "thin_poo", "fragment level, thinned sums, parent-of-origin test")
gate(m, "thin_dense_poo", "fragment level, thinned sums, SNP-dense genes, parent-of-origin test", 0.07, 0.10)
cat(sprintf("TASK5_INPUT rmd02_unthinned=%.4f rmd02_thinned=%.4f rmd02_dense_unthinned=%.4f rmd02_dense_thinned=%.4f unthinned_sum_poo=%.4f\n",
            mean(m[, "rmd02_unthinned"]), mean(m[, "rmd02_thinned"]), mean(m[, "rmd02_dense_unthinned"]), mean(m[, "rmd02_dense_thinned"]), mean(m[, "unthin_poo"])))

# ---- 9. runtime (one core): projected time for 60,000 units x 12 rows must fit the 4 h run-script limit with margin
set.seed(900); ub <- sim_units(5000, rep(c(1, -1), 6), 0.02)
tt <- system.time(rb <- ase_glm_test(ub, list(strain = "(Intercept)", parent_of_origin = "d"), RHO_MIN))[["elapsed"]]
proj <- tt * 60000 / 5000 / 3600
ok(proj < 3, sprintf("runtime: 5000 units x 12 rows, 2 tests, one core: %.0f s; projected 60,000 units %.2f h (< 3 h)", tt, proj))

quit(status = fail)
```

- [ ] **Step 3: Write the test runner**

Create `ase-pipeline/tests/r/run_stats_tests.sh`:

```bash
#!/bin/bash
#SBATCH -J ase_stats_tests
#SBATCH -p bcc
#SBATCH -n 8 --mem=16G -t 3:00:00
# Usage, from the repository root (the log goes where -o says, never into the repository):
#   sbatch -p bcc -o <scratch>/logs/ase_stats_tests_%j.out ase-pipeline/tests/r/run_stats_tests.sh [stage1|stage2|all]
set -uo pipefail
module add singularity/3.10.4 || exit 1
SIF=/net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif
WHAT=${1:-all}; rc=0
cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
if [ "$WHAT" = stage1 ] || [ "$WHAT" = all ]; then
  singularity exec --bind /net/bmc-lab3 "$SIF" Rscript ase-pipeline/tests/r/test_ase_stats.R ase-pipeline/ase-pipeline.md || rc=1
fi
if [ "$WHAT" = stage2 ] || [ "$WHAT" = all ]; then
  singularity exec --bind /net/bmc-lab3 "$SIF" Rscript ase-pipeline/tests/r/test_ase_stats_stage2.R ase-pipeline/ase-pipeline.md || rc=1
fi
echo "EXIT $rc"; exit $rc
```

- [ ] **Step 4: Run the Stage 2 tests to verify they fail**

Run: `sbatch -p bcc -o $SCR/logs/t1_red_%j.out ase-pipeline/tests/r/run_stats_tests.sh stage2`
Expected: the log shows `FAIL: missing from the stats block: chrom_class, thin_snps, ...` and `EXIT 1` (the `const("THIN_BP")` line is never reached).

- [ ] **Step 5: Add the functions to the stats block and the `THIN_BP` constant**

In Step 8's table, add this row directly after the `RHO_MIN` row:

```
| `THIN_BP` | 500 | gene-level tests of Rmd 03 and 04: SNPs closer than this many exonic bases can be counted in the same read pair (ASEReadCounter counts a fragment at every SNP it covers), so only the deepest SNP per window is used; set it to about the longest fragment (insert) length of the library |
```

In Step 14, insert this code directly above the line `# --- ase-stats-end`:

```r
chrom_class <- function(contig) {   # "autosome", "X", "Y" or "MT" (Ensembl and chr-prefixed names); Rmd 03/04 test autosomes only
  s <- toupper(sub("^chr", "", as.character(contig), ignore.case = TRUE))
  ifelse(s == "X", "X", ifelse(s == "Y", "Y", ifelse(s %in% c("MT", "M"), "MT", "autosome")))
}
thin_snps <- function(pos, depth, window) {   # logical keep: deepest SNP first (ties: lower position); no two kept SNPs closer than window
  keep <- rep(FALSE, length(pos))
  for (i in order(-depth, pos)) if (!any(keep & abs(pos - pos[i]) < window)) keep[i] <- TRUE
  keep
}
bb_ll_rows <- function(x, n, p, rho) {   # per-row beta-binomial log-likelihood, p may differ by row; rho < 1e-8 -> binomial
  p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
  if (rho < 1e-8) return(dbinom(x, n, p, log = TRUE))
  a <- p * (1 - rho) / rho; b <- (1 - p) * (1 - rho) / rho
  lchoose(n, x) + lbeta(x + a, n - x + b) - lbeta(a, b)
}
bb_glm_fit <- function(y, n, X, rho, bound = 15) {   # beta-binomial GLM, logit link, dispersion rho FIXED; coefficients boxed to [-bound, bound]
  X <- as.matrix(X)
  if (ncol(X) == 0) return(list(beta = numeric(0), loglik = sum(bb_ll_rows(y, n, rep(0.5, length(y)), rho)), conv = 0L, at_bound = FALSE))
  nll <- function(b) -sum(bb_ll_rows(y, n, plogis(drop(X %*% b)), rho))
  gr <- function(b) {   # analytic gradient (finite differences stall near the bounds of monoallelic genes)
    p <- pmin(pmax(plogis(drop(X %*% b)), 1e-10), 1 - 1e-10)
    s <- if (rho < 1e-8) y - n * p else {
      k <- (1 - rho) / rho
      k * p * (1 - p) * (digamma(y + p * k) - digamma(n - y + (1 - p) * k) - digamma(p * k) + digamma((1 - p) * k))
    }
    -drop(crossprod(X, s))
  }
  st <- tryCatch(as.numeric(qr.solve(X, qlogis((y + 0.5) / (n + 1)))), error = function(e) rep(0, ncol(X)))
  st[!is.finite(st)] <- 0; st <- pmin(pmax(st, -bound + 1), bound - 1)
  o <- optim(st, nll, gr, method = "L-BFGS-B", lower = -bound, upper = bound)
  list(beta = setNames(o$par, colnames(X)), loglik = -o$value, conv = o$convergence, at_bound = any(abs(o$par) > bound - 1e-3))
}
bb_glm_lrt <- function(y, n, X, drop_cols, rho, bound = 15) {   # LRT of H0: coefficients of drop_cols = 0, rho fixed; df = length(drop_cols)
  X <- as.matrix(X)
  f1 <- bb_glm_fit(y, n, X, rho, bound)
  f0 <- bb_glm_fit(y, n, X[, setdiff(colnames(X), drop_cols), drop = FALSE], rho, bound)
  stat <- max(0, 2 * (f1$loglik - f0$loglik))
  list(beta = f1$beta, stat = stat, df = length(drop_cols), p = pchisq(stat, length(drop_cols), lower.tail = FALSE),
       conv = max(f1$conv, f0$conv), at_bound = f1$at_bound)
}
bb_moment_phi <- function(units, max_iter = 25, tol = 1e-6, bound = 15) {
  # pooled moment (Williams-type) estimate of the beta-binomial dispersion over units (each list(y, n, X), rows with n > 0):
  # solves sum(Pearson residual^2 / (1 + (n - 1) phi)) = residual df, the means refitted at each phi. Dividing by the residual
  # df (rows minus coefficients) is what avoids the Neyman-Scott shrinkage of a maximum-likelihood fit with free means.
  phi <- 0
  for (it in seq_len(max_iter)) {
    parts <- lapply(units, function(u) tryCatch({
      f <- bb_glm_fit(u$y, u$n, u$X, phi, bound)
      p <- pmin(pmax(plogis(drop(as.matrix(u$X) %*% f$beta)), 1e-6), 1 - 1e-6)
      r2 <- (u$y - u$n * p)^2 / (u$n * p * (1 - p)); dfu <- as.numeric(length(u$y) - qr(as.matrix(u$X))$rank)
      if (dfu > 0 && all(is.finite(r2))) list(r2 = r2, n = u$n, df = dfu) else NULL
    }, error = function(e) NULL))
    parts <- parts[!vapply(parts, is.null, logical(1))]
    if (length(parts) == 0) return(NA_real_)
    r2 <- unlist(lapply(parts, `[[`, "r2")); nn <- unlist(lapply(parts, `[[`, "n")); df <- sum(vapply(parts, `[[`, numeric(1), "df"))
    g <- function(ph) sum(r2 / (1 + (nn - 1) * ph)) - df
    new <- if (g(0) <= 0) 0 else if (g(0.99) > 0) 0.99 else uniroot(g, c(0, 0.99), tol = 1e-12)$root
    if (abs(new - phi) < tol) return(new)
    phi <- new
  }
  warning("dispersion estimate did not converge in ", max_iter, " iterations; the last value is used")
  phi
}
ase_glm_test <- function(units, tests, rho_min, min_df_unit = 2, bound = 15) {
  # units: named list of list(y, n, X), one row per sample (rows with n > 0 only), X with named columns
  # tests: named list, test name -> the X columns set to 0 under H0 (likelihood-ratio test, df = number of columns)
  # dispersion: phi_common = pooled moment estimate over every estimable unit; per unit phi_used = max(phi_common, rho_min,
  # phi_unit when the unit has at least min_df_unit residual df): the larger of the shared and the unit's own estimate
  nr <- vapply(units, function(u) length(u$y), integer(1))
  rk <- vapply(units, function(u) if (length(u$y) == 0) 0L else qr(as.matrix(u$X))$rank, integer(1))
  est <- nr > 0 & rk == vapply(units, function(u) ncol(as.matrix(u$X)), integer(1)) & nr > rk
  common <- if (any(est)) bb_moment_phi(units[est], bound = bound) else NA_real_
  cn <- unique(unlist(lapply(units, function(u) colnames(u$X))))
  one <- function(i) {
    u <- units[[i]]
    row <- c(list(unit = names(units)[i], n_rows = nr[i], df_resid = nr[i] - rk[i], depth = sum(u$n),
                  phi_common = common, phi_unit = NA_real_, phi_used = NA_real_, status = "not_estimable"),
             setNames(as.list(rep(NA_real_, length(cn))), paste0("beta_", cn)),
             setNames(as.list(rep(NA_real_, length(tests))), paste0("stat_", names(tests))),
             setNames(as.list(rep(NA_real_, length(tests))), paste0("p_", names(tests))))
    if (est[i] && !is.na(common)) {
      res <- tryCatch({
        pu <- if (nr[i] - rk[i] >= min_df_unit) bb_moment_phi(list(u), bound = bound) else NA_real_
        ph <- max(c(common, rho_min, pu), na.rm = TRUE)
        list(pu = pu, ph = ph, tt = lapply(tests, function(cols) bb_glm_lrt(u$y, u$n, u$X, cols, ph, bound)))
      }, error = function(e) NULL)
      if (is.null(res)) row$status <- "fit_error" else {
        row$phi_unit <- res$pu; row$phi_used <- res$ph
        b <- res$tt[[1]]$beta; row[paste0("beta_", names(b))] <- as.list(unname(b))
        conv <- max(vapply(res$tt, function(t) t$conv, numeric(1)))
        for (nm in names(tests)) {
          row[[paste0("stat_", nm)]] <- res$tt[[nm]]$stat
          row[[paste0("p_", nm)]] <- if (conv == 0) res$tt[[nm]]$p else NA_real_
        }
        row$status <- if (conv != 0) "not_converged" else if (any(vapply(res$tt, function(t) t$at_bound, logical(1)))) "ok_at_bound" else "ok"
      }
    }
    as.data.frame(row, check.names = FALSE, stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, lapply(seq_along(units), one)); rownames(out) <- NULL; out
}
ase_paired_test <- function(d, rho_min, min_individuals = 2, bound = 15) {
  # d: data.frame(snp, individual, cond, y, n); cond 0 = reference condition, 1 = tested condition; y = REF count, n = REF + ALT
  # 1. per individual: pair dispersion phi_i, moment estimate with a free REF fraction per SNP, floored at rho_min
  # 2. per SNP and individual: LRT of one shared REF fraction against one per condition (1 df) with phi_i
  # 3. per SNP: the statistics summed over the informative individuals (df = their number). Direction-free, because
  #    which allele carries a regulatory variant differs between individuals (a shared slope would cancel real changes)
  empty <- data.frame(snp = character(0), n_individuals = integer(0), n_failed = integer(0), stat = numeric(0), df = integer(0),
                      p = numeric(0), mean_delta = numeric(0), max_abs_delta = numeric(0), n_up = integer(0), n_down = integer(0),
                      status = character(0), stringsAsFactors = FALSE)
  d <- d[is.finite(d$y) & is.finite(d$n) & d$n > 0, c("snp", "individual", "cond", "y", "n"), drop = FALSE]
  key <- paste(d$individual, d$snp, sep = "\t")
  paired <- tapply(d$cond, key, function(v) any(v == 0) & any(v == 1))
  d <- d[as.logical(paired[key]), , drop = FALSE]
  if (nrow(d) == 0) return(list(snp = empty, phi = data.frame(individual = character(0), phi_pair = numeric(0))))
  key <- paste(d$individual, d$snp, sep = "\t")
  cells <- split(seq_len(nrow(d)), key); first <- vapply(cells, `[`, integer(1), 1)
  cell_ind <- d$individual[first]; cell_snp <- d$snp[first]
  inds <- sort(unique(cell_ind))
  phi <- vapply(inds, function(ind) {
    u <- lapply(cells[cell_ind == ind], function(i) list(y = d$y[i], n = d$n[i], X = matrix(1, length(i), 1, dimnames = list(NULL, "(Intercept)"))))
    ph <- bb_moment_phi(u, bound = bound)
    if (is.na(ph)) NA_real_ else max(ph, rho_min)
  }, numeric(1))
  names(phi) <- inds
  cr <- do.call(rbind, lapply(seq_along(cells), function(j) {
    i <- cells[[j]]; ph <- phi[[cell_ind[j]]]
    r <- if (is.na(ph)) NULL else tryCatch(bb_glm_lrt(d$y[i], d$n[i], cbind("(Intercept)" = 1, cond = d$cond[i]), "cond", ph, bound),
                                           error = function(e) NULL)
    if (is.null(r) || r$conv != 0) return(data.frame(snp = cell_snp[j], stat = NA_real_, delta = NA_real_, stringsAsFactors = FALSE))
    data.frame(snp = cell_snp[j], stat = r$stat, delta = plogis(sum(r$beta)) - plogis(r$beta[[1]]), stringsAsFactors = FALSE)
  }))
  snp <- do.call(rbind, lapply(split(cr, cr$snp), function(s) {
    k <- sum(is.finite(s$stat)); so <- s[is.finite(s$stat), , drop = FALSE]; enough <- k >= min_individuals
    data.frame(snp = s$snp[1], n_individuals = k, n_failed = nrow(s) - k,
               stat = if (enough) sum(so$stat) else NA_real_, df = k,
               p = if (enough) pchisq(sum(so$stat), k, lower.tail = FALSE) else NA_real_,
               mean_delta = if (k > 0) mean(so$delta) else NA_real_, max_abs_delta = if (k > 0) max(abs(so$delta)) else NA_real_,
               n_up = sum(so$delta > 0), n_down = sum(so$delta < 0),
               status = if (enough) "ok" else "too_few_individuals", stringsAsFactors = FALSE)
  }))
  rownames(snp) <- NULL
  list(snp = snp, phi = data.frame(individual = inds, phi_pair = unname(phi), stringsAsFactors = FALSE))
}
```

Change the Step 14 intro sentence "These base-R functions are pasted verbatim into Rmd 02 (`{TODAY}_{WD_NAME}_02_imbalance.Rmd`)" to "These base-R functions are pasted verbatim into Rmd 02, Rmd 03 and Rmd 04 (`{TODAY}_{WD_NAME}_02_imbalance.Rmd`, `..._03_reciprocal.Rmd`, `..._04_differential.Rmd`)".

- [ ] **Step 6: Run both test files and handle the gates**

Run: `sbatch -p bcc -o $SCR/logs/t1_green_%j.out ase-pipeline/tests/r/run_stats_tests.sh all`
Expected: the Stage 1 file still prints only `ok` lines (62 ok). The Stage 2 file prints only `ok` lines, the `report`/`EVIDENCE` lines, and one `TASK5_INPUT` line. The log ends with `EXIT 0`.

Code defects: if a test fails because of a coding defect in a new function (an R error, a wrong column name, a status not produced), fix the function and re-run.

Gate failures: if a simulation gate fails (a size, power or sign line), change nothing and report BLOCKED with the full log. The controller rules; the options are a documented method change or a gate ruling by the user.

Runtime: the whole job is expected to take under 2 hours. If `projected 60,000 units` is at or above 3 h, report BLOCKED as well (the run scripts of Tasks 3 and 4 are single-threaded with `-t 4:00:00`).

Copy the log to `$SCR/logs/t1_final.log`.

- [ ] **Step 7: Document the measured behaviour in Step 14**

Append these bullets after the last Step 14 bullet. Fill every number from `$SCR/logs/t1_final.log`, copying the printed values; nothing is estimated by hand:

- **Stage 2 models (Rmd 03 and 04): `ase_glm_test`.** Each gene (or SNP) is one unit, with one row per sample. For F1 genes, a row is the sum of strain-A and total counts over the gene's SNPs after `thin_snps` (one SNP per `THIN_BP` window in exon coordinates, deepest first), so a read pair is not counted twice. The beta-binomial GLM (`bb_glm_fit`: logit link, dispersion fixed, coefficients boxed to ±15 so that monoallelic genes give a finite likelihood-ratio statistic) is tested by `bb_glm_lrt`. The dispersion is `phi_common` from `bb_moment_phi`: the spread between replicate samples, divided by the residual df, so there is no Neyman-Scott shrinkage. Each unit uses `phi_used = max(phi_common, RHO_MIN, phi_unit)`, and `phi_unit` counts only with at least 2 residual df. Measured, 10 seeds x 2000 null genes per design: the reciprocal null sizes (strain / parent of origin) are `<numbers from the "reciprocal ..." report lines for 3+3, 2+2, 3+2 and 4+2>`, and with heterogeneous dispersion `<numbers>`. Power for 3+3 is `<numbers>`. The F1 differential null size is `<numbers>`, including imprinted genes with the cross-direction term. Without the term, the imprinted-gene size in an unbalanced design is `<EVIDENCE number>`.
- **Outbred differential: `ase_paired_test`.** A per-individual pair dispersion (moment estimate, free REF fraction per SNP, floored at `RHO_MIN`) is used for a 1-df LRT per SNP and individual. The statistics are summed over individuals (df = number of individuals), which is direction-free because the allele carrying a regulatory variant differs between individuals. Genes use `acat`. Measured: the SNP and gene null sizes are `<numbers for I = 2, 3, 4, 6>`, and power is `<phase-heterogeneous>` / `<consistent>` for a 0.5 -> 0.8 change with 4 individuals.
- **SNPs sharing read pairs.** In a fragment-level simulation (ASEReadCounter-style counting, 30 percent SNP-dense genes), the thinned gene tests have null sizes of `<thin_strain>` / `<thin_poo>` (dense genes `<thin_dense_poo>`), against `<unthin_poo>` without thinning. The Stage 1 Rmd 02 F1 gene LRT has a size of `<rmd02_unthinned>` unthinned and `<rmd02_thinned>` thinned.

- [ ] **Step 8: Checker lines, proven RED, then GREEN**

Append to `ase-pipeline/tests/check_skill.sh` before the final `[ $fail -eq 0 ]` line:

```bash
# --- Stage 2 Task 1 (statistics); each string is absent from the Stage 1 skill (152dc20)
need "chrom_class <- function(contig)"
need "thin_snps <- function(pos, depth, window)"
need "bb_glm_fit <- function(y, n, X, rho, bound = 15)"
need "bb_glm_lrt <- function(y, n, X, drop_cols, rho, bound = 15)"
need "bb_moment_phi <- function(units, max_iter = 25, tol = 1e-6, bound = 15)"
need "ase_glm_test <- function(units, tests, rho_min, min_df_unit = 2, bound = 15)"
need "ase_paired_test <- function(d, rho_min, min_individuals = 2, bound = 15)"
need '| `THIN_BP` | 500 |'
need "pasted verbatim into Rmd 02, Rmd 03 and Rmd 04"
need "**Stage 2 models (Rmd 03 and 04): \`ase_glm_test\`.**"
need "**Outbred differential: \`ase_paired_test\`.**"
need "**SNPs sharing read pairs.**"
[ -s "$HERE/r/test_ase_stats_stage2.R" ] || { echo "FAIL: missing tests/r/test_ase_stats_stage2.R"; fail=1; }
# --- end Stage 2 Task 1
```
The Task 1 block uses `$HERE`, which the checker defines in the Task 7 block (`HERE=$(dirname "$0")`), so the new block must sit after that line. Proof: run the checker on `git show HEAD:ase-pipeline/ase-pipeline.md > $SCR/pre.md`. Expected: 12 new `FAIL` lines. Then run it on the working copy. Expected: `PASS`.

- [ ] **Step 9: Amend the spec**

In `docs/superpowers/specs/2026-09-29-ase-pipeline-design.md`, make three replacements:

- Replace the paragraph that starts `**Rmd 03 — Reciprocal F1` with:

  "**Rmd 03 — Reciprocal F1 (F1 mode, at least 2 samples in each cross direction).** Per gene, a beta-binomial GLM `logit(p) = β0 + β1·d` (plus sum-to-zero condition terms when several conditions exist). Here p is the strain-A allele fraction, and d = +1 for `AxB` (strain A is the mother) and −1 for `BxA`. **β0 is the strain effect (cis-regulatory divergence); β1 is the parent-of-origin effect (imprinting).** Each sample contributes one row per gene: the counts summed over the gene's SNPs after thinning to one SNP per `THIN_BP` exon-coordinate bases (so no read pair is counted twice). The dispersion is estimated from the spread between replicate animals (pooled moment estimate, per-gene maximum rule, floored at `RHO_MIN`), and β0 and β1 are tested by likelihood-ratio tests. X, Y and MT are excluded (hemizygous X and maternal MT mimic imprinting)."
- Replace the paragraph that starts `**Rmd 04 — Differential ASE` with:

  "**Rmd 04 — Differential ASE (at least two conditions with replicates).** Each non-reference condition is compared with a user-chosen reference condition.
  - **F1:** per gene (thinned, summed rows as in Rmd 03) and per SNP, `logit(p) = β0 + β1·condition` (+ the cross-direction term when both directions are present; the Rmd stops when condition and direction are confounded), tested by a likelihood-ratio test with the Rmd 03 dispersion rules.
  - **Outbred (paired individuals):** per SNP and individual, a beta-binomial LRT of one REF fraction against one per condition, with a per-individual pair dispersion. The statistics are summed over individuals, so the test is direction-free, because the allele carrying a regulatory variant differs between individuals. Genes use ACAT and are labelled 'unphased, no direction'."
- Replace the sentence that starts `The beta-binomial GLMs in Rmd 03 and 04 use` (up to "no image rebuild is needed.") with:

  "The beta-binomial GLMs in Rmd 03 and 04 are fitted by base-R code in the shipped statistics block (bounded L-BFGS-B, dispersion fixed at the estimate above). `aod::betabin` (aod 1.3.3, in the `bulkrnaseq` image) cross-checks that likelihood in the unit tests. A free per-gene aod dispersion from 4-12 rows is biased low, and its unbounded fit fails on monoallelic genes. `lme4` is not used: a random individual effect is not estimable with 2-6 individuals, and a shared condition slope cancels phase-dependent changes."

- [ ] **Step 10: Commit**

```bash
git status && git diff --stat
git add ase-pipeline/ase-pipeline.md ase-pipeline/tests/r/test_ase_stats_stage2.R ase-pipeline/tests/r/run_stats_tests.sh ase-pipeline/tests/check_skill.sh docs/superpowers/specs/2026-09-29-ase-pipeline-design.md
git commit -m "ase-pipeline: Stage 2 statistics (beta-binomial GLM, moment dispersion, paired outbred test) with simulation tests

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Stage 2 synthetic ground truth (reciprocal x two-condition F1, paired two-condition outbred)

**Files:**
- Create: `ase-pipeline/tests/synthetic/simulate_ase_stage2.R`, `ase-pipeline/tests/synthetic/check_simulation_stage2.sh`, `ase-pipeline/tests/synthetic/check_biology_stage2.R`, `ase-pipeline/tests/synthetic/prove_checks_stage2.sh`
- Do not modify `simulate_ase_data.R`: Stage 1 data must stay byte-identical for its seed.

**Interfaces:**
- Produces: `Rscript simulate_ase_stage2.R <outdir> <seed>` writes the files below. Tasks 3, 4 and 8 depend on these names.

  | Path under `<outdir>` | Content |
  |---|---|
  | `seed.txt` | the seed |
  | `genome/genome.fa`, `genome/genome.gtf` | chr1, 1,000,000 bp, 60 plus-strand genes of 4-5 exons of 150-300 bp (spliced length ≥ 600) |
  | `f1_recip/samples.csv` | `sample,fastq_1,fastq_2,condition,cross_direction,individual`; 12 rows `f1r_s01`..`f1r_s12`: `AxB`/`ctrl` x3, `AxB`/`treat` x3, `BxA`/`ctrl` x3, `BxA`/`treat` x3; `individual` = `sample` |
  | `f1_recip/{sample}_1.fastq.gz`, `_2.fastq.gz` | paired reads, length 100 |
  | `f1_recip/parental_snps.vcf` | sites-only, REF = strain A allele, ALT = strain B allele |
  | `f1_recip/truth_genes.tsv` | `gene_id n_snps class b0 b1 bc p_A_AxB_ctrl p_A_AxB_treat p_A_BxA_ctrl p_A_BxA_treat expect_strain expect_poo expect_diff` |
  | `f1_recip/truth_snps.tsv` | `chrom pos ref alt gene_id tx_off` |
  | `f1_recip/truth_fragments.tsv` | `sample gene_id frag_A frag_B` (realised fragment counts per allele) |
  | `f1_recip/direct_counts/{sample}.table` | ASEReadCounter columns, counted from the simulated fragments |
  | `outbred_diff/samples.csv` | 8 rows `ob_s1`..`ob_s8`: individuals `ind1`..`ind4`, each once in `ctrl` and once in `treat`; `cross_direction` = `NA` |
  | `outbred_diff/{sample}_1/_2.fastq.gz` | as above |
  | `outbred_diff/{ind}.all.vcf`, `{ind}.het.vcf` | as in Stage 1 |
  | `outbred_diff/truth_genes.tsv` | `gene_id n_snps class f_ctrl f_treat consistent expect_diff` |
  | `outbred_diff/truth_snps.tsv` | as for F1 |
  | `outbred_diff/truth_individual_genes.tsv` | `individual gene_id phase n_het p_alt_ctrl p_alt_treat` |
  | `outbred_diff/truth_fragments.tsv` | as for F1 |
  | `outbred_diff/direct_counts/{sample}.table`, `.unfiltered.table`, `.wasp_stats.tsv` | the unfiltered table is a copy of the table; the WASP stats have header `vW vA n` and rows `1 1 <REF fragments>`, `1 2 <ALT fragments>`, `none - 0` |

  Direct counts ignore sequencing errors. They test the Rmds, not the aligner. Each row has `lowMAPQDepth`, `lowBaseQDepth`, `otherBases` and `improperPairs` = 0 and `rawDepth` = `totalCount`, and is written only when `totalCount` ≥ 1.

- [ ] **Step 1: Write the check script first**

Create `ase-pipeline/tests/synthetic/check_simulation_stage2.sh`:

```bash
#!/bin/bash
# Usage: check_simulation_stage2.sh <outdir>   (after simulate_ase_stage2.R); exact quotas, designs and consistency
set -u; D=${1:?usage: check_simulation_stage2.sh <outdir>}; fail=0
bad() { echo "FAIL: $*"; fail=1; }
chk() { [ -s "$1" ] || bad "missing/empty $1"; }
export LC_ALL=C
F=$D/f1_recip; O=$D/outbred_diff
chk $D/seed.txt; chk $D/genome/genome.fa; chk $D/genome/genome.gtf; chk $F/parental_snps.vcf; chk $F/truth_fragments.tsv; chk $O/truth_fragments.tsv
for m in $F $O; do
  chk $m/samples.csv; chk $m/truth_genes.tsv; chk $m/truth_snps.tsv
  [ "$(head -1 $m/samples.csv 2>/dev/null)" = "sample,fastq_1,fastq_2,condition,cross_direction,individual" ] || bad "$m/samples.csv header"
  for s in $(tail -n +2 $m/samples.csv 2>/dev/null | cut -d, -f1); do
    chk $m/${s}_1.fastq.gz; chk $m/${s}_2.fastq.gz; chk $m/direct_counts/$s.table
    n1=$(zcat $m/${s}_1.fastq.gz 2>/dev/null | awk 'END{print NR/4}'); n2=$(zcat $m/${s}_2.fastq.gz 2>/dev/null | awk 'END{print NR/4}')
    nf=$(awk -F'\t' -v s=$s '$1==s{t+=$3+$4} END{print t+0}' $m/truth_fragments.tsv)
    [ "$n1" = "$n2" ] && [ "$n1" = "$nf" ] || bad "$s read pairs ($n1 / $n2) differ from truth_fragments ($nf)"
    for g in $m/${s}_1.fastq.gz $m/${s}_2.fastq.gz; do
      [ "$(zcat $g 2>/dev/null | awk 'NR%4==2{print length($0)}' | sort -u)" = 100 ] || bad "$g read length not exactly 100"
    done
    awk -F'\t' 'NR==1{ if ($1!="contig"||$6!="refCount"||$7!="altCount"||$8!="totalCount") exit 1; next }
                { if ($6+$7!=$8 || $11!=$8 || $8<1) exit 1 }' $m/direct_counts/$s.table || bad "$m/direct_counts/$s.table header or ref+alt/total/rawDepth"
  done
done
# F1 design: exactly 3 samples per direction x condition; individual = sample
got=$(tail -n +2 $F/samples.csv | awk -F, '{print $5"_"$4}' | sort | uniq -c | awk '{print $2":"$1}' | paste -sd' ')
[ "$got" = "AxB_ctrl:3 AxB_treat:3 BxA_ctrl:3 BxA_treat:3" ] || bad "f1_recip design '$got'"
tail -n +2 $F/samples.csv | awk -F, '$6!=$1{e=1} END{exit e}' || bad "f1_recip: individual must equal sample"
# outbred design: 4 individuals, each exactly once in ctrl and once in treat, cross_direction NA
tail -n +2 $O/samples.csv | awk -F, '{c[$6","$4]++; i[$6]++; if ($5!="NA") e=1}
  END{ n=0; for (k in i) { n++; if (c[k",ctrl"]!=1 || c[k",treat"]!=1) e=1 } if (n!=4 || NR!=8) e=1; exit e }' || bad "outbred_diff design (4 paired individuals, 8 rows, cross_direction NA)"
# F1 truth: exact planted classes and values; cell fractions = plogis(b0 + b1 d + bc t)
keys=$(awk -F'\t' 'function pl(x){return 1/(1+exp(-x))}
  NR>1{ print $3":"sprintf("%.2f",pl($4))":"sprintf("%.2f",pl($5))":"sprintf("%.2f",pl($6))
        for (j=0;j<4;j++){ d=(j<2)?1:-1; t=j%2; p=pl($4+$5*d+$6*t); if ((p-$(7+j))^2>1e-8) print "BADCELL:"$1 } }' $F/truth_genes.tsv | sort | uniq -c | awk '{print $2"="$1}' | paste -sd' ')
want="bias:0.50:0.50:0.50=1 dense_null:0.50:0.50:0.50=2 diff_down:0.50:0.50:0.25=1 diff_on_strain:0.70:0.50:0.30=1 diff_up:0.50:0.50:0.75=2 maternal:0.50:0.85:0.50=1 maternal:0.50:0.95:0.50=1 null:0.50:0.50:0.50=44 paternal:0.50:0.10:0.50=1 paternal:0.50:0.15:0.50=1 strain_A_high:0.70:0.50:0.50=1 strain_A_high:0.75:0.50:0.50=1 strain_B_high:0.25:0.50:0.50=1 strain_B_high:0.30:0.50:0.50=1 strain_and_maternal:0.70:0.80:0.50=1"
[ "$keys" = "$want" ] || bad "f1_recip truth classes/values: got '$keys'"
awk -F'\t' 'NR>1 && (($3=="dense_null" && $2!=8) || ($3=="bias" && $2!=5) || ($3!="dense_null" && $3!="bias" && $2!=4)){e=1} END{exit e}' $F/truth_genes.tsv || bad "f1_recip n_snps (8 dense, 5 bias, 4 otherwise)"
# dense_null genes: 8 SNPs within 300 transcript bases; bias gene: 5 SNPs within 60
awk -F'\t' 'NR==FNR{ if (FNR>1) c[$1]=$3; next } FNR>1{ k=$5; if (!(k in lo) || $6<lo[k]) lo[k]=$6; if (!(k in hi) || $6>hi[k]) hi[k]=$6 }
  END{ for (k in c) { if (c[k]=="dense_null" && hi[k]-lo[k]>=300) e=1; if (c[k]=="bias" && hi[k]-lo[k]>=60) e=1 } exit e }' $F/truth_genes.tsv $F/truth_snps.tsv || bad "dense/bias SNP spans"
# outbred truth: classes and values
keys=$(awk -F'\t' 'NR>1{print $3":"sprintf("%.2f",$4)":"sprintf("%.2f",$5)":"$6}' $O/truth_genes.tsv | sort | uniq -c | awk '{print $2"="$1}' | paste -sd' ')
want="base_imbalanced:0.75:0.75:FALSE=3 bias:0.50:0.50:FALSE=1 diff_consistent:0.50:0.75:TRUE=2 diff_phase:0.50:0.80:FALSE=3 null:0.50:0.50:FALSE=51"
[ "$keys" = "$want" ] || bad "outbred_diff truth classes/values: got '$keys'"
# outbred per-individual truth: p_alt follows phase; every planted gene has >= 3 individuals with a het SNP; diff_phase phases mixed
awk -F'\t' 'NR==FNR{ if (FNR>1){ c[$1]=$3; fc[$1]=$4; ft[$1]=$5; cons[$1]=$6 } next }
  FNR>1{ g=$2; if (cons[g]=="TRUE") { ec=fc[g]; et=ft[g] } else if ($3==1) { ec=fc[g]; et=ft[g] } else { ec=1-fc[g]; et=1-ft[g] }
         if ($4>0 && ((($5-ec)^2>1e-8) || (($6-et)^2>1e-8))) e=1
         if ($4>0) { nh[g]++; if ($3==1) up[g]=1; else dn[g]=1 } }
  END{ for (g in c) { if (c[g]!="null" && c[g]!="bias" && nh[g]<3) e=1; if (c[g]=="diff_phase" && !(up[g] && dn[g])) e=1 } exit e }' \
  $O/truth_genes.tsv $O/truth_individual_genes.tsv || bad "outbred_diff truth_individual_genes (phase, p_alt, >= 3 het individuals, mixed phases)"
# VCFs: REF = genome base; parental VCF = truth SNPs; het VCF = exactly the 0/1 rows of the all VCF
awk 'NR==FNR{ if(/^>/) next; g=g $0; next } !/^#/{ if(substr(g,$2,1)!=$4) e=1 } END{exit e}' $D/genome/genome.fa $F/parental_snps.vcf || bad "parental_snps.vcf REF differs from genome"
vcf_rows() { grep -v '^#' "$1" | awk -F'\t' '{print $1"\t"$2"\t"$4"\t"$5}' | sort; }
truth_rows() { tail -n +2 "$1" | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4}' | sort; }
diff <(vcf_rows $F/parental_snps.vcf) <(truth_rows $F/truth_snps.tsv) >/dev/null || bad "parental_snps.vcf differs from truth_snps.tsv"
for i in ind1 ind2 ind3 ind4; do
  chk $O/$i.all.vcf; chk $O/$i.het.vcf
  diff <(vcf_rows $O/$i.all.vcf) <(truth_rows $O/truth_snps.tsv) >/dev/null || bad "$i.all.vcf differs from truth_snps.tsv"
  diff <(grep -v '^#' $O/$i.all.vcf | awk -F'\t' '$10=="0/1"') <(grep -v '^#' $O/$i.het.vcf) >/dev/null || bad "$i.het.vcf is not the 0/1 rows of $i.all.vcf"
done
# direct counts: positions are truth SNPs (outbred: that individual's het SNPs); F1 cell fractions near the planted values
for s in $(tail -n +2 $O/samples.csv | cut -d, -f1); do
  ind=$(awk -F, -v s=$s '$1==s{print $6}' $O/samples.csv)
  awk -F'\t' 'NR==FNR{ if (!/^#/ && $10=="0/1") h[$2]=1; next } FNR>1 && !($2 in h){e=1} END{exit e}' $O/$ind.het.vcf $O/direct_counts/$s.table || bad "$s direct counts at non-het positions"
  cmp -s $O/direct_counts/$s.table $O/direct_counts/$s.unfiltered.table || bad "$s unfiltered.table must equal the table"
  [ "$(head -1 $O/direct_counts/$s.wasp_stats.tsv 2>/dev/null)" = "$(printf 'vW\tvA\tn')" ] || bad "$s wasp_stats.tsv header"
done
awk -F'\t' 'FNR==1{ f++; next } f==1{ gene[$2]=$5; next }                       # truth_snps: pos -> gene
  f==2{ split($0,a,","); if (a[1]!="sample") cell[a[1]]=a[5]"_"a[4]; next }     # samples.csv read as tab file: whole line in $0
  f==3{ p[$1"|AxB_ctrl"]=$7; p[$1"|AxB_treat"]=$8; p[$1"|BxA_ctrl"]=$9; p[$1"|BxA_treat"]=$10; next }
  { s=FILENAME; sub(/.*\//,"",s); sub(/\.table$/,"",s); k=gene[$2]"|"cell[s]; A[k]+=$6; N[k]+=$8 }
  END{ for (k in N) if (N[k]>=200 && (A[k]/N[k]-p[k])^2>0.08^2) { print "  " k, A[k]/N[k], p[k]; e=1 } exit e }' \
  $F/truth_snps.tsv $F/samples.csv $F/truth_genes.tsv $F/direct_counts/*.table || bad "F1 direct-count strain-A fractions differ from truth by more than 0.08"
[ $fail -eq 0 ] && echo PASS || exit 1
```

Run it on an empty directory (`mkdir -p $SCR/empty && bash ase-pipeline/tests/synthetic/check_simulation_stage2.sh $SCR/empty`). Expected: `FAIL` lines and exit 1.

- [ ] **Step 2: Write the generator**

`simulate_ase_stage2.R` (R, `Biostrings` and base R; about 300 lines) follows this algorithm, seeded by the second argument (written to `seed.txt`):

1. **Constants:** `READLEN <- 100L; ERR <- 0.002; LAMBDA <- 450; PHI_BIO <- 0.005; NGENES <- 60L; GLEN <- 1000000L`.
2. **Genome:** the Stage 1 synthetic-genome block (`simulate_ase_data.R`, the `if (SYNTH) { ... }` block), with `GLEN`, `NGENES` genes, and `ne <- sample(4:5, 1); el <- sample(150:300, ne, replace = TRUE)`. Every gene then has a spliced length of at least 600.
3. **Copied functions:** copy verbatim from `simulate_ase_data.R` the functions `alt_of`, `spliced`, `revcomp`, `haplo`, `add_err`, `write_fq` and `vcf_hdr`, and the gene-list construction (`genes[[g]]` with `start`, `end`, `strand`, `len`).
4. **Classes:** shuffle the 60 gene ids with `sample()`. Assign in order: 1 `bias` gene, 2 `dense_null` genes, then the F1 classes of the table below, then `null` for the rest.

   ```r
   plan <- data.frame(
     class = c("strain_A_high", "strain_A_high", "strain_B_high", "strain_B_high", "maternal", "maternal", "paternal", "paternal",
               "strain_and_maternal", "diff_up", "diff_up", "diff_down", "diff_on_strain"),
     fA = c(0.75, 0.70, 0.25, 0.30, 0.5, 0.5, 0.5, 0.5, 0.70, 0.5, 0.5, 0.5, 0.70),       # strain-A fraction (b0 = qlogis(fA))
     fM = c(0.5, 0.5, 0.5, 0.5, 0.85, 0.95, 0.15, 0.10, 0.80, 0.5, 0.5, 0.5, 0.5),        # maternal fraction (b1 = qlogis(fM))
     fC = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.75, 0.75, 0.25, 0.30),        # condition effect (bc = qlogis(fC))
     expect_strain = c("A_higher", "A_higher", "B_higher", "B_higher", "none", "none", "none", "none", "A_higher", "any", "any", "any", "any"),
     expect_poo = c("none", "none", "none", "none", "maternal", "maternal", "paternal", "paternal", "maternal", "none", "none", "none", "none"),
     expect_diff = c(rep("none", 9), "up", "up", "down", "down"), stringsAsFactors = FALSE)
   # bias, dense_null and null genes: b0 = b1 = bc = 0, expect_* = "none"
   ```

   Write `b0`, `b1` and `bc` with `format(x, digits = 12)`, and the four cell fractions `p_A_<dir>_<cond> = plogis(b0 + b1 * d + bc * t)` (d = +1 `AxB`, −1 `BxA`; t = 1 `treat`).
5. **SNPs (transcript offsets):**
   - ordinary genes: 4 offsets in `[1, len]`, pairwise ≥ 40 apart (`repeat` until true);
   - `bias` gene: 5 offsets within 60 consecutive bases;
   - `dense_null` genes: 8 offsets within 300 consecutive bases, pairwise ≥ 20 apart.

   Map each offset to the genome with `i <- max(which(cum < off)); pos <- start[i] + off - cum[i] - 1`, where `cum <- c(0, cumsum(end - start + 1))`. The REF base is the genome base, and the ALT base is `alt_of(REF)`. Write `truth_snps.tsv` sorted by `pos`, with `tx_off`.
6. **Fragments that also give direct counts:** use this function instead of Stage 1's `frags`:

   ```r
   frags2 <- function(hseq, n, snp_offs) {   # read pairs plus, per fragment, the SNP offsets either mate covers (ASEReadCounter rule)
     L <- nchar(hseq)
     if (n == 0) return(list(r1 = character(0), r2 = character(0), cover = matrix(FALSE, 0, length(snp_offs))))
     fl <- pmin(pmax(round(rnorm(n, 250, 30)), 150), 350, L); st <- floor(runif(n) * (L - fl + 1)) + 1
     r1 <- substring(hseq, st, st + READLEN - 1); r2 <- revcomp(substring(hseq, st + fl - READLEN, st + fl - 1))
     cover <- (outer(st, snp_offs, "<=") & outer(st + READLEN - 1, snp_offs, ">=")) |
              (outer(st + fl - READLEN, snp_offs, "<=") & outer(st + fl - 1, snp_offs, ">="))
     list(r1 = r1, r2 = r2, cover = cover)
   }
   ```

   All genes are on the plus strand, so transcript offsets equal spliced offsets. Per sample and gene, draw `n <- rpois(1, LAMBDA)`, then `nA <- rbinom(1, n, p_sg)` fragments from the REF-side haplotype (strain A / REF) and `n - nA` from the ALT-side haplotype. Add `colSums(cover)` of the REF-side fragments to `refCount` and of the ALT-side fragments to `altCount` at those SNPs. Append `(sample, gene_id, nA, n - nA)` to `truth_fragments.tsv` (for outbred, `frag_A` = REF-side fragments). Apply `add_err` to the reads only.
7. **F1:** 12 samples in the design of the interface table. Per sample and gene, `p_sg <- rbeta(1, p * (1 - PHI_BIO) / PHI_BIO, (1 - p) * (1 - PHI_BIO) / PHI_BIO)`, with p the cell's strain-A fraction. The F1 haplotypes are all-REF and all-ALT (as in Stage 1). Write `parental_snps.vcf` exactly as Stage 1 does.
8. **Outbred:**
   - Shuffle the genes again. Assign the Stage 1-style `bias` gene (the same 5-SNP gene as in F1), 3 `base_imbalanced` (`f_ctrl = f_treat = 0.75`, `consistent = FALSE`), 3 `diff_phase` (0.50 → 0.80, `FALSE`) and 2 `diff_consistent` (0.50 → 0.75, `TRUE`); the rest are `null` (0.5, 0.5, `FALSE`).
   - Genotypes per individual and SNP: `sample(c("0/1", "0/0", "1/1"), prob = c(0.6, 0.2, 0.2))`. For every non-null, non-bias gene, redraw that gene's genotypes of all individuals until at least 3 individuals have at least one `0/1` SNP in it.
   - Phase per individual and gene: `sample(c(1, -1), 1)`. For `diff_phase` genes it is fixed to `c(1, 1, -1, -1)` for `ind1`..`ind4`.
   - ALT fraction per cell: if `consistent`, `p_alt = f`; else if `phase == 1`, `p_alt = f`; else `p_alt = 1 - f`. For genes where the individual has no het SNP, `p_alt = 0.5`.
   - Per sample, `p_sg` is drawn with `PHI_BIO` around `p_alt`. Haplotypes, VCFs and `truth_individual_genes.tsv` follow Stage 1's outbred block. Direct counts are written only at the individual's `0/1` SNPs.
   - Also write `.unfiltered.table` as a copy and `.wasp_stats.tsv` as in the interface.
9. **Samples sheets:** write both `samples.csv` files with absolute FASTQ paths (`normalizePath(outdir)`).

- [ ] **Step 3: Write the FASTQ-level biology check**

Create `check_biology_stage2.R` by copying `check_biology.R` and changing only these parts:
- `D` points at `<outdir>/f1_recip` and reads `../genome/genome.fa` and `../genome/genome.gtf`.
- `ts` is `truth_snps.tsv`, with the window built from `tx_off` directly (`off <- s$tx_off`).
- The "planted" value per gene is the realised REF-fragment fraction pooled over all 12 samples: `sum(frag_A) / sum(frag_A + frag_B)` from `truth_fragments.tsv`.
- The pass rule is `abs(observed_alt_fraction - (1 - realised_A_fraction)) <= 0.05` for every gene.

It prints `BIOLOGY PASS` or `BIOLOGY FAIL` and exits 0 or 1.

- [ ] **Step 4: Write the injected-defect proofs**

Create `ase-pipeline/tests/synthetic/prove_checks_stage2.sh`:

```bash
#!/bin/bash
# Usage: prove_checks_stage2.sh <good outdir> <scratch dir>: the check must PASS on the good data and FAIL on each injected defect
set -u; G=${1:?}; W=${2:?}; HERE=$(dirname "$0"); fail=0
bash $HERE/check_simulation_stage2.sh $G >/dev/null || { echo "FAIL: check does not pass on the good data"; exit 1; }
defect() {   # $1 = name, $2 = shell command run inside the copy ($C)
  rm -rf $W/copy; cp -r $G $W/copy; C=$W/copy; eval "$2"
  if bash $HERE/check_simulation_stage2.sh $C >/dev/null 2>&1; then echo "FAIL: defect '$1' not detected"; fail=1; else echo "ok   defect '$1' detected"; fi
}
defect "AxB sample relabelled BxA"     'sed -i "0,/,AxB,/s//,BxA,/" $C/f1_recip/samples.csv'
defect "maternal b1 set to 0"          'awk -F"\t" "BEGIN{OFS=FS} NR>1 && \$3==\"maternal\" && !x {\$5=0; x=1} 1" $C/f1_recip/truth_genes.tsv > $C/t && mv $C/t $C/f1_recip/truth_genes.tsv'
defect "REF/ALT swapped in parental VCF" 'awk -F"\t" "BEGIN{OFS=FS} !/^#/ && !x {t=\$4; \$4=\$5; \$5=t; x=1} 1" $C/f1_recip/parental_snps.vcf > $C/t && mv $C/t $C/f1_recip/parental_snps.vcf'
defect "one read pair dropped from mate 1" 's=$(tail -n +2 $C/f1_recip/samples.csv | head -1 | cut -d, -f1); zcat $C/f1_recip/${s}_1.fastq.gz | tail -n +5 | gzip > $C/t && mv $C/t $C/f1_recip/${s}_1.fastq.gz'
defect "direct count ref+alt != total"  's=$(tail -n +2 $C/f1_recip/samples.csv | head -1 | cut -d, -f1); awk -F"\t" "BEGIN{OFS=FS} NR==2{\$6=\$6+5} 1" $C/f1_recip/direct_counts/$s.table > $C/t && mv $C/t $C/f1_recip/direct_counts/$s.table'
defect "outbred individual unpaired"   'sed -i "0,/,ind1$/s//,ind9/" $C/outbred_diff/samples.csv'
defect "cell fraction inconsistent"     'awk -F"\t" "BEGIN{OFS=FS} NR>1 && \$3==\"strain_A_high\" && !x {\$7=0.5; x=1} 1" $C/f1_recip/truth_genes.tsv > $C/t && mv $C/t $C/f1_recip/truth_genes.tsv'
defect "diff_phase phases not mixed"   'awk -F"\t" "BEGIN{OFS=FS} NR==FNR{ if (\$3==\"diff_phase\") g[\$1]=1; next } FNR>1 && (\$2 in g) {\$3=1} 1" $C/outbred_diff/truth_genes.tsv $C/outbred_diff/truth_individual_genes.tsv > $C/t && mv $C/t $C/outbred_diff/truth_individual_genes.tsv'
defect "direct F1 fractions altered"    'for f in $C/f1_recip/direct_counts/*.table; do awk -F"\t" "BEGIN{OFS=FS} NR>1{\$6=int(\$8*0.9); \$7=\$8-\$6} 1" $f > $C/t && mv $C/t $f; done'
defect "seed.txt removed"              'rm -f $C/seed.txt'
rm -rf $W/copy
[ $fail -eq 0 ] && echo "PROOFS PASS" || exit 1
```
The "maternal b1 set to 0" defect breaks both the value quota and the cell consistency, so either rule catches it. The "cell fraction inconsistent" defect breaks only the cell rule.

- [ ] **Step 5: Generate, check, prove, and test determinism on the cluster**

```bash
cat > $SCR/t2_gen.sh <<'EOF'
#!/bin/bash
#SBATCH -p bcc
#SBATCH -n 1 --mem=8G -t 1:00:00
set -uo pipefail
module add singularity/3.10.4 || exit 1
SIF=/net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif; T=ase-pipeline/tests/synthetic
cd /net/bmc-lab3/data/bcc/yannvrb/Github_repos/bioinformatics-claude-skills || exit 1
singularity exec --bind /net/bmc-lab3 $SIF Rscript $T/simulate_ase_stage2.R /net/bmc-lab3/data/bcc/yannvrb/ase_synthetic_stage2 20260930 || exit 1
singularity exec --bind /net/bmc-lab3 $SIF Rscript $T/simulate_ase_stage2.R /net/bmc-lab3/data/bcc/yannvrb/ase_stage2_scratch/regen 20260930 || exit 1
singularity exec --bind /net/bmc-lab3 $SIF Rscript $T/check_biology_stage2.R /net/bmc-lab3/data/bcc/yannvrb/ase_synthetic_stage2 || exit 1
EOF
sbatch -p bcc -o $SCR/logs/t2_gen_%j.out $SCR/t2_gen.sh
```
After the job ends:
- `bash ase-pipeline/tests/synthetic/check_simulation_stage2.sh /net/bmc-lab3/data/bcc/yannvrb/ase_synthetic_stage2` prints `PASS`.
- `bash ase-pipeline/tests/synthetic/prove_checks_stage2.sh /net/bmc-lab3/data/bcc/yannvrb/ase_synthetic_stage2 $SCR/proofs` prints 10 `ok` lines and `PROOFS PASS`.
- `diff <(cd /net/bmc-lab3/data/bcc/yannvrb/ase_synthetic_stage2 && find . -type f ! -name samples.csv -exec md5sum {} + | sort -k2) <(cd $SCR/regen && find . -type f ! -name samples.csv -exec md5sum {} + | sort -k2)` prints nothing. The samples sheets differ only by their absolute paths.
- The biology log ends with `BIOLOGY PASS`.

The two checks are shell (zcat/awk over about 50 MB of FASTQ). If they take more than a few seconds on the login node, run them inside the same kind of `sbatch` job instead.

- [ ] **Step 6: Commit (scripts only, never the data)**

```bash
git status && git diff --stat
git add ase-pipeline/tests/synthetic/simulate_ase_stage2.R ase-pipeline/tests/synthetic/check_simulation_stage2.sh ase-pipeline/tests/synthetic/check_biology_stage2.R ase-pipeline/tests/synthetic/prove_checks_stage2.sh
git commit -m "ase-pipeline: Stage 2 synthetic ground truth (reciprocal x condition F1, paired outbred) with exact checks and defect proofs

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Rmd 03 — reciprocal F1 (Step 16), run script, render tool and evaluator

**Files:**
- Modify: `ase-pipeline/ase-pipeline.md`: new `## Step 16 — Rmd 03: reciprocal F1 (strain and parent-of-origin effects)`, inserted between the end of Step 15 (the `---` line before `## Notes for the assistant`) and the Notes.
- Create: `ase-pipeline/tests/synthetic/render_from_skill.sh`, `ase-pipeline/tests/synthetic/evaluate_stage2.R`
- Modify: `ase-pipeline/tests/check_skill.sh`

**Interfaces:**
- Consumes:
  - from Task 1: `chrom_class`, `thin_snps`, `ase_glm_test` and `THIN_BP`;
  - from Rmd 01: `ase_checkpoint.rds` (`sites` with `sample, contig, position, SNP, ref_n, alt_n, total, condition, cross_direction, individual`; `samples`; `constants`; `bias`);
  - from Task 2: `f1_recip/direct_counts`, `truth_genes.tsv`.
- Produces:
  - `{CWD}/{TODAY}_{WD_NAME}_03_reciprocal.Rmd`, which renders to `{RESULTS_DIR}/{TODAY}_{WD_NAME}_03_reciprocal.html`;
  - `{RESULTS_DIR}/{TODAY}_{WD_NAME}_ASE_reciprocal.xlsx` (sheets `Gene`, `SNPs_used`, `Design`, `Excluded`, `Summary`) and `{TODAY}_{WD_NAME}_ASE_reciprocal.pdf`;
  - `{RESULTS_DIR}/ase_reciprocal_checkpoint.rds` = `list(gene, not_tested, snps_used, design, excluded_chrom, summary, constants)`. `gene` has the columns of `ase_glm_test` plus `gene_id, n_snps_used, n_samples_AxB, n_samples_BxA, frac_A_AxB, frac_A_BxA, b0_strain, b1_parent_of_origin, frac_A, maternal_frac, padj_strain, padj_parent_of_origin, sig_strain, sig_parent_of_origin, direction_strain, direction_parent_of_origin, category`;
  - `{RESULTS_DIR}/summary_numbers_reciprocal.tsv` with columns `test, genes_tested, genes_not_tested, sig_genes, direction_pos, n_pos, direction_neg, n_neg, phi_common, samples_AxB, samples_BxA` (rows `strain`, `parent_of_origin`);
  - `{RESULTS_DIR}/scripts/run_03_reciprocal.sh`;
  - `render_from_skill.sh <skill.md> <step number> <values.tsv> <out.Rmd>`;
  - `Rscript evaluate_stage2.R reciprocal <RESULTS_DIR> <truth dir>`, which prints `PASS`/`FAIL` lines and exits 0/1.

- [ ] **Step 1: Failing checker lines**

Append to `check_skill.sh`:

```bash
# --- Stage 2 Task 3 (Rmd 03, Step 16)
need "## Step 16 — Rmd 03: reciprocal F1"
need "_03_reciprocal.Rmd"
need "_ASE_reciprocal.xlsx"
need "ase_reciprocal_checkpoint.rds"
need "summary_numbers_reciprocal.tsv"
need "run_03_reciprocal.sh"
need 'THIN_BP     <- {THIN_BP}'
need 'dirs$d <- ifelse(dirs$cross_direction == "AxB", 1, -1)'
need 'ase_glm_test(units, list(strain = "(Intercept)", parent_of_origin = "d"), RHO_MIN)'
need "cross_direction must be AxB"
need 'dplyr::mutate(keep = thin_snps(exon_pos, depth, THIN_BP))'
need 'chrom == "autosome"'
need "b1 is not affected by a constant mapping bias"
# --- end Stage 2 Task 3
```
Run the checker. Expected: FAIL (13 lines).

- [ ] **Step 2: Write Step 16 into the skill**

Insert the following as the new Step 16 (between Step 15's closing `---` and `## Notes for the assistant`, followed by its own `---`):

`````markdown
## Step 16 — Rmd 03: reciprocal F1 (strain and parent-of-origin effects)

Only when "Reciprocal F1 analysis" was selected in Step 7 (`{MODE}` = `f1`). Write `{CWD}/{TODAY}_{WD_NAME}_03_reciprocal.Rmd` with the same conventions as Step 12, and the same `{AUTHOR}` and `{PROJECT_TITLE}`. It loads `ase_checkpoint.rds` (Rmd 01), checks the constants against it, pastes the statistics functions of Step 14 into the chunk marked below, and writes these files to `{RESULTS_DIR}`: `{TODAY}_{WD_NAME}_ASE_reciprocal.xlsx` (sheets `Gene`, `SNPs_used`, `Design`, `Excluded`, `Summary`), `{TODAY}_{WD_NAME}_ASE_reciprocal.pdf`, `ase_reciprocal_checkpoint.rds` and `summary_numbers_reciprocal.tsv`.

- **Model.** Per gene, `logit(p) = b0 + b1 * d`, where p is the **strain-A (`{STRAIN_A}`, REF) fraction** and d = +1 for `AxB` (strain A is the mother) and −1 for `BxA`. When several conditions exist, sum-to-zero condition terms are added, so that b0 and b1 are averages over the conditions. **b0 is the strain effect** (cis-regulatory divergence): b0 > 0 means "`{STRAIN_A}` higher". **b1 is the parent-of-origin effect** (imprinting): b1 > 0 means "maternal higher", and `maternal_frac = plogis(b1)`.
- **Rows.** One row per sample and gene: the strain-A and total counts summed over the gene's SNPs after `thin_snps` (one SNP per `THIN_BP` window in exon coordinates, deepest SNP first, the same SNPs in every sample). ASEReadCounter counts a read pair at every SNP it covers, so the gene test uses no SNP pair that one read pair can span. A gene with a single SNP is tested normally, because the replication comes from the samples. SNPs that lie in two genes are left out.
- **Dispersion and tests.** `ase_glm_test` (Step 14): the pooled dispersion between replicate animals (`phi_common`), each gene's own estimate when it is larger and the gene has at least 2 residual df, floored at `RHO_MIN`. b0 and b1 are tested by likelihood-ratio tests (1 df each), with BH within each test. `sig_strain` requires `padj_strain < FDR_SIG` and `|frac_A - 0.5| >= ABS_DEV_SIG`; `sig_parent_of_origin` requires `padj_parent_of_origin < FDR_SIG` and `|maternal_frac - 0.5| >= ABS_DEV_SIG`. Monoallelic genes (complete imprinting) have their coefficients at the ±15 bound (`status` `ok_at_bound`) with a valid p-value.
- **Requirements.** Every sample's `cross_direction` must be exactly `AxB` or `BxA`, and each direction needs at least 2 samples; otherwise the Rmd stops. Per gene, at least 2 samples with coverage in each direction; genes below that are listed in `not_tested`.
- **Excluded chromosomes.** X, Y and MT (and `chrX`, `chrY`, `chrM`) are not tested. A hemizygous X in males and the maternally inherited MT would look like maternal imprinting. They are counted in the `Excluded` sheet.
- **Reference bias.** b0 is exactly what a mapping bias toward strain A would imitate, so the Rmd 01 bias flags are printed first and repeated next to the strain results. b1 is not affected by a constant mapping bias toward one strain (it cancels between the two directions).

````rmd
---
title: "{PROJECT_TITLE} - ASE reciprocal F1"
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
library(GenomicRanges); library(rtracklayer); library(openxlsx)   # Bioconductor first
library(tidyverse)
```

## Constants and checkpoint

```{r constants}
MODE        <- "{MODE}"
STRAIN_A    <- "{STRAIN_A}"
STRAIN_B    <- "{STRAIN_B}"
RESULTS_DIR <- "{RESULTS_DIR}"
GTF_PATH    <- "{GTF_PATH}"
DATE_TAG    <- "{TODAY}_{WD_NAME}"
MIN_DEPTH   <- {MIN_DEPTH}
FDR_SIG     <- {FDR_SIG}
ABS_DEV_SIG <- {ABS_DEV_SIG}
RHO_MIN     <- {RHO_MIN}             # floor for the dispersion used by the tests
BIAS_TOL    <- {BIAS_TOL}
THIN_BP     <- {THIN_BP}             # SNPs closer than this (exon coordinates) can share a read pair: one kept per window
if (MODE != "f1") stop("Rmd 03 (reciprocal F1) needs MODE = f1", call. = FALSE)
ck <- readRDS(file.path(RESULTS_DIR, "ase_checkpoint.rds"))
same <- c(MODE = identical(ck$constants$MODE, MODE),
          vapply(c("MIN_DEPTH", "FDR_SIG", "ABS_DEV_SIG", "BIAS_TOL"),
                 function(k) isTRUE(all.equal(ck$constants[[k]], get(k))), logical(1)))
if (!all(same)) stop("constants differ from the Rmd 01 checkpoint (", paste(names(same)[!same], collapse = ", "),
                     "); re-render Rmd 01 with the same values", call. = FALSE)
sites <- ck$sites; bias <- ck$bias
knitr::kable(bias, digits = 4, caption = paste0("Reference-bias diagnostic from Rmd 01 (REF = ", STRAIN_A, "); a screen, not a test. ",
                                                 "The strain effect b0 is sensitive to it; the parent-of-origin effect b1 is not"))
if (any(bias$flagged)) cat("FLAGGED samples:", paste(bias$sample[bias$flagged], collapse = ", "), "- read the strain effects with caution\n")
```

## Statistics functions

```{r stats}
# <<< paste here the code of Step 14 between the two marker lines (the marker lines themselves are not pasted) >>>
```

## Design

```{r design}
dirs <- unique(ck$samples[, c("sample", "condition", "cross_direction")])
bad <- dirs$sample[!dirs$cross_direction %in% c("AxB", "BxA")]
if (length(bad) > 0) stop("cross_direction must be AxB (", STRAIN_A, " mother) or BxA (", STRAIN_B, " mother); other values for: ",
                          paste(bad, collapse = ", "), call. = FALSE)
n_dir <- table(factor(dirs$cross_direction, levels = c("AxB", "BxA")))
if (any(n_dir < 2)) stop("the reciprocal analysis needs at least 2 samples in each cross direction (AxB: ", n_dir[["AxB"]],
                         ", BxA: ", n_dir[["BxA"]], ")", call. = FALSE)
dirs$d <- ifelse(dirs$cross_direction == "AxB", 1, -1)   # d = +1: strain A is the mother (maternal allele = STRAIN_A)
knitr::kable(dirs, caption = "Samples, condition and cross direction (d = +1: AxB, strain A maternal)")
```

## Gene map in exon coordinates and SNP thinning

```{r genemap}
gtf <- rtracklayer::import(GTF_PATH)
if (is.null(gtf$gene_id)) stop("the GTF has no gene_id attribute", call. = FALSE)
ex <- gtf[gtf$type == "exon" & !is.na(gtf$gene_id)]
exd <- as.data.frame(GenomicRanges::reduce(split(ex, ex$gene_id)))              # exon union per gene: group_name = gene_id
exd <- exd[order(exd$group_name, exd$start), ]
exd$offset <- ave(exd$width, exd$group_name, FUN = function(w) cumsum(w) - w)    # exonic bases of the gene before this exon
pos <- unique(sites[, c("contig", "position")])
gr <- GenomicRanges::GRanges(pos$contig, IRanges::IRanges(pos$position, width = 1))
exr <- GenomicRanges::GRanges(as.character(exd$seqnames), IRanges::IRanges(exd$start, exd$end))
hits <- GenomicRanges::findOverlaps(gr, exr); qh <- S4Vectors::queryHits(hits); sh <- S4Vectors::subjectHits(hits)
if (length(qh) == 0) stop("no SNP lies in an exon of the GTF (contig names: counts ", paste(head(unique(pos$contig), 3), collapse = ", "),
                          "; GTF ", paste(head(unique(as.character(exd$seqnames)), 3), collapse = ", "), ")", call. = FALSE)
snp_gene <- data.frame(contig = pos$contig[qh], position = pos$position[qh], gene_id = exd$group_name[sh],
                       exon_pos = exd$offset[sh] + pos$position[qh] - exd$start[sh] + 1, stringsAsFactors = FALSE)
ng <- table(paste(snp_gene$contig, snp_gene$position))
snp_gene$multi_gene <- as.vector(ng[paste(snp_gene$contig, snp_gene$position)]) > 1
use <- sites %>% dplyr::inner_join(dplyr::filter(snp_gene, !multi_gene), by = c("contig", "position")) %>%
  dplyr::mutate(chrom = chrom_class(contig))
excluded_chrom <- dplyr::count(dplyr::filter(use, chrom != "autosome"), chrom, name = "sample_site_rows")
use <- dplyr::filter(use, chrom == "autosome")
kept <- use %>% dplyr::group_by(gene_id, contig, position, exon_pos) %>% dplyr::summarise(depth = sum(total), .groups = "drop") %>%
  dplyr::group_by(gene_id) %>% dplyr::mutate(keep = thin_snps(exon_pos, depth, THIN_BP)) %>% dplyr::ungroup() %>% dplyr::filter(keep)
gene_rows <- use %>% dplyr::semi_join(kept, by = c("gene_id", "contig", "position")) %>%
  dplyr::group_by(gene_id, sample, condition, cross_direction) %>%
  dplyr::summarise(y = sum(ref_n), n = sum(total), .groups = "drop") %>% dplyr::filter(n > 0)
cat(sum(snp_gene$multi_gene), "SNP-gene pairs in overlapping genes left out;", nrow(kept), "SNPs kept after thinning in",
    length(unique(kept$gene_id)), "genes\n")
if (nrow(excluded_chrom) > 0) knitr::kable(excluded_chrom, caption = "Rows on X, Y or MT: not tested (hemizygous X and maternal MT mimic imprinting)")
```

## Per-gene model and tests

```{r test}
build <- lapply(split(gene_rows, gene_rows$gene_id), function(r) {
  nA <- sum(r$cross_direction == "AxB"); nB <- sum(r$cross_direction == "BxA")
  if (nA < 2 || nB < 2) return(list(unit = NULL, why = sprintf("fewer than 2 samples with coverage in a direction (AxB %d, BxA %d)", nA, nB)))
  X <- cbind("(Intercept)" = 1, d = ifelse(r$cross_direction == "AxB", 1, -1))
  cl <- sort(unique(r$condition))
  if (length(cl) > 1) {   # condition adjustment, sum-to-zero coding: b0 and b1 are averages over the conditions
    C <- contr.sum(length(cl)); colnames(C) <- paste0("cond_", cl[-length(cl)])
    X <- cbind(X, C[match(r$condition, cl), , drop = FALSE])
  }
  list(unit = list(y = r$y, n = r$n, X = X), why = NA_character_)
})
units <- lapply(Filter(function(b) !is.null(b$unit), build), `[[`, "unit")
not_tested <- data.frame(gene_id = names(build), reason = vapply(build, `[[`, character(1), "why"), stringsAsFactors = FALSE)
not_tested <- not_tested[!is.na(not_tested$reason), , drop = FALSE]
if (length(units) == 0) stop("no gene has at least 2 samples with coverage in each cross direction", call. = FALSE)
res <- ase_glm_test(units, list(strain = "(Intercept)", parent_of_origin = "d"), RHO_MIN)
obs <- gene_rows %>% dplyr::group_by(gene_id) %>%
  dplyr::summarise(n_samples_AxB = sum(cross_direction == "AxB"), n_samples_BxA = sum(cross_direction == "BxA"),
                   frac_A_AxB = sum(y[cross_direction == "AxB"]) / sum(n[cross_direction == "AxB"]),
                   frac_A_BxA = sum(y[cross_direction == "BxA"]) / sum(n[cross_direction == "BxA"]), .groups = "drop")
gene <- res %>% dplyr::rename(gene_id = unit) %>%
  dplyr::left_join(dplyr::count(kept, gene_id, name = "n_snps_used"), by = "gene_id") %>%
  dplyr::left_join(obs, by = "gene_id") %>%
  dplyr::mutate(b0_strain = `beta_(Intercept)`, b1_parent_of_origin = beta_d,
                frac_A = plogis(b0_strain), maternal_frac = plogis(b1_parent_of_origin),
                padj_strain = p.adjust(p_strain, "BH"), padj_parent_of_origin = p.adjust(p_parent_of_origin, "BH"),
                sig_strain = !is.na(padj_strain) & padj_strain < FDR_SIG & abs(frac_A - 0.5) >= ABS_DEV_SIG,
                sig_parent_of_origin = !is.na(padj_parent_of_origin) & padj_parent_of_origin < FDR_SIG & abs(maternal_frac - 0.5) >= ABS_DEV_SIG,
                direction_strain = dplyr::case_when(!sig_strain ~ "none", b0_strain > 0 ~ paste(STRAIN_A, "higher"), TRUE ~ paste(STRAIN_B, "higher")),
                direction_parent_of_origin = dplyr::case_when(!sig_parent_of_origin ~ "none", b1_parent_of_origin > 0 ~ "maternal higher",
                                                              TRUE ~ "paternal higher"),
                category = dplyr::case_when(sig_strain & sig_parent_of_origin ~ "strain and parent-of-origin", sig_strain ~ "strain (cis)",
                                            sig_parent_of_origin ~ "parent-of-origin", TRUE ~ "none")) %>%
  dplyr::arrange(gene_id)
cat("Dispersion between replicate animals (phi_common):", signif(gene$phi_common[1], 4), "; floor RHO_MIN =", RHO_MIN, "\n")
knitr::kable(dplyr::count(gene, status), caption = "Fit status (ok_at_bound = monoallelic gene, coefficient at the bound, p valid)")
knitr::kable(head(dplyr::filter(gene, category != "none"), 40), digits = 4,
             caption = paste0("Significant genes (first 40): strain effect = ", STRAIN_A, " fraction averaged over directions; ",
                              "parent-of-origin effect = maternal fraction", if (any(bias$flagged)) "; some samples flagged for reference bias" else ""))
```

## Summary

```{r summary}
row_for <- function(test, sig, b, pos_lab, neg_lab) {
  p <- gene[[paste0("p_", test)]]
  data.frame(test = test, genes_tested = sum(!is.na(p)), genes_not_tested = nrow(not_tested) + sum(is.na(p)), sig_genes = sum(sig),
             direction_pos = pos_lab, n_pos = sum(sig & b > 0, na.rm = TRUE), direction_neg = neg_lab, n_neg = sum(sig & b < 0, na.rm = TRUE),
             phi_common = gene$phi_common[1], samples_AxB = n_dir[["AxB"]], samples_BxA = n_dir[["BxA"]], stringsAsFactors = FALSE)
}
summary_rec <- rbind(row_for("strain", gene$sig_strain, gene$b0_strain, paste(STRAIN_A, "higher"), paste(STRAIN_B, "higher")),
                     row_for("parent_of_origin", gene$sig_parent_of_origin, gene$b1_parent_of_origin, "maternal higher", "paternal higher"))
knitr::kable(summary_rec, digits = 4, caption = "Summary (FDR_SIG, ABS_DEV_SIG as in the constants block)")
```

## Figure

```{r figure}
gene_plot <- dplyr::filter(gene, !is.na(p_strain))
stopifnot(sum(gene_plot$sig_strain) == summary_rec$sig_genes[1], sum(gene_plot$sig_parent_of_origin) == summary_rec$sig_genes[2])
cat("Figure sig counts equal Summary counts: TRUE\n")
p1 <- ggplot(gene_plot, aes(frac_A, maternal_frac, colour = category)) + geom_hline(yintercept = 0.5, linetype = 2) +
  geom_vline(xintercept = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.8) +
  labs(x = paste0(STRAIN_A, " fraction (strain effect, plogis(b0))"), y = "maternal fraction (parent-of-origin effect, plogis(b1))",
       colour = NULL, caption = paste0("sig: BH < ", FDR_SIG, " and |fraction - 0.5| >= ", ABS_DEV_SIG,
                                       if (any(bias$flagged)) "; some samples flagged for reference bias" else "")) + theme_bw()
print(p1)
ggsave(file.path(RESULTS_DIR, paste0(DATE_TAG, "_ASE_reciprocal.pdf")), p1, width = 9, height = 7)
```

## Export

```{r export}
out_names <- function(d) {
  names(d) <- sub("^frac_A$", paste0("frac_", STRAIN_A), names(d))
  names(d) <- sub("^frac_A_(AxB|BxA)$", paste0("frac_", STRAIN_A, "_\\1"), names(d))
  d
}
xlsx_file <- file.path(RESULTS_DIR, paste0(DATE_TAG, "_ASE_reciprocal.xlsx"))
wb <- openxlsx::createWorkbook()
for (nm in c("Gene", "SNPs_used", "Design", "Excluded", "Summary")) openxlsx::addWorksheet(wb, nm)
openxlsx::writeData(wb, "Gene", out_names(gene))
openxlsx::writeData(wb, "SNPs_used", kept)
openxlsx::writeData(wb, "Design", dirs)
openxlsx::writeData(wb, "Excluded", rbind(data.frame(what = sprintf("chromosome %s", excluded_chrom$chrom), n = excluded_chrom$sample_site_rows),
                                           data.frame(what = sprintf("not tested: %s", not_tested$reason), n = rep(1, nrow(not_tested)))) %>%
                                       dplyr::group_by(what) %>% dplyr::summarise(n = sum(n), .groups = "drop"))
openxlsx::writeData(wb, "Summary", summary_rec)
openxlsx::saveWorkbook(wb, xlsx_file, overwrite = TRUE)
saveRDS(list(gene = gene, not_tested = not_tested, snps_used = kept, design = dirs, excluded_chrom = excluded_chrom, summary = summary_rec,
             constants = c(ck$constants, list(RHO_MIN = RHO_MIN, THIN_BP = THIN_BP))),
        file.path(RESULTS_DIR, "ase_reciprocal_checkpoint.rds"))
write.table(summary_rec, file.path(RESULTS_DIR, "summary_numbers_reciprocal.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", xlsx_file, "and summary_numbers_reciprocal.tsv\n")
sessionInfo()
```
````

Render it with `{RESULTS_DIR}/scripts/run_03_reciprocal.sh`. This is the same script as `run_01_import_qc.sh` (Step 12), with job name `ase_03_reciprocal`, log `run_03_reciprocal_%j.out`, the Rmd 03 file name, the same `--bind` rule and `-n 1 --mem=16G -t 4:00:00`. The Rmd runs single-threaded. The time covers genome-wide gene sets (the Step 14 benchmark projects under 3 h for 60,000 genes). If the job stops with "cross_direction must be AxB", correct `{SAMPLES_CSV}` (Step 5), re-render Rmd 01, then Rmd 03. Submission order: Step 15.
`````

- [ ] **Step 3: Write the render tool and the evaluator**

Create `ase-pipeline/tests/synthetic/render_from_skill.sh`:

```bash
#!/bin/bash
# Usage: render_from_skill.sh <skill.md> <step number, e.g. 12|13|16|17> <values.tsv> <out.Rmd>
# Cuts the ````rmd block of "## Step <N> " out of the skill, splices the Step 14 statistics block (between the marker lines) into
# the "<<< paste here the code of Step 14" line, and substitutes {NAME} from values.tsv (NAME<TAB>value); stops on leftovers.
set -uo pipefail
SKILL=${1:?}; N=${2:?}; VALS=${3:?}; OUT=${4:?}
awk -v n="$N" '
  index($0, "## Step " n " ") == 1 { in_step = 1; next }
  in_step && /^## Step / { in_step = 0 }
  in_step && /^````rmd$/ { grab = 1; next }
  grab && /^````$/ { grab = 0; in_step = 0; next }
  grab { print }' "$SKILL" > "$OUT.tmp" || exit 1
[ -s "$OUT.tmp" ] || { echo "ERROR: no rmd block in Step $N" >&2; exit 1; }
awk '/^# --- ase-stats-begin/ {f = 1; next} /^# --- ase-stats-end/ {f = 0} f' "$SKILL" > "$OUT.stats" || exit 1
awk -v sf="$OUT.stats" '/<<< paste here the code of Step 14/ { while ((getline l < sf) > 0) print l; next } { print }' "$OUT.tmp" > "$OUT.tmp2" || exit 1
awk -F'\t' 'NR == FNR { v[$1] = $2; next } { for (k in v) { t = "{" k "}"; while ((i = index($0, t)) > 0) $0 = substr($0, 1, i - 1) v[k] substr($0, i + length(t)) } print }' \
  "$VALS" "$OUT.tmp2" > "$OUT" || exit 1
rm -f "$OUT.tmp" "$OUT.tmp2" "$OUT.stats"
if grep -nE '\{[A-Z][A-Z_0-9]*\}' "$OUT"; then echo "ERROR: unsubstituted placeholders above" >&2; exit 1; fi
echo "wrote $OUT"
```

Create `ase-pipeline/tests/synthetic/evaluate_stage2.R` (the Task 4 differential parts are added in Task 4):

```r
# Usage: Rscript evaluate_stage2.R <reciprocal|differential_f1|differential_outbred> <RESULTS_DIR> <truth dir (f1_recip or outbred_diff)>
args <- commandArgs(trailingOnly = TRUE); what <- args[1]; RES <- args[2]; TRUTH <- args[3]
fail <- 0
crit <- function(cond, msg) { cat(if (isTRUE(cond)) "PASS" else "FAIL", msg, "\n"); if (!isTRUE(cond)) fail <<- 1 }
tg <- read.delim(file.path(TRUTH, "truth_genes.tsv"), stringsAsFactors = FALSE)
if (what == "reciprocal") {
  g <- readRDS(file.path(RES, "ase_reciprocal_checkpoint.rds"))$gene
  m <- merge(tg, g, by = "gene_id", all.x = TRUE)
  for (k in c("sig_strain", "sig_parent_of_origin")) m[[k]][is.na(m[[k]])] <- FALSE
  print(m[m$class != "null", c("gene_id", "class", "expect_strain", "expect_poo", "frac_A", "maternal_frac", "padj_strain",
                               "padj_parent_of_origin", "direction_strain", "direction_parent_of_origin", "status")], row.names = FALSE)
  exS <- m$expect_strain %in% c("A_higher", "B_higher")
  okS <- m$sig_strain & ((m$expect_strain == "A_higher" & m$b0_strain > 0) | (m$expect_strain == "B_higher" & m$b0_strain < 0))
  crit(all(okS[exS]), sprintf("strain effect recovered with the correct sign: %d of %d", sum(okS[exS]), sum(exS)))
  exP <- m$expect_poo %in% c("maternal", "paternal")
  okP <- m$sig_parent_of_origin & ((m$expect_poo == "maternal" & m$b1_parent_of_origin > 0) | (m$expect_poo == "paternal" & m$b1_parent_of_origin < 0))
  crit(all(okP[exP]), sprintf("parent-of-origin effect recovered with the correct sign: %d of %d", sum(okP[exP]), sum(exP)))
  so <- m$class %in% c("strain_A_high", "strain_B_high"); po <- m$class %in% c("maternal", "paternal")
  crit(sum(m$sig_parent_of_origin[so]) == 0, sprintf("no parent-of-origin call among strain-only genes: %d of %d", sum(m$sig_parent_of_origin[so]), sum(so)))
  crit(sum(m$sig_strain[po]) == 0, sprintf("no strain call among parent-of-origin-only genes: %d of %d", sum(m$sig_strain[po]), sum(po)))
  nul <- m$class %in% c("null", "dense_null", "bias")
  crit(sum(m$sig_strain[nul]) <= 2 && sum(m$sig_parent_of_origin[nul]) <= 2,
       sprintf("null genes: strain false positives %d, parent-of-origin false positives %d, of %d (limit 2 each)",
               sum(m$sig_strain[nul]), sum(m$sig_parent_of_origin[nul]), sum(nul)))
  crit(all(m$status[nul | exS | exP] %in% c("ok", "ok_at_bound")), "every planted and null gene was tested")
  su <- readRDS(file.path(RES, "ase_reciprocal_checkpoint.rds"))$snps_used
  dn <- su[su$gene_id %in% m$gene_id[m$class == "dense_null"], ]
  crit(all(tapply(dn$exon_pos, dn$gene_id, function(p) length(p) < 2 || min(diff(sort(p))) >= 500)),
       "dense_null genes: kept SNPs at least THIN_BP (500) apart")
}
quit(status = fail)
```

- [ ] **Step 4: Run the checker. Expected: PASS**

- [ ] **Step 5: Render on the direct counts (F1) and evaluate against truth**

```bash
P=$SCR/t3_f1; RD=$P/results/2026-09-30_t3_f1; mkdir -p $RD/ase_counts $RD/logs
cp /net/bmc-lab3/data/bcc/yannvrb/ase_synthetic_stage2/f1_recip/direct_counts/*.table $RD/ase_counts/
cp /net/bmc-lab3/data/bcc/yannvrb/ase_synthetic_stage2/f1_recip/samples.csv $P/t3_f1_samples.csv
cp -r /net/bmc-lab3/data/bcc/yannvrb/ase_synthetic_stage2/genome $P/genome
printf '%s\t%s\n' MODE f1 STRAIN_A C57BL_6NJ STRAIN_B A_J SAMPLES_CSV $P/t3_f1_samples.csv RESULTS_DIR $RD GTF_PATH $P/genome/genome.gtf \
  MIN_DEPTH 10 FDR_SIG 0.05 ABS_DEV_SIG 0.1 BIAS_TOL 0.03 RHO_MIN 0.01 THIN_BP 500 TODAY 2026-09-30 WD_NAME t3_f1 \
  PROJECT_TITLE "Stage 2 test" AUTHOR tester REF_CONDITION ctrl > $P/values.tsv
T=ase-pipeline/tests/synthetic
bash $T/render_from_skill.sh ase-pipeline/ase-pipeline.md 12 $P/values.tsv $P/2026-09-30_t3_f1_01_import_qc.Rmd
bash $T/render_from_skill.sh ase-pipeline/ase-pipeline.md 16 $P/values.tsv $P/2026-09-30_t3_f1_03_reciprocal.Rmd
cat > $P/render.sh <<EOF
#!/bin/bash
#SBATCH -p bcc
#SBATCH -n 1 --mem=16G -t 1:00:00
set -uo pipefail
module add singularity/3.10.4 || exit 1
R() { singularity exec --bind /net/bmc-lab3 /net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif "\$@"; }
R Rscript -e "rmarkdown::render('$P/2026-09-30_t3_f1_01_import_qc.Rmd', output_dir = '$RD', knit_root_dir = '$P')" || exit 1
R Rscript -e "rmarkdown::render('$P/2026-09-30_t3_f1_03_reciprocal.Rmd', output_dir = '$RD', knit_root_dir = '$P')" || exit 1
R Rscript /net/bmc-lab3/data/bcc/yannvrb/Github_repos/bioinformatics-claude-skills/$T/evaluate_stage2.R reciprocal $RD /net/bmc-lab3/data/bcc/yannvrb/ase_synthetic_stage2/f1_recip
EOF
sbatch -p bcc -o $P/render_%j.out $P/render.sh
```
Expected: both renders exit 0, `summary_numbers_reciprocal.tsv` has 2 data rows, the xlsx has 5 sheets, and every evaluator line is `PASS`. The job exits 0. This is the same list of criteria that Task 8 applies to the aligned data. If an evaluator line fails, first check the sign convention and the thinning in the Rmd, never the truth.

- [ ] **Step 6: Negative tests (Review Focus 1 and 3)**

1. **Free-text direction:** copy `$P` to `$SCR/t3_bad`. In the copied samples CSV, replace every `AxB` by `B6xAJ` (`sed -i 's/,AxB,/,B6xAJ,/'`). Update `values.tsv`, regenerate both Rmds, and render. Expected: Rmd 01 succeeds, and Rmd 03 stops with the message starting `cross_direction must be AxB`; the job log contains that text.
2. **X-linked gene:** copy `$P` to `$SCR/t3_chrx`. In the copied direct-count tables and the copied GTF, rename the contig of one `maternal` gene's SNPs and of its GTF lines to `chrX`. Use `awk` keyed on that gene's positions from `truth_snps.tsv` and its `gene_id` in the GTF. Render. Expected: the gene is absent from `gene`, and the `Excluded` sheet has a `chromosome X` row.

- [ ] **Step 7: Commit**

```bash
git status && git diff --stat
git add ase-pipeline/ase-pipeline.md ase-pipeline/tests/check_skill.sh ase-pipeline/tests/synthetic/render_from_skill.sh ase-pipeline/tests/synthetic/evaluate_stage2.R
git commit -m "ase-pipeline: Rmd 03 reciprocal F1 (Step 16), render tool and truth evaluator

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Rmd 04 — differential ASE (Step 17)

**Files:**
- Modify: `ase-pipeline/ase-pipeline.md`: new `## Step 17 — Rmd 04: differential ASE between conditions`, inserted after Step 16 and before `## Notes for the assistant`.
- Modify: `ase-pipeline/tests/synthetic/evaluate_stage2.R` (two new branches), `ase-pipeline/tests/check_skill.sh`

**Interfaces:**
- Consumes: Task 1 functions (`ase_glm_test`, `ase_paired_test`, `thin_snps`, `chrom_class`, `acat`), the Rmd 01 checkpoint, the genemap and thinning code of Step 16 (repeated verbatim, because the Rmds are self-contained), and the placeholder `{REF_CONDITION}` (set in Step 7 by Task 6; the test values file of Task 3 already contains it).
- Produces:
  - `{CWD}/{TODAY}_{WD_NAME}_04_differential.Rmd`, which renders to `{RESULTS_DIR}/{TODAY}_{WD_NAME}_04_differential.html`;
  - `{TODAY}_{WD_NAME}_ASE_differential.xlsx` (sheets `SNP`, `Gene`, `Design`, `Excluded`, `Summary`) and `{TODAY}_{WD_NAME}_ASE_differential.pdf`;
  - `ase_differential_checkpoint.rds` = `list(snp, gene, not_tested, design, excluded_chrom, snp_gene, dispersion, summary, constants)`;
  - `summary_numbers_differential.tsv` with columns `contrast, mode, level, tested, not_tested, sig, n_up, n_down, dispersion`;
  - `run_04_differential.sh`.

  Column contract:
  - F1 `gene` has `contrast, gene_id, b0, b_condition, frac_A_ref, frac_A_test, delta_frac, p_condition, padj, sig, direction, status, phi_common, phi_unit, phi_used, n_rows`. F1 `snp` has the same columns with `SNP` in place of `gene_id`.
  - Outbred `snp` has `contrast, SNP, n_individuals, n_failed, stat, df, p, padj, mean_delta, max_abs_delta, n_up, n_down, sig, direction, status`.
  - Outbred `gene` has `contrast, gene_id, n_snps, acat_p, padj, max_abs_delta, sig, note`.
  - The contrast label is `paste(L, "vs", REF_CONDITION)`, for example `treat vs ctrl`.

- [ ] **Step 1: Failing checker lines**

```bash
# --- Stage 2 Task 4 (Rmd 04, Step 17)
need "## Step 17 — Rmd 04: differential ASE between conditions"
need "_04_differential.Rmd"
need "_ASE_differential.xlsx"
need "ase_differential_checkpoint.rds"
need "summary_numbers_differential.tsv"
need "run_04_differential.sh"
need 'REF_CONDITION <- "{REF_CONDITION}"'
need "is confounded with the cross direction"
need 'ase_glm_test(gu$units, list(condition = "cond"), RHO_MIN)'
need "r <- ase_paired_test(as.data.frame(dp), RHO_MIN)"
need '"mixed (phase differs)"'
need "acat_p = acat(p), max_abs_delta = max(max_abs_delta)"
# --- end Stage 2 Task 4
```
Run the checker. Expected: 12 new FAIL lines.

- [ ] **Step 2: Write Step 17 into the skill**

Insert:

`````markdown
## Step 17 — Rmd 04: differential ASE between conditions

Only when "Differential ASE between conditions" was selected in Step 7. Write `{CWD}/{TODAY}_{WD_NAME}_04_differential.Rmd` with the conventions of Step 12. It writes these files to `{RESULTS_DIR}`: `{TODAY}_{WD_NAME}_ASE_differential.xlsx` (sheets `SNP`, `Gene`, `Design`, `Excluded`, `Summary`), `{TODAY}_{WD_NAME}_ASE_differential.pdf`, `ase_differential_checkpoint.rds` and `summary_numbers_differential.tsv`. Every condition other than `{REF_CONDITION}` (Step 7) is compared with it, one contrast at a time (1-df likelihood-ratio test), with BH within each contrast and level (SNP or gene).

- **F1.** Per gene (rows = samples, strain-A and total counts summed over the thinned SNPs exactly as in Rmd 03) and per SNP (rows = samples), `logit(p) = b0 + b_condition * cond`. Here p is the strain-A fraction and cond = 1 for the tested condition. When both cross directions are among the contrast's samples, the cross-direction term `d` is added, so that parent-of-origin genes are not mistaken for condition effects. If condition and direction are confounded (every sample of one condition has one direction and every sample of the other the other direction), the Rmd stops, because the two effects cannot be separated. The dispersion and the tests are those of Rmd 03 (`ase_glm_test`). `delta_frac = plogis(b0 + b_condition) - plogis(b0)`, and `sig` requires `padj < FDR_SIG` and `|delta_frac| >= ABS_DEV_SIG`.
- **Outbred.** Only individuals sampled in both conditions of a contrast are used (the others are listed as unpaired). `ase_paired_test` works per individual. First, a pair dispersion is estimated from all its SNPs (each SNP keeps its own REF fraction, floored at `RHO_MIN`). Then each SNP gets a likelihood-ratio test of one REF fraction against one per condition. The per-SNP statistic is the sum over the informative individuals (df = their number), and at least 2 individuals are needed. The test has no direction, because the allele that carries a regulatory variant differs between individuals: in one person the REF allele of a SNP can rise while in another it falls. `mean_delta` (REF fraction change), `n_up` and `n_down` describe the individual changes. A significant SNP is labelled "REF higher (or lower) in the tested condition (all individuals)" only when all individuals agree, and "mixed (phase differs)" otherwise. Genes use `acat` over their SNP p-values and are labelled "unphased, no direction". `sig` requires `padj < FDR_SIG` and `max_abs_delta >= ABS_DEV_SIG`.
- X, Y and MT are excluded in both modes, as in Rmd 03. The Rmd 01 bias flags are printed first.

````rmd
---
title: "{PROJECT_TITLE} - ASE differential between conditions"
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
library(GenomicRanges); library(rtracklayer); library(openxlsx)   # Bioconductor first
library(tidyverse)
```

## Constants and checkpoint

```{r constants}
MODE          <- "{MODE}"
STRAIN_A      <- "{STRAIN_A}"
STRAIN_B      <- "{STRAIN_B}"
RESULTS_DIR   <- "{RESULTS_DIR}"
GTF_PATH      <- "{GTF_PATH}"
DATE_TAG      <- "{TODAY}_{WD_NAME}"
MIN_DEPTH     <- {MIN_DEPTH}
FDR_SIG       <- {FDR_SIG}
ABS_DEV_SIG   <- {ABS_DEV_SIG}
RHO_MIN       <- {RHO_MIN}
BIAS_TOL      <- {BIAS_TOL}
THIN_BP       <- {THIN_BP}
REF_CONDITION <- "{REF_CONDITION}"
stopifnot(MODE %in% c("f1", "outbred"))
ck <- readRDS(file.path(RESULTS_DIR, "ase_checkpoint.rds"))
same <- c(MODE = identical(ck$constants$MODE, MODE),
          vapply(c("MIN_DEPTH", "FDR_SIG", "ABS_DEV_SIG", "BIAS_TOL"),
                 function(k) isTRUE(all.equal(ck$constants[[k]], get(k))), logical(1)))
if (!all(same)) stop("constants differ from the Rmd 01 checkpoint (", paste(names(same)[!same], collapse = ", "),
                     "); re-render Rmd 01 with the same values", call. = FALSE)
sites <- ck$sites; bias <- ck$bias
knitr::kable(bias, digits = 4, caption = paste0("Reference-bias diagnostic from Rmd 01 (REF = ", STRAIN_A, "); a screen, not a test"))
if (any(bias$flagged)) cat("FLAGGED samples:", paste(bias$sample[bias$flagged], collapse = ", "), "\n")
```

## Statistics functions

```{r stats}
# <<< paste here the code of Step 14 between the two marker lines (the marker lines themselves are not pasted) >>>
```

## Design and contrasts

```{r design}
smp <- unique(ck$samples[, c("sample", "condition", "cross_direction", "individual")])
if (!REF_CONDITION %in% smp$condition) stop("reference condition '", REF_CONDITION, "' is not a condition in samples.csv", call. = FALSE)
if (MODE == "f1") {
  bad <- smp$sample[!smp$cross_direction %in% c("AxB", "BxA")]
  if (length(bad) > 0) stop("cross_direction must be AxB (", STRAIN_A, " mother) or BxA (", STRAIN_B, " mother); other values for: ",
                            paste(bad, collapse = ", "), call. = FALSE)
}
levels_test <- setdiff(sort(unique(smp$condition)), REF_CONDITION)
if (length(levels_test) == 0) stop("differential ASE needs at least two conditions", call. = FALSE)
design <- list(); skipped <- character()
for (L in levels_test) {
  s <- smp[smp$condition %in% c(REF_CONDITION, L), ]
  if (MODE == "f1") {
    nc <- table(factor(s$condition, levels = c(REF_CONDITION, L)))
    if (any(nc < 2)) { skipped[L] <- sprintf("fewer than 2 samples in a condition (%s %d, %s %d)", REF_CONDITION, nc[[1]], L, nc[[2]]); next }
    use_d <- all(c("AxB", "BxA") %in% s$cross_direction)
    if (use_d && qr(cbind(1, s$condition == L, s$cross_direction == "AxB"))$rank < 3)
      stop("condition '", L, "' versus '", REF_CONDITION, "' is confounded with the cross direction (each condition has only one ",
           "direction), so a condition effect cannot be separated from a parent-of-origin effect. Add samples of both directions to a ",
           "condition, or restrict samples.csv to one direction and re-run Rmd 01", call. = FALSE)
    design[[L]] <- list(samples = s, use_d = use_d)
  } else {
    pairs <- intersect(s$individual[s$condition == REF_CONDITION], s$individual[s$condition == L])
    if (length(pairs) < 2) { skipped[L] <- sprintf("fewer than 2 individuals sampled in both %s and %s", REF_CONDITION, L); next }
    design[[L]] <- list(samples = s[s$individual %in% pairs, ], unpaired = sort(setdiff(s$individual, pairs)))
  }
}
if (length(skipped) > 0) cat("Contrasts not tested:", paste(names(skipped), skipped, sep = ": ", collapse = "; "), "\n")
if (length(design) == 0) stop("no contrast can be tested: ", paste(names(skipped), skipped, sep = ": ", collapse = "; "), call. = FALSE)
design_tbl <- dplyr::bind_rows(lapply(names(design), function(L) dplyr::mutate(design[[L]]$samples, contrast = paste(L, "vs", REF_CONDITION),
  cross_direction_term = if (MODE == "f1") design[[L]]$use_d else NA,
  unpaired_individuals = if (MODE == "outbred") paste(design[[L]]$unpaired, collapse = ";") else NA_character_)))
knitr::kable(design_tbl, caption = "Samples per contrast (outbred: paired individuals only; unpaired ones listed)")
```

## Gene map in exon coordinates and SNP thinning

```{r genemap}
gtf <- rtracklayer::import(GTF_PATH)
if (is.null(gtf$gene_id)) stop("the GTF has no gene_id attribute", call. = FALSE)
ex <- gtf[gtf$type == "exon" & !is.na(gtf$gene_id)]
exd <- as.data.frame(GenomicRanges::reduce(split(ex, ex$gene_id)))
exd <- exd[order(exd$group_name, exd$start), ]
exd$offset <- ave(exd$width, exd$group_name, FUN = function(w) cumsum(w) - w)
pos <- unique(sites[, c("contig", "position")])
gr <- GenomicRanges::GRanges(pos$contig, IRanges::IRanges(pos$position, width = 1))
exr <- GenomicRanges::GRanges(as.character(exd$seqnames), IRanges::IRanges(exd$start, exd$end))
hits <- GenomicRanges::findOverlaps(gr, exr); qh <- S4Vectors::queryHits(hits); sh <- S4Vectors::subjectHits(hits)
if (length(qh) == 0) stop("no SNP lies in an exon of the GTF (contig names differ?)", call. = FALSE)
snp_gene <- data.frame(contig = pos$contig[qh], position = pos$position[qh], gene_id = exd$group_name[sh],
                       exon_pos = exd$offset[sh] + pos$position[qh] - exd$start[sh] + 1, stringsAsFactors = FALSE)
ng <- table(paste(snp_gene$contig, snp_gene$position))
snp_gene$multi_gene <- as.vector(ng[paste(snp_gene$contig, snp_gene$position)]) > 1
auto <- dplyr::mutate(sites, chrom = chrom_class(contig))
excluded_chrom <- dplyr::count(dplyr::filter(auto, chrom != "autosome"), chrom, name = "sample_site_rows")
auto <- dplyr::filter(auto, chrom == "autosome")
use <- dplyr::inner_join(auto, dplyr::filter(snp_gene, !multi_gene), by = c("contig", "position"))
kept <- use %>% dplyr::group_by(gene_id, contig, position, exon_pos) %>% dplyr::summarise(depth = sum(total), .groups = "drop") %>%
  dplyr::group_by(gene_id) %>% dplyr::mutate(keep = thin_snps(exon_pos, depth, THIN_BP)) %>% dplyr::ungroup() %>% dplyr::filter(keep)
gene_rows <- use %>% dplyr::semi_join(kept, by = c("gene_id", "contig", "position")) %>%
  dplyr::group_by(gene_id, sample, condition, cross_direction) %>%
  dplyr::summarise(y = sum(ref_n), n = sum(total), .groups = "drop") %>% dplyr::filter(n > 0)
snp_to_gene <- dplyr::distinct(use, SNP, gene_id)
if (nrow(excluded_chrom) > 0) knitr::kable(excluded_chrom, caption = "Rows on X, Y or MT: not tested")
```

## F1: per-gene and per-SNP condition effect

```{r f1, eval = (MODE == "f1")}
f1_units <- function(rows, id, L, use_d) {
  u <- lapply(split(rows, rows[[id]]), function(r) {
    cc <- as.numeric(r$condition == L)
    if (sum(cc == 1) < 2 || sum(cc == 0) < 2) return(NULL)
    X <- cbind("(Intercept)" = 1, cond = cc)
    if (use_d && length(unique(r$cross_direction)) == 2) X <- cbind(X, d = ifelse(r$cross_direction == "AxB", 1, -1))
    list(y = r$y, n = r$n, X = X)
  })
  nul <- vapply(u, is.null, logical(1))
  list(units = u[!nul], not_tested = names(u)[nul])
}
f1_table <- function(res, L, idname) {
  res %>% dplyr::rename(!!idname := unit) %>%
    dplyr::mutate(contrast = paste(L, "vs", REF_CONDITION), b0 = `beta_(Intercept)`, b_condition = beta_cond,
                  frac_A_ref = plogis(b0), frac_A_test = plogis(b0 + b_condition), delta_frac = frac_A_test - frac_A_ref,
                  padj = p.adjust(p_condition, "BH"), sig = !is.na(padj) & padj < FDR_SIG & abs(delta_frac) >= ABS_DEV_SIG,
                  direction = dplyr::case_when(!sig ~ "none", b_condition > 0 ~ paste0(STRAIN_A, " fraction higher in ", L),
                                               TRUE ~ paste0(STRAIN_A, " fraction lower in ", L)))
}
snp_rows <- dplyr::transmute(auto, SNP, sample, condition, cross_direction, y = ref_n, n = total) %>% dplyr::filter(n > 0)
gene_l <- list(); snp_l <- list(); nt_l <- list()
for (L in names(design)) {
  ks <- design[[L]]$samples$sample; ct <- paste(L, "vs", REF_CONDITION)
  gu <- f1_units(dplyr::filter(gene_rows, sample %in% ks), "gene_id", L, design[[L]]$use_d)
  su <- f1_units(dplyr::filter(snp_rows, sample %in% ks), "SNP", L, design[[L]]$use_d)
  if (length(gu$units) > 0) gene_l[[L]] <- f1_table(ase_glm_test(gu$units, list(condition = "cond"), RHO_MIN), L, "gene_id")
  if (length(su$units) > 0) snp_l[[L]] <- f1_table(ase_glm_test(su$units, list(condition = "cond"), RHO_MIN), L, "SNP")
  nt_l[[L]] <- rbind(data.frame(contrast = rep(ct, length(gu$not_tested)), level = rep("gene", length(gu$not_tested)), id = gu$not_tested,
                                reason = rep("fewer than 2 samples with coverage in a condition", length(gu$not_tested))),
                     data.frame(contrast = rep(ct, length(su$not_tested)), level = rep("SNP", length(su$not_tested)), id = su$not_tested,
                                reason = rep("fewer than 2 samples with coverage in a condition", length(su$not_tested))))
}
gene <- dplyr::bind_rows(gene_l); snp <- dplyr::bind_rows(snp_l); not_tested <- dplyr::bind_rows(nt_l)
if (nrow(gene) == 0) stop("no gene could be tested in any contrast", call. = FALSE)
dispersion <- dplyr::distinct(gene, contrast, phi_common)
knitr::kable(head(dplyr::filter(gene, sig), 40), digits = 4, caption = paste0("Significant genes (first 40); delta_frac = change of the ", STRAIN_A, " fraction"))
```

## Outbred: paired per-SNP test and gene combination

```{r outbred, eval = (MODE == "outbred")}
snp_l <- list(); gene_l <- list(); phi_l <- list(); nt_l <- list()
for (L in names(design)) {
  s <- design[[L]]$samples; ct <- paste(L, "vs", REF_CONDITION)
  dp <- auto %>% dplyr::filter(sample %in% s$sample) %>%
    dplyr::transmute(snp = SNP, individual, cond = as.numeric(condition == L), y = ref_n, n = total)
  r <- ase_paired_test(as.data.frame(dp), RHO_MIN)
  st <- r$snp %>% dplyr::rename(SNP = snp) %>%
    dplyr::mutate(contrast = ct, padj = p.adjust(p, "BH"), sig = !is.na(padj) & padj < FDR_SIG & max_abs_delta >= ABS_DEV_SIG,
                  direction = dplyr::case_when(!sig ~ "none", n_down == 0 ~ paste0("REF higher in ", L, " (all individuals)"),
                                               n_up == 0 ~ paste0("REF lower in ", L, " (all individuals)"), TRUE ~ "mixed (phase differs)"))
  gt <- st %>% dplyr::inner_join(snp_to_gene, by = "SNP") %>% dplyr::filter(!is.na(p)) %>% dplyr::group_by(gene_id) %>%
    dplyr::summarise(n_snps = dplyr::n(), acat_p = acat(p), max_abs_delta = max(max_abs_delta), .groups = "drop") %>%
    dplyr::mutate(contrast = ct, padj = p.adjust(acat_p, "BH"), sig = !is.na(padj) & padj < FDR_SIG & max_abs_delta >= ABS_DEV_SIG,
                  note = "unphased, no direction")
  snp_l[[L]] <- st; gene_l[[L]] <- gt; phi_l[[L]] <- dplyr::mutate(r$phi, contrast = ct)
  bad <- st$status != "ok"; nt_l[[L]] <- data.frame(contrast = rep(ct, sum(bad)), level = rep("SNP", sum(bad)), id = st$SNP[bad], reason = st$status[bad])
}
snp <- dplyr::bind_rows(snp_l); gene <- dplyr::bind_rows(gene_l); not_tested <- dplyr::bind_rows(nt_l)
dispersion <- dplyr::bind_rows(phi_l)
knitr::kable(dispersion, digits = 4, caption = "Pair dispersion per individual (floored at RHO_MIN)")
knitr::kable(head(dplyr::filter(snp, sig), 40), digits = 4, caption = "Significant SNPs (first 40); unphased: direction only when all individuals agree")
```

## Summary

```{r summary}
cnt <- function(tab, lev) {
  if (nrow(tab) == 0) return(NULL)
  dplyr::bind_rows(lapply(split(tab, tab$contrast), function(t) {
    pcol <- if ("p_condition" %in% names(t)) t$p_condition else if ("acat_p" %in% names(t)) t$acat_p else t$p
    up <- if (MODE == "f1") sum(t$sig & t$delta_frac > 0) else if (lev == "SNP") sum(t$sig & t$n_down == 0) else NA_integer_
    dn <- if (MODE == "f1") sum(t$sig & t$delta_frac < 0) else if (lev == "SNP") sum(t$sig & t$n_up == 0) else NA_integer_
    disp <- if (MODE == "f1") sprintf("phi_common %.4f", t$phi_common[1]) else {
      ph <- dispersion$phi_pair[dispersion$contrast == t$contrast[1]]; sprintf("phi_pair %.4f-%.4f", min(ph), max(ph)) }
    data.frame(contrast = t$contrast[1], mode = MODE, level = lev, tested = sum(!is.na(pcol)),
               not_tested = sum(not_tested$contrast == t$contrast[1] & not_tested$level == lev) + sum(is.na(pcol)),
               sig = sum(t$sig), n_up = up, n_down = dn, dispersion = disp, stringsAsFactors = FALSE)
  }))
}
summary_diff <- dplyr::bind_rows(cnt(gene, "gene"), cnt(snp, "SNP"))
knitr::kable(summary_diff, caption = paste0("Summary (n_up / n_down: F1 = ", STRAIN_A, " fraction higher / lower in the tested condition; ",
                                            "outbred SNPs = all individuals agree; outbred genes: no direction)"))
```

## Figure

```{r figure}
fig <- if (MODE == "f1") dplyr::filter(gene, !is.na(p_condition)) %>% dplyr::transmute(contrast, x = delta_frac, p = p_condition, sig) else
  dplyr::filter(snp, !is.na(p)) %>% dplyr::transmute(contrast, x = mean_delta, p = p, sig)
lev <- if (MODE == "f1") "gene" else "SNP"
fig_n <- tapply(fig$sig, fig$contrast, sum); sum_n <- setNames(summary_diff$sig[summary_diff$level == lev], summary_diff$contrast[summary_diff$level == lev])
stopifnot(identical(as.integer(fig_n[names(sum_n)]), as.integer(sum_n)))
cat("Figure sig counts equal Summary counts: TRUE\n")
p1 <- ggplot(fig, aes(x, -log10(p), colour = sig)) + geom_vline(xintercept = 0, linetype = 2) + geom_point(alpha = 0.7) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick")) + facet_wrap(~contrast) +
  labs(x = if (MODE == "f1") paste0("change of the ", STRAIN_A, " fraction (tested - reference)") else "mean change of the REF fraction (unphased)",
       y = "-log10 p", caption = paste0("sig: BH < ", FDR_SIG, " and |change| >= ", ABS_DEV_SIG)) + theme_bw()
print(p1)
ggsave(file.path(RESULTS_DIR, paste0(DATE_TAG, "_ASE_differential.pdf")), p1, width = 10, height = 7)
```

## Export

```{r export}
xlsx_file <- file.path(RESULTS_DIR, paste0(DATE_TAG, "_ASE_differential.xlsx"))
wb <- openxlsx::createWorkbook()
for (nm in c("SNP", "Gene", "Design", "Excluded", "Summary")) openxlsx::addWorksheet(wb, nm)
openxlsx::writeData(wb, "SNP", snp); openxlsx::writeData(wb, "Gene", gene); openxlsx::writeData(wb, "Design", design_tbl)
openxlsx::writeData(wb, "Excluded", dplyr::bind_rows(dplyr::transmute(excluded_chrom, what = sprintf("chromosome %s", chrom), n = sample_site_rows),
                                                      if (nrow(not_tested) > 0) dplyr::count(not_tested, contrast, level, reason, name = "n") %>%
                                                        dplyr::transmute(what = sprintf("%s %s %s", contrast, level, reason), n) else NULL))
openxlsx::writeData(wb, "Summary", summary_diff)
openxlsx::saveWorkbook(wb, xlsx_file, overwrite = TRUE)
saveRDS(list(snp = snp, gene = gene, not_tested = not_tested, design = design_tbl, excluded_chrom = excluded_chrom, snp_gene = snp_to_gene,
             dispersion = dispersion, summary = summary_diff, constants = c(ck$constants, list(RHO_MIN = RHO_MIN, THIN_BP = THIN_BP, REF_CONDITION = REF_CONDITION))),
        file.path(RESULTS_DIR, "ase_differential_checkpoint.rds"))
write.table(summary_diff, file.path(RESULTS_DIR, "summary_numbers_differential.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote", xlsx_file, "and summary_numbers_differential.tsv\n")
sessionInfo()
```
````

Render it with `{RESULTS_DIR}/scripts/run_04_differential.sh`. This is the same script as `run_01_import_qc.sh`, with job name `ase_04_differential`, log `run_04_differential_%j.out`, the Rmd 04 file name and `-n 1 --mem=16G -t 4:00:00` (single-threaded; F1 tests every gene and every SNP). If the job stops with "is confounded with the cross direction", show the message to the user. The design cannot answer the question; never edit the Rmd. Submission order: Step 15.
`````

- [ ] **Step 3: Extend the evaluator**

In `evaluate_stage2.R`, insert before `quit(status = fail)`:

```r
if (what == "differential_f1") {
  ck <- readRDS(file.path(RES, "ase_differential_checkpoint.rds")); g <- ck$gene[ck$gene$contrast == "treat vs ctrl", ]
  m <- merge(tg, g, by = "gene_id", all.x = TRUE); m$sig[is.na(m$sig)] <- FALSE
  print(m[m$class != "null", c("gene_id", "class", "expect_diff", "frac_A_ref", "frac_A_test", "delta_frac", "padj", "direction", "status")], row.names = FALSE)
  ex <- m$expect_diff %in% c("up", "down")
  okD <- m$sig & ((m$expect_diff == "up" & m$delta_frac > 0) | (m$expect_diff == "down" & m$delta_frac < 0))
  crit(all(okD[ex]), sprintf("F1 differential genes recovered with the correct sign: %d of %d", sum(okD[ex]), sum(ex)))
  sp <- m$class %in% c("strain_A_high", "strain_B_high", "maternal", "paternal", "strain_and_maternal")
  crit(sum(m$sig[sp]) == 0, sprintf("no differential call among strain / parent-of-origin genes (no condition effect): %d of %d", sum(m$sig[sp]), sum(sp)))
  nul <- m$class %in% c("null", "dense_null", "bias")
  crit(sum(m$sig[nul]) <= 2, sprintf("null genes: %d false positives of %d (limit 2)", sum(m$sig[nul]), sum(nul)))
  crit(isTRUE(ck$design$cross_direction_term[1]), "the cross-direction term was used (both directions present)")
}
if (what == "differential_outbred") {
  ck <- readRDS(file.path(RES, "ase_differential_checkpoint.rds")); g <- ck$gene[ck$gene$contrast == "treat vs ctrl", ]
  m <- merge(tg, g, by = "gene_id", all.x = TRUE); m$sig[is.na(m$sig)] <- FALSE
  print(m[m$class != "null", c("gene_id", "class", "n_snps", "acat_p", "padj", "max_abs_delta", "sig")], row.names = FALSE)
  ex <- m$class %in% c("diff_phase", "diff_consistent")
  crit(sum(m$sig[ex]) >= 4, sprintf("outbred differential genes detected (gene level, unphased): %d of %d (need >= 4)", sum(m$sig[ex]), sum(ex)))
  crit(sum(m$sig[m$class == "base_imbalanced"]) == 0, "no call among imbalanced but unchanged genes")
  nul <- m$class %in% c("null", "bias")
  crit(sum(m$sig[nul]) <= 2, sprintf("null genes: %d false positives of %d (limit 2)", sum(m$sig[nul]), sum(nul)))
  s <- merge(ck$snp[ck$snp$contrast == "treat vs ctrl" & ck$snp$sig, ], ck$snp_gene, by = "SNP")
  s <- merge(s, tg[, c("gene_id", "class")], by = "gene_id")
  crit(all(s$direction[s$class == "diff_consistent"] == "REF lower in treat (all individuals)") && any(s$class == "diff_consistent"),
       "consistent genes: significant SNPs say 'REF lower in treat (all individuals)'")
  crit(any(s$direction[s$class == "diff_phase"] == "mixed (phase differs)"), "phase-heterogeneous genes: at least one SNP labelled 'mixed (phase differs)'")
}
```

- [ ] **Step 4: Run the checker. Expected: PASS**

- [ ] **Step 5: Render on the direct counts (F1 and outbred) and evaluate**

- **F1:** reuse `$SCR/t3_f1` (Task 3 Step 5). Add `render_from_skill.sh ... 17 ... $P/2026-09-30_t3_f1_04_differential.Rmd` and a render job that runs only Rmd 04 (the Rmd 01 checkpoint exists) and then `evaluate_stage2.R differential_f1 $RD .../f1_recip`.
- **Outbred:** use a new `$SCR/t4_ob` built as in Task 3 Step 5 from `outbred_diff`. Copy `direct_counts/*` (tables, `.unfiltered.table`, `.wasp_stats.tsv`) into `results/2026-09-30_t4_ob/ase_counts/`, with `MODE outbred`, `STRAIN_A REF`, `STRAIN_B ALT`, `WD_NAME t4_ob`. Render Rmd 01 and Rmd 04, then run `evaluate_stage2.R differential_outbred $RD .../outbred_diff`.

Expected: all renders exit 0, `summary_numbers_differential.tsv` has 2 data rows per project (gene and SNP for `treat vs ctrl`), and every evaluator line is `PASS`.

- [ ] **Step 6: Negative tests (Review Focus 2 and 4)**

1. **Confounded F1 design:** copy `$SCR/t3_f1` to `$SCR/t4_conf`. Remove the `AxB`/`treat` and `BxA`/`ctrl` rows from the samples CSV and their tables from `ase_counts`. Re-render Rmd 01, then Rmd 04. Expected: Rmd 04 stops with the message containing `is confounded with the cross direction`.
2. **Unpaired outbred individual:** copy `$SCR/t4_ob` to `$SCR/t4_unp`. Add a ninth row `ob_s9` for a new individual `ind9` in `ctrl` only, reusing a copy of `ob_s1`'s three files under the name `ob_s9`. Re-render Rmd 01 and Rmd 04. Expected: the `Design` table lists `ind9` under `unpaired_individuals`, the SNP results for the other individuals are identical to `$SCR/t4_ob` (compare the `snp` tables with `all.equal` in a small R job), and the evaluator still PASSes.

- [ ] **Step 7: Commit**

```bash
git status && git diff --stat
git add ase-pipeline/ase-pipeline.md ase-pipeline/tests/check_skill.sh ase-pipeline/tests/synthetic/evaluate_stage2.R
git commit -m "ase-pipeline: Rmd 04 differential ASE (Step 17): F1 GLM with direction term, paired direction-free outbred test

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Rmd 02 F1 gene test on thinned SNPs (conditional)

**Gate (controller checks before dispatch):** take the `TASK5_INPUT` line of `$SCR/logs/t1_final.log`. Run this task only if `rmd02_unthinned > 0.065` (or `rmd02_dense_unthinned > 0.08`) and `rmd02_thinned <= 0.065`. Otherwise skip it. In that case the Stage 1 limitation text stays, and Task 7 adds the measured sizes to the README limitation. Record the ruling in the ledger. This task also depends on planner decision 13.

**Files:**
- Modify: `ase-pipeline/ase-pipeline.md` Step 13:
  - the constants chunk: add `THIN_BP     <- {THIN_BP}` after the `RHO_MIN` line;
  - the `genemap` chunk: add `exon_pos` exactly as in Step 16;
  - the `snp` chunk F1 branch: pass the thinned gene labels to `bb_estimate_rho_gene`;
  - the `gene` chunk F1 branch: thinned SNPs only, plus `n_snps_used`;
  - the note texts.
- Modify: Step 14 bullet "**Independence assumption (F1 gene LRT and `bb_estimate_rho_gene`).**" (replaced), and `check_skill.sh` (the `need` lines of the final-review fix wave that pin the old note text are replaced by the lines below).

**Interfaces:**
- Consumes: `thin_snps` and `THIN_BP` (Task 1).
- Produces: the Rmd 02 F1 `Gene` sheet gains `n_snps_used`, and its `note` column becomes "thinned: at most one SNP per THIN_BP window (exon coordinates), so no read pair is counted twice; approximate (ignores alternative splicing)". No other output changes.

- [ ] **Step 1: Failing checker lines (and retire the old ones)**

In `check_skill.sh`, delete these four lines of the Stage 1 final-review block (they pin the text this task replaces):
- `need "**Independence assumption (F1 gene LRT and \`bb_estimate_rho_gene\`).**"`
- `need "come from simulations with **independent SNPs**"`
- `need 'note = "SNP counts treated as independent: may be anti-conservative in SNP-dense genes at moderate depth"'`
- `need "**Limitation (F1):** SNP counts that share read pairs are treated as independent"`

The README checks for the Stage 1 limitation are handled in Task 7. Append:

```bash
# --- Stage 2 Task 5 (Rmd 02 F1 gene test on thinned SNPs)
need 'note = "thinned: at most one SNP per THIN_BP window (exon coordinates), so no read pair is counted twice; approximate (ignores alternative splicing)"'
need 'snp_gene$thin_keep'
need "**Thinned SNPs (F1 gene LRT and \`bb_estimate_rho_gene\`).**"
forbid "Planned for Stage 2 (not implemented): thin each gene's SNPs"
# --- end Stage 2 Task 5
```
Run the checker. Expected: FAIL (3 needs, 1 forbid).

- [ ] **Step 2: Edit Rmd 02 (Step 13)**

- In the `genemap` chunk, compute `exon_pos` with the Step 16 `exd`/`offset` code (replace the `hits` construction on `ex` with the one on `exr`, keeping `snp_gene` columns `contig, position, gene_id` plus `exon_pos`).
- Then add:

  ```r
  thin_depth <- dplyr::inner_join(sites, snp_gene, by = c("contig", "position")) %>%
    dplyr::group_by(gene_id, contig, position, exon_pos) %>% dplyr::summarise(depth = sum(total), .groups = "drop") %>%
    dplyr::group_by(gene_id) %>% dplyr::mutate(thin_keep = thin_snps(exon_pos, depth, THIN_BP)) %>% dplyr::ungroup()
  snp_gene <- dplyr::left_join(snp_gene, dplyr::select(thin_depth, gene_id, contig, position, thin_keep), by = c("gene_id", "contig", "position"))
  snp_gene$thin_keep[is.na(snp_gene$thin_keep)] <- FALSE
  ```
- In the `snp` chunk F1 branch, change the gene labels passed to `bb_estimate_rho_gene` to the thinned SNPs only:

  ```r
  sg_keep <- snp_gene[snp_gene$thin_keep, ]
  g <- sg_keep$gene_id[match(paste(d$contig, d$position), paste(sg_keep$contig, sg_keep$position))]
  ```
- In the `gene` chunk F1 branch, build `snp_by_gene` from the thinned SNPs only:

  ```r
  snp_by_gene_f1 <- dplyr::inner_join(snp, dplyr::filter(snp_gene, thin_keep), by = c("contig", "position"))
  ```

  Use it in place of `snp_by_gene` for F1. Add `n_snps_used = nrow(g)` and the new `note` text. Replace the NOTE `cat(...)` and the table caption with: "F1 gene test: at most one SNP per THIN_BP window (exon coordinates) is used, so no read pair is counted twice (approximate: alternative splicing can bring distant exons into one fragment)."
- In the Step 13 bullets and in the "## Gene level" prose of Rmd 02, replace every sentence about SNP counts being treated as independent (the Step 13 gene-level bullet's `note` clause and the "**Limitation (F1):**" sentence) with the thinned wording above.
- Replace the Step 14 bullet "**Independence assumption ...**" with "**Thinned SNPs (F1 gene LRT and `bb_estimate_rho_gene`).**", followed by the Task 1 measured sizes (`rmd02_unthinned` → `rmd02_thinned`, dense genes likewise) and the statement that the SNP-level tests are unchanged.

- [ ] **Step 3: Re-render Rmd 02 on the Stage 1 acceptance tables and the Stage 2 direct counts**

1. **Stage 1 tables:** copy `/net/bmc-lab3/data/bcc/yannvrb/ase_accept2/f1/` to `$SCR/t5_s1`. Change `RESULTS_DIR` and `CWD` in a `values.tsv` to the copy (`TODAY 2026-09-29`, `WD_NAME f1`, `GTF_PATH` = the copied genome GTF, `THIN_BP 500`). Regenerate Rmd 02 with `render_from_skill.sh ... 13 ...`, and render it in an `sbatch` job (the Rmd 01 checkpoint exists).
   - Evaluate with a small R job: the four planted genes are `sig` with the correct direction in every sample, and null genes are not.
   - Expected: 24 of 24 planted and 0 of 42 null, as in Stage 1. Report any change with numbers; a drop below 20 of 24 is BLOCKED.
2. **Stage 2 direct counts:** render Rmd 02 on `$SCR/t3_f1`. Expected: exit 0, and `n_snps_used` ≤ 2 for every `dense_null` gene.

- [ ] **Step 4: Run the checker (expected PASS) and the unit tests (Stage 1 file only; nothing in the stats block changed); commit**

```bash
bash ase-pipeline/tests/check_skill.sh ase-pipeline/ase-pipeline.md
git status && git diff --stat
git add ase-pipeline/ase-pipeline.md ase-pipeline/tests/check_skill.sh
git commit -m "ase-pipeline: Rmd 02 F1 gene test on thinned SNPs (no read pair counted twice)

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Wizard, submission chain and summary page integration

**Files:**
- Modify: `ase-pipeline/ase-pipeline.md`: Step 5 (F1 `cross_direction` rule and validation), Step 7 (menu text, reference condition, design checks), Step 15 (submission block, run scripts, summary page), Notes (Stage bullet, statistics-block bullet)
- Modify: `ase-pipeline/tests/check_skill.sh`

**Interfaces:**
- Consumes:
  - the file names of Tasks 3 and 4: `run_03_reciprocal.sh`, `run_04_differential.sh`, `summary_numbers_reciprocal.tsv`, `summary_numbers_differential.tsv`, `{TODAY}_{WD_NAME}_03_reciprocal.html`, `_ASE_reciprocal.xlsx`, `_ASE_reciprocal.pdf`, `_04_differential.html`, `_ASE_differential.xlsx`, `_ASE_differential.pdf`;
  - the `{REF_CONDITION}` placeholder.
- Produces: `{REF_CONDITION}` (set in Step 7), and the `{ANALYSES}` values `per_sample`, `reciprocal`, `differential` (and `phaser`, still "available in a later stage").

- [ ] **Step 1: Failing checker lines (and update the replaced Stage 1 lines)**

In `check_skill.sh`, replace `need "offered only when \`{MODE}\` = \`f1\` and both \`cross_direction\` values"` with the first line below, and replace `need "offered only when at least two conditions each have replicates"` with the second line. Then append the rest:

```bash
need "offered only when \`{MODE}\` = \`f1\` and at least 2 samples have \`cross_direction\` \`AxB\` and at least 2 have \`BxA\`"
need "offered only when at least two conditions each have replicates and"
# --- Stage 2 Task 6 (wizard, chain, summary)
need "\`cross_direction\` is exactly \`AxB\` or \`BxA\`"
need "Which condition is the reference (baseline)?"
need "Store it as \`{REF_CONDITION}\`"
need "condition and cross direction are confounded"
need "at least 2 individuals sampled in both"
need 'R3=$(sbatch -p bcc --parsable --dependency=afterok:$R2 $S/run_03_reciprocal.sh)'
need 'R4=$(sbatch -p bcc --parsable --dependency=afterok:$R2 $S/run_04_differential.sh)'
need "### Reciprocal F1 section"
need "### Differential ASE section"
need "offered only when \`{MODE}\` = \`outbred\`; available in a later stage"
need "pasted verbatim into Rmd 02, 03 and 04"
forbid "Stage 1 implements only the always-on analysis"
forbid "**Stage 1 only.**"
# --- end Stage 2 Task 6
```
Run the checker. Expected: FAIL (the new lines and the two replaced lines).

- [ ] **Step 2: Step 5 — cross direction rule**

In Step 5's F1 bullet, replace "and, per sample, `cross_direction` (for example `AxB` or `BxA`, maternal strain first)" with:

"and, per sample, `cross_direction`: `cross_direction` is exactly `AxB` (mother `{STRAIN_A}`, father `{STRAIN_B}`) or `BxA` (mother `{STRAIN_B}`). The letters are literally A and B, not strain names, because the parent-of-origin sign of Rmd 03 depends on it. Ask per group of samples (numbered: 1. `AxB` · 2. `BxA`), never as free text."

Add to the Step 5 validation rules: "in F1 mode every `cross_direction` is `AxB` or `BxA` (stop and ask again otherwise)".

- [ ] **Step 3: Step 7 — menu and design checks**

Replace the Step 7 list and its closing paragraph with:

```markdown
1. **Per-sample allelic imbalance** (always on, cannot be deselected): reference-bias diagnostic plus per-gene and per-site allelic ratios for every sample (Rmd 01 and Rmd 02).
2. **Reciprocal F1 analysis** (strain effect versus parent-of-origin effect, Rmd 03, Step 16): offered only when `{MODE}` = `f1` and at least 2 samples have `cross_direction` `AxB` and at least 2 have `BxA`.
3. **Differential ASE between conditions** (Rmd 04, Step 17): offered only when at least two conditions each have replicates and the design allows a paired or adjusted comparison. In F1 mode, condition and cross direction must not be confounded: when both directions are present, it must not be true that every sample of one condition has one direction and every sample of the other condition has the other. In outbred mode, at least 2 individuals sampled in both the reference condition and each other condition are needed (the test pairs each individual's samples).
4. **phASER haplotype phasing**: offered only when `{MODE}` = `outbred`; available in a later stage.

Show only the options whose preconditions hold, and say why any other is hidden, for example "Differential ASE is hidden: condition and cross direction are confounded (all ctrl samples are AxB and all treat samples BxA), so a condition effect cannot be told apart from a parent-of-origin effect." Store the selection as `{ANALYSES}`. If the user chooses phASER, say "available in a later stage" and continue.

If differential ASE is selected, ask: "Which condition is the reference (baseline)?" (numbered list of the conditions in `{SAMPLES_CSV}`, in the order they first appear). Store it as `{REF_CONDITION}`. Every other condition is compared with it, one contrast at a time. For each contrast in outbred mode, list the individuals that are sampled in only one of the two conditions: they are left out of that contrast.
```

- [ ] **Step 4: Step 15 — submission chain and summary page**

In the Step 15 submission block, insert after the `R2=` line:

```bash
# 5. ONLY IF reciprocal F1 was selected (Step 7): Rmd 03 after Rmd 02
R3=$(sbatch -p bcc --parsable --dependency=afterok:$R2 $S/run_03_reciprocal.sh)
# 6. ONLY IF differential ASE was selected (Step 7): Rmd 04 after Rmd 02 (runs in parallel with Rmd 03)
R4=$(sbatch -p bcc --parsable --dependency=afterok:$R2 $S/run_04_differential.sh)
```

Change the `echo` line to `echo "prep $P, array $A, Rmd01 $R1, Rmd02 $R2, Rmd03 ${R3:-none}, Rmd04 ${R4:-none}"`. Add this sentence after the block: "Delete the lines of an analysis that was not selected. Wait for the last job submitted (`R4`, else `R3`, else `R2`) before writing the summary page; if Rmd 03 or Rmd 04 fails, the per-sample results stay valid and the page says which analysis failed and why."

In "Dropping a failed sample", append: "Re-run Rmd 03 and Rmd 04 after Rmd 02 in the same way (a dropped sample can change whether a design is still valid, so re-check Step 7's preconditions first)."

Add two subsections at the end of the summary-report part, before "**Verify before finishing.**":

```markdown
### Reciprocal F1 section (only when Rmd 03 ran)

Read `{RESULTS_DIR}/summary_numbers_reciprocal.tsv` with `awk -F'\t'` (columns `test, genes_tested, genes_not_tested, sig_genes, direction_pos, n_pos, direction_neg, n_neg, phi_common, samples_AxB, samples_BxA`) and show one table row per test ("strain effect" and "parent-of-origin effect") with the tested and significant genes and the counts per direction. Add the sentences: "p is the `{STRAIN_A}` fraction; d = +1 when `{STRAIN_A}` is the mother (AxB). The strain effect is sensitive to reference mapping bias (see the flags above); the parent-of-origin effect is not. X, Y and MT are not tested." Add link cards to `{TODAY}_{WD_NAME}_03_reciprocal.html`, `{TODAY}_{WD_NAME}_ASE_reciprocal.xlsx` and `{TODAY}_{WD_NAME}_ASE_reciprocal.pdf`.

### Differential ASE section (only when Rmd 04 ran)

Read `{RESULTS_DIR}/summary_numbers_differential.tsv` with `awk -F'\t'` (columns `contrast, mode, level, tested, not_tested, sig, n_up, n_down, dispersion`) and show one row per contrast and level. Add, in F1 mode: "n_up / n_down: `{STRAIN_A}` fraction higher / lower in the tested condition than in `{REF_CONDITION}`". In outbred mode: "unphased: SNP directions are counted only when all individuals agree; gene calls have no direction". Add link cards to `{TODAY}_{WD_NAME}_04_differential.html`, `{TODAY}_{WD_NAME}_ASE_differential.xlsx` and `{TODAY}_{WD_NAME}_ASE_differential.pdf`.
```

In the "Verify before finishing" paragraph, extend the list of reported paths with "and, when they ran, the Rmd 03/04 HTML files, xlsx files and summary TSVs".

- [ ] **Step 5: Notes for the assistant**

Replace the bullet `- **Stage 1 only.** ...` with:

"- **Analyses available.** Per-sample imbalance (Rmd 01/02), reciprocal F1 (Rmd 03) and differential ASE (Rmd 04) are implemented. phASER is a later stage: say "available in a later stage" and continue. Never add a model or package that the Step 14 block does not contain. `lme4` and `aod` are not used by the Rmds (Step 14 explains why)."

Replace "is pasted verbatim into Rmd 02; never edit it inside a generated Rmd" with "is pasted verbatim into Rmd 02, 03 and 04; never edit it inside a generated Rmd".

- [ ] **Step 6: Test the summary-page instructions on the Task 3/4 outputs**

In `$SCR/t3_f1/results/2026-09-30_t3_f1/`, run Rmd 02 via `render_from_skill.sh ... 13` (needed for `summary_numbers.tsv`). Then write `t3_f1_summary_report.html` exactly as Step 15 instructs: bash plus awk, reading the three TSVs, no R. Run the Step 15 href loop and the no-`http`/no-`/` check. Expected: both print nothing. Do the same for `$SCR/t4_ob` (outbred; differential section only). Keep both pages for the Task 8 comparison.

- [ ] **Step 7: Run the checker (expected PASS); commit**

```bash
bash ase-pipeline/tests/check_skill.sh ase-pipeline/ase-pipeline.md
git status && git diff --stat
git add ase-pipeline/ase-pipeline.md ase-pipeline/tests/check_skill.sh
git commit -m "ase-pipeline: wizard menu, reference condition, Rmd 03/04 run chain and summary sections

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 7: README, root README and final checker consolidation

**Files:**
- Modify: `ase-pipeline/README.md`, root `README.md` (the `/ase-pipeline` row), `ase-pipeline/tests/check_skill.sh`

**Interfaces:**
- Consumes: the measured numbers of `$SCR/logs/t1_final.log` and the Task 3/4 evaluator logs, and the Task 5 ruling (applied or skipped).
- Produces: README text that Task 8 fills with the acceptance result, using the marker line `Stage 2 synthetic acceptance run: PENDING`.

- [ ] **Step 1: Failing checker lines**

```bash
# --- Stage 2 Task 7 (README)
RD="$HERE/../README.md"
grep -qF "(Stages 1 and 2)" "$RD" || { echo "FAIL: README title must say Stages 1 and 2"; fail=1; }
grep -qF "Reciprocal F1 (Rmd 03)" "$RD" && grep -qF "Differential ASE (Rmd 04)" "$RD" || { echo "FAIL: README must describe Rmd 03 and Rmd 04"; fail=1; }
grep -qF "summary_numbers_reciprocal.tsv" "$RD" && grep -qF "summary_numbers_differential.tsv" "$RD" || { echo "FAIL: README outputs must list the Stage 2 TSVs"; fail=1; }
grep -qF "direction-free" "$RD" || { echo "FAIL: README must explain the direction-free outbred test"; fail=1; }
grep -qF "X, Y and MT" "$RD" || { echo "FAIL: README must state the X/Y/MT exclusion"; fail=1; }
grep -qE "Stage 2 synthetic acceptance run: (PENDING|DONE)" "$RD" || { echo "FAIL: README must state the Stage 2 acceptance status"; fail=1; }
! grep -qF "Stage 1 covers **per-sample allelic imbalance only**" "$RD" || { echo "FAIL: README still says Stage 1 only"; fail=1; }
grep -q "/ase-pipeline.*reciprocal F1" "$HERE/../../README.md" || { echo "FAIL: root README row must mention reciprocal F1"; fail=1; }
# --- end Stage 2 Task 7
```
If Task 5 was applied, also delete the Stage 1 README checks `grep -qF "SNPs that share read pairs are treated as independent"` and `grep -qF "Future work (Stage 2): thin each gene's SNPs"`, and add `grep -qF "thinned to one SNP per" "$RD" || { echo "FAIL: README must describe the thinned F1 gene test"; fail=1; }`. If Task 5 was skipped, keep them. Run the checker. Expected: FAIL.

- [ ] **Step 2: Edit `ase-pipeline/README.md`**

1. **Title:** "# `/ase-pipeline` — Allele-Specific Expression Skill (Stages 1 and 2)". Replace the "Stage 1 covers ..." paragraph with: "Stages 1 and 2 cover per-sample allelic imbalance, the reciprocal F1 analysis (strain versus parent-of-origin effects) and differential ASE between conditions. phASER phasing (Stage 3) is not included."
2. **Wizard steps table:** change row 7 to "Analysis menu: per-sample (always), reciprocal F1, differential ASE (with reference condition); phASER later". Add rows `16 | Rmd 03: reciprocal F1` and `17 | Rmd 04: differential ASE`.
3. **Outputs table:** add rows for `run_03_reciprocal.sh`/`run_04_differential.sh`, the two HTML files, the two xlsx files (with their sheets), the two PDFs, the two checkpoints and the two `summary_numbers_*.tsv` files. Update the chain sentence to "prep -> array -> Rmd 01 -> Rmd 02 -> Rmd 03 and Rmd 04 (each `afterok` on Rmd 02, only when selected)".
4. **Key design points:** add three bullets.
   - "**Reciprocal F1 (Rmd 03)**: ..." Summarise Step 16: p = strain-A fraction, d = +1 for AxB, b0 strain, b1 parent of origin, rows = samples with thinned SNP sums, dispersion between replicate animals with the maximum rule, LRTs, X/Y/MT excluded, b1 robust to constant mapping bias.
   - "**Differential ASE (Rmd 04)**: ..." Summarise Step 17: F1 GLM with the direction term and the confounding stop; outbred paired, direction-free per-SNP test summed over individuals, ACAT genes, no lme4 and why.
   - "**Base-R models**: ..." Cover planner decision 3.
5. **Validation status:** add a "Stage 2" subsection.
   - Unit tests: copy the measured sizes, power and sign rates from `t1_final.log`, citing `tests/r/test_ase_stats_stage2.R` and the job needing at least 8 CPUs.
   - The render and evaluator results on direct counts (Tasks 3 and 4), the negative tests, and (if run) the Task 5 re-render numbers.
   - The line `Stage 2 synthetic acceptance run: PENDING`.
   - "What was not exercised": real data, more than two conditions, unbalanced real designs, multi-transcript genes (thinning uses the exon union), and genome-scale runtime (projected only, from the benchmark).
6. **Known limitations:** add these items.
   - Thinning uses exon-union coordinates, so alternative splicing can bring two distant exons into one fragment. THIN_BP should match the longest fragment.
   - One dispersion per analysis plus each gene's own estimate when larger; with few replicates the per-gene estimate is noisy (conservative).
   - Outbred differential: at least 2 paired individuals; the pair dispersion includes truly changed SNPs (conservative); gene calls have no direction.
   - X, Y and MT are not analysed in Rmd 03/04.
   - A reciprocal design needs both directions with replicates; condition-by-direction interactions are not modelled.

   Replace "Stage 1 covers only per-sample imbalance; ..." with "phASER (Stage 3) is not implemented." Apply the Task 5 ruling: if Task 5 was applied, replace the SNP-dependence limitation with the thinned-test description and its measured sizes; if it was skipped, append the measured Stage 1 gene-test sizes from `TASK5_INPUT` to that limitation.
7. **Root `README.md`:** in the `/ase-pipeline` row, replace "Stage 1: per-sample imbalance only" with "per-sample imbalance, reciprocal F1 (strain and parent-of-origin effects) and differential ASE between conditions; phASER not included".

- [ ] **Step 3: Run the checker (expected PASS); commit**

```bash
bash ase-pipeline/tests/check_skill.sh ase-pipeline/ase-pipeline.md
git status && git diff --stat
git add ase-pipeline/README.md README.md ase-pipeline/tests/check_skill.sh
git commit -m "ase-pipeline: README and registration for Stage 2

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
```

---

### Task 8: Synthetic acceptance on the cluster (controller-run, not a subagent)

Run by the controller, because it needs cluster jobs, judgment and the wizard followed by script. It is the same method as Stage 1 Task 8b: scripts and Rmds are generated from the branch skill's code blocks by extraction and placeholder substitution (`render_from_skill.sh` for the Rmds, and a scratch `gen.sh` that cuts each script block for the bash scripts). The installed copy `~/.claude/commands/ase-pipeline.md` is not overwritten without the user's approval: the acceptance reads the branch file.

**Files:** Modify `ase-pipeline/README.md` (the acceptance line and results). The evaluation uses the committed `evaluate_stage2.R`.

- [ ] **Step 1: Data and projects**

The Stage 2 synthetic data is `/net/bmc-lab3/data/bcc/yannvrb/ase_synthetic_stage2` (Task 2; re-check with `check_simulation_stage2.sh`: PASS). Create `/net/bmc-lab3/data/bcc/yannvrb/ase_accept_stage2/{f1r,obd}`, each with a writable copy of `genome/`.

The wizard answers, following Stage 1 8b:
- email `yannvrb@mit.edu`;
- F1: `STRAIN_A` C57BL_6NJ (REF), `STRAIN_B` A_J, parental VCF `f1_recip/parental_snps.vcf` (sites-only, `MIN_PARENTAL_SITES=20`), sample sheet = `f1_recip/samples.csv` (already in the Step 5 format), analyses `1,2,3`, reference condition `ctrl`;
- outbred: VCFs `ind*.all.vcf`, sheet = `outbred_diff/samples.csv`, analyses `1,3`, reference `ctrl`;
- constants at their defaults (`THIN_BP` 500), and the small-genome resource tier.

Run Step 2's `check_r_packages.sh` job in each directory.

- [ ] **Step 2: Submit the Step 15 chain in each project**

Submit prep → array → Rmd 01 → Rmd 02 → (F1: Rmd 03 and Rmd 04; outbred: Rmd 04), all `sbatch -p bcc --parsable --dependency=afterok`. Record the job ids in `results/*/logs/submission_ids.txt`. Check the dependencies and the commands with `scontrol show job` (`sacct` is unreachable).

- [ ] **Step 3: Evaluate against truth**

Run an `sbatch` job:
- `evaluate_stage2.R reciprocal <f1r RESULTS_DIR> .../f1_recip`
- `evaluate_stage2.R differential_f1 <f1r RESULTS_DIR> .../f1_recip`
- `evaluate_stage2.R differential_outbred <obd RESULTS_DIR> .../outbred_diff`

Acceptance criteria, all stated with numbers in the report:
- **F1 reciprocal:** the 5 strain genes (2 `strain_A_high`, 2 `strain_B_high`, `strain_and_maternal`) are significant with the correct sign. The 5 parent-of-origin genes (2 maternal, including the 0.95 near-monoallelic one; 2 paternal; `strain_and_maternal`) are significant with the correct sign. There are 0 parent-of-origin calls among the 4 strain-only genes, 0 strain calls among the 4 parent-of-origin-only genes, and at most 2 false positives per test among the 47 null, dense and bias genes. The false-positive rate is stated.
- **F1 differential:** the 4 differential genes (`diff_up` ×2, `diff_down`, `diff_on_strain`) are significant with the correct sign. There are 0 calls among the 9 strain and parent-of-origin genes, at most 2 of 47 null genes, and the direction term was used.
- **Outbred differential:** at least 4 of the 5 changed genes are detected (gene level). There are 0 of 3 `base_imbalanced` calls and at most 2 of 52 null calls. The consistent genes' significant SNPs say "REF lower in treat (all individuals)", and at least one phase-heterogeneous gene shows "mixed (phase differs)".
- **Every Rmd** (01-04 in F1; 01, 02 and 04 in outbred) builds without hand edits, the figure/summary assertions hold, and every `ase_counts/*.table` is non-empty.
- **Rmd 02 on the aligned Stage 2 F1 data** (per-sample): the planted strain genes are detected in each sample, direction-consistent with the sample's cross direction. It is reported, not gated.
- **Aligned versus direct counts:** compare the aligned results with the Task 3/4 direct-count results (same genes called), and state any difference and its cause (mapping, masking, sequencing errors).

- [ ] **Step 4: Summary pages**

Write both `{WD_NAME}_summary_report.html` pages by Step 15 (awk from the TSVs, no R). The href loop and the no-`http` check print nothing.

- [ ] **Step 5: Defects, README, ledger**

Defects go through the process: one fix subagent per defect wave, a scoped re-review, and a re-run of the affected steps. When all criteria pass:
- replace `Stage 2 synthetic acceptance run: PENDING` with `Stage 2 synthetic acceptance run: DONE (2026-..)` and the numbers;
- run the checker;
- commit (`git status && git diff --stat` first) with message `ase-pipeline: record Stage 2 synthetic acceptance results` plus the Co-Authored-By trailer.

Then comes the whole-branch final review. Install to `~/.claude/commands/` and merge only with the user's approval. Never push.

---

## Self-review notes

- **Spec coverage:**
  - Rmd 03: the model, the sign conventions and the replicate requirement (Tasks 1 and 3).
  - Rmd 04: the F1 per-gene and per-SNP models, the individual term for outbred and the paired per-SNP tests, re-specified as the direction-free paired test with the reason recorded in the planner decisions and the spec amendment (Tasks 1 and 4).
  - aod/lme4 usage re-specified (Task 1 Step 9).
  - Wizard menu constraints (Task 6).
  - Run scripts and the `afterok` chain (Tasks 3, 4 and 6).
  - Summary page (Task 6).
  - Synthetic reciprocal design and ground truth (Task 2).
  - Unit tests with simulation (Task 1).
  - Checker (every task).
  - README validation status (Tasks 7 and 8).
  - Live acceptance (Task 8).
  - The spec's "Every Rmd builds from the real checkpoint files without hand edits" is checked in Tasks 3, 4 and 8.
  - Stage 3 is excluded.
- **Lessons carried from Stage 1:**
  - Statistics come before any Rmd text, each with simulation gates: null size across 10 seeds with pooled and worst-seed limits, power, sign, dispersion recovery, heterogeneous dispersion, unbalanced designs, single-SNP genes, the `RHO_MIN` floor, fit failures and complete separation.
  - Summing counts is used only after thinning, with a replicate-level dispersion, and a fragment-level simulation of shared read pairs gates it.
  - Neyman-Scott bias is avoided by the df-divided moment estimator (section 4 checks it with 2 coefficients per 6 rows).
  - Imbalance-inflated dispersion is avoided because the full model carries the effects.
  - The acceptance data exercises every new branch, because Stage 1's acceptance found defects that the unit tests missed.
- **Interface consistency:**
  - The function names and output columns defined in Task 1 are the ones used in Tasks 3 and 4 (`beta_(Intercept)`, `beta_d`, `beta_cond`, `p_strain`, `p_parent_of_origin`, `p_condition`, `phi_common`, `status`), and in the evaluator (`b0_strain`, `b1_parent_of_origin`, `delta_frac`, `sig`, `direction`, `snp_gene`).
  - The TSV columns of Tasks 3 and 4 match Task 6's summary sections.
  - `THIN_BP` and `REF_CONDITION` are defined once (Step 8 and Step 7) and substituted in every Rmd.
- **Review Focus:** each of the 5 items has a test in its owning task: sign in Tasks 1 and 3; confounding in Tasks 1 and 4; monoallelic genes and X in Tasks 1 and 3; unpaired individuals in Tasks 1 and 4; SNP-dense genes in Tasks 1 and 3.
- **Known risk left to execution:** the R code in this plan has not been run. The binding parts are the tests and gates, not the exact code lines. Coding defects are fixed in the function. Statistical gate failures go to the controller (BLOCKED), never to a weaker gate.
