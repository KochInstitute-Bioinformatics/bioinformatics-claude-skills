#!/bin/bash
# Static checks for ase-pipeline.md. Usage: check_skill.sh <skill.md>
set -u
SKILL=${1:?usage: check_skill.sh <skill.md>}
FIX=$(dirname "$0")/fixtures
fail=0

if [ ! -s "$SKILL" ]; then echo "FAIL: skill file missing or empty: $SKILL"; exit 1; fi
for f in star_help.txt gatk_ASEReadCounter_help.txt; do
  [ -s "$FIX/$f" ] || { echo "cannot read fixture $FIX/$f"; exit 2; }
done

need()   { grep -qF -- "$1" "$SKILL" || { echo "FAIL: missing required text: $1"; fail=1; }; }
forbid() { ! grep -qF -- "$1" "$SKILL" || { echo "FAIL: forbidden text present: $1"; fail=1; }; }

# (a) every --flag in the skill is a STAR/GATK option from the recorded help fixtures,
#     or an allowlisted flag of another tool (never add a STAR or GATK flag here).
# GATK: options are listed as "--name,-x <Type>" at line start.
# STAR: parameters are "name  default" at line start, followed by an indented "type: ..." line.
tool_flags=$( { grep -oE '^--[A-Za-z][A-Za-z0-9_-]*' "$FIX/gatk_ASEReadCounter_help.txt"
                awk 'prev ~ /^[A-Za-z][A-Za-z0-9_]*( |$)/ && $0 ~ /^ +(string|int|double|uint|bool|-|[A-Za-z0-9]+\(?s?\)?:)/ {split(prev,a," "); print a[1]} {prev=$0}' "$FIX/star_help.txt"
              } | sed 's/^--//' | sort -u )
# Non-STAR/non-GATK flags: sbatch/singularity (bind mem mail-type mail-user array dependency parsable),
# bcftools (regions samples genotype types min-alleles max-alleles rename-chrs).
allow="bind mem mail-type mail-user array dependency parsable regions samples genotype types min-alleles max-alleles rename-chrs"
[ -n "${SHOW_TOOL_FLAGS:-}" ] && echo "$tool_flags"
for flag in $(grep -oE '(^|[ `(=])--[A-Za-z][A-Za-z0-9_-]*' "$SKILL" | sed -E 's/^[^-]*--//' | sort -u); do
  echo "$tool_flags $allow" | tr ' ' '\n' | grep -qx -- "$flag" || { echo "FAIL: flag not in fixtures or allowlist: --$flag"; fail=1; }
done

# (b) content checks, appended per task -----------------------------------
# --- Task 3 (Steps 0-3)
need "# ase-pipeline — allele-specific expression"
need "## Step 0"
need "## Step 3"
forbid "module add htslib"
forbid "module add bcftools"
need "singularity exec --bind"
need "|| { echo"
need "#SBATCH -p bcc"
need "sbatch -p bcc"
forbid "module add samtools"
forbid "module add star"
forbid "module add gatk"
forbid "module add picard"
forbid "module load htslib"
# --- end Task 3
# --- Task 4 (Steps 4-9); strings chosen so only Task 4 text satisfies them
for n in 4 5 6 7 8 9; do need "## Step $n"; done
need "star_ase_sjdb{SJDB_OVERHANG}"
need "star_ase_masked_{STRAIN_A}_{STRAIN_B}_sjdb{SJDB_OVERHANG}"   # final-review I2: F1 index keyed on strain pair + mask
need "genomeSAindexNbases"
need "max(4,"
need "{SAMPLES_CSV}"
need "sample,fastq_1,fastq_2,condition,cross_direction,individual"
need "MIN_DEPTH"
need "FDR_SIG"
need "ABS_DEV_SIG"
need "BIAS_TOL"
need "offered only when \`{MODE}\` = \`f1\` and at least 2 samples have \`cross_direction\` \`AxB\` and at least 2 have \`BxA\`"
need "offered only when at least two conditions each have replicates and"
need "offered only when \`{MODE}\` = \`outbred\`"
need "replicates = at least 2 samples in the condition"
need "#SBATCH --array=1-{ARRAY_N}"
need "-n 8 --mem=48G -t 4:00:00"
need "available in a later stage"
need "f1_het_sites.vcf.gz"
need "RNA-derived genotypes are circular"
need "Step 5, part 2"
need "first \`fastq_1\` of each of up to 5 distinct samples"
forbid "strandedness"
# --- end Task 4
# --- Task 5 (Steps 10-11); strings chosen so only Task 5 text satisfies them
for n in 10 11; do need "## Step $n"; done
need 'singularity exec --bind "$BIND" "$SIF"'
need "bcftools consensus"
need "third allele"
need "--waspOutputMode SAMtag"
need "--varVCFfile"
need "vW"
need "MarkDuplicates"
need "ASEReadCounter"
need "--min-mapping-quality"
need "--min-base-quality"
need "--count-overlap-reads-handling"
need "SLURM_ARRAY_TASK_ID"
need "set -uo pipefail"
need "needs non-empty"
need "ase_counts/\$SAMPLE.table"
need "--outSAMattrRGline ID:"
need "--outSAMattributes NH HI AS nM vA vG vW"
need "wasp_stats.tsv"
need "prep_f1_reference.sh"
need "prep_genotypes.sh"
need "align_count_f1.sh"
need "align_wasp_count.sh"
need "extract_mgp_parental_vcf.sh"
need "bcftools concat"
forbid "--outSAMattributes All"
forbid "module add bcftools"
need "concat_mgp_parental_vcf.sh"
need "rm -f \"\$TABLE\" \"\$STATS\""
need "faidx -n 60"
need "soft-masked"
# --- end Task 5
# --- Task 6 (Steps 12-14); strings chosen so only Task 6 text satisfies them
for n in 12 13 14; do need "## Step $n"; done
need "# --- ase-stats-begin"
need "# --- ase-stats-end"
need "bb_pvalue <- function"
need "acat <- function"
need "ase_checkpoint.rds"
need "_ASE_imbalance.xlsx"
need "<- {MIN_DEPTH}"
need "never by globbing"
need "GenomicRanges::findOverlaps"
need "unphased, no direction"
need "_01_import_qc.Rmd"
need "_02_imbalance.Rmd"
need 'knitr::opts_chunk$set(cache = FALSE, echo = TRUE'
grep -qx "options(scipen = 9)" "$SKILL" || { echo "FAIL: missing Rmd line options(scipen = 9)"; fail=1; }
forbid "source("
[ "$(grep -c '^# --- ase-stats-begin' "$SKILL")" -eq 1 ] || { echo "FAIL: need exactly one ase-stats-begin marker"; fail=1; }
[ "$(grep -c '^# --- ase-stats-end' "$SKILL")" -eq 1 ] || { echo "FAIL: need exactly one ase-stats-end marker"; fail=1; }
# --- end Task 6
# --- Task 6 fix (F1 gene level: free-mean rho and gene LRT)
need "bb_gene_lrt"
need "bb_estimate_rho_gene"
need "counts are not summed before testing"
need "RHO_MIN     <- {RHO_MIN}"
need "rho_corrected"
need "rho <- max(rho_corrected, RHO_MIN)"
need '| `RHO_MIN` | 0.01 |'
need "rho estimate at the upper boundary"
need "frac_\", alt_label"
need "acat <- function(p) {   # Cauchy combination with equal weights; p is capped"
# --- end Task 6 fix
# --- Task 7 (Step 15, run scripts, submission chain, summary page, README, registration); each string is absent from the pre-Task-7 skill
# --dependency and --parsable are allowlisted above: they are sbatch options (job chaining), not STAR/GATK options.
need "## Step 15"
need "_summary_report.html"
need "relative links"
need '$S/run_01_import_qc.sh'
need '$S/run_02_imbalance.sh'
need '--dependency=afterok:$R1'
need "-n 1 --mem=16G -t 1:00:00"
forbid "-n 8 --mem=16G"
need "# ONLY IF the mouse helper is used"
need "awk -F'\t'"
need "while read -r f; do [ -s"
need "knit_root_dir = '{CWD}'"
need "summary_numbers.tsv"
need "summary_numbers <- summary_tbl"
need "bias_flag = ref_bias_flagged"
need "Verify before finishing"
need "shown prominently"
need "Alignments whose vW tag is 2-7 were removed"
need "## Notes for the assistant"
forbid "-n 2 --mem=8G -t 0:30:00"
forbid "-n 2 --mem=16G"
# README files are repo files, not the skill text
HERE=$(dirname "$0")
grep -qF "Validation status" "$HERE/../README.md" 2>/dev/null || { echo "FAIL: ase-pipeline/README.md missing or lacks: Validation status"; fail=1; }
grep -qF "acceptance run" "$HERE/../README.md" 2>/dev/null && grep -qF "(2026-09-29, DONE)" "$HERE/../README.md" || { echo "FAIL: ase-pipeline/README.md must record the end-to-end acceptance run as done"; fail=1; }
grep -q "/ase-pipeline" "$HERE/../../README.md" 2>/dev/null || { echo "FAIL: root README.md has no /ase-pipeline row"; fail=1; }
grep -qF "Nothing on real biological data" "$HERE/../README.md" 2>/dev/null || { echo "FAIL: ase-pipeline/README.md must state: Nothing on real biological data"; fail=1; }
grep -qF "not exercised" "$HERE/../README.md" 2>/dev/null && grep -qF "extract_mgp_parental_vcf.sh" "$HERE/../README.md" || { echo "FAIL: ase-pipeline/README.md must state what was not exercised (incl. the mouse helper)"; fail=1; }
# --- end Task 7
# --- Task 8 fix (acceptance-run defects D1-D7); each string is absent from the pre-fix skill (f20ba79)
# D1: outbred rho from the central sites with a truncation-corrected likelihood; the naive H0 fit is a diagnostic only
need "bb_estimate_rho_trim <- function"
need "truncated likelihood"
need "rho <- if (is.na(rho_trim)) NA_real_ else max(rho_trim, RHO_MIN)"
need "rho_robust = dplyr::first(if (MODE == \"f1\") rho_corrected else rho_trim)"
need "rho_used, rho_robust, rho_h0_naive, mean_ref_frac, ref_frac_before_wasp, ref_frac_after_wasp,"
need "\`rho_used\`, \`rho_robust\`, \`rho_h0_naive\`, \`mean_ref_frac\`, \`ref_frac_before_wasp\`, \`ref_frac_after_wasp\`, \`bias_flag\`"
need '| `RHO_MIN` | 0.01 | both modes'
need "Outbred mode uses \`bb_estimate_rho_trim\`"
forbid "**Outbred** uses \`rho_h0\`"
forbid "never inflates false positives (conservative)"
forbid "so the estimate is conservative"
forbid "rho_used, rho_h0, mean_ref_frac"
forbid "TRIM_RESULTS_PLACEHOLDER"
# D2: results directory uses YYYY-MM-DD (user data-safety rule); one date placeholder everywhere
need "formatted as \`YYYY-MM-DD\`"
need "\`{RESULTS_DIR}\` = \`{CWD}/results/{TODAY}_{WD_NAME}\`"
need "DATE_TAG    <- \"{TODAY}_{WD_NAME}\""
need "\`{TODAY}_{WD_NAME}_01_import_qc.html\`"
forbid "TODAY_YYMMDD"
forbid "YYMMDD"
# D3: the prep job caches all five containers, Picard included (the array job only checks)
need 'fetch_sif "$PICARD_SIF" "https://depot.galaxyproject.org/singularity/picard:3.1.1--hdfd78af_0"'
need "fetches all five containers"
forbid "block C (STAR, GATK, BCFTOOLS, SAMTOOLS; plus"
forbid "block C (STAR, GATK, BCFTOOLS, SAMTOOLS plus"
# D4: unfiltered (pre-WASP) counts, before/after REF fraction, information only; tables read per sample, never globbed
need 'UNF_TABLE="$R/ase_counts/$SAMPLE.unfiltered.table"'
need 'rm -f "$UNF_TABLE"'
need 'check_table "$TABLE"'
need 'check_table "$UNF_TABLE"'
need "ref_frac_before_wasp"
need "ref_frac_after_wasp"
need "reported for information, not a gate"
need 'stopifnot(setequal(unique(sites_raw$sample), samples$sample))'
need '!grepl("\\.unfiltered\\.table$", tbl_files)'
forbid 'list.files(COUNTS_DIR, pattern = "\\.table$")), samples$sample)'
# D5: array log pattern
need "#SBATCH -o {RESULTS_DIR}/logs/align_count_f1_%A_%a.out"
need "#SBATCH -o {RESULTS_DIR}/logs/align_wasp_count_%A_%a.out"
# D6: ready-made --bind lines for the Rmd run scripts
need '`{GTF_DIR}` inside `{CWD}`: `singularity exec --bind {CWD} {R_SIF} \`'
need '`{GTF_DIR}` outside `{CWD}`: `singularity exec --bind {CWD},{GTF_DIR} {R_SIF} \`'
# D7: plain parental VCF accepted; explicit fallback strain names
need "a plain \`.vcf\` is accepted"
need "\`{STRAIN_A}\` = \`C57BL_6NJ\` and \`{STRAIN_B}\` = \`A_J\`"
forbid "biallelic SNP VCF (bgzipped and tabix-indexed) where"
# README (repo file) follows the same fixes
RD="$HERE/../README.md"
grep -qF "central sites" "$RD" || { echo "FAIL: ase-pipeline/README.md must describe the trimmed outbred rho (central sites)"; fail=1; }
! grep -qF "Outbred uses the H0-based estimate" "$RD" || { echo "FAIL: ase-pipeline/README.md still says outbred uses the H0-based estimate"; fail=1; }
! grep -qF "YYMMDD" "$RD" || { echo "FAIL: ase-pipeline/README.md still uses YYMMDD"; fail=1; }
grep -qF "unfiltered.table" "$RD" || { echo "FAIL: ase-pipeline/README.md must list the outbred unfiltered.table"; fail=1; }
# --- end Task 8 fix
# --- Final-review fix wave (I1-I6 + cheap minors); each string is absent from the reviewed skill/README/spec (e6e3da9)
# I1: F1 prep checks that the strains differ; strain A must be the reference-consistent one
need "A_ALT_SITES=stop"
need "MIN_PARENTAL_SITES=1000"
need 'if (ga == "0/0" && gb == "1/1") c = "keep"; else if (ga == "1/1" && gb == "0/0") c = "a_alt"'
need 'dropped $N_SAME (same genotype in both strains), $N_HET (heterozygous or missing), $N_AALT'
need 'sites have $A = 1/1 and $B = 0/0: strain A carries the non-reference allele there'
need '[ "$N_SITES" -ge "$MIN_PARENTAL_SITES" ]'
need "**Sites-only VCF** (no sample columns): every biallelic SNP is taken on trust"
need "F1 mode needs one parent to be the reference strain"   # Task 6 fix M1: was "F1 mode in Stage 1 needs ..."
need "the VCF was not called against this assembly"
forbid "F1 mode needs strain A = the reference strain, REF = strain A, ALT = strain B (and matching contig names)"
# I2: masked STAR index keyed on the mask; reused only when index_key.txt matches; array jobs read the recorded path
need "MASK_KEY=\$(grep -v '^#' \"\$REF_DIR/masked_sites.vcf\" | cut -f1,2,4,5 | md5sum | cut -c1-8)"
need 'STAR_INDEX="{STAR_INDEX}_$MASK_KEY"; INDEX_KEY="$MASK_KEY"'
need '[ "$(cat "$STAR_INDEX/index_key.txt" 2>/dev/null)" = "$INDEX_KEY" ]'
need 'echo "$INDEX_KEY" > "$STAR_INDEX/index_key.txt"'
need 'STAR_INDEX=$(cat "$R/reference/star_index_path.txt" 2>/dev/null)'
need '--genomeDir "$STAR_INDEX"'
forbid '--genomeDir "{STAR_INDEX}"'
# I3: the mouse helper never reads the FASTA or .fai; the FASTA contig guard is in the prep job
need "**The helper needs neither the FASTA nor its \`.fai\`**"
need "sets \`{N_CHROM}\` = 20"
need "the \`add_bind\` line is reduced to \`add_bind \"{CWD}\"\`, because the helper reads no genome file"
need "right after block C, before any command that reads the FASTA or the GTF"
# masked genome: genotyped mask VCF + consensus -s MASK (a sites-only VCF with -H A was applied in only 4 of 10 identical runs)
need 'bcftools consensus -s MASK -f "$FASTA" "$REF_DIR/masked_sites.vcf.gz"'
need 'third(toupper($4),toupper($5)),".",".",".","GT","1/1"'
forbid "bcftools consensus -H A"
need "which downloads/decompresses the FASTA, writes the .fai and compares the parental VCF's contigs with it"
need 'contig mismatch: $N_OFF of $N_BI parental SNPs lie on contigs that are not in $FASTA'
forbid "grep -m1 '^>' \"{FASTA_PATH}\""
forbid "taken from the \`.fai\`"
# I4: outbred contig guard in prep_genotypes.sh before any alignment
need 'contig mismatch for individual $IND: $N_OFF of $N_ALL heterozygous sites lie on contigs that are not in $FASTA'
need "then **block R first** (the contig guard below needs \`\$FASTA.fai\`)"
need "**Contig names must match the FASTA.**"
# I5: count tables read with explicit column classes (contig as character)
need 'colClasses = c(contig = "character", variantID = "character",'
need "**Every table is read with explicit column classes**"
forbid 'read.delim(f, stringsAsFactors = FALSE, check.names = FALSE), error'
# I6: SNP-dependence limitation (Step 13, Step 14, Rmd 02 gene table, README): superseded by the thinned F1 gene test, see Stage 2 Task 5 below
# minors: GTF bind, fetch_sif calls, wrapper sets, paths with spaces, dropping a failed sample, README wording, spec
need 'add_bind "$(dirname "{FASTA_PATH}")"; add_bind "$(dirname "{GTF_PATH}")"'
need "must **call** it once per container it uses"
need 'call** it for each container the script uses, on the lines right after block C (see "Per-script fetch lines" in Step 10): `fetch_sif "$STAR_SIF"; fetch_sif "$GATK_SIF"; fetch_sif "$SAMTOOLS_SIF"; fetch_sif "$PICARD_SIF"`'
need '| `align_count_f1.sh`, `align_wasp_count.sh` | `star`, `gatk`, `samtools`, `picard` |'
need "choose them by the function name at the start of the line (never by position)"
need "**Paths with spaces or commas are not supported.**"
need "no space or comma in any \`fastq_1\` / \`fastq_2\` path"
need "**Dropping a failed sample.**"
need "remove that sample's row from \`{SAMPLES_CSV}\`"
grep -qF "the interactive dialogue of Steps 0-9" "$RD" || { echo "FAIL: ase-pipeline/README.md must say the interactive dialogue was not exercised"; fail=1; }
! grep -qF "a fresh run of the installed skill" "$RD" || { echo "FAIL: ase-pipeline/README.md still overstates the acceptance run"; fail=1; }
SPEC="$HERE/../../docs/superpowers/specs/2026-09-29-ase-pipeline-design.md"
[ -s "$SPEC" ] || { echo "FAIL: spec missing: $SPEC"; fail=1; }
! grep -qF "YYMMDD" "$SPEC" || { echo "FAIL: spec still uses YYMMDD"; fail=1; }
! grep -qF "gene counts are exact sums" "$SPEC" || { echo "FAIL: spec still describes the F1 gene test as summed counts"; fail=1; }
grep -qF "likelihood-ratio test (H0 p = 0.5 vs p free, 1 df)" "$SPEC" || { echo "FAIL: spec must describe the F1 gene LRT"; fail=1; }
# --- end final-review fix wave


# --- Stage 2 statistics block; each string is absent from the Stage 1 skill (c073dfb)
need "chrom_class <- function(contig)"
need "thin_snps <- function(pos, depth, window)"
need "bb_glm_fit <- function(y, n, X, rho, bound = 15)"
need "bb_glm_lrt <- function(y, n, X, drop_cols, rho, bound = 15)"
need "bb_moment_phi <- function(units, max_iter = 25, tol = 1e-6, bound = 15)"
need "ase_glm_test <- function(units, tests, rho_min, min_df_unit = 2, bound = 15)"
need "ase_paired_test <- function(d, rho_min, min_individuals = 2, bound = 15)"
need "bb_pair_phi_trim <- function(y, n, cell, keep = 0.9, max_iter = 50)"
need '| `THIN_BP` | 500 |'
need "**Stage 2 models (Rmd 03 and 04): \`ase_glm_test\`.**"
need "**Outbred differential: \`ase_paired_test\`.**"
need "**SNPs sharing read pairs.**"
[ -s "$HERE/r/test_ase_stats_stage2.R" ] || { echo "FAIL: missing tests/r/test_ase_stats_stage2.R"; fail=1; }
# --- end Stage 2 statistics block

# --- Stage 2 statistics review fixes; each check fails on the 5fdcae4 skill and tests
need 'thin_snps: ", sum(bad), " SNP(s) with a missing position or depth'
need 'bb_glm_lrt: tested column(s) not in the design: '
need 'for (p in c("openxlsx","tidyverse","GenomicRanges","rtracklayer"))'
forbid 'c("aod","lme4","openxlsx"'
need "The paired SNP sizes and power are over the tested SNPs only"
need "excludes the SNP-dense genes"
RT="$HERE/r/test_ase_stats_stage2.R"
grep -qF 'tst <- x$status == "ok"; sg <- tst & x$p < 0.05' "$RT" || { echo "FAIL: paired sizes/power must be computed over tested SNPs only"; fail=1; }
grep -qF 'UNTHINNED_F1_GENE_SIZE' "$RT" && ! grep -qF 'TASK5_INPUT' "$RT" || { echo "FAIL: the unit test must print UNTHINNED_F1_GENE_SIZE (no process token)"; fail=1; }
grep -qF '#SBATCH -N 1 -n 1 -c 8 --mem=16G' "$HERE/r/run_stats_tests.sh" || { echo "FAIL: run_stats_tests.sh must request -N 1 -n 1 -c 8"; fail=1; }
! grep -qF '`aod`, `lme4`, `openxlsx`' "$RD" || { echo "FAIL: ase-pipeline/README.md still lists aod and lme4 as required"; fail=1; }
! grep -qF 'require `aod` and `lme4`' "$SPEC" || { echo "FAIL: spec still requires aod and lme4"; fail=1; }
# --- end Stage 2 statistics review fixes

# --- Stage 2 Task 3 (Rmd 03, Step 16); each string is absent from the 3aa9286 skill
need "## Step 16 — Rmd 03: reciprocal F1"
need "{CWD}/{TODAY}_{WD_NAME}_03_reciprocal.Rmd"   # "_03_reciprocal.Rmd" alone already occurs in Step 14
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
# chrom_class: only contigs of the GTF gene map are tested; alt, random and unplaced contigs count as autosomes (said in the Rmd text)
need "Only the contigs present in the GTF gene map are tested"
need "alt, random and unplaced contigs are treated as autosomes"
# --- end Stage 2 Task 3

# --- Stage 2 Task 3 review fixes; each check fails on the 70ba7ef skill and tools
# I1: condition confounded with cross direction stops; per-gene non-estimable designs are listed; no testable gene stops
need 'if (qr(Xd)$rank < ncol(Xd))'
need "condition is confounded with cross direction"
need "design not estimable (condition confounded with direction in the covered samples)"
need "no gene could be tested (fit status: "
# M1: strain names checked against the Rmd 01 checkpoint
need 'STRAIN_A = identical(ck$constants$STRAIN_A, STRAIN_A), STRAIN_B = identical(ck$constants$STRAIN_B, STRAIN_B)'
# M8: b1 is one effect common to all conditions
need "b1 is one parent-of-origin effect common to all conditions"
forbid "b0 and b1 are averages over the conditions"
# M9: nothing left after the gene map names the gene map, not the coverage
need "no autosomal SNP in a single GTF gene is left after the gene map"
# M2/M3 (evaluator) and M4 (render tool) are repository files
SYN="${SYN_DIR:-$HERE/synthetic}"   # SYN_DIR: only to run these checks on other copies of the tools
grep -qF 'summary_numbers_reciprocal.tsv' "$SYN/evaluate_stage2.R" || { echo "FAIL: evaluate_stage2.R must check summary_numbers_reciprocal.tsv"; fail=1; }
grep -qF 'direction strings match the signs' "$SYN/evaluate_stage2.R" || { echo "FAIL: evaluate_stage2.R must check the direction strings"; fail=1; }
! grep -qF '>= 500)' "$SYN/evaluate_stage2.R" || { echo "FAIL: evaluate_stage2.R must read THIN_BP from the values file"; fail=1; }
grep -qF "trap 'rm -f" "$SYN/render_from_skill.sh" || { echo "FAIL: render_from_skill.sh must remove its temporary files on failure"; fail=1; }
grep -qF 'no statistics block between the ase-stats markers' "$SYN/render_from_skill.sh" || { echo "FAIL: render_from_skill.sh must stop without the stats markers"; fail=1; }
# --- end Stage 2 Task 3 review fixes

# --- Stage 2 Task 4 (Rmd 04, Step 17); each check fails on the f083a98 skill and tools
need "## Step 17 — Rmd 04: differential ASE between conditions"
need "{CWD}/{TODAY}_{WD_NAME}_04_differential.Rmd"   # "_04_differential.Rmd" alone already occurs in Step 14
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
# Task 3 review fixes repeated in Rmd 04 (the Rmds are self-contained): these strings must occur in Step 16 AND Step 17
need_n() { [ "$(grep -cF -- "$2" "$SKILL")" -ge "$1" ] || { echo "FAIL: fewer than $1 lines with: $2"; fail=1; }; }
need_n 2 'STRAIN_A = identical(ck$constants$STRAIN_A, STRAIN_A), STRAIN_B = identical(ck$constants$STRAIN_B, STRAIN_B)'
need_n 2 "no autosomal SNP in a single GTF gene is left after the gene map"
need_n 2 "design not estimable (condition confounded with direction in the covered samples)"
need_n 3 "alt, random and unplaced contigs are treated as autosomes"
need "no gene could be tested in any contrast (fit status: "
need "no SNP could be tested in any contrast"
# the outbred not_tested count must not count a SNP twice (listed in not_tested AND with p NA)
need "not_tested = length(unique(c(nt_ids, ids[is.na(pcol)])))"
# evaluator: the differential branches exist and check the summary file
grep -qF 'what == "differential_outbred"' "$SYN/evaluate_stage2.R" || { echo "FAIL: evaluate_stage2.R lacks the differential_outbred branch"; fail=1; }
grep -qF 'summary_numbers_differential.tsv' "$SYN/evaluate_stage2.R" || { echo "FAIL: evaluate_stage2.R must check summary_numbers_differential.tsv"; fail=1; }
# --- end Stage 2 Task 4

# --- Stage 2 Task 4 review fixes; each check fails on the 14d43ef skill and tools
# I1: a contrast with nothing tested gets its own message, one Summary row per contrast x level, figure counts with a zero default
need "nothing tested in contrast"
need "for (ct in contrasts)"
need 'fig_n <- vapply(names(sum_n), function(ct) sum(fig$sig[fig$contrast == ct]), integer(1))'
need '"phi_pair NA"'
# m1: the outbred label means all TESTED individuals (skill, evaluator)
need "(all tested individuals)"
forbid "(all individuals)"
! grep -qF '(all individuals)' "$SYN/evaluate_stage2.R" || { echo "FAIL: evaluate_stage2.R still uses the label '(all individuals)'"; fail=1; }
# m2: with the d term the fractions are direction-averaged
need "frac_A_ref and frac_A_test are direction-averaged"
# m4: the two Rmds count excluded chromosome rows differently, and say so
need "Rmd 03 counts the sample-site rows of SNPs in exons of a single GTF gene"
need "Rmd 04 counts every sample-site row on X, Y or MT"
# m5: a wrong REF_CONDITION needs Rmd 04 written again, not Rmd 01
need "Rmd 01 need not be re-rendered"
# --- end Stage 2 Task 4 review fixes

# --- Stage 2 Task 5 (Rmd 02 F1 gene test on thinned SNPs); each check fails on the da225c8 skill, README and spec
need 'note = "thinned: at most one SNP per THIN_BP window (exon coordinates), so no read pair is counted twice; approximate (ignores alternative splicing)"'
need 'snp_gene$thin_keep'
need "**Thinned SNPs (F1 gene LRT and \`bb_estimate_rho_gene\`).**"
forbid "Planned for Stage 2 (not implemented): thin each gene's SNPs"
need 'n_snps_used = nrow(g)'
need "F1 gene test: at most one SNP per THIN_BP window (exon coordinates) is used, so no read pair is counted twice (approximate: alternative splicing can bring distant exons into one fragment)."
need 'g_thin <- sg_keep$gene_id[match(paste(d$contig, d$position), paste(sg_keep$contig, sg_keep$position))]'
forbid 'g <- sg_keep$gene_id'   # the per-SNP rho stays on all SNPs (controller ruling 1)
# controller ruling 2: the gene LRT has its own rho_gene: thinned -> unthinned free-mean -> H0-based
need 'rt <- collect_warnings(bb_estimate_rho_gene(d$alt_n, d$total, g_thin))'
need 'gsrc <- "unthinned SNPs, free-mean per gene, bias-corrected (fallback: fewer than 5 genes with 2 or more thinned SNPs)"'
need 'rho_gene_source = dplyr::first(rho_gene_source)'
need "The F1 gene LRT has its own \`rho_gene\`, with this fallback chain"
need "**Fallback chain of \`rho_gene\`.**"
need "the Stage 1 figure of 24 of 24 planted gene x sample tests no longer applies"
grep -qF "frag_fallback_run <- function(seed" "$RT" && grep -qF "RHO_GENE_FALLBACK" "$RT" || { echo "FAIL: test_ase_stats_stage2.R must simulate the rho_gene fallback"; fail=1; }
grep -qF "costs power in genes whose SNPs are close together" "$RD" || { echo "FAIL: README must state the power cost of the thinned F1 gene test"; fail=1; }
need 'snp_by_gene_f1 <- dplyr::inner_join(snp, dplyr::filter(snp_gene, thin_keep), by = c("contig", "position"))'
need "the depth is pooled over all samples of the project"
need '| `THIN_BP` | 500 | F1 gene-level test of Rmd 02 and the gene-level tests of Rmd 03 and 04'
forbid "SNP counts treated as independent"
forbid "**Limitation (F1):** SNP counts that share read pairs are treated as independent"
grep -qF "thinned to one SNP per" "$RD" || { echo "FAIL: README must describe the thinned F1 gene test"; fail=1; }
! grep -qF "SNPs that share read pairs are treated as independent" "$RD" || { echo "FAIL: README still says the F1 gene test treats SNP counts as independent"; fail=1; }
! grep -qF "treated as independent, a documented limitation" "$SPEC" || { echo "FAIL: spec still says the F1 gene test treats SNP counts as independent"; fail=1; }
# --- Stage 2 Task 5 review fixes
# I1: the per-SNP rho labels and snp_annot keep the Stage 1 findOverlaps-on-exons table (GTF exon order, distinct); the
#     exon-coordinate table for thinning is separate (sg_exon) and only adds the thin_keep column
need 'hits <- GenomicRanges::findOverlaps(gr, ex)'
need 'gene_id = ex$gene_id[S4Vectors::subjectHits(hits)], stringsAsFactors = FALSE) %>% dplyr::distinct()'
need 'thin_depth <- dplyr::inner_join(sites, sg_exon, by = c("contig", "position"))'
# I2: the fallback caveat, in Step 13 (rho_gene_source), Step 14, the run-time NOTE, and the README
need_n 3 "the fallback with the unthinned dispersion has not been verified when the true overdispersion is above RHO_MIN; because SNPs that share reads push that estimate low, the gene test may be anti-conservative in that case"
grep -qF "the fallback with the unthinned dispersion has not been verified when the true overdispersion is above RHO_MIN" "$RD" || { echo "FAIL: README must state that the rho_gene fallback is not verified above RHO_MIN"; fail=1; }
# I3: the USE sites of the behaviour (these pass on 58805da; each one is proved by a mutation that removes the behaviour)
need 'gene <- dplyr::bind_rows(lapply(split(snp_by_gene_f1, snp_by_gene_f1$sample), function(d) {'
need 'r <- if (is.na(g$rho_gene[1])) list(p = NA_real_, phat = NA_real_) else bb_gene_lrt(g$alt_n, g$total, g$rho_gene[1])'
need 'd$rho_gene <- rho_gene   # F1 gene LRT only'
need '} else if (!is.na(rho_corrected)) {'
need 'dplyr::group_by(gene_id) %>% dplyr::mutate(thin_keep = thin_snps(exon_pos, depth, THIN_BP)) %>% dplyr::ungroup()'
need 'g <- snp_gene$gene_id[match(paste(d$contig, d$position), paste(snp_gene$contig, snp_gene$position))]'
need 'r1 <- collect_warnings(bb_estimate_rho_gene(d$alt_n, d$total, g))'
need 'd$p <- mapply(bb_p_safe, d$alt_n, d$total, MoreArgs = list(rho = rho))'
forbid 'dplyr::mutate(thin_keep = TRUE'
# minors: Step 15 rho_used, gene_level label, clear stop, outbred Summary, Step 13 wording, README validation line
need '`rho_used` is the overdispersion the SNP tests used (the F1 gene tests use `rho_gene`, see Step 13)'
need '"per-gene LRT on thinned SNP-level counts (at most one SNP per THIN_BP window; shared strain fraction)"'
need "have SNPs in GTF genes but none of the SNPs kept by the thinning"
need 'if (MODE != "f1") summary_tbl$rho_gene_source <- NULL'
forbid "saying that the test is thinned. F1 gene test:"
! grep -qF "F1: 24 of 24 planted gene x sample tests significant with the correct direction, 0 of 42 null false positives." "$RD" || { echo "FAIL: README validation line still gives the Stage 1 24 of 24 as current"; fail=1; }
grep -qF "the three lost tests are genes thinned to one SNP" "$RD" || { echo "FAIL: README validation line must give the thinned F1 result"; fail=1; }
# --- end Stage 2 Task 5 review fixes
# --- end Stage 2 Task 5
# --- Stage 2 Task 6 (wizard, chain, summary)
need "\`cross_direction\` is exactly \`AxB\` or \`BxA\`"
need "Which condition is the reference (baseline)?"
need "Store it as \`{REF_CONDITION}\`"
need "condition and cross direction are confounded"
need "at least 2 individuals sampled in both"
need 'R3=$(sbatch -p bcc --parsable --dependency=afterok:$R1 $S/run_03_reciprocal.sh)'
need 'R4=$(sbatch -p bcc --parsable --dependency=afterok:$R1 $S/run_04_differential.sh)'
need "### Reciprocal F1 section"
need "### Differential ASE section"
need "offered only when \`{MODE}\` = \`outbred\`; available in a later stage"
need "pasted verbatim into Rmd 02, 03 and 04"
forbid "Stage 1 implements only the always-on analysis"
forbid "**Stage 1 only.**"
# additions beyond the brief (review context of Tasks 3-5): Rmd 03 confounded stop, per-contrast checks,
# missing-script guard, chain echo, nothing-tested contrast rows, thinning depends on the sample set, X/Y/MT in Rmd 04
need "at least one condition must hold samples of both directions"
need "checked for every contrast against \`{REF_CONDITION}\`"
need 'MISSING=$(for f in $NEED; do [ -s "$S/$f" ] || echo "$f"; done)'
need "a contrast with \`tested\` = 0"
need "removing a sample can change which SNP a dense gene keeps"
need "X, Y and MT are not tested in either mode"
forbid "The summary page is written from \`summary_numbers.tsv\`, not from R."
# --- end Stage 2 Task 6
# --- Stage 2 Task 6 review fixes; each check fails on the ec81a61 skill/README
# ruling: Rmd 03/04 depend on Rmd 01 (they read only ase_checkpoint.rds), not on Rmd 02
forbid '--dependency=afterok:$R2 $S/run_03_reciprocal.sh'
forbid '--dependency=afterok:$R2 $S/run_04_differential.sh'
need "# 5. Rmd 03 only when reciprocal F1 was selected (Step 7): after Rmd 01"
grep -qF 'Rmd 03 / Rmd 04 (`afterok` on Rmd 01' "$RD" || { echo "FAIL: README chain must put Rmd 03 / Rmd 04 after Rmd 01"; fail=1; }
# I1: wait on every submitted id, stop on never-satisfiable dependencies, bounded
forbid 'squeue -h -j $LAST'
need 'JOBS=$(cut -f2 "$IDS" | paste -sd, -)'
need "Q=\$(squeue -h -j \"\$JOBS\" -o '%i %r' 2>/dev/null)"
need "grep -qv ' DependencyNeverSatisfied\$'"
need 'scancel $STUCK'
need '[ $n -ge 240 ]'
# M3: every sbatch id checked; "not selected" distinct from a failed submission
need 'got() { [ -n "$2" ] || die'
need 'R3="not selected"; R4="not selected"'
forbid 'Rmd03 ${R3:-none}'
# I2: outbred menu precondition over any pair of conditions; per-contrast check after the reference is chosen
need "at least one pair of conditions with at least 2 individuals sampled in both"
# M1, M2, M4, M5
forbid "not supported in Stage 1"
need "xargs -r scancel"
need "control, ctrl, untreated, wt or vehicle"
need "choose another reference condition first"
need "grep -o '## Contrasts not tested:[^<]*'"
# --- end Stage 2 Task 6 review fixes
# --- Stage 2 acceptance fix wave (T1-T3); each check fails on the 32033d4 skill
# T1: block C code fence holds no fetch_sif call (a literal copy must work for every script); calls live in a labelled list
need "**Per-script fetch lines**"
forbid "# ... one fetch_sif line per container this script uses"
blockC_calls=$(awk '/^### Shared block C/ {b=1} b && /^```bash$/ {f=1; next} f && /^```$/ {exit} f && /^fetch_sif "\$/ {n++} END {print n+0}' "$SKILL")
[ "$blockC_calls" = 0 ] || { echo "FAIL: block C code fence contains $blockC_calls fetch_sif call line(s)"; fail=1; }
blockC_def=$(awk '/^### Shared block C/ {b=1} b && /^```bash$/ {f=1; next} f && /^```$/ {exit} f && /^fetch_sif\(\) \{/ {n++} END {print n+0}' "$SKILL")
[ "$blockC_def" = 1 ] || { echo "FAIL: block C code fence must hold the fetch_sif definition"; fail=1; }
# T2: Step 11 points to Step 15, no duplicate submission block
forbid 'P=$(sbatch -p bcc --parsable {RESULTS_DIR}/scripts/prep_f1_reference.sh)'
need "do not submit these scripts by hand"
# T3: mode-neutral submit_chain.sh
need 'f1)      PREP_SCRIPT=prep_f1_reference.sh; ARRAY_SCRIPT=align_count_f1.sh ;;'
need 'outbred) PREP_SCRIPT=prep_genotypes.sh;    ARRAY_SCRIPT=align_wasp_count.sh ;;'
need 'P=$(sbatch -p bcc --parsable $DEP $S/$PREP_SCRIPT); got $PREP_SCRIPT "$P"'
need 'A=$(sbatch -p bcc --parsable --dependency=afterok:$P $S/$ARRAY_SCRIPT); got $ARRAY_SCRIPT "$A"'
need 'case " $ANALYSES " in *reciprocal*)'
need 'case " $ANALYSES " in *differential*)'
need 'NEED="$PREP_SCRIPT $ARRAY_SCRIPT run_01_import_qc.sh run_02_imbalance.sh"'
need "Do not edit the block except for the mouse helper"
forbid 'got prep_f1_reference.sh'
forbid "Delete line 5 unless"
forbid "F1 = prep_f1_reference.sh, outbred = prep_genotypes.sh"
# --- end Stage 2 acceptance fix wave
# --- Stage 2 README and registration; each check fails on the 212361b README and root README
RD="$HERE/../README.md"
rneed() { grep -qF -- "$1" "$RD" || { echo "FAIL: README must contain: $1"; fail=1; }; }
rforbid() { ! grep -qF -- "$1" "$RD" || { echo "FAIL: README must not contain: $1"; fail=1; }; }
grep -qF "(Stages 1 and 2)" "$RD" || { echo "FAIL: README title must say Stages 1 and 2"; fail=1; }
grep -qF "Reciprocal F1 (Rmd 03)" "$RD" && grep -qF "Differential ASE (Rmd 04)" "$RD" || { echo "FAIL: README must describe Rmd 03 and Rmd 04"; fail=1; }
grep -qF "summary_numbers_reciprocal.tsv" "$RD" && grep -qF "summary_numbers_differential.tsv" "$RD" || { echo "FAIL: README outputs must list the Stage 2 TSVs"; fail=1; }
grep -qF "direction-free" "$RD" || { echo "FAIL: README must explain the direction-free outbred test"; fail=1; }
grep -qF "X, Y and MT" "$RD" || { echo "FAIL: README must state the X/Y/MT exclusion"; fail=1; }
grep -qF "Stage 2 synthetic acceptance run: DONE (2026-09-30)" "$RD" && ! grep -qF "synthetic acceptance run: PENDING" "$RD" || { echo "FAIL: README must state Stage 2 acceptance DONE (2026-09-30)"; fail=1; }
! grep -qF "Stage 1 covers **per-sample allelic imbalance only**" "$RD" || { echo "FAIL: README still says Stage 1 only"; fail=1; }
grep -q "/ase-pipeline.*reciprocal F1" "$HERE/../../README.md" || { echo "FAIL: root README row must mention reciprocal F1"; fail=1; }
! grep -q "/ase-pipeline.*Stage 1: per-sample imbalance only" "$HERE/../../README.md" || { echo "FAIL: root README row still says Stage 1 only"; fail=1; }
rforbid "Stage 1 covers only per-sample imbalance"
rforbid "(only the per-sample analysis is implemented in Stage 1)"
rforbid "not supported in Stage 1"
# what Stage 2 adds and its requirements
rneed "**What Stage 2 adds**"
rneed "must be exactly \`AxB\` or \`BxA\`"
rneed "reference condition"
rneed "| 16 | Rmd 03: reciprocal F1 |"
rneed "| 17 | Rmd 04: differential ASE |"
rneed "ASE_reciprocal.xlsx"
rneed "ASE_differential.xlsx"
rneed "No mixed model (lme4) is used"
# measured validation numbers (source: the final unit-test log and the evaluator runs)
rneed "173 checks ok, 0 failed"
rneed "0.0370 (strain) / 0.0367 (parent of origin)"
rneed "tested fraction 0.4244 / 0.6961 / 0.8434 / 0.9646"
rneed "at least 0.70 and at least 90 percent of the oracle power"
rneed "0.1450 (0.2084 in SNP-dense genes), the thinned test that Rmd 02 now uses 0.0497"
rneed "5 of 5 planted strain genes and 5 of 5 planted parent-of-origin genes found with the correct sign, 0 of 47 null genes"
rneed "4 of 4 planted genes found with the correct sign"
rneed "5 of 5 planted genes found, 0 of 3 imbalanced-but-unchanged genes and 0 of 52 null genes"
rneed "21 of 24 planted, 0 of 42 null, 65 significant SNPs unchanged"
# limitations stated plainly
rneed "with 30 percent changed (phase-heterogeneous) it is 0.6698 against an oracle of 0.9163"
rneed "null sizes are about 0.036 at a true dispersion of 0.02"
rneed "phASER (Stage 3) is not implemented"
rneed "unplaced contigs are treated as autosomes"
rneed "adding samples or conditions can change which SNP is kept"
rneed "\`sacct\` is unavailable on the test cluster"
rneed "The mouse helper was only partly exercised"
# no internal process labels in shipped text
! grep -qE "Task [0-9]" "$RD" "$HERE/../ase-pipeline.md" || { echo "FAIL: README or skill contains an internal 'Task N' label"; fail=1; }
# Stage 2 acceptance results recorded in the README (fail on the 32033d4 README)
rneed "0 of 47 null genes called in each test"
rneed "0 of 9 calls among strain / parent-of-origin genes"
rneed "28 of 28 count tables non-empty"
rneed "66 s in F1 mode, 77 s in outbred mode"
rneed "planted strain genes detected in 50 of 60 gene x sample tests"
rneed "null false positives 1 of 564"
rneed "SYN000047"
rneed "the cause was not isolated"
rneed "say nothing about real data"
# --- end Stage 2 README and registration


# --- Stage 2 final-review fix wave (fail on 6416df4)
# I1: one spelling of the analysis words, menu mapping in Step 7, guard in submit_chain.sh
forbid "per_sample"
need 'menu number 1 is the word `per-sample`, 2 is `reciprocal`, 3 is `differential`'
need 'holds the words, separated by spaces, never the menu numbers'
need 'case "$ANALYSES" in "") die "ANALYSES is empty: write the Step 7 words'
need 'for w in $ANALYSES; do case "$w" in per-sample|reciprocal|differential) ;; *) die "ANALYSES contains '"'"'$w'"'"'; the accepted words are per-sample, reciprocal and differential'
# M1: the Rmd 03 / Rmd 04 submission must sit directly under its case line (the NEED line cannot satisfy this)
pair_need() { awk -v a="$1" -v b="$2" 'index($0,b) && p {f=1} {p = index($0,a) > 0} END {exit !f}' "$SKILL" || { echo "FAIL: line not directly under its case line: $2"; fail=1; }; }
pair_need 'case " $ANALYSES " in *reciprocal*)' '  R3=$(sbatch -p bcc --parsable --dependency=afterok:$R1 $S/run_03_reciprocal.sh); got run_03_reciprocal.sh "$R3" ;;'
pair_need 'case " $ANALYSES " in *differential*)' '  R4=$(sbatch -p bcc --parsable --dependency=afterok:$R1 $S/run_04_differential.sh); got run_04_differential.sh "$R4" ;;'
# M3: one sample per animal in F1 mode, duplicated individual stops
need "one sample per animal"
need "merge the FASTQ files of the technical replicates of an animal before the wizard"
need "in F1 mode every \`individual\` value must be unique"
need "F1 mode needs one sample per animal; these individual values occur more than once:"
need "awk -F, 'NR > 1 {c[\$6]++} END {for (i in c) if (c[i] > 1) print i}' {SAMPLES_CSV}"
# M5: no overwrite of an existing chain_job_ids.tsv; re-run behaviour stated
need '[ ! -e "$IDS" ] || die "a chain has been submitted from this results directory'
need "check squeue, then remove or rename"
need "mv {RESULTS_DIR}/logs/chain_job_ids.tsv {RESULTS_DIR}/logs/chain_job_ids_previous.tsv"
need "never delete an existing results directory"
# M4: what to do after the wait timeout
need "On TIMEOUT (exit 2)"
need "do not cancel them"
# M6: X is not set aside in Rmd 01 / Rmd 02
need "Rmd 01 and Rmd 02 do not set X aside"
need "male F1 animals carry one X"
# M7: runtime wording from the measurement
forbid "projects under 3 h"
need "6.4 ms per unit (32 s for 5000 units"
need "projects 0.11 h for 60,000 genes"
# M9
forbid "only when all individuals agree"
need "only when all tested individuals agree"
# M10
forbid "(it rewrites \`chain_job_ids.tsv\`"
need "the new \`submit_rmds.sh\` writes a new \`chain_job_ids.tsv\`"
# M11: process traces
forbid "debug job"
forbid "Stage 1 synthetic acceptance"
forbid "the Stage 2 beta-binomial"
! grep -qE 'job [0-9]{6,}' "$SKILL" || { echo "FAIL: skill mentions a job id"; fail=1; }
! grep -qE 'Singularity\. [a-z]' "$HERE/../../README.md" || { echo "FAIL: root README row: lower-case sentence start"; fail=1; }
# M13
need "PAT='error|halted|due to time limit|cancelled|out of memory|oom-kill'"
need 'grep -h -m1 -iE "$PAT"'
forbid "-iE 'error|halted'"
# README: M6, M8
rneed "do not set X aside"
rneed "verified with a stub scheduler, not live"
rneed "manually adapted chain"
rforbid "the skill's own \`submit_chain.sh\` and \`wait_chain.sh\` ran"
# --- end Stage 2 final-review fix wave

# --- Stage 3 statistics (hap_gene_test); each string is absent from the Stage 2 skill (de43444)
need "hap_gene_test <- function(a, b, sample, rho_min, min_genes = 20)"
need "per sample rho_used = max(rho_cohort, rho_own, rho_min)"
need "pasted verbatim into Rmd 02, Rmd 03, Rmd 04 and Rmd 05"
need "**phASER gene test (Rmd 05): \`hap_gene_test\`.**"
[ -s "$HERE/r/test_ase_stats_stage3.R" ] || { echo "FAIL: missing tests/r/test_ase_stats_stage3.R"; fail=1; }
grep -qF 'test_ase_stats_stage3.R' "$HERE/r/run_stats_tests.sh" || { echo "FAIL: run_stats_tests.sh does not run the Stage 3 tests"; fail=1; }
# --- end Stage 3 statistics

[ $fail -eq 0 ] && echo "PASS" || exit 1
