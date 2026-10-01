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
for f in phaser_help.txt phaser_gene_ae_help.txt; do
  [ -s "$FIX/$f" ] && tool_flags="$tool_flags
$(grep -oE -- '--[a-z][a-z0-9_]*' "$FIX/$f" | sed 's/^--//' | sort -u)"
done
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
need "pasted verbatim into Rmd 02, 03, 04 and 05"   # Stage 3 Task 4 (M11): Rmd 05 added to the Notes line
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
grep -qF "Reciprocal F1 (Rmd 03)" "$RD" && grep -qF "Differential ASE (Rmd 04)" "$RD" || { echo "FAIL: README must describe Rmd 03 and Rmd 04"; fail=1; }
grep -qF "summary_numbers_reciprocal.tsv" "$RD" && grep -qF "summary_numbers_differential.tsv" "$RD" || { echo "FAIL: README outputs must list the Stage 2 TSVs"; fail=1; }
grep -qF "direction-free" "$RD" || { echo "FAIL: README must explain the direction-free outbred test"; fail=1; }
grep -qF "X, Y and MT" "$RD" || { echo "FAIL: README must state the X/Y/MT exclusion"; fail=1; }
grep -qF "Stage 2 synthetic acceptance run: DONE (2026-09-30)" "$RD" && ! grep -qF "Stage 2 synthetic acceptance run: PENDING" "$RD" || { echo "FAIL: README must state Stage 2 acceptance DONE (2026-09-30)"; fail=1; }
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
need "PAT='error|exception|halted|due to time limit|cancelled|out of memory|oom-kill'"   # Stage 3: 'exception' added for phASER (Python) tracebacks
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
grep -qE '^[^#]*Rscript ase-pipeline/tests/r/test_ase_stats_stage3\.R' "$HERE/r/run_stats_tests.sh" || { echo "FAIL: run_stats_tests.sh does not run the Stage 3 tests"; fail=1; }
# --- end Stage 3 statistics

# --- Stage 3 Task 1 review fixes (I1, M1-M8); the code lines pin the rule itself, not only its comment
need '  ru <- if (is.na(coh)) rep(NA_real_, length(smp)) else pmax(coh, ifelse(is.na(own["rho", ]), 0, own["rho", ]), rho_min)'
need '  if (anyNA(sample)) stop("hap_gene_test: sample labels must not be NA")'
need '  use <- is.finite(a) & is.finite(b) & a >= 0 & b >= 0 & n > 0 & a == round(a) & b == round(b)'
need "**One block per gene without genome-wide phasing:**"
need "it keeps only the most-covered block (or single variant) of the gene"
need "\`hap_gene_test\` does not apply \`MIN_DEPTH\`; Rmd 05 sets rows below \`MIN_DEPTH\` to \`NA\` before the call"
need "adding or removing samples changes \`rho_cohort\`"
need "**Limit: many imbalanced genes inflate the dispersion, which is conservative but costs power.**"
need "\`rho_cohort\` 0.1324, power 0.1460 against an oracle 0.5961"
need "so the haplotype-test power and tested fractions above are optimistic for that case"
forbid "so the pooled counts are not inflated by SNPs that share reads"
grep -qF 'non_integer = c(2.5, 7.5)' "$HERE/r/test_ase_stats_stage3.R" || { echo "FAIL: Stage 3 tests do not cover unusable count rows one by one"; fail=1; }
# --- end Stage 3 Task 1 review fixes

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

# --- Stage 3 Task 3 review fixes (I1-I4, M1-M3, M5, M8); structural checks run on the block cut from Step 18
# I1: no pip --user wording, and the flag allowlist is exactly the reviewed list (a new entry must be added here on purpose)
forbid "pip --user"
[ "$allow" = "bind mem mail-type mail-user array dependency parsable regions samples genotype types min-alleles max-alleles rename-chrs" ] ||
  { echo "FAIL: the flag allowlist of section (a) differs from the reviewed list"; fail=1; }
SETUP=$(bash "$HERE/synthetic/cut_block.sh" "$SKILL" '### `setup_phaser_env.sh`' 2>/dev/null)
[ -n "$SETUP" ] || { echo "FAIL: no code block under ### \`setup_phaser_env.sh\`"; fail=1; }
sline() { printf '%s\n' "$SETUP" | grep -qxF -- "$1" || { echo "FAIL: setup_phaser_env.sh lacks the line: $1"; fail=1; }; }
sline 'set -uo pipefail'
sline 'PHASER_COMMIT=aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301   # no releases or tags upstream: pinned by SHA (Step 18, bumping)'
sline 'case "$PHASER_HOME" in /*) ;; *) die "PHASER_HOME must be an absolute path, not '"'"'$PHASER_HOME'"'"'" ;; esac'
sline 'ENV_SHA=$(printf '"'"'%s\n'"'"' "$ENV_YML" | sha256sum | cut -d '"'"' '"'"' -f 1)   # line 2 of install_ok.txt: a changed package list reinstalls'
sline 'if [ "$(sed -n 1p "$PHASER_HOME/install_ok.txt" 2>/dev/null)" = "commit $PHASER_COMMIT" ] && [ "$(sed -n 2p "$PHASER_HOME/install_ok.txt" 2>/dev/null)" = "env_sha256 $ENV_SHA" ]; then'
sline '  echo "phASER already installed in $PHASER_HOME"; exit 0'
sline 'printf '"'"'%s\n'"'"' "$ENV_YML" > "$PHASER_HOME/environment.phaser.yml" || die "cannot write $PHASER_HOME/environment.phaser.yml"'
sline 'L=$(ls -A "$PHASER_HOME") || die "cannot list $PHASER_HOME"'
sline '[ -e "$PHASER_HOME/$MARK" ] || [ -z "$L" ] || die "$PHASER_HOME is not empty and has no $MARK (it was not created by setup_phaser_env.sh): choose another, new or empty directory for phASER"'
sline 'touch "$PHASER_HOME/$MARK" || die "cannot write $PHASER_HOME/$MARK"'
sline 'mkdir "$PHASER_HOME/.installing" 2>/dev/null || die "another installation is running in $PHASER_HOME (if not, remove $PHASER_HOME/.installing)"'
sline '[ "$(git -C "$PHASER_HOME/src/phaser" rev-parse HEAD)" = "$PHASER_COMMIT" ] || die "the checked-out commit differs from $PHASER_COMMIT"'
sline 'module add miniconda3/v4 || die "cannot load the miniconda3/v4 module"'
sline 'unset PYTHONPATH PYTHONHOME    # an inherited Python search path would shadow the environment the same way'
sline '"$PY" -c "import numpy, scipy, pysam, pandas, intervaltree" || die "the environment'"'"'s Python packages do not import"'
sline 'for t in python samtools bcftools bgzip tabix bedtools; do [ -x "$PHASER_HOME/env/bin/$t" ] || die "$t is missing from $PHASER_HOME/env/bin"; done'
sline '{ echo "commit $PHASER_COMMIT"; echo "env_sha256 $ENV_SHA"; "$PY" -c '"'"'import sys; print("python", sys.version.split()[0])'"'"'; date -u; } > "$PHASER_HOME/install_ok.txt" \'
printf '%s\n' "$SETUP" | grep -qE '^export PYTHONNOUSERSITE=1( |$)' || { echo "FAIL: setup_phaser_env.sh does not export PYTHONNOUSERSITE=1"; fail=1; }
[ "$(printf '%s\n' "$SETUP" | grep -c '^module add miniconda3/v4')" = 1 ] || { echo "FAIL: the miniconda3/v4 line must be inside setup_phaser_env.sh"; fail=1; }
# every pinned package of phaser_env.txt is pinned to the same version in the block (python first, then the others)
while read -r p v _; do sline "  - $p=$v"; done < <(grep -E '^(python|numpy|scipy|pandas|pysam|intervaltree|samtools|bcftools|htslib|bedtools) ' "$FIX/phaser_env.txt")
# M1: the commit named in the prose is the commit pinned in the block
pc=$(printf '%s\n' "$SETUP" | sed -n 's/^PHASER_COMMIT=\([0-9a-f]\{40\}\).*/\1/p')
[ -n "$pc" ] && grep -qF "pinned to commit \`$pc\`" "$SKILL" || { echo "FAIL: the prose commit differs from PHASER_COMMIT of setup_phaser_env.sh"; fail=1; }
# M3, M5, M8 prose
need "exports \`PYTHONNOUSERSITE=1\` and unsets \`PYTHONPATH\` and \`PYTHONHOME\`"
need "about 2.6 GB free"
need "line 1 is exactly \`commit <PHASER_COMMIT>\` and line 2 is exactly \`env_sha256 <sha256 of the package list>\`"
need "a changed package list changes this hash"
# --- end Stage 3 Task 3 review fixes

# --- Stage 3 Task 4 (Rmd 05, Step 19); each string is absent from the Task 3 skill
need "## Step 19 — Rmd 05: phASER gene-level haplotype imbalance (outbred, optional)"
need "{CWD}/{TODAY}_{WD_NAME}_05_phaser.Rmd"   # "_05_phaser.Rmd" alone is already in Step 14 (Task 1)
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
# carried items: M9 (dispersion warnings per sample), M10 (runtime scale), M11 (lists of Rmds), M8 (cohort median in the Rmd text)
need 'ht$samples$rho_warning <- ""'
need '    w <- collect_warnings(bb_estimate_rho_trim(genes$aCount[u], genes$totalCount[u]))$warn'
need "its time grows in proportion to a sample's haplotype reads: about 4.6 s per million haplotype reads (projected)"
need "is pasted verbatim into Rmd 02, 03, 04 and 05"
forbid "is pasted verbatim into Rmd 02, 03 and 04;"
need "Rmd 02, 03, 04 and 05 read its checkpoint"
forbid "Rmd 02, 03 and 04 read its checkpoint"
need "so adding or removing samples can change the p-values of the samples already analysed"
need 'if (anyDuplicated(genes[, c("sample", "gene_id")]))'
# --- end Stage 3 Task 4

# --- Stage 3 Task 4 review fixes (I1-I3, M1, M3-M6); each skill string is absent from the 84a9c8e skill
need "{TODAY}_{WD_NAME}_ASE_phaser_gene.tsv.gz"
need "{TODAY}_{WD_NAME}_ASE_phaser_comparison.tsv.gz"
need 'XL_MAX <- 1048575   # data rows per worksheet: Excel holds 1,048,576 rows including the header'
need "share the same Excel limit of 1,048,576 rows per sheet"
need "about 4.6 s per million haplotype reads"
forbid "Raise \`-t\` for more than about 200 samples of that size"
need 'stop("no gene with haplotype reads in the phASER table of sample(s) "'
need 'stop("missing counts (aCount, bCount or totalCount) in the phASER gene tables: "'
need '"; expected (columns of phaser_gene_ae at the commit pinned in Step 18): "'
need "no SNP of the gene was tested in Rmd 02 (none reached \`MIN_DEPTH\`, or none lies in its exons after the Rmd 01 filters)"
forbid '"phASER only (no SNP tested unphased)": no SNP of the gene reached `MIN_DEPTH` in Rmd 02;'
# committed proofs (I2): the evaluator proofs, the per-sample warning test, the stop tests, the orientation-invariance test
for f in prove_evaluate_stage3.R test_rmd05_warnings.R test_rmd05_stops.R check_orientation_invariance_stage3.R; do
  [ -s "$HERE/synthetic/$f" ] || { echo "FAIL: missing tests/synthetic/$f"; fail=1; }
done
# every evaluator gate tag "[name]" has a tamper proof that names it (a new gate without a proof is noticed here)
EV3="$HERE/synthetic/evaluate_stage3.R"; PR3="$HERE/synthetic/prove_evaluate_stage3.R"
tags=$(grep -oE '"\[[a-z0-9_]+\] ' "$EV3" 2>/dev/null | grep -oE '\[[a-z0-9_]+\]' | sort -u)
[ "$(printf '%s\n' "$tags" | grep -c .)" -ge 12 ] || { echo "FAIL: evaluate_stage3.R gates carry fewer than 12 [tag]s"; fail=1; }
for t in $tags; do grep -qF -- "\"$t" "$PR3" "$HERE/synthetic/prove_evaluate_stage3_phaser.R" 2>/dev/null || { echo "FAIL: prove_evaluate_stage3.R (rmd mode) or prove_evaluate_stage3_phaser.R (phaser mode) has no proof for evaluator gate $t"; fail=1; }; done
# I1: the oracle and the completeness gate read the input phASER tables, never only the Rmd's own gene table; M3: true rho parsed from the generator
grep -qF 'file.path(RES, "phaser", paste0(sm$sample, ".gene_ae.txt"))' "$EV3" 2>/dev/null || { echo "FAIL: evaluate_stage3.R does not read the input phASER tables"; fail=1; }
grep -qF 'PHI_BIO <-' "$EV3" 2>/dev/null && ! grep -qE '^ *PHI_TRUE <- [0-9]' "$EV3" || { echo "FAIL: evaluate_stage3.R must parse PHI_BIO from the generator, not hard-code it"; fail=1; }
# --- end Stage 3 Task 4 review fixes

# --- Stage 3 Task 5 (phaser_count.sh, gene spans, assembler); each string is absent from the Task 4 skill
need "### Gene spans for phASER (end of \`prep_genotypes.sh\`, before block I)"
need 'SPAN_BED="{RESULTS_DIR}/reference/genes_span.bed"'
need "### \`phaser_count.sh\`"
need '--paired_end 1 --mapq 255 --baseq {MIN_BASEQ}'
need '--pass_only 0 --unique_ids 1 --gw_phase_vcf "$PHASED_GT"'
need '--gw_phase_vcf "$PHASED_GT" --threads 1 --temp_dir "$TMPD"'   # Python 3.14 workers do not inherit phASER's global args (NameError above 1 thread)
# every --threads on a phaser.py / phaser_gene_ae.py command line (continuation lines joined) is exactly "--threads 1", and there is one
th_vals=$(awk '{ if (sub(/\\$/, "")) { buf = buf $0 " "; next } print buf $0; buf = "" }' "$SKILL" | grep -E 'phaser/phaser\.py|phaser_gene_ae/phaser_gene_ae\.py' |
  grep -oE -- '--threads[ =][^ ]+' | sed -E 's/^--threads[ =]//' | sort | uniq -c)
[ -n "$th_vals" ] && ! printf '%s\n' "$th_vals" | awk '{print $2}' | grep -qvx '1' ||
  { echo "FAIL: phASER commands must use --threads 1 (found: $(echo $th_vals)); above 1 thread phaser.py stops with NameError under Python 3.14"; fail=1; }
need "phaser.py\` with more than one thread stops with \`NameError: name 'args' is not defined\`"
need "has never been tested. **Consequence:** request one core"
forbid 'rm -f "$OUT".*'   # review I1: a glob s1.* also removes the outputs of a sample s1.redo; explicit suffixes only
need "phaser_gene_ae splits variant IDs on '_'"
need "PHASED_GT=1 but only"
need 'N_COV=$(awk -F'"'"'\t'"'"' '"'"'NR > 1 && $7 > 0'"'"' "$OUT.gene_ae.txt.part"'
need "(ASEReadCounter and phASER; a plain VCF fails there). Block R has already run before the loop; then the gene-spans block below; then block I"
[ "$(grep -c 'PHASER_COMMIT=aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301' "$SKILL")" = 3 ] || { echo "FAIL: the pinned commit must be set in setup_phaser_env.sh, phaser_count.sh and submit_chain.sh (Step 18: all three places)"; fail=1; }
[ -s "$HERE/synthetic/assemble_outbred_from_skill.sh" ] || { echo "FAIL: missing tests/synthetic/assemble_outbred_from_skill.sh"; fail=1; }
# carried items (Task 3 ruling) and guards: lines of the block cut from Step 18, so text only in prose or in setup_phaser_env.sh never satisfies them
PHC=$(bash "$HERE/synthetic/cut_block.sh" "$SKILL" '### `phaser_count.sh`' 2>/dev/null)
[ -n "$PHC" ] || { echo "FAIL: no code block under ### \`phaser_count.sh\`"; fail=1; }
pline() { printf '%s\n' "$PHC" | grep -qxF -- "$1" || { echo "FAIL: phaser_count.sh lacks the line: $1"; fail=1; }; }
pline 'set -uo pipefail'
pline 'unset PYTHONPATH PYTHONHOME    # an inherited Python search path would shadow the environment the same way'
pline 'export PYTHONNOUSERSITE=1      # a per-user site-packages directory of the same Python version would shadow the environment'
pline 'export PATH="$PHASER_HOME/env/bin:$PATH"; PY="$PHASER_HOME/env/bin/python"; SRC="$PHASER_HOME/src/phaser"'
pline '[ "$(sed -n 1p "$PHASER_HOME/install_ok.txt" 2>/dev/null)" = "commit $PHASER_COMMIT" ] && sed -n 2p "$PHASER_HOME/install_ok.txt" 2>/dev/null | grep -qE '"'"'^env_sha256 [0-9a-f]{64}$'"'"' && [ -x "$PY" ] \'
pline '[ -n "$SAMPLE" ] && [ -n "$INDIVIDUAL" ] || die "empty sample or individual in row $SLURM_ARRAY_TASK_ID of {SAMPLES_CSV}"'
pline 'clean_out                      # never keep outputs of an earlier run: a failed sample must have none'
pline '[ "$N_OFF" -eq 0 ] || die "$N_OFF heterozygous sites of $VCF lie on contigs that are not in the BAM header"'
pline '[ "$N_BLK" -gt 0 ] || { clean_out; die "phASER wrote no haplotype block (check the BAM, the VCF sample $VS and the contig names; log $LOG)"; }'
pline '  || { tail -5 "$LOG" >&2; clean_out; die "phaser.py failed (log $LOG)"; }'
pline '  || { tail -5 "$LOG" >&2; clean_out; die "phaser_gene_ae.py failed (log $LOG)"; }'
pstart() { printf '%s\n' "$PHC" | awk -v p="$1" 'index($0, p) == 1 {f = 1} END {exit !f}' || { echo "FAIL: phaser_count.sh lacks a line starting with: $1"; fail=1; }; }
pstart '[ "$NS" -eq 1 ] || die "$VCF has $NS sample columns'
pstart '[ "$N_US" -eq 0 ] || die "$N_US heterozygous sites lie on contigs whose names contain '"'_'"
pstart '[ "$N_BED" -gt 0 ] || die "no gene of $BED lies on a contig of the BAM'
pstart '  [ $((N_PH * 10)) -ge $((N_HET * 9)) ] || die "PHASED_GT=1 but only'
pstart '[ "$N_COV" -gt 0 ] || { clean_out; die "no gene has haplotype counts (totalCount > 0)'
if printf '%s\n' "$PHC" | grep -qE '^ *(module|conda|source) '; then echo "FAIL: phaser_count.sh must load no module and activate no conda environment"; fail=1; fi
SPC=$(bash "$HERE/synthetic/cut_block.sh" "$SKILL" '### Gene spans for phASER' 2>/dev/null)
printf '%s\n' "$SPC" | grep -qF '$3 == "exon" && match($9, /gene_id "[^"]+"/)' && printf '%s\n' "$SPC" | grep -qF 'print c[1] "\t" lo[g] - 1 "\t" hi[g] "\t" g' ||
  { echo "FAIL: the gene-spans block must build one 0-based row per gene from the exons"; fail=1; }
# the phaser-mode gates of the evaluator have their own tamper proofs
PR3P="$HERE/synthetic/prove_evaluate_stage3_phaser.R"
[ -s "$PR3P" ] || { echo "FAIL: missing tests/synthetic/prove_evaluate_stage3_phaser.R"; fail=1; }
for t in $(grep -oE '"\[ph_[a-z0-9_]+\] ' "$EV3" 2>/dev/null | grep -oE '\[[a-z0-9_]+\]' | sort -u); do
  grep -qF -- "\"$t" "$PR3P" 2>/dev/null || { echo "FAIL: prove_evaluate_stage3_phaser.R has no proof for evaluator gate $t"; fail=1; }; done
[ "$(grep -oE '"\[ph_[a-z0-9_]+\] ' "$EV3" 2>/dev/null | sort -u | grep -c .)" -ge 6 ] || { echo "FAIL: evaluate_stage3.R has fewer than 6 phaser-mode [ph_*] gates"; fail=1; }
# --- end Stage 3 Task 5

# --- Stage 3 Task 5 review fixes (I1, I2, M1-M7, M9-M11); each line fails on the 68e4bb8 skill
# I1: this sample's outputs by explicit suffix, never a glob
pline 'PH_SFX="haplotypic_counts.txt haplotypes.txt allelic_counts.txt allele_config.txt variant_connections.txt vcf.gz vcf.gz.tbi gene_ae.txt gene_ae.txt.part"'
pline 'clean_out() { local x; for x in $PH_SFX; do rm -f "$OUT.$x"; done; }'
if printf '%s\n' "$PHC" | grep -E 'rm ' | grep -vxF -e 'clean_out() { local x; for x in $PH_SFX; do rm -f "$OUT.$x"; done; }' -e "trap 'rm -rf \"\$TMPD\"' EXIT     # this sample's temporary files go on success and on failure" | grep -q .; then
  echo "FAIL: phaser_count.sh may remove files only through clean_out and the TMPD trap"; fail=1; fi
# M1: pins of the temp directory, the contig computations, the input guards and the --features path; the stale-output removal
#     comes before the first input guard
pline 'OUT="$R/phaser/$SAMPLE"; LOG="$R/logs/phaser_$SAMPLE.log"; TMPD="$R/tmp/phaser_$SAMPLE"'
pline "trap 'rm -rf \"\$TMPD\"' EXIT     # this sample's temporary files go on success and on failure"
pline 'N_OFF=$(bcftools query -f '"'"'%CHROM\n'"'"' "$VCF" | awk '"'"'NR == FNR {c[$1] = 1; next} !($1 in c)'"'"' "$BAM_CTG" - | wc -l)'
pline 'N_BED=$(awk '"'"'NR == FNR {c[$1] = 1; next} ($1 in c)'"'"' "$BAM_CTG" "$BED" | wc -l)'
pline '[ -s "$BAM" ] && [ -s "$BAM.bai" ] || die "missing $BAM or its index (run the per-sample array job first)"'
pline '[ -s "$VCF" ] && [ -s "$VCF.tbi" ] || die "missing $VCF or its index (run prep_genotypes.sh first)"'
pline '[ -s "$BED" ] || die "missing $BED (written by prep_genotypes.sh)"'
pline '"$PY" "$SRC/phaser_gene_ae/phaser_gene_ae.py" --haplotypic_counts "$OUT.haplotypic_counts.txt" --features "$BED" --o "$OUT.gene_ae.txt.part" >> "$LOG" 2>&1 \'
pline 'mv "$OUT.gene_ae.txt.part" "$OUT.gene_ae.txt" || { clean_out; die "cannot rename $OUT.gene_ae.txt.part"; }'
pline ': > "$LOG" || die "cannot write $LOG"   # the log describes this run only'
ln_of() { printf '%s\n' "$PHC" | awk -v p="$1" 'index($0, p) == 1 {print NR; exit}'; }
l_clean=$(ln_of 'clean_out    '); l_bam=$(ln_of '[ -s "$BAM" ]'); l_trap=$(ln_of "trap 'rm -rf")
[ -n "$l_clean" ] && [ -n "$l_bam" ] && [ -n "$l_trap" ] && [ "$l_trap" -lt "$l_clean" ] && [ "$l_clean" -lt "$l_bam" ] ||
  { echo "FAIL: phaser_count.sh must set the TMPD trap, then remove stale outputs, before the first input guard (lines trap $l_trap, clean_out $l_clean, BAM guard $l_bam)"; fail=1; }
printf '%s\n' "$SPC" | grep -qxF '    else { if ($4 < lo[g]) lo[g] = $4; if ($5 > hi[g]) hi[g] = $5 } }' ||
  { echo "FAIL: the gene-spans block must take the minimum start and maximum end over the gene's exons"; fail=1; }
# M2: no path characters in the sample name (script and Step 5)
pline 'case "$SAMPLE" in */*|*..*|.*|-*|*[[:space:]]*) die "sample name '"'"'$SAMPLE'"'"' must not contain '"'"'/'"'"', '"'"'..'"'"' or whitespace, nor start with '"'"'.'"'"' or '"'"'-'"'"' (Step 5)" ;; esac'
need "no \`/\` or \`..\` in \`sample\`, and \`sample\` does not start with \`.\`"
# M5, M6: PHASED_GT is 0 or 1; the BAM holds paired reads
pline 'case "$PHASED_GT" in 0|1) ;; *) die "PHASED_GT must be 0 or 1, not '"'"'$PHASED_GT'"'"' (Step 7)" ;; esac'
pline '[ -n "$(samtools view -f 1 "$BAM" 2>/dev/null | head -n 1)" ] || die "$BAM holds no paired reads; phASER is run with --paired_end 1 and is offered only for paired-end data (Step 7)"'
# M7: honest --mapq wording
need "phASER also drops reads whose alignment score falls below its own quantile cutoff"
forbid "the same reads as ASEReadCounter's \`--min-mapping-quality\` 10 keeps"
# I2: contig names with '_' are found before alignment: prep warns (non-fatal), the wizard refuses phASER for such a reference
printf '%s\n' "$SPC" | grep -qF 'N_US_CTG=$(cut -f1 "$FASTA.fai" | grep -c _)' ||
  { echo "FAIL: the gene-spans block of prep_genotypes.sh must count FASTA contig names with '_' (N_US_CTG)"; fail=1; }
printf '%s\n' "$SPC" | grep -qF "echo \"WARNING: \$N_US_CTG contig names of the FASTA contain '_' (for example \$(cut -f1 \"\$FASTA.fai\" | grep _ | head -3 | paste -sd' ')): phASER cannot be used" ||
  { echo "FAIL: prep_genotypes.sh must warn about contig names with '_'"; fail=1; }
need "**phASER and contig names with \`_\`.**"
need "phASER cannot be used with this reference:"
# M9: the negative tests of phaser_count.sh are committed
[ -s "$HERE/synthetic/test_phaser_count_guards.sh" ] || { echo "FAIL: missing tests/synthetic/test_phaser_count_guards.sh"; fail=1; }
for m in "lies on a contig of the BAM" "PHASED_GT=1 but only 0 of" "contigs whose names contain '_'" "no complete phASER installation" "must not contain '/'" "PHASED_GT must be 0 or 1" "holds no paired reads" "s1.redo outputs survive"; do
  grep -qF -- "$m" "$HERE/synthetic/test_phaser_count_guards.sh" 2>/dev/null || { echo "FAIL: test_phaser_count_guards.sh does not test: $m"; fail=1; }; done
# M10: the two_block oracle matches contig and position
grep -qF 'paste(act$contig, act$position)' "$EV3" 2>/dev/null || { echo "FAIL: the two_block oracle of evaluate_stage3.R must match contig and position"; fail=1; }
# M11: Step 10 points to Step 19's note on spans
need "The span includes introns: every heterozygous SNP inside it counts for the gene, including SNPs of genes nested in its introns or overlapping it (Step 19, Gene spans)."
# --- end Stage 3 Task 5 review fixes

# --- Stage 3 chain integration (wizard Steps 2/7/9, submit_chain.sh, wait_chain.sh, summary page, Notes); each check fails on ceb4667
forbid "available in a later stage"
forbid "phASER is a later stage"
forbid "Menu number 4 (phASER) runs nothing yet"
forbid "the only module ever loaded is"
forbid "the job just before the first cancelled one in chain order"
forbid "so the new hash must also replace the old one wherever the chain compares it"
need "No conda environment is needed for this skill, except for the optional phASER analysis (Steps 18-19)"
need "offered only when \`{MODE}\` = \`outbred\` and the data are paired-end"
need 'menu number 1 is the word `per-sample`, 2 is `reciprocal`, 3 is `differential`, 4 is `phaser`'
need "Are the genotype VCFs phased?"
need "Store it as \`{PHASED_GT}\`"
need "Where should phASER be installed?"
need "set \`{PHASED_GT}\` to 0 and \`{PHASER_HOME}\` to the word \`none\`"
need "\`{PHASER_RESOURCES}\` for \`phaser_count.sh\` (phASER only; one array task per sample, always one core):"
need "- otherwise: \`-n 1 --mem=16G -t 2:00:00\`."
need "these figures say nothing about real data"
need "### phASER section (only when Rmd 05 ran)"
need "\`{TODAY}_{WD_NAME}_ASE_phaser_gene.tsv.gz\` and \`{TODAY}_{WD_NAME}_ASE_phaser_comparison.tsv.gz\` (the full tables"
need "a direction (haplotype A or B higher) is given only for genome-wide phased genes"
need "run_05_phaser_<jobid>.out"
need "the one exception is \`module add miniconda3/v4\`"
need "- **phASER.** Always the repository https://github.com/secastel/phaser at the pinned commit"
# {PHASER_RESOURCES}: one core (phaser.py runs with --threads 1); no tier of that list asks for more
! grep -qE -- '`-n ([2-9]|[1-9][0-9]+) --mem=' <(sed -n '/^`{PHASER_RESOURCES}` for/,/^$/p' "$SKILL") ||
  { echo "FAIL: {PHASER_RESOURCES} must request one core (-n 1)"; fail=1; }
[ -n "$(sed -n '/^`{PHASER_RESOURCES}` for/,/^$/p' "$SKILL")" ] || { echo "FAIL: no {PHASER_RESOURCES} list in Step 9"; fail=1; }
# submit_chain.sh and wait_chain.sh: exact lines inside the cut blocks (comments or prose cannot satisfy them)
SUBC=$(bash "$HERE/synthetic/cut_block.sh" "$SKILL" "### Submission order (one block" 2>/dev/null)
WAITC=$(bash "$HERE/synthetic/cut_block.sh" "$SKILL" "**Waiting.**" 2>/dev/null)
cline() { printf '%s\n' "$SUBC" | grep -qxF -- "$1" || { echo "FAIL: submit_chain.sh lacks the line: $1"; fail=1; }; }
wline() { printf '%s\n' "$WAITC" | grep -qxF -- "$1" || { echo "FAIL: wait_chain.sh lacks the line: $1"; fail=1; }; }
cline 'MODE="{MODE}"; ANALYSES="{ANALYSES}"; PHASER_HOME="{PHASER_HOME}"   # Step 7: MODE f1 or outbred; ANALYSES the words, e.g. "per-sample phaser"; PHASER_HOME a path or none'
cline 'for w in $ANALYSES; do case "$w" in per-sample|reciprocal|differential|phaser) ;; *) die "ANALYSES contains '"'"'$w'"'"'; the accepted words are per-sample, reciprocal, differential and phaser, separated by spaces (not menu numbers or other spellings)" ;; esac; done'
cline 'case " $ANALYSES " in *" phaser "*) [ "$MODE" = outbred ] || die "phaser analysis is outbred only"; NEED="$NEED setup_phaser_env.sh phaser_count.sh run_05_phaser.sh" ;; esac'
cline 'PHASER_COMMIT=aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301; INSTALL=no'
cline '  ESHA=$(awk '"'"'index($0, "ENV_YML=$(cat <<") == 1 {f = 1; next} f && $0 == "YML" {exit} f'"'"' "$S/setup_phaser_env.sh" | sha256sum | cut -d '"'"' '"'"' -f 1)'
cline '  [ "$(sed -n 1p "$PHASER_HOME/install_ok.txt" 2>/dev/null)" = "commit $PHASER_COMMIT" ] && [ "$(sed -n 2p "$PHASER_HOME/install_ok.txt" 2>/dev/null)" = "env_sha256 $ESHA" ] || INSTALL=yes ;;'
cline 'DEP=""; R3="not selected"; R4="not selected"; R5="not selected"; PS=""   # DEP stays empty unless the mouse helper below is used'
cline 'echo "Rmd 03: $R3; Rmd 04: $R4; Rmd 05: $R5"'
pair_need '  case "$INSTALL" in yes) PS=$(sbatch -p bcc --parsable $S/setup_phaser_env.sh); got setup_phaser_env.sh "$PS" ;; esac' \
          '  PH=$(sbatch -p bcc --parsable --dependency=afterok:$A${PS:+:$PS} $S/phaser_count.sh); got phaser_count.sh "$PH"'
pair_need '  PH=$(sbatch -p bcc --parsable --dependency=afterok:$A${PS:+:$PS} $S/phaser_count.sh); got phaser_count.sh "$PH"' \
          '  R5=$(sbatch -p bcc --parsable --dependency=afterok:$PH:$R2 $S/run_05_phaser.sh); got run_05_phaser.sh "$R5" ;;'
pair_need 'case " $ANALYSES " in *" phaser "*)' '  case "$INSTALL" in yes) PS=$(sbatch -p bcc --parsable $S/setup_phaser_env.sh); got setup_phaser_env.sh "$PS" ;; esac'
# the commit of the install test in submit_chain.sh is the commit pinned by setup_phaser_env.sh ($pc, Task 3 review block)
sc=$(printf '%s\n' "$SUBC" | sed -n 's/^PHASER_COMMIT=\([0-9a-f]\{40\}\); INSTALL=no$/\1/p')
[ -n "$sc" ] && [ "$sc" = "${pc:-}" ] || { echo "FAIL: submit_chain.sh PHASER_COMMIT (${sc:-none}) differs from setup_phaser_env.sh (${pc:-none})"; fail=1; }
wline "PAT='error|exception|halted|due to time limit|cancelled|out of memory|oom-kill'   # first matching line of a job log (a time-limit kill has no \"error\" word)"
wline '                         run_05_*) echo summary_numbers_phaser.tsv ;; esac; }'
wline '                       run_05_*) echo run_02_imbalance.sh phaser_count.sh ;; esac; }'
wline '    echo "$name $id CANCELLED: never satisfiable, upstream failure: ${UP[$name]}"; bad=1'
# the stub-scheduler dry run of both scripts (bash and stubs only, under env -i; see the script header) must pass
[ -s "$HERE/chain/dry_run_chain.sh" ] || { echo "FAIL: missing tests/chain/dry_run_chain.sh"; fail=1; }
DRYD=$(mktemp -d "${TMPDIR:-/tmp}/ase_dry.XXXXXX") && {
  dry_out=$(env -i PATH=/usr/bin:/bin HOME="${HOME:-/nonexistent}" /bin/bash --noprofile --norc "$HERE/chain/dry_run_chain.sh" "$SKILL" "$DRYD/w" 2>&1)
  printf '%s\n' "$dry_out" | tail -1 | grep -qx "DRY RUN PASS" || { echo "FAIL: tests/chain/dry_run_chain.sh did not pass:"; printf '%s\n' "$dry_out" | grep -E '^FAIL|^  ' | head -20; fail=1; }
  rm -rf "$DRYD"; }
# --- end Stage 3 chain integration

# --- Stage 3 chain integration review fixes (I1, M1-M4, M6-M8, M11); each check fails on f49d308
# I1: exact phASER re-submission recipes (whole fenced block compared), and the old inexact wording is gone
forbid "\`phaser_count.sh\` without its dependency on the old array job once the BAMs exist"
need "**Re-submitting the phASER part (exact recipes).**"
need "under \`set -u\` an unset variable inside \`\$(...)\` aborts only that substitution"
rblock() { local got; got=$(bash "$HERE/synthetic/cut_block.sh" "$SKILL" "**Recipe $1 " 2>/dev/null)
  [ "$got" = "$2" ] || { echo "FAIL: recipe $1 block differs from the tested lines"; fail=1; }; }
L_PS='case "$INSTALL" in yes) PS=$(sbatch -p bcc --parsable $S/setup_phaser_env.sh); got setup_phaser_env.sh "$PS" ;; esac'
L_PH='PH=$(sbatch -p bcc --parsable ${PS:+--dependency=afterok:$PS} $S/phaser_count.sh); got phaser_count.sh "$PH"'
rblock A 'R5=$(sbatch -p bcc --parsable $S/run_05_phaser.sh); got run_05_phaser.sh "$R5"'
rblock B 'R2=$(sbatch -p bcc --parsable $S/run_02_imbalance.sh); got run_02_imbalance.sh "$R2"
R5=$(sbatch -p bcc --parsable --dependency=afterok:$R2 $S/run_05_phaser.sh); got run_05_phaser.sh "$R5"'
rblock C "$L_PS
$L_PH"'
R5=$(sbatch -p bcc --parsable --dependency=afterok:$PH $S/run_05_phaser.sh); got run_05_phaser.sh "$R5"'
rblock D 'R1=$(sbatch -p bcc --parsable $S/run_01_import_qc.sh); got run_01_import_qc.sh "$R1"
R2=$(sbatch -p bcc --parsable --dependency=afterok:$R1 $S/run_02_imbalance.sh); got run_02_imbalance.sh "$R2"
case " $ANALYSES " in *differential*)
  R4=$(sbatch -p bcc --parsable --dependency=afterok:$R1 $S/run_04_differential.sh); got run_04_differential.sh "$R4" ;;
esac'"
$L_PS
$L_PH"'
R5=$(sbatch -p bcc --parsable --dependency=afterok:$PH:$R2 $S/run_05_phaser.sh); got run_05_phaser.sh "$R5"'
for r in A B C D; do grep -q "RECIPE=$r scenario" "$HERE/chain/dry_run_chain.sh" || { echo "FAIL: the dry run does not run recipe $r"; fail=1; }; done
# M1, M2: runtime citation and the large-tier memory
forbid "phASER's authors report"
need "(\`docs/benchmarks/runtime_benchmark_report.md\` at the pinned commit) reports 89-466 s per sample"
need "The synthetic timings above say nothing for real data."
forbid "- otherwise: \`-n 1 --mem=16G -t 4:00:00\`."
# M3, M4: the Step 7 copy of the commit and the dry run's pinned hash are in the bump procedure
grep -qF "its first line is \`commit ${pc:-none}\`" "$SKILL" || { echo "FAIL: the Step 7 install test must name the commit of setup_phaser_env.sh"; fail=1; }
need "All four must agree"
[ "$(grep -c "${pc:-none}" "$SKILL")" = 5 ] || { echo "FAIL: the pinned commit must occur on exactly 5 lines (Step 7 install test, submit_chain.sh, Step 18 prose, setup_phaser_env.sh, phaser_count.sh); found $(grep -c "${pc:-none}" "$SKILL")"; fail=1; }
need "The stub dry run \`ase-pipeline/tests/chain/dry_run_chain.sh\` deliberately pins the hash of the tested installation (\`ENV_SHA_INSTALLED\`"
grep -qE '^ENV_SHA_INSTALLED=[0-9a-f]{64}$' "$HERE/chain/dry_run_chain.sh" || { echo "FAIL: dry_run_chain.sh must pin ENV_SHA_INSTALLED"; fail=1; }
# M6: commit match anchored at its end; the phASER scripts use the PHASER_HOME of submit_chain.sh
cline '  for f in setup_phaser_env.sh phaser_count.sh; do grep -q "^PHASER_\(HOME=.*; PHASER_\)\{0,1\}COMMIT=$PHASER_COMMIT\($\|[ ;]\)" "$S/$f" || die "$f does not pin phASER commit $PHASER_COMMIT (Step 18: every copy must agree) - nothing submitted"'
cline '    grep -qF "PHASER_HOME=\"$PHASER_HOME\"" "$S/$f" || die "$f uses another PHASER_HOME than '"'"'$PHASER_HOME'"'"' (substitute the Step 7 path in every script) - nothing submitted"; done'
# M7: concurrent installations into one PHASER_HOME
need "the job that starts second fails at once with \"another installation is running\" (it never waits)"
# M8: nothing-tested note of the phASER section
need "Mark a sample with \`genes_tested\` = 0 as \"nothing tested\""
# M11: the dry run deletes its scratch argument only under the shared scratch area
[ "$(grep -c 'rm -rf "$W"' "$HERE/chain/dry_run_chain.sh")" = 1 ] && grep -qF '  /net/bmc-lab3/data/bcc/ase_scratch/?*) rm -rf "$W" && mkdir -p "$W" || exit 1 ;;' "$HERE/chain/dry_run_chain.sh" ||
  { echo "FAIL: dry_run_chain.sh may delete its scratch argument only under /net/bmc-lab3/data/bcc/ase_scratch/"; fail=1; }
# --- end Stage 3 chain integration review fixes
# --- Stage 3 README and registration; each check fails on the 11162ea README and root README
RD="$HERE/../README.md"
RR="$HERE/../../README.md"
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
grep -qE "Stage 3 synthetic acceptance run: DONE \(2026-[0-9-]+\)" "$RD" && ! grep -qF "Stage 3 synthetic acceptance run: PENDING" "$RD" || { echo "FAIL: README must state Stage 3 synthetic acceptance DONE"; fail=1; }
rneed "102 evaluator checks PASS and 0 FAIL"
rneed "13 with unphased genotypes, while the unphased ACAT test of Rmd 02 detects 14"
rneed "phASER was therefore not better than ACAT on moderate genes"
grep -qE "^Stage 3 real-data smoke test: DONE \(2026-[0-9-]+\)$" "$RD" && ! grep -qF "Stage 3 real-data smoke test: PENDING" "$RD" || { echo "FAIL: README must state the Stage 3 real-data smoke test DONE (whole line, nothing trailing)"; fail=1; }
grep -q "/ase-pipeline.*phASER haplotype counts" "$RR" || { echo "FAIL: root README row must mention phASER haplotype counts"; fail=1; }
! grep -q "/ase-pipeline.*phASER not included" "$RR" || { echo "FAIL: root README row still says phASER not included"; fail=1; }
# honesty statements: each is absent from the 11162ea README
rneed "single-threaded"
rneed "NameError: name 'args' is not defined"
rneed "a Python 3.13 or older pin was never tested"
rneed "31 of the 50 genes"
rneed "power gain measured by the fragment-level simulation is optimistic"
rneed "contig names containing \`_\`"
rneed "analyses word \`phaser\`"
rneed "mode and paired-end data only"
rneed "no differential ASE on phASER counts"
rneed "median of the samples' own dispersions"
rneed "1,048,575"
rneed "ASE_phaser_gene.tsv.gz"
rneed "network access on a compute node"
rneed "Nothing on real biological data has been run"
# measured numbers (sources in the Task 7 report): each is absent from the 11162ea README
rneed "213 checks ok, 0 failed"
rneed "0.8560 against an oracle of 0.8751"
rneed "0.1460 against an oracle of 0.5961"
rneed "0.6330 / 0.6757 / 0.7789"
rneed "0.7877 against 0.5873"
rneed "hap_strong 24 of 24 gene x sample cells"
rneed "16 of 24 (unphased ACAT 12)"
rneed "15 of 18"
rneed "11 of 12"
rneed "2 of 389 (phased) and 1 of 389 (unphased), unphased ACAT 0 of 402"
rneed "66 of 66"
rneed "\`n_variants\` 5,5"
rneed "transcript offset 2"
rneed "631 s"
rneed "--threads 1"
rneed "41 checks"
# no stale phASER status text and no process labels
! grep -qiE "phASER[^.|]*(not implemented|not included|later stage)" "$RD" || { echo "FAIL: README still says phASER is not implemented, not included or a later stage"; fail=1; }
! grep -qE "phASER[^|]*not included" "$RR" || { echo "FAIL: root README still says phASER is not included"; fail=1; }
! grep -qE "(Task [0-9]|task-[0-9]|\.superpowers)" "$RD" "$RR" || { echo "FAIL: README or root README contains a process label (Task N, task-N, .superpowers)"; fail=1; }
# the 'What was not exercised' list names the Stage 3 gaps
rneed "single-end data"
rneed "phaser_pop"
# --- end Stage 3 README and registration

# --- Stage 3 real-data fix wave (gene spans, phASER resources, BLAS threads, unphased loss, README record); each check fails on 2d7b325
# (1) gene spans: only genes on contigs of the FASTA are written; the '_' warning uses the FASTA's contig names only (as Step 7)
sline() { printf '%s\n' "$SPC" | grep -qxF -- "$1" || { echo "FAIL: the gene-spans block lacks the line: $1"; fail=1; }; }
sline 'SPAN_ALL="$SPAN_BED.all.tmp"   # every gene of the GTF; only the rows on contigs of $FASTA.fai go into $SPAN_BED'
sline '  | sort -k1,1 -k2,2n > "$SPAN_ALL" || { echo "ERROR: cannot write $SPAN_ALL" >&2; exit 1; }'
sline "awk -F'\\t' 'NR == FNR {c[\$1] = 1; next} (\$1 in c)' \"\$FASTA.fai\" \"\$SPAN_ALL\" > \"\$SPAN_BED\" || { echo \"ERROR: cannot write \$SPAN_BED\" >&2; exit 1; }"
sline 'N_SPAN=$(wc -l < "$SPAN_ALL"); N_SPAN_ON=$(wc -l < "$SPAN_BED")'
sline 'N_US_CTG=$(cut -f1 "$FASTA.fai" | grep -c _)'
sline 'if [ "$N_US_CTG" -gt 0 ]; then'
[ "$(printf '%s\n' "$SPC" | grep -c '> "$SPAN_BED"')" = 1 ] || { echo "FAIL: the gene-spans block must write \$SPAN_BED only once (the FASTA-contig filter)"; fail=1; }
! printf '%s\n' "$SPC" | grep -qF 'N_US_SPAN' || { echo "FAIL: the '_' warning must not count gene-span contigs (N_US_SPAN): GTF-only patch contigs gave a false alarm"; fail=1; }
need "awk '/^>/ {print substr(\$1, 2)}' {FASTA_PATH} | grep -c _"
forbid "otherwise \`grep '^>' {FASTA_PATH} | grep -c _\`"
need "Only the FASTA's contig names count, the same input and rule as the prep job's warning (Step 10, Gene spans)"
need "Contigs with \`_\` that exist only in the GTF do not trigger it"
[ -s "$HERE/synthetic/test_gene_spans_contigs.sh" ] || { echo "FAIL: missing tests/synthetic/test_gene_spans_contigs.sh"; fail=1; }
for m in '### Gene spans for phASER' 'HSCHR6_MHC_COX' 'chrUn_test' '(for example )' 'ALL GENE-SPAN TESTS OK'; do
  grep -qF -- "$m" "$HERE/synthetic/test_gene_spans_contigs.sh" 2>/dev/null || { echo "FAIL: test_gene_spans_contigs.sh does not test: $m"; fail=1; }; done
# (2) {PHASER_RESOURCES} large tier and the measured / projected statement of Step 9
forbid "- otherwise: \`-n 1 --mem=32G -t 4:00:00\`."
need "- small genome (under 100 Mb): \`-n 1 --mem=4G -t 0:30:00\`;"
forbid "Neither the time nor the memory of a real human sample has been measured"
need "\`phaser.py\` took 2-3.6 s and about 117 MB per sample"
need "took 32 s and 1.59 GB, so memory grows mainly with the number of heterozygous sites"
need "**Projected, not measured:** a whole-genome sample like these (about 30 million read pairs) needs about 3-5 min and 2-3 GB"
need "**assumed** genome-wide heterozygous-site count"
need "it agrees with the upstream benchmark's 88 s for NA06986"
# (3) one BLAS/OpenMP thread in phaser_count.sh, on the line right after PYTHONNOUSERSITE
pline 'export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1   # nodes do not bind a -n 1 job to one core: numpy used 71 s CPU vs 6.6 s, no speed gain'
l_pnus=$(ln_of 'export PYTHONNOUSERSITE=1'); l_omp=$(ln_of 'export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1')
[ -n "$l_pnus" ] && [ -n "$l_omp" ] && [ "$l_omp" -eq $((l_pnus + 1)) ] ||
  { echo "FAIL: phaser_count.sh must export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 on the line after PYTHONNOUSERSITE (lines ${l_pnus:-none}, ${l_omp:-none})"; fail=1; }
need "the script alone used 71 s of CPU time against 6.6 s with these variables, with no gain in wall time"
# (4) Step 7: measured loss with unphased genotypes; phASER recommended mainly with phased genotypes; unphased still allowed
need "phASER is recommended mainly with phased genotypes."
need "that rule lost 76-77% of the covered SNPs of the covered genes, against 25-35% on the synthetic data"
need "You may still run phASER with unphased genotypes."
# (5) README record of the real-data smoke test (the whole-line DONE status is checked above)
rneed "1,099 of 1,140 directly connected SNP pairs (96.4%)"
rneed "96.5% / 97.1% / 98.9% agree for pairs linked by at least 2 / 5 / 10 reads"
rneed "lost 76-77% of the covered SNPs of the covered genes (1,692 of 2,218 and 1,914 of 2,478)"
rneed "The only flag change: \`--mapq 91\` was chosen for GEM's MAPQ scale, while the skill's \`--mapq 255\` is for STAR."
rneed "renamed from \`chr22\` to \`22\`"
rneed "marked them with \`samtools markdup\`"
rneed "of the 13 pairs more than 10 kb apart, 9 agree. The cause of these discordant pairs is unexplained"
rneed "uses an **assumed** genome-wide count of heterozygous sites"
rneed "no whole-genome run was made"
rneed "a median 0.74 of the phased totals"
rforbid "the real-data smoke test is pending and optional"
rforbid "nothing was measured with this pinned environment"
rforbid "-n 1 --mem=32G -t 4:00:00"
# --- end Stage 3 real-data fix wave

[ $fail -eq 0 ] && echo "PASS" || exit 1
