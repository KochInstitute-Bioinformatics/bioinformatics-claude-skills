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
# bcftools (regions samples genotype types min-alleles max-alleles).
allow="bind mem mail-type mail-user array dependency parsable regions samples genotype types min-alleles max-alleles"
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
need "star_ase_masked_sjdb{SJDB_OVERHANG}"
need "genomeSAindexNbases"
need "max(4,"
need "{SAMPLES_CSV}"
need "sample,fastq_1,fastq_2,condition,cross_direction,individual"
need "MIN_DEPTH"
need "FDR_SIG"
need "ABS_DEV_SIG"
need "BIAS_TOL"
need "offered only when \`{MODE}\` = \`f1\` and both \`cross_direction\` values"
need "offered only when at least two conditions each have replicates"
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
grep -qF "PENDING" "$HERE/../README.md" 2>/dev/null || { echo "FAIL: ase-pipeline/README.md must mark the end-to-end acceptance run as PENDING"; fail=1; }
grep -q "/ase-pipeline" "$HERE/../../README.md" 2>/dev/null || { echo "FAIL: root README.md has no /ase-pipeline row"; fail=1; }
grep -qF "Nothing on real biological data" "$HERE/../README.md" 2>/dev/null || { echo "FAIL: ase-pipeline/README.md must state: Nothing on real biological data"; fail=1; }
grep -qF "never submitted" "$HERE/../README.md" 2>/dev/null || { echo "FAIL: ase-pipeline/README.md must state which scripts were never submitted"; fail=1; }
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


[ $fail -eq 0 ] && echo "PASS" || exit 1
