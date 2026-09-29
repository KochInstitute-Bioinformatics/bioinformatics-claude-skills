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
need "bb_pvalue"
need "acat"
need "ase_checkpoint.rds"
need "_ASE_imbalance.xlsx"
need "MIN_DEPTH"
need "never by globbing"
need "GenomicRanges::findOverlaps"
need "unphased, no direction"
need "_01_import_qc.Rmd"
need "_02_imbalance.Rmd"
need "cache = FALSE"
need "options(scipen = 9)"
forbid "source("
[ "$(grep -c '^# --- ase-stats-begin' "$SKILL")" -eq 1 ] || { echo "FAIL: need exactly one ase-stats-begin marker"; fail=1; }
[ "$(grep -c '^# --- ase-stats-end' "$SKILL")" -eq 1 ] || { echo "FAIL: need exactly one ase-stats-end marker"; fail=1; }
# --- end Task 6



[ $fail -eq 0 ] && echo "PASS" || exit 1
