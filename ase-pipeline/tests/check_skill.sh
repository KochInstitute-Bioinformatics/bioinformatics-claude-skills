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
need "F1 mode in Stage 1 needs one parent to be the reference strain"
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
# I6: SNP-dependence limitation documented (Step 13, Step 14, Rmd 02 gene table, README)
need "**Independence assumption (F1 gene LRT and \`bb_estimate_rho_gene\`).**"
need "come from simulations with **independent SNPs**"
need 'note = "SNP counts treated as independent: may be anti-conservative in SNP-dense genes at moderate depth"'
need "**Limitation (F1):** SNP counts that share read pairs are treated as independent"
grep -qF "SNPs that share read pairs are treated as independent" "$RD" || { echo "FAIL: ase-pipeline/README.md must state the SNP-dependence limitation of the F1 gene test"; fail=1; }
grep -qF "Future work (Stage 2): thin each gene's SNPs" "$RD" || { echo "FAIL: ase-pipeline/README.md must list SNP thinning as Stage 2 future work"; fail=1; }
# minors: GTF bind, fetch_sif calls, wrapper sets, paths with spaces, dropping a failed sample, README wording, spec
need 'add_bind "$(dirname "{FASTA_PATH}")"; add_bind "$(dirname "{GTF_PATH}")"'
need "must **call** it once per container it uses"
need 'call** it for each container the script uses, on the line after the definition: `fetch_sif "$STAR_SIF"; fetch_sif "$GATK_SIF"; fetch_sif "$SAMTOOLS_SIF"; fetch_sif "$PICARD_SIF"`'
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


# --- Stage 2 Task 1 (statistics); each string is absent from the Stage 1 skill (152dc20)
need "chrom_class <- function(contig)"
need "thin_snps <- function(pos, depth, window)"
need "bb_glm_fit <- function(y, n, X, rho, bound = 15)"
need "bb_glm_lrt <- function(y, n, X, drop_cols, rho, bound = 15)"
need "bb_moment_phi <- function(units, max_iter = 25, tol = 1e-6, bound = 15)"
need "ase_glm_test <- function(units, tests, rho_min, min_df_unit = 2, bound = 15)"
need "ase_paired_test <- function(d, rho_min, min_individuals = 2, bound = 15)"
need "bb_pair_phi_trim <- function(y, n, cell, keep = 0.9, max_iter = 50)"
need '| `THIN_BP` | 500 |'
need "pasted verbatim into Rmd 02, Rmd 03 and Rmd 04"
need "**Stage 2 models (Rmd 03 and 04): \`ase_glm_test\`.**"
need "**Outbred differential: \`ase_paired_test\`.**"
need "**SNPs sharing read pairs.**"
[ -s "$HERE/r/test_ase_stats_stage2.R" ] || { echo "FAIL: missing tests/r/test_ase_stats_stage2.R"; fail=1; }
# --- end Stage 2 Task 1

[ $fail -eq 0 ] && echo "PASS" || exit 1
