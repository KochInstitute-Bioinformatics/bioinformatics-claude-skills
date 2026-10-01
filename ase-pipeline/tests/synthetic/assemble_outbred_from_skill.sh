#!/bin/bash
# Usage: assemble_outbred_from_skill.sh <skill.md> <values.tsv> <genotype_map.tsv>
# Builds every outbred-mode script and Rmd 01/02/05 of a project from the skill's code blocks (cut by headings or anchor lines,
# never by line numbers), as the wizard would: placeholders substituted from values.tsv (NAME<TAB>value; must define CWD,
# RESULTS_DIR, TODAY and WD_NAME); the genotype map (individual<TAB>vcf<TAB>sample in vcf) written into prep_genotypes.sh; the test-data
# settings FILTER_EXPR='' and MIN_SITES_WARN=20; the GTF bind dropped from the run scripts (the test genome lies inside CWD); the
# four mouse-helper lines of submit_chain.sh deleted. Stops on any unsubstituted placeholder.
# Every block is cut inside a pipeline, so each pipeline's status is checked (a missing block stops the assembly instead of
# leaving a truncated script); text only, no tool is run, so it may run on the login node.
set -uo pipefail
SKILL=${1:?}; VALS=${2:?}; MAP=${3:?}; HERE=$(cd "$(dirname "$0")" && pwd)
die() { echo "ERROR: $*" >&2; exit 1; }
[ -s "$SKILL" ] && [ -s "$VALS" ] && [ -s "$MAP" ] || die "skill, values or genotype map missing or empty"
val() { awk -F'\t' -v k="$1" '$1 == k {print $2; f = 1} END {exit !f}' "$VALS" || die "values file has no $1"; }
CWD=$(val CWD) || exit 1; RD=$(val RESULTS_DIR) || exit 1; TD=$(val TODAY) || exit 1; WN=$(val WD_NAME) || exit 1
TAG="${TD}_$WN"; S=$RD/scripts
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
  blk "### Shared block I"; echo 'echo "prep_genotypes done"'; } | sub > "$S/prep_genotypes.sh" || die "cannot assemble prep_genotypes.sh"
# align_wasp_count.sh: header, block C with the check-only fetch_sif and its four calls, common start, WASP block, empty-table guard
{ hdr 'align_wasp_count_%A_%a.out' '{ARRAY_RESOURCES}' array
  blk "### Shared block C"
  echo 'fetch_sif() { [ -s "$1" ] || { echo "ERROR: missing container $1 (run the prep job first)" >&2; exit 1; }; }'
  echo 'fetch_sif "$STAR_SIF"; fetch_sif "$GATK_SIF"; fetch_sif "$SAMTOOLS_SIF"; fetch_sif "$PICARD_SIF"'
  blk "**Common start of both scripts**"; echo '[ -n "$INDIVIDUAL" ] || die "empty individual"'
  echo 'HET_PLAIN="$R/genotypes/$INDIVIDUAL.het.vcf"; HET_GZ="$R/genotypes/$INDIVIDUAL.het.vcf.gz"'
  blk '### `align_wasp_count.sh` (outbred mode)'; blk "### Empty-table guard"; } | sub > "$S/align_wasp_count.sh" || die "cannot assemble align_wasp_count.sh"
# Rmd run scripts from the run_01 template (first bind case: the test genome lies inside CWD)
R01=$(blk "The template below shows the second case:") || exit 1
for r in 01_import_qc 02_imbalance 05_phaser; do
  printf '%s\n' "$R01" | sed -e "s/ase_01_import_qc/ase_$r/; s/run_01_import_qc_/run_${r}_/; s/_01_import_qc\.Rmd/_$r.Rmd/; s/Rmd 01 failed/Rmd ${r%%_*} failed/" \
    -e 's/--bind {CWD},{GTF_DIR} /--bind {CWD} /' | sub > "$S/run_$r.sh" || die "cannot assemble run_$r.sh"
done
sed -i 's/-t 1:00:00/-t 4:00:00/' "$S/run_05_phaser.sh" || die "cannot edit run_05_phaser.sh"     # Step 19: -n 1 --mem=16G -t 4:00:00
blk '### `setup_phaser_env.sh`' | sub > "$S/setup_phaser_env.sh" || die "cannot assemble setup_phaser_env.sh"
blk '### `phaser_count.sh`' | sub > "$S/phaser_count.sh" || die "cannot assemble phaser_count.sh"
blk "### Submission order (one block" | awk '/^# Without the helper, delete the next four lines/ {print; skip = 4; next} skip > 0 {skip--; next} {print}' | sub > "$S/submit_chain.sh" \
  || die "cannot assemble submit_chain.sh"
blk "**Waiting.**" | sub > "$S/wait_chain.sh" || die "cannot assemble wait_chain.sh"
for n in 12:01_import_qc 13:02_imbalance 19:05_phaser; do
  bash "$HERE/render_from_skill.sh" "$SKILL" "${n%%:*}" "$VALS" "$CWD/${TAG}_${n#*:}.Rmd" > /dev/null || die "render of Step ${n%%:*} failed"
done
for f in "$S"/*.sh; do [ -s "$f" ] || die "empty script $f"; done
if grep -nE '\{[A-Z][A-Z_0-9]*\}' "$S"/*.sh; then die "unsubstituted placeholders above"; fi
for f in "$S"/*.sh; do bash -n "$f" || die "syntax error in $f"; done
echo "assembled $(ls "$S" | wc -l) scripts in $S and 3 Rmds in $CWD"
