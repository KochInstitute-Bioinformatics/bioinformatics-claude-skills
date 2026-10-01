#!/bin/bash
# Usage: bash test_gene_spans_contigs.sh <skill.md> <genome.gtf> <genome.fa.fai> <work dir (must not exist)> [expected BED]
# Tests the gene-span block of prep_genotypes.sh (Step 10, cut from the skill) and the wizard's contig check of Step 7 (the two
# commands cut from the "phASER and contig names with `_`" paragraph) on the synthetic genome, with three variants:
#   A. the genome as it is: every gene is written, no warning, Step 7 counts 0 (and, with [expected BED], the BED is unchanged)
#   B. one extra gene on a patch contig "HSCHR6_MHC_COX" that exists only in the GTF: no warning, the contig is not in the BED,
#      the log counts it among the genes of the GTF only; Step 7 counts 0 (phASER offered)
#   C. one extra contig "chrUn_test" in the FASTA (.fai and header) with a gene on it: the warning fires and names chrUn_test;
#      Step 7 counts 1 with the .fai and with the FASTA header lines (phASER refused)
# The FASTA used for the Step 7 header command holds header lines only (one 'N' per contig), with a description that contains '_'
# so the first-word rule is tested. Text tools only (awk, sort, grep, bash): it may run on the login node on these small files.
# Ends with "ALL GENE-SPAN TESTS OK" (exit 0) or "SOME GENE-SPAN TESTS BAD" (exit 1).
set -uo pipefail
SKILL=${1:?usage}; GTF=${2:?usage}; FAI=${3:?usage}; WORK=${4:?usage}; EXP_BED=${5:-}; HERE=$(cd "$(dirname "$0")" && pwd)
die() { echo "ERROR: $*" >&2; exit 1; }
case "$WORK" in /*) ;; *) die "work dir must be an absolute path" ;; esac
[ ! -e "$WORK" ] || die "$WORK exists; give a new work dir"
[ -s "$SKILL" ] && [ -s "$GTF" ] && [ -s "$FAI" ] || die "skill, GTF or .fai missing or empty"
mkdir -p "$WORK" || die "cannot create $WORK"
BLOCK=$(bash "$HERE/cut_block.sh" "$SKILL" "### Gene spans for phASER") || die "no gene-span block in $SKILL"
PARA=$(grep -F '**phASER and contig names with `_`.**' "$SKILL") || die "no Step 7 contig paragraph in $SKILL"
C_FAI=$(printf '%s\n' "$PARA" | grep -oE '`cut -f1 \{FASTA_PATH\}\.fai[^`]*`' | head -1 | tr -d '`')
C_HDR=$(printf '%s\n' "$PARA" | grep -oE "\`awk '/\^>/[^\`]*\{FASTA_PATH\}[^\`]*\`" | head -1 | tr -d '`')
[ -n "$C_FAI" ] && [ -n "$C_HDR" ] || die "cannot cut the two Step 7 commands (got: '$C_FAI' / '$C_HDR')"
N_GTF_GENES=$(awk -F'\t' '$3 == "exon"' "$GTF" | grep -o 'gene_id "[^"]*"' | sort -u | wc -l)
bad=0
ok()  { echo "ok   $*"; }
nok() { echo "BAD  $*"; bad=1; }
# run_case <name> <gtf> <fai>: builds the case dir, runs the block, prints nothing; sets LOGF, BEDF, FA
run_case() {
  local d="$WORK/$1"; mkdir -p "$d/genome" || die "mkdir $d"
  cp "$2" "$d/genome/genome.gtf" && cp "$3" "$d/genome/genome.fa.fai" || die "copy for case $1"
  FA="$d/genome/genome.fa"
  awk -F'\t' '{print ">" $1 " dna:chromosome note_with_underscore"; print "N"}' "$d/genome/genome.fa.fai" > "$FA" || die "fasta $1"
  printf '%s\n' "$BLOCK" | awk -v r="$d/results" -v g="$d/genome/genome.gtf" '{gsub(/\{RESULTS_DIR\}/, r); gsub(/\{GTF_PATH\}/, g); print}' > "$d/span_block.sh"
  LOGF="$d/span.log"; BEDF="$d/results/reference/genes_span.bed"
  FASTA="$FA" bash "$d/span_block.sh" > "$LOGF" 2>&1; echo "rc $?" >> "$LOGF"
}
step7() { local c=${1//\{FASTA_PATH\}/$FA}; bash -c "$c"; }
N_FAI=$(wc -l < "$FAI")

# A. the genome as it is
run_case A "$GTF" "$FAI"
grep -qxF "rc 0" "$LOGF" && ok "A: block exit 0" || nok "A: block failed: $(tail -3 "$LOGF")"
grep -qxF "Gene spans for phASER: $N_GTF_GENES genes, $N_GTF_GENES on contigs of the FASTA (only these are written to $BEDF)" "$LOGF" \
  && ok "A: $N_GTF_GENES/$N_GTF_GENES gene spans" || nok "A: span line: $(grep 'Gene spans' "$LOGF")"
grep -q WARNING "$LOGF" && nok "A: unexpected warning" || ok "A: no warning"
[ "$(wc -l < "$BEDF")" = "$N_GTF_GENES" ] && ok "A: BED has $N_GTF_GENES rows" || nok "A: BED rows $(wc -l < "$BEDF")"
[ ! -e "$BEDF.all.tmp" ] && ok "A: temporary span file removed" || nok "A: $BEDF.all.tmp left"
if [ -n "$EXP_BED" ]; then cmp -s "$BEDF" "$EXP_BED" && ok "A: BED identical to $EXP_BED" || nok "A: BED differs from $EXP_BED"; fi
[ "$(step7 "$C_FAI")" = 0 ] && [ "$(step7 "$C_HDR")" = 0 ] && ok "A: Step 7 counts 0 (.fai and header; '_' in descriptions ignored)" \
  || nok "A: Step 7 counts $(step7 "$C_FAI") / $(step7 "$C_HDR")"

# B. a patch contig with '_' only in the GTF
GB="$WORK/gtf_patch.gtf"
{ cat "$GTF"; printf 'HSCHR6_MHC_COX\tsynthetic\texon\t100\t900\t.\t+\t.\tgene_id "PATCH0001"; transcript_id "PATCH0001.1";\n'; } > "$GB"
run_case B "$GB" "$FAI"
grep -qxF "rc 0" "$LOGF" && ok "B: block exit 0" || nok "B: block failed: $(tail -3 "$LOGF")"
grep -qxF "Gene spans for phASER: $((N_GTF_GENES + 1)) genes, $N_GTF_GENES on contigs of the FASTA (only these are written to $BEDF)" "$LOGF" \
  && ok "B: $((N_GTF_GENES + 1)) genes in the GTF, $N_GTF_GENES written" || nok "B: span line: $(grep 'Gene spans' "$LOGF")"
grep -q WARNING "$LOGF" && nok "B: warning fired for a GTF-only '_' contig: $(grep WARNING "$LOGF" | cut -c1-160)" || ok "B: no warning"
cut -f1 "$BEDF" | grep -q _ && nok "B: a '_' contig is in the BED" || ok "B: no '_' contig in the BED"
grep -qF PATCH0001 "$BEDF" && nok "B: the patch gene is in the BED" || ok "B: patch gene not in the BED"
[ "$(step7 "$C_FAI")" = 0 ] && [ "$(step7 "$C_HDR")" = 0 ] && ok "B: Step 7 counts 0 (phASER offered)" \
  || nok "B: Step 7 counts $(step7 "$C_FAI") / $(step7 "$C_HDR")"

# C. a contig with '_' in the FASTA (.fai) with a gene on it
GC="$WORK/gtf_un.gtf"; FC="$WORK/fai_un.fai"
{ cat "$GTF"; printf 'chrUn_test\tsynthetic\texon\t100\t900\t.\t+\t.\tgene_id "UN0001"; transcript_id "UN0001.1";\n'; } > "$GC"
{ cat "$FAI"; printf 'chrUn_test\t5000\t0\t60\t61\n'; } > "$FC"
run_case C "$GC" "$FC"
grep -qxF "rc 0" "$LOGF" && ok "C: block exit 0 (warning only)" || nok "C: block failed: $(tail -3 "$LOGF")"
grep -qF "WARNING: 1 contig names of the FASTA contain '_' (for example chrUn_test): phASER cannot be used with this reference" "$LOGF" \
  && ok "C: warning fires and names chrUn_test" || nok "C: warning: $(grep WARNING "$LOGF" | cut -c1-200)"
grep -q "(for example )" "$LOGF" && nok "C: empty example list" || ok "C: example list not empty"
[ "$(step7 "$C_FAI")" = 1 ] && [ "$(step7 "$C_HDR")" = 1 ] && ok "C: Step 7 counts 1 with the .fai and with the header lines (phASER refused)" \
  || nok "C: Step 7 counts $(step7 "$C_FAI") / $(step7 "$C_HDR")"

[ "$bad" = 0 ] && { echo "ALL GENE-SPAN TESTS OK"; exit 0; }
echo "SOME GENE-SPAN TESTS BAD"; exit 1
