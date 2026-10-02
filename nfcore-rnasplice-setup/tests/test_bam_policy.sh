#!/bin/bash
# Usage: test_bam_policy.sh <skill.md>
# Runs bam_input_allowed from the skill (block after "**BAM input rule.**") for every strandedness x read type and compares
# each answer with the expectation computed here from fixtures/gate_values.tsv. Prints "BAM POLICY PASS" or exits 1.
set -u
source "$(cd "$(dirname "$0")" && pwd)/test_env.sh" "$@"
SKILL=${1:?usage: test_bam_policy.sh <skill.md>}
HERE=$(cd "$(dirname "$0")" && pwd)
gv() { awk -F'\t' -v k="$1" '$1 == k {print $2; exit}' "$HERE/fixtures/gate_values.tsv"; }
CODE=$(bash "$HERE/cut_block.sh" "$SKILL" "**BAM input rule.**") || { echo "FAIL: BAM input rule block not found"; exit 1; }
eval "$CODE"
declare -F bam_input_allowed >/dev/null || { echo "FAIL: bam_input_allowed is not defined by the block"; exit 1; }
G_LT=$(gv BAM_RMATS_LIBTYPE); G_RT=$(gv BAM_RMATS_READTYPE); G_DS=$(gv BAM_DEXSEQ_STRAND); G_FC=$(gv BAM_FC_STRAND)
fail=0
for lib in unstranded forward reverse; do
  case $lib in unstranded) lt=fr-unstranded ds=no fc=0 ;; forward) lt=fr-secondstrand ds=yes fc=1 ;; reverse) lt=fr-firststrand ds=reverse fc=2 ;; esac
  for lay in paired single; do
    if [ "$G_LT" = none ]; then exp=REFUSED
    elif [ "$G_LT" = from_sheet ]; then exp=ALLOWED
    elif [ "$lt" = "$G_LT" ] && [ "$lay" = "$G_RT" ] && { [ "$G_DS" = n/a ] || [ "$ds" = "$G_DS" ]; } && { [ "$G_FC" = n/a ] || [ "$fc" = "$G_FC" ]; }; then exp=ALLOWED
    else exp=REFUSED; fi
    out=$(bam_input_allowed "$lib" "$lay"); rc=$?
    if [ "$exp" = ALLOWED ]; then
      [ $rc -eq 0 ] && [[ $out == ALLOWED* ]] || { echo "FAIL: case $lib/$lay: expected ALLOWED, got rc=$rc '$out'"; fail=1; }
    else
      [ $rc -ne 0 ] && [[ $out == REFUSED* ]] && [[ $out == *"Start from FASTQ"* ]] || { echo "FAIL: case $lib/$lay: expected a REFUSED answer pointing to FASTQ, got rc=$rc '$out'"; fail=1; }
    fi
  done
done
out=$(bam_input_allowed auto paired); rc=$?; [ $rc -ne 0 ] && [[ $out == REFUSED* ]] || { echo "FAIL: case invalid_strandedness"; fail=1; }
out=$(bam_input_allowed unstranded both); rc=$?; [ $rc -ne 0 ] && [[ $out == REFUSED* ]] || { echo "FAIL: case invalid_layout"; fail=1; }
forbidden_ran && fail=1
[ $fail -eq 0 ] && echo "BAM POLICY PASS" || exit 1
