#!/bin/bash
# Usage: test_strandedness.sh <skill.md>
# Runs infer_strandedness from the skill (block after "**Strandedness from an existing nf-core/rnaseq run.**") on RSeQC
# infer_experiment.txt files written here. Prints "STRANDEDNESS PASS" or exits 1.
set -u
source "$(cd "$(dirname "$0")" && pwd)/test_env.sh" "$@"
SKILL=${1:?usage: test_strandedness.sh <skill.md>}
HERE=$(cd "$(dirname "$0")" && pwd)
CODE=$(bash "$HERE/cut_block.sh" "$SKILL" "**Strandedness from an existing nf-core/rnaseq run.**") || { echo "FAIL: strandedness block not found"; exit 1; }
eval "$CODE"
declare -F infer_strandedness >/dev/null || { echo "FAIL: infer_strandedness is not defined by the block"; exit 1; }
mk() { # mk <name> <PE|SE> <fraction explained by sense pattern> <fraction explained by antisense pattern>
  if [ "$2" = PE ]; then
    printf '\n\nThis is PairEnd Data\nFraction of reads failed to determine: 0.0100\nFraction of reads explained by "1++,1--,2+-,2-+": %s\nFraction of reads explained by "1+-,1-+,2++,2--": %s\n' "$3" "$4" > "$TEST_TMP/$1"
  else
    printf '\n\nThis is SingleEnd Data\nFraction of reads failed to determine: 0.0100\nFraction of reads explained by "++,--": %s\nFraction of reads explained by "+-,-+": %s\n' "$3" "$4" > "$TEST_TMP/$1"
  fi; }
fail=0
check() { local got; got=$(infer_strandedness "$TEST_TMP/$1"); [ "$got" = "$2" ] || { echo "FAIL: case $1: expected $2, got '$got'"; fail=1; }; }
mk pe_forward PE 0.9512 0.0388;    check pe_forward forward
mk pe_reverse PE 0.0210 0.9690;    check pe_reverse reverse
mk pe_unstranded PE 0.4950 0.4950; check pe_unstranded unstranded
mk se_reverse SE 0.0500 0.9400;    check se_reverse reverse
mk se_forward SE 0.8800 0.1100;    check se_forward forward
mk pe_ambiguous PE 0.7000 0.2900;  check pe_ambiguous unclear
: > "$TEST_TMP/empty";             check empty unclear
printf 'Fraction of reads explained by "1++,1--,2+-,2-+": 0.95\n' > "$TEST_TMP/one_line"; check one_line unclear
# nf-core/rnaseq rule: forward if f >= 0.8, reverse if r >= 0.8, unstranded if |f - r| <= 0.1, else unclear (boundaries)
mk pe_f_at_080 PE 0.8000 0.1900;   check pe_f_at_080 forward
mk pe_f_below PE 0.7999 0.1900;    check pe_f_below unclear
mk se_r_at_080 SE 0.1900 0.8000;   check se_r_at_080 reverse
mk pe_diff_010 PE 0.4500 0.5500;   check pe_diff_010 unstranded
mk pe_diff_011 PE 0.4400 0.5500;   check pe_diff_011 unclear
mk pe_70_30 PE 0.7000 0.3000;      check pe_70_30 unclear
# a missing file, and a file with repeated lines (two reports concatenated), are unclear
check no_such_file unclear
mk rep_a PE 0.0200 0.9700; mk rep_b PE 0.9700 0.0200; cat "$TEST_TMP/rep_a" "$TEST_TMP/rep_b" > "$TEST_TMP/repeated"; check repeated unclear
forbidden_ran && fail=1
[ $fail -eq 0 ] && echo "STRANDEDNESS PASS" || exit 1
