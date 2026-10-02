#!/bin/bash
# Usage: test_sample_names.sh <skill.md>
# Runs the functions of the skill's block after "**Sample names and read pairs.**" (fastq_r2_name, fastq_sample_name,
# sanitize_sample_name, sample_name_collisions, input_path_ok) on fixed cases, and checks every sanitised name against the
# sample rule of the pinned samplesheet schemas (fixtures/schema_input.json and schema_input_genome_bam.json, copied from
# nf-core/rnasplice 1b447239488097651d8eac44bca2c1556865eb0f): the pattern and the reserved-word list are read from the
# fixtures, so a schema change makes this test fail. Prints "SAMPLE NAMES PASS" or exits 1.
set -u
source "$(cd "$(dirname "$0")" && pwd)/test_env.sh" "$@"
SKILL=${1:?usage: test_sample_names.sh <skill.md>}
HERE=$(cd "$(dirname "$0")" && pwd); FIX="$HERE/fixtures"
fail=0
# The schema rule this test translates: pattern (as a POSIX ERE below) + "not enum" list. Both schemas must carry it.
PAT_JSON='^(?!\\.\\.\\d+)(?!\\.$)[a-zA-Z.]([a-zA-Z0-9._]*)?$'
for s in schema_input.json schema_input_genome_bam.json; do
  [ -s "$FIX/$s" ] || { echo "FAIL: fixture $s missing"; exit 1; }
  got=$(awk '/"sample": \{/ {f = 1} f && /"pattern"/ {sub(/^[ \t]*"pattern": "/, ""); sub(/",[ \t]*$/, ""); print; exit}' "$FIX/$s")
  [ "$got" = "$PAT_JSON" ] || { echo "FAIL: $s sample pattern changed: '$got' (update the skill rule and this test)"; fail=1; }
done
# A name that starts with a letter needs no lookahead: ^[a-zA-Z.] then [a-zA-Z0-9._]*.
PAT_ERE='^[a-zA-Z.][a-zA-Z0-9._]*$'
RESERVED=$(awk '/"sample": \{/ {f = 1} f && /"enum": \[/ {e = 1; next} e && /\]/ {exit} e {gsub(/[ \t",]/, ""); print}' "$FIX/schema_input.json")
[ "$(echo "$RESERVED" | grep -c .)" -ge 20 ] || { echo "FAIL: reserved-word list not read from schema_input.json"; exit 1; }
schema_ok() { [[ $1 =~ $PAT_ERE ]] && ! echo "$RESERVED" | grep -qxF -- "$1"; }

CODE=$(bash "$HERE/cut_block.sh" "$SKILL" "**Sample names and read pairs.**") || { echo "FAIL: sample-name block not found"; exit 1; }
eval "$CODE"
for fn in fastq_r2_name fastq_sample_name sanitize_sample_name sample_name_collisions input_path_ok; do
  declare -F "$fn" >/dev/null || { echo "FAIL: $fn is not defined by the block"; exit 1; }
done

# 1. R1 -> R2 file name (paired-end detection and pairing)
r2() { local got rc; got=$(fastq_r2_name "$1"); rc=$?
  if [ -n "$2" ]; then [ $rc -eq 0 ] && [ "$got" = "$2" ] || { echo "FAIL: r2 case $1: expected '$2', got rc=$rc '$got'"; fail=1; }
  else [ $rc -ne 0 ] && [ -z "$got" ] || { echo "FAIL: r2 case $1: expected no R2 (not an R1 file), got rc=$rc '$got'"; fail=1; }; fi; }
r2 X_S1_L001_R1_001.fastq.gz X_S1_L001_R2_001.fastq.gz
r2 ctrl1_R1.fastq.gz ctrl1_R2.fastq.gz
r2 data/ctrl1_R1.fq.gz data/ctrl1_R2.fq.gz
r2 s_1.fastq.gz s_2.fastq.gz
r2 s_1.fq.gz s_2.fq.gz
r2 s_1_sequence.fastq.gz s_2_sequence.fastq.gz
r2 A_R1_rep_R1_001.fastq.gz A_R1_rep_R2_001.fastq.gz
r2 ctrl1_R2.fastq.gz ""
r2 ctrl.fastq.gz ""
r2 R1sample_2.fastq.gz ""

# 2. Sample name cut (before sanitisation)
nm() { local got; got=$(fastq_sample_name "$1"); [ "$got" = "$2" ] || { echo "FAIL: name case $1: expected '$2', got '$got'"; fail=1; }; }
nm X_S1_L001_R1_001.fastq.gz X
nm X_S1_L002_R2_001.fastq.gz X
nm ctrl1_R1.fastq.gz ctrl1
nm data/ctrl1_R2.fq.gz ctrl1
nm s_1.fq.gz s
nm s_2.fastq.gz s
nm s_1_sequence.fastq.gz s
nm ctrl.fastq.gz ctrl
nm ctrl.fq.gz ctrl
nm data/WT-1_R1.fastq.gz WT-1
nm KO_S12.fastq.gz KO

# 3. Sanitisation: expected value, and every result satisfies the pinned schema rule
sn() { local got; got=$(sanitize_sample_name "$1")
  [ "$got" = "$2" ] || { echo "FAIL: sanitise case '$1': expected '$2', got '$got'"; fail=1; }
  schema_ok "$got" || { echo "FAIL: sanitise case '$1': '$got' violates the pinned schema rule"; fail=1; }; }
sn 1abc S1abc
sn _x S_x
sn -x S_x
sn NA NA_S
sn TRUE TRUE_S
sn in in_S
sn NaN NaN_S
sn NA_integer_ NA_integer__S
sn function function_S
sn WT-1 WT_1
sn "WT 1(a)" WT_1_a_
sn a.b a_b
sn ctrl1 ctrl1
sn Na Na
sn If If
sn "" S
# every reserved word of the schema list that can survive the character rule
while IFS= read -r w; do
  case "$w" in *.*) continue ;; esac
  got=$(sanitize_sample_name "$w"); schema_ok "$got" || { echo "FAIL: reserved word '$w' -> '$got' violates the pinned schema rule"; fail=1; }
done <<< "$RESERVED"

# 4. Collisions after sanitisation: lanes (same name before sanitisation) vs collisions created by sanitisation
out=$(printf '%s\n' X X WT-1 WT_1 ctrl 1a S1a | sample_name_collisions)
echo "$out" | grep -qxF "X: same name before sanitisation (lanes or technical replicates?)" || { echo "FAIL: collision case lanes: got '$out'"; fail=1; }
echo "$out" | grep -qxF "WT_1: created by sanitisation from WT-1, WT_1 (different samples?)" || { echo "FAIL: collision case sanitisation: got '$out'"; fail=1; }
echo "$out" | grep -qxF "S1a: created by sanitisation from 1a, S1a (different samples?)" || { echo "FAIL: collision case digit prefix: got '$out'"; fail=1; }
echo "$out" | grep -q "^ctrl:" && { echo "FAIL: collision case unique: ctrl reported: '$out'"; fail=1; }
[ "$(echo "$out" | grep -c .)" -eq 3 ] || { echo "FAIL: collision case count: expected 3 lines, got '$out'"; fail=1; }
out=$(printf '%s\n' a b c | sample_name_collisions); [ -z "$out" ] || { echo "FAIL: collision case none: got '$out'"; fail=1; }

# 5. Paths: whitespace or a comma stops (the schemas reject whitespace; a comma breaks the CSV)
po() { local got rc; got=$(input_path_ok "$1"); rc=$?
  if [ "$2" = ok ]; then [ $rc -eq 0 ] || { echo "FAIL: path case '$1': expected ok, got rc=$rc '$got'"; fail=1; }
  else [ $rc -ne 0 ] && [[ $got == STOP* ]] || { echo "FAIL: path case '$1': expected STOP, got rc=$rc '$got'"; fail=1; }; fi; }
po data/ok_R1.fastq.gz ok
po /abs/dir/s.bam ok
po "data/a b_R1.fastq.gz" stop
po "my data/s.bam" stop
po data/a,b_R1.fastq.gz stop
po "$(printf 'data/a\tb.bam')" stop

forbidden_ran && fail=1
[ $fail -eq 0 ] && echo "SAMPLE NAMES PASS" || exit 1
