#!/bin/bash
# Usage: test_read_length.sh <skill.md>
# Runs the read-length block of the skill (block after "**Read-length detection.**": detect_read_length, read_length_report)
# and the run lines (block after "**Running the read-length check.**", placeholders filled in) on generated gzipped FASTQ.
# Prints "READ LENGTH PASS" or exits 1.
set -u
source "$(cd "$(dirname "$0")" && pwd)/test_env.sh" "$@"
SKILL=${1:?usage: test_read_length.sh <skill.md>}
HERE=$(cd "$(dirname "$0")" && pwd)
CODE=$(bash "$HERE/cut_block.sh" "$SKILL" "**Read-length detection.**") || { echo "FAIL: read-length block not found"; exit 1; }
RUN=$(bash "$HERE/cut_block.sh" "$SKILL" "**Running the read-length check.**") || { echo "FAIL: read-length run block not found"; exit 1; }
eval "$CODE"
declare -F detect_read_length >/dev/null || { echo "FAIL: detect_read_length is not defined by the block"; exit 1; }
declare -F read_length_report >/dev/null || { echo "FAIL: read_length_report is not defined by the block"; exit 1; }
reads() { awk -v L="$1" -v n="$2" 'BEGIN { s = ""; for (i = 0; i < L; i++) s = s "A"; q = s; gsub(/A/, "I", q); for (r = 1; r <= n; r++) printf "@r%d\n%s\n+\n%s\n", r, s, q }'; }
fq() { reads "$2" "$3" | gzip -c > "$TEST_TMP/$1"; }
fail=0
check() { local name=$1 want=$2 got; shift 2; got=$(detect_read_length "$@"); [ "$got" = "$want" ] || { echo "FAIL: case $name: expected '$want', got '$got'"; fail=1; }; }
# check_err <name> <file>...: detect_read_length must print one ERROR line and return 1
check_err() { local name=$1 got rc; shift; got=$(detect_read_length "$@"); rc=$?
  [ $rc -ne 0 ] && [[ $got == ERROR:* ]] && [ "$(printf '%s\n' "$got" | wc -l)" -eq 1 ] || { echo "FAIL: case $name: expected one ERROR line and a non-zero exit, got '$got' (exit $rc)"; fail=1; }; }
# report <name> <rows file> <expected line>...: read_length_report must print every expected line (exact) and exit 0
report() { local name=$1 rows=$2 out w; shift 2; out=$(read_length_report < "$rows") || { echo "FAIL: case $name: report failed: $out"; fail=1; return; }
  for w in "$@"; do printf '%s\n' "$out" | grep -qxF -- "$w" || { echo "FAIL: case $name: missing line '$w' in:"; printf '%s\n' "$out"; fail=1; }; done; }
# verdict <name> <rows file> <OK|WARNING|CONFOUNDED>: exactly one verdict line, and it is the expected one
verdict() { local out v; out=$(read_length_report < "$2"); v=$(printf '%s\n' "$out" | grep -E '^(OK|WARNING|CONFOUNDED):' | cut -d: -f1 | tr '\n' ' ')
  [ "$v" = "$3 " ] || { echo "FAIL: case $1: verdict '$v', expected '$3'"; printf '%s\n' "$out"; fail=1; }; }
row() { printf '%s\t%s\t%s\t%s\n' "$1" "$2" "$TEST_TMP/$3" "${4:+$TEST_TMP/$4}"; }

# detect_read_length on files
fq a.fastq.gz 101 30; fq b.fastq.gz 101 30; fq c.fastq.gz 76 50; fq d.fq.gz 151 10; fq e.fastq.gz 76 30
{ reads 60 1000; reads 90 2000; } | gzip -c > "$TEST_TMP/big.fastq.gz"
check single "101 30" "$TEST_TMP/a.fastq.gz"
check mixed "101 60" "$TEST_TMP/a.fastq.gz" "$TEST_TMP/b.fastq.gz" "$TEST_TMP/c.fastq.gz"
check tie "101 30" "$TEST_TMP/a.fastq.gz" "$TEST_TMP/e.fastq.gz"
check fq_gz "151 10" "$TEST_TMP/d.fq.gz"
check first1000 "60 1000" "$TEST_TMP/big.fastq.gz"
reads 75 20 | sed 's/$/\r/' | gzip -c > "$TEST_TMP/crlf.fastq.gz"
check crlf "75 20" "$TEST_TMP/crlf.fastq.gz"
: > "$TEST_TMP/empty.fastq.gz"; printf '' | gzip -c > "$TEST_TMP/emptygz.fastq.gz"
reads 101 200 | gzip -c > "$TEST_TMP/full.gz"; head -c "$(( $(wc -c < "$TEST_TMP/full.gz") * 2 / 3 ))" "$TEST_TMP/full.gz" > "$TEST_TMP/trunc.fastq.gz"
printf 'this is not\na fastq file\n' | gzip -c > "$TEST_TMP/text.fastq.gz"
check_err empty "$TEST_TMP/empty.fastq.gz"
check_err empty_gz "$TEST_TMP/emptygz.fastq.gz"
check_err missing "$TEST_TMP/nonexistent.fastq.gz"
check_err truncated "$TEST_TMP/trunc.fastq.gz"
check_err not_fastq "$TEST_TMP/text.fastq.gz"
check_err no_args
check_err one_bad_of_two "$TEST_TMP/a.fastq.gz" "$TEST_TMP/empty.fastq.gz"

# read_length_report on samplesheet rows
fq L150a.fastq.gz 150 30; fq L150b.fastq.gz 150 30; fq L100a.fastq.gz 100 30; fq L100b.fastq.gz 100 30
fq L150c.fastq.gz 150 30; fq L150d.fastq.gz 150 30; fq L100c.fastq.gz 100 30; fq L100d.fastq.gz 100 30
{ row s1 WT L150a.fastq.gz; row s2 WT L150b.fastq.gz; row s3 KO L150c.fastq.gz; row s4 KO L150d.fastq.gz; } > "$TEST_TMP/equal.tsv"
report all_equal "$TEST_TMP/equal.tsv" "SAMPLE s1 WT R1 150" "READ_LENGTH 150" "OK: every file has 150 bp reads"
verdict all_equal "$TEST_TMP/equal.tsv" OK
# two lengths, each condition has both; the first row alone would give 100, all rows together give 150 (tie -> longer)
{ row s1 WT L100a.fastq.gz; row s2 WT L150a.fastq.gz; row s3 KO L150b.fastq.gz; row s4 KO L100b.fastq.gz; } > "$TEST_TMP/unconf.tsv"
report unconfounded "$TEST_TMP/unconf.tsv" "SAMPLE s1 WT R1 100" "SAMPLE s2 WT R1 150" "READ_LENGTH 150"
verdict unconfounded "$TEST_TMP/unconf.tsv" WARNING
# condition WT at 150, KO at 100
{ row s1 WT L150a.fastq.gz; row s2 WT L150b.fastq.gz; row s3 KO L100a.fastq.gz; row s4 KO L100b.fastq.gz; } > "$TEST_TMP/conf.tsv"
report confounded "$TEST_TMP/conf.tsv" "SAMPLE s3 KO R1 100" "READ_LENGTH 150" "CONFOUNDED: read lengths differ between conditions (WT: 150 bp; KO: 100 bp)"
verdict confounded "$TEST_TMP/conf.tsv" CONFOUNDED
# partly confounded: KO has no 150 bp sample
{ row s1 WT L150a.fastq.gz; row s2 WT L100a.fastq.gz; row s3 KO L100b.fastq.gz; row s4 KO L100c.fastq.gz; } > "$TEST_TMP/partconf.tsv"
verdict partly_confounded "$TEST_TMP/partconf.tsv" CONFOUNDED
# paired-end, R1 150 and R2 100 in every sample: lengths differ, not confounded, READ_LENGTH from R1
{ row s1 WT L150a.fastq.gz L100a.fastq.gz; row s2 WT L150b.fastq.gz L100b.fastq.gz; row s3 KO L150c.fastq.gz L100c.fastq.gz; row s4 KO L150d.fastq.gz L100d.fastq.gz; } > "$TEST_TMP/pairedr2.tsv"
report paired_r2 "$TEST_TMP/pairedr2.tsv" "SAMPLE s1 WT R1 150 R2 100" "READ_LENGTH 150"
verdict paired_r2 "$TEST_TMP/pairedr2.tsv" WARNING
# paired-end, only KO's R2 is shorter: confounded
{ row s1 WT L150a.fastq.gz L150b.fastq.gz; row s2 WT L150c.fastq.gz L150d.fastq.gz; row s3 KO L150a.fastq.gz L100a.fastq.gz; row s4 KO L150b.fastq.gz L100b.fastq.gz; } > "$TEST_TMP/pairedconf.tsv"
verdict paired_r2_confounded "$TEST_TMP/pairedconf.tsv" CONFOUNDED
# an unreadable file: ERROR lines, no READ_LENGTH, non-zero exit
for bad in empty.fastq.gz trunc.fastq.gz; do
  { row s1 WT L150a.fastq.gz; row s2 WT "$bad"; row s3 KO L150c.fastq.gz; row s4 KO L150d.fastq.gz; } > "$TEST_TMP/bad.tsv"
  out=$(read_length_report < "$TEST_TMP/bad.tsv"); rc=$?
  { [ $rc -ne 0 ] && printf '%s\n' "$out" | grep -q '^ERROR: read length not set' && ! printf '%s\n' "$out" | grep -q '^READ_LENGTH'; } \
    || { echo "FAIL: case report_error_$bad: expected ERROR, no READ_LENGTH and a non-zero exit, got (exit $rc):"; printf '%s\n' "$out"; fail=1; }
done
# paired R2 unreadable
{ row s1 WT L150a.fastq.gz empty.fastq.gz; row s2 WT L150b.fastq.gz L150c.fastq.gz; } > "$TEST_TMP/badr2.tsv"
out=$(read_length_report < "$TEST_TMP/badr2.tsv"); rc=$?
{ [ $rc -ne 0 ] && ! printf '%s\n' "$out" | grep -q '^READ_LENGTH'; } || { echo "FAIL: case report_error_r2: expected a non-zero exit and no READ_LENGTH (exit $rc)"; fail=1; }

# the run lines read the samplesheet columns by name (CRLF and an empty line tolerated)
mkdir -p "$TEST_TMP/wd/fq" && cp "$TEST_TMP/L150a.fastq.gz" "$TEST_TMP/L100a.fastq.gz" "$TEST_TMP/L150b.fastq.gz" "$TEST_TMP/L100b.fastq.gz" "$TEST_TMP/wd/fq/"
printf 'sample,fastq_1,fastq_2,strandedness,condition\r\nA1,fq/L150a.fastq.gz,fq/L100a.fastq.gz,reverse,WT\r\nB1,fq/L150b.fastq.gz,fq/L100b.fastq.gz,reverse,KO\r\n\r\n' > "$TEST_TMP/wd/x_samplesheet.csv"
RUNF=$(printf '%s\n' "$RUN" | sed -e "s|{CWD}|$TEST_TMP/wd|g" -e 's|{SAMPLESHEET_CSV}|x_samplesheet.csv|g')
out=$( eval "$RUNF" )
for w in "SAMPLE A1 WT R1 150 R2 100" "SAMPLE B1 KO R1 150 R2 100" "READ_LENGTH 150"; do
  printf '%s\n' "$out" | grep -qxF -- "$w" || { echo "FAIL: case run_lines: missing line '$w' in:"; printf '%s\n' "$out"; fail=1; }
done
forbidden_ran && fail=1
[ $fail -eq 0 ] && echo "READ LENGTH PASS" || exit 1
