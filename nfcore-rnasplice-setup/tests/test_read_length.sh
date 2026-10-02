#!/bin/bash
# Usage: test_read_length.sh <skill.md>
# Runs detect_read_length from the skill (block after "**Read-length detection.**") on generated gzipped FASTQ.
# Prints "READ LENGTH PASS" or exits 1.
set -u
source "$(cd "$(dirname "$0")" && pwd)/test_env.sh" "$@"
SKILL=${1:?usage: test_read_length.sh <skill.md>}
HERE=$(cd "$(dirname "$0")" && pwd)
CODE=$(bash "$HERE/cut_block.sh" "$SKILL" "**Read-length detection.**") || { echo "FAIL: read-length block not found"; exit 1; }
eval "$CODE"
declare -F detect_read_length >/dev/null || { echo "FAIL: detect_read_length is not defined by the block"; exit 1; }
reads() { awk -v L="$1" -v n="$2" 'BEGIN { s = ""; for (i = 0; i < L; i++) s = s "A"; q = s; gsub(/A/, "I", q); for (r = 1; r <= n; r++) printf "@r%d\n%s\n+\n%s\n", r, s, q }'; }
fq() { reads "$2" "$3" | gzip -c > "$TEST_TMP/$1"; }
fail=0
check() { local name=$1 want=$2 got; shift 2; got=$(detect_read_length "$@"); [ "$got" = "$want" ] || { echo "FAIL: case $name: expected '$want', got '$got'"; fail=1; }; }
fq a.fastq.gz 101 30; fq b.fastq.gz 101 30; fq c.fastq.gz 76 50; fq d.fq.gz 151 10; fq e.fastq.gz 76 30
{ reads 60 1000; reads 90 2000; } | gzip -c > "$TEST_TMP/big.fastq.gz"
check single "101 30" "$TEST_TMP/a.fastq.gz"
check mixed "101 60" "$TEST_TMP/a.fastq.gz" "$TEST_TMP/b.fastq.gz" "$TEST_TMP/c.fastq.gz"
check tie "101 30" "$TEST_TMP/a.fastq.gz" "$TEST_TMP/e.fastq.gz"
check fq_gz "151 10" "$TEST_TMP/d.fq.gz"
check first1000 "60 1000" "$TEST_TMP/big.fastq.gz"
forbidden_ran && fail=1
[ $fail -eq 0 ] && echo "READ LENGTH PASS" || exit 1
