#!/bin/bash
# Self-test of the proof tools: check_skill.sh fixture check, prove_red.sh guards, prove_mutations.sh key check.
# Usage: test_proof_tools.sh   (no arguments; uses bash, git, awk, grep, sed only; works in a mktemp -d copy under /tmp)
# Builds a two-commit scratch repository (c0: tests without a skill; c1: + the current skill) from the working tree files.
# Needs git 1.8.5 or newer (git -C); see run_all_tests.sh for RNASPLICE_REPORT_TEST_GIT_DIR.
# Prints "PROOF TOOLS PASS (<n>)" or exits 1.
set -u
HERE=$(cd "$(dirname "$0")" && pwd); SRC=$(dirname "$HERE")
T=$(mktemp -d /tmp/tmp.XXXXXXXXXX) || exit 1
trap 'case "$T" in /tmp/tmp.??????????) rm -rf -- "$T" ;; esac' EXIT
case "$T" in /tmp/tmp.??????????) ;; *) echo "unexpected temp dir: $T"; exit 1 ;; esac
R="$T/repo"; D="$R/nfcore-rnasplice-report"; G=(git -C "$R" -c user.name=selftest -c user.email=selftest@invalid -c commit.gpgsign=false)
mkdir -p "$D/tests/fixtures" "$T/tmpdir"
# Every script of tests/ is copied: run_all_tests.sh and the checker check that each is named in a list, and
# prove_mutations.sh runs the runners named in mutations.tsv on the unmutated skill first.
cp "$HERE"/*.sh "$HERE"/mutations.tsv "$D/tests/"
cp "$HERE"/fixtures/*.tsv "$D/tests/fixtures/"
git init -q "$R" || exit 1
"${G[@]}" add -A && "${G[@]}" commit -q -m c0 || exit 1
# The checker also checks the README next to the skill.
cp "$SRC/nfcore-rnasplice-report.md" "$SRC/README.md" "$D/"
cp "$SRC/../README.md" "$R/README.md"
"${G[@]}" add -A && "${G[@]}" commit -q -m c1 || exit 1
C0=$("${G[@]}" rev-parse HEAD~1); C1=$("${G[@]}" rev-parse HEAD)
cp "$D/tests/check_skill.sh" "$T/check.orig"
n=0; bad=0
# expect <exit code: 0, 1, 2 or nz (any non-zero)> <text expected in the output> <command...>
expect() { local want=$1 text=$2 rc; shift 2; n=$((n + 1))
  "$@" > "$T/out" 2>&1; rc=$?
  case "$want" in nz) [ $rc -ne 0 ] ;; *) [ $rc -eq "$want" ] ;; esac \
    || { echo "CASE $n: exit $rc, expected $want: $*"; tail -5 "$T/out"; bad=1; return; }
  grep -qF -- "$text" "$T/out" || { echo "CASE $n: output lacks '$text': $*"; tail -5 "$T/out"; bad=1; }; }
check() { bash "$D/tests/check_skill.sh" "$D/nfcore-rnasplice-report.md"; }

# 1. check_skill.sh: the unchanged tree passes; a missing fixture exits 2; a script that no list names fails the checker.
expect 0 "PASS" check
mv "$D/tests/fixtures/sign_rows.tsv" "$T/sign_rows.tsv"
expect 2 "cannot read fixture" check
mv "$T/sign_rows.tsv" "$D/tests/fixtures/sign_rows.tsv"
printf '#!/bin/bash\necho selftest\n' > "$D/tests/test_selftest_unnamed.sh"
expect 1 "FAIL: tests/test_selftest_unnamed.sh must be named in exactly one list" check
rm -f "$D/tests/test_selftest_unnamed.sh"
expect 0 "PASS" check

# 2. prove_red.sh: a bad base is refused; a base without the skill passes; no added needs pass.
expect 1 "is not a commit" bash "$D/tests/prove_red.sh" not_a_commit
expect 1 "is not a commit" bash "$D/tests/prove_red.sh" 1ccfd9x
expect 1 "is not a commit" bash "$D/tests/prove_red.sh" 0123456789abcdef0123456789abcdef01234567
expect 0 "RED PASS (no skill file at $C0" bash "$D/tests/prove_red.sh" "$C0"
expect 0 "RED PASS (0 new needs" bash "$D/tests/prove_red.sh" "$C1"
# A genuinely new need (text added to the skill) is RED.
sed -i 's/^# --- end skeleton checks$/need "SELFTEST_NEW_TEXT"\n&/' "$D/tests/check_skill.sh"; echo "SELFTEST_NEW_TEXT" >> "$D/nfcore-rnasplice-report.md"
expect 0 "RED PASS (1 new needs fail on $C1)" bash "$D/tests/prove_red.sh" "$C1"
git -C "$R" checkout -q -- nfcore-rnasplice-report
# A new need whose text was already in the skill at the base (a phrase of the title line) is NOT RED.
sed -i 's/^# --- end skeleton checks$/need "Report Skill"\n&/' "$D/tests/check_skill.sh"
expect 1 "NOT RED on $C1" bash "$D/tests/prove_red.sh" "$C1"
cp "$T/check.orig" "$D/tests/check_skill.sh"
# Added need lines of which none runs: refused.
sed -i 's/^# --- end skeleton checks$/[ 1 = 2 ] \&\& need "SELFTEST_NEVER_RUN"\n&/' "$D/tests/check_skill.sh"
expect 1 "the checker executed none of them" bash "$D/tests/prove_red.sh" "$C1"
cp "$T/check.orig" "$D/tests/check_skill.sh"
# The checker breaks (fixture missing): refused, never RED PASS.
mv "$D/tests/fixtures/sign_rows.tsv" "$T/sign_rows.tsv"
expect 1 "checker broke" bash "$D/tests/prove_red.sh" "$C1"
expect 1 "checker broke" bash "$D/tests/prove_red.sh" "$C0"
mv "$T/sign_rows.tsv" "$D/tests/fixtures/sign_rows.tsv"

# 3. prove_mutations.sh: an unknown @KEY@ fails its row; TMPDIR does not move or leak its temp dir.
cp "$D/tests/mutations.tsv" "$T/mut.orig"
printf '%s\t%s\t%s\t%s\t%s\n' T-typo check skill 's/^Stub:/Stub: --x/' 'missing required text: @NO_SUCH_KEY@' >> "$D/tests/mutations.tsv"
expect 1 "MUTATION T-typo: unknown mutation key" bash "$D/tests/prove_mutations.sh" "$D/nfcore-rnasplice-report.md"
cp "$T/mut.orig" "$D/tests/mutations.tsv"
# the default mode skips the rtest rows (R tests through sbatch) and says how many; TMPDIR does not move or leak its temp dir
nr=$(awk -F'\t' '$2 == "rtest"' "$D/tests/mutations.tsv" | grep -c .)
expect 0 "; $nr rtest rows skipped: run_all_tests.sh --full runs them)" env TMPDIR="$T/tmpdir" bash "$D/tests/prove_mutations.sh" "$D/nfcore-rnasplice-report.md"
n=$((n + 1)); [ -z "$(ls -A "$T/tmpdir")" ] || { echo "CASE $n: prove_mutations.sh left files in TMPDIR"; bad=1; }
# an unknown mode is rejected
expect 1 "usage: prove_mutations.sh <skill.md> [static|rtest]" bash "$D/tests/prove_mutations.sh" "$D/nfcore-rnasplice-report.md" bogus

[ $bad -eq 0 ] && echo "PROOF TOOLS PASS ($n)" || exit 1
