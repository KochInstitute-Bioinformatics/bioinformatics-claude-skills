#!/bin/bash
# Usage: test_rtest_runner.sh [<skill.md>]   (the argument is ignored; the test builds its own tree)
# Tests the R test runner without R and without SLURM: test_r.sh and prove_mutations.sh (modes static and rtest) run under
# env -i against a stub sbatch (type-checked to resolve to the stub) that writes a fake job log where -o says and fails when the
# skill it is given contains STUB_MAKE_FAIL. R, Singularity and the other heavy commands are failing stubs (test_env.sh).
# Checks: test_r.sh submits run_r_tests.sh with sbatch --wait on partition bcc, its log inside a fresh scratch directory that is
# removed afterwards, passes and fails with the job, and refuses to run without sbatch; prove_mutations.sh skips rtest rows in
# the default mode, runs them (more than 4: two batches) in mode rtest on copies inside the shared scratch, reports a row whose
# job passes, refuses mode rtest without sbatch, and rejects an unknown mode.
# Prints "RTEST RUNNER PASS (<n>)" or exits 1.
set -u
source "$(cd "$(dirname "$0")" && pwd)/test_env.sh" "$@"
HERE=$(cd "$(dirname "$0")" && pwd); SRC=$(dirname "$HERE")
BASE="$TEST_TMP/shared"; CALLS="$TEST_TMP/calls"; mkdir -p "$BASE" "$TEST_TMP/okbin" "$TEST_TMP/minbin" || exit 1
# the stub sbatch: records its arguments, writes the job log at the -o path (%j = 4242), fails on STUB_MAKE_FAIL
cat > "$TEST_TMP/okbin/sbatch" <<EOF
#!/bin/bash
echo "\$*" >> "$CALLS"
out=""; rest=()
while [ \$# -gt 0 ]; do case "\$1" in -o) out=\$2; shift 2 ;; -p) shift 2 ;; --wait) shift ;; *) rest+=("\$1"); shift ;; esac; done
log=\${out//%j/4242}
[ -n "\$out" ] && [ \${#rest[@]} -eq 3 ] || { echo "stub sbatch: unexpected arguments"; exit 2; }
if grep -q STUB_MAKE_FAIL "\${rest[1]}"; then printf 'FAIL: stub failure\nR TESTS FAIL (1/1)\nEXIT 1\n' > "\$log"; exit 1; fi
printf 'ok   stub\nR TESTS PASS (1)\nEXIT 0\n' > "\$log"; exit 0
EOF
chmod +x "$TEST_TMP/okbin/sbatch"
for c in bash dirname; do ln -s "$(type -P "$c")" "$TEST_TMP/minbin/$c"; done
OKPATH="$TEST_TMP/okbin:$STUBS:/usr/bin:/bin"; MINPATH="$TEST_TMP/minbin"
[ "$(env -i PATH="$OKPATH" /bin/bash --noprofile --norc -c 'type -P sbatch')" = "$TEST_TMP/okbin/sbatch" ] || { echo "FAIL: the stub sbatch does not resolve first"; exit 1; }
[ -z "$(env -i PATH="$MINPATH" /bin/bash --noprofile --norc -c 'type -P sbatch')" ] || { echo "FAIL: sbatch resolves on the minimal PATH"; exit 1; }
n=0; bad=0
# runit <expected exit> <PATH> <command...>: under env -i with RSR_SHARED_TMP=$BASE; output in $TEST_TMP/out
runit() { local want=$1 p=$2 rc; shift 2; n=$((n + 1)); : > "$CALLS"
  env -i HOME="$HOME" PATH="$p" RSR_SHARED_TMP="$BASE" /bin/bash --noprofile --norc "$@" > "$TEST_TMP/out" 2>&1; rc=$?
  [ $rc -eq "$want" ] || { echo "FAIL: case $n ($*): exit $rc, expected $want"; tail -8 "$TEST_TMP/out"; bad=1; return 1; }; }
has() { grep -qF -- "$1" "$TEST_TMP/out" || { echo "FAIL: case $n: output lacks: $1"; tail -8 "$TEST_TMP/out"; bad=1; }; }
clean() { [ -z "$(ls -A "$BASE")" ] || { echo "FAIL: case $n: scratch left behind in $BASE: $(ls -A "$BASE")"; bad=1; }; }
ncalls() { grep -c . "$CALLS"; }

# a tree with copies of the skill and the tests (the real fixtures are not needed: the stub never runs R)
D="$TEST_TMP/tree/nfcore-rnasplice-report"; mkdir -p "$D/tests/fixtures" "$D/tests/r" || exit 1
cp "$SRC/nfcore-rnasplice-report.md" "$SRC/README.md" "$D/"; cp "$HERE"/*.sh "$HERE/mutations.tsv" "$D/tests/"
# the checker also checks the row of the skill in the repository root README (one level above the skill directory)
cp "$SRC/../README.md" "$TEST_TMP/tree/README.md" || exit 1
cp "$HERE"/fixtures/*.tsv "$D/tests/fixtures/"; cp "$HERE/r/run_r_tests.sh" "$D/tests/r/"
SK="$D/nfcore-rnasplice-report.md"; FAILSK="$TEST_TMP/fail_skill.md"; { cat "$SK"; echo STUB_MAKE_FAIL; } > "$FAILSK"

# 1. test_r.sh: passes with the job; sbatch --wait -p bcc, the log in a fresh scratch directory, the job script and its arguments
runit 0 "$OKPATH" "$D/tests/test_r.sh" "$SK" && { has "R SUITE PASS"; has "R TESTS PASS (1)"
  grep -qE "^--wait -p bcc -o $BASE/rsr\.[A-Za-z0-9]{10}/rsr_r_%j\.out $D/tests/r/run_r_tests\.sh $SK $D/tests\$" "$CALLS" \
    || { echo "FAIL: case $n: sbatch was not called as expected: $(cat "$CALLS")"; bad=1; }; }
clean
# 2. test_r.sh: fails with the job and prints its log
runit 1 "$OKPATH" "$D/tests/test_r.sh" "$FAILSK" && { has "R SUITE FAIL (sbatch exit 1)"; has "FAIL: stub failure"; }
clean
# 3. test_r.sh refuses without sbatch, and calls nothing
runit 1 "$MINPATH" "$D/tests/test_r.sh" "$SK" && { has "sbatch not found on the PATH"; [ "$(ncalls)" -eq 0 ] || { echo "FAIL: case $n: sbatch called"; bad=1; }; }

# prove_mutations.sh on the tree, with its own mutations.tsv: one check row and five rtest rows that turn the stub job red
M="$D/tests/mutations.tsv"
printf '# id\trunner\ttarget\tsed program\texpected text\n' > "$M"
printf 'T-check\tcheck\tskill\t1d\tmissing required text: # nf-core/rnasplice Report Skill\n' >> "$M"
for i in 1 2 3 4 5; do printf 'T-red%s\trtest\tskill\t$s/$/\\nSTUB_MAKE_FAIL %s/\tFAIL: stub failure\n' "$i" "$i" >> "$M"; done
# 4. default mode: the check row runs, the rtest rows are skipped and counted, sbatch is never called
runit 0 "$OKPATH" "$D/tests/prove_mutations.sh" "$SK" && { has "MUTATIONS PASS (1; 5 rtest rows skipped: run_all_tests.sh --full runs them)"
  [ "$(ncalls)" -eq 0 ] || { echo "FAIL: case $n: sbatch called in the default mode"; bad=1; }; }
# 5. mode rtest: the unmutated skill first, then the five rows (two batches); every mutated copy lies in the shared scratch
runit 0 "$OKPATH" "$D/tests/prove_mutations.sh" "$SK" rtest && { has "MUTATIONS PASS (5 rtest rows)"; has "MUTATION T-red5: red as expected"
  [ "$(ncalls)" -eq 6 ] || { echo "FAIL: case $n: expected 6 sbatch calls, got $(ncalls)"; bad=1; }
  [ "$(grep -cE " $BASE/rsrmut\.[A-Za-z0-9]{10}/row[1-5]/nfcore-rnasplice-report\.md " "$CALLS")" -eq 5 ] \
    || { echo "FAIL: case $n: the mutated copies were not under $BASE/rsrmut.*/row<k>/: $(cat "$CALLS")"; bad=1; }; }
clean
# 6. mode rtest: a row whose job passes is reported and fails the run
printf 'T-green\trtest\tskill\t$s/$/\\nharmless line/\tFAIL: stub failure\n' >> "$M"
runit 1 "$OKPATH" "$D/tests/prove_mutations.sh" "$SK" rtest && has "MUTATION T-green: test_r.sh passed (it must fail)"
clean
# 7. mode rtest refuses without sbatch; an unknown mode is rejected
runit 1 "$MINPATH" "$D/tests/prove_mutations.sh" "$SK" rtest && has "sbatch not found on the PATH"
runit 1 "$OKPATH" "$D/tests/prove_mutations.sh" "$SK" bogus && has "usage: prove_mutations.sh <skill.md> [static|rtest]"

# 8. run_r_tests.sh with an unknown third argument fails with a clear message and runs nothing
n=$((n + 1)); env -i HOME="$HOME" PATH="$OKPATH" SLURM_JOB_ID=1 /bin/bash --noprofile --norc "$D/tests/r/run_r_tests.sh" "$SK" "$D/tests" nosuch.R > "$TEST_TMP/out" 2>&1; rc=$?
[ $rc -ne 0 ] && has "unknown test script: nosuch.R" || { echo "FAIL: case $n: unknown script argument: exit $rc"; bad=1; }

forbidden_ran && bad=1
[ $bad -eq 0 ] && echo "RTEST RUNNER PASS ($n)" || exit 1
