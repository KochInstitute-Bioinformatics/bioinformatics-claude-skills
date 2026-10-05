#!/bin/bash
# Usage: test_run_all.sh <skill.md>
# Tests tests/run_all_tests.sh itself (under a second) on a fake tree: a copy of the runner, fake test scripts for every name in
# its TESTS, FULL_TESTS and NOT_TESTS lists (each logs its name, its argument, its PATH and whether a planted variable leaked,
# and fails when a flag file says so), fake git executables of chosen versions, and a fake root README.
# Checks: the exact PASS lines and exit 0; a failing suite gives exit 1, "SOME TESTS FAILED" and no PASS line, and the later
# suites still run; every test of TESTS runs once, in order, with the skill path, under env -i with PATH=[git dir:]/usr/bin:/bin;
# --full runs FULL_TESTS (default does not, and says SKIPPED); the git version comparison (1.8.3.1, 1.8.4 stop; 1.8.5, 1.10.0,
# 2.49.0 pass); RNASPLICE_TEST_GIT_DIR without git; the usage check; an unnamed tests/*.sh; the root README row and internet-clause checks.
# CHECK_RUN_ALL (test hook of prove_mutations.sh, target runall) replaces the runner under test; default: tests/run_all_tests.sh.
# Prints "RUN ALL PASS (<n>)" or exits 1.
set -u
[ "${TEST_CLEAN_ENV:-}" = 1 ] || [ -z "${CHECK_RUN_ALL:-}" ] || set -- "${1:-}" "$CHECK_RUN_ALL"
source "$(cd "$(dirname "$0")" && pwd)/test_env.sh" "$@"
SKILL=${1:?usage: test_run_all.sh <skill.md>}
HERE=$(cd "$(dirname "$0")" && pwd)
RA=${2:-$HERE/run_all_tests.sh}
[ -s "$RA" ] || { echo "FAIL: runner $RA is missing or empty"; exit 1; }
F="$TEST_TMP/tree"; D="$F/nfcore-rnasplice-setup"; L="$TEST_TMP/log"; FL="$TEST_TMP/flags"
mkdir -p "$D/tests" "$FL" || exit 1
cp "$RA" "$D/tests/run_all_tests.sh" || exit 1
echo "fake skill" > "$D/nfcore-rnasplice-setup.md"
ROW='| nf-core/rnasplice setup | `/nfcore-rnasplice-setup` | fake row |'
INET='- Internet access ...; `/nfcore-rnasplice-setup` needs internet on the compute nodes: its job downloads the pinned nf-core/rnasplice revision and the containers'
printf '# fake root README\n\n%s\n%s\n' "$ROW" "$INET" > "$F/README.md"
lst() { sed -n "s/^$1=\"\\(.*\\)\"\$/\\1/p" "$D/tests/run_all_tests.sh"; }
TESTS=$(lst TESTS); FULL=$(lst FULL_TESTS); NOT=$(lst NOT_TESTS)
[ -n "$TESTS" ] && [ -n "$FULL" ] || { echo "FAIL: cannot read TESTS or FULL_TESTS from the runner"; exit 1; }
for t in $TESTS $FULL $NOT; do
  [ "$t" = run_all_tests.sh ] && continue
  printf '#!/bin/bash\necho "%s|${1:-}|$PATH|${LEAK_PROBE:-unset}" >> "%s"\n[ -e "%s/%s" ] && exit 1\nexit 0\n' "$t" "$L" "$FL" "$t" > "$D/tests/$t"
done
mkgit() { mkdir -p "$TEST_TMP/git_$1"; printf '#!/bin/bash\necho "git version %s"\n' "$1" > "$TEST_TMP/git_$1/git"; chmod +x "$TEST_TMP/git_$1/git"; }
for v in 1.8.3.1 1.8.4 1.8.5 1.10.0 2.49.0; do mkgit "$v"; done
n=0; bad=0
# runit <expected exit> [VAR=value ...] -- [runner args]: runs the fake runner under env -i; output in $TEST_TMP/out
runit() { local want=$1 rc; shift; local ev=(); while [ "$1" != -- ]; do ev+=("$1"); shift; done; shift
  : > "$L"; n=$((n + 1))
  env -i HOME="$HOME" PATH=/usr/bin:/bin LEAK_PROBE=leaked ${ev[@]+"${ev[@]}"} /bin/bash --noprofile --norc "$D/tests/run_all_tests.sh" "$@" > "$TEST_TMP/out" 2>&1; rc=$?
  [ $rc -eq "$want" ] || { echo "FAIL: case $n (runner args '$*'${ev[*]+, ${ev[*]}}): exit $rc, expected $want"; tail -5 "$TEST_TMP/out"; bad=1; return 1; }; }
has() { grep -qxF -- "$1" "$TEST_TMP/out" || { echo "FAIL: case $n: output lacks the line: $1"; tail -5 "$TEST_TMP/out"; bad=1; }; }
hasnot() { ! grep -qF -- "$1" "$TEST_TMP/out" || { echo "FAIL: case $n: output must not contain: $1"; bad=1; }; }
ran() { cut -d'|' -f1 "$L" | tr '\n' ' ' | sed 's/ $//'; }
PASS_DEFAULT="ALL TESTS PASS (default run: test_proof_tools.sh was SKIPPED; run with --full to include it)"
PASS_FULL="ALL TESTS PASS (full run, proof-tool self-test included)"
want_list=$(echo $TESTS)

# 1. all pass, default run
runit 0 -- && { has "$PASS_DEFAULT"; has "== SKIPPED: $FULL (proof-tool self-test, about 100 minutes; run with --full)"
  [ "$(ran)" = "$want_list" ] || { echo "FAIL: case $n: tests run: '$(ran)', expected every test of TESTS once, in order: '$want_list'"; bad=1; }
  awk -F'|' -v sk="$D/tests/../nfcore-rnasplice-setup.md" '$2 != sk || $3 != "/usr/bin:/bin" || $4 != "unset"' "$L" | grep -q . \
    && { echo "FAIL: case $n: a test got another argument, another PATH or the caller's environment (not env -i):"; awk -F'|' -v sk="$D/tests/../nfcore-rnasplice-setup.md" '$2 != sk || $3 != "/usr/bin:/bin" || $4 != "unset"' "$L" | head -3; bad=1; }; }
# 2. each suite failing in turn: exit 1, no PASS line, the later suites still run
for t in $TESTS; do
  touch "$FL/$t"
  runit 1 -- && { has "FAIL: $t"; has "SOME TESTS FAILED"; hasnot "ALL TESTS PASS"
    [ "$(ran)" = "$want_list" ] || { echo "FAIL: case $n: with $t failing, tests run: '$(ran)', expected '$want_list'"; bad=1; }; }
  rm -f "$FL/$t"
done
# 3. --full with a new git: the extra suite runs after the others, exact full PASS line; the git dir is first on the PATH of each test
runit 0 RNASPLICE_TEST_GIT_DIR="$TEST_TMP/git_2.49.0" -- --full && { has "$PASS_FULL"; has "git 2.49.0 ($TEST_TMP/git_2.49.0/git)"
  [ "$(ran)" = "$want_list $FULL" ] || { echo "FAIL: case $n: --full ran '$(ran)', expected '$want_list $FULL'"; bad=1; }
  grep -q "^$FULL||$TEST_TMP/git_2.49.0:/usr/bin:/bin|unset\$" "$L" || { echo "FAIL: case $n: $FULL did not run without arguments under env -i with the git dir first on PATH"; bad=1; }; }
# 4. --full, the extra suite fails
touch "$FL/$FULL"
runit 1 RNASPLICE_TEST_GIT_DIR="$TEST_TMP/git_2.49.0" -- --full && { has "FAIL: $FULL"; hasnot "ALL TESTS PASS"; }
rm -f "$FL/$FULL"
# 5. the git version comparison: too old stops before any test; 1.8.5 and newer (version order, not text order) pass
for v in 1.8.3.1 1.8.4; do
  runit 1 RNASPLICE_TEST_GIT_DIR="$TEST_TMP/git_$v" -- --full && { has "ERROR: the proof-tool self-test needs git 1.8.5 or newer (git -C); found '$v' at '$TEST_TMP/git_$v/git'."; hasnot "ALL TESTS PASS"
    [ -z "$(ran)" ] || { echo "FAIL: case $n: tests ran although git $v is too old: $(ran)"; bad=1; }; }
done
for v in 1.8.5 1.10.0; do runit 0 RNASPLICE_TEST_GIT_DIR="$TEST_TMP/git_$v" -- --full && has "$PASS_FULL"; done
# 6. RNASPLICE_TEST_GIT_DIR without an executable git
runit 1 RNASPLICE_TEST_GIT_DIR="$TEST_TMP/nogit" -- --full && { has "ERROR: RNASPLICE_TEST_GIT_DIR=$TEST_TMP/nogit has no executable git"; [ -z "$(ran)" ] || { echo "FAIL: case $n: tests ran"; bad=1; }; }
# 7. usage
runit 1 -- --bogus && { has "usage: run_all_tests.sh [--full]"; [ -z "$(ran)" ] || { echo "FAIL: case $n: tests ran"; bad=1; }; }
# 8. a script in tests/ that no list names
printf '#!/bin/bash\n' > "$D/tests/test_unnamed_fake.sh"
runit 1 -- && { has "FAIL: tests/test_unnamed_fake.sh is not named in run_all_tests.sh (TESTS, FULL_TESTS or NOT_TESTS)"; hasnot "ALL TESTS PASS"; }
rm -f "$D/tests/test_unnamed_fake.sh"
# 9. root README: row missing; row still marked in development; the compute-node internet clause missing
printf '# fake root README\n\n%s\n' "$INET" > "$F/README.md"
runit 1 -- && { has "FAIL: root README row missing"; hasnot "ALL TESTS PASS"; }
printf '# fake root README\n\n%s\n%s\n%s\n' "$ROW" "$INET" "| x | y | (in development, not yet validated) z |" > "$F/README.md"
runit 1 -- && { has "FAIL: root README row still marked in development"; hasnot "ALL TESTS PASS"; }
printf '# fake root README\n\n%s\n' "$ROW" > "$F/README.md"
runit 1 -- && { has "FAIL: root README Requirements must say that /nfcore-rnasplice-setup needs internet on the compute nodes"; hasnot "ALL TESTS PASS"; }
printf '# fake root README\n\n%s\n%s\n' "$ROW" "$INET" > "$F/README.md"
runit 0 -- && has "$PASS_DEFAULT"

forbidden_ran && bad=1
[ $bad -eq 0 ] && echo "RUN ALL PASS ($n)" || exit 1
