#!/bin/bash
# Usage: test_r.sh <skill.md>
# Runs the R tests of the skill (the five scripts tests/r/test_*.R, see tests/r/run_r_tests.sh) in the R image on a compute node: submits tests/r/run_r_tests.sh
# with sbatch --wait on the site partition (RSR_PARTITION, else SLURM_PARTITION of ~/.claude/bioinfo-skills.conf, else none) and prints the job's log. This script itself only calls sbatch; R never runs here.
# Refuses to run when sbatch is not on the PATH. The skill file must be on the shared file system (the job reads it).
# The job log goes to a fresh directory under RSR_SHARED_TMP (required: a scratch directory visible from the compute nodes),
# removed at exit.
# Part of run_all_tests.sh --full only (it takes minutes); prove_mutations.sh uses it as the runner of rtest rows.
# Prints the log and "R SUITE PASS" (exit 0), or "R SUITE FAIL ..." (exit 1).
set -u
SKILL=${1:?usage: test_r.sh <skill.md>}
case "$SKILL" in /*) ;; *) SKILL="$PWD/$SKILL" ;; esac
HERE=$(cd "$(dirname "$0")" && pwd)
type -P sbatch > /dev/null || { echo "R SUITE FAIL: sbatch not found on the PATH; the R tests run only through sbatch, never on a login node"; exit 1; }
[ -s "$SKILL" ] || { echo "R SUITE FAIL: skill file missing or empty: $SKILL"; exit 1; }
BASE=${RSR_SHARED_TMP:-}
[ -n "$BASE" ] || { echo "R SUITE FAIL: set RSR_SHARED_TMP to a scratch directory that the compute nodes can see"; exit 1; }
PART=${RSR_PARTITION:-$(sed -n 's/^SLURM_PARTITION=//p' "${BIOINFO_SKILLS_CONF:-$HOME/.claude/bioinfo-skills.conf}" 2>/dev/null | tail -n 1)}
mkdir -p "$BASE" || exit 1
T=$(mktemp -d "$BASE/rsr.XXXXXXXXXX") || exit 1
trap 'rm -f "$T"/rsr_r_*.out "$T/submit.out"; rmdir "$T"' EXIT
case "$T" in "$BASE"/rsr.??????????) ;; *) echo "unexpected temp dir: $T"; exit 1 ;; esac
sbatch --wait ${PART:+-p "$PART"} -o "$T/rsr_r_%j.out" "$HERE/r/run_r_tests.sh" "$SKILL" "$HERE" > "$T/submit.out" 2>&1; rc=$?
cat "$T/submit.out"
logs=$(ls "$T"/rsr_r_*.out 2> /dev/null)
[ "$(echo "$logs" | grep -c .)" -eq 1 ] || { echo "R SUITE FAIL: expected one job log in $T, found: ${logs:-none}"; exit 1; }
cat "$logs"
if [ $rc -eq 0 ] && grep -q '^R TESTS PASS' "$logs"; then echo "R SUITE PASS"; exit 0; fi
echo "R SUITE FAIL (sbatch exit $rc)"; exit 1
