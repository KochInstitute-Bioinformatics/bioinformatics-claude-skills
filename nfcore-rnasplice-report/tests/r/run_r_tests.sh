#!/bin/bash
#SBATCH -J rsr_r_tests
#SBATCH -N 1 -n 1 -c 2 --mem=8G -t 2:00:00
# Usage (the log goes where -o says, never into the repository):
#   sbatch [-p <partition>] -o <scratch>/logs/rsr_r_%j.out <tests dir>/r/run_r_tests.sh <skill.md> <tests dir>  [<one test script>]
# The optional third argument (for example test_enrich_summary.R) runs that script only; without it all scripts run (the walltime is 2 h for the full set).
# Both paths absolute and on the shared file system (the job runs on a compute node). tests/test_r.sh submits this script with
# sbatch and waits for it; prove_mutations.sh (mode rtest) does the same on a mutated copy of the skill.
# Runs tests/r/test_harmonise.R, tests/r/test_events_report.R, tests/r/test_usage_cross.R, tests/r/test_enrich_summary.R and tests/r/test_real_hipscm.R in the R image, one after the other. It refuses to run outside a SLURM job: R never runs on a login node.
# The log holds one "R TESTS PASS (<n>)" or "R TESTS FAIL (...)" line per script, then "EXIT <code>".
set -uo pipefail
[ -n "${SLURM_JOB_ID:-}" ] || { echo "ERROR: run_r_tests.sh runs R: submit it with sbatch, never run it on a login node"; exit 1; }
SKILL=${1:?usage: run_r_tests.sh <skill.md> <tests dir> [<one test script>]}; TESTS=${2:?usage: run_r_tests.sh <skill.md> <tests dir> [<one test script>]}
case "${3:-}" in ""|test_harmonise.R|test_events_report.R|test_usage_cross.R|test_enrich_summary.R|test_real_hipscm.R) ;; *) echo "ERROR: unknown test script: ${3} (one of test_harmonise.R test_events_report.R test_usage_cross.R test_enrich_summary.R test_real_hipscm.R, or none for all)"; exit 1 ;; esac
# A job submitted from a clean environment (env -i) has no module function: load Lmod's bash init (the profile script of the
# cluster does nothing inside a SLURM job).
CONF=${BIOINFO_SKILLS_CONF:-$HOME/.claude/bioinfo-skills.conf}
conf() { sed -n "s/^$1=//p" "$CONF" 2>/dev/null | tail -n 1; }
if ! type module > /dev/null 2>&1; then
  [ -z "$(conf MODULEPATH)" ] || export MODULEPATH=$(conf MODULEPATH)
  [ -z "$(conf MODULE_INIT)" ] || . "$(conf MODULE_INIT)" || exit 1
fi
eval "$(conf SINGULARITY_LOAD)" || exit 1   # the user's own site config; empty = Singularity is already on the PATH
SIF=${RSR_SIF:-$(conf SIF_DIR)/bulkrnaseq_latest.sif}
BIND="$(dirname "$SKILL"),$TESTS"; [ -z "$(conf SINGULARITY_BIND)" ] || BIND="$(conf SINGULARITY_BIND),$BIND"
rc=0
for t in test_harmonise.R test_events_report.R test_usage_cross.R test_enrich_summary.R test_real_hipscm.R; do
  [ -z "${3:-}" ] || [ "$t" = "$3" ] || continue
  echo "== $t"
  singularity exec --bind "$BIND" "$SIF" Rscript "$TESTS/r/$t" "$SKILL" "$TESTS/fixtures" || rc=1
done
echo "EXIT $rc"; exit $rc
