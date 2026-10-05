#!/bin/bash
#SBATCH -J rsr_r_tests
#SBATCH -p bcc
#SBATCH -N 1 -n 1 -c 2 --mem=8G -t 0:30:00
# Usage (the log goes where -o says, never into the repository):
#   sbatch -p bcc -o <scratch>/logs/rsr_r_%j.out <tests dir>/r/run_r_tests.sh <skill.md> <tests dir>
# Both paths absolute and on the shared file system (the job runs on a compute node). tests/test_r.sh submits this script with
# sbatch and waits for it; prove_mutations.sh (mode rtest) does the same on a mutated copy of the skill.
# Runs tests/r/test_harmonise.R in the R image. It refuses to run outside a SLURM job: R never runs on a login node.
# The log ends with "R TESTS PASS (<n>)" or "R TESTS FAIL (...)", then "EXIT <code>".
set -uo pipefail
[ -n "${SLURM_JOB_ID:-}" ] || { echo "ERROR: run_r_tests.sh runs R: submit it with sbatch, never run it on a login node"; exit 1; }
SKILL=${1:?usage: run_r_tests.sh <skill.md> <tests dir>}; TESTS=${2:?usage: run_r_tests.sh <skill.md> <tests dir>}
# A job submitted from a clean environment (env -i) has no module function: load Lmod's bash init (the profile script of the
# cluster does nothing inside a SLURM job).
if ! type module > /dev/null 2>&1; then
  export MODULEPATH=${MODULEPATH:-/opt/ohpc/pub/modulefiles:/home/software/modulefiles}
  . /opt/ohpc/admin/lmod/lmod/init/bash || exit 1
fi
module add singularity/3.10.4 || exit 1
SIF=/net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif
rc=0
singularity exec --bind /net/bmc-lab3 "$SIF" Rscript "$TESTS/r/test_harmonise.R" "$SKILL" "$TESTS/fixtures" || rc=1
echo "EXIT $rc"; exit $rc
