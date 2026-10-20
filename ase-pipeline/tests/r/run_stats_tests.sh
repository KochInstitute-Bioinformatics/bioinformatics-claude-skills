#!/bin/bash
#SBATCH -J ase_stats_tests
#SBATCH -N 1 -n 1 -c 8 --mem=16G -t 3:00:00
# Usage, from the repository root (the log goes where -o says, never into the repository):
#   sbatch [-p <partition>] -o <scratch>/logs/ase_stats_tests_%j.out ase-pipeline/tests/r/run_stats_tests.sh [stage1|stage2|stage3|all]
set -uo pipefail
# Site settings come from ~/.claude/bioinfo-skills.conf (SINGULARITY_LOAD, SIF_DIR, SINGULARITY_BIND, MODULE_INIT, MODULEPATH).
CONF=${BIOINFO_SKILLS_CONF:-$HOME/.claude/bioinfo-skills.conf}
conf() { sed -n "s/^$1=//p" "$CONF" 2>/dev/null | tail -n 1; }
if ! type module > /dev/null 2>&1; then
  [ -z "$(conf MODULEPATH)" ] || export MODULEPATH=$(conf MODULEPATH)
  [ -z "$(conf MODULE_INIT)" ] || . "$(conf MODULE_INIT)" || exit 1
fi
eval "$(conf SINGULARITY_LOAD)" || exit 1   # the user's own site config; empty = Singularity is already on the PATH
SIF=${ASE_R_SIF:-$(conf SIF_DIR)/bulkrnaseq_latest.sif}
BIND=$(pwd -P); [ -z "$(conf SINGULARITY_BIND)" ] || BIND="$(conf SINGULARITY_BIND),$BIND"
WHAT=${1:-all}; rc=0
cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
if [ "$WHAT" = stage1 ] || [ "$WHAT" = all ]; then
  singularity exec --bind "$BIND" "$SIF" Rscript ase-pipeline/tests/r/test_ase_stats.R ase-pipeline/ase-pipeline.md || rc=1
fi
if [ "$WHAT" = stage2 ] || [ "$WHAT" = all ]; then
  singularity exec --bind "$BIND" "$SIF" Rscript ase-pipeline/tests/r/test_ase_stats_stage2.R ase-pipeline/ase-pipeline.md || rc=1
fi
if [ "$WHAT" = stage3 ] || [ "$WHAT" = all ]; then
  singularity exec --bind "$BIND" "$SIF" Rscript ase-pipeline/tests/r/test_ase_stats_stage3.R ase-pipeline/ase-pipeline.md || rc=1
fi
echo "EXIT $rc"; exit $rc
