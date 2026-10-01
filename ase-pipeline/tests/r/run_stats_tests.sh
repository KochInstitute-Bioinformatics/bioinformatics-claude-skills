#!/bin/bash
#SBATCH -J ase_stats_tests
#SBATCH -p bcc
#SBATCH -N 1 -n 1 -c 8 --mem=16G -t 3:00:00
# Usage, from the repository root (the log goes where -o says, never into the repository):
#   sbatch -p bcc -o <scratch>/logs/ase_stats_tests_%j.out ase-pipeline/tests/r/run_stats_tests.sh [stage1|stage2|stage3|all]
set -uo pipefail
module add singularity/3.10.4 || exit 1
SIF=/net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif
WHAT=${1:-all}; rc=0
cd "${SLURM_SUBMIT_DIR:-.}" || exit 1
if [ "$WHAT" = stage1 ] || [ "$WHAT" = all ]; then
  singularity exec --bind /net/bmc-lab3 "$SIF" Rscript ase-pipeline/tests/r/test_ase_stats.R ase-pipeline/ase-pipeline.md || rc=1
fi
if [ "$WHAT" = stage2 ] || [ "$WHAT" = all ]; then
  singularity exec --bind /net/bmc-lab3 "$SIF" Rscript ase-pipeline/tests/r/test_ase_stats_stage2.R ase-pipeline/ase-pipeline.md || rc=1
fi
if [ "$WHAT" = stage3 ] || [ "$WHAT" = all ]; then
  singularity exec --bind /net/bmc-lab3 "$SIF" Rscript ase-pipeline/tests/r/test_ase_stats_stage3.R ase-pipeline/ase-pipeline.md || rc=1
fi
echo "EXIT $rc"; exit $rc
