#!/bin/bash
# Usage: sbatch --export=ALL,LITMASTER_SCRIPTS=<dir of the scripts> run_stage.sbatch <script.py | /abs/path.py> <args...>
# Runs one litmaster stage on a compute node (never on the login node). Slurm only. Add your partition on the command line (sbatch -p <partition> ...).
# LITMASTER_PY = a python that has ToolUniverse installed (required).
#SBATCH -J litmaster
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 00:40:00
#SBATCH -o logs/litmaster_%j.out
#SBATCH -e logs/litmaster_%j.err
set -euo pipefail
: "${LITMASTER_PY:?set LITMASTER_PY to a python with ToolUniverse installed}"
PY="$LITMASTER_PY"
: "${LITMASTER_SCRIPTS:?set LITMASTER_SCRIPTS to the scripts directory}"
case "$1" in /*) S="$1";; *) S="$LITMASTER_SCRIPTS/$1";; esac
echo "Job $SLURM_JOB_ID on $(hostname) $(date) | $S ${*:2}"
"$PY" -I "$S" "${@:2}"
echo "Finished $(date)"
