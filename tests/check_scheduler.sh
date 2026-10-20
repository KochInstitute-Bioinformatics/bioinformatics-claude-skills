#!/bin/bash
# Usage: tests/check_scheduler.sh
# (a) the scheduler rules and table in every job-writing skill equal docs/scheduler-table.md;
# (c) bioinfo-skills.conf.example documents the scheduler keys.
# Prints "SCHEDULER CHECK PASS" or exits 1.
set -u
ROOT=$(cd "$(dirname "$0")/.." && pwd); fail=0
SKILLS="nfcore-rnaseq-setup/nfcore-rnaseq-setup.md nfcore-scrnaseq-setup/nfcore-scrnaseq-setup.md nfcore-rnavar-setup/nfcore-rnavar-setup.md
        nfcore-rnasplice-setup/nfcore-rnasplice-setup.md nfcore-rnasplice-report/nfcore-rnasplice-report.md
        bulk-rnaseq-pipeline/bulk-rnaseq-pipeline.md seurat-scrna-pipeline/seurat-scrna-pipeline.md"
for s in $SKILLS; do
  [ -f "$ROOT/$s" ] || { echo "FAIL: missing skill $s"; fail=1; continue; }
  grep -q '^<!-- SCHEDULER:BEGIN rules -->$' "$ROOT/$s" || { echo "FAIL: $s has no scheduler markers (run tools/embed-scheduler.sh)"; fail=1; continue; }
  bash "$ROOT/tools/embed-scheduler.sh" --check "$ROOT/$s" > /dev/null || { echo "FAIL: $s differs from docs/scheduler-table.md (run tools/embed-scheduler.sh)"; fail=1; }
done
for k in SCHEDULER QUEUE QUEUE_ON_CLI MAX_TIME MAX_MEM_GB MAX_CPUS PBS_FLAVOR; do
  grep -q "^$k=" "$ROOT/bioinfo-skills.conf.example" || { echo "FAIL: bioinfo-skills.conf.example does not document $k"; fail=1; }
done
[ $fail -eq 0 ] && echo "SCHEDULER CHECK PASS" || exit 1
