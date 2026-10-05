#!/bin/bash
# Usage: run_all_tests.sh [--full]
# Runs every static test of the skill, each under `env -i` with an explicit PATH; prints "ALL TESTS PASS ..." or exits 1.
#   default  the suites in TESTS below (about 50 minutes, bash, awk and sed only); the proof-tool self-test is SKIPPED, and the
#            last line says so
#   --full   also runs the proof-tool self-test (FULL_TESTS, about 100 minutes)
# The proof-tool self-test needs git 1.8.5 or newer (git -C); /usr/bin/git on the cluster is 1.8.3. To use a newer git, set
# RNASPLICE_TEST_GIT_DIR to the directory that holds it (for example the bin directory of a conda environment with git):
#   RNASPLICE_TEST_GIT_DIR=$HOME/.conda/envs/git-new/bin bash tests/run_all_tests.sh --full
# That directory is put first on the PATH of every test; nothing else is taken from the caller's environment except HOME.
# Every *.sh file in tests/ must be named in exactly one list below (TESTS, FULL_TESTS or NOT_TESTS); check_skill.sh and
# this script both check it, so a new test cannot be left out silently.
# test_run_all.sh tests this script on a fake tree (exit codes, PASS lines, every test run, the git version comparison).
# CHECK_RUN_ALL (read by check_skill.sh and test_run_all.sh, not by this script) is a test hook of prove_mutations.sh only: it
# points them at a mutated copy of this file.
#
# Live acceptance on the nf-core test data (not run by this script: by hand, on the cluster, the pipeline only through sbatch;
# needed after any change to a generated file, for example the nextflow.config template, changed after the real-data test):
#   1. In an empty project directory, follow the skill with the nf-core test FASTQ files (the 4 paired-end human chrX samples of
#      tests/fixtures/test_samplesheet.csv, downloaded to that directory) and the test FASTA and GTF of nf-core/test-datasets
#      (branch rnasplice; see tests/fixtures/gate_report.md, G0b) as a custom reference; choose all five analyses.
#   2. Before submitting: nextflow.config must be byte-identical to the Step 10 template
#      (bash tests/cut_block.sh <skill> "**nextflow.config template.**" | cmp - nextflow.config), and the params file must equal
#      the output of tests/render_params.sh with the same values.
#   3. sbatch the generated script. Afterwards: 'Pipeline completed successfully' in its log; every task COMPLETED in
#      {OUTDIR}/pipeline_info/execution_trace.txt; its cpus, memory and time columns equal the withName selectors for every
#      process they match; the STAR_ALIGN .command.sh files hold no --quantMode or --quantTranscriptomeSAMoutput; every path of
#      the Step 12 hand-off note exists.
#   4. Record the run in the README (Validation status) and pin it in check_skill.sh, with its mutation rows.
set -u
TESTS="check_skill.sh test_strandedness.sh test_strand_salmon.sh test_bam_policy.sh test_sample_names.sh test_validate_sheets.sh test_read_length.sh test_render_params.sh dry_run_submit.sh dry_run_helpers.sh test_run_all.sh prove_mutations.sh"
FULL_TESTS="test_proof_tools.sh"
NOT_TESTS="cut_block.sh test_env.sh render_params.sh prove_red.sh run_all_tests.sh"
HERE=$(cd "$(dirname "$0")" && pwd); SK="$HERE/../nfcore-rnasplice-setup.md"; ROOT_README="$HERE/../../README.md"
full=0
case "${1:-}" in
  "") ;;
  --full) full=1 ;;
  *) echo "usage: run_all_tests.sh [--full]"; exit 1 ;;
esac
rc=0
# Every script in tests/ is named in one of the lists.
for f in "$HERE"/*.sh; do
  b=${f##*/}
  case " $TESTS $FULL_TESTS $NOT_TESTS " in *" $b "*) ;; *) echo "FAIL: tests/$b is not named in run_all_tests.sh (TESTS, FULL_TESTS or NOT_TESTS)"; rc=1 ;; esac
done
# The PATH of every test: an optional newer git first, then the system directories.
P=/usr/bin:/bin
if [ -n "${RNASPLICE_TEST_GIT_DIR:-}" ]; then
  [ -x "$RNASPLICE_TEST_GIT_DIR/git" ] || { echo "ERROR: RNASPLICE_TEST_GIT_DIR=$RNASPLICE_TEST_GIT_DIR has no executable git"; exit 1; }
  P="$RNASPLICE_TEST_GIT_DIR:$P"
  g=$(env -i PATH="$P" /bin/bash --noprofile --norc -c 'type -P git')
  [ "$g" = "$RNASPLICE_TEST_GIT_DIR/git" ] || { echo "ERROR: git resolves to '$g', not to $RNASPLICE_TEST_GIT_DIR/git"; exit 1; }
fi
if [ $full -eq 1 ]; then
  gver=$(env -i PATH="$P" /bin/bash --noprofile --norc -c 'git --version' 2>/dev/null | awk '{print $3}')
  gpath=$(env -i PATH="$P" /bin/bash --noprofile --norc -c 'type -P git')
  if [ -z "$gver" ] || [ "$(printf '%s\n%s\n' "$gver" 1.8.5 | sort -V | head -n1)" != 1.8.5 ]; then
    echo "ERROR: the proof-tool self-test needs git 1.8.5 or newer (git -C); found '${gver:-no git}' at '${gpath:-nowhere}'."
    echo "       Set RNASPLICE_TEST_GIT_DIR to a directory holding a newer git, for example:"
    echo "       RNASPLICE_TEST_GIT_DIR=\$HOME/.conda/envs/git-new/bin bash tests/run_all_tests.sh --full"
    exit 1
  fi
  echo "git $gver ($gpath)"
fi
run() { echo "== $1"; env -i HOME="${HOME:-/nonexistent}" PATH="$P" /bin/bash --noprofile --norc "$HERE/$1" "${@:2}" || { echo "FAIL: $1"; rc=1; }; }
for t in $TESTS; do run "$t" "$SK"; done
if [ $full -eq 1 ]; then
  for t in $FULL_TESTS; do run "$t"; done
else
  echo "== SKIPPED: $FULL_TESTS (proof-tool self-test, about 100 minutes; run with --full)"
fi
grep -qF '| nf-core/rnasplice setup | `/nfcore-rnasplice-setup` |' "$ROOT_README" || { echo "FAIL: root README row missing"; rc=1; }
! grep -qF '(in development, not yet validated)' "$ROOT_README" || { echo "FAIL: root README row still marked in development"; rc=1; }
grep -qF '`/nfcore-rnasplice-setup` needs internet on the compute nodes: its job downloads the pinned nf-core/rnasplice revision and the containers' "$ROOT_README" \
  || { echo "FAIL: root README Requirements must say that /nfcore-rnasplice-setup needs internet on the compute nodes"; rc=1; }
[ $rc -eq 0 ] || { echo "SOME TESTS FAILED"; exit 1; }
if [ $full -eq 1 ]; then
  echo "ALL TESTS PASS (full run, proof-tool self-test included)"
else
  echo "ALL TESTS PASS (default run: test_proof_tools.sh was SKIPPED; run with --full to include it)"
fi
