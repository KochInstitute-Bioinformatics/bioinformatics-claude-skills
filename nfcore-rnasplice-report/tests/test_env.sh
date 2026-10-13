# Sourced by the test scripts (source "$HERE/test_env.sh" "$@"): re-runs the calling script under env -i with a clean PATH,
# creates a stub directory in which every command that must never run during a test (Nextflow, conda, Java, Python, R,
# Singularity, SLURM, Lmod, STAR, Salmon, wget, curl) is a stub that fails and leaves a trace, and type-checks the stubs.
# Provides: $TEST_TMP (a mktemp -d directory, removed at exit), $STUBS, forbidden_ran (true when a forbidden stub ran).
if [ "${TEST_CLEAN_ENV:-}" != 1 ]; then
  exec env -i TEST_CLEAN_ENV=1 HOME="${HOME:-/nonexistent}" PATH=/usr/bin:/bin /bin/bash --noprofile --norc "$0" "$@"
fi
TEST_TMP=$(mktemp -d) || exit 1
case "$TEST_TMP" in /tmp/tmp.?*) ;; *) echo "unexpected temp dir: $TEST_TMP"; exit 1 ;; esac
trap 'rm -rf -- "$TEST_TMP"' EXIT
STUBS="$TEST_TMP/stubs"; mkdir -p "$STUBS" || exit 1
FORBIDDEN_CMDS="nextflow conda java python python3 R Rscript singularity sbatch srun squeue module STAR salmon wget curl"
for c in $FORBIDDEN_CMDS; do
  printf '#!/bin/bash\necho "FORBIDDEN %s $*" >> "%s/forbidden"\nexit 99\n' "$c" "$TEST_TMP" > "$STUBS/$c"; chmod +x "$STUBS/$c"
done
export PATH="$STUBS:/usr/bin:/bin"
for c in $FORBIDDEN_CMDS; do [ "$(type -P "$c")" = "$STUBS/$c" ] || { echo "stub $c does not resolve to $STUBS/$c"; exit 1; }; done
forbidden_ran() { if [ -s "$TEST_TMP/forbidden" ]; then echo "FAIL: a forbidden command ran:"; cat "$TEST_TMP/forbidden"; return 0; fi; return 1; }
