#!/bin/bash
# Usage: dry_run_submit.sh <skill.md>
# Cuts the Step 11 submission script (block after "**Submission script.**"), fills its placeholders as the wizard would, and
# runs it in scenarios against stubs of module, conda, singularity and nextflow (the real condainit is never sourced; the
# environment is clean (test_env.sh: env -i), so no exported Lmod `module` function can shadow the stub).
# The stub nextflow records every call: a passing scenario has exactly one `run` with the exact launch arguments and no `pull`.
# Prints "DRY RUN SUBMIT PASS" or exits 1.
set -u
source "$(cd "$(dirname "$0")" && pwd)/test_env.sh" "$@"
SKILL=${1:?usage: dry_run_submit.sh <skill.md>}
HERE=$(cd "$(dirname "$0")" && pwd)
gv() { awk -F'\t' -v k="$1" '$1 == k {print $2; exit}' "$HERE/fixtures/gate_values.tsv"; }
REV=$(gv PIPELINE_REVISION); TAG=$(gv VERSION_TAG); NF_T=$(gv NEXTFLOW_TESTED); NF_MIN=$(gv NEXTFLOW_MIN); NF_MAX=$(gv NEXTFLOW_MAX_EXCL)
S="$TEST_TMP/scen"; R="$TEST_TMP/run"; mkdir -p "$S" "$R" "$TEST_TMP/nonf" "$TEST_TMP/home" || exit 1
export STUB_STATE="$TEST_TMP/state"; mkdir -p "$STUB_STATE" || exit 1
cat > "$S/module" <<'EOF'
#!/bin/bash
echo "module $*" >> "$STUB_STATE/calls"
grep -qxF -- "$2" "$STUB_STATE/module_fail" 2>/dev/null && exit 1
exit 0
EOF
cat > "$S/conda" <<'EOF'
#!/bin/bash
echo "conda $*" >> "$STUB_STATE/calls"
[ -e "$STUB_STATE/conda_fail" ] && exit 1
exit 0
EOF
cat > "$S/singularity" <<'EOF'
#!/bin/bash
echo "singularity $*" >> "$STUB_STATE/calls"; exit 0
EOF
cat > "$S/nextflow" <<'EOF'
#!/bin/bash
case "$1" in
  -version) echo "nextflow -version" >> "$STUB_STATE/calls"
            printf '\n      N E X T F L O W\n      version %s build 5940\n      created 01-01-2026 00:00 UTC\n' "$(cat "$STUB_STATE/nf_version")" ;;
  run) echo "nextflow $* CACHE=${NXF_SINGULARITY_CACHEDIR:-unset}" >> "$STUB_STATE/calls"; [ -e "$STUB_STATE/run_fail" ] && exit 3 ;;
  *) echo "nextflow $*" >> "$STUB_STATE/calls" ;;
esac
exit 0
EOF
: > "$S/condainit"
chmod +x "$S/module" "$S/conda" "$S/singularity" "$S/nextflow"
cp "$S/module" "$S/conda" "$S/singularity" "$TEST_TMP/nonf/"
# No stub bypass: the stubs resolve first, no shell function (an exported Lmod `module`, a condainit `conda`) is inherited.
for c in module conda singularity nextflow; do
  [ "$(PATH="$S:$PATH" type -P "$c")" = "$S/$c" ] || { echo "stub $c does not resolve to $S/$c"; exit 1; }
  [ "$(PATH="$S:$PATH" bash -c "type -t $c")" = file ] || { echo "$c is not the stub file in a child shell (an inherited function would bypass it)"; exit 1; }
done
[ -z "$(env | grep '^BASH_FUNC_')" ] || { echo "exported shell functions in the test environment: $(env | grep -o '^BASH_FUNC_[^=]*')"; exit 1; }
[ -z "$(PATH="$TEST_TMP/nonf:/usr/bin:/bin" type -P nextflow)" ] || { echo "a real nextflow is on /usr/bin:/bin; the missing-nextflow case cannot be tested"; exit 1; }
bash "$HERE/cut_block.sh" "$SKILL" "**Submission script.**" > "$TEST_TMP/tmpl.sh" || { echo "FAIL: submission script block not found"; exit 1; }
[ "$(grep -c '^source /home/software/conda/miniconda3/bin/condainit ' "$TEST_TMP/tmpl.sh")" -eq 1 ] || { echo "FAIL: expected exactly one condainit source line"; exit 1; }
[ "$(grep -c '^[[:space:]]*\(source\|\.\) ' "$TEST_TMP/tmpl.sh")" -eq 1 ] || { echo "FAIL: unexpected extra source line"; exit 1; }
sed -e "s#^source /home/software/conda/miniconda3/bin/condainit #source \"$S/condainit\" #" \
    -e "s#{USER_EMAIL}#test@example.org#g; s#{CONDA_ENV}#nf-test#g; s#{VERSION_TAG}#$TAG#g; s#{VERSION}#$REV#g; s#{PARAMS_YAML}#T_params.yaml#g" \
    "$TEST_TMP/tmpl.sh" > "$R/submit.sh"
! grep -nE '(^|[^$])\{[A-Z_]+\}' "$R/submit.sh" || { echo "FAIL: placeholder left in the submission script"; exit 1; }
! grep -q '/home/software/conda' "$R/submit.sh" || { echo "FAIL: the real condainit would be sourced"; exit 1; }
LAUNCH="nextflow run nf-core/rnasplice -r $REV -c nextflow.config -profile slurm,singularity -params-file T_params.yaml"
fail=0
scenario() { # scenario <name> <nextflow version> <expected exit 0|1|nz> <expect one run: yes|no> [module_fail=<module> conda_fail no_nextflow run_fail cache=<dir>]
  local name=$1 ver=$2 erc=$3 erun=$4 p="$S:$PATH" w rc runs cache=""; shift 4
  rm -f "$STUB_STATE/calls" "$STUB_STATE/module_fail" "$STUB_STATE/conda_fail" "$STUB_STATE/run_fail"
  printf '%s\n' "$ver" > "$STUB_STATE/nf_version"
  for w in "$@"; do case "$w" in
    module_fail=*) echo "${w#module_fail=}" >> "$STUB_STATE/module_fail" ;;
    conda_fail) : > "$STUB_STATE/conda_fail" ;;
    run_fail) : > "$STUB_STATE/run_fail" ;;
    no_nextflow) p="$TEST_TMP/nonf:/usr/bin:/bin" ;;
    cache=*) cache=${w#cache=} ;;
  esac; done
  if [ -n "$cache" ]; then
    ( cd "$R" && HOME="$TEST_TMP/home" NXF_SINGULARITY_CACHEDIR="$cache" PATH="$p" bash ./submit.sh ) > "$TEST_TMP/out" 2>&1; rc=$?
  else
    ( cd "$R" && HOME="$TEST_TMP/home" PATH="$p" bash ./submit.sh ) > "$TEST_TMP/out" 2>&1; rc=$?
  fi
  case "$erc" in nz) [ "$rc" -ne 0 ] ;; *) [ "$rc" -eq "$erc" ] ;; esac \
    || { echo "FAIL: case $name: exit $rc, expected $erc"; sed 's/^/    /' "$TEST_TMP/out"; fail=1; }
  runs=$(grep -c '^nextflow run ' "$STUB_STATE/calls" 2>/dev/null); runs=${runs:-0}
  if [ "$erun" = yes ]; then
    [ "$runs" -eq 1 ] && grep -q "^$LAUNCH CACHE=" "$STUB_STATE/calls" \
      || { echo "FAIL: case $name: expected exactly one 'nextflow run' with the exact arguments"; cat "$STUB_STATE/calls" 2>/dev/null; fail=1; }
  else
    [ "$runs" -eq 0 ] || { echo "FAIL: case $name: nextflow run was called"; fail=1; }
  fi
  ! grep -q '^nextflow pull' "$STUB_STATE/calls" 2>/dev/null || { echo "FAIL: case $name: nextflow pull was called (a pinned revision is fetched on first use)"; fail=1; }
  ! grep -v '^nextflow \(run\|-version\)' "$STUB_STATE/calls" 2>/dev/null | grep -q '^nextflow' \
    || { echo "FAIL: case $name: nextflow was called with something other than -version or run"; fail=1; }
}
scenario tested "$NF_T" 0 yes
awk '/^module add miniconda3\/v4$/ {a = NR} /^conda activate nf-test$/ {b = NR} /^module add singularity\/3.10.4$/ {c = NR}
     /^nextflow -version$/ && !d {d = NR} /^nextflow run / {e = NR}
     END {exit !(a && b && c && d && e && a < b && b < c && c < d && d < e)}' "$STUB_STATE/calls" \
  || { echo "FAIL: case order: module miniconda3/v4, conda activate, module singularity/3.10.4, nextflow -version, nextflow run"; cat "$STUB_STATE/calls"; fail=1; }
grep -q "^$LAUNCH CACHE=$TEST_TMP/home/.singularity/cache\$" "$STUB_STATE/calls" && [ -d "$TEST_TMP/home/.singularity/cache" ] \
  || { echo "FAIL: case cache_default: NXF_SINGULARITY_CACHEDIR must default to \$HOME/.singularity/cache (created)"; fail=1; }
scenario cache_preset "$NF_T" 0 yes "cache=$TEST_TMP/mycache"
grep -q "^$LAUNCH CACHE=$TEST_TMP/mycache\$" "$STUB_STATE/calls" || { echo "FAIL: case cache_preset: a preset NXF_SINGULARITY_CACHEDIR must be kept"; fail=1; }
scenario at_min "$NF_MIN" 0 yes
scenario below_min 0.0.1 1 no
grep -q "older than $NF_MIN" "$TEST_TMP/out" || { echo "FAIL: case below_min: message"; fail=1; }
JB=$(echo "$NF_MIN" | awk -F. '{ if ($3 > 0) print $1 "." $2 "." $3 - 1; else if ($2 > 0) print $1 "." $2 - 1 ".99"; else print $1 - 1 ".99.99" }')
scenario just_below_min "$JB" 1 no
if [ "$NF_MAX" = none ]; then scenario far_future 99.0.0 0 yes
else scenario at_max "$NF_MAX" 1 no; scenario far_future 99.0.0 1 no; fi
scenario no_version garbage 1 no
grep -q "could not read the Nextflow version" "$TEST_TMP/out" || { echo "FAIL: case no_version: message"; fail=1; }
scenario run_fails "$NF_T" nz yes run_fail
scenario singularity_module_fails "$NF_T" 1 no module_fail=singularity/3.10.4
! grep -q '^nextflow' "$STUB_STATE/calls" || { echo "FAIL: case singularity_module_fails: nextflow was called"; fail=1; }
scenario miniconda_module_fails "$NF_T" 1 no module_fail=miniconda3/v4
! grep -q '^conda' "$STUB_STATE/calls" || { echo "FAIL: case miniconda_module_fails: conda was called"; fail=1; }
scenario conda_fails "$NF_T" 1 no conda_fail
! grep -q '^module add singularity' "$STUB_STATE/calls" || { echo "FAIL: case conda_fails: the script went on after conda activate failed"; fail=1; }
scenario nextflow_missing "$NF_T" 1 no no_nextflow
grep -q "nextflow is not in conda environment" "$TEST_TMP/out" || { echo "FAIL: case nextflow_missing: message"; fail=1; }
forbidden_ran && fail=1
[ $fail -eq 0 ] && echo "DRY RUN SUBMIT PASS" || exit 1
