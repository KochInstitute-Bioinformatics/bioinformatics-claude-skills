#!/bin/bash
# Stub test of the STAR index block of Step 12 (build_star_index_rnavar_{REF_TAG}.sh): the module check,
# --outTmpDir next to the index, the literal-path removal of that directory before and after STAR, and the
# exit code. No real tool runs: STAR and module are stub executables on a private PATH, and the block runs
# under env -i (no exported shell function can override a stub).
# Usage: test_star_tmp.sh <skill.md> <work_dir>   (<work_dir> must not exist; it is removed file by file at the end)
set -u
SKILL=${1:?usage: test_star_tmp.sh <skill.md> <work_dir>}
W=${2:?usage: test_star_tmp.sh <skill.md> <work_dir>}
case "$W" in /*) ;; *) echo "work_dir must be an absolute path"; exit 2 ;; esac
mkdir "$W" || { echo "cannot create $W (it must not exist)"; exit 2; }
fail=0
GD="$W/genome"; IDX="$GD/index"; TMPD="$IDX/_STARtmp_rnavar_sjdb150"
cleanup() {
  rm -f "$W/bin/STAR" "$W/bin/module" "$W/helper.sh" "$W/module_rc" "$W/star_rc" "$W/star.log" "$W/out" \
        "$IDX/keep.txt" "$IDX/star_rnavar_sjdb150/keep.txt" "$TMPD/chunk" "$TMPD/old_chunk" "$W/cwd/_STARtmp/chunk"
  rmdir "$TMPD" "$W/cwd/_STARtmp" "$IDX/star_rnavar_sjdb150" "$IDX" "$GD" "$W/cwd" "$W/bin" "$W" 2>/dev/null
}
trap cleanup EXIT

block=$(awk '/^```bash/{f=1;b="";next} /^```/{if(f&&b~/module add star\//)print b; f=0;next} f{b=b $0 "\n"}' "$SKILL")
[ -n "$block" ] || { echo "FAIL: no STAR helper block (with 'module add star/') found in $SKILL"; exit 1; }
mkdir -p "$W/bin" "$W/cwd"
printf '%s' "$block" | sed "s|{GENOME_DIR}|$GD|g; s|{SJDB_OVERHANG}|150|g; s|{GTF_PATH}|$GD/g.gtf|g" > "$W/helper.sh"

cat > "$W/bin/module" <<'EOF'
#!/bin/bash
exit "$(cat "$STUB_DIR/module_rc")"
EOF
# STAR stub: logs its arguments; like STAR, refuses an existing --outTmpDir; without the option it uses ./_STARtmp
cat > "$W/bin/STAR" <<'EOF'
#!/bin/bash
echo "$*" >> "$STUB_DIR/star.log"
T=./_STARtmp
while [ $# -gt 0 ]; do [ "$1" = "--outTmpDir" ] && T=$2; shift; done
if [ -e "$T" ]; then echo "EXITING because of fatal ERROR: could not make temporary directory: $T" >&2; exit 102; fi
mkdir "$T" && : > "$T/chunk"
exit "$(cat "$STUB_DIR/star_rc")"
EOF
chmod +x "$W/bin/module" "$W/bin/STAR"

run() {  # run <module_rc> <star_rc>; sets rc
  echo "$1" > "$W/module_rc"; echo "$2" > "$W/star_rc"; rm -f "$W/star.log"
  env -i PATH="$W/bin:/usr/bin:/bin" STUB_DIR="$W" FASTA="$GD/g.fa" SA_INDEX_NBASES=6 SLURM_NTASKS=2 \
    bash -c 'cd "$1" || exit 99
             [ "$(type -P STAR)" = "$2/bin/STAR" ] && [ "$(type -t STAR)" = file ] || { echo "STUB CHECK FAILED: STAR"; exit 98; }
             [ "$(type -P module)" = "$2/bin/module" ] && [ "$(type -t module)" = file ] || { echo "STUB CHECK FAILED: module"; exit 98; }
             bash "$2/helper.sh"' _ "$W/cwd" "$W" > "$W/out" 2>&1
  rc=$?
  grep -q "STUB CHECK FAILED" "$W/out" && { echo "FAIL: stub check failed, stopping"; cat "$W/out"; exit 1; }
}
ok()  { echo "ok: $1"; }
bad() { echo "FAIL: $1"; fail=1; }

# 1. clean build: exit 0, STAR got --outTmpDir next to the index, nothing left in the working dir or the index dir
run 0 0
[ "$rc" -eq 0 ] && ok "clean build exits 0" || bad "clean build exit $rc: $(cat "$W/out")"
grep -qF -- "--outTmpDir $TMPD" "$W/star.log" 2>/dev/null && ok "STAR got --outTmpDir $TMPD" || bad "STAR did not get --outTmpDir $TMPD"
[ ! -e "$W/cwd/_STARtmp" ] && ok "no _STARtmp in the working directory" || bad "_STARtmp left in the working directory"
[ ! -e "$TMPD" ] && ok "temp dir removed after the build" || bad "temp dir left: $TMPD"
# 2. leftover temp dir of an interrupted run, plus siblings that must survive
: > "$IDX/keep.txt"; : > "$IDX/star_rnavar_sjdb150/keep.txt"; mkdir -p "$TMPD"; : > "$TMPD/old_chunk"
run 0 0
[ "$rc" -eq 0 ] && ok "re-run with a leftover temp dir exits 0" || bad "re-run with a leftover temp dir exit $rc: $(cat "$W/out")"
[ ! -e "$TMPD" ] && ok "leftover temp dir removed" || bad "leftover temp dir still there"
[ -e "$IDX/keep.txt" ] && [ -e "$IDX/star_rnavar_sjdb150/keep.txt" ] && ok "sibling files untouched" || bad "a sibling file was removed"
# 3. STAR fails: non-zero exit, clear message, temp dir removed
run 0 1
[ "$rc" -eq 1 ] && ok "STAR failure gives exit 1" || bad "STAR failure gave exit $rc"
grep -qF "ERROR: STAR genomeGenerate failed (exit 1)" "$W/out" && ok "STAR failure message" || bad "no STAR failure message: $(cat "$W/out")"
[ ! -e "$TMPD" ] && ok "temp dir removed after a failure" || bad "temp dir left after a failure"
# 4. module add fails: exit 1 before STAR runs
run 1 0
[ "$rc" -eq 1 ] && ok "module failure gives exit 1" || bad "module failure gave exit $rc"
[ ! -e "$W/star.log" ] && ok "STAR not run after a module failure" || bad "STAR ran after a module failure"
grep -qF "ERROR: cannot load star" "$W/out" && ok "module failure message" || bad "no module failure message"

[ $fail -eq 0 ] && echo "PASS" || exit 1
