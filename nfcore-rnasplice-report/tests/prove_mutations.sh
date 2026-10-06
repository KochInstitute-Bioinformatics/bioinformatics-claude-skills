#!/bin/bash
# Usage: prove_mutations.sh <skill.md> [static|rtest]   (default static)
# Each row of mutations.tsv (id<TAB>runner<TAB>target<TAB>sed program<TAB>expected text) is applied to a copy of the skill, its
# README and tests/run_all_tests.sh (target: skill, readme, rootreadme (the repository root README.md) or runall; the runner sees the copy of run_all_tests.sh through
# CHECK_RUN_ALL); @KEY@ in the program and the text is replaced by the mutation key KEY. A row passes when
# the sed program changed the target, the runner exits non-zero, and the runner printed the expected text.
# An @KEY@ that is not a non-empty mutation key fails its row (it would otherwise become an empty string).
# A sed program that contains | uses another delimiter (s/.../.../ or s#...#...#), and the reverse.
# Runners: check runtest rtest (see runner() below).
# Modes:
#   static  the rows of the runners check and runtest (bash, awk and sed only; seconds). Rows of runner rtest are skipped and
#           counted in the last line. This is the mode of the default run_all_tests.sh.
#   rtest   only the rows of runner rtest, which are exercised only by run_all_tests.sh --full (or by hand): each row's mutated
#           skill is written to a fresh directory under RSR_SHARED_TMP (default
#           /net/bmc-lab3/data/bcc/ase_scratch/rnasplice_report/tmp, shared with the compute nodes) and test_r.sh submits the R
#           tests on it with sbatch --wait (R runs only on a compute node). At most 4 jobs run at a time. Refuses to run without
#           sbatch. Minutes per row.
# Prints "MUTATIONS PASS (<n>...)" or exits 1.
set -u
SKILL=${1:?usage: prove_mutations.sh <skill.md> [static|rtest]}
MODE=${2:-static}
case "$MODE" in static|rtest) ;; *) echo "usage: prove_mutations.sh <skill.md> [static|rtest]"; exit 1 ;; esac
case "$SKILL" in /*) ;; *) SKILL="$PWD/$SKILL" ;; esac
HERE=$(cd "$(dirname "$0")" && pwd); README="$(dirname "$SKILL")/README.md"; ROOT_README="$(dirname "$SKILL")/../README.md"
gv() { awk -F'\t' -v k="$1" '$1 == k {print $2; exit}' "$HERE/fixtures/mutation_keys.tsv"; }
# subst: replace every @KEY@ by its mutation key; an unknown or empty KEY is an error (return 1, message on stderr).
subst() { local s=$1 k v; while [[ $s =~ @([A-Z_0-9]+)@ ]]; do k=${BASH_REMATCH[1]}; v=$(gv "$k")
  [ -n "$v" ] || { echo "unknown or empty mutation key @$k@" >&2; return 1; }; s=${s//@$k@/$v}; done; printf '%s' "$s"; }
runner() { case "$1" in check) echo check_skill.sh ;; runtest) echo test_run_all.sh ;; rtest) echo test_r.sh ;; *) echo "" ;; esac; }
in_mode() { if [ "$MODE" = rtest ]; then [ "$1" = rtest ]; else [ "$1" != rtest ]; fi; }

if [ "$MODE" = rtest ]; then
  type -P sbatch > /dev/null || { echo "ERROR: sbatch not found on the PATH; rtest rows run the R tests only through sbatch"; exit 1; }
  BASE=${RSR_SHARED_TMP:-/net/bmc-lab3/data/bcc/ase_scratch/rnasplice_report/tmp}
  mkdir -p "$BASE" || exit 1
  T=$(mktemp -d "$BASE/rsrmut.XXXXXXXXXX") || exit 1
  trap 'for d in "$T"/row[0-9]*; do [ -d "$d" ] || continue; rm -f "$d/nfcore-rnasplice-report.md" "$d/README.md" "$d/ROOT_README.md" "$d/run_all_tests.sh" "$d/out" "$d/rc"; rmdir "$d"; done; rm -f "$T/out"; rmdir "$T"' EXIT
  case "$T" in "$BASE"/rsrmut.??????????) ;; *) echo "unexpected temp dir: $T"; exit 1 ;; esac
else
  T=$(mktemp -d /tmp/tmp.XXXXXXXXXX) || exit 1
  trap 'rm -f "$T/nfcore-rnasplice-report.md" "$T/README.md" "$T/ROOT_README.md" "$T/run_all_tests.sh" "$T/out"; rmdir "$T"' EXIT
  case "$T" in /tmp/tmp.??????????) ;; *) echo "unexpected temp dir: $T"; exit 1 ;; esac
fi
# every runner of the mode passes on the unmutated skill first
for r in $(awk -F'\t' '$1 !~ /^#/ && NF >= 5 {print $2}' "$HERE/mutations.tsv" | sort -u); do
  s=$(runner "$r"); [ -n "$s" ] || { echo "unknown runner: $r"; exit 1; }
  in_mode "$r" || continue
  bash "$HERE/$s" "$SKILL" > "$T/out" 2>&1 || { echo "$s fails on the unmutated skill:"; tail -20 "$T/out"; exit 1; }
done
# one mutated copy per row: static rows run at once; rtest rows are collected and run below
n=0; bad=0; skipped=0; k=0; ids=(); dirs=(); expects=()
while IFS=$'\t' read -r id run target prog expect; do
  case "$id" in ''|'#'*) continue ;; esac
  if ! in_mode "$run"; then skipped=$((skipped + 1)); continue; fi
  n=$((n + 1)); prog=$(subst "$prog") && expect=$(subst "$expect") || { echo "MUTATION $id: unknown mutation key"; bad=1; continue; }
  if [ "$MODE" = rtest ]; then k=$((k + 1)); D="$T/row$k"; mkdir "$D" || exit 1; else D=$T; fi
  cp "$SKILL" "$D/nfcore-rnasplice-report.md"; cp "$README" "$D/README.md" 2>/dev/null || : > "$D/README.md"
  cp "$HERE/run_all_tests.sh" "$D/run_all_tests.sh" 2>/dev/null || : > "$D/run_all_tests.sh"
  cp "$ROOT_README" "$D/ROOT_README.md" 2>/dev/null || : > "$D/ROOT_README.md"
  case "$target" in
    skill) f="$D/nfcore-rnasplice-report.md"; orig=$SKILL ;;
    readme) f="$D/README.md"; orig=$README ;;
    runall) f="$D/run_all_tests.sh"; orig="$HERE/run_all_tests.sh" ;;
    rootreadme) f="$D/ROOT_README.md"; orig=$ROOT_README ;;
    *) echo "MUTATION $id: bad target $target"; bad=1; continue ;;
  esac
  sed -i -e "$prog" "$f" || { echo "MUTATION $id: sed error"; bad=1; continue; }
  if cmp -s "$orig" "$f"; then echo "MUTATION $id: the sed program changed nothing"; bad=1; continue; fi
  s=$(runner "$run")
  if [ "$MODE" = rtest ]; then ids+=("$id"); dirs+=("$D"); expects+=("$expect"); continue; fi
  if CHECK_RUN_ALL="$D/run_all_tests.sh" CHECK_ROOT_README="$D/ROOT_README.md" bash "$HERE/$s" "$D/nfcore-rnasplice-report.md" > "$D/out" 2>&1; then echo "MUTATION $id: $s passed (it must fail)"; bad=1; continue; fi
  grep -qF -- "$expect" "$D/out" || { echo "MUTATION $id: $s failed without printing '$expect':"; tail -20 "$D/out"; bad=1; }
done < "$HERE/mutations.tsv"
# rtest rows: test_r.sh (sbatch --wait) on each mutated copy, at most 4 at a time
if [ "$MODE" = rtest ]; then
  i=0
  while [ $i -lt ${#ids[@]} ]; do
    for j in $(seq $i $((i + 3))); do
      [ $j -lt ${#ids[@]} ] || break
      ( bash "$HERE/test_r.sh" "${dirs[$j]}/nfcore-rnasplice-report.md" > "${dirs[$j]}/out" 2>&1; echo $? > "${dirs[$j]}/rc" ) &
    done
    wait; i=$((i + 4))
  done
  for j in "${!ids[@]}"; do
    rc=$(cat "${dirs[$j]}/rc" 2> /dev/null || echo missing)
    if [ "$rc" = 0 ]; then echo "MUTATION ${ids[$j]}: test_r.sh passed (it must fail)"; bad=1; continue; fi
    grep -qF -- "${expects[$j]}" "${dirs[$j]}/out" || { echo "MUTATION ${ids[$j]}: test_r.sh failed (exit $rc) without printing '${expects[$j]}':"; tail -20 "${dirs[$j]}/out"; bad=1; continue; }
    echo "MUTATION ${ids[$j]}: red as expected"
  done
fi
[ $bad -eq 0 ] || exit 1
if [ "$MODE" = rtest ]; then echo "MUTATIONS PASS ($n rtest rows)"
else echo "MUTATIONS PASS ($n; $skipped rtest rows skipped: run_all_tests.sh --full runs them)"; fi
