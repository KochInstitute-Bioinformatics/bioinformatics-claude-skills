#!/bin/bash
# Usage: prove_mutations.sh <skill.md>
# Each row of mutations.tsv (id<TAB>runner<TAB>target<TAB>sed program<TAB>expected text) is applied to a copy of the skill and
# its README (target: skill or readme); @KEY@ in the program and the text is replaced by the gate value KEY. A row passes when
# the sed program changed the target, the runner exits non-zero, and the runner printed the expected text.
# Runners: check strand bam validate readlen render submit helper (see runner() below).
# Prints "MUTATIONS PASS (<n>)" or exits 1.
set -u
SKILL=${1:?usage: prove_mutations.sh <skill.md>}
case "$SKILL" in /*) ;; *) SKILL="$PWD/$SKILL" ;; esac
HERE=$(cd "$(dirname "$0")" && pwd); README="$(dirname "$SKILL")/README.md"
T=$(mktemp -d) || exit 1
case "$T" in /tmp/tmp.*) ;; *) echo "unexpected temp dir: $T"; exit 1 ;; esac
trap 'rm -f "$T/nfcore-rnasplice-setup.md" "$T/README.md" "$T/out"; rmdir "$T"' EXIT
gv() { awk -F'\t' -v k="$1" '$1 == k {print $2; exit}' "$HERE/fixtures/gate_values.tsv"; }
subst() { local s=$1 k; while [[ $s =~ @([A-Z_0-9]+)@ ]]; do k=${BASH_REMATCH[1]}; s=${s//@$k@/$(gv "$k")}; done; printf '%s' "$s"; }
runner() { case "$1" in
  check) echo check_skill.sh ;; strand) echo test_strandedness.sh ;; bam) echo test_bam_policy.sh ;;
  validate) echo test_validate_sheets.sh ;; readlen) echo test_read_length.sh ;; render) echo test_render_params.sh ;;
  submit) echo dry_run_submit.sh ;; helper) echo dry_run_helpers.sh ;; *) echo "" ;; esac; }
for r in $(awk -F'\t' '$1 !~ /^#/ && NF >= 5 {print $2}' "$HERE/mutations.tsv" | sort -u); do
  s=$(runner "$r"); [ -n "$s" ] || { echo "unknown runner: $r"; exit 1; }
  bash "$HERE/$s" "$SKILL" > "$T/out" 2>&1 || { echo "$s fails on the unmutated skill:"; tail -20 "$T/out"; exit 1; }
done
n=0; bad=0
while IFS=$'\t' read -r id run target prog expect; do
  case "$id" in ''|'#'*) continue ;; esac
  n=$((n + 1)); prog=$(subst "$prog"); expect=$(subst "$expect")
  cp "$SKILL" "$T/nfcore-rnasplice-setup.md"; cp "$README" "$T/README.md" 2>/dev/null || : > "$T/README.md"
  case "$target" in
    skill) f="$T/nfcore-rnasplice-setup.md"; orig=$SKILL ;;
    readme) f="$T/README.md"; orig=$README ;;
    *) echo "MUTATION $id: bad target $target"; bad=1; continue ;;
  esac
  sed -i -e "$prog" "$f" || { echo "MUTATION $id: sed error"; bad=1; continue; }
  if cmp -s "$orig" "$f"; then echo "MUTATION $id: the sed program changed nothing"; bad=1; continue; fi
  s=$(runner "$run")
  if bash "$HERE/$s" "$T/nfcore-rnasplice-setup.md" > "$T/out" 2>&1; then echo "MUTATION $id: $s passed (it must fail)"; bad=1; continue; fi
  grep -qF -- "$expect" "$T/out" || { echo "MUTATION $id: $s failed without printing '$expect':"; tail -20 "$T/out"; bad=1; }
done < "$HERE/mutations.tsv"
[ $bad -eq 0 ] && echo "MUTATIONS PASS ($n)" || exit 1
