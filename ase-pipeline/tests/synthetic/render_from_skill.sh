#!/bin/bash
# Usage: render_from_skill.sh <skill.md> <step number, e.g. 12|13|16|17> <values.tsv> <out.Rmd>
# Cuts the ````rmd block of "## Step <N> " out of the skill, splices the Step 14 statistics block (between the marker lines) into
# the "<<< paste here the code of Step 14" line, and substitutes {NAME} from values.tsv (NAME<TAB>value); stops on leftovers.
set -uo pipefail
SKILL=${1:?}; N=${2:?}; VALS=${3:?}; OUT=${4:?}
awk -v n="$N" '
  index($0, "## Step " n " ") == 1 { in_step = 1; next }
  in_step && /^## Step / { in_step = 0 }
  in_step && /^````rmd$/ { grab = 1; next }
  grab && /^````$/ { grab = 0; in_step = 0; next }
  grab { print }' "$SKILL" > "$OUT.tmp" || exit 1
[ -s "$OUT.tmp" ] || { echo "ERROR: no rmd block in Step $N" >&2; exit 1; }
awk '/^# --- ase-stats-begin/ {f = 1; next} /^# --- ase-stats-end/ {f = 0} f' "$SKILL" > "$OUT.stats" || exit 1
awk -v sf="$OUT.stats" '/<<< paste here the code of Step 14/ { while ((getline l < sf) > 0) print l; next } { print }' "$OUT.tmp" > "$OUT.tmp2" || exit 1
awk -F'\t' 'NR == FNR { v[$1] = $2; next } { for (k in v) { t = "{" k "}"; while ((i = index($0, t)) > 0) $0 = substr($0, 1, i - 1) v[k] substr($0, i + length(t)) } print }' \
  "$VALS" "$OUT.tmp2" > "$OUT" || exit 1
rm -f "$OUT.tmp" "$OUT.tmp2" "$OUT.stats"
if grep -nE '\{[A-Z][A-Z_0-9]*\}' "$OUT"; then echo "ERROR: unsubstituted placeholders above" >&2; exit 1; fi
echo "wrote $OUT"
