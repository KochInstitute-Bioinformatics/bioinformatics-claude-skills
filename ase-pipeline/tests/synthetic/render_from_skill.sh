#!/bin/bash
# Usage: render_from_skill.sh <skill.md> <step number, e.g. 12|13|16|17> <values.tsv> <out.Rmd>
# Cuts the ````rmd block of "## Step <N> " out of the skill, splices the Step 14 statistics block (between the marker lines) into
# the "<<< paste here the code of Step 14" line, and substitutes {NAME} from values.tsv (NAME<TAB>value); stops on leftovers.
# Exits non-zero (and removes its temporary files and any partial out.Rmd) on: no rmd block, no statistics block between the
# markers, an empty values file, a paste line left in the output, an empty output, or unsubstituted placeholders.
set -uo pipefail
SKILL=${1:?}; N=${2:?}; VALS=${3:?}; OUT=${4:?}
trap 'rm -f "$OUT.tmp" "$OUT.tmp2" "$OUT.stats"' EXIT
die() { echo "ERROR: $*" >&2; rm -f "$OUT"; exit 1; }
[ -s "$SKILL" ] || die "skill file missing or empty: $SKILL"
[ -s "$VALS" ] || die "values file missing or empty: $VALS"
awk -v n="$N" '
  index($0, "## Step " n " ") == 1 { in_step = 1; next }
  in_step && /^## Step / { in_step = 0 }
  in_step && /^````rmd$/ { grab = 1; next }
  grab && /^````$/ { grab = 0; in_step = 0; next }
  grab { print }' "$SKILL" > "$OUT.tmp" || die "awk failed"
[ -s "$OUT.tmp" ] || die "no rmd block in Step $N"
for m in begin end; do
  [ "$(grep -c "^# --- ase-stats-$m" "$SKILL")" -eq 1 ] || die "no statistics block between the ase-stats markers in $SKILL (need exactly one '# --- ase-stats-$m' line)"
done
awk '/^# --- ase-stats-begin/ {f = 1; next} /^# --- ase-stats-end/ {f = 0} f' "$SKILL" > "$OUT.stats" || die "awk failed"
[ -s "$OUT.stats" ] || die "no statistics block between the ase-stats markers in $SKILL"
awk -v sf="$OUT.stats" '/<<< paste here the code of Step 14/ { while ((getline l < sf) > 0) print l; next } { print }' "$OUT.tmp" > "$OUT.tmp2" || die "awk failed"
awk -F'\t' 'NR == FNR { sub(/\r$/, ""); if (NF >= 2) v[$1] = $2; next }
  { for (k in v) { t = "{" k "}"; while ((i = index($0, t)) > 0) $0 = substr($0, 1, i - 1) v[k] substr($0, i + length(t)) } print }' \
  "$VALS" "$OUT.tmp2" > "$OUT" || die "awk failed"
[ -s "$OUT" ] || die "empty output $OUT"
if grep -qF '<<< paste here' "$OUT"; then die "the statistics paste line is still in $OUT"; fi
if grep -nE '\{[A-Z][A-Z_0-9]*\}' "$OUT"; then die "unsubstituted placeholders above"; fi
echo "wrote $OUT"
