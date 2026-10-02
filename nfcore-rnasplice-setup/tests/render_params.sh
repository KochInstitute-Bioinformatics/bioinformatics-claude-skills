#!/bin/bash
# Usage: render_params.sh <skill.md> <values.tsv> <out.yaml>
# Builds the params file exactly as Step 11 describes: the Step 11 template (block after "**Params file template.**") with its
# {MODULE_PARAMS} line replaced by the Step 8 block (after "**Module and option keys (written to the params file).**"), then
# every {NAME} replaced by its value from values.tsv (NAME<TAB>VALUE). The star_index line is deleted when its value is empty
# (BAM input, or no reusable STAR index); any other missing or empty value, a value containing {...}, or a repeated key is an
# error, and then no output file is written. Prints "RENDERED <out.yaml>" or exits 1.
set -u
SKILL=${1:?usage: render_params.sh <skill.md> <values.tsv> <out.yaml>}
VALS=${2:?usage: render_params.sh <skill.md> <values.tsv> <out.yaml>}
OUT=${3:?usage: render_params.sh <skill.md> <values.tsv> <out.yaml>}
HERE=$(cd "$(dirname "$0")" && pwd)
T=$(mktemp -d) || exit 1
trap 'rm -f "$T/base" "$T/mod" "$T/joined" "$T/out"; rmdir "$T"' EXIT
case "$T" in /tmp/tmp.?*) ;; *) echo "unexpected temp dir: $T"; exit 1 ;; esac
bash "$HERE/cut_block.sh" "$SKILL" "**Params file template.**" > "$T/base" || { echo "ERROR: params file template not found"; exit 1; }
bash "$HERE/cut_block.sh" "$SKILL" "**Module and option keys (written to the params file).**" > "$T/mod" || { echo "ERROR: module block not found"; exit 1; }
[ "$(grep -cx '{MODULE_PARAMS}' "$T/base")" -eq 1 ] || { echo "ERROR: the template must contain exactly one {MODULE_PARAMS} line"; exit 1; }
awk -v mod="$T/mod" '$0 == "{MODULE_PARAMS}" { while ((getline l < mod) > 0) print l; next } { print }' "$T/base" > "$T/joined"
awk -F'\t' '
  NR == FNR {
    if ($1 != "") { has[$1] = 1; v[$1] = $2; if ($2 ~ /\{[A-Z_]+\}/) { print "ERROR: value of " $1 " contains a placeholder" > "/dev/stderr"; bad = 1 } }
    next
  }
  $0 == "star_index: \"{STAR_INDEX}\"" && (!("STAR_INDEX" in has) || v["STAR_INDEX"] == "") { next }
  { line = $0
    while (match(line, /\{[A-Z_]+\}/)) {
      k = substr(line, RSTART + 1, RLENGTH - 2)
      if (!(k in has) || v[k] == "") { print "ERROR: no value for {" k "}" > "/dev/stderr"; bad = 1; line = substr(line, 1, RSTART - 1) substr(line, RSTART + RLENGTH); continue }
      line = substr(line, 1, RSTART - 1) v[k] substr(line, RSTART + RLENGTH)
    }
    if (match(line, /^[A-Za-z_0-9]+:/)) { key = substr(line, 1, RLENGTH - 1); if (key in seen) { print "ERROR: repeated key " key > "/dev/stderr"; bad = 1 }; seen[key] = 1 }
    print line
  }
  END { exit bad }' "$VALS" "$T/joined" > "$T/out" || exit 1
mv "$T/out" "$OUT" && echo "RENDERED $OUT"
