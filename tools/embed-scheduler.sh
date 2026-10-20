#!/bin/bash
# Usage: tools/embed-scheduler.sh [--check] <skill.md>...
# Copies the "rules" and "table" regions of docs/scheduler-table.md (the lines from <!-- SCHEDULER:BEGIN x --> to
# <!-- SCHEDULER:END x -->, markers included) into each skill. A skill that has the markers gets its regions replaced; a skill
# without them gets both regions inserted before its first line that starts with "**Safety notice". With --check nothing is
# written and the exit status is 1 when any file would change. Exit 2: a file has neither markers nor that anchor.
set -u
SRC="$(cd "$(dirname "$0")" && pwd)/../docs/scheduler-table.md"
check=0; [ "${1:-}" = "--check" ] && { check=1; shift; }
[ $# -gt 0 ] || { echo "usage: embed-scheduler.sh [--check] <skill.md>..."; exit 2; }
[ -s "$SRC" ] || { echo "ERROR: $SRC is missing or empty" >&2; exit 2; }
T=$(mktemp -d) || exit 2
trap 'rm -f "$T"/*; rmdir "$T"' EXIT
for r in rules table; do
  awk -v b="<!-- SCHEDULER:BEGIN $r -->" -v e="<!-- SCHEDULER:END $r -->" '$0 == b {p = 1} p {print} $0 == e {p = 0}' "$SRC" > "$T/$r"
  [ -s "$T/$r" ] || { echo "ERROR: region $r not found in $SRC" >&2; exit 2; }
done
status=0
for f in "$@"; do
  [ -f "$f" ] || { echo "ERROR: no such file: $f" >&2; exit 2; }
  if grep -q '^<!-- SCHEDULER:BEGIN rules -->$' "$f"; then mode=replace
  elif grep -q '^\*\*Safety notice' "$f"; then mode=insert
  else echo "ERROR: $f has no scheduler markers and no '**Safety notice' line" >&2; exit 2; fi
  awk -v mode="$mode" -v rules="$T/rules" -v table="$T/table" '
    function emit(file,   l) { while ((getline l < file) > 0) print l; close(file) }
    mode == "insert" && !done && /^\*\*Safety notice/ { emit(rules); print ""; emit(table); print ""; done = 1 }
    mode == "replace" && $0 == "<!-- SCHEDULER:BEGIN rules -->" { emit(rules); skip = "rules"; next }
    mode == "replace" && $0 == "<!-- SCHEDULER:BEGIN table -->" { emit(table); skip = "table"; next }
    skip != "" { if ($0 == "<!-- SCHEDULER:END " skip " -->") skip = ""; next }
    { print }' "$f" > "$T/new" || exit 2
  if ! cmp -s "$f" "$T/new"; then
    if [ $check -eq 1 ]; then echo "DIFFERS: $f"; status=1; else cp "$T/new" "$f"; echo "updated: $f"; fi
  fi
done
exit $status
