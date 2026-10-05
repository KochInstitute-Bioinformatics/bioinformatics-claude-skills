#!/bin/bash
# Usage: cut_block.sh <skill.md> <anchor>
# Prints the body of the first fenced code block (``` or ```<lang>, closed by ```) that opens after the first line containing
# <anchor> (fixed string). Exits 1 when the anchor or the block is missing. Four-backtick (````rmd) fences never match.
set -uo pipefail
SKILL=${1:?usage: cut_block.sh <skill.md> <anchor>}; ANCHOR=${2:?usage: cut_block.sh <skill.md> <anchor>}
awk -v a="$ANCHOR" '!seen && index($0, a) { seen = 1; next }
  seen && !inb && /^```[a-z]*$/ { inb = 1; next }
  inb && /^```$/ { found = 1; exit }
  inb { print }
  END { exit !found }' "$SKILL"
