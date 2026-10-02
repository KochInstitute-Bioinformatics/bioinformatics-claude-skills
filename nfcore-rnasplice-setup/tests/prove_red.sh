#!/bin/bash
# Usage: prove_red.sh <base commit>   (run from anywhere inside the repository, before committing)
# Every need/needr line added to tests/check_skill.sh since <base> (working tree included) must fail on the skill and the
# README as they were at <base>. Prints "RED PASS (...)" or exits 1.
set -u
BASE=${1:?usage: prove_red.sh <base commit>}
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(git -C "$HERE" rev-parse --show-toplevel) || exit 1
D=nfcore-rnasplice-setup
T=$(mktemp -d) || exit 1
case "$T" in /tmp/tmp.*) ;; *) echo "unexpected temp dir: $T"; exit 1 ;; esac
trap 'rm -f "$T/nfcore-rnasplice-setup.md" "$T/README.md" "$T/out" "$T/list"; rmdir "$T"' EXIT
git -C "$ROOT" show "$BASE:$D/nfcore-rnasplice-setup.md" > "$T/nfcore-rnasplice-setup.md" 2>/dev/null || : > "$T/nfcore-rnasplice-setup.md"
git -C "$ROOT" show "$BASE:$D/README.md" > "$T/README.md" 2>/dev/null || : > "$T/README.md"
if [ ! -s "$T/nfcore-rnasplice-setup.md" ]; then echo "RED PASS (no skill file at $BASE: every need fails)"; exit 0; fi
bash "$HERE/check_skill.sh" "$T/nfcore-rnasplice-setup.md" > "$T/out" 2>&1
CHECK_LIST_NEEDS=1 bash "$HERE/check_skill.sh" "$ROOT/$D/nfcore-rnasplice-setup.md" 2>/dev/null | grep '^NEED ' > "$T/list"
added=$(git -C "$ROOT" diff -U0 "$BASE" -- "$D/tests/check_skill.sh" \
  | awk '/^@@/ { split($3, p, ","); s = substr(p[1], 2); n = (p[2] == "") ? 1 : p[2]; for (i = 0; i < n; i++) print s + i }')
n=0; bad=0
while IFS= read -r line; do
  ln=${line#NEED }; ln=${ln%% *}; text=${line#NEED $ln }
  echo "$added" | grep -qx -- "$ln" || continue
  n=$((n + 1))
  grep -qxF -- "FAIL: missing required text: $text" "$T/out" || grep -qxF -- "FAIL: README missing required text: $text" "$T/out" \
    || { echo "NOT RED on $BASE (checker line $ln): $text"; bad=1; }
done < "$T/list"
[ $bad -eq 0 ] && echo "RED PASS ($n new needs fail on $BASE)" || exit 1
