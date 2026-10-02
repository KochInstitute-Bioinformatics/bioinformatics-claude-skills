#!/bin/bash
# Usage: prove_red.sh <base commit>   (run from anywhere inside the repository, before committing)
# Every need/needr line added to tests/check_skill.sh since <base> (working tree included) must fail on the skill and the
# README as they were at <base>. Prints "RED PASS (...)" or exits 1.
# Exits 1 also when: <base> is not a commit; the checker exits with anything but 0 or 1 (a broken fixture or extraction);
# need/needr lines were added but none of them was executed. Limitation: a need called from inside a helper function
# reports the helper's line, not an added line, so it is not proven here (see the header of check_skill.sh).
set -u
BASE=${1:?usage: prove_red.sh <base commit>}
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(git -C "$HERE" rev-parse --show-toplevel) || exit 1
D=nfcore-rnasplice-setup
git -C "$ROOT" rev-parse --verify -q "$BASE^{commit}" > /dev/null || { echo "ERROR: base '$BASE' is not a commit of this repository"; exit 1; }
T=$(mktemp -d /tmp/tmp.XXXXXXXXXX) || exit 1
trap 'rm -f "$T/nfcore-rnasplice-setup.md" "$T/README.md" "$T/out" "$T/list" "$T/cur"; rmdir "$T"' EXIT
case "$T" in /tmp/tmp.??????????) ;; *) echo "unexpected temp dir: $T"; exit 1 ;; esac
# The checker must run to completion (exit 0 or 1) on the current skill.
CHECK_LIST_NEEDS=1 bash "$HERE/check_skill.sh" "$ROOT/$D/nfcore-rnasplice-setup.md" > "$T/cur" 2>&1; rc=$?
[ $rc -le 1 ] || { echo "ERROR: check_skill.sh exited $rc on the current skill (checker broke):"; grep -v '^NEED ' "$T/cur" | tail -5; exit 1; }
grep '^NEED ' "$T/cur" > "$T/list"
# A file absent at <base> (git cat-file -e fails on a valid commit) is checked as an empty file.
for f in nfcore-rnasplice-setup.md README.md; do
  if git -C "$ROOT" cat-file -e "$BASE:$D/$f" 2> /dev/null; then
    git -C "$ROOT" show "$BASE:$D/$f" > "$T/$f" || { echo "ERROR: cannot read $D/$f at $BASE"; exit 1; }
  else
    : > "$T/$f"
  fi
done
if [ ! -s "$T/nfcore-rnasplice-setup.md" ]; then echo "RED PASS (no skill file at $BASE: every need fails)"; exit 0; fi
bash "$HERE/check_skill.sh" "$T/nfcore-rnasplice-setup.md" > "$T/out" 2>&1; rc=$?
[ $rc -le 1 ] || { echo "ERROR: check_skill.sh exited $rc on the skill at $BASE (checker broke):"; tail -5 "$T/out"; exit 1; }
diff=$(git -C "$ROOT" diff -U0 "$BASE" -- "$D/tests/check_skill.sh")
added=$(echo "$diff" | awk '/^@@/ { split($3, p, ","); s = substr(p[1], 2); n = (p[2] == "") ? 1 : p[2]; for (i = 0; i < n; i++) print s + i }')
added_needs=$(echo "$diff" | grep -v '^+++' | grep '^+' | grep -cE '(^\+|[^A-Za-z0-9_])needr?[[:space:]]+"')
n=0; bad=0
while IFS= read -r line; do
  ln=${line#NEED }; ln=${ln%% *}; text=${line#NEED $ln }
  echo "$added" | grep -qx -- "$ln" || continue
  n=$((n + 1))
  grep -qxF -- "FAIL: missing required text: $text" "$T/out" || grep -qxF -- "FAIL: README missing required text: $text" "$T/out" \
    || { echo "NOT RED on $BASE (checker line $ln): $text"; bad=1; }
done < "$T/list"
if [ "$added_needs" -gt 0 ] && [ $n -eq 0 ]; then
  echo "ERROR: $added_needs need/needr lines added since $BASE, but the checker executed none of them"; bad=1
fi
[ $bad -eq 0 ] && echo "RED PASS ($n new needs fail on $BASE)" || exit 1
