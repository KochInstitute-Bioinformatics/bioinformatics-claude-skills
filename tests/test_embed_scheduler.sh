#!/bin/bash
# Usage: tests/test_embed_scheduler.sh
# Tests tools/embed-scheduler.sh on fake skills in a mktemp directory: first insertion before the "**Safety notice" line,
# idempotence, drift detection and repair, and a file without the anchor. Prints "EMBED TEST PASS (4)" or exits 1.
set -u
if [ "${EMBED_TEST_CLEAN:-}" != 1 ]; then exec env -i EMBED_TEST_CLEAN=1 HOME="${HOME:-/nonexistent}" PATH=/usr/bin:/bin /bin/bash --noprofile --norc "$0" "$@"; fi
HERE=$(cd "$(dirname "$0")" && pwd); GEN="$HERE/../tools/embed-scheduler.sh"; SRC="$HERE/../docs/scheduler-table.md"
T=$(mktemp -d /tmp/tmp.XXXXXXXXXX) || exit 1
case "$T" in /tmp/tmp.??????????) ;; *) echo "unexpected temp dir: $T"; exit 1 ;; esac
trap 'rm -f "$T"/*; rmdir "$T"' EXIT
fail=0; n=0
ok()  { n=$((n + 1)); echo "ok   $1"; }
bad() { echo "FAIL $1"; fail=1; }
[ -x "$GEN" ] || { echo "FAIL: $GEN is missing or not executable"; exit 1; }
[ -s "$SRC" ] || { echo "FAIL: $SRC is missing or empty"; exit 1; }
printf '# Fake skill\n\n## Site configuration (read first)\n\nKeys.\n\n**Safety notice — tell the user once.**\n\n---\n\n## Step 0\n' > "$T/a.md"
cp "$T/a.md" "$T/a.orig"
bash "$GEN" "$T/a.md" > "$T/out" 2>&1 || bad "first run failed: $(cat "$T/out")"
b=$(grep -n 'SCHEDULER:BEGIN rules' "$T/a.md" | head -1 | cut -d: -f1); s=$(grep -n '^\*\*Safety notice' "$T/a.md" | cut -d: -f1)
if [ -n "$b" ] && [ "$b" -lt "$s" ] && grep -q 'SCHEDULER:END table' "$T/a.md" && grep -q 'SCHEDULER:BEGIN table' "$T/a.md"; then ok "inserted before the safety notice"; else bad "inserted: blocks missing or after the safety notice"; fi
cp "$T/a.md" "$T/a.once"; bash "$GEN" "$T/a.md" > /dev/null 2>&1
if cmp -s "$T/a.md" "$T/a.once"; then ok "idempotent"; else bad "idempotent: second run changed the file"; fi
sed -i 's/^\(.*SCHEDULER:BEGIN table.*\)$/\1\nDRIFT LINE/' "$T/a.md"
if bash "$GEN" --check "$T/a.md" > /dev/null 2>&1; then bad "drift not detected by --check"; else
  bash "$GEN" "$T/a.md" > /dev/null 2>&1; if cmp -s "$T/a.md" "$T/a.once" && ! grep -q 'DRIFT LINE' "$T/a.md"; then ok "drift detected and repaired"; else bad "drift: not repaired"; fi; fi
printf '# No anchor here\n' > "$T/b.md"; bash "$GEN" "$T/b.md" > /dev/null 2>&1; rc=$?
if [ "$rc" -eq 2 ]; then ok "no anchor exits 2"; else bad "no anchor: exit $rc, expected 2"; fi
[ $fail -eq 0 ] && echo "EMBED TEST PASS ($n)" || exit 1
