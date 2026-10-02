#!/bin/bash
# Self-test of the proof tools: check_skill.sh gate-value validation, prove_red.sh guards, prove_mutations.sh key check.
# Usage: test_proof_tools.sh   (no arguments; uses bash, git, awk, grep, sed only; works in a mktemp -d copy under /tmp)
# Builds a two-commit scratch repository (c0: tests without a skill; c1: + the current skill) from the working tree files.
# Prints "PROOF TOOLS PASS (<n>)" or exits 1.
set -u
HERE=$(cd "$(dirname "$0")" && pwd); SRC=$(dirname "$HERE")
T=$(mktemp -d /tmp/tmp.XXXXXXXXXX) || exit 1
trap 'case "$T" in /tmp/tmp.??????????) rm -rf -- "$T" ;; esac' EXIT
case "$T" in /tmp/tmp.??????????) ;; *) echo "unexpected temp dir: $T"; exit 1 ;; esac
R="$T/repo"; D="$R/nfcore-rnasplice-setup"; G=(git -C "$R" -c user.name=selftest -c user.email=selftest@invalid -c commit.gpgsign=false)
mkdir -p "$D/tests/fixtures" "$T/tmpdir"
cp "$HERE"/check_skill.sh "$HERE"/cut_block.sh "$HERE"/prove_red.sh "$HERE"/prove_mutations.sh "$HERE"/mutations.tsv "$D/tests/"
# The runners named in mutations.tsv (prove_mutations.sh runs each on the unmutated skill first).
cp "$HERE"/test_env.sh "$HERE"/test_strandedness.sh "$HERE"/test_bam_policy.sh "$HERE"/test_sample_names.sh "$HERE"/test_validate_sheets.sh "$HERE"/test_read_length.sh \
   "$HERE"/render_params.sh "$HERE"/test_render_params.sh "$HERE"/dry_run_submit.sh "$HERE"/dry_run_helpers.sh "$D/tests/"
cp "$HERE"/fixtures/* "$D/tests/fixtures/"
git init -q "$R" || exit 1
"${G[@]}" add -A && "${G[@]}" commit -q -m c0 || exit 1
cp "$SRC/nfcore-rnasplice-setup.md" "$D/"
"${G[@]}" add -A && "${G[@]}" commit -q -m c1 || exit 1
C0=$("${G[@]}" rev-parse HEAD~1); C1=$("${G[@]}" rev-parse HEAD)
GV="$D/tests/fixtures/gate_values.tsv"; cp "$GV" "$T/gv.orig"; cp "$D/tests/check_skill.sh" "$T/check.orig"
n=0; bad=0
# expect <exit code: 0, 1, 2 or nz (any non-zero)> <text expected in the output> <command...>
expect() { local want=$1 text=$2 rc; shift 2; n=$((n + 1))
  "$@" > "$T/out" 2>&1; rc=$?
  case "$want" in nz) [ $rc -ne 0 ] ;; *) [ $rc -eq "$want" ] ;; esac \
    || { echo "CASE $n: exit $rc, expected $want: $*"; tail -5 "$T/out"; bad=1; return; }
  grep -qF -- "$text" "$T/out" || { echo "CASE $n: output lacks '$text': $*"; tail -5 "$T/out"; bad=1; }; }
check() { bash "$D/tests/check_skill.sh" "$D/nfcore-rnasplice-setup.md"; }
setgv() { cp "$T/gv.orig" "$GV"; awk -F'\t' -v OFS='\t' -v k="$1" -v v="$2" '$1 == k {$2 = v} {print}' "$T/gv.orig" > "$GV"; }

# 1. check_skill.sh: the unchanged fixtures pass; each invalid gate value exits 2 naming the key.
expect 0 "PASS" check
for kv in "GATE_OUTCOME=X" "GATE_OUTCOME=STOP" "NEXTFLOW_MAX_EXCL=None" "NEXTFLOW_MAX_EXCL=26.04.0" "PIPELINE_REVISION=1b44723" \
          "PIPELINE_REVISION=1.0.4" "VERSION_TAG=dev-1b44724" "HAS_MAX_PARAMS=No" "HAS_RESOURCE_LIMITS=true" "SALMON_ROUTE=pseudo-only" \
          "PSEUDO_OFF_LINE=pseudo_aligner: false" "BAM_SHEET_HEADER=sample,genome_bam" "BAM_RMATS_LIBTYPE=fr-firstStrand" \
          "BAM_RMATS_READTYPE=pe" "BAM_DEXSEQ_STRAND=2" "BAM_FC_STRAND=3" "BAM_NEEDS_BAI=maybe" "RMATS_BAMLIST_ORDER=sorted" \
          "RMATS_B1_GROUP=Treatment" "STAR_VERSION_GENOME=2.7" "SALMON_INDEX_VERSION=v5" "DTU_FILTER_SCOPE=all" \
          "TEST_CONTRAST=GBR-YRI" "TEST_READ_LENGTH=0" "EDGER_DEU_FUNCTION=diffSpliceDGE glmFitt" "NEXTFLOW_TESTED=26.04" \
          "NEXTFLOW_MIN=>=26.04.0" "CONDA_ENV_TESTED=nf env"; do
  setgv "${kv%%=*}" "${kv#*=}"; expect 2 "gate value ${kv%%=*}=" check
done
cp "$T/gv.orig" "$GV"
expect 0 "PASS" check

# 2. prove_red.sh: a bad base is refused; a base without the skill passes; no added needs pass.
expect 1 "is not a commit" bash "$D/tests/prove_red.sh" not_a_commit
expect 1 "is not a commit" bash "$D/tests/prove_red.sh" 1ccfd9x
expect 1 "is not a commit" bash "$D/tests/prove_red.sh" 0123456789abcdef0123456789abcdef01234567
expect 0 "RED PASS (no skill file at $C0" bash "$D/tests/prove_red.sh" "$C0"
expect 0 "RED PASS (0 new needs" bash "$D/tests/prove_red.sh" "$C1"
# A genuinely new need (text added to the skill) is RED.
sed -i 's/^# --- end Task 1$/need "SELFTEST_NEW_TEXT"\n&/' "$D/tests/check_skill.sh"; echo "SELFTEST_NEW_TEXT" >> "$D/nfcore-rnasplice-setup.md"
expect 0 "RED PASS (1 new needs fail on $C1)" bash "$D/tests/prove_red.sh" "$C1"
git -C "$R" checkout -q -- nfcore-rnasplice-setup
# A new need whose text was already in the skill at the base is NOT RED.
sed -i 's/^# --- end Task 1$/need "Do not run \\`nextflow\\` here"\n&/' "$D/tests/check_skill.sh"
expect 1 "NOT RED on $C1" bash "$D/tests/prove_red.sh" "$C1"
cp "$T/check.orig" "$D/tests/check_skill.sh"
# Added need lines of which none runs: refused.
sed -i 's/^# --- end Task 1$/[ 1 = 2 ] \&\& need "SELFTEST_NEVER_RUN"\n&/' "$D/tests/check_skill.sh"
expect 1 "the checker executed none of them" bash "$D/tests/prove_red.sh" "$C1"
cp "$T/check.orig" "$D/tests/check_skill.sh"
# The checker breaks (fixture missing, or an invalid gate value): refused, never RED PASS.
mv "$D/tests/fixtures/output_tree.txt" "$T/output_tree.txt"
expect 1 "checker broke" bash "$D/tests/prove_red.sh" "$C1"
expect 1 "checker broke" bash "$D/tests/prove_red.sh" "$C0"
mv "$T/output_tree.txt" "$D/tests/fixtures/output_tree.txt"
setgv GATE_OUTCOME X
expect 1 "checker broke" bash "$D/tests/prove_red.sh" "$C1"
cp "$T/gv.orig" "$GV"

# 3. prove_mutations.sh: an unknown @KEY@ fails its row; TMPDIR does not move or leak its temp dir.
cp "$D/tests/mutations.tsv" "$T/mut.orig"
printf '%s\t%s\t%s\t%s\t%s\n' T-typo check skill 's/^Extract `tag_name`/Extract --x `tag_name`/' 'missing required text: @NO_SUCH_KEY@' >> "$D/tests/mutations.tsv"
expect 1 "MUTATION T-typo: unknown gate key" bash "$D/tests/prove_mutations.sh" "$D/nfcore-rnasplice-setup.md"
cp "$T/mut.orig" "$D/tests/mutations.tsv"
expect 0 "MUTATIONS PASS" env TMPDIR="$T/tmpdir" bash "$D/tests/prove_mutations.sh" "$D/nfcore-rnasplice-setup.md"
n=$((n + 1)); [ -z "$(ls -A "$T/tmpdir")" ] || { echo "CASE $n: prove_mutations.sh left files in TMPDIR"; bad=1; }

[ $bad -eq 0 ] && echo "PROOF TOOLS PASS ($n)" || exit 1
