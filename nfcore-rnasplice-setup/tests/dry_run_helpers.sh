#!/bin/bash
# Usage: dry_run_helpers.sh <skill.md>
# Runs the genome download helper (block after "**Genome download helper.**") against a wget stub (no network) and a logging
# gunzip stub that hands over to the system gzip, with the system gzip for `gzip -t`, in fresh, re-run, corrupt, truncated,
# failing and relative-path scenarios. Every wget call must resume (-c); every gunzip call must write to stdout (-c, never -k).
# Prints "DRY RUN HELPERS PASS" or exits 1.
set -u
source "$(cd "$(dirname "$0")" && pwd)/test_env.sh" "$@"
SKILL=${1:?usage: dry_run_helpers.sh <skill.md>}
HERE=$(cd "$(dirname "$0")" && pwd)
S="$TEST_TMP/scen"; mkdir -p "$S" || exit 1
export STUB_STATE="$TEST_TMP/state"; mkdir -p "$STUB_STATE" || exit 1
case "$(type -P gzip)" in /usr/bin/gzip|/bin/gzip) SYSGZIP=$(type -P gzip) ;; *) echo "gzip is not the system gzip"; exit 1 ;; esac
cat > "$S/wget" <<'EOF'
#!/bin/bash
# wget stub: wget -c -O <out> <url>; per-URL mode in $STUB_STATE/wget_mode ("<url> fail" or "<url> truncated")
echo "wget $*" >> "$STUB_STATE/calls"
out=""; url=""
while [ $# -gt 0 ]; do case "$1" in -O) out=$2; shift 2 ;; -*) shift ;; *) url=$1; shift ;; esac; done
mode=$(awk -v u="$url" '$1 == u {print $2}' "$STUB_STATE/wget_mode" 2>/dev/null)
[ "$mode" = fail ] && exit 4
case "$url" in
  *.fa.gz) printf '>1\nACGTACGT\n' ;;
  *.gtf.gz) printf '1\ttest\texon\t1\t8\t.\t+\t.\tgene_id "G1";\n' ;;
esac | gzip -c > "$out"
if [ "$mode" = truncated ]; then head -c 10 "$out" > "$out.t" && mv "$out.t" "$out"; fi
exit 0
EOF
printf '#!/bin/bash\necho "gunzip $*" >> "$STUB_STATE/calls"\nexec %s -d "$@"\n' "$SYSGZIP" > "$S/gunzip"
chmod +x "$S/wget" "$S/gunzip"; export PATH="$S:$PATH"
for c in wget gunzip; do
  [ "$(type -P "$c")" = "$S/$c" ] || { echo "stub $c does not resolve to $S/$c"; exit 1; }
  [ "$(bash -c "type -t $c")" = file ] || { echo "$c is not the stub file in a child shell (an inherited function would bypass it)"; exit 1; }
done
[ -z "$(env | grep '^BASH_FUNC_')" ] || { echo "exported shell functions in the test environment"; exit 1; }
bash "$HERE/cut_block.sh" "$SKILL" "**Genome download helper.**" > "$TEST_TMP/helper.tmpl" || { echo "FAIL: download helper block not found"; exit 1; }
G="$TEST_TMP/genome"; FU=https://example.invalid/genome.fa.gz; GU=https://example.invalid/genes.gtf.gz
mk_helper() { # mk_helper <GENOME_DIR as written into the script>
  sed -e "s#{GENOME_DIR}#$1#g; s#{FASTA_PATH}#$G/genome.fa#g; s#{GTF_PATH}#$G/genes.gtf#g" \
      -e "s#{ENSEMBL_FASTA_URL}#$FU#g; s#{ENSEMBL_GTF_URL}#$GU#g; s#{USER_EMAIL}#test@example.org#g; s#{REF_TAG}#test#g" \
      "$TEST_TMP/helper.tmpl" > "$TEST_TMP/helper.sh"
  ! grep -nE '(^|[^$])\{[A-Z_]+\}' "$TEST_TMP/helper.sh" || { echo "FAIL: placeholder left in the helper"; exit 1; }
}
fail=0
runh() { rm -f "$STUB_STATE/calls"; ( cd "$TEST_TMP" && bash ./helper.sh ) > "$TEST_TMP/out" 2>&1; local rc=$?
  ! grep '^wget ' "$STUB_STATE/calls" 2>/dev/null | grep -vq '^wget -c ' || { echo "FAIL: wget was called without -c (no resume)"; fail=1; }
  ! grep '^gunzip ' "$STUB_STATE/calls" 2>/dev/null | grep -vq '^gunzip -c ' || { echo "FAIL: gunzip was called without -c"; fail=1; }
  return $rc; }
nwget() { local n; n=$(grep -c '^wget ' "$STUB_STATE/calls" 2>/dev/null); echo "${n:-0}"; }
mk_helper "$G"; mkdir -p "$G"; echo keep > "$G/sentinel"
runh; rc=$?
{ [ $rc -eq 0 ] && grep -qx '>1' "$G/genome.fa" && grep -q 'gene_id "G1"' "$G/genes.gtf" && [ -s "$G/genome.fa.gz" ] && [ -s "$G/genes.gtf.gz" ] \
  && [ "$(nwget)" -eq 2 ] && [ -z "$(ls "$G" | grep '\.part$')" ] && [ -s "$G/sentinel" ]; } || { echo "FAIL: case fresh"; cat "$TEST_TMP/out"; fail=1; }
runh; rc=$?; { [ $rc -eq 0 ] && [ "$(nwget)" -eq 0 ]; } || { echo "FAIL: case rerun"; fail=1; }
rm -f "$G/genome.fa"; runh; rc=$?; { [ $rc -eq 0 ] && [ "$(nwget)" -eq 0 ] && grep -qx '>1' "$G/genome.fa"; } || { echo "FAIL: case gz_present"; fail=1; }
rm -f "$G/genes.gtf"; printf 'garbage' > "$G/genes.gtf.gz"; runh; rc=$?
{ [ $rc -eq 0 ] && [ "$(nwget)" -eq 1 ] && grep -q 'gene_id "G1"' "$G/genes.gtf"; } || { echo "FAIL: case corrupt_gz"; fail=1; }
rm -f "$G/genes.gtf" "$G/genes.gtf.gz"; echo "$GU truncated" > "$STUB_STATE/wget_mode"; runh; rc=$?
{ [ $rc -ne 0 ] && grep -q 'is not a valid gzip file' "$TEST_TMP/out" && [ ! -e "$G/genes.gtf" ] && [ ! -e "$G/genes.gtf.gz" ] \
  && [ ! -e "$G/genes.gtf.gz.part" ] && [ -e "$G/genes.gtf.gz.invalid" ]; } \
  || { echo "FAIL: case truncated_download"; cat "$TEST_TMP/out"; fail=1; }
rm -f "$STUB_STATE/wget_mode"; runh; rc=$?
{ [ $rc -eq 0 ] && grep -q 'gene_id "G1"' "$G/genes.gtf"; } || { echo "FAIL: case rerun_after_truncated"; cat "$TEST_TMP/out"; fail=1; }
rm -f "$G/genes.gtf" "$G/genes.gtf.gz" "$G/genes.gtf.gz.invalid"; echo "$GU fail" > "$STUB_STATE/wget_mode"; runh; rc=$?
{ [ $rc -ne 0 ] && grep -q "download failed for $GU" "$TEST_TMP/out" && [ ! -e "$G/genes.gtf" ] && grep -qx '>1' "$G/genome.fa"; } || { echo "FAIL: case wget_fails"; fail=1; }
rm -f "$STUB_STATE/wget_mode"
mk_helper genome_rel; runh; rc=$?
{ [ $rc -ne 0 ] && grep -q 'must be an absolute path' "$TEST_TMP/out" && [ "$(nwget)" -eq 0 ] && [ ! -e "$TEST_TMP/genome_rel" ]; } || { echo "FAIL: case relative_dir"; fail=1; }
[ -s "$G/sentinel" ] || { echo "FAIL: a file the helper did not create was removed"; fail=1; }
forbidden_ran && fail=1
[ $fail -eq 0 ] && echo "DRY RUN HELPERS PASS" || exit 1
