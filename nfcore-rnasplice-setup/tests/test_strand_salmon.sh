#!/bin/bash
# Usage: test_strand_salmon.sh <skill.md>
# Strandedness from a Salmon subsample (Step 4), in two parts:
#  (1) salmon_strandedness (block after "**Reading the Salmon results.**") on the real lib_format_counts.json and meta_info.json of
#      the real-data test (fixtures realdata_strand_*: ISR 722,724 vs ISF 68 and ISR 745,906 vs ISF 54, both reverse) and on
#      written cases: forward, unstranded, single-end SR/SF/U, the 0.8 and 0.1 boundaries of the RSeQC rule, disagreement with
#      Salmon's expected_format, a low compatible-fragment ratio, missing or repeated counts, other orientations, no counts,
#      missing files;
#  (2) a stub dry run of the helper script (block after "**Strandedness helper script.**"), run from another directory under
#      env -i (test_env.sh): module, singularity (a fake Salmon that writes the real fixture outputs), wget and zcat (system gzip)
#      are logging stubs; every other forbidden command (salmon itself, sbatch, nextflow, conda, ...) stays a failing stub.
# Prints "STRAND SALMON PASS" or exits 1.
set -u
source "$(cd "$(dirname "$0")" && pwd)/test_env.sh" "$@"
SKILL=${1:?usage: test_strand_salmon.sh <skill.md>}
HERE=$(cd "$(dirname "$0")" && pwd); FIX="$HERE/fixtures"
fail=0

# ---- (1) the parser
CODE=$(bash "$HERE/cut_block.sh" "$SKILL" "**Reading the Salmon results.**") || { echo "FAIL: the Salmon results block was not found"; exit 1; }
eval "$CODE"
declare -F salmon_strandedness >/dev/null || { echo "FAIL: salmon_strandedness is not defined by the block"; exit 1; }
P="$TEST_TMP/p"; mkdir -p "$P" || exit 1
real() { # real <run accession>: the real Salmon output of the real-data test
  mkdir -p "$P/$1/aux_info" && cp "$FIX/realdata_strand_$1_lib_format_counts.json" "$P/$1/lib_format_counts.json" \
    && cp "$FIX/realdata_strand_$1_meta_info.json" "$P/$1/aux_info/meta_info.json" || exit 1; }
mk() { # mk <name> <expected_format> <compatible_fragment_ratio> [KEY=count ...]: a lib_format_counts.json (other counts 0)
  local d="$P/$1" ef=$2 cfr=$3 k v kv; shift 3; mkdir -p "$d/aux_info" || exit 1
  { printf '{\n    "read_files": "[ a_1.fq.gz, a_2.fq.gz]",\n    "expected_format": "%s",\n    "compatible_fragment_ratio": %s,\n' "$ef" "$cfr"
    printf '    "num_compatible_fragments": 1000,\n'
    for k in MSF OSF ISF MSR OSR ISR SF SR MU OU IU U; do
      v=0; for kv in "$@"; do [ "${kv%%=*}" = "$k" ] && v=${kv#*=}; done; printf '    "%s": %s,\n' "$k" "$v"
    done
    printf '    "strand_mapping_bias": 0.5\n}\n'; } > "$d/lib_format_counts.json"
  printf '{\n    "salmon_version": "1.10.3",\n    "percent_mapped": 80.5,\n    "start_time": "Fri Oct 02 13:36:58 2026"\n}\n' > "$d/aux_info/meta_info.json"
}
check() { local got; got=$(salmon_strandedness "$P/$1"); [ "${got%% *}" = "$2" ] || { echo "FAIL: case $1: expected $2, got '$got'"; fail=1; }; }
# the real-data test (report: ISR 722,724 vs ISF 68 and ISR 745,906 vs ISF 54; compatible_fragment_ratio 1.0)
real SRR20021259; real SRR20021257
check SRR20021259 reverse; check SRR20021257 reverse
line=$(salmon_strandedness "$P/SRR20021259")
[ "$line" = "reverse $P/SRR20021259 expected_format=ISR forward=0.0001 reverse=0.9999 compatible_fragment_ratio=1.0 percent_mapped=75.3196" ] \
  || { echo "FAIL: case real_line: '$line'"; fail=1; }
# paired-end
mk pe_forward ISF 1.0 ISF=900 ISR=50;            check pe_forward forward
mk pe_unstranded IU 1.0 ISF=500 ISR=500;         check pe_unstranded unstranded
mk pe_iu_only IU 1.0 IU=1000;                    check pe_iu_only unstranded  # edge case: real unstranded data show ISF ~ ISR, IU 0
# live re-run on the unstranded nf-core test data (2026-10-03): ISF 22822, ISR 22215, SF 960, SR 775, IU 0, expected_format IU
mk live_unstranded IU 1.0 ISF=22822 ISR=22215 SF=960 SR=775; check live_unstranded unstranded
# the RSeQC rule's boundaries: >= 0.8 for a strand, |forward - reverse| <= 0.1 for unstranded
mk pe_r_at_080 ISR 1.0 ISR=80 ISF=20;            check pe_r_at_080 reverse
mk pe_f_at_080 ISF 1.0 ISF=80 ISR=20;            check pe_f_at_080 forward
mk pe_r_below ISR 1.0 ISR=7999 ISF=2001;         check pe_r_below unclear
mk pe_diff_010 IU 1.0 ISF=45 ISR=55;             check pe_diff_010 unstranded
mk pe_diff_011 IU 1.0 ISF=445 ISR=555;           check pe_diff_011 unclear
mk pe_iu_dilutes ISR 1.0 ISR=790 ISF=10 IU=200;  check pe_iu_dilutes unclear
# disagreement with Salmon's own detection, compatible-fragment ratio
mk pe_disagree ISF 1.0 ISR=900 ISF=100;          check pe_disagree unclear
mk pe_disagree_u IU 1.0 ISR=950 ISF=50;          check pe_disagree_u unclear
mk pe_unstr_vs_isr ISR 1.0 ISF=500 ISR=500;      check pe_unstr_vs_isr unclear
mk pe_lowcfr ISR 0.5 ISR=900;                    check pe_lowcfr unclear
mk pe_cfr_at_080 ISR 0.8 ISR=900;                check pe_cfr_at_080 reverse
mk pe_cfr_below ISR 0.7999 ISR=900;              check pe_cfr_below unclear
# single-end
mk se_reverse SR 1.0 SR=900 SF=20;               check se_reverse reverse
mk se_forward SF 1.0 SF=900 SR=20;               check se_forward forward
mk se_unstranded U 1.0 SF=480 SR=520;            check se_unstranded unstranded
mk se_paired_keys SR 1.0 ISR=900;                check se_paired_keys unclear
# other orientations, no counts, broken files
mk other_orient MSR 1.0 MSR=900;                 check other_orient unclear
mk no_counts ISR 1.0;                            check no_counts unclear
mk missing_key ISR 1.0 ISR=900; sed -i '/"ISF"/d' "$P/missing_key/lib_format_counts.json"; check missing_key unclear
mk dup_key ISR 1.0 ISR=900; sed -i 's/^    "ISR": 900,$/&\n    "ISR": 5,/' "$P/dup_key/lib_format_counts.json"; check dup_key unclear
mk no_cfr ISR 1.0 ISR=900; sed -i '/"compatible_fragment_ratio"/d' "$P/no_cfr/lib_format_counts.json"; check no_cfr unclear
mk bad_count ISR 1.0 ISR=900; sed -i 's/^    "ISF": 0,$/    "ISF": x,/' "$P/bad_count/lib_format_counts.json"; check bad_count unclear
check no_such_dir unclear
mk no_meta ISR 1.0 ISR=900; rm -f "$P/no_meta/aux_info/meta_info.json"; check no_meta reverse
case "$(salmon_strandedness "$P/no_meta")" in *percent_mapped=NA) ;; *) echo "FAIL: case no_meta: percent_mapped must be NA"; fail=1 ;; esac
case "$(salmon_strandedness "$P/pe_forward")" in *percent_mapped=80.5) ;; *) echo "FAIL: case pe_forward: percent_mapped from meta_info.json"; fail=1 ;; esac

# ---- (2) stub dry run of the helper script
S="$TEST_TMP/scen"; R="$TEST_TMP/project"; ELSE="$TEST_TMP/elsewhere"; IDX="$TEST_TMP/genomes/human/hg38_ens116/index/salmon"
mkdir -p "$S" "$R/data" "$ELSE" "$IDX" "$TEST_TMP/home" "$TEST_TMP/nosing" || exit 1
export STUB_STATE="$TEST_TMP/state" STUB_FIX="$FIX"; mkdir -p "$STUB_STATE" || exit 1
case "$(type -P gzip)" in /usr/bin/gzip|/bin/gzip) SYSGZIP=$(type -P gzip) ;; *) echo "gzip is not the system gzip"; exit 1 ;; esac
cat > "$S/module" <<'EOF'
#!/bin/bash
echo "module $*" >> "$STUB_STATE/calls"
[ -e "$STUB_STATE/module_fail" ] && exit 1
exit 0
EOF
cat > "$S/singularity" <<'EOF'
#!/bin/bash
# fake Salmon in a container: writes the real Salmon outputs of the real-data test (sample SRR20021259) into the -o directory
echo "singularity $*" >> "$STUB_STATE/calls"
[ -e "$STUB_STATE/salmon_fail" ] && exit 1
out=""; while [ $# -gt 0 ]; do [ "$1" = -o ] && out=$2; shift; done
[ -n "$out" ] || exit 2
mkdir -p "$out/aux_info" || exit 2
[ -e "$STUB_STATE/salmon_nojson" ] && exit 0
cp "$STUB_FIX/realdata_strand_SRR20021259_lib_format_counts.json" "$out/lib_format_counts.json" && cp "$STUB_FIX/realdata_strand_SRR20021259_meta_info.json" "$out/aux_info/meta_info.json"
EOF
cat > "$S/wget" <<'EOF'
#!/bin/bash
echo "wget $*" >> "$STUB_STATE/calls"
[ -e "$STUB_STATE/wget_fail" ] && exit 4
out=""; while [ $# -gt 0 ]; do [ "$1" = -O ] && out=$2; shift; done
echo "fake image" > "$out"
EOF
printf '#!/bin/bash\necho "zcat $*" >> "$STUB_STATE/calls"\nexec %s -dc "$@"\n' "$SYSGZIP" > "$S/zcat"
chmod +x "$S/module" "$S/singularity" "$S/wget" "$S/zcat"
for c in $FORBIDDEN_CMDS; do [ "$c" = singularity ] || cp "$STUBS/$c" "$TEST_TMP/nosing/"; done
cp "$S/module" "$S/wget" "$S/zcat" "$TEST_TMP/nosing/"
for c in module singularity wget zcat; do
  [ "$(PATH="$S:$PATH" type -P "$c")" = "$S/$c" ] || { echo "stub $c does not resolve to $S/$c"; exit 1; }
  [ "$(PATH="$S:$PATH" bash -c "type -t $c")" = file ] || { echo "$c is not the stub file in a child shell (an inherited function would bypass it)"; exit 1; }
done
[ -z "$(env | grep '^BASH_FUNC_')" ] || { echo "exported shell functions in the test environment"; exit 1; }
for c in salmon sbatch srun nextflow conda python python3 R Rscript STAR java curl; do
  case "$(PATH="$S:$PATH" type -P "$c")" in "$STUBS/$c") ;; *) echo "forbidden command $c is not a stub on the scenario PATH"; exit 1 ;; esac
done
[ -z "$(PATH="$TEST_TMP/nosing:/usr/bin:/bin" type -P singularity)" ] || { echo "a real singularity is on /usr/bin:/bin; the missing-singularity case cannot be tested"; exit 1; }
for f in info.json pos.bin seq.bin mphf.bin ctable.bin; do echo x > "$IDX/$f"; done
fq() { # fq <file.fq.gz> <reads>
  awk -v n="$2" 'BEGIN { for (i = 1; i <= n; i++) printf "@r%d\nACGTACGT\n+\nIIIIIIII\n", i }' | "$SYSGZIP" -1 > "$1"; }
fq "$R/data/A_1.fq.gz" 1000005; fq "$R/data/A_2.fq.gz" 1000005
fq "$R/data/B_1.fq.gz" 10; fq "$R/data/B_2.fq.gz" 10
fq "$R/data/C.fq.gz" 10
fq "$R/data/D_1.fq.gz" 10; fq "$R/data/D_2.fq.gz" 9
: > "$R/data/E_1.fq.gz"; fq "$R/data/E_2.fq.gz" 10
bash "$HERE/cut_block.sh" "$SKILL" "**Strandedness helper script.**" > "$TEST_TMP/strand.tmpl" || { echo "FAIL: the strandedness helper block was not found"; exit 1; }
[ "$(grep -cx 'run_sample "{SAMPLE}" "{FASTQ_1}" "{FASTQ_2}"' "$TEST_TMP/strand.tmpl")" -eq 1 ] || { echo "FAIL: the helper must hold one template run_sample line"; exit 1; }
OUT="$R/strandedness_salmon"; SIF="$TEST_TMP/home/.singularity/cache/depot.galaxyproject.org-singularity-salmon-1.10.3--h6dccd9a_2.img"
URL=https://depot.galaxyproject.org/singularity/salmon:1.10.3--h6dccd9a_2
mk_helper() { # mk_helper <run_sample lines...>: fills the placeholders as the wizard would
  local l lines=""; for l in "$@"; do lines+="$l"$'\n'; done
  awk -v lines="$lines" '$0 == "run_sample \"{SAMPLE}\" \"{FASTQ_1}\" \"{FASTQ_2}\"" { printf "%s", lines; next } { print }' "$TEST_TMP/strand.tmpl" \
    | sed -e "s#{USER_EMAIL}#test@example.org#g; s#{WD_NAME}#project#g; s#{CWD}#$R#g; s#{STRAND_INDEX}#$IDX#g" > "$TEST_TMP/helper.sh"
  ! grep -nE '(^|[^$])\{[A-Z_0-9]+\}' "$TEST_TMP/helper.sh" || { echo "FAIL: placeholder left in the helper"; exit 1; }
}
runh() { # runh <name> <expected exit 0|1> [module_fail salmon_fail salmon_nojson wget_fail no_singularity no_sif]
  local name=$1 erc=$2 p="$S:$PATH" w rc; shift 2
  local keep=0; for w in "$@"; do [ "$w" = keep_outputs ] && keep=1; done
  # the helper refuses an existing sample directory (it never deletes): the test clears them between cases
  if [ $keep -eq 0 ]; then for w in "$OUT"/*/; do [ -d "$w" ] && rm -rf -- "$w"; done; fi
  rm -f "$STUB_STATE/calls" "$STUB_STATE/module_fail" "$STUB_STATE/salmon_fail" "$STUB_STATE/salmon_nojson" "$STUB_STATE/wget_fail"
  mkdir -p "$(dirname "$SIF")" && echo "cached image" > "$SIF" || exit 1
  for w in "$@"; do case "$w" in
    no_singularity) p="$TEST_TMP/nosing:/usr/bin:/bin" ;;
    no_sif) mv "$SIF" "$TEST_TMP/held_sif" || exit 1 ;;
    keep_outputs) ;;
    *) : > "$STUB_STATE/$w" ;;
  esac; done
  ( cd "$ELSE" && HOME="$TEST_TMP/home" PATH="$p" bash "$TEST_TMP/helper.sh" ) > "$TEST_TMP/out" 2>&1; rc=$?
  [ "$rc" -eq "$erc" ] || { echo "FAIL: case $name: exit $rc, expected $erc"; sed 's/^/    /' "$TEST_TMP/out"; fail=1; }
  rm -f "$TEST_TMP/held_sif"
}
lines_of() { "$SYSGZIP" -dc "$1" | wc -l; }
ncalls() { local n; n=$(grep -c "^$1 " "$STUB_STATE/calls" 2>/dev/null); echo "${n:-0}"; }
mkdir -p "$OUT" && echo keep > "$OUT/sentinel" || exit 1
# paired-end, relative and absolute paths, image in the cache; then the parser reads the outputs
mk_helper 'run_sample "A" "data/A_1.fq.gz" "data/A_2.fq.gz"' "run_sample \"B\" \"$R/data/B_1.fq.gz\" \"$R/data/B_2.fq.gz\""
runh paired 0
{ grep -q '^ALL DONE: results in ' "$TEST_TMP/out" && [ "$(lines_of "$OUT/A_sub_1.fq.gz")" -eq 4000000 ] && [ "$(lines_of "$OUT/A_sub_2.fq.gz")" -eq 4000000 ] \
  && [ "$(lines_of "$OUT/B_sub_1.fq.gz")" -eq 40 ] && grep -qx 'DONE A (1000000 reads)' "$TEST_TMP/out" && grep -qx 'DONE B (10 reads)' "$TEST_TMP/out"; } \
  || { echo "FAIL: case paired: subsamples or messages"; sed 's/^/    /' "$TEST_TMP/out"; fail=1; }
for s in A B; do
  grep -qxF "singularity exec -B $R,$IDX $SIF salmon quant -i $IDX -l A -1 $OUT/${s}_sub_1.fq.gz -2 $OUT/${s}_sub_2.fq.gz -p 8 -o $OUT/$s" "$STUB_STATE/calls" \
    || { echo "FAIL: case paired: Salmon call for $s:"; cat "$STUB_STATE/calls"; fail=1; }
  [ "$(salmon_strandedness "$OUT/$s" | cut -d' ' -f1)" = reverse ] || { echo "FAIL: case paired: the parser on the helper output of $s"; fail=1; }
done
{ [ "$(ncalls singularity)" -eq 2 ] && [ "$(ncalls wget)" -eq 0 ] && grep -qx 'module add singularity/3.10.4' "$STUB_STATE/calls"; } \
  || { echo "FAIL: case paired: calls"; cat "$STUB_STATE/calls"; fail=1; }
awk '/^module add singularity/ && !m {m = NR} /^singularity / && !s {s = NR} END {exit !(m && s && m < s)}' "$STUB_STATE/calls" || { echo "FAIL: case paired: module add singularity must come first"; fail=1; }
[ ! -e "$ELSE/strandedness_salmon" ] || { echo "FAIL: case paired: output written in the submit directory, not in {CWD}"; fail=1; }
# a rerun with an existing sample directory is refused (old results are never read as new ones)
mk_helper 'run_sample "B" "data/B_1.fq.gz" "data/B_2.fq.gz"'
runh stale_output 1 keep_outputs
{ grep -q "B: $OUT/B exists from an earlier run" "$TEST_TMP/out" && [ "$(ncalls singularity)" -eq 0 ] && [ "$(ncalls zcat)" -eq 0 ]; } || { echo "FAIL: case stale_output"; sed 's/^/    /' "$TEST_TMP/out"; fail=1; }
# single-end: salmon -r, no -1/-2
mk_helper 'run_sample "C" "data/C.fq.gz"'
runh single 0
grep -qxF "singularity exec -B $R,$IDX $SIF salmon quant -i $IDX -l A -r $OUT/C_sub_1.fq.gz -p 8 -o $OUT/C" "$STUB_STATE/calls" \
  || { echo "FAIL: case single: Salmon call:"; cat "$STUB_STATE/calls"; fail=1; }
# the image is downloaded (resumable) only when it is not in the cache
mk_helper 'run_sample "B" "data/B_1.fq.gz" "data/B_2.fq.gz"'
runh sif_download 0 no_sif
{ grep -qxF "wget -c -O $SIF.part $URL" "$STUB_STATE/calls" && [ -s "$SIF" ] && [ "$(ncalls singularity)" -eq 1 ]; } || { echo "FAIL: case sif_download"; cat "$STUB_STATE/calls"; fail=1; }
runh wget_fails 1 no_sif wget_fail
{ grep -q 'could not download the Salmon container' "$TEST_TMP/out" && [ "$(ncalls singularity)" -eq 0 ]; } || { echo "FAIL: case wget_fails"; fail=1; }
runh module_fails 1 module_fail
[ "$(ncalls singularity)" -eq 0 ] && [ "$(ncalls zcat)" -eq 0 ] || { echo "FAIL: case module_fails: the helper went on"; fail=1; }
runh singularity_missing 1 no_singularity
grep -q 'singularity is not on PATH' "$TEST_TMP/out" || { echo "FAIL: case singularity_missing: message"; fail=1; }
runh salmon_fails 1 salmon_fail
grep -q 'B: Salmon failed' "$TEST_TMP/out" || { echo "FAIL: case salmon_fails: message"; fail=1; }
mk_helper 'run_sample "N" "data/B_1.fq.gz" "data/B_2.fq.gz"'
runh no_json 1 salmon_nojson
grep -q 'N: Salmon wrote no lib_format_counts.json' "$TEST_TMP/out" || { echo "FAIL: case no_json: message"; fail=1; }
# input problems: one failing sample does not hide the others, and the job ends with exit 1
mk_helper 'run_sample "X" "data/missing_1.fq.gz" "data/missing_2.fq.gz"' 'run_sample "B" "data/B_1.fq.gz" "data/B_2.fq.gz"'
runh missing_fastq 1
{ grep -q 'X: data/missing_1.fq.gz is missing or empty' "$TEST_TMP/out" && grep -qx 'DONE B (10 reads)' "$TEST_TMP/out" && grep -q 'at least one sample failed' "$TEST_TMP/out"; } \
  || { echo "FAIL: case missing_fastq"; sed 's/^/    /' "$TEST_TMP/out"; fail=1; }
mk_helper 'run_sample "D" "data/D_1.fq.gz" "data/D_2.fq.gz"'
runh r1_r2_differ 1
{ grep -q 'D: the R1 and R2 subsamples differ (40 and 36 lines)' "$TEST_TMP/out" && [ "$(ncalls singularity)" -eq 0 ]; } || { echo "FAIL: case r1_r2_differ"; sed 's/^/    /' "$TEST_TMP/out"; fail=1; }
mk_helper 'run_sample "E" "data/E_1.fq.gz" "data/E_2.fq.gz"'
runh empty_fastq 1
grep -q 'E: data/E_1.fq.gz is missing or empty' "$TEST_TMP/out" || { echo "FAIL: case empty_fastq"; sed 's/^/    /' "$TEST_TMP/out"; fail=1; }
# an incomplete index stops the job before anything else
rm -f "$IDX/pos.bin"
mk_helper 'run_sample "B" "data/B_1.fq.gz" "data/B_2.fq.gz"'
runh bad_index 1
{ grep -q 'pos.bin is missing or empty: not a usable Salmon index' "$TEST_TMP/out" && [ ! -s "$STUB_STATE/calls" ]; } || { echo "FAIL: case bad_index"; fail=1; }
[ "$(cat "$OUT/sentinel")" = keep ] || { echo "FAIL: a file the helper did not create was changed or removed"; fail=1; }
forbidden_ran && fail=1
[ $fail -eq 0 ] && echo "STRAND SALMON PASS" || exit 1
