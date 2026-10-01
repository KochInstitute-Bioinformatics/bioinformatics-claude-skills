#!/bin/bash
# Usage: dry_run_chain.sh <skill.md> <scratch dir>  (deleted and re-created only under /net/bmc-lab3/data/bcc/ase_scratch/;
# anywhere else the run uses a new mktemp -d directory inside it and removes only that)
# Cuts submit_chain.sh and wait_chain.sh out of Step 15 (cut_block.sh), edits them as the wizard would (placeholders; the four
# mouse-helper lines deleted unless a scenario keeps them), and runs them against stub sbatch / squeue / scancel / sleep. Each
# submission scenario compares the submitted jobs and their afterok dependencies (as script names) with the expected list, or
# the stop message; each wait scenario compares the per-job report, the squeue job list and the scancel calls.
# Prints "DRY RUN PASS" or exits 1.
#
# Safety: only bash and stubs run. The script re-executes itself under `env -i` with PATH = <stub dir>:/usr/bin:/bin, so no
# exported shell function (for example the Lmod `module` function) can shadow a stub, and it checks that every stub resolves
# to the stub file before any scenario runs. Commands that must never run here (conda, python, phaser, R, module, singularity,
# STAR) are stubs that fail and leave a trace, which fails the run.
set -u
if [ "${DRY_RUN_CLEAN_ENV:-}" != 1 ]; then
  exec env -i DRY_RUN_CLEAN_ENV=1 HOME="${HOME:-/nonexistent}" PATH=/usr/bin:/bin /bin/bash --noprofile --norc "$0" "$@"
fi
SKILL=${1:?usage: dry_run_chain.sh <skill.md> <scratch dir>}; W=${2:?usage: dry_run_chain.sh <skill.md> <scratch dir>}
case "$SKILL" in /*) ;; *) SKILL="$PWD/$SKILL" ;; esac
case "$W" in /*) ;; *) W="$PWD/$W" ;; esac
HERE=$(cd "$(dirname "$0")" && pwd); CUT="$HERE/../synthetic/cut_block.sh"; fail=0
COMMIT=aa1f8ec5fe1cc676e37cfa6f6a0bce6b09070301
# line 2 of install_ok.txt written by the real installation of the pinned package list (cluster run of setup_phaser_env.sh);
# submit_chain.sh must compute the same hash from the script, otherwise it would reinstall on every run
ENV_SHA_INSTALLED=ef6e8597e66a4e6c23c27fee56181c0ae37810c1e94370a375de44974f6507a2
# the scratch argument is deleted only under the shared scratch area; anywhere else the test works in a new mktemp -d directory
# inside it and removes only that directory at the end (a project directory given by mistake is never deleted)
case "$W" in */../*|*/..|*/./*) echo "FAIL: scratch dir '$W' must be a plain absolute path"; exit 1 ;; esac
W=${W%/}; OWN_TMP=""
case "$W" in
  /net/bmc-lab3/data/bcc/ase_scratch/?*) rm -rf "$W" && mkdir -p "$W" || exit 1 ;;
  ""|/) echo "FAIL: refusing scratch dir '/'"; exit 1 ;;
  *) mkdir -p "$W" && OWN_TMP=$(mktemp -d "$W/dry_run_chain.XXXXXX") || exit 1; W=$OWN_TMP ;;
esac
mkdir -p "$W/bin" || exit 1
B="$W/bin"
cat > "$B/sbatch" <<'EOF'
#!/bin/bash
# stub: prints the next job id and logs "<id> <dependency or -> <script>"; refuses a missing script and an empty id in a
# dependency; prints no id for the script named in $STUB_DIR/empty_for
n=$(( $(cat "$STUB_DIR/n" 2>/dev/null || echo 100) + 1 )); echo "$n" > "$STUB_DIR/n"
dep=-; script=""
for a in "$@"; do case "$a" in --dependency=*) dep=${a#--dependency=} ;; -p|bcc|--parsable) ;; *) script=$a ;; esac; done
case "$dep" in *afterok:|*::*|*:) echo "stub sbatch: empty id in dependency '$dep'" >&2; exit 1 ;; esac
[ -s "$script" ] || { echo "stub sbatch: no such script $script" >&2; exit 1; }
[ "$(cat "$STUB_DIR/empty_for" 2>/dev/null)" = "$(basename "$script")" ] && exit 0
echo "$n $dep $(basename "$script")" >> "$STUB_DIR/calls"; echo "$n"
EOF
cat > "$B/squeue" <<'EOF'
#!/bin/bash
# stub: logs the -j list; call k prints $STUB_DIR/q.k (nothing when absent, that is: the queue is empty)
k=$(( $(cat "$STUB_DIR/sq_n" 2>/dev/null || echo 0) + 1 )); echo "$k" > "$STUB_DIR/sq_n"
prev=""; for a in "$@"; do [ "$prev" = -j ] && echo "$a" >> "$STUB_DIR/squeue_j"; prev=$a; done
[ -f "$STUB_DIR/q.$k" ] && cat "$STUB_DIR/q.$k"; exit 0
EOF
printf '#!/bin/bash\necho "$*" >> "$STUB_DIR/scancel"\n' > "$B/scancel"
printf '#!/bin/bash\necho "$*" >> "$STUB_DIR/sleep"\n' > "$B/sleep"
FORBIDDEN="conda python python3 phaser phaser.py Rscript R module singularity STAR pip sacct"
for c in $FORBIDDEN; do printf '#!/bin/bash\necho "FORBIDDEN %s $*" >> "%s/forbidden"; exit 99\n' "$c" "$W" > "$B/$c"; done
chmod +x "$B/"*
# every stub must resolve to its file before anything runs (type -t: file, not a function, alias or builtin)
for c in sbatch squeue scancel sleep $FORBIDDEN; do
  t=$(PATH="$B:/usr/bin:/bin" type -t "$c"); p=$(PATH="$B:/usr/bin:/bin" type -P "$c")
  [ "$t" = file ] && [ "$p" = "$B/$c" ] || { echo "FAIL: stub $c does not resolve to $B/$c (type: ${t:-none}, path: ${p:-none}); nothing run"; exit 1; }
done
echo "stubs resolve: sbatch squeue scancel sleep and $(echo $FORBIDDEN | wc -w) forbidden-command stubs"
SUB=$(bash "$CUT" "$SKILL" "### Submission order (one block") || { echo "FAIL: no submit_chain.sh block"; exit 1; }
WAIT=$(bash "$CUT" "$SKILL" "**Waiting.**") || { echo "FAIL: no wait_chain.sh block"; exit 1; }
SETUP=$(bash "$CUT" "$SKILL" '### `setup_phaser_env.sh`') || { echo "FAIL: no setup_phaser_env.sh block"; exit 1; }
PHC=$(bash "$CUT" "$SKILL" '### `phaser_count.sh`') || { echo "FAIL: no phaser_count.sh block"; exit 1; }
run_stub() {   # $1 dir, $2 script: runs it with the stubs only (clean environment, stub dir first on PATH)
  (cd "$1" && env -i HOME="$HOME" STUB_DIR="$1" PATH="$B:/usr/bin:/bin" /bin/bash --noprofile --norc "$2" 2>&1)
}
named() {   # calls file -> "script afterok:<names>" lines
  awk '{ id[$1] = $3; d = $2
         if (d != "-") { sub(/^afterok:/, "", d); n = split(d, x, ":"); d = "afterok:"; for (i = 1; i <= n; i++) d = d (i > 1 ? ":" : "") id[x[i]] }
         print $3, d }' "$1"; }
ALL="prep_genotypes.sh align_wasp_count.sh prep_f1_reference.sh align_count_f1.sh run_01_import_qc.sh run_02_imbalance.sh
     run_03_reciprocal.sh run_04_differential.sh setup_phaser_env.sh phaser_count.sh run_05_phaser.sh
     extract_mgp_parental_vcf.sh concat_mgp_parental_vcf.sh"
# scenario: $1 name, $2 MODE, $3 ANALYSES, $4 install_ok.txt: yes (both lines of the real installation) | stale (line 2 differs)
# | no, $5 expected lines or "DIE:<text>", $6 scripts left out; optional env: HOMEVAL (PHASER_HOME), KEEP_HELPER=1,
# EMPTY_FOR (script that gets no job id), PRE_IDS=1 (chain_job_ids.tsv exists), BAD_COMMIT (script with another pinned commit)
scenario() {
  local D="$W/$1" got out rc s home
  mkdir -p "$D/res/scripts" "$D/res/logs" "$D/home"
  for s in $ALL; do
    case " $6 " in *" $s "*) continue ;; esac
    case $s in setup_phaser_env.sh) printf '%s\n' "$SETUP" > "$D/res/scripts/$s" ;;
               phaser_count.sh) printf '%s\n' "$PHC" > "$D/res/scripts/$s" ;;
               *) echo '#!/bin/bash' > "$D/res/scripts/$s" ;; esac
  done
  home=${HOMEVAL:-$D/home}
  for s in setup_phaser_env.sh phaser_count.sh; do [ -f "$D/res/scripts/$s" ] && sed -i "s|{PHASER_HOME}|$home|g" "$D/res/scripts/$s"; done
  [ -n "${BAD_COMMIT:-}" ] && sed -i "s/$COMMIT/${COMMIT_REPL:-0000000000000000000000000000000000000000}/" "$D/res/scripts/$BAD_COMMIT"
  [ -n "${OTHER_HOME_IN:-}" ] && sed -i "s|PHASER_HOME=\"$home\"|PHASER_HOME=\"$D/other_home\"|" "$D/res/scripts/$OTHER_HOME_IN"
  case "$4" in yes)   printf 'commit %s\nenv_sha256 %s\npython 3.14.7\n' "$COMMIT" "$ENV_SHA_INSTALLED" > "$D/home/install_ok.txt" ;;
               stale) printf 'commit %s\nenv_sha256 %s\npython 3.14.7\n' "$COMMIT" "$(printf '%064d' 0)" > "$D/home/install_ok.txt" ;; esac
  [ -n "${EMPTY_FOR:-}" ] && echo "$EMPTY_FOR" > "$D/empty_for"
  [ -n "${PRE_IDS:-}" ] && printf 'old.sh\t1\n' > "$D/res/logs/chain_job_ids.tsv"
  if [ -n "${KEEP_HELPER:-}" ]; then printf '%s\n' "$SUB"
  else printf '%s\n' "$SUB" | awk '/^# Without the helper, delete the next four lines/ {print; skip = 4; next} skip > 0 {skip--; next} {print}'
  fi | sed -e "s|{RESULTS_DIR}|$D/res|g" -e "s|{CWD}|$D|g" -e "s|{MODE}|$2|g" -e "s|{ANALYSES}|$3|g" -e "s|{PHASER_HOME}|$home|g" > "$D/full_chain.sh"
  case "${RECIPE:-}" in
    "") cp "$D/full_chain.sh" "$D/submit_chain.sh" ;;
    naive) grep -v "^  PH=" "$D/full_chain.sh" > "$D/submit_chain.sh" ;;   # the copy-and-delete hazard the recipes avoid
    *) { awk '{print} /^DEP=""; R3="not selected"; R4="not selected"; R5="not selected"; PS=""/ {exit}' "$D/full_chain.sh"
         bash "$CUT" "$SKILL" "**Recipe $RECIPE " || echo "echo MISSING RECIPE $RECIPE; exit 3"; tail -n 2 "$D/full_chain.sh"; } > "$D/submit_chain.sh" ;;
  esac
  grep -q '{[A-Z_]*}' "$D/submit_chain.sh" && { echo "FAIL: $1: placeholder left in submit_chain.sh"; fail=1; return; }
  out=$(run_stub "$D" "$D/submit_chain.sh"); rc=$?
  case "${RECIPE:-}" in
    naive) printf '%s' "$out" | grep -q 'PH: unbound variable' || { echo "FAIL: $1: expected the unbound-variable hazard: $out"; fail=1; } ;;
    ?*) ! printf '%s' "$out" | grep -q 'unbound variable' || { echo "FAIL: $1: a recipe line uses an unset variable: $out"; fail=1; } ;;
  esac
  case "$5" in
    DIE:*) if [ $rc -ne 0 ] && printf '%s' "$out" | grep -qF -- "${5#DIE:}" && [ ! -s "$D/calls" ]; then echo "ok   $1: stops ('${5#DIE:}'), nothing submitted"
           else echo "FAIL: $1: expected a stop with '${5#DIE:}' and no submission (rc $rc): $out"; fail=1; fi
           if [ -n "${PRE_IDS:-}" ] && [ "$(cat "$D/res/logs/chain_job_ids.tsv")" != "$(printf 'old.sh\t1')" ]; then
             echo "FAIL: $1: the existing chain_job_ids.tsv was changed"; fail=1; fi ;;
    STOP:*) got=$([ -s "$D/calls" ] && named "$D/calls")
           if [ $rc -ne 0 ] && printf '%s' "$out" | grep -qF -- "${5#STOP:}" && [ "$got" = "$EXPECT_PART" ] &&
              [ "$(cut -f1 "$D/res/logs/chain_job_ids.tsv" | paste -sd' ' -)" = "$(printf '%s\n' "$EXPECT_PART" | cut -d' ' -f1 | paste -sd' ' -)" ]
           then echo "ok   $1: stops ('${5#STOP:}') after $(printf '%s\n' "$got" | wc -l) jobs, all of them recorded in chain_job_ids.tsv"
           else echo "FAIL: $1 (rc $rc)"; printf '  expected:\n%s\n  got:\n%s\n  output: %s\n' "$EXPECT_PART" "$got" "$out"; fail=1; fi ;;
    *) got=$([ -s "$D/calls" ] && named "$D/calls")
       if [ $rc -eq 0 ] && [ "$got" = "$5" ] && [ "$(cut -f2 "$D/res/logs/chain_job_ids.tsv" | paste -sd' ' -)" = "$(cut -d' ' -f1 "$D/calls" | paste -sd' ' -)" ]
       then echo "ok   $1: $(printf '%s\n' "$got" | wc -l) jobs with the expected dependencies; $(printf '%s\n' "$out" | tail -1)"
       else echo "FAIL: $1 (rc $rc)"; printf '  expected:\n%s\n  got:\n%s\n  output: %s\n' "$5" "$got" "$out"; fail=1; fi ;;
  esac
}
BASE_OB='prep_genotypes.sh -
align_wasp_count.sh afterok:prep_genotypes.sh
run_01_import_qc.sh afterok:align_wasp_count.sh
run_02_imbalance.sh afterok:run_01_import_qc.sh'
PH_INST='phaser_count.sh afterok:align_wasp_count.sh
run_05_phaser.sh afterok:phaser_count.sh:run_02_imbalance.sh'
PH_NEW='setup_phaser_env.sh -
phaser_count.sh afterok:align_wasp_count.sh:setup_phaser_env.sh
run_05_phaser.sh afterok:phaser_count.sh:run_02_imbalance.sh'
BASE_F1='prep_f1_reference.sh -
align_count_f1.sh afterok:prep_f1_reference.sh
run_01_import_qc.sh afterok:align_count_f1.sh
run_02_imbalance.sh afterok:run_01_import_qc.sh'
echo "--- submit_chain.sh"
scenario ob_phaser_installed outbred "per-sample phaser" yes "$BASE_OB
$PH_INST" ""
scenario ob_phaser_install outbred "per-sample phaser" no "$BASE_OB
$PH_NEW" ""
scenario ob_phaser_stale_env outbred "per-sample phaser" stale "$BASE_OB
$PH_NEW" ""
scenario ob_diff_phaser outbred "per-sample differential phaser" yes "$BASE_OB
run_04_differential.sh afterok:run_01_import_qc.sh
$PH_INST" ""
scenario ob_diff_phaser_install outbred "per-sample differential phaser" no "$BASE_OB
run_04_differential.sh afterok:run_01_import_qc.sh
$PH_NEW" ""
scenario ob_per_sample_only outbred "per-sample" no "$BASE_OB" "setup_phaser_env.sh phaser_count.sh run_05_phaser.sh"
HOMEVAL=none scenario ob_no_phaser_home_none outbred "per-sample differential" no "$BASE_OB
run_04_differential.sh afterok:run_01_import_qc.sh" "setup_phaser_env.sh phaser_count.sh run_05_phaser.sh"
scenario f1_all f1 "per-sample reciprocal differential" no "$BASE_F1
run_03_reciprocal.sh afterok:run_01_import_qc.sh
run_04_differential.sh afterok:run_01_import_qc.sh" ""
KEEP_HELPER=1 scenario f1_mouse_helper f1 "per-sample" no "extract_mgp_parental_vcf.sh -
concat_mgp_parental_vcf.sh afterok:extract_mgp_parental_vcf.sh
prep_f1_reference.sh afterok:concat_mgp_parental_vcf.sh
align_count_f1.sh afterok:prep_f1_reference.sh
run_01_import_qc.sh afterok:align_count_f1.sh
run_02_imbalance.sh afterok:run_01_import_qc.sh" ""
scenario f1_phaser f1 "per-sample phaser" yes "DIE:phaser analysis is outbred only" ""
scenario ob_reciprocal outbred "per-sample reciprocal" yes "DIE:reciprocal analysis is F1 only" ""
scenario menu_numbers outbred "per-sample 4" yes "DIE:ANALYSES contains '4'" ""
scenario unknown_word outbred "per-sample phASER" yes "DIE:ANALYSES contains 'phASER'" ""
scenario empty_analyses outbred "" yes "DIE:ANALYSES is empty" ""
scenario missing_rmd05 outbred "per-sample phaser" yes "DIE:MISSING in" "run_05_phaser.sh"
scenario missing_phaser_count outbred "per-sample phaser" yes "DIE:MISSING in" "phaser_count.sh"
scenario missing_setup outbred "per-sample phaser" yes "DIE:MISSING in" "setup_phaser_env.sh"
HOMEVAL=none scenario phaser_home_none outbred "per-sample phaser" no "DIE:PHASER_HOME must be the absolute path" ""
HOMEVAL=rel/home scenario phaser_home_relative outbred "per-sample phaser" no "DIE:PHASER_HOME must be the absolute path" ""
BAD_COMMIT=phaser_count.sh scenario commit_mismatch_count outbred "per-sample phaser" yes "DIE:phaser_count.sh does not pin phASER commit" ""
BAD_COMMIT=setup_phaser_env.sh scenario commit_mismatch_setup outbred "per-sample phaser" yes "DIE:setup_phaser_env.sh does not pin phASER commit" ""
PRE_IDS=1 scenario existing_ids outbred "per-sample phaser" yes "DIE:a chain has been submitted from this results directory" ""
EXPECT_PART="$BASE_OB" EMPTY_FOR=phaser_count.sh scenario empty_id_phaser outbred "per-sample phaser" yes "STOP:sbatch returned no job id for phaser_count.sh" ""
EXPECT_PART="$BASE_OB
setup_phaser_env.sh -
phaser_count.sh afterok:align_wasp_count.sh:setup_phaser_env.sh" EMPTY_FOR=run_05_phaser.sh scenario empty_id_rmd05 outbred "per-sample phaser" no "STOP:sbatch returned no job id for run_05_phaser.sh" ""
BAD_COMMIT=phaser_count.sh COMMIT_REPL="${COMMIT}0" scenario commit_with_suffix outbred "per-sample phaser" yes "DIE:phaser_count.sh does not pin phASER commit" ""
OTHER_HOME_IN=phaser_count.sh scenario home_mismatch_count outbred "per-sample phaser" yes "DIE:phaser_count.sh uses another PHASER_HOME" ""
OTHER_HOME_IN=setup_phaser_env.sh scenario home_mismatch_setup outbred "per-sample phaser" yes "DIE:setup_phaser_env.sh uses another PHASER_HOME" ""
echo "--- re-submission recipes of Step 15 (part 1 of submit_chain.sh + one recipe + the two echo lines)"
RECIPE=A scenario recipe_A_rmd05_alone outbred "per-sample phaser" yes "run_05_phaser.sh -" ""
RECIPE=B scenario recipe_B_rmd02_then_rmd05 outbred "per-sample phaser" yes "run_02_imbalance.sh -
run_05_phaser.sh afterok:run_02_imbalance.sh" ""
RECIPE=C scenario recipe_C_phaser_installed outbred "per-sample phaser" yes "phaser_count.sh -
run_05_phaser.sh afterok:phaser_count.sh" ""
RECIPE=C scenario recipe_C_phaser_with_install outbred "per-sample phaser" no "setup_phaser_env.sh -
phaser_count.sh afterok:setup_phaser_env.sh
run_05_phaser.sh afterok:phaser_count.sh" ""
RECIPE=D scenario recipe_D_rmd01_on_diff outbred "per-sample differential phaser" yes "run_01_import_qc.sh -
run_02_imbalance.sh afterok:run_01_import_qc.sh
run_04_differential.sh afterok:run_01_import_qc.sh
phaser_count.sh -
run_05_phaser.sh afterok:phaser_count.sh:run_02_imbalance.sh" ""
RECIPE=D scenario recipe_D_rmd01_on_install outbred "per-sample phaser" no "run_01_import_qc.sh -
run_02_imbalance.sh afterok:run_01_import_qc.sh
setup_phaser_env.sh -
phaser_count.sh afterok:setup_phaser_env.sh
run_05_phaser.sh afterok:phaser_count.sh:run_02_imbalance.sh" ""
# recipe E (a sample dropped): no phaser_count.sh and no installation job, even when the installation test says "install"
RECIPE=E scenario recipe_E_drop_sample_diff outbred "per-sample differential phaser" yes "run_01_import_qc.sh -
run_02_imbalance.sh afterok:run_01_import_qc.sh
run_04_differential.sh afterok:run_01_import_qc.sh
run_05_phaser.sh afterok:run_02_imbalance.sh" ""
RECIPE=E scenario recipe_E_drop_sample_not_installed outbred "per-sample phaser" no "run_01_import_qc.sh -
run_02_imbalance.sh afterok:run_01_import_qc.sh
run_05_phaser.sh afterok:run_02_imbalance.sh" ""
EXPECT_PART="$BASE_OB" RECIPE=naive scenario naive_copy_without_PH outbred "per-sample phaser" yes "STOP:sbatch returned no job id for run_05_phaser.sh" ""

# wait scenarios: $1 name, $2 job list (script<TAB>id lines), $3 result files written after the job list, $4 logs ("file|line"
# entries separated by newlines), $5 squeue outputs per call (blocks separated by a line "--"), $6 expected report lines (each
# must appear exactly), $7 expected scancel calls (one line per call, or empty), $8 expected exit status
wait_case() {
  local D="$W/$1" out rc k blk l f
  mkdir -p "$D/res/logs"
  printf '%b' "$2" > "$D/res/logs/chain_job_ids.tsv"; touch -d '2 minutes ago' "$D/res/logs/chain_job_ids.tsv"
  for f in $3; do touch "$D/res/$f"; done
  while IFS='|' read -r f l; do [ -n "$f" ] && printf '%s\n' "$l" >> "$D/res/logs/$f"; done <<< "$4"
  k=1; blk=""
  while IFS= read -r l; do
    if [ "$l" = "--" ]; then printf '%s' "$blk" > "$D/q.$k"; k=$((k + 1)); blk=""; else blk="$blk$l"$'\n'; fi
  done <<< "$5"
  [ -n "$blk" ] && printf '%s' "$blk" > "$D/q.$k"
  printf '%s\n' "$WAIT" | sed -e "s|{RESULTS_DIR}|$D/res|g" > "$D/wait_chain.sh"
  out=$(run_stub "$D" "$D/wait_chain.sh"); rc=$?
  local ok=1 want_j
  while IFS= read -r l; do [ -z "$l" ] || printf '%s\n' "$out" | grep -qxF -- "$l" || { ok=0; echo "  missing report line: $l"; }; done <<< "$6"
  [ "$(cat "$D/scancel" 2>/dev/null)" = "$7" ] || { ok=0; echo "  scancel calls: '$(cat "$D/scancel" 2>/dev/null)', expected '$7'"; }
  want_j=$(printf '%b' "$2" | cut -f2 | paste -sd, -)
  [ -z "$(sort -u "$D/squeue_j" | grep -vxF -- "$want_j")" ] && [ -s "$D/squeue_j" ] || { ok=0; echo "  squeue was not asked for every job id ($want_j): $(sort -u "$D/squeue_j" | paste -sd' ' -)"; }
  [ $rc -eq "$8" ] || { ok=0; echo "  exit status $rc, expected $8"; }
  if [ $ok = 1 ]; then echo "ok   $1: $(printf '%s\n' "$out" | grep -c .) report lines as expected, scancel '${7:-none}', exit $rc"
  else echo "FAIL: $1"; printf '%s\n' "$out" | sed 's/^/  | /'; fail=1; fi
}
echo "--- wait_chain.sh"
J7='prep_genotypes.sh\t101\nalign_wasp_count.sh\t102\nrun_01_import_qc.sh\t103\nrun_02_imbalance.sh\t104\nrun_04_differential.sh\t105\nphaser_count.sh\t106\nrun_05_phaser.sh\t107\n'
OKRES="ase_checkpoint.rds summary_numbers.tsv summary_numbers_differential.tsv"
wait_case wait_rmd05_ok "$J7" "$OKRES summary_numbers_phaser.tsv" "" "" "run_05_phaser.sh 107 ok
phaser_count.sh 106 finished
run_04_differential.sh 105 ok" "" 0
wait_case wait_rmd05_failed "$J7" "$OKRES" "run_05_phaser_107.out|Error: phASER gene counts missing for sample(s) s2" "" \
  "run_05_phaser.sh 107 FAILED: Error: phASER gene counts missing for sample(s) s2
phaser_count.sh 106 finished" "" 1
wait_case rmd02_blocks_rmd05 "$J7" "ase_checkpoint.rds summary_numbers_differential.tsv" \
  "run_02_imbalance_104.out|Error in bb_fit(x): no sites" "105 None
107 DependencyNeverSatisfied
--
107 DependencyNeverSatisfied" "run_02_imbalance.sh 104 FAILED: Error in bb_fit(x): no sites
run_04_differential.sh 105 ok
phaser_count.sh 106 finished
run_05_phaser.sh 107 CANCELLED: never satisfiable, upstream failure: run_02_imbalance.sh 104" "107" 1
wait_case phaser_array_fails "$J7" "$OKRES" "phaser_count_106_1.out|Sample s1: 50 haplotype blocks
phaser_count_106_2.out|Traceback (most recent call last):
phaser_count_106_2.out|  File \"phaser.py\", line 812, in <module>
phaser_count_106_2.out|NameError: name 'args' is not defined
phaser_count_106_2.out|ERROR: sample s2: phaser.py failed (log res/logs/phaser_s2.log)" "107 DependencyNeverSatisfied" \
  "run_02_imbalance.sh 104 ok
run_04_differential.sh 105 ok
phaser_count.sh 106 finished; first error line of its log: NameError: name 'args' is not defined
run_05_phaser.sh 107 CANCELLED: never satisfiable, upstream failure: phaser_count.sh 106" "107" 1
wait_case phaser_and_rmd02_fail "$J7" "ase_checkpoint.rds summary_numbers_differential.tsv" "run_02_imbalance_104.out|Error: stop
phaser_count_106_1.out|FATAL ERROR - no reads found" "107 DependencyNeverSatisfied" \
  "run_05_phaser.sh 107 CANCELLED: never satisfiable, upstream failure: run_02_imbalance.sh 104; phaser_count.sh 106
phaser_count.sh 106 finished; first error line of its log: FATAL ERROR - no reads found" "107" 1
J8='prep_genotypes.sh\t101\nalign_wasp_count.sh\t102\nrun_01_import_qc.sh\t103\nrun_02_imbalance.sh\t104\nsetup_phaser_env.sh\t105\nphaser_count.sh\t106\nrun_05_phaser.sh\t107\n'
wait_case install_fails "$J8" "ase_checkpoint.rds summary_numbers.tsv" "setup_phaser_env_105.out|ERROR: conda env create failed" \
  "106_[1-3] DependencyNeverSatisfied
107 DependencyNeverSatisfied" "setup_phaser_env.sh 105 finished; first error line of its log: ERROR: conda env create failed
phaser_count.sh 106 CANCELLED: never satisfiable, upstream failure: setup_phaser_env.sh 105
run_05_phaser.sh 107 CANCELLED: never satisfiable, upstream failure: setup_phaser_env.sh 105
run_02_imbalance.sh 104 ok" "106 107" 1
wait_case array_fails "$J8" "" "align_wasp_count_102_3.out|ERROR: STAR failed for sample s3" \
  "103 DependencyNeverSatisfied
104 DependencyNeverSatisfied
106_[1-3] DependencyNeverSatisfied
107 DependencyNeverSatisfied" "align_wasp_count.sh 102 finished; first error line of its log: ERROR: STAR failed for sample s3
run_01_import_qc.sh 103 CANCELLED: never satisfiable, upstream failure: align_wasp_count.sh 102
run_02_imbalance.sh 104 CANCELLED: never satisfiable, upstream failure: align_wasp_count.sh 102
phaser_count.sh 106 CANCELLED: never satisfiable, upstream failure: align_wasp_count.sh 102
run_05_phaser.sh 107 CANCELLED: never satisfiable, upstream failure: align_wasp_count.sh 102" "103 104 106 107" 1

if [ -s "$W/forbidden" ]; then echo "FAIL: a forbidden command was called:"; cat "$W/forbidden"; fail=1; fi
n_sleep=$(cat "$W"/*/sleep 2>/dev/null | grep -c '^60$')   # the stub itself is bin/sleep and holds no such line
[ "$n_sleep" = 1 ] && echo "sleep stub: called once (rmd02_blocks_rmd05: one polling round while Rmd 04 still ran), no real wait" ||
  { echo "FAIL: expected exactly one 'sleep 60' (rmd02_blocks_rmd05), got $n_sleep"; fail=1; }
[ -n "$OWN_TMP" ] && rm -rf "$OWN_TMP"   # only the directory this run created with mktemp -d
[ $fail -eq 0 ] && echo "DRY RUN PASS" || exit 1
