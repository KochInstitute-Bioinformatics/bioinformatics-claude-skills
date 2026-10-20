#!/bin/bash
# Test of the SLURM generator of the skill (Step 10). The generator is a few short Bash blocks and two template files: the blocks
# and the templates are cut out by their headings (the templates are written to templates/ as the Write tool would write them), the
# placeholders of the blocks are replaced by test values, and each block is run under env -i on a scratch tree, in order, as the
# assistant of a user would. Asserts the header of every script, the R call, the jobs that exist for each set of modules and
# sections, the chain in run_all.sh (run against a stub sbatch that only records its arguments), the refusals, and that nothing is
# written to the pipeline output. Also tests the Step 3 block that makes the output directory absolute and the Step 6 block (exit
# status without a pairs file), and that every Bash block of Step 10 is short and has no heredoc. No real sbatch, singularity or R.
set -u
here=$(cd "$(dirname "$0")" && pwd)
skill=${1:-$here/../nfcore-rnasplice-report.md}
tmp=$(mktemp -d "${TMPDIR:-/tmp}/rsr_slurm.XXXXXX") || exit 1
case "$tmp" in /*/rsr_slurm.??????) ;; *) echo "unexpected temp dir: $tmp"; exit 1 ;; esac
cleanup() { chmod -R u+w "$tmp" 2>/dev/null; rm -rf "$tmp"; }   # the directory was made by mktemp above and checked against its pattern
trap cleanup EXIT
fail=0
ok()    { echo "ok   $1"; }
bad()   { echo "FAIL $1"; fail=1; }
check() { if [ "$2" = "$3" ]; then ok "$1"; else bad "$1: expected '$2' got '$3'"; fi; }
has()   { if grep -qF -- "$2" "$3" 2>/dev/null; then ok "$1"; else bad "$1: '$2' not found in $3"; fi; }
hasnt() { if grep -qF -- "$2" "$3" 2>/dev/null; then bad "$1: '$2' found in $3"; else ok "$1"; fi; }
SIFP=/srv/site_root/images/bulkrnaseq_latest.sif
SYS=/usr/bin:/bin
PART=testpart; SINGL='module add singularity/3.10.4'; SBIND=/srv/site_root; BIMG=docker://example/bulkrnaseq:latest
blk() { bash "$here/cut_block.sh" "$skill" "$1"; }
B_CHECK=$(blk '**Checks and job list.**') && B_SET=$(blk '**Settings.**') && B_REND=$(blk '**Rendering the scripts.**') \
  && blk '**Job script template.**' > "$tmp/run_job.sbatch.tpl" && blk '**Submit script template.**' > "$tmp/run_all.tpl" \
  || { bad "a block or template was not found after its heading (Checks and job list, Job script template, Submit script template, Settings, Rendering the scripts)"; echo "SLURM RENDER FAIL"; exit 1; }
[ -s "$tmp/run_job.sbatch.tpl" ] && [ -s "$tmp/run_all.tpl" ] && [ -n "$B_CHECK$B_SET$B_REND" ] || { bad "a block or template is empty"; echo "SLURM RENDER FAIL"; exit 1; }
# A fake pipeline output (read-only afterwards) and a fake report directory under results/.
mkdir -p "$tmp/wd/out/star" "$tmp/wd/results"
echo x > "$tmp/wd/out/star/keep.txt"
chmod -R a-w "$tmp/wd/out"
echo "a,b" > "$tmp/wd/contrasts.csv"; : > "$tmp/wd/samples.csv"
# Stub tools, type-checked: the generated scripts are never executed, only run_all.sh against the sbatch stub.
mkdir -p "$tmp/stub"
cat > "$tmp/stub/sbatch" <<'STUB'
#!/bin/bash
n=$(( $(wc -l < "$STUBLOG") + 1 ))
echo "sbatch $* -> $((1000 + n))" >> "$STUBLOG"
echo $((1000 + n))
STUB
chmod +x "$tmp/stub/sbatch"
[ "$(env -i PATH="$tmp/stub:$SYS" /bin/bash -c 'type -P sbatch')" = "$tmp/stub/sbatch" ] || { bad "stub sbatch is not first on the PATH"; exit 1; }
# subst <text>: replaces the placeholders of the skill by the test values (the variables of render).
PROJP=wd; CWDP=$tmp/wd; GTFP=$tmp/wd/genes.gtf; CSVP=$tmp/wd/contrasts.csv
se() { printf '%s' "$1" | sed 's/[\\&|]/\\&/g'; }
subst() {
  printf '%s\n' "$1" | sed \
    -e "s|{REPORT_DIR}|$RD|g" -e "s|{OUTDIR_RNASPLICE}|${OUTX:-$tmp/wd/out}|g" -e "s|{CWD}|$(se "$CWDP")|g" \
    -e "s|{SIF}|$(se "$SIFP")|g" -e "s|{EMAIL}|$(se "$EMAIL")|g" \
    -e "s|{SLURM_PARTITION}|$PART|g" -e "s|{SINGULARITY_LOAD}|$SINGL|g" -e "s|{SINGULARITY_BIND}|$SBIND|g" -e "s|{BULK_IMAGE}|$BIMG|g" \
    -e "s|{PROJECT}|$(se "$PROJP")|g" -e "s|{MEM}|$MEM|g" -e "s|{TIME}|$TIME|g" -e "s|{MODULES}|$mods|g" -e "s|{SECTIONS}|$secs|g" \
    -e "s|{FDR}|0.05|g" -e "s|{DPSI}|0.1|g" -e "s|{RMATS_TYPE}|JC|g" -e "s|{MIN_READS}|10|g" -e "s|{ORGANISM}|human|g" \
    -e "s|{ENRICH_DBS}|GO KEGG|g" -e "s|{GTF}|$(se "$GTFP")|g" -e "s|{CONTRASTS_CSV}|$(se "$CSVP")|g" \
    -e "s|{SAMPLESHEET}|$tmp/wd/samples.csv|g" -e "s|{PAIRED_RMATS}|false|g" -e "s|{PAIRED_SUPPA}|true|g"
}
runb() { env -i PATH="$SYS" /bin/bash --noprofile --norc -c "$1"; }
# render <name> <modules> <sections> [var=value ...]: creates the empty Rmd files (except the one named by DROP), runs the checks block,
# writes the templates, runs the settings and rendering blocks; the output of all blocks is in $tmp/<name>.out, the exit status of
# the first failing block in $tmp/<name>.rc, the files in $tmp/<name>/. Later blocks are not run after an ERROR, as the skill says.
render() {
  local name=$1 mods=$2 secs=$3; shift 3
  local R="$tmp/wd/results/2026-10-05_wd_report_$name"
  local MEM=32G TIME=4:00:00 EMAIL=x@y.z RD=$R DROP="" OUTX=$tmp/wd/out kv f prog rc=0 PROJP=$PROJP CWDP=$tmp/wd GTFP=$tmp/wd/genes.gtf CSVP=$tmp/wd/contrasts.csv SIFP=$SIFP PART=$PART SINGL="$SINGL" SBIND=$SBIND
  for kv in "$@"; do printf -v "${kv%%=*}" %s "${kv#*=}"; done
  for prog in "$B_CHECK" "$B_SET" "$B_REND"; do
    subst "$prog" | grep -o "[^$]{[A-Z][A-Z_]*}" | sort -u
  done > "$tmp/$name.unreplaced"
  mkdir -p "$RD" 2>/dev/null
  for f in 01_import 02_rmats 03_suppa 04_usage 05_cross_tool 06_enrichment 07_summary; do { : > "$RD/$f.Rmd"; } 2>/dev/null; done
  if [ -n "$DROP" ]; then rm -f "$RD/$DROP.Rmd"; fi
  : > "$tmp/$name.out"
  runb "$(subst "$B_CHECK")" >> "$tmp/$name.out" 2>&1; rc=$?
  if [ $rc -eq 0 ]; then
    cp "$tmp/run_job.sbatch.tpl" "$tmp/run_all.tpl" "$RD/templates/"
    runb "$(subst "$B_SET")" >> "$tmp/$name.out" 2>&1; rc=$?
  fi
  if [ $rc -eq 0 ]; then runb "$(subst "$B_REND")" >> "$tmp/$name.out" 2>&1; rc=$?; fi
  echo "$rc" > "$tmp/$name.rc"
  ln -sfn "$RD" "$tmp/$name"
}
jobs_of() { sed -n 's/^JOBS=//p' "$tmp/$1.out" | tail -n 1; }
ALL="rmats suppa dexseq_exon edger_exon dtu"
# --- all modules, events only
render ev "$ALL" "events"
check "ev: exit status" 0 "$(cat "$tmp/ev.rc")"
check "ev: jobs" "01 02 03 04 07" "$(jobs_of ev)"
# --- all modules, all four sections
render all "$ALL" "events cross topgene enrich"
check "all: exit status" 0 "$(cat "$tmp/all.rc")"
check "all: jobs" "01 02 03 04 05 06 07" "$(jobs_of all)"
# the job lists differ exactly by the Rmd 05 and 06 jobs
for n in 05 06; do
  [ -f "$tmp/all/run_$n.sbatch" ] && ok "all: run_$n.sbatch exists" || bad "all: run_$n.sbatch is missing"
  [ ! -e "$tmp/ev/run_$n.sbatch" ] && ok "ev: no run_$n.sbatch" || bad "ev: run_$n.sbatch must not exist"
done
# --- header and R call of every script of the full set
for n in 01 02 03 04 05 06 07; do
  f="$tmp/all/run_$n.sbatch"
  has "run_$n: shebang" '#!/bin/bash' "$f"
  has "run_$n: job name" "#SBATCH -J wd_report_$n" "$f"
  has "run_$n: partition" '#SBATCH -p testpart' "$f"
  has "run_$n: one node, one task, four cpus" '#SBATCH -N 1' "$f"; has "run_$n: tasks" '#SBATCH -n 1' "$f"; has "run_$n: cpus" '#SBATCH -c 4' "$f"
  has "run_$n: memory" '#SBATCH --mem=32G' "$f"
  has "run_$n: time" '#SBATCH -t 4:00:00' "$f"
  has "run_$n: log under the report directory" "#SBATCH -o $tmp/wd/results/2026-10-05_wd_report_all/logs/run_$n." "$f"
  has "run_$n: mail user" '#SBATCH --mail-user=x@y.z' "$f"
  has "run_$n: mail type" '#SBATCH --mail-type=END,FAIL' "$f"
  has "run_$n: module guarded" 'module add singularity/3.10.4 || exit 1' "$f"
  has "run_$n: container" "singularity exec --bind '/srv/site_root" "$f"
  has "run_$n: Rscript render" "Rscript -e 'rmarkdown::render(" "$f"
  has "run_$n: output_dir under results" "output_dir = \"$tmp/wd/results/2026-10-05_wd_report_all\"" "$f"
  bash -n "$f" && ok "run_$n: bash -n" || bad "run_$n: bash -n"
done
for pair in 01:01_import 02:02_rmats 03:03_suppa 04:04_usage 05:05_cross_tool 06:06_enrichment 07:07_summary; do
  has "run_${pair%%:*}: its Rmd" "${pair#*:}.Rmd" "$tmp/all/run_${pair%%:*}.sbatch"
done
has "run_07: the summary page is written to report_summary.html" "output_file = \"report_summary.html\"" "$tmp/all/run_07.sbatch"
hasnt "run_01: no output_file" "output_file" "$tmp/all/run_01.sbatch"
has "params: fdr" 'fdr = 0.05' "$tmp/all/run_02.sbatch"
has "params: dpsi" 'dpsi = 0.1' "$tmp/all/run_02.sbatch"
has "params: rmats type" "rmats_type = \"JC\"" "$tmp/all/run_02.sbatch"
has "params: min reads" 'min_reads = 10' "$tmp/all/run_02.sbatch"
has "params: events on" 'events = TRUE' "$tmp/all/run_02.sbatch"
has "params: topgene on" 'topgene = TRUE' "$tmp/all/run_02.sbatch"
has "params: events only: topgene off" 'topgene = FALSE' "$tmp/ev/run_02.sbatch"
has "params: events only: events on" 'events = TRUE' "$tmp/ev/run_02.sbatch"
has "params: pairing" 'paired_suppa = TRUE' "$tmp/all/run_03.sbatch"
# --- no mail lines without an address
render nomail "rmats" "events" EMAIL=none
hasnt "nomail: no mail-user" 'mail-user' "$tmp/nomail/run_01.sbatch"
hasnt "nomail: no mail-type" 'mail-type' "$tmp/nomail/run_01.sbatch"
# --- no partition line without a partition, and an empty Singularity load becomes `true`
render nopart "rmats" "events" PART=
check "nopart: exit status" 0 "$(cat "$tmp/nopart.rc")"
hasnt "nopart: no partition line" '#SBATCH -p' "$tmp/nopart/run_01.sbatch"
has "nopart: the job script is otherwise complete" '#SBATCH -N 1' "$tmp/nopart/run_01.sbatch"
render nosing "rmats" "events" SINGL=
check "nosing: exit status" 0 "$(cat "$tmp/nosing.rc")"
has "nosing: empty SINGULARITY_LOAD becomes true" 'true || exit 1' "$tmp/nosing/run_01.sbatch"
# --- other module sets
render rm "rmats" "events cross topgene enrich"
check "rmats only: jobs (one tool: no cross-tool job)" "01 02 06 07" "$(jobs_of rm)"
render rs "rmats suppa" "events cross"
check "rmats and suppa: cross-tool job" "01 02 03 05 07" "$(jobs_of rs)"
render dtu "dtu" "events cross enrich"
check "dtu only: jobs" "01 04 06 07" "$(jobs_of dtu)"
render ex "dexseq_exon edger_exon" "events cross"
check "two exon tools count as two tools" "01 04 05 07" "$(jobs_of ex)"
render none "rmats suppa" ""
check "no sections: jobs" "01 02 03 07" "$(jobs_of none)"
# --- run_all.sh: chain (02-06 afterok on 01 and cancelled when 01 fails, 07 afterany on 01-06), light only
f="$tmp/all/run_all.sh"
has "run_all: dependency afterok" 'dependency=afterok:' "$f"
has "run_all: 01 has no dependency" 'JOB01=$(sbatch --parsable run_01.sbatch)' "$f"
has "run_all: 02 after 01" 'JOB02=$(sbatch --parsable --dependency=afterok:$JOB01 --kill-on-invalid-dep=yes run_02.sbatch)' "$f"
has "run_all: 06 after 01" 'JOB06=$(sbatch --parsable --dependency=afterok:$JOB01 --kill-on-invalid-dep=yes run_06.sbatch)' "$f"
has "run_all: summary after all others (afterany)" 'JOB07=$(sbatch --parsable --dependency=afterany:$JOB01:$JOB02:$JOB03:$JOB04:$JOB05:$JOB06 run_07.sbatch)' "$f"
hasnt "run_all: 07 not afterok" 'afterok:$JOB01:' "$f"
hasnt "run_all: no Rscript" 'Rscript' "$f"
hasnt "run_all: no singularity" 'singularity' "$f"
hasnt "run_all: no R" 'rmarkdown' "$f"
hasnt "ev: run_all has no 05" 'run_05' "$tmp/ev/run_all.sh"
has "ev: summary after 01 02 03 04" 'dependency=afterany:$JOB01:$JOB02:$JOB03:$JOB04 run_07.sbatch' "$tmp/ev/run_all.sh"
bash -n "$f" && ok "run_all: bash -n" || bad "run_all: bash -n"
# run it against the stub: seven submissions, the ids chained
STUBLOG="$tmp/stub.log"; : > "$STUBLOG"
env -i PATH="$tmp/stub:$SYS" STUBLOG="$STUBLOG" /bin/bash --noprofile --norc "$f" > "$tmp/runall.out" 2>&1
check "run_all against the stub: seven sbatch calls" 7 "$(wc -l < "$STUBLOG")"
has "stub: 02 depends on the id of 01" 'sbatch --parsable --dependency=afterok:1001 --kill-on-invalid-dep=yes run_02.sbatch' "$STUBLOG"
has "stub: 07 depends on the ids of 01-06" 'sbatch --parsable --dependency=afterany:1001:1002:1003:1004:1005:1006 run_07.sbatch' "$STUBLOG"
# --- only 01 and 07 (no module): 01 has no dependency, 07 waits for 01
render only "" ""
check "no modules: exit status" 0 "$(cat "$tmp/only.rc")"
check "no modules: jobs" "01 07" "$(jobs_of only)"
has "no modules: 01 without dependency" 'JOB01=$(sbatch --parsable run_01.sbatch)' "$tmp/only/run_all.sh"
has "no modules: 07 after 01 (afterany)" 'JOB07=$(sbatch --parsable --dependency=afterany:$JOB01 run_07.sbatch)' "$tmp/only/run_all.sh"
# --- all 18 parameters of the Rmd contract are in the render call of every job
for p in "outdir = \"$tmp/wd/out\"" "report_dir = \"$tmp/wd/results/2026-10-05_wd_report_all\"" "contrasts = \"$tmp/wd/contrasts.csv\"" \
  "samplesheet = \"$tmp/wd/samples.csv\"" "modules = \"$ALL\"" 'fdr = 0.05' 'dpsi = 0.1' "rmats_type = \"JC\"" 'min_reads = 10' \
  'events = TRUE' 'cross = TRUE' 'topgene = TRUE' 'enrich = TRUE' "organism = \"human\"" "enrich_dbs = \"GO KEGG\"" \
  "gtf = \"$tmp/wd/genes.gtf\"" 'paired_rmats = FALSE' 'paired_suppa = TRUE'; do
  for n in 01 04 07; do has "param $p in run_$n" "$p" "$tmp/all/run_$n.sbatch"; done
done
# --- a missing Rmd stops the checks before anything is written (only the directories are created)
render norm "rmats" "events" DROP=02_rmats
check "missing Rmd: exit status" 1 "$(cat "$tmp/norm.rc")"
has "missing Rmd: message" "02_rmats.Rmd does not exist: write the Rmd files first" "$tmp/norm.out"
check "missing Rmd: no run_* file written" "" "$(ls "$tmp/wd/results/2026-10-05_wd_report_norm" | grep '^run_' | tr '\n' ' ')"
check "missing Rmd: no job list written" "" "$(ls "$tmp/wd/results/2026-10-05_wd_report_norm/templates" 2>/dev/null | grep jobs.txt)"
# --- the rendering block refuses to run without the earlier steps
mkdir -p "$tmp/wd/results/2026-10-05_wd_report_early"
RD="$tmp/wd/results/2026-10-05_wd_report_early" mods=rmats secs=events MEM=32G TIME=4:00:00 EMAIL=none
runb "$(subst "$B_REND")" > "$tmp/early.out" 2>&1; check "rendering without the earlier steps is refused" 1 "$?"
check "rendering without the earlier steps writes no script" "" "$(ls "$tmp/wd/results/2026-10-05_wd_report_early" | tr '\n' ' ')"
# --- resources stay within the default limits, and the generated scripts hold only the requested values
for n in 01 02 03 04 05 06 07; do
  m=$(sed -n 's/^#SBATCH --mem=\([0-9]*\)G$/\1/p' "$tmp/all/run_$n.sbatch"); t=$(sed -n 's/^#SBATCH -t \([0-9]*\):.*/\1/p' "$tmp/all/run_$n.sbatch")
  [ -n "$m" ] && [ "$m" -le 64 ] && [ -n "$t" ] && [ "$t" -le 4 ] && ok "run_$n: at most 64G and 4h" || bad "run_$n: memory '$m' G or time '$t' h out of the default limits"
done
# --- refusals
render badmem "rmats" "events" MEM='32G;rm'
check "bad memory value is refused" 1 "$(cat "$tmp/badmem.rc")"
render badtime "rmats" "events" TIME='4h'
check "bad time value is refused" 1 "$(cat "$tmp/badtime.rc")"
render big "rmats" "events" MEM=96G TIME=8:00:00
check "values above the default limits are written when given" 0 "$(cat "$tmp/big.rc")"
has "big: the note about the limits" 'above the default limits' "$tmp/big.out"
render outside "rmats" "events" RD="$tmp/wd/elsewhere"
check "report directory outside results/ is refused" 1 "$(cat "$tmp/outside.rc")"
has "outside: the message" 'not under a results/ directory' "$tmp/outside.out"
render inside "rmats" "events" RD="$tmp/wd/out/results/r"
check "report directory inside the pipeline output is refused" 1 "$(cat "$tmp/inside.rc")"
has "inside: the message" "is inside the pipeline output $tmp/wd/out" "$tmp/inside.out"
# --- a path outside /srv/site_root is bound into the container: its directory is in the --bind list of every script
render bindx "rmats" "events" GTFP=/opt/rsr_ext/genes.gtf
check "bind outside /srv/site_root: exit status" 0 "$(cat "$tmp/bindx.rc")"
for n in 01 02 07; do
  bl=$(sed -n "s/.*--bind '\([^']*\)'.*/\1/p" "$tmp/bindx/run_$n.sbatch")
  case ",$bl," in *",/opt/rsr_ext,"*) ok "run_$n: the directory of a gtf outside /srv/site_root is bound ($bl)" ;; *) bad "run_$n: /opt/rsr_ext is not in the --bind list '$bl'" ;; esac
  case "$bl" in /srv/site_root,*|/srv/site_root) ok "run_$n: /srv/site_root is bound first" ;; *) bad "run_$n: the --bind list '$bl' does not start with /srv/site_root" ;; esac
done
# --- the final line names the one submit command, and the pipeline output is unchanged
has "ev: one submit command" "SUBMIT=bash $tmp/wd/results/2026-10-05_wd_report_ev/run_all.sh" "$tmp/ev.out"
check "nothing written to the pipeline output" "$tmp/wd/out/star/keep.txt" "$(find "$tmp/wd/out" -type f | tr '\n' ' ' | sed 's/ $//')"
# --- special characters in values: either the scripts hold the value exactly, or the generator fails loudly (never a broken script)
render amp "rmats" "events" EMAIL='a&b@c.d' GTFP="$tmp/wd/g&h.gtf"
check "ampersand in the e-mail and the gtf path: exit status" 0 "$(cat "$tmp/amp.rc")"
has "ampersand in the e-mail" '#SBATCH --mail-user=a&b@c.d' "$tmp/amp/run_01.sbatch"
has "ampersand in the gtf path" "gtf = \"$tmp/wd/g&h.gtf\"" "$tmp/amp/run_07.sbatch"
render pipe "rmats" "events" CSVP="$tmp/wd/c|d.csv"
check "pipe in a path: exit status" 0 "$(cat "$tmp/pipe.rc")"
has "pipe in the contrasts path" "contrasts = \"$tmp/wd/c|d.csv\"" "$tmp/pipe/run_01.sbatch"
render bsl "rmats" "events" CSVP="$tmp/wd/c\\d.csv"
check "backslash in a value is refused: exit status" 1 "$(cat "$tmp/bsl.rc")"
has "backslash in a value: message" "ERROR" "$tmp/bsl.out"
check "backslash in a value: no run script written" "" "$(ls "$tmp/wd/results/2026-10-05_wd_report_bsl" | grep '^run_' | tr '\n' ' ')"
render dqv "rmats" "events" GTFP="$tmp/wd/g\"h.gtf"
check "double quote in a value is refused: exit status" 1 "$(cat "$tmp/dqv.rc")"
check "double quote in a value: no run script written" "" "$(ls "$tmp/wd/results/2026-10-05_wd_report_dqv" | grep '^run_' | tr '\n' ' ')"
render spc "rmats" "events" CWDP="$tmp/wd sp"
check "space in the working directory (it would be bound): exit status" 1 "$(cat "$tmp/spc.rc")"
has "space in the working directory: message" "ERROR" "$tmp/spc.out"
check "space in the working directory: no run script written" "" "$(ls "$tmp/wd/results/2026-10-05_wd_report_spc" | grep '^run_' | tr '\n' ' ')"
# --- shell syntax in values is taken literally: no expansion, no command runs (canary files), backslashes kept
render lit1 "rmats" "events" EMAIL='a$HOME@b.c' GTFP="$tmp/wd/g\`>$tmp/CANARY1\`.gtf" CSVP="$tmp/wd/\$(>$tmp/CANARY2).csv"
check "dollar, backtick and command substitution in values: exit status" 0 "$(cat "$tmp/lit1.rc")"
has "dollar sign kept in the e-mail" '#SBATCH --mail-user=a$HOME@b.c' "$tmp/lit1/run_01.sbatch"
has "backtick command kept in the gtf path" "gtf = \"$tmp/wd/g\`>$tmp/CANARY1\`.gtf\"" "$tmp/lit1/run_01.sbatch"
has "command substitution kept in the contrasts path" "contrasts = \"$tmp/wd/\$(>$tmp/CANARY2).csv\"" "$tmp/lit1/run_01.sbatch"
[ ! -e "$tmp/CANARY1" ] && [ ! -e "$tmp/CANARY2" ] && ok "no command of a value was run" || bad "a command inside a value was run (canary file exists)"
render lit2 "rmats" "events" EMAIL="a\$(>$tmp/CANARY3)@b.c"
has "command substitution kept in the e-mail" "--mail-user=a\$(>$tmp/CANARY3)@b.c" "$tmp/lit2/run_01.sbatch"
[ ! -e "$tmp/CANARY3" ] && ok "no command of the e-mail was run" || bad "the e-mail command was run"
# --- the generated scripts are run (stubbed sbatch, module and singularity): values with $HOME, a backtick command and $(...) in the
# report directory, the image, the gtf and the contrasts path reach the Rscript command line byte for byte, and no command of a value runs
printf '#!/bin/bash\nexit 0\n' > "$tmp/stub/module"
printf '#!/bin/bash\nprintf "%%s\\n" "$@" > "$SLOG"\n' > "$tmp/stub/singularity"
chmod +x "$tmp/stub/module" "$tmp/stub/singularity"
for t in module singularity sbatch; do [ "$(env -i PATH="$tmp/stub:$SYS" /bin/bash -c "type -P $t")" = "$tmp/stub/$t" ] || bad "stub $t is not first on the PATH"; done
XRD="$tmp/wd/results/r\$HOME\`>$tmp/CANARY4\`"
XSIF="/srv/site_root/s\$HOME\$(>$tmp/CANARY5).sif"
XGTF="$tmp/wd/g\`>$tmp/CANARY6\`&|.gtf"
XCSV="$tmp/wd/\$(>$tmp/CANARY7).csv"
render exec "rmats suppa" "events" RD="$XRD" SIFP="$XSIF" GTFP="$XGTF" CSVP="$XCSV"
check "values with shell syntax: exit status" 0 "$(cat "$tmp/exec.rc")"
STUBLOG="$tmp/exec.stub.log"; : > "$STUBLOG"
env -i PATH="$tmp/stub:$SYS" STUBLOG="$STUBLOG" /bin/bash --noprofile --norc "$tmp/exec/run_all.sh" > "$tmp/exec.runall.out" 2>&1
check "run_all.sh with shell syntax in the report directory: exit status" 0 "$?"
check "run_all.sh submitted the four jobs" 4 "$(wc -l < "$STUBLOG")"
SLOG="$tmp/exec.sing.log"; : > "$SLOG"
env -i PATH="$tmp/stub:$SYS" SLOG="$SLOG" /bin/bash --noprofile --norc "$tmp/exec/run_01.sbatch" > "$tmp/exec.job.out" 2>&1
check "run_01.sbatch with shell syntax in the values: exit status" 0 "$?"
has "image argument byte-exact" "$XSIF" "$SLOG"
has "report directory byte-exact in the render call" "rmarkdown::render(\"$XRD/01_import.Rmd\"," "$SLOG"
has "gtf byte-exact in the params" "gtf = \"$XGTF\"" "$SLOG"
has "contrasts path byte-exact in the params" "contrasts = \"$XCSV\"" "$SLOG"
for c in CANARY4 CANARY5 CANARY6 CANARY7; do [ ! -e "$tmp/$c" ] && ok "no command of a value ran ($c)" || bad "a command inside a value ran: $c exists"; done
# a value that would break a script, whatever the cause, stops the rendering: an empty script or a leftover marker is an error
render emptyv "rmats" "events"
printf 's|__R__|x|g\n' > "$tmp/wd/results/2026-10-05_wd_report_emptyv/templates/settings.sed"
RD=$tmp/wd/results/2026-10-05_wd_report_emptyv mods=rmats secs=events MEM=32G TIME=4:00:00 EMAIL=x@y.z
runb "$(subst "$B_REND")" > "$tmp/emptyv2.out" 2>&1; check "leftover marker in a rendered script: exit status" 1 "$?"
has "leftover marker: message" "ERROR" "$tmp/emptyv2.out"
# --- a project name or path with two underscores and capitals is a legitimate value, not a leftover marker
render legit "rmats" "events" PROJP='my__PROJECT__x' CSVP="$tmp/wd/a__B__c.csv"
check "value with __CAPITALS__: exit status" 0 "$(cat "$tmp/legit.rc")"
has "value with __CAPITALS__: job name" '#SBATCH -J my__PROJECT__x_report_01' "$tmp/legit/run_01.sbatch"
has "value with __CAPITALS__: contrasts path" "contrasts = \"$tmp/wd/a__B__c.csv\"" "$tmp/legit/run_01.sbatch"
# --- Step 2: the scripts carry the mail lines only when mail is wanted
for em in none x@y.z; do
  for pair in "**Image missing.**:pull_image" "**Package check.**:check_packages"; do
    anc=${pair%%:*}; nm=${pair#*:}
    blkx=$(blk "$anc") || { bad "Step 2 block after $anc not found"; continue; }
    RD="$tmp/wd/results/2026-10-05_wd_report_s2_${nm}_$em"; mods=rmats; secs=events; MEM=32G; TIME=4:00:00; EMAIL=$em
    STUBLOG="$tmp/s2.stub.log"; : > "$STUBLOG"
    env -i PATH="$tmp/stub:$SYS" STUBLOG="$STUBLOG" /bin/bash --noprofile --norc -c "$(subst "$blkx")" > "$tmp/s2.out" 2>&1; rc=$?
    check "Step 2 $nm with mail=$em: exit status" 0 "$rc"
    if [ "$em" = none ]; then
      hasnt "Step 2 $nm without mail: no mail-user line" 'mail-user' "$RD/image_check/$nm.sbatch"
      hasnt "Step 2 $nm without mail: no mail-type line" 'mail-type' "$RD/image_check/$nm.sbatch"
    else
      has "Step 2 $nm with mail: mail-user line" "#SBATCH --mail-user=$em" "$RD/image_check/$nm.sbatch"
    fi
  done
done
# --- Extracting the Rmd files: the seven blocks of the skill are cut out into the report directory, byte for byte
ex=$(blk '**Extracting the Rmd files.**') || bad "block after the heading Extracting the Rmd files not found"
mkdir -p "$tmp/xtmp"
# the extraction block deletes and moves nothing, and makes no temporary directory (validate, then write)
for w in rm rmdir mv mktemp trap unlink; do
  if printf '%s\n' "$ex" | grep -qw -- "$w"; then bad "the extraction block uses $w"; else ok "the extraction block does not use $w"; fi
done
runex() { env -i PATH="$SYS" TMPDIR="$tmp/xtmp" /bin/bash --noprofile --norc -c "$(printf '%s\n' "$ex" | sed -e "s|{REPORT_DIR}|$1|g" -e "s|{SKILL_FILE}|$2|g" -e "s|{OUTDIR_RNASPLICE}|${3:-$tmp/wd/out}|g")"; }
refrmd() { sed -n "/^\*\*File \`$1\.Rmd\`\.\*\*\$/,/^\`\`\`\`\$/p" "$skill" | sed '1,2d;$d'; }
EX="$tmp/wd/results/2026-10-05_wd_report_rmdx"
runex "$EX" "$skill" > "$tmp/ex.out" 2>&1; check "extraction: exit status" 0 "$?"
check "extraction: seven Rmd files and nothing else" "01_import.Rmd 02_rmats.Rmd 03_suppa.Rmd 04_usage.Rmd 05_cross_tool.Rmd 06_enrichment.Rmd 07_summary.Rmd" "$(ls -A "$EX" | tr '\n' ' ' | sed 's/ $//')"
for n in 01_import 02_rmats 03_suppa 04_usage 05_cross_tool 06_enrichment 07_summary; do
  refrmd "$n" > "$tmp/ref.$n"
  if cmp -s "$tmp/ref.$n" "$EX/$n.Rmd" && [ -s "$tmp/ref.$n" ]; then ok "extraction: $n.Rmd identical to the block"; else bad "extraction: $n.Rmd differs from the block"; fi
  has "extraction: line count of $n.Rmd printed" "$n.Rmd: $(wc -l < "$tmp/ref.$n") lines" "$tmp/ex.out"
done
mkdir -p "$tmp/wd/results/2026-10-05_wd_report_ex_nofile"; runex "$tmp/wd/results/2026-10-05_wd_report_ex_nofile" "$tmp/no_such_skill.md" > "$tmp/ex2.out" 2>&1; check "extraction: missing skill file: exit status" 1 "$?"
has "extraction: missing skill file: the path is named" "$tmp/no_such_skill.md" "$tmp/ex2.out"
check "extraction: missing skill file: nothing written" "" "$(ls -A "$tmp/wd/results/2026-10-05_wd_report_ex_nofile" 2>/dev/null)"
awk '/^\*\*File `07_summary.Rmd`\.\*\*$/ { exit } { print }' "$skill" > "$tmp/skill6.md"
mkdir -p "$tmp/wd/results/2026-10-05_wd_report_ex6"; runex "$tmp/wd/results/2026-10-05_wd_report_ex6" "$tmp/skill6.md" > "$tmp/ex3.out" 2>&1; check "extraction: a skill file with six Rmd blocks: exit status" 1 "$?"
has "extraction: six blocks: the message" "ERROR" "$tmp/ex3.out"
check "extraction: six blocks: nothing written" "" "$(ls -A "$tmp/wd/results/2026-10-05_wd_report_ex6" 2>/dev/null)"
sed '/^\*\*File `03_suppa.Rmd`\.\*\*$/{n;n;s/^---$/--/}' "$skill" > "$tmp/skillh.md"
mkdir -p "$tmp/wd/results/2026-10-05_wd_report_exh"; runex "$tmp/wd/results/2026-10-05_wd_report_exh" "$tmp/skillh.md" > "$tmp/ex4.out" 2>&1; check "extraction: a block that does not start with --- : exit status" 1 "$?"
check "extraction: bad header: nothing written" "" "$(ls -A "$tmp/wd/results/2026-10-05_wd_report_exh" 2>/dev/null)"
awk '/^\*\*File `04_usage.Rmd`\.\*\*$/ { print; getline; print; getline; skip = 1 } skip && /^````$/ { skip = 0 } !skip { print }' "$skill" > "$tmp/skille.md"
mkdir -p "$tmp/wd/results/2026-10-05_wd_report_exe"; runex "$tmp/wd/results/2026-10-05_wd_report_exe" "$tmp/skille.md" > "$tmp/ex5.out" 2>&1; check "extraction: an empty block: exit status" 1 "$?"
check "extraction: empty block: nothing written" "" "$(ls -A "$tmp/wd/results/2026-10-05_wd_report_exe" 2>/dev/null)"
sed 's/$/\r/' "$skill" > "$tmp/skillcr.md"
mkdir -p "$tmp/wd/results/2026-10-05_wd_report_excr"; runex "$tmp/wd/results/2026-10-05_wd_report_excr" "$tmp/skillcr.md" > "$tmp/ex6.out" 2>&1; check "extraction: a skill file with CRLF line endings: exit status" 1 "$?"
has "extraction: CRLF: the cause is named" "line endings are CRLF" "$tmp/ex6.out"
check "extraction: CRLF: nothing written" "" "$(ls -A "$tmp/wd/results/2026-10-05_wd_report_excr" 2>/dev/null)"
check "extraction: no temporary directory is left behind (success and every failure)" "" "$(ls -A "$tmp/xtmp")"
# a block whose closing fence is missing, or that holds another opener, is refused by name before anything is written
awk '/^\*\*File `03_suppa.Rmd`\.\*\*$/ { inb = 1 } inb && /^````$/ && !done { done = 1; next } { print }' "$skill" > "$tmp/skillc3.md"
last=$(grep -n '^````$' "$skill" | tail -n 1 | cut -d: -f1); sed "${last}d" "$skill" > "$tmp/skillc7.md"
awk '/^\*\*File `02_rmats.Rmd`\.\*\*$/ { h = 1 } h == 1 && /^````rmd$/ { h = 2; print; next } h == 2 { print; print "````rmd"; h = 3; next } { print }' "$skill" > "$tmp/skillop.md"
for pair in "c3:03_suppa" "c7:07_summary" "op:02_rmats"; do
  k=${pair%%:*}; nm=${pair#*:}; d="$tmp/wd/results/2026-10-05_wd_report_fence_$k"; mkdir -p "$d"
  runex "$d" "$tmp/skill$k.md" > "$tmp/exf$k.out" 2>&1; check "extraction: fence problem $k: exit status" 1 "$?"
  has "extraction: fence problem $k: the block is named" "$nm.Rmd" "$tmp/exf$k.out"
  has "extraction: fence problem $k: ERROR" "ERROR" "$tmp/exf$k.out"
  check "extraction: fence problem $k: nothing written" "" "$(ls -A "$d")"
done
check "extraction: no temporary directory after the fence cases" "" "$(ls -A "$tmp/xtmp")"
# the opening line of a block must be the line right after its heading: a missing opener (03 would otherwise get the content of 04), a
# text line or a blank line between the heading and the opener are refused, naming the block, before anything is written
awk '{ print } /^\*\*File `03_suppa.Rmd`\.\*\*$/ { getline }' "$skill" > "$tmp/skillno.md"
awk '{ print } /^\*\*File `03_suppa.Rmd`\.\*\*$/ { print "Some text." }' "$skill" > "$tmp/skilltx.md"
awk '{ print } /^\*\*File `03_suppa.Rmd`\.\*\*$/ { print "" }' "$skill" > "$tmp/skillbl.md"
check "fence fixture: the opener of 03 is deleted" 1 "$(( $(grep -c '^````rmd$' "$skill") - $(grep -c '^````rmd$' "$tmp/skillno.md") ))"
for k in no tx bl; do
  d="$tmp/wd/results/2026-10-05_wd_report_open_$k"; mkdir -p "$d"
  runex "$d" "$tmp/skill$k.md" > "$tmp/exo$k.out" 2>&1; check "extraction: opener problem $k: exit status" 1 "$?"
  has "extraction: opener problem $k: the block is named in the message" "ERROR: the line right after the heading of 03_suppa.Rmd is not the opening line "'````rmd'" of its block; nothing was written to $d" "$tmp/exo$k.out"
  check "extraction: opener problem $k: nothing written" "" "$(ls -A "$d")"
done
# --- the report directory is checked before anything is written into it: the block of Step 3, the extraction block and the checks
# block refuse a report directory inside the pipeline output (here a writable one, so that a missing check would write files)
mkdir -p "$tmp/wd/out2/star"
rdc=$(blk '**Report directory check.** Before anything') || bad "the block after the heading Report directory check not found"
runrd() { env -i PATH="$SYS" /bin/bash --noprofile --norc -c "$(printf '%s\n' "$rdc" | sed -e "s|{REPORT_DIR}|$1|g" -e "s|{OUTDIR_RNASPLICE}|$2|g")"; }
runrd "$tmp/wd/out2/results/r" "$tmp/wd/out2" > "$tmp/rdc1.out" 2>&1; check "report directory check: inside the pipeline output: exit status" 1 "$?"
has "report directory check: inside the pipeline output: message" "is inside the pipeline output $tmp/wd/out2" "$tmp/rdc1.out"
runrd "$tmp/wd/out2" "$tmp/wd/out2" > "$tmp/rdc2.out" 2>&1; check "report directory check: the pipeline output itself: exit status" 1 "$?"
runrd "$tmp/wd/elsewhere2" "$tmp/wd/out2" > "$tmp/rdc3.out" 2>&1; check "report directory check: outside results/: exit status" 1 "$?"
has "report directory check: outside results/: message" "not under a results/ directory" "$tmp/rdc3.out"
runrd "$tmp/wd/results/2026-10-05_wd_report_ok" "$tmp/wd/out2" > "$tmp/rdc4.out" 2>&1; check "report directory check: a directory next to the output: exit status" 0 "$?"
has "report directory check: OK line" "REPORT_DIR_OK=$tmp/wd/results/2026-10-05_wd_report_ok" "$tmp/rdc4.out"
check "report directory check: writes nothing" "" "$(ls -A "$tmp/wd/out2/star")"
runex "$tmp/wd/out2/results/r" "$skill" "$tmp/wd/out2" > "$tmp/exin.out" 2>&1; check "extraction: report directory inside the pipeline output: exit status" 1 "$?"
has "extraction: report directory inside the pipeline output: message" "is inside the pipeline output $tmp/wd/out2" "$tmp/exin.out"
[ ! -e "$tmp/wd/out2/results" ] && ok "extraction: nothing created inside the pipeline output" || bad "extraction: $tmp/wd/out2/results was created inside the pipeline output"
render inside2 "rmats" "events" RD="$tmp/wd/out2/results/r2" OUTX="$tmp/wd/out2"
check "checks block: report directory inside a writable pipeline output: exit status" 1 "$(cat "$tmp/inside2.rc")"
has "checks block: report directory inside the pipeline output: message" "is inside the pipeline output $tmp/wd/out2" "$tmp/inside2.out"
check "checks block: inside the pipeline output: no logs, templates or scripts written" "" "$(ls "$tmp/wd/out2/results/r2" | grep -v '\.Rmd$' | tr '\n' ' ')"
# --- the marker list of the rendering block names every marker of the templates
mk=$(printf '%s\n' "$B_REND" | sed -n "s/^MK='\(.*\)'\$/\1/p; s/^MK=\"\$MK|\(.*\)\"\$/|\1/p" | tr -d '\n')
tplm=$(cat "$tmp/run_job.sbatch.tpl" "$tmp/run_all.tpl" | grep -oE '__[A-Z]+(_[A-Z]+)*__' | sort -u)
check "the templates hold markers (the list is not empty)" 31 "$(printf "%s\n" $tplm | wc -l)"
miss=""; for m in $tplm; do case "|$mk|" in *"|$m|"*) ;; *) miss="$miss $m" ;; esac; done
check "every template marker is in the marker list of the rendering block" "" "$miss"
# --- no __PLACEHOLDER__ of the templates is left in a generated script
left=$(cat "$tmp"/all/run_*.sbatch "$tmp"/all/run_all.sh "$tmp"/only/run_all.sh | grep -o '__[A-Z_]*__' | sort -u | tr '\n' ' ')
check "no __PLACEHOLDER__ left in the generated scripts" "" "$left"
check "templates are kept in the report directory" "jobs.txt run_all.tpl run_job.sbatch.tpl settings.sed" "$(ls "$tmp/all/templates" | tr '\n' ' ' | sed 's/ $//')"
# --- every Bash block of Step 10 is short, has no function together with a heredoc, and at most two heredocs
s10=$(awk '/^## Step 10 /{p=1; next} /^## /{p=0} p' "$skill")
printf '%s\n' "$s10" | awk '
  /^```bash$/ { inb = 1; n = 0; next }
  inb && /^```$/ { inb = 0; print n; next }
  inb { n++ }' > "$tmp/s10.lens"
check "Step 10 has four Bash blocks" 4 "$(wc -l < "$tmp/s10.lens")"
check "no Bash block of Step 10 has more than 40 lines" 0 "$(awk '$1 > 40' "$tmp/s10.lens" | wc -l)"
check "no heredoc in a Bash block of Step 10" 0 "$(printf '%s\n' "$s10" | awk '/^```bash$/{inb=1; next} /^```$/{inb=0} inb && /<</' | wc -l)"
# --- Step 3: the output directory is made absolute. The default of the params file is relative to the directory that holds the
# params file (the run directory, {OUTDIR_BASE} = {PARAMS_DIR}), not to the directory the wizard runs in: the block runs from another
# directory ($tmp/elsewhere, which has a results_run/out of its own that must not be taken) and never reads {CWD}.
abs=$(blk '**Absolute output directory.**') || bad "Step 3 block with the heading Absolute output directory not found"
mkdir -p "$tmp/wd/results_run/out" "$tmp/elsewhere/results_run/out"
check "the Step 3 block resolves from {OUTDIR_BASE}, not from {CWD}" "0 1" "$(printf '%s\n' "$abs" | grep -c '{CWD}') $(printf '%s\n' "$abs" | grep -c "^cd '{OUTDIR_BASE}' ")"
runc() { (cd "$tmp/elsewhere" && env -i PATH="$SYS" /bin/bash --noprofile --norc -c "$1"); }
for given in "results_run/out" "$tmp/wd/results_run/out"; do
  out=$(runc "$(printf '%s\n' "$abs" | sed -e "s|{OUTDIR_BASE}|$tmp/wd|g" -e "s|{OUTDIR_GIVEN}|$given|g")" 2>&1 | sed -n 's/^OUTDIR_RNASPLICE=//p')
  check "absolute output directory from '$given' (run directory $tmp/wd, shell in $tmp/elsewhere)" "$tmp/wd/results_run/out" "$out"
done
runc "$(printf '%s\n' "$abs" | sed -e "s|{OUTDIR_BASE}|$tmp/wd|g" -e "s|{OUTDIR_GIVEN}|nowhere|g")" > "$tmp/abs.out" 2>&1; check "missing output directory: exit status" 1 "$?"
has "missing output directory: message" "ERROR: output directory nowhere not found from $tmp/wd" "$tmp/abs.out"
# --- derived values with a quote, a backslash or a line break stop the block that derives them, loudly and before anything is written
mkdir -p "$tmp/wd/o\"q"
runc "$(printf '%s\n' "$abs" | sed -e "s|{OUTDIR_BASE}|$tmp/wd|g" -e "s|{OUTDIR_GIVEN}|o\"q|g")" > "$tmp/absq.out" 2>&1; check "output directory with a double quote: exit status" 1 "$?"
has "output directory with a double quote: message" "ERROR: the output directory $tmp/wd/o\"q holds a single quote, a double quote, a backslash or a line break" "$tmp/absq.out"
hasnt "output directory with a double quote: no OUTDIR_RNASPLICE line" "OUTDIR_RNASPLICE=" "$tmp/absq.out"
s0=$(blk 'Record the current working directory') || bad "the Step 0 block after 'Record the current working directory' not found"
mkdir -p "$tmp/home" "$tmp/s0/wd'q" "$tmp/s0/wd\"d" "$tmp/s0/wd\\b" "$tmp/s0/ok" "$tmp/s0/ho'me"
runs0() { (cd "$1" && env -i PATH="$SYS" HOME="$2" /bin/bash --noprofile --norc -c "$s0"); }
runs0 "$tmp/s0/ok" "$tmp/home" > "$tmp/s0ok.out" 2>&1; check "Step 0 block: plain paths: exit status" 0 "$?"
has "Step 0 block: CWD" "CWD=$tmp/s0/ok" "$tmp/s0ok.out"; has "Step 0 block: WD_NAME" "WD_NAME=ok" "$tmp/s0ok.out"
has "Step 0 block: SKILL_FILE" "SKILL_FILE=$tmp/home/.claude/commands/nfcore-rnasplice-report.md" "$tmp/s0ok.out"
grep -qE '^TODAY_ISO=[0-9]{4}-[0-9]{2}-[0-9]{2}$' "$tmp/s0ok.out" && ok "Step 0 block: TODAY_ISO" || bad "Step 0 block: no TODAY_ISO=YYYY-MM-DD line"
for w in "wd'q" "wd\"d" "wd\\b"; do
  runs0 "$tmp/s0/$w" "$tmp/home" > "$tmp/s0bad.out" 2>&1; check "Step 0 block: working directory $w: exit status" 1 "$?"
  has "Step 0 block: working directory $w: message" "ERROR: the path $tmp/s0/$w holds a single quote" "$tmp/s0bad.out"
  hasnt "Step 0 block: working directory $w: no CWD line" "CWD=" "$tmp/s0bad.out"
  check "Step 0 block: working directory $w: nothing written" "" "$(ls -A "$tmp/s0/$w")"
done
runs0 "$tmp/s0/ok" "$tmp/s0/ho'me" > "$tmp/s0home.out" 2>&1; check "Step 0 block: a home directory with a quote (skill file path): exit status" 1 "$?"
has "Step 0 block: home with a quote: message" "ERROR: the path $tmp/s0/ho'me/.claude/commands/nfcore-rnasplice-report.md holds" "$tmp/s0home.out"
check "Step 0 blocks wrote nothing" "" "$(ls -A "$tmp/s0/ok")$(ls -A "$tmp/home")"
# the params block stops on a quote or a backslash inside a path value (the YAML quotes around a value are fine)
pp=$(blk 'First look for the params file, in one Bash call') || bad "the params block of Step 3 not found"
mkdir -p "$tmp/pq1" "$tmp/pq2" "$tmp/pq3" "$tmp/pq4"
printf 'input: "s.csv"\ncontrasts: '"'"'c.csv'"'"'\noutdir: results/x\ngtf: "/g/genes.gtf"\n' > "$tmp/pq1/x_params.yaml"
printf 'input: "s.csv"\noutdir: "results/it'"'"'s"\n' > "$tmp/pq2/x_params.yaml"
printf 'gtf: "/g/ge\\\\nes.gtf"\n' > "$tmp/pq3/x_params.yaml"
printf 'contrasts: "c\\"d.csv"\n' > "$tmp/pq4/x_params.yaml"
for k in 1 2 3 4; do runb "$(printf '%s\n' "$pp" | sed -e "s|{CWD}|$tmp/pq$k|g")" > "$tmp/pq$k.out" 2>&1; echo $? > "$tmp/pq$k.rc"; done
check "params block: quoted and plain values: exit status" 0 "$(cat "$tmp/pq1.rc")"
has "params block: the values are printed" 'outdir: results/x' "$tmp/pq1.out"
for k in 2:outdir 3:gtf 4:contrasts; do
  check "params block: a quote or backslash in $k: exit status" 1 "$(cat "$tmp/pq${k%%:*}.rc")"
  has "params block: a quote or backslash in $k: message" "ERROR: the value of ${k#*:} in the params file holds a quote or a backslash" "$tmp/pq${k%%:*}.out"
done
# Step 0: a report directory that exists gets the first numbered directory that does not exist yet
rdx=$(blk '**Report directory.** Check whether') || bad "the Step 0 block after Report directory not found"
FR="$tmp/wd/results/2026-10-05_wd_report_fr"; mkdir -p "$FR" "${FR}_2"
runb "$(printf '%s\n' "$rdx" | sed -e "s|{REPORT_DIR}|$FR|g")" > "$tmp/fr.out" 2>&1; check "Step 0 report directory block: exit status" 0 "$?"
has "Step 0: the first free numbered directory (_2 exists)" "FREE ${FR}_3" "$tmp/fr.out"
hasnt "Step 0: _2 is not offered when it exists" "FREE ${FR}_2" "$tmp/fr.out"
runb "$(printf '%s\n' "$rdx" | sed -e "s|{REPORT_DIR}|$tmp/wd/results/2026-10-05_wd_report_fresh|g")" > "$tmp/fr2.out" 2>&1
has "Step 0: a new report directory" "NEW $tmp/wd/results/2026-10-05_wd_report_fresh" "$tmp/fr2.out"
check "Step 0 report directory block: nothing created" "" "$(ls -A "$FR"; [ -e "${FR}_3" ] && echo created)"
# --- Step 6: no pairs file is a normal run (exit status 0, one clear line)
pb=$(blk 'The pairing is read from the run, not asked for') || bad "Step 6 block not found"
mkdir -p "$tmp/nopairs"
echo "rmats_paired_stats: false" > "$tmp/nopairs/x_params.yaml"
runb "$(printf '%s\n' "$pb" | sed -e "s|{CWD}|$tmp/nopairs|g")" > "$tmp/pairs.out" 2>&1; check "Step 6 without a pairs file: exit status" 0 "$?"
has "Step 6 without a pairs file: the line" "no pairs file in $tmp/nopairs" "$tmp/pairs.out"
has "Step 6: the switch is printed" "rmats_paired_stats: false" "$tmp/pairs.out"
printf 'sample,subject\ns1,p1\n' > "$tmp/nopairs/x_pairs.csv"
runb "$(printf '%s\n' "$pb" | sed -e "s|{CWD}|$tmp/nopairs|g")" > "$tmp/pairs2.out" 2>&1; check "Step 6 with a pairs file: exit status" 0 "$?"
has "Step 6 with a pairs file: its content" "s1,p1" "$tmp/pairs2.out"
check "no placeholder left unreplaced in the block" "" "$(cat "$tmp"/*.unreplaced | tr "\n" " ")"
chmod -R u+w "$tmp/wd/out"
[ $fail -eq 0 ] && echo "SLURM RENDER PASS" || { echo "SLURM RENDER FAIL"; exit 1; }
