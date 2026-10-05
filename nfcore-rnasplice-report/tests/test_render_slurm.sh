#!/bin/bash
# Test of the SLURM generator of the skill (Step 10): the block after the anchor "**SLURM scripts.**" is cut out, its
# placeholders are replaced by test values, and it is run under env -i on a scratch tree. Asserts the header of every script, the R
# call, the jobs that exist for each set of modules and sections, the chain in run_all.sh (run against a stub sbatch that
# only records its arguments), and that nothing is written to the pipeline output. No real sbatch, singularity or R is started.
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
SYS=/usr/bin:/bin
block=$(bash "$here/cut_block.sh" "$skill" '**SLURM scripts.**') || { bad "template block not found after the anchor **SLURM scripts.**"; echo "SLURM RENDER FAIL"; exit 1; }
[ -n "$block" ] || { bad "template block is empty"; echo "SLURM RENDER FAIL"; exit 1; }
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
# render <name> <modules> <sections> [var=value ...]: runs the block, leaves the output in $tmp/<name>.out and the files in $tmp/<name>/
render() {
  local name=$1 mods=$2 secs=$3; shift 3
  local R="$tmp/wd/results/2026-10-05_wd_report_$name"
  local MEM=32G TIME=4:00:00 EMAIL=x@y.z RD=$R kv
  for kv in "$@"; do printf -v "${kv%%=*}" %s "${kv#*=}"; done
  local prog
  prog=$(printf '%s\n' "$block" | sed \
    -e "s|{REPORT_DIR}|$RD|g" -e "s|{OUTDIR_RNASPLICE}|$tmp/wd/out|g" -e "s|{CWD}|$tmp/wd|g" \
    -e "s|{SIF}|/net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif|g" -e "s|{EMAIL}|$EMAIL|g" \
    -e "s|{PROJECT}|wd|g" -e "s|{MEM}|$MEM|g" -e "s|{TIME}|$TIME|g" -e "s|{MODULES}|$mods|g" -e "s|{SECTIONS}|$secs|g" \
    -e "s|{FDR}|0.05|g" -e "s|{DPSI}|0.1|g" -e "s|{RMATS_TYPE}|JC|g" -e "s|{MIN_READS}|10|g" -e "s|{ORGANISM}|human|g" \
    -e "s|{ENRICH_DBS}|GO KEGG|g" -e "s|{GTF}|$tmp/wd/genes.gtf|g" -e "s|{CONTRASTS_CSV}|$tmp/wd/contrasts.csv|g" \
    -e "s|{SAMPLESHEET}|$tmp/wd/samples.csv|g" -e "s|{PAIRED_RMATS}|false|g" -e "s|{PAIRED_SUPPA}|true|g")
  printf "%s\n" "$prog" | grep -o "[^$]{[A-Z][A-Z_]*}" | sort -u > "$tmp/$name.unreplaced"
  mkdir -p "$RD" 2>/dev/null
  local f
  for f in 01_import 02_rmats 03_suppa 04_usage 05_cross_tool 06_enrichment 07_summary; do { : > "$RD/$f.Rmd"; } 2>/dev/null; done
  env -i PATH="$SYS" /bin/bash --noprofile --norc -c "$prog" > "$tmp/$name.out" 2>&1
  echo "$?" > "$tmp/$name.rc"
  ln -sfn "$RD" "$tmp/$name"
}
jobs_of() { sed -n 's/^JOBS=//p' "$tmp/$1.out"; }
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
  has "run_$n: partition" '#SBATCH -p bcc' "$f"
  has "run_$n: one node, one task, four cpus" '#SBATCH -N 1' "$f"; has "run_$n: tasks" '#SBATCH -n 1' "$f"; has "run_$n: cpus" '#SBATCH -c 4' "$f"
  has "run_$n: memory" '#SBATCH --mem=32G' "$f"
  has "run_$n: time" '#SBATCH -t 4:00:00' "$f"
  has "run_$n: log under the report directory" "#SBATCH -o $tmp/wd/results/2026-10-05_wd_report_all/logs/run_$n." "$f"
  has "run_$n: mail user" '#SBATCH --mail-user=x@y.z' "$f"
  has "run_$n: mail type" '#SBATCH --mail-type=END,FAIL' "$f"
  has "run_$n: module guarded" 'module add singularity/3.10.4 || exit 1' "$f"
  has "run_$n: container" 'singularity exec --bind /net/bmc-lab3' "$f"
  has "run_$n: Rscript render" 'Rscript -e "rmarkdown::render(' "$f"
  has "run_$n: output_dir under results" "output_dir = '$tmp/wd/results/2026-10-05_wd_report_all'" "$f"
  bash -n "$f" && ok "run_$n: bash -n" || bad "run_$n: bash -n"
done
for pair in 01:01_import 02:02_rmats 03:03_suppa 04:04_usage 05:05_cross_tool 06:06_enrichment 07:07_summary; do
  has "run_${pair%%:*}: its Rmd" "${pair#*:}.Rmd" "$tmp/all/run_${pair%%:*}.sbatch"
done
has "params: fdr" 'fdr = 0.05' "$tmp/all/run_02.sbatch"
has "params: dpsi" 'dpsi = 0.1' "$tmp/all/run_02.sbatch"
has "params: rmats type" "rmats_type = 'JC'" "$tmp/all/run_02.sbatch"
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
# --- run_all.sh: chain with afterok, light only
f="$tmp/all/run_all.sh"
has "run_all: dependency afterok" 'dependency=afterok:' "$f"
has "run_all: 01 has no dependency" 'JOB01=$(sbatch --parsable run_01.sbatch)' "$f"
has "run_all: 02 after 01" 'JOB02=$(sbatch --parsable --dependency=afterok:$JOB01 run_02.sbatch)' "$f"
has "run_all: 06 after 01" 'JOB06=$(sbatch --parsable --dependency=afterok:$JOB01 run_06.sbatch)' "$f"
has "run_all: summary after all others" 'JOB07=$(sbatch --parsable --dependency=afterok:$JOB02:$JOB03:$JOB04:$JOB05:$JOB06 run_07.sbatch)' "$f"
hasnt "run_all: no Rscript" 'Rscript' "$f"
hasnt "run_all: no singularity" 'singularity' "$f"
hasnt "run_all: no R" 'rmarkdown' "$f"
hasnt "ev: run_all has no 05" 'run_05' "$tmp/ev/run_all.sh"
has "ev: summary after 02 03 04" 'dependency=afterok:$JOB02:$JOB03:$JOB04 run_07.sbatch' "$tmp/ev/run_all.sh"
bash -n "$f" && ok "run_all: bash -n" || bad "run_all: bash -n"
# run it against the stub: seven submissions, the ids chained
STUBLOG="$tmp/stub.log"; : > "$STUBLOG"
env -i PATH="$tmp/stub:$SYS" STUBLOG="$STUBLOG" /bin/bash --noprofile --norc "$f" > "$tmp/runall.out" 2>&1
check "run_all against the stub: seven sbatch calls" 7 "$(wc -l < "$STUBLOG")"
has "stub: 02 depends on the id of 01" 'sbatch --parsable --dependency=afterok:1001 run_02.sbatch' "$STUBLOG"
has "stub: 07 depends on the ids of 02-06" 'sbatch --parsable --dependency=afterok:1002:1003:1004:1005:1006 run_07.sbatch' "$STUBLOG"
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
# --- the final line names the one submit command, and the pipeline output is unchanged
has "ev: one submit command" "SUBMIT=bash $tmp/wd/results/2026-10-05_wd_report_ev/run_all.sh" "$tmp/ev.out"
check "nothing written to the pipeline output" "$tmp/wd/out/star/keep.txt" "$(find "$tmp/wd/out" -type f | tr '\n' ' ' | sed 's/ $//')"
check "no placeholder left unreplaced in the block" "" "$(cat "$tmp"/*.unreplaced | tr "\n" " ")"
chmod -R u+w "$tmp/wd/out"
[ $fail -eq 0 ] && echo "SLURM RENDER PASS" || { echo "SLURM RENDER FAIL"; exit 1; }
