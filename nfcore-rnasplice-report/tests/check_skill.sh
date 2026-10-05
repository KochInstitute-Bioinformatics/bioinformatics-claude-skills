#!/bin/bash
# Static checks for nfcore-rnasplice-report.md. Usage: check_skill.sh <skill.md>
# Prints PASS (exit 0) or one FAIL line per failed check (exit 1). Uses only bash, grep, awk and sed; no network.
# The README checked is the README.md next to the skill file. CHECK_LIST_NEEDS=1 also prints "NEED <line> <text>" per need.
# CHECK_RUN_ALL=<file> checks that file instead of tests/run_all_tests.sh: a test hook of prove_mutations.sh (target runall) only.
# Exit codes: 0 PASS, 1 a check failed, 2 a fixture is missing.
# Call need/needr at top level (also inside a loop, `case` or `||`), never from inside a helper function: BASH_LINENO[0]
# would then be the line inside the helper, and prove_red.sh could not map the need to an added checker line.
set -u
SKILL=${1:?usage: check_skill.sh <skill.md>}
HERE=$(cd "$(dirname "$0")" && pwd); FIX="$HERE/fixtures"
README="$(dirname "$SKILL")/README.md"
fail=0
[ -s "$SKILL" ] || { echo "FAIL: skill file missing or empty: $SKILL"; exit 1; }
for f in sign_rows.tsv; do
  [ -s "$FIX/$f" ] || { echo "cannot read fixture $FIX/$f"; exit 2; }
done

need()   { [ -n "${CHECK_LIST_NEEDS:-}" ] && echo "NEED ${BASH_LINENO[0]} $1"; grep -qF -- "$1" "$SKILL" || { echo "FAIL: missing required text: $1"; fail=1; }; }
needr()  { [ -n "${CHECK_LIST_NEEDS:-}" ] && echo "NEED ${BASH_LINENO[0]} $1"; grep -qF -- "$1" "$README" 2>/dev/null || { echo "FAIL: README missing required text: $1"; fail=1; }; }
forbid() { ! grep -qF -- "$1" "$SKILL" || { echo "FAIL: forbidden text present: $1"; fail=1; }; }
# forbid_re <extended regex> <label>: no line of the skill may match the regex.
forbid_re() { ! grep -qE -- "$1" "$SKILL" || { echo "FAIL: forbidden pattern present: $2"; fail=1; }; }

# --- skeleton checks
# The skill starts with its title line (the first line, not any line).
[ "$(head -n1 "$SKILL")" = "# nf-core/rnasplice Report Skill" ] || { echo "FAIL: missing required text: # nf-core/rnasplice Report Skill"; fail=1; }
# Command-line flags (--name): prose outside ``` fences never has one. Inside a fence a flag is allowed only in a command segment that
# starts with singularity exec, sbatch, Rscript or srun (also as the capture jid=$(sbatch ...)), in a line that continues such a
# command (the previous line ends with a backslash), and on #SBATCH lines. A line is split on ; && || and | into segments and
# the rule is applied to each segment, so a flag after the first command of a line does not pass.
# A fence is closed by a line of the same number of backticks as the one that opened it (four-backtick Rmd fences hold ``` lines).
if awk '
  function cmd(s) { return s ~ /^(singularity exec|sbatch|Rscript|srun)([[:space:]]|$)/ }
  function segok(s,   c) { c = s; sub(/^[A-Za-z_][A-Za-z0-9_${}]*=\\?\$\(/, "", c); return cmd(s) || cmd(c) }
  { t = $0; sub(/^[[:space:]]+/, "", t); n = 0; while (substr(t, n + 1, 1) == "`") n++ }
  !infence && n >= 3 { infence = 1; flen = n; allowed = 0; next }
  infence && n == flen && n == length(t) { infence = 0; allowed = 0; next }
  !infence { if ($0 ~ /(^|[[:space:]`(])--[A-Za-z]/) found = 1; next }
  t ~ /^#SBATCH / { allowed = ($0 ~ /\\$/); next }
  { ns = split(t, seg, /&&|\|\||;|\|/); lastok = 0
    for (i = 1; i <= ns; i++) {
      s = seg[i]; sub(/^[[:space:]]+/, "", s)
      lastok = segok(s) || (i == 1 && allowed)
      if (s ~ /(^|[[:space:]`(])--[A-Za-z]/ && !lastok) found = 1
    }
    allowed = (lastok && ($0 ~ /\\$/)) }
  END { exit !found }' "$SKILL"; then
  echo "FAIL: a --flag appears in prose or outside the allowed commands (singularity exec, sbatch, Rscript, srun, #SBATCH lines)"; fail=1
fi
# No process words in the skill or the README. "Task [0-9]" is case-sensitive (a lower-case "task 3" is ordinary prose).
for re in 'Task [0-9]' 'ledger' 'fix round' 'reviewer'; do
  if [ "$re" = 'Task [0-9]' ]; then gi=-E; else gi=-iE; fi
  ! grep -q $gi -- "$re" "$SKILL" || { echo "FAIL: forbidden process word in the skill: $re"; fail=1; }
  ! grep -q $gi -- "$re" "$README" 2>/dev/null || { echo "FAIL: forbidden process word in the README: $re"; fail=1; }
done
# --- end skeleton checks

# --- wizard checks (Steps 0-6, the inputs)
# Every step is numbered and states a default.
for s in 0 1 2 3 4 5 6; do
  need "## Step $s — "
  awk -v s="$s" '/^## Step / { insec = ($3 == s); next } insec && /[Dd]efault/ { f = 1 } END { exit !f }' "$SKILL" || { echo "FAIL: Step $s states no default"; fail=1; }
done
need 'ask one question at a time, state the default with each question'
need 'Shell variables and functions do not persist between Bash tool calls'
need 'one Bash call per procedure'
need 'R, Nextflow, Java and every other heavy tool run only on compute nodes'
need 'never on a login node'
need 'Default SLURM requests are at most 64 GB and 4 hours'
need 'Every `module add` is guarded with `|| exit 1`'
need '`{WD_NAME}` = basename of `{CWD}`'
need '`{TODAY_ISO}` = today'\''s date as `YYYY-MM-DD`'
need '`{REPORT_DIR}` = `{CWD}/results/{TODAY_ISO}_{WD_NAME}_report`'
need 'already exists (a previous report?)'
need 'Reuse it: files of the same name are replaced and the others are kept'
need 'A previous report is never overwritten silently'
need 'Do not pre-fill any address.'
need 'Store `{EMAIL}`'
need 'Store `{SIF}`'
need 'The image is pulled on a compute node, never on the login node'
need '#SBATCH -p bcc'
need 'module add singularity/3.10.4 || exit 1'
need '"DEXSeq", "edgeR", "stageR", "ggplot2", "dplyr", "tidyr", "readr", "knitr", "rmarkdown"'
need '"ReactomePA", "KEGG.db", "reactome.db"'
need 'the image must be rebuilt by them'
need 'packages are never installed at run time'
need 'KEGG/Reactome require an image rebuild'
need 'store `{ENRICH_DBS}` = `GO KEGG Reactome`'
need 'is the default for `{OUTDIR_RNASPLICE}`'
need 'has more than a header line'
need 'an empty file or a file with only a header is not found'
need 'star/rmats/{CONTRAST}/{SE,A5SS,A3SS,MXE,RI}.MATS.{JC,JCEC}.txt'
need 'star/rmats/{CONTRAST}_paired/'
need 'RMATS_DIR_PAIRED=star/rmats/${c}_paired'
need 'salmon/suppa/diffsplice/per_local_event/local_{T}-{C}.dpsi'
need 'star/dexseq_exon/results/DEXSeqResults.{CONTRAST}.csv'
need 'star/edger/contrast_{CONTRAST}.usage.exon.csv'
need 'salmon/dexseq_dtu/results/dexseq/DEXSeqResults.{T}-{C}.tsv'
need 'never by splitting file names'
need 'The rMATS tables have the column `ID` twice'
need 'has no `salmon/` directory, so SUPPA2 and DTU are simply not found'
need 'Store the modules as `{MODULES}`'
need 'echo "MODULES=${mods# }"'
need 'sample,condition,genome_bam,strandedness,single_end'
need 'the FASTQ or BAM files the sheet points to are never opened'
need 'Contrast `{CONTRAST}`: is `{T}` the treatment and `{C}` the control'
need 'This fixes every sign'
need 'Store `{CONTRASTS}` as the confirmed list'
need '**Names built from the contrasts file.**'
need 'contrasts are never mixed'
need 'Store the GTF path as `{GTF}`'
need 'ENSMUSG'
need 'Store `{ORGANISM}` ='
need 'Store `{PAIRED_RMATS}` and `{PAIRED_SUPPA}`'
need '`rmats_paired_stats` (rMATS) and `diffsplice_paired` (SUPPA2)'
need 'never claims unpaired statistics for a paired run'
# Names of files and groups are built from the contrasts file, never recovered by splitting on - or _.
forbid_re 'awk -F[-_]|cut -d[-_]|IFS=[-_]' 'file or group names split on - or _'
# R is started only through singularity exec (sbatch), and packages are never installed at run time.
forbid_re '^(Rscript|R -e|R -f|R --)' 'R started outside singularity exec'
forbid_re 'install\.packages|BiocManager|install_github|conda install|pip install' 'a package is installed at run time'
# The FASTQ or BAM files of the sample sheet are never opened.
forbid_re 'zcat[^|]*(fastq|fq)' 'a FASTQ file is opened'
# Every module add is guarded.
if grep 'module add' "$SKILL" | grep -qv '|| exit 1'; then echo "FAIL: a module add without || exit 1"; fail=1; fi
# --- end wizard checks
# --- wizard checks (Steps 7-10, thresholds, menu, SLURM)
for s in 7 8 9 10; do
  need "## Step $s — "
  awk -v s="$s" '/^## Step / { insec = ($3 == s); next } insec && /[Dd]efault/ { f = 1 } END { exit !f }' "$SKILL" || { echo "FAIL: Step $s states no default"; fail=1; }
done
need 'Default: every job that applies.'
need '(FDR) cutoff for a significant event, gene or transcript. Default: 0.05.'
need 'Default: 0.1. It is not applied to exon usage or transcript usage'
need '`{RMATS_TYPE}`: which rMATS count type to read, `JC` (junction counts only) or `JCEC`'
need 'Default: `JC`.'
need 'minimum read support of an rMATS event: for each group, the mean over its samples of `IJC` plus `SJC`, required in both groups'
need 'Default: 10. Events below it are flagged (`low_support`), kept in the tables and excluded from the significance counts of the reports, never deleted'
need 'the rMATS FDR is not recomputed (it stays as rMATS wrote it)'
need 'On very shallow data (like the nf-core test data) most events fall below 10; the user may lower it.'
forbid 'dropped from the tests'
need 'applies the Benjamini-Hochberg correction to it and says so'
need 'a space-separated subset of `events cross topgene enrich`'
need 'two or more of `rmats`, `suppa`, `dexseq_exon`, `edger_exon` and `dtu` are in `{MODULES}` (each counts as a separate tool)'
need 'drawn from the pipeline'\''s own tables only, never from BAM files'
need 'Offered only when `{ORGANISM}` is `human` or `mouse`'
need '`events` and `topgene` are both passed to the module reports as parameters'
need 'Memory `{MEM}`: default 32G (at most 64G).'
need 'Time `{TIME}`: default 4:00:00 (at most 4:00:00).'
need 'the name used in the job names. Default: `{WD_NAME}`.'
need 'one node, one task and four CPUs (`-N 1 -n 1 -c 4`)'
need 'is written only when the user asks for it by name'
need '| `run_01.sbatch` | `01_import.Rmd` | always |'
need '| `run_02.sbatch` | `02_rmats.Rmd` | `rmats` in `{MODULES}` |'
need '| `run_03.sbatch` | `03_suppa.Rmd` | `suppa` in `{MODULES}` |'
need '| `run_04.sbatch` | `04_usage.Rmd` (DEXSeq and edgeR exon usage, and DTU) | any of `dexseq_exon`, `edger_exon`, `dtu` in `{MODULES}` |'
need '| `run_05.sbatch` | `05_cross_tool.Rmd` | `cross` in `{SECTIONS}` and two or more tools in `{MODULES}` |'
need '| `run_06.sbatch` | `06_enrichment.Rmd` | `enrich` in `{SECTIONS}` |'
need '| `run_07.sbatch` | `07_summary.Rmd` (always last) | always |'
need 'Jobs 02 to 06 start when job 01 has finished successfully (`dependency=afterok`, with the option that cancels them when job 01 fails)'
need 'job 07 starts when jobs 01 to 06 (those that exist) have ended, whatever their outcome (`dependency=afterany`)'
need '**SLURM scripts.**'
need 'case "$R" in */results/*) ;; *) fail "report directory $R is not under a results/ directory" ;; esac'
need 'is inside the pipeline output $OUT'
need 'are above the default limits (64G, 4:00:00), as requested'
need '#SBATCH -N 1'
need '#SBATCH -n 1'
need '#SBATCH -c 4'
need '#SBATCH --mem=$MEM'
need '#SBATCH -t $TIME'
need 'singularity exec --bind $BIND "$SIF" Rscript -e "rmarkdown::render('\''$R/$2'\'','
need "output_dir = '\$R', knit_root_dir = '\$R',"
need 'events = $EV, cross = $CR, topgene = $TG, enrich = $EN'
need 'if sec cross && [ $ntools -ge 2 ]; then JOBS="$JOBS 05:05_cross_tool.Rmd"; fi'
need 'if sec enrich; then JOBS="$JOBS 06:06_enrichment.Rmd"; fi'
need 'JOBS="$JOBS 07:07_summary.Rmd"'
need 'JOB$n=\$(sbatch --parsable --dependency=afterany:$ALLDEP run_$n.sbatch) || exit 1'
need 'JOB$n=\$(sbatch --parsable --dependency=afterok:\$JOB01 --kill-on-invalid-dep=yes run_$n.sbatch) || exit 1'
need '    outdir = '\''$OUT'\'', report_dir = '\''$R'\'','
need '    contrasts = '\''$CONTRASTS_CSV'\'', samplesheet = '\''$SAMPLESHEET'\'','
need '    modules = '\''$MODULES'\'', rmats_type = '\''$RMATS_TYPE'\'','
need '    fdr = $FDR, dpsi = $DPSI, min_reads = $MIN_READS,'
need '    organism = '\''$ORGANISM'\'', gtf = '\''$GTF'\'','
need '    enrich_dbs = '\''$ENRICH_DBS'\'','
need '    paired_rmats = $PR, paired_suppa = $PS))'
need '`run_all.sh`, which only calls sbatch'
need 'echo "SUBMIT=bash $R/run_all.sh"'
need '**Final instructions.**'
need 'Run this one command:'
need 'the output of the run is read-only for the skill and is never written to'
# The final instructions print exactly one submit command (a bash or sbatch line naming a run script of the report directory).
n_submit=$(grep -cE '^[[:space:]]*(bash|sbatch)[[:space:]].*\{REPORT_DIR\}.*run_' "$SKILL")
[ "$n_submit" -eq 1 ] || { echo "FAIL: the final instructions must print exactly one submit command (bash {REPORT_DIR}/run_all.sh), found $n_submit"; fail=1; }
# Defaults of the resources are at most 64G and 4h; so is every literal request of a script written in the skill.
mem_def=$(sed -n 's/.*Memory `{MEM}`: default \([0-9]*\)G.*/\1/p' "$SKILL" | head -n1)
time_def=$(sed -n 's/.*Time `{TIME}`: default \([0-9]*\):\([0-9][0-9]\):\([0-9][0-9]\).*/\1\2\3/p' "$SKILL" | head -n1)
[ -n "$mem_def" ] && [ "$mem_def" -le 64 ] || { echo "FAIL: default memory is missing or above 64G: '$mem_def'"; fail=1; }
[ -n "$time_def" ] && [ "$time_def" -le 40000 ] || { echo "FAIL: default time is missing or above 4:00:00: '$time_def'"; fail=1; }
if grep -oE '^#SBATCH --mem=[0-9]+[GM]' "$SKILL" | sed 's/.*=//' | awk '/G$/ { v = $0; sub(/G/, "", v); if (v + 0 > 64) bad = 1 } /M$/ { v = $0; sub(/M/, "", v); if (v + 0 > 65536) bad = 1 } END { exit !bad }'; then
  echo "FAIL: a literal memory request above 64G in a script of the skill"; fail=1
fi
if grep -oE '^#SBATCH -t [0-9]+(-[0-9]+)?:[0-9]+:[0-9]+' "$SKILL" | sed 's/^#SBATCH -t //' | awk -F'[-:]' '{ if (NF == 4) s = ($1 * 24 + $2) * 3600 + $3 * 60 + $4; else s = $1 * 3600 + $2 * 60 + $3; if (s > 14400) bad = 1 } END { exit !bad }'; then
  echo "FAIL: a literal time request above 4:00:00 in a script of the skill"; fail=1
fi
# --- end wizard checks (Steps 7-10)

# --- Rmd 01 (import and harmonise)
need '**Writing the Rmd files.**'
need 'Write that block with the Write tool to `{REPORT_DIR}/<name>.Rmd` byte for byte, without the two fence lines'
need '**File `01_import.Rmd`.**'
# the five sign rules of harmonise_effect, exactly (treatment minus control)
need '    rmats = value,          # IncLevelDifference: b1 = treatment, already treatment - control'
need '    suppa = -value,         # SUPPA2 dPSI is control - treatment'
need '    exon_dexseq = -value,   # log2fold_{CONTROL}_{TREATMENT}'
need '    dtu = -value,           # log2fold_{CONTROL}_{TREATMENT}'
need '    exon_edger = value,     # logFC is treatment - control'
need '  ifelse(is.na(effect) | effect == 0, "none", ifelse(effect > 0, "treatment_up", "control_up"))'
# a sign flip happens only in those lines: -value nowhere else, and no bare -effect (also - x$effect) anywhere
if grep -E '(^|[^A-Za-z0-9_])-[[:space:]]*value([^A-Za-z0-9_]|$)' "$SKILL" | grep -vxF -e '    suppa = -value,         # SUPPA2 dPSI is control - treatment' \
     -e '    exon_dexseq = -value,   # log2fold_{CONTROL}_{TREATMENT}' -e '    dtu = -value,           # log2fold_{CONTROL}_{TREATMENT}' | grep -q .; then
  echo "FAIL: a sign flip (-value) outside harmonise_effect"; fail=1
fi
forbid_re '(^|[^A-Za-z0-9_])-[[:space:]]*([A-Za-z_.][A-Za-z0-9_.]*\$)?effect([^A-Za-z0-9_]|$)' 'a bare -effect (a sign flip outside harmonise_effect)'
need 'This report applies the Benjamini-Hochberg correction itself, separately for each contrast and for each file (local events, isoforms)'
need 'harmonised_table(contrast, "suppa", gene, NA_character_, ftype, event, dpsi, bh_adjust(p), NA_real_)'
need 'qc$suppa_untested[["local"]] <- sum(is.na(parts$local$effect))'
need 'qc$suppa_untested[["isoform"]] <- sum(is.na(parts$isoform$effect))'
need 'strip_quotes <- function(x) gsub("^\"|\"$", "", x)'
need 'names(tab) <- dedupe_names(names(tab))'
need '# support of a group = mean over its samples of (IJC + SJC); n_support = min(support of group 1, support of group 2).'
# low-support rMATS events are flagged, never deleted, and excluded from the significant counts
need '  low <- is.na(ns) | ns < min_reads'
need '"padj", "direction", "n_support", "low_support")'
need '    !tab$low_support'
grep -qxF '      out$rmats <- tab' "$SKILL" || { echo "FAIL: the rMATS table must be kept whole (out\$rmats <- tab): flagged rows are never deleted"; fail=1; }
need 'are flagged in the column `low_support`: they stay in the table, are excluded from the significant counts and are counted below; nothing is deleted, and the rMATS FDR is not recomputed.'
# contract of a missing module file
need 'no RDS is written for that tool and contrast and the file is listed in `qc$files_missing`; Rmds 02 to 06 check that each RDS exists before they read it'
need 'effect_scale, padj, direction, n_support, low_support`'
need 'if (length(miss)) stop("file ", path, ": column not found: ", paste(miss, collapse = ", "),'
# no network and no installation inside the Rmd (the image is used as it is)
forbid_re 'download\.file|url\(|install\.packages' 'network access or package installation in an Rmd'
# --- end Rmd 01

# --- Rmd 02 (rMATS) and Rmd 03 (SUPPA2)
need '**File `02_rmats.Rmd`.**'
need '**File `03_suppa.Rmd`.**'
need 'DIRECTION_WORDING <- c(treatment_up = "more included in treatment", control_up = "more skipped in treatment")'
need 'none_message <- function(fdr, dpsi) paste0("none at FDR ", format(fdr), ", |dPSI| ", format(dpsi))'
need '  !is.na(tab$padj) & tab$padj < fdr & !is.na(tab$effect) & abs(tab$effect) >= dpsi & !tab$low_support'
need '  !is.na(tab$padj) & tab$padj < fdr & !is.na(tab$effect) & abs(tab$effect) >= dpsi'
need '  cn <- split_counts(tab, params$fdr, params$dpsi, RMATS_TYPES)'
need '  cn <- split_counts(tab, params$fdr, params$dpsi, types)'
need '  g1 <- group_samples(ss, treatment, bam_treatment)'
need '  g2 <- group_samples(ss, control, bam_control)'
need '  s1 <- g1$samples'
need '  s2 <- g2$samples'
need '  f <- file.path(outdir, "star", "rmats", "bamlist", paste0(condition, "_bamlist.txt"))'
need '  list(samples = bamlist, from = "rMATS sample list", warning = NULL)'
need 'WARNING: no rMATS sample list found for '
need 'psi_legend <- function(treatment, control) c(paste0("treatment (", treatment, ")"), paste0("control (", control, ")"))'
need '  ev <- sub("^[^;]*;", "", id)'
need 'The adjusted p-value used here is the Benjamini-Hochberg correction that `01_import.Rmd` applied'
need '`rds/{CONTRAST}.rmats.rds`'
need '`rds/{CONTRAST}.suppa.rds`'
need '`tables/{CONTRAST}.rmats.{type}.significant.csv`'
need '`figures/{CONTRAST}.rmats.volcano.{type}.png`'
need '`figures/{CONTRAST}.rmats.psi_top.png`'
need '`tables/{CONTRAST}.suppa.{type}.significant.csv`'
need '`figures/{CONTRAST}.suppa.volcano.{type}.png`'
rmd_block() { awk -v a="**File \`$1\`.**" '$0 == a { seen = 1; next } seen && !inb && $0 == "````rmd" { inb = 1; next } inb && $0 == "````" { exit } inb { print }' "$SKILL"; }
for n in 02_rmats.Rmd 03_suppa.Rmd; do
  b=$(rmd_block "$n")
  [ -n "$b" ] || { echo "FAIL: $n: no rmd block under its heading"; fail=1; continue; }
  # the sign sentence stands in the text before the first code chunk with functions
  ls=$(printf '%s\n' "$b" | grep -nF 'Every effect is stated as **treatment minus control**: a positive effect is higher in the treatment group' | head -n1 | cut -d: -f1)
  lf=$(printf '%s\n' "$b" | grep -n '^```{r functions' | head -n1 | cut -d: -f1)
  { [ -n "$ls" ] && [ -n "$lf" ] && [ "$ls" -lt "$lf" ]; } || { echo "FAIL: $n: the sign sentence (Every effect is stated as treatment minus control) must stand before the functions chunk"; fail=1; }
  # the parameters the report switches on are read (the 18 params of the contract; no other name)
  case "$n" in 02_rmats.Rmd) extra="samplesheet outdir min_reads rmats_type paired_rmats" ;; *) extra="paired_suppa" ;; esac
  for v in fdr dpsi events topgene report_dir contrasts modules $extra; do
    printf '%s\n' "$b" | grep -qF "params\$$v" || { echo "FAIL: $n: the parameter params\$$v is not used"; fail=1; }
  done
  for v in $(printf '%s\n' "$b" | grep -o 'params\$[A-Za-z_]*' | sort -u | sed 's/params\$//'); do
    case " outdir report_dir contrasts samplesheet modules fdr dpsi rmats_type min_reads events cross topgene enrich organism enrich_dbs gtf paired_rmats paired_suppa " in
      *" $v "*) ;; *) echo "FAIL: $n: params\$$v is not one of the 18 parameters of the contract"; fail=1 ;; esac
  done
  # the thresholds come from the parameters: no literal 0.05 or 0.1 in the body of the Rmd (the front matter holds the defaults)
  if printf '%s\n' "$b" | awk 'f >= 2 { print } /^---$/ { f++ }' | grep -qE '(^|[^0-9.])0\.(05|1)([^0-9]|$)'; then
    echo "FAIL: $n: a literal threshold (0.05 or 0.1) in the body of the Rmd instead of params\$fdr and params\$dpsi"; fail=1
  fi
  printf '%s\n' "$b" | grep -qF 'paste0("none at FDR "' || { echo "FAIL: $n: the empty-set message none at FDR is missing"; fail=1; }
  printf '%s\n' "$b" | grep -qF 'three nearest events by adjusted p-value' || { echo "FAIL: $n: the three nearest events are not printed under the empty-set message"; fail=1; }
  printf '%s\n' "$b" | grep -qF 'sig_flag(st, params$fdr, params$dpsi)' || { echo "FAIL: $n: the significant events are not selected with params\$fdr and params\$dpsi"; fail=1; }
done
b=$(rmd_block 02_rmats.Rmd)
printf '%s\n' "$b" | grep -qE '^\s+s1 <- g1\$samples$' || { echo "FAIL: 02_rmats.Rmd: the PSI plot must take the sample order from rMATS's own sample lists (s1 <- g1\$samples, group 1 = the treatment)"; fail=1; }
printf '%s\n' "$b" | grep -qE '^\s+s2 <- g2\$samples$' || { echo "FAIL: 02_rmats.Rmd: group 2 (the control) must take the sample order of its own sample list (s2 <- g2\$samples)"; fail=1; }
printf '%s\n' "$b" | grep -qF 'psi_long(inc, ss, Tr, Co, bam_t, bam_c)' || { echo "FAIL: 02_rmats.Rmd: the PSI plot is not built from the rMATS sample lists (psi_long(inc, ss, Tr, Co, bam_t, bam_c))"; fail=1; }
printf '%s\n' "$b" | grep -qF 'bam_t <- read_bamlist(params$outdir, Tr); bam_c <- read_bamlist(params$outdir, Co)' || { echo "FAIL: 02_rmats.Rmd: the sample lists of the run are not read (read_bamlist)"; fail=1; }
if printf '%s\n' "$b" | grep -qF 'order of the sample sheet.'; then echo "FAIL: 02_rmats.Rmd: the PSI sample order must not be described as the order of the sample sheet"; fail=1; fi
printf '%s\n' "$b" | grep -qF 'low_support_report(tab, params$fdr, params$dpsi)' || { echo "FAIL: 02_rmats.Rmd: the low_support count is not reported"; fail=1; }
# --- end Rmd 02 and 03

# tests/run_all_tests.sh names every script of tests/ in exactly one list (TESTS, FULL_TESTS or NOT_TESTS), runs each under
# env -i, and says when the proof-tool self-test was skipped. CHECK_RUN_ALL overrides the file (used by prove_mutations.sh).
RUN_ALL=${CHECK_RUN_ALL:-$HERE/run_all_tests.sh}
if [ -s "$RUN_ALL" ]; then
  ra_tests=$(sed -n 's/^TESTS="\(.*\)"$/\1/p' "$RUN_ALL"); ra_full=$(sed -n 's/^FULL_TESTS="\(.*\)"$/\1/p' "$RUN_ALL"); ra_not=$(sed -n 's/^NOT_TESTS="\(.*\)"$/\1/p' "$RUN_ALL")
  for f in "$HERE"/*.sh; do
    b=${f##*/}; c=0
    for l in "$ra_tests" "$ra_full" "$ra_not"; do case " $l " in *" $b "*) c=$((c + 1)) ;; esac; done
    [ $c -eq 1 ] || { echo "FAIL: tests/$b must be named in exactly one list of run_all_tests.sh (TESTS, FULL_TESTS or NOT_TESTS), found in $c"; fail=1; }
  done
  [ "$ra_not" = "cut_block.sh test_env.sh prove_red.sh run_all_tests.sh run_r_tests.sh" ] \
    || { echo "FAIL: run_all_tests.sh NOT_TESTS must be exactly the helpers (cut_block.sh test_env.sh prove_red.sh run_all_tests.sh run_r_tests.sh), found: $ra_not"; fail=1; }
  case " $ra_full " in *" test_proof_tools.sh "*) ;; *) echo "FAIL: run_all_tests.sh FULL_TESTS must include test_proof_tools.sh, found: $ra_full"; fail=1 ;; esac
  # the R tests and the rtest mutation rows run only in the full mode (through sbatch); the default suite stays static
  for t in test_r.sh test_rtest_runner.sh; do
    case " $ra_full " in *" $t "*) ;; *) echo "FAIL: run_all_tests.sh FULL_TESTS must include $t, found: $ra_full"; fail=1 ;; esac
  done
  grep -qxF '  run prove_mutations.sh "$SK" rtest' "$RUN_ALL" || { echo "FAIL: run_all_tests.sh --full must run prove_mutations.sh in mode rtest"; fail=1; }
  case " $ra_tests " in *" test_proof_tools.sh "*) echo "FAIL: run_all_tests.sh TESTS must not include test_proof_tools.sh (it belongs to FULL_TESTS)"; fail=1 ;; esac
  grep -qF 'env -i HOME="${HOME:-/nonexistent}" PATH="$P" /bin/bash --noprofile --norc "$HERE/$1"' "$RUN_ALL" \
    || { echo "FAIL: run_all_tests.sh must run every test under env -i with an explicit PATH"; fail=1; }
  grep -qF 'test_proof_tools.sh was SKIPPED' "$RUN_ALL" || { echo "FAIL: run_all_tests.sh must say when the proof-tool self-test was skipped"; fail=1; }
  grep -qF "1.8.5 | sort -V | head -n1)\" != 1.8.5 ]; then" "$RUN_ALL" || { echo "FAIL: run_all_tests.sh must check for git 1.8.5 or newer before the proof-tool self-test"; fail=1; }
  case " $ra_tests " in *" test_run_all.sh "*) ;; *) echo "FAIL: run_all_tests.sh TESTS must run test_run_all.sh (the test of the runner itself)"; fail=1 ;; esac
else
  echo "FAIL: tests/run_all_tests.sh is missing or empty"; fail=1
fi
[ $fail -eq 0 ] && echo "PASS" || exit 1
