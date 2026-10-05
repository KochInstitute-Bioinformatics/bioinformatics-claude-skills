# nf-core/rnasplice Report Skill

You are helping the user turn the output of an nf-core/rnasplice run (set up with `/nfcore-rnasplice-setup`) into R Markdown reports and SLURM scripts on an HPC cluster using SLURM and Singularity. The run's tables come from rMATS, SUPPA2, DEXSeq and edgeR exon usage, and DEXSeq transcript usage; the report adapts to the modules that produced output. Every direction in every table and plot is stated as treatment minus control, and the treatment and control of each contrast are fixed once, in Step 4. This skill reads the pipeline's output and never writes into it. Walk through the steps below in order, ask one question at a time, state the default with each question, and perform automated steps silently.

Rules that hold for every step:
- R, Nextflow, Java and every other heavy tool run only on compute nodes, through sbatch (or salloc and srun for debugging), never on a login node. The login node is used only for light commands (`pwd`, `ls`, `test`, `head`, `grep`, `awk`, `sed`, `wc`, `mkdir`, `squeue`).
- Default SLURM requests are at most 64 GB and 4 hours.
- The files of the run are read-only. Everything the skill writes goes to the report directory `{REPORT_DIR}` (Step 0), under `results/`.
- **Bash tool calls.** Shell variables and functions do not persist between Bash tool calls: one Bash call per procedure, and every block below is pasted whole into one call, with its placeholders (`{...}`) replaced by the literal values from the earlier steps. A value needed in a later call is written out literally there, never taken from a shell variable of an earlier call.
- Every `module add` is guarded with `|| exit 1` (under Lmod a failed `module add` silently breaks every later one).

---

## Step 0 — Working directory and constants

Run `pwd` to record the current working directory (`{CWD}`). The sample sheet, the contrasts sheet and the params file of the run are expected here (the files `/nfcore-rnasplice-setup` wrote). Tell the user: "Working directory: `{CWD}`"

Also derive:
- `{WD_NAME}` = basename of `{CWD}`
- `{TODAY_ISO}` = today's date as `YYYY-MM-DD`
- `{REPORT_DIR}` = `{CWD}/results/{TODAY_ISO}_{WD_NAME}_report` (the reports, tables, figures and scripts go here; the name differs from the pipeline's own `results/{TODAY_ISO}_{WD_NAME}` so that the report never lands inside the pipeline output)

**Report directory.** Check whether `{REPORT_DIR}` exists, in one Bash call:
```bash
if [ -e "{REPORT_DIR}" ]; then echo "EXISTS {REPORT_DIR}"; ls -la "{REPORT_DIR}" | head -n 20; else echo "NEW {REPORT_DIR}"; fi
```
When it exists, ask before using it (numbered; default 1): "`{REPORT_DIR}` already exists (a previous report?). 1. Use a new directory with a numbered suffix, `{REPORT_DIR}_2` (default) · 2. Reuse it: files of the same name are replaced and the others are kept · 3. Give another path." A previous report is never overwritten silently: nothing is written into an existing directory before the user has answered. The answer sets `{REPORT_DIR}`; Step 3 checks that it is not inside the pipeline output directory.

Tell the user: "Working directory: `{CWD}`; report directory: `{REPORT_DIR}`".

---

## Step 1 — Email address

Ask the user: "What email address should SLURM use for notifications (completion and failure of the report jobs)? Default: none, no mail is sent."

Do not pre-fill any address. Wait for the answer before Step 2. Store `{EMAIL}`. When the user wants no mail, the two `#SBATCH` mail lines are left out of every script.

---

## Step 2 — R container image

The reports run R inside the Singularity image `{SIF}`. Ask (default 1): "Which R image? 1. `/net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif` (default) · 2. another image path." Store `{SIF}`.

**Image present.** Check, in one Bash call: `ls -la "{SIF}"`.

**Image missing.** The image is pulled on a compute node, never on the login node. Write and submit this job (replace `{SIF}`, `{REPORT_DIR}` and `{EMAIL}`; the job stops with an error when the module or the pull fails), tell the user the job id, wait until the job has left the queue (`squeue -u $USER`), then continue:
```bash
mkdir -p "{REPORT_DIR}/image_check"
cat > "{REPORT_DIR}/image_check/pull_image.sbatch" <<'EOF'
#!/bin/bash
#SBATCH -p bcc
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -t 1:00:00
#SBATCH --mail-user={EMAIL}
#SBATCH --mail-type=END,FAIL
#SBATCH -o {REPORT_DIR}/image_check/pull_image.%j.log
module add singularity/3.10.4 || exit 1
cd "$(dirname "{SIF}")" || exit 1
singularity pull "$(basename "{SIF}")" docker://yannvrb56/bulkrnaseq:latest || exit 1
ls -la "{SIF}"
EOF
sbatch "{REPORT_DIR}/image_check/pull_image.sbatch"
```

**Package check.** The image must hold the packages the reports load. The check runs on a compute node (R is never started on the login node). Write and submit it, wait for the job to leave the queue, then read its log:
```bash
mkdir -p "{REPORT_DIR}/image_check"
cat > "{REPORT_DIR}/image_check/check_packages.R" <<'EOF'
need <- c("DEXSeq", "edgeR", "stageR", "ggplot2", "dplyr", "tidyr", "readr", "knitr", "rmarkdown")
opt  <- c("UpSetR", "ComplexUpset", "clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db", "enrichplot", "ggrepel", "scales",
          "jsonlite", "yaml", "msigdbr", "KEGGREST", "GO.db", "AnnotationDbi", "ReactomePA", "KEGG.db", "reactome.db")
for (p in need) cat(sprintf("MANDATORY %s %s\n", p, requireNamespace(p, quietly = TRUE)))
for (p in opt)  cat(sprintf("OPTIONAL %s %s\n", p, requireNamespace(p, quietly = TRUE)))
EOF
cat > "{REPORT_DIR}/image_check/check_packages.sbatch" <<'EOF'
#!/bin/bash
#SBATCH -p bcc
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 0:20:00
#SBATCH --mail-user={EMAIL}
#SBATCH --mail-type=END,FAIL
#SBATCH -o {REPORT_DIR}/image_check/check_packages.%j.log
module add singularity/3.10.4 || exit 1
singularity exec --bind /net/bmc-lab3 "{SIF}" Rscript "{REPORT_DIR}/image_check/check_packages.R"
EOF
sbatch "{REPORT_DIR}/image_check/check_packages.sbatch"
```
Read the lines `MANDATORY <package> TRUE|FALSE` and `OPTIONAL <package> TRUE|FALSE` of the log:
- A mandatory package is FALSE: stop. Tell the user which packages are missing and that the image must be rebuilt by them; packages are never installed at run time (not in a job, not by the wizard), because the reports must run from the image as it is.
- ReactomePA, KEGG.db or reactome.db is FALSE (and no mandatory package is): enrichment is offered with GO only; tell the user "KEGG/Reactome require an image rebuild". Store `{ENRICH_DBS}` = `GO`.
- All three are TRUE: store `{ENRICH_DBS}` = `GO KEGG Reactome`.
- Also note which of UpSetR, ComplexUpset and the organism packages (org.Hs.eg.db, org.Mm.eg.db) are FALSE; the cross-tool comparison and the enrichment for that organism are not offered without them.

Default: the shared image, which holds all nine mandatory packages. Tell the user the result in one line.

---

## Step 3 — Locate the run and detect the modules

Ask (default: the `outdir` of the params file): "Where is the output directory of the rnasplice run?" First look for the params file, in one Bash call:
```bash
ls "{CWD}"/*_params.yaml 2>/dev/null
grep -E '^(input|contrasts|source|outdir|gtf|genome|rmats|suppa|dexseq_exon|edger_exon|dexseq_dtu|rmats_paired_stats|diffsplice_paired):' "{CWD}"/*_params.yaml 2>/dev/null
```
One params file found: its `outdir` (relative to `{CWD}` when it is not absolute) is the default for `{OUTDIR_RNASPLICE}`, its `contrasts` file is `{CONTRASTS_CSV}` and its `input` file is `{SAMPLESHEET}` (relative paths are relative to `{CWD}`; write the absolute paths). Several or none: ask which params file, or for the paths of the three files. The output directory must exist, and `{REPORT_DIR}` must not lie inside it.

**Module detection.** A module counts as found only when its main result file exists and has more than a header line; an empty file or a file with only a header is not found. The files are named from the contrasts file (treatment `{T}`, control `{C}`, contrast name `{CONTRAST}`), never by splitting file names: rMATS `star/rmats/{CONTRAST}/{SE,A5SS,A3SS,MXE,RI}.MATS.{JC,JCEC}.txt`, or `star/rmats/{CONTRAST}_paired/` with the same file names when the run used `rmats_paired_stats: true` (either directory counts); SUPPA2 `salmon/suppa/diffsplice/per_local_event/local_{T}-{C}.dpsi`; DEXSeq exon `star/dexseq_exon/results/DEXSeqResults.{CONTRAST}.csv`; edgeR exon `star/edger/contrast_{CONTRAST}.usage.exon.csv`; DTU `salmon/dexseq_dtu/results/dexseq/DEXSeqResults.{T}-{C}.tsv`. SUPPA2 and DTU file names use `{T}-{C}` (treatment, hyphen, control); the other three use the contrast name. A BAM-input run has no `salmon/` directory, so SUPPA2 and DTU are simply not found. The rMATS tables have the column `ID` twice (first column, and again after the coordinate columns); detection only counts lines and never reads columns, so this does not matter here, but every reader of these tables must account for it. Run, in one Bash call (replace `{OUTDIR_RNASPLICE}` and `{CONTRASTS_CSV}` with the literal paths):
```bash
OUT="{OUTDIR_RNASPLICE}"
CSV="{CONTRASTS_CSV}"
[ -d "$OUT" ] || { echo "ERROR: output directory $OUT does not exist"; exit 1; }
[ -s "$CSV" ] || { echo "ERROR: contrasts file $CSV is missing or empty"; exit 1; }
[ "$(head -n1 "$CSV" | tr -d '\r')" = "contrast,treatment,control" ] || { echo "ERROR: header of $CSV is not contrast,treatment,control"; exit 1; }
# a file counts only when it exists and has more than a header line (two or more non-blank lines)
has_rows() { [ -f "$1" ] && awk 'NF { n++; if (n > 1) exit } END { exit !(n > 1) }' "$1"; }
found=""
add() { case " $found " in *" $1 "*) ;; *) found="$found $1" ;; esac; }
chk() { if has_rows "$3"; then echo "found $1 $2 $3"; add "$1"; else echo "missing $1 $2 $3"; fi; }
n=0
while IFS=, read -r c t k _ || [ -n "$c" ]; do
  k=${k%$'\r'}; n=$((n + 1)); [ $n -eq 1 ] && continue; [ -n "$c" ] || continue
  hit=""
  for d in "$c" "${c}_paired"; do for ty in SE A5SS A3SS MXE RI; do for m in JC JCEC; do has_rows "$OUT/star/rmats/$d/$ty.MATS.$m.txt" && hit="$hit $d/$ty.$m"; done; done; done
  if [ -n "$hit" ]; then echo "found rmats $c$hit"; add rmats; else echo "missing rmats $c $OUT/star/rmats/$c/ and $OUT/star/rmats/${c}_paired/ (no MATS table with data)"; fi
  chk suppa "$c" "$OUT/salmon/suppa/diffsplice/per_local_event/local_${t}-${k}.dpsi"
  chk dexseq_exon "$c" "$OUT/star/dexseq_exon/results/DEXSeqResults.${c}.csv"
  chk edger_exon "$c" "$OUT/star/edger/contrast_${c}.usage.exon.csv"
  chk dtu "$c" "$OUT/salmon/dexseq_dtu/results/dexseq/DEXSeqResults.${t}-${k}.tsv"
done < "$CSV"
[ -d "$OUT/salmon" ] || echo "note: no salmon directory in $OUT (BAM-input run, or the Salmon branch did not run): SUPPA2 and DTU are not found"
mods=""
for m in rmats suppa dexseq_exon edger_exon dtu; do case " $found " in *" $m "*) mods="$mods $m" ;; esac; done
[ -n "$mods" ] || echo "nothing found: compare the missing paths above with the output directory and the contrast names"
echo "MODULES=${mods# }"
```
The block prints one `found` or `missing` line per module and contrast (a missing line names the path it looked for) and one last line `MODULES=` followed by the modules found, in the order `rmats suppa dexseq_exon edger_exon dtu`. Store the modules as `{MODULES}`. Show the user what was found and what is missing, and which sections of the report that removes: a module that did not run is not offered later. Cross-check against the params file: a module the params switched on but that was not found is reported as a problem of the run (name its missing path), and a module found although the params had it off is mentioned. If `{MODULES}` is empty, stop and say so.

Default: all modules found are reported on. Ask (default: yes): "Report on these modules: `{MODULES}`? You can drop modules from the list."

---

## Step 4 — Sample sheet and contrasts (this is where every sign is fixed)

Read the pipeline's own sample sheet `{SAMPLESHEET}` and contrasts file `{CONTRASTS_CSV}` and print them. This reads the two sheets only; the FASTQ or BAM files the sheet points to are never opened. The sheet has either the FASTQ schema (`sample,fastq_1,fastq_2,strandedness,condition`) or the genome BAM schema (`sample,condition,genome_bam,strandedness,single_end`); the block finds the `condition` column by its header name. Run, in one Bash call:
```bash
SS="{SAMPLESHEET}"
CSV="{CONTRASTS_CSV}"
[ -s "$SS" ] || { echo "ERROR: sample sheet $SS is missing or empty"; exit 1; }
[ -s "$CSV" ] || { echo "ERROR: contrasts file $CSV is missing or empty"; exit 1; }
echo "== sample sheet $SS"; cat "$SS"
echo "== contrasts $CSV"; cat "$CSV"
echo "== samples per condition, and the groups of the contrasts file in the sheet"
awk -F, -v ss="$SS" '{ sub(/\r$/, "") }
  FNR == 1 { if (FILENAME == ss) for (i = 1; i <= NF; i++) if ($i == "condition") col = i; next }
  FILENAME == ss { if (col) { n[$col]++ } ; next }
  NF >= 3 { want[$2]; want[$3] }
  END { if (!col) print "ERROR: no condition column in the sample sheet"
        for (g in n) print "condition", g, n[g]
        for (g in want) print (g in n) ? "ok " g : "MISSING in the sample sheet: " g }' "$SS" "$CSV"
```
A group of the contrasts file that is absent from the `condition` column of the sheet is an error: stop and tell the user.

**Confirm treatment and control.** For each contrast ask (numbered; default 1): "Contrast `{CONTRAST}`: is `{T}` the treatment and `{C}` the control, so that every effect in every table and plot of the report is `{T}` minus `{C}` (positive means higher in `{T}`)? This fixes every sign: the report states the direction once from this answer, and every tool's own sign convention is converted to it. 1. Yes (default) · 2. No." The tools computed their tables with the treatment and control of the contrasts file, so the report cannot relabel them: on 2, stop and tell the user that the contrasts file must be corrected and the pipeline run again. Store `{CONTRASTS}` as the confirmed list of `contrast treatment control` triples.

**Names built from the contrasts file.** The expected file and column names of each contrast are built from its treatment and control; file names are never split to recover them (a group name may contain `-` and `_`). Print them, in one Bash call (replace `{CONTRASTS_CSV}`):
```bash
CSV="{CONTRASTS_CSV}"
[ -s "$CSV" ] || { echo "ERROR: contrasts file $CSV is missing or empty"; exit 1; }
[ "$(head -n1 "$CSV" | tr -d '\r')" = "contrast,treatment,control" ] || { echo "ERROR: header of $CSV is not contrast,treatment,control"; exit 1; }
n=0
while IFS=, read -r c t k _ || [ -n "$c" ]; do
  k=${k%$'\r'}; n=$((n + 1)); [ $n -eq 1 ] && continue; [ -n "$c" ] || continue
  echo "CONTRAST=$c"
  echo "TREATMENT=$t"
  echo "CONTROL=$k"
  echo "RMATS_DIR=star/rmats/$c"
  echo "RMATS_DIR_PAIRED=star/rmats/${c}_paired"
  echo "DEXSEQ_EXON_FILE=star/dexseq_exon/results/DEXSeqResults.$c.csv"
  echo "EDGER_EXON_FILE=star/edger/contrast_$c.usage.exon.csv"
  echo "SUPPA_FILE=salmon/suppa/diffsplice/per_local_event/local_${t}-${k}.dpsi"
  echo "DTU_FILE=salmon/dexseq_dtu/results/dexseq/DEXSeqResults.${t}-${k}.tsv"
  echo "DEXSEQ_COLUMN=log2fold_${k}_${t}"
  echo "SUPPA_COLUMN=local_${t}-local_${k}_dPSI"
done < "$CSV"
[ $n -gt 1 ] || echo "ERROR: no contrast row in $CSV"
```
The paths are relative to `{OUTDIR_RNASPLICE}`. Several contrasts give one block each; the report is per contrast, with the contrast named in every output file, and contrasts are never mixed.

---

## Step 5 — Annotation and organism

The GTF the run used gives the gene identifiers and the organism. Default: the `gtf` key of the params file. Run, in one Bash call (replace `{CWD}`):
```bash
grep -E '^(gtf|genome|fasta):' "{CWD}"/*_params.yaml 2>/dev/null
```
Store the GTF path as `{GTF}` (ask for it when the params file has none, for example a run that used `genome`). Then guess the organism from the identifiers in the GTF, in one Bash call (light work; replace `{GTF}`):
```bash
G="{GTF}"
[ -f "$G" ] || { echo "ERROR: $G is not a file"; exit 1; }
case "$G" in *.gz) rd=zcat ;; *) rd=cat ;; esac
$rd "$G" | head -n 5000 | grep -m1 -o 'ENS[A-Z]*G[0-9]\{6,\}' || echo "no Ensembl gene identifier found in the first 5000 lines"
```
`ENSG` followed by digits means human (`org.Hs.eg.db`), `ENSMUSG` mouse (`org.Mm.eg.db`). Ask (default: the guess): "Organism: 1. human (default when ENSG) · 2. mouse (default when ENSMUSG) · 3. other." Store `{ORGANISM}` = `human`, `mouse` or `other`. For `other`, enrichment is not offered (the image holds annotation packages for human and mouse only) and the reports say so.

---

## Step 6 — Pairing

The report states the pairing of the run and never claims unpaired statistics for a paired run. The pairing is read from the run, not asked for: the params file holds the switches `rmats_paired_stats` (rMATS) and `diffsplice_paired` (SUPPA2), and a paired design has a pairs file (header `sample,subject`) next to the sheets. Run, in one Bash call (replace `{CWD}`):
```bash
grep -E '^(rmats_paired_stats|diffsplice_paired):' "{CWD}"/*_params.yaml 2>/dev/null
ls "{CWD}"/*_pairs.csv 2>/dev/null && cat "{CWD}"/*_pairs.csv
```
Store `{PAIRED_RMATS}` and `{PAIRED_SUPPA}` (`true` or `false`) from the two switches. A switch that the params file does not name has the pinned pipeline's default: `rmats_paired_stats` false, `diffsplice_paired` true; say that the value is a default. Print the pairing to the user in one line per tool ("rMATS: unpaired statistics; SUPPA2: paired statistics" or the reverse) and, when a pairs file exists, its sample and subject columns. The two switches can differ; each report section states the pairing of its own tool. The DEXSeq and edgeR exon and DTU analyses use the sample sheet design of the run as it is; the report does not change it.

Default: the pairing as read. Ask (default: yes): "Is this pairing what you intended for the run? If not, the run has to be set up again; the report cannot change it."


---

## Step 7 — Thresholds

Ask one question per value, state the default with each, and store the answers:
- `{FDR}`: adjusted p-value (FDR) cutoff for a significant event, gene or transcript. Default: 0.05.
- `{DPSI}`: minimum absolute change in percent spliced in (dPSI) for an rMATS or SUPPA2 event to count as significant, in addition to `{FDR}`. Default: 0.1. It is not applied to exon usage or transcript usage, which are cut by `{FDR}` only; the report states this.
- `{RMATS_TYPE}`: which rMATS count type to read, `JC` (junction counts only) or `JCEC` (junction counts plus reads on the target exon). Default: `JC`. Asked only when `rmats` is in `{MODULES}`; a type whose files the run did not produce is not offered.
- `{MIN_READS}`: minimum read support of an rMATS event: for each group, the mean over its samples of `IJC` plus `SJC`, required in both groups (the smaller of the two group means must reach it). Default: 10. Events below it are dropped from the tests, and their number is reported. Asked only when `rmats` is in `{MODULES}`.

SUPPA2 reports a p-value but no adjusted p-value; the report applies the Benjamini-Hochberg correction to it and says so, and `{FDR}` is applied to the corrected value. Print the thresholds in one line and ask (default: yes): "Use these thresholds?"

---

## Step 8 — Analyses menu

Offer only the sections the run supports; a section that cannot be produced is not offered, and the report says why. Ask one question per section (numbered; default as stated) and store the sections switched on as `{SECTIONS}`, a space-separated subset of `events cross topgene enrich`:
- `events`: event and result tables with a volcano plot for each module in `{MODULES}`. Default: yes. This is the content of the rMATS, SUPPA2 and exon and transcript usage reports.
- `cross`: cross-tool comparison (UpSet plot of the significant genes of each tool and the direction concordance). Offered only when two or more of `rmats`, `suppa`, `dexseq_exon`, `edger_exon` and `dtu` are in `{MODULES}` (each counts as a separate tool) and UpSetR or ComplexUpset is in the image (Step 2). Default: yes when offered.
- `topgene`: top-gene plots (per-sample PSI of the top events, exon usage and isoform proportions), drawn from the pipeline's own tables only, never from BAM files. Default: yes.
- `enrich`: GO enrichment, plus KEGG and Reactome as far as `{ENRICH_DBS}` (Step 2) allows, of the significant genes, with the genes that could be tested as the background. Offered only when `{ORGANISM}` is `human` or `mouse`. Default: yes when offered.

Which Rmd and job each switch controls is fixed in Step 10. `events` and `topgene` are both passed to the module reports as parameters. Ask (default: yes): "Generate these sections: `{SECTIONS}`?"

---

## Step 9 — SLURM resources and project name

Ask one question per value and state the default:
- Memory `{MEM}`: default 32G (at most 64G).
- Time `{TIME}`: default 4:00:00 (at most 4:00:00).
- `{PROJECT}`: the name used in the job names. Default: `{WD_NAME}`.

Every generated script uses partition `bcc`, one node, one task and four CPUs (`-N 1 -n 1 -c 4`). A value above 64G or 4:00:00 is written only when the user asks for it by name; the generator then prints a note that the value is above the default limits. Tell the user the values in one line.

---

## Step 10 — Generate the reports and the SLURM scripts

All files go to `{REPORT_DIR}` (Step 0), which is under `results/`. Nothing is written to `{OUTDIR_RNASPLICE}`: the output of the run is read-only for the skill and is never written to.

The Rmd files have fixed names, and there is one job per Rmd that applies:

| Job script | Rmd file | When the job exists |
|---|---|---|
| `run_01.sbatch` | `01_import.Rmd` | always |
| `run_02.sbatch` | `02_rmats.Rmd` | `rmats` in `{MODULES}` |
| `run_03.sbatch` | `03_suppa.Rmd` | `suppa` in `{MODULES}` |
| `run_04.sbatch` | `04_usage.Rmd` (DEXSeq and edgeR exon usage, and DTU) | any of `dexseq_exon`, `edger_exon`, `dtu` in `{MODULES}` |
| `run_05.sbatch` | `05_cross_tool.Rmd` | `cross` in `{SECTIONS}` and two or more tools in `{MODULES}` |
| `run_06.sbatch` | `06_enrichment.Rmd` | `enrich` in `{SECTIONS}` |
| `run_07.sbatch` | `07_summary.Rmd` (always last) | always |

Default: every job that applies. Write the Rmd files that apply to `{REPORT_DIR}`, then generate the job scripts and `run_all.sh`. Jobs 02 to 06 start when job 01 has finished successfully (`dependency=afterok`, with the option that cancels them when job 01 fails) and run side by side; job 07 starts when jobs 01 to 06 (those that exist) have ended, whatever their outcome (`dependency=afterany`), so the summary also lists the reports that failed.

**Writing the Rmd files.** The content of each Rmd file is the block under its heading **File `<name>.Rmd`.** in the section "Report files" below (a fence of four backticks, language `rmd`). Write that block with the Write tool to `{REPORT_DIR}/<name>.Rmd` byte for byte, without the two fence lines. The Rmd files hold no placeholders: every value of the earlier steps reaches them as a render parameter of the job script (the `params = list(...)` of the block below).

**SLURM scripts.** Paste this block whole into one Bash call, with the placeholders replaced by the literal values of the earlier steps (`{EMAIL}` is `none` when no mail is wanted). It checks the values, refuses a report directory outside `results/` or inside the output of the run, writes one `run_NN.sbatch` per job that applies, and writes `run_all.sh`, which only calls sbatch:
```bash
R="{REPORT_DIR}"
OUT="{OUTDIR_RNASPLICE}"
WORK="{CWD}"
SIF="{SIF}"
EMAIL="{EMAIL}"
JOBNAME="{PROJECT}_report"
MEM="{MEM}"
TIME="{TIME}"
MODULES="{MODULES}"
SECTIONS="{SECTIONS}"
FDR="{FDR}"
DPSI="{DPSI}"
RMATS_TYPE="{RMATS_TYPE}"
MIN_READS="{MIN_READS}"
ORGANISM="{ORGANISM}"
ENRICH_DBS="{ENRICH_DBS}"
GTF="{GTF}"
CONTRASTS_CSV="{CONTRASTS_CSV}"
SAMPLESHEET="{SAMPLESHEET}"
PAIRED_RMATS="{PAIRED_RMATS}"
PAIRED_SUPPA="{PAIRED_SUPPA}"
fail() { echo "ERROR: $1"; exit 1; }
case "$R" in */results/*) ;; *) fail "report directory $R is not under a results/ directory" ;; esac
case "$R/" in "$OUT"/*) fail "report directory $R is inside the pipeline output $OUT" ;; esac
[[ "$MEM" =~ ^[0-9]+[GM]$ ]] || fail "memory $MEM is not of the form 32G"
[[ "$TIME" =~ ^([0-9]+-)?[0-9]+:[0-9][0-9]:[0-9][0-9]$ ]] || fail "time $TIME is not of the form 4:00:00"
[[ "$FDR" =~ ^[0-9]*\.?[0-9]+$ && "$DPSI" =~ ^[0-9]*\.?[0-9]+$ && "$MIN_READS" =~ ^[0-9]+$ ]] || fail "FDR, DPSI or MIN_READS is not a number"
case "$RMATS_TYPE" in JC|JCEC) ;; *) fail "RMATS_TYPE $RMATS_TYPE is not JC or JCEC" ;; esac
# memory in MB and time in seconds, for the note about the default limits
case "$MEM" in *G) mb=$(( ${MEM%G} * 1024 )) ;; *) mb=${MEM%M} ;; esac
tt=$TIME; d=0; case "$tt" in *-*) d=${tt%%-*}; tt=${tt#*-} ;; esac
hh=${tt%%:*}; mm=${tt#*:}; mm=${mm%%:*}; ss=${tt##*:}
secs=$(( ((10#$d * 24 + 10#$hh) * 60 + 10#$mm) * 60 + 10#$ss ))
if [ "$mb" -gt 65536 ] || [ "$secs" -gt 14400 ]; then echo "NOTE: $MEM and $TIME are above the default limits (64G, 4:00:00), as requested"; fi
has() { case " $MODULES " in *" $1 "*) return 0 ;; esac; return 1; }
sec() { case " $SECTIONS " in *" $1 "*) return 0 ;; esac; return 1; }
tf() { if "$@"; then echo TRUE; else echo FALSE; fi; }
EV=$(tf sec events); CR=$(tf sec cross); TG=$(tf sec topgene); EN=$(tf sec enrich)
PR=$(tf test "$PAIRED_RMATS" = true); PS=$(tf test "$PAIRED_SUPPA" = true)
ntools=0
for m in rmats suppa dexseq_exon edger_exon dtu; do if has $m; then ntools=$((ntools + 1)); fi; done
# the jobs that apply, as number:Rmd
JOBS="01:01_import.Rmd"
if has rmats; then JOBS="$JOBS 02:02_rmats.Rmd"; fi
if has suppa; then JOBS="$JOBS 03:03_suppa.Rmd"; fi
if has dexseq_exon || has edger_exon || has dtu; then JOBS="$JOBS 04:04_usage.Rmd"; fi
if sec cross && [ $ntools -ge 2 ]; then JOBS="$JOBS 05:05_cross_tool.Rmd"; fi
if sec enrich; then JOBS="$JOBS 06:06_enrichment.Rmd"; fi
JOBS="$JOBS 07:07_summary.Rmd"
for j in $JOBS; do [ -f "$R/${j#*:}" ] || fail "$R/${j#*:} does not exist: write the Rmd files first"; done
mkdir -p "$R/logs" || fail "cannot create $R/logs"
# bind the cluster share, and the directory of any path that lies outside it
BIND=/net/bmc-lab3
for p in "$R" "$OUT" "$WORK" "$SIF" "$GTF" "$CONTRASTS_CSV" "$SAMPLESHEET"; do
  case "$p" in /net/bmc-lab3/*) ;; /*) q=$(dirname "$p"); case ",$BIND," in *",$q,"*) ;; *) BIND="$BIND,$q" ;; esac ;; esac
done
write_job() {
  f="$R/run_$1.sbatch"
  cat > "$f" <<EOF
#!/bin/bash
#SBATCH -J ${JOBNAME}_$1
#SBATCH -p bcc
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=$MEM
#SBATCH -t $TIME
#SBATCH -o $R/logs/run_$1.%j.log
EOF
  if [ -n "$EMAIL" ] && [ "$EMAIL" != none ]; then
    cat >> "$f" <<EOF
#SBATCH --mail-user=$EMAIL
#SBATCH --mail-type=END,FAIL
EOF
  fi
  cat >> "$f" <<EOF
module add singularity/3.10.4 || exit 1
cd "$R" || exit 1
singularity exec --bind $BIND "$SIF" Rscript -e "rmarkdown::render('$R/$2',
  output_dir = '$R', knit_root_dir = '$R',
  params = list(
    outdir = '$OUT', report_dir = '$R',
    contrasts = '$CONTRASTS_CSV', samplesheet = '$SAMPLESHEET',
    modules = '$MODULES', rmats_type = '$RMATS_TYPE',
    fdr = $FDR, dpsi = $DPSI, min_reads = $MIN_READS,
    events = $EV, cross = $CR, topgene = $TG, enrich = $EN,
    organism = '$ORGANISM', gtf = '$GTF',
    enrich_dbs = '$ENRICH_DBS',
    paired_rmats = $PR, paired_suppa = $PS))"
EOF
}
cat > "$R/run_all.sh" <<EOF
#!/bin/bash
# Submits the report jobs of $R in order; light work only, every job runs R on a compute node.
cd "$R" || exit 1
EOF
ALLDEP=""; ALLJ=""; NUMS=""
for j in $JOBS; do
  n=${j%%:*}; rmd=${j#*:}
  write_job "$n" "$rmd"
  case "$n" in
    01)
      cat >> "$R/run_all.sh" <<EOF
JOB$n=\$(sbatch --parsable run_$n.sbatch) || exit 1
EOF
      ;;
    07)
      cat >> "$R/run_all.sh" <<EOF
JOB$n=\$(sbatch --parsable --dependency=afterany:$ALLDEP run_$n.sbatch) || exit 1
EOF
      ;;
    *)
      cat >> "$R/run_all.sh" <<EOF
JOB$n=\$(sbatch --parsable --dependency=afterok:\$JOB01 --kill-on-invalid-dep=yes run_$n.sbatch) || exit 1
EOF
      ;;
  esac
  case "$n" in 07) ;; *) ALLDEP="${ALLDEP:+$ALLDEP:}\$JOB$n" ;; esac
  ALLJ="$ALLJ \$JOB$n"; NUMS="$NUMS $n"
done
cat >> "$R/run_all.sh" <<EOF
echo "submitted:$ALLJ"
EOF
ls -la "$R"
echo "JOBS=${NUMS# }"
echo "SUBMIT=bash $R/run_all.sh"
```

The block prints the files written, `JOBS=` with the job numbers, and `SUBMIT=` with the one command that submits everything.

**Final instructions.** Print this to the user, with `{REPORT_DIR}` replaced, and nothing else to run:
```
Reports are ready to submit. Run this one command:
  bash {REPORT_DIR}/run_all.sh
It submits the jobs in order (each waits for the one it needs). The reports and the summary page appear in {REPORT_DIR}; the logs are in {REPORT_DIR}/logs.
Check progress with: squeue -u $USER
```
`run_all.sh` calls only sbatch; R runs inside the jobs, on compute nodes. When a job fails, read its log, fix the cause, and submit the jobs that did not run with sbatch from `{REPORT_DIR}`.

---

## Report files

Each block below is the whole content of one Rmd file of `{REPORT_DIR}` (see **Writing the Rmd files.** in Step 10). The values of the wizard reach the Rmd as render parameters; the parameter names are those of the job scripts: `outdir`, `report_dir`, `contrasts`, `samplesheet`, `modules`, `fdr`, `dpsi`, `rmats_type`, `min_reads`, `events`, `cross`, `topgene`, `enrich`, `organism`, `enrich_dbs`, `gtf`, `paired_rmats`, `paired_suppa`.

The import report `01_import.Rmd` reads the table of every module in `{MODULES}` for every contrast, converts each effect to treatment minus control, and writes one RDS per contrast and tool (`rds/{CONTRAST}.{tool}.rds`, tool one of `rmats`, `suppa`, `exon_dexseq`, `exon_edger`, `dtu`) with the columns `contrast, tool, gene_id, gene_symbol, feature_type, feature_id, effect, effect_scale, padj, direction, n_support`, plus `rds/{CONTRAST}.qc.rds` with the counts the report prints. A missing or renamed column stops the import with a message that names the file and the column.

**File `01_import.Rmd`.**
````rmd
---
title: "Import and harmonise the nf-core/rnasplice tables"
output:
  html_document:
    toc: true
params:
  outdir: ""
  report_dir: "."
  contrasts: ""
  samplesheet: ""
  modules: ""
  fdr: 0.05
  dpsi: 0.1
  rmats_type: "JC"
  min_reads: 10
  events: true
  cross: true
  topgene: true
  enrich: true
  organism: "human"
  enrich_dbs: "GO"
  gtf: ""
  paired_rmats: false
  paired_suppa: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)
```

This report reads the result tables of the nf-core/rnasplice run in `r params$outdir` (read only) and writes one table per contrast and tool to `r file.path(params$report_dir, "rds")`. Every effect is stated as **treatment minus control**: a positive effect is higher in the treatment group (more inclusion for rMATS and SUPPA2, higher usage for the exon and transcript usage tools).

## Sign rules

| Tool | Column read | Convention of the tool | Conversion to treatment minus control |
|---|---|---|---|
| rMATS | `IncLevelDifference` | IncLevel1 minus IncLevel2, group 1 (b1) is the treatment: treatment minus control | kept |
| SUPPA2 | `{kind}_{T}-{kind}_{C}_dPSI` | control minus treatment | sign flipped |
| DEXSeq exon usage | `log2fold_{C}_{T}` | control versus treatment | sign flipped |
| DEXSeq transcript usage (DTU) | `log2fold_{C}_{T}` | control versus treatment | sign flipped |
| edgeR exon usage | `logFC` | treatment minus control | kept |

`direction` is `treatment_up` for a positive effect, `control_up` for a negative one and `none` for zero or a missing effect.

SUPPA2 reports a p-value but no adjusted p-value. This report applies the Benjamini-Hochberg correction itself, separately for each contrast and for each file (local events, isoforms). An event whose dPSI is `nan` was not tested: it is kept with a missing effect, counted as untested below, and left out of the correction (SUPPA2 writes a p-value of 1 for it, which is not a test). rMATS events whose read support is below `r params$min_reads` (for each group the mean over its samples of IJC plus SJC; the smaller of the two group means) are removed from the rMATS table, and their number is reported.

```{r functions}
TOOL_OF_MODULE <- c(rmats = "rmats", suppa = "suppa", dexseq_exon = "exon_dexseq", edger_exon = "exon_edger", dtu = "dtu")
CONTRACT_COLUMNS <- c("contrast", "tool", "gene_id", "gene_symbol", "feature_type", "feature_id", "effect", "effect_scale",
                      "padj", "direction", "n_support")
RMATS_TYPES <- c("SE", "A5SS", "A3SS", "MXE", "RI")

harmonise_effect <- function(tool, value) {
  # returns the effect as treatment minus control
  switch(tool,
    rmats = value,          # IncLevelDifference: b1 = treatment, already treatment - control
    suppa = -value,         # SUPPA2 dPSI is control - treatment
    exon_dexseq = -value,   # log2fold_{CONTROL}_{TREATMENT}
    dtu = -value,           # log2fold_{CONTROL}_{TREATMENT}
    exon_edger = value,     # logFC is treatment - control
    stop("unknown tool: ", tool))
}
direction_of <- function(effect) {
  ifelse(is.na(effect) | effect == 0, "none", ifelse(effect > 0, "treatment_up", "control_up"))
}

# rMATS writes GeneID and geneSymbol in double quotes; the quotes are removed from every text field read here
strip_quotes <- function(x) gsub("^\"|\"$", "", x)

# a repeated column name gets the suffix _2, _3, ...: rMATS writes the column ID twice (first, and after the coordinates)
dedupe_names <- function(h) {
  k <- stats::ave(seq_along(h), h, FUN = seq_along)
  ifelse(k > 1, paste0(h, "_", k), h)
}

# reads a delimited file as text columns; a file with a header only gives zero rows with the header's columns.
# row_label = TRUE: the rows have one field more than the header (SUPPA2 writes the row label without a header cell)
read_delim_chr <- function(path, sep, row_label = FALSE) {
  if (!file.exists(path)) stop("file ", path, " not found: none of its columns can be read")
  l <- sub("\r$", "", readLines(path, warn = FALSE))
  l <- l[nzchar(l)]
  if (length(l) == 0) stop("file ", path, " is empty: no header line, so no column can be found")
  nm <- strsplit(l[1], sep, fixed = TRUE)[[1]]
  if (row_label) nm <- c("row_label", nm)
  if (length(l) == 1) {
    tab <- as.data.frame(stats::setNames(replicate(length(nm), character(0), simplify = FALSE), nm),
                         check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    tab <- tryCatch(utils::read.table(text = l[-1], sep = sep, quote = "", header = FALSE, colClasses = "character",
                                      na.strings = character(0), comment.char = "", col.names = paste0("V", seq_along(nm)),
                                      fill = FALSE, stringsAsFactors = FALSE),
                    error = function(e) stop("file ", path, ": the rows do not match the ", length(nm),
                                             " columns of the header (", paste(nm, collapse = ", "), "): ", conditionMessage(e)))
    if (ncol(tab) != length(nm)) stop("file ", path, ": ", ncol(tab), " fields per row, but the header has the ", length(nm),
                                      " columns ", paste(nm, collapse = ", "))
    names(tab) <- nm
  }
  tab
}

# stops naming the file and every missing column
need_cols <- function(tab, cols, path) {
  miss <- setdiff(cols, names(tab))
  if (length(miss)) stop("file ", path, ": column not found: ", paste(miss, collapse = ", "),
                         " (columns present: ", paste(names(tab), collapse = ", "), ")")
}

# text to number; NA, nan and empty fields are missing values, any other non-number stops naming the file and the column
to_num <- function(x, path, col) {
  x <- trimws(strip_quotes(x))
  miss <- x %in% c("", "NA", "nan", "NaN", "-nan")
  v <- suppressWarnings(as.numeric(ifelse(miss, NA, x)))
  bad <- is.na(v) & !miss
  if (any(bad)) stop("file ", path, ": column ", col, " has a value that is not a number: '", x[bad][1], "'")
  v
}

# Benjamini-Hochberg over the non-missing p-values; a missing p-value stays missing
bh_adjust <- function(p) {
  out <- rep(NA_real_, length(p))
  ok <- !is.na(p)
  out[ok] <- stats::p.adjust(p[ok], method = "BH")
  out
}

harmonised_table <- function(contrast, tool, gene_id, gene_symbol, feature_type, feature_id, value, padj, n_support) {
  n <- length(gene_id)
  if (n == 0) feature_id <- character(0)   # paste0 over an empty column returns one string, not none
  effect <- harmonise_effect(tool, value)
  data.frame(contrast = rep(contrast, n), tool = rep(tool, n), gene_id = as.character(gene_id),
             gene_symbol = rep_len(as.character(gene_symbol), n), feature_type = rep_len(as.character(feature_type), n),
             feature_id = as.character(feature_id), effect = as.numeric(effect),
             effect_scale = rep(if (tool %in% c("rmats", "suppa")) "dPSI" else "log2FC", n), padj = as.numeric(padj),
             direction = as.character(direction_of(effect)), n_support = rep_len(as.numeric(n_support), n),
             stringsAsFactors = FALSE)
}

read_contrasts <- function(path) {
  tab <- read_delim_chr(path, ",")
  need_cols(tab, c("contrast", "treatment", "control"), path)
  tab <- tab[nzchar(tab$contrast), c("contrast", "treatment", "control")]
  if (nrow(tab) == 0) stop("file ", path, ": no row with a value in the column contrast")
  rownames(tab) <- NULL
  tab
}

# rMATS: per event one count per sample, as a comma-separated list: IJC_SAMPLE_1 and SJC_SAMPLE_1 for group 1 (b1, the
# treatment), IJC_SAMPLE_2 and SJC_SAMPLE_2 for group 2 (the control).
# support of a group = mean over its samples of (IJC + SJC); n_support = min(support of group 1, support of group 2).
# An event is kept for the tests when n_support >= min_reads.
rmats_support <- function(tab, path) {
  grp <- function(g) {
    ci <- paste0("IJC_SAMPLE_", g); cs <- paste0("SJC_SAMPLE_", g)
    vapply(seq_len(nrow(tab)), function(i) {
      a <- to_num(strsplit(tab[[ci]][i], ",", fixed = TRUE)[[1]], path, ci)
      b <- to_num(strsplit(tab[[cs]][i], ",", fixed = TRUE)[[1]], path, cs)
      if (length(a) != length(b)) stop("file ", path, ": columns ", ci, " and ", cs, " list different numbers of samples in row ", i)
      mean(a + b)
    }, numeric(1))
  }
  pmin(grp(1), grp(2))
}

read_rmats_raw <- function(path) {
  tab <- read_delim_chr(path, "\t")
  names(tab) <- dedupe_names(names(tab))
  tab
}

read_rmats <- function(path, type, contrast) {
  tab <- read_rmats_raw(path)
  need_cols(tab, c("ID", "ID_2", "GeneID", "geneSymbol", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2",
                   "FDR", "IncLevelDifference"), path)
  if (any(tab$ID != tab$ID_2)) stop("file ", path, ": the two ID columns (ID and ID_2, the second ID) differ")
  harmonised_table(contrast, "rmats", strip_quotes(tab$GeneID), strip_quotes(tab$geneSymbol), type, paste0(type, ":", tab$ID),
                   to_num(tab$IncLevelDifference, path, "IncLevelDifference"), to_num(tab$FDR, path, "FDR"),
                   rmats_support(tab, path))
}

# SUPPA2 column names, built from the contrasts file (kind is local or transcript), never split from file names
suppa_columns <- function(kind, treatment, control) paste0(kind, "_", treatment, "-", kind, "_", control, c("_dPSI", "_p-val"))

read_suppa <- function(path, kind, treatment, control, contrast) {
  tab <- read_delim_chr(path, "\t", row_label = TRUE)
  cols <- suppa_columns(kind, treatment, control)
  need_cols(tab, cols, path)
  lab <- strip_quotes(tab$row_label)
  if (!all(grepl(";", lab, fixed = TRUE)))
    stop("file ", path, ": the row label column (GENE;EVENT) has a label without ';': ", lab[!grepl(";", lab, fixed = TRUE)][1])
  gene <- sub(";.*$", "", lab)      # the text before the first ; only
  event <- sub("^[^;]*;", "", lab)  # the rest, which may hold further ; characters
  dpsi <- to_num(tab[[cols[1]]], path, cols[1])
  p <- to_num(tab[[cols[2]]], path, cols[2])
  p[is.na(dpsi)] <- NA              # dPSI nan: not tested, left out of the correction
  ftype <- if (kind == "local") sub(":.*$", "", event) else rep("isoform", length(event))
  harmonised_table(contrast, "suppa", gene, NA_character_, ftype, event, dpsi, bh_adjust(p), NA_real_)
}

# DEXSeq column of the fold change, built from the contrasts file: log2fold_{CONTROL}_{TREATMENT}
dexseq_fold_column <- function(treatment, control) paste0("log2fold_", control, "_", treatment)

read_dexseq <- function(path, tool, treatment, control, contrast) {
  fc <- dexseq_fold_column(treatment, control)
  if (!tool %in% c("exon_dexseq", "dtu")) stop("file ", path, ": read_dexseq reads exon_dexseq or dtu tables (column ", fc, "), not ", tool)
  tab <- read_delim_chr(path, if (grepl("\\.csv$", path)) "," else "\t")
  need_cols(tab, c("groupID", "featureID", "padj", fc), path)
  gid <- strip_quotes(tab$groupID)
  fid <- strip_quotes(tab$featureID)
  feat <- if (tool == "exon_dexseq") paste0(gid, ":", fid) else fid
  harmonised_table(contrast, tool, gid, NA_character_, if (tool == "exon_dexseq") "exon_bin" else "transcript", feat,
                   to_num(tab[[fc]], path, fc), to_num(tab$padj, path, "padj"), NA_real_)
}

read_edger_exon <- function(path, contrast) {
  tab <- read_delim_chr(path, ",")
  need_cols(tab, c("Geneid", "Chr", "Start", "End", "logFC", "FDR"), path)
  gid <- strip_quotes(tab$Geneid)
  feat <- paste0(gid, ":", strip_quotes(tab$Chr), ":", tab$Start, "-", tab$End)
  harmonised_table(contrast, "exon_edger", gid, NA_character_, "exon", feat, to_num(tab$logFC, path, "logFC"),
                   to_num(tab$FDR, path, "FDR"), NA_real_)
}

# counts per tool: rows, tested (padj present), significant (padj < fdr, and abs(effect) >= dpsi on the dPSI scale)
count_table <- function(tab, fdr, dpsi) {
  sig <- !is.na(tab$padj) & tab$padj < fdr & !is.na(tab$effect) & (tab$effect_scale != "dPSI" | abs(tab$effect) >= dpsi)
  data.frame(rows = nrow(tab), tested = sum(!is.na(tab$padj)), effect_missing = sum(is.na(tab$effect)),
             treatment_up = sum(tab$direction == "treatment_up"), control_up = sum(tab$direction == "control_up"),
             none = sum(tab$direction == "none"), significant = sum(sig),
             significant_treatment_up = sum(sig & tab$direction == "treatment_up"),
             significant_control_up = sum(sig & tab$direction == "control_up"))
}
```

```{r import}
modules <- strsplit(trimws(params$modules), "[[:space:]]+")[[1]]
modules <- modules[nzchar(modules)]
unknown <- setdiff(modules, names(TOOL_OF_MODULE))
if (length(unknown)) stop("parameter modules: unknown module ", unknown[1], " (known: ", paste(names(TOOL_OF_MODULE), collapse = ", "), ")")
if (!params$rmats_type %in% c("JC", "JCEC")) stop("parameter rmats_type: ", params$rmats_type, " is not JC or JCEC")
rds_dir <- file.path(params$report_dir, "rds")
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
contrasts <- read_contrasts(params$contrasts)
path_of <- function(...) file.path(params$outdir, ...)
results <- list()
for (i in seq_len(nrow(contrasts))) {
  C <- contrasts$contrast[i]; Tr <- contrasts$treatment[i]; Co <- contrasts$control[i]
  qc <- list(contrast = C, treatment = Tr, control = Co,
             thresholds = c(fdr = params$fdr, dpsi = params$dpsi, min_reads = params$min_reads),
             rmats_type = params$rmats_type, paired_rmats = params$paired_rmats, paired_suppa = params$paired_suppa,
             files_read = character(0), files_missing = character(0), rmats_dir = NA_character_,
             rmats_low_support = NA_integer_, suppa_untested = c(local = NA_integer_, isoform = NA_integer_))
  out <- list()
  if ("rmats" %in% modules) {
    # the directory of the contrast or, for a run with rmats_paired_stats true, the one with the suffix _paired (as in the
    # module detection); when both hold tables, the one that matches the run's pairing switch is read
    dirs <- file.path("star", "rmats", c(C, paste0(C, "_paired")))
    if (isTRUE(params$paired_rmats)) dirs <- rev(dirs)
    files_in <- function(d) path_of(d, sprintf("%s.MATS.%s.txt", RMATS_TYPES, params$rmats_type))
    hit <- dirs[vapply(dirs, function(d) any(file.exists(files_in(d))), logical(1))]
    if (length(hit) == 0) {
      qc$files_missing <- c(qc$files_missing, files_in(dirs[1]))
    } else {
      qc$rmats_dir <- path_of(hit[1])
      parts <- list()
      for (ty in RMATS_TYPES) {
        f <- path_of(hit[1], sprintf("%s.MATS.%s.txt", ty, params$rmats_type))
        if (file.exists(f)) { parts[[ty]] <- read_rmats(f, ty, C); qc$files_read <- c(qc$files_read, f) }
        else qc$files_missing <- c(qc$files_missing, f)
      }
      tab <- do.call(rbind, unname(parts))
      keep <- !is.na(tab$n_support) & tab$n_support >= params$min_reads
      qc$rmats_low_support <- sum(!keep)
      out$rmats <- tab[keep, , drop = FALSE]
    }
  }
  if ("suppa" %in% modules) {
    lf <- path_of("salmon", "suppa", "diffsplice", "per_local_event", sprintf("local_%s-%s.dpsi", Tr, Co))
    tf <- path_of("salmon", "suppa", "diffsplice", "per_isoform", sprintf("transcript_%s-%s.dpsi", Tr, Co))
    parts <- list()
    if (file.exists(lf)) {
      parts$local <- read_suppa(lf, "local", Tr, Co, C); qc$files_read <- c(qc$files_read, lf)
      qc$suppa_untested[["local"]] <- sum(is.na(parts$local$effect))
    } else qc$files_missing <- c(qc$files_missing, lf)
    if (file.exists(tf)) {
      parts$isoform <- read_suppa(tf, "transcript", Tr, Co, C); qc$files_read <- c(qc$files_read, tf)
      qc$suppa_untested[["isoform"]] <- sum(is.na(parts$isoform$effect))
    } else qc$files_missing <- c(qc$files_missing, tf)
    if (length(parts)) out$suppa <- do.call(rbind, unname(parts))
  }
  if ("dexseq_exon" %in% modules) {
    f <- path_of("star", "dexseq_exon", "results", sprintf("DEXSeqResults.%s.csv", C))
    if (file.exists(f)) { out$exon_dexseq <- read_dexseq(f, "exon_dexseq", Tr, Co, C); qc$files_read <- c(qc$files_read, f) }
    else qc$files_missing <- c(qc$files_missing, f)
  }
  if ("edger_exon" %in% modules) {
    f <- path_of("star", "edger", sprintf("contrast_%s.usage.exon.csv", C))
    if (file.exists(f)) { out$exon_edger <- read_edger_exon(f, C); qc$files_read <- c(qc$files_read, f) }
    else qc$files_missing <- c(qc$files_missing, f)
  }
  if ("dtu" %in% modules) {
    f <- path_of("salmon", "dexseq_dtu", "results", "dexseq", sprintf("DEXSeqResults.%s-%s.tsv", Tr, Co))
    if (file.exists(f)) { out$dtu <- read_dexseq(f, "dtu", Tr, Co, C); qc$files_read <- c(qc$files_read, f) }
    else qc$files_missing <- c(qc$files_missing, f)
  }
  for (tool in names(out)) {
    rownames(out[[tool]]) <- NULL
    stopifnot(identical(names(out[[tool]]), CONTRACT_COLUMNS))
    saveRDS(out[[tool]], file.path(rds_dir, sprintf("%s.%s.rds", C, tool)))
  }
  qc$counts <- if (length(out)) cbind(tool = names(out), do.call(rbind, lapply(out, count_table, fdr = params$fdr, dpsi = params$dpsi)))
               else NULL
  if (!is.null(qc$counts)) rownames(qc$counts) <- NULL
  saveRDS(qc, file.path(rds_dir, sprintf("%s.qc.rds", C)))
  results[[C]] <- qc
}
```

## Samples and pairing

```{r samples, results="asis"}
ss <- read_delim_chr(params$samplesheet, ",")
need_cols(ss, c("sample", "condition"), params$samplesheet)
cat("Samples per condition in", params$samplesheet, "\n\n")
cat(knitr::kable(as.data.frame(table(condition = ss$condition), responseName = "samples"), format = "markdown"), sep = "\n")
cat("\n\nrMATS:", if (isTRUE(params$paired_rmats)) "paired statistics" else "unpaired statistics",
    "; SUPPA2:", if (isTRUE(params$paired_suppa)) "paired statistics" else "unpaired statistics",
    "(read from the run's parameters). The exon and transcript usage tools use the design of the sample sheet as it is.\n")
```

## Per contrast: files, sign and counts

```{r qc, results="asis"}
for (qc in results) {
  cat("\n\n### Contrast ", qc$contrast, "\n\n", sep = "")
  cat("Treatment `", qc$treatment, "`, control `", qc$control, "`: every effect below is `", qc$treatment, "` minus `",
      qc$control, "`. SUPPA2 dPSI is control minus treatment and DEXSeq `", dexseq_fold_column(qc$treatment, qc$control),
      "` is control versus treatment: both were sign flipped; rMATS and edgeR were kept.\n\n", sep = "")
  cat("Thresholds: FDR ", params$fdr, ", abs(dPSI) ", params$dpsi, " (rMATS and SUPPA2 only), rMATS minimum read support ",
      params$min_reads, ", rMATS counts ", params$rmats_type, ".\n\n", sep = "")
  if (length(qc$files_read)) cat("Files read:\n\n", paste0("- `", qc$files_read, "`\n"), "\n", sep = "")
  if (length(qc$files_missing)) cat("Files not found (their section is left out):\n\n", paste0("- `", qc$files_missing, "`\n"), "\n", sep = "")
  if (!is.na(qc$rmats_dir)) cat("rMATS tables read from `", qc$rmats_dir, "`; events removed for read support below ",
                                params$min_reads, ": ", qc$rmats_low_support, ".\n\n", sep = "")
  if (!all(is.na(qc$suppa_untested))) cat("SUPPA2 events not tested (dPSI nan; kept with a missing effect, left out of the Benjamini-Hochberg correction): local events ",
                                          qc$suppa_untested[["local"]], ", isoforms ", qc$suppa_untested[["isoform"]], ".\n\n", sep = "")
  if (!is.null(qc$counts)) cat(knitr::kable(qc$counts, format = "markdown"), sep = "\n")
  else cat("No table was read for this contrast.\n")
}
```
````
