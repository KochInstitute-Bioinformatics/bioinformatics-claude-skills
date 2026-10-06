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
- All three are TRUE: store `{ENRICH_DBS}` = `GO KEGG Reactome` (the enrichment report runs GO only; Step 8 says so).
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
- `{MIN_READS}`: minimum read support of an rMATS event: for each group, the mean over its samples of `IJC` plus `SJC`, required in both groups (the smaller of the two group means must reach it). Default: 10. Events below it are flagged (`low_support`), kept in the tables and excluded from the significance counts of the reports, never deleted; their number is reported, and the rMATS FDR is not recomputed (it stays as rMATS wrote it). On very shallow data (like the nf-core test data) most events fall below 10; the user may lower it. Asked only when `rmats` is in `{MODULES}`.

SUPPA2 reports a p-value but no adjusted p-value; the report applies the Benjamini-Hochberg correction to it and says so, and `{FDR}` is applied to the corrected value. Print the thresholds in one line and ask (default: yes): "Use these thresholds?"

---

## Step 8 — Analyses menu

Offer only the sections the run supports; a section that cannot be produced is not offered, and the report says why. Ask one question per section (numbered; default as stated) and store the sections switched on as `{SECTIONS}`, a space-separated subset of `events cross topgene enrich`:
- `events`: event and result tables with a volcano plot for each module in `{MODULES}`. Default: yes. This is the content of the rMATS, SUPPA2 and exon and transcript usage reports.
- `cross`: cross-tool comparison (UpSet plot of the significant genes of each tool and the direction concordance). Offered only when two or more of `rmats`, `suppa`, `dexseq_exon`, `edger_exon` and `dtu` are in `{MODULES}` (each counts as a separate tool) and UpSetR or ComplexUpset is in the image (Step 2). Default: yes when offered.
- `topgene`: top-gene plots (per-sample PSI of the top events, exon usage and isoform proportions), drawn from the pipeline's own tables only, never from BAM files. Default: yes.
- `enrich`: GO enrichment (biological process, molecular function, cellular component) of the significant genes of each tool, in the groups treatment_up, control_up and all significant genes (a group needs at least 10 genes), with the genes that were testable in that tool as the background. Only GO is run: KEGG and Reactome are never run by this report, and when they are asked for in `{ENRICH_DBS}` the report prints "KEGG/Reactome require an image rebuild". Offered only when `{ORGANISM}` is `human` or `mouse` (any other organism is answered in the report with "organism not supported"). Default: yes when offered.

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
  OUTF=""; if [ "$1" = 07 ]; then OUTF="output_file = 'report_summary.html', "; fi
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
  output_dir = '$R', ${OUTF}knit_root_dir = '$R',
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
  Start with the summary page, the index of every report, figure and table: {REPORT_DIR}/report_summary.html
Check progress with: squeue -u $USER
```
`run_all.sh` calls only sbatch; R runs inside the jobs, on compute nodes. When a job fails, read its log, fix the cause, and submit the jobs that did not run with sbatch from `{REPORT_DIR}`.

---

## Report files

Each block below is the whole content of one Rmd file of `{REPORT_DIR}` (see **Writing the Rmd files.** in Step 10). The values of the wizard reach the Rmd as render parameters; the parameter names are those of the job scripts: `outdir`, `report_dir`, `contrasts`, `samplesheet`, `modules`, `fdr`, `dpsi`, `rmats_type`, `min_reads`, `events`, `cross`, `topgene`, `enrich`, `organism`, `enrich_dbs`, `gtf`, `paired_rmats`, `paired_suppa`.

The import report `01_import.Rmd` reads the table of every module in `{MODULES}` for every contrast, converts each effect to treatment minus control, and writes one RDS per contrast and tool (`rds/{CONTRAST}.{tool}.rds`, tool one of `rmats`, `suppa`, `exon_dexseq`, `exon_edger`, `dtu`) with the columns `contrast, tool, gene_id, gene_symbol, feature_type, feature_id, effect, effect_scale, padj, direction, n_support, low_support`, plus `rds/{CONTRAST}.qc.rds` with the counts the report prints. A missing or renamed column stops the import with a message that names the file and the column. When the file of a module is missing for a contrast, no RDS is written for that tool and contrast and the file is listed in `qc$files_missing`; Rmds 02 to 06 check that each RDS exists before they read it. `low_support` is TRUE for an rMATS event whose read support is below `{MIN_READS}` and FALSE for every other row and tool; flagged rows stay in the RDS, and the reports exclude them from the significant calls but count and report them.

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

SUPPA2 reports a p-value but no adjusted p-value. This report applies the Benjamini-Hochberg correction itself, separately for each contrast and for each file (local events, isoforms). An event whose dPSI is `nan` was not tested: it is kept with a missing effect, counted as untested below, and left out of the correction (SUPPA2 writes a p-value of 1 for it, which is not a test). rMATS events whose read support is below `r params$min_reads` (for each group the mean over its samples of IJC plus SJC; the smaller of the two group means) are flagged in the column `low_support`: they stay in the table, are excluded from the significant counts and are counted below; nothing is deleted, and the rMATS FDR is not recomputed.

```{r functions}
TOOL_OF_MODULE <- c(rmats = "rmats", suppa = "suppa", dexseq_exon = "exon_dexseq", edger_exon = "exon_edger", dtu = "dtu")
CONTRACT_COLUMNS <- c("contrast", "tool", "gene_id", "gene_symbol", "feature_type", "feature_id", "effect", "effect_scale",
                      "padj", "direction", "n_support", "low_support")
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

harmonised_table <- function(contrast, tool, gene_id, gene_symbol, feature_type, feature_id, value, padj, n_support,
                             low_support = FALSE) {
  n <- length(gene_id)
  if (n == 0) feature_id <- character(0)   # paste0 over an empty column returns one string, not none
  effect <- harmonise_effect(tool, value)
  data.frame(contrast = rep(contrast, n), tool = rep(tool, n), gene_id = as.character(gene_id),
             gene_symbol = rep_len(as.character(gene_symbol), n), feature_type = rep_len(as.character(feature_type), n),
             feature_id = as.character(feature_id), effect = as.numeric(effect),
             effect_scale = rep(if (tool %in% c("rmats", "suppa")) "dPSI" else "log2FC", n), padj = as.numeric(padj),
             direction = as.character(direction_of(effect)), n_support = rep_len(as.numeric(n_support), n),
             low_support = rep_len(as.logical(low_support), n),
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
# An event is flagged low_support when n_support < min_reads (or missing); it stays in the table, never deleted.
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

read_rmats <- function(path, type, contrast, min_reads) {
  tab <- read_rmats_raw(path)
  need_cols(tab, c("ID", "ID_2", "GeneID", "geneSymbol", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2",
                   "FDR", "IncLevelDifference"), path)
  if (any(tab$ID != tab$ID_2)) stop("file ", path, ": the two ID columns (ID and ID_2, the second ID) differ")
  ns <- rmats_support(tab, path)
  low <- is.na(ns) | ns < min_reads
  harmonised_table(contrast, "rmats", strip_quotes(tab$GeneID), strip_quotes(tab$geneSymbol), type, paste0(type, ":", tab$ID),
                   to_num(tab$IncLevelDifference, path, "IncLevelDifference"), to_num(tab$FDR, path, "FDR"), ns, low)
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

# counts per tool: rows, tested (padj present), low_support (flagged, rMATS only), significant (padj < fdr, abs(effect) >= dpsi
# on the dPSI scale, and not low_support)
count_table <- function(tab, fdr, dpsi) {
  sig <- !is.na(tab$padj) & tab$padj < fdr & !is.na(tab$effect) & (tab$effect_scale != "dPSI" | abs(tab$effect) >= dpsi) &
    !tab$low_support
  data.frame(rows = nrow(tab), tested = sum(!is.na(tab$padj)), low_support = sum(tab$low_support), effect_missing = sum(is.na(tab$effect)),
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
        if (file.exists(f)) { parts[[ty]] <- read_rmats(f, ty, C, params$min_reads); qc$files_read <- c(qc$files_read, f) }
        else qc$files_missing <- c(qc$files_missing, f)
      }
      tab <- do.call(rbind, unname(parts))
      qc$rmats_low_support <- sum(tab$low_support)
      out$rmats <- tab
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
  if (!is.na(qc$rmats_dir)) cat("rMATS tables read from `", qc$rmats_dir, "`; events flagged low_support (read support below ",
                                params$min_reads, "; kept in the table, excluded from the significant counts, FDR as written by rMATS): ",
                                qc$rmats_low_support, ".\n\n", sep = "")
  if (!all(is.na(qc$suppa_untested))) cat("SUPPA2 events not tested (dPSI nan; kept with a missing effect, left out of the Benjamini-Hochberg correction): local events ",
                                          qc$suppa_untested[["local"]], ", isoforms ", qc$suppa_untested[["isoform"]], ".\n\n", sep = "")
  if (!is.null(qc$counts)) cat(knitr::kable(qc$counts, format = "markdown"), sep = "\n")
  else cat("No table was read for this contrast.\n")
}
```
````

The rMATS report `02_rmats.Rmd` reads `rds/{CONTRAST}.rmats.rds` (written by `01_import.Rmd`) for every contrast and writes, per contrast, `tables/{CONTRAST}.rmats.counts.csv`, `tables/{CONTRAST}.rmats.{type}.significant.csv`, `figures/{CONTRAST}.rmats.volcano.{type}.png`, `figures/{CONTRAST}.rmats.effect_hist.png` and `figures/{CONTRAST}.rmats.psi_top.png` (`type` is one of `SE`, `A5SS`, `A3SS`, `MXE`, `RI`). An event is significant when its adjusted p-value is below `{FDR}`, its absolute dPSI is at least `{DPSI}`, and it is not flagged `low_support`. The per-sample PSI of the top events is not in the RDS; the report reads `IncLevel1` and `IncLevel2` again from the rMATS table of the same directory the import used (the contrast's own, or the one with the suffix `_paired`). The sample order of the two lists is that of rMATS's own sample lists, `star/rmats/bamlist/{CONDITION}_bamlist.txt` (treatment is group 1, control group 2); when a list is missing or its names are not the samples of the sample sheet, the report uses the order of the sample sheet and prints a warning. A contrast whose RDS is missing is reported as missing and the report still renders.

**File `02_rmats.Rmd`.**
````rmd
---
title: "rMATS splicing events"
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
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE, dev = "png", dpi = 150, fig.width = 7, fig.height = 5)
```

This report shows the rMATS skipped exon (SE), alternative 5' and 3' splice site (A5SS, A3SS), mutually exclusive exons (MXE) and retained intron (RI) events of each contrast, read from the tables that `01_import.Rmd` wrote to `r file.path(params$report_dir, "rds")`. Every effect is stated as **treatment minus control**: a positive effect is higher in the treatment group (more inclusion for rMATS and SUPPA2, higher usage for the exon and transcript usage tools). For an event, a positive dPSI means it is more included in treatment, a negative dPSI means it is more skipped in treatment. The effect is the harmonised column of the RDS and is used as it is.

An event is significant when its adjusted p-value (the FDR written by rMATS) is below `r params$fdr`, its absolute dPSI is at least `r params$dpsi`, and it is not flagged `low_support` (read support below `r params$min_reads`: for each group the mean over its samples of IJC plus SJC, the smaller of the two group means). Flagged events stay in the tables and plots, are counted, and are left out of the significant counts.

```{r functions}
suppressPackageStartupMessages(library(ggplot2))
RMATS_TYPES <- c("SE", "A5SS", "A3SS", "MXE", "RI")
DIRECTION_WORDING <- c(treatment_up = "more included in treatment", control_up = "more skipped in treatment")
direction_wording <- function(direction) unname(DIRECTION_WORDING[direction])
none_message <- function(fdr, dpsi) paste0("none at FDR ", format(fdr), ", |dPSI| ", format(dpsi))

# significant: padj below the FDR, absolute effect at least the dPSI threshold, and not flagged low_support
sig_flag <- function(tab, fdr, dpsi) {
  !is.na(tab$padj) & tab$padj < fdr & !is.na(tab$effect) & abs(tab$effect) >= dpsi & !tab$low_support
}

# one row per type: events, tested (padj present), flagged low_support, significant, and the significant ones split by direction
split_counts <- function(tab, fdr, dpsi, types) {
  sig <- sig_flag(tab, fdr, dpsi)
  rows <- lapply(types, function(ty) {
    i <- tab$feature_type == ty
    r <- data.frame(type = ty, events = sum(i), tested = sum(i & !is.na(tab$padj)), low_support = sum(i & tab$low_support),
                    significant = sum(i & sig), stringsAsFactors = FALSE)
    r[[direction_wording("treatment_up")]] <- sum(i & sig & tab$direction == "treatment_up")
    r[[direction_wording("control_up")]] <- sum(i & sig & tab$direction == "control_up")
    r
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# flagged events, and how many of them would be significant if they were not flagged
low_support_report <- function(tab, fdr, dpsi) {
  would <- sig_flag(transform(tab, low_support = FALSE), fdr, dpsi) & tab$low_support
  c(flagged = sum(tab$low_support), would_be_significant = sum(would))
}

# events not flagged and with an adjusted p-value, nearest first (smallest padj, then the larger change)
nearest_events <- function(tab, n = 3) {
  ok <- tab[!tab$low_support & !is.na(tab$padj), , drop = FALSE]
  ok <- ok[order(ok$padj, 1 - abs(ok$effect), ok$feature_id), , drop = FALSE]
  head(ok, n)
}
# the significant events first, then the nearest ones
top_events <- function(tab, fdr, dpsi, n) {
  ok <- tab[!tab$low_support & !is.na(tab$padj), , drop = FALSE]
  ok <- ok[order(!sig_flag(ok, fdr, dpsi), ok$padj, 1 - abs(ok$effect), ok$feature_id), , drop = FALSE]
  head(ok, n)
}
show_events <- function(d) {
  w <- direction_wording(d$direction)
  data.frame(gene = d$gene_symbol, event = d$feature_id, effect = round(d$effect, 3), padj = signif(d$padj, 3),
             direction = ifelse(is.na(w), "no change", w), stringsAsFactors = FALSE)
}

# a small reader for the rMATS tables (the column ID appears twice; the second one becomes ID_2) and for the sheets
strip_quotes <- function(x) gsub("^\"|\"$", "", x)
dedupe_names <- function(h) {
  k <- stats::ave(seq_along(h), h, FUN = seq_along)
  ifelse(k > 1, paste0(h, "_", k), h)
}
read_delim_chr <- function(path, sep) {
  if (!file.exists(path)) stop("file ", path, " not found: none of its columns can be read")
  l <- sub("\r$", "", readLines(path, warn = FALSE))
  l <- l[nzchar(l)]
  if (length(l) == 0) stop("file ", path, " is empty: no header line, so no column can be found")
  nm <- strsplit(l[1], sep, fixed = TRUE)[[1]]
  if (length(l) == 1) {
    return(as.data.frame(stats::setNames(replicate(length(nm), character(0), simplify = FALSE), nm),
                         check.names = FALSE, stringsAsFactors = FALSE))
  }
  tab <- utils::read.table(text = l[-1], sep = sep, quote = "", header = FALSE, colClasses = "character", na.strings = character(0),
                           comment.char = "", col.names = paste0("V", seq_along(nm)), fill = FALSE, stringsAsFactors = FALSE)
  names(tab) <- nm
  tab
}
need_cols <- function(tab, cols, path) {
  miss <- setdiff(cols, names(tab))
  if (length(miss)) stop("file ", path, ": column not found: ", paste(miss, collapse = ", "),
                         " (columns present: ", paste(names(tab), collapse = ", "), ")")
}
read_rmats_raw <- function(path) {
  tab <- read_delim_chr(path, "\t")
  names(tab) <- dedupe_names(names(tab))
  tab
}
read_contrasts <- function(path) {
  tab <- read_delim_chr(path, ",")
  need_cols(tab, c("contrast", "treatment", "control"), path)
  tab <- tab[nzchar(tab$contrast), c("contrast", "treatment", "control")]
  rownames(tab) <- NULL
  tab
}

# the directory of the rMATS tables of a contrast: its own, or the one with the suffix _paired (a run with rmats_paired_stats
# true); when both hold tables of the type, the one that matches the run's pairing switch comes first
find_rmats_dir <- function(outdir, contrast, paired, type) {
  dirs <- file.path(outdir, "star", "rmats", c(contrast, paste0(contrast, "_paired")))
  if (isTRUE(paired)) dirs <- rev(dirs)
  has <- vapply(dirs, function(d) any(file.exists(file.path(d, sprintf("%s.MATS.%s.txt", RMATS_TYPES, type)))), logical(1))
  if (any(has)) dirs[has][1] else NA_character_
}

# IncLevel1 and IncLevel2 (comma lists, one value per sample) of the events with these feature_ids (TYPE:ID), re-read from the tables
read_inclevels <- function(dir, feature_ids, rmats_type) {
  out <- list()
  for (ty in RMATS_TYPES) {
    want <- sub("^[^:]*:", "", feature_ids[startsWith(feature_ids, paste0(ty, ":"))])
    f <- file.path(dir, sprintf("%s.MATS.%s.txt", ty, rmats_type))
    if (length(want) == 0 || !file.exists(f)) next
    tab <- read_rmats_raw(f)
    need_cols(tab, c("ID", "ID_2", "IncLevel1", "IncLevel2"), f)
    k <- tab[match(want, tab$ID), , drop = FALSE]
    k <- k[!is.na(k$ID), , drop = FALSE]
    out[[ty]] <- data.frame(feature_id = paste0(ty, ":", k$ID), IncLevel1 = k$IncLevel1, IncLevel2 = k$IncLevel2, stringsAsFactors = FALSE)
  }
  if (length(out) == 0) return(data.frame(feature_id = character(0), IncLevel1 = character(0), IncLevel2 = character(0), stringsAsFactors = FALSE))
  res <- do.call(rbind, unname(out))
  idx <- match(feature_ids, res$feature_id)
  res[idx[!is.na(idx)], , drop = FALSE]
}

psi_values <- function(x) suppressWarnings(as.numeric(strsplit(x, ",", fixed = TRUE)[[1]]))
# rMATS lists its samples in the sample lists it was started with: star/rmats/bamlist/{condition}_bamlist.txt, BAM names separated by
# commas (the sample name followed by _sorted.bam). IncLevel1 and the counts of group 1 follow the list of the treatment (b1),
# IncLevel2 the list of the control (b2). The order of the sample sheet is not guaranteed to be the same.
read_bamlist <- function(outdir, condition) {
  f <- file.path(outdir, "star", "rmats", "bamlist", paste0(condition, "_bamlist.txt"))
  if (!file.exists(f)) return(NULL)
  x <- trimws(unlist(strsplit(readLines(f, warn = FALSE), ",", fixed = TRUE)))
  sub("(_sorted)?\\.bam$", "", basename(x[nzchar(x)]))
}
# the samples of a group in the order of its sample list; without a list, or when its names are not the sample names of the group
# in the sample sheet, the order of the sample sheet with a warning
group_samples <- function(ss, condition, bamlist) {
  sheet <- unique(ss$sample[ss$condition == condition])
  if (is.null(bamlist))
    return(list(samples = sheet, from = "sample sheet",
                warning = paste0("WARNING: no rMATS sample list found for ", condition, ": the order of the sample sheet is used, which may differ from the order in IncLevel1 and IncLevel2")))
  if (anyDuplicated(bamlist) || !setequal(bamlist, sheet))
    return(list(samples = sheet, from = "sample sheet",
                warning = paste0("WARNING: the names of the rMATS sample list of ", condition, " (", paste(bamlist, collapse = ", "),
                                 ") are not the samples of this condition in the sample sheet: the order of the sample sheet is used, which may differ from the order in IncLevel1 and IncLevel2")))
  list(samples = bamlist, from = "rMATS sample list", warning = NULL)
}
# one row per event and sample. IncLevel1 lists the samples of group 1 (the treatment), IncLevel2 those of group 2 (the control)
psi_long <- function(inc, ss, treatment, control, bam_treatment = NULL, bam_control = NULL) {
  g1 <- group_samples(ss, treatment, bam_treatment)
  g2 <- group_samples(ss, control, bam_control)
  s1 <- g1$samples
  s2 <- g2$samples
  rows <- lapply(seq_len(nrow(inc)), function(i) {
    v1 <- psi_values(inc$IncLevel1[i]); v2 <- psi_values(inc$IncLevel2[i])
    if (length(v1) != length(s1) || length(v2) != length(s2))
      stop("event ", inc$feature_id[i], ": IncLevel1 and IncLevel2 list ", length(v1), " and ", length(v2), " samples, the sample order has ",
           length(s1), " samples of ", treatment, " and ", length(s2), " samples of ", control)
    data.frame(feature_id = inc$feature_id[i], sample = c(s1, s2), condition = rep(c(treatment, control), c(length(s1), length(s2))),
               role = rep(c("treatment", "control"), c(length(s1), length(s2))), psi = c(v1, v2), stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, rows)
  attr(res, "order_from") <- c(treatment = g1$from, control = g2$from)
  attr(res, "order_warnings") <- c(g1$warning, g2$warning)
  res
}
psi_legend <- function(treatment, control) c(paste0("treatment (", treatment, ")"), paste0("control (", control, ")"))
psi_plot <- function(long, labels, treatment, control, contrast) {
  long$sample <- factor(long$sample, levels = unique(long$sample))
  long$role <- factor(long$role, levels = c("treatment", "control"), labels = psi_legend(treatment, control))
  gs <- labels$gene_symbol[match(long$feature_id, labels$feature_id)]
  long$event <- ifelse(is.na(gs) | !nzchar(gs), long$feature_id, paste0(gs, " (", long$feature_id, ")"))
  ggplot(long, aes(x = sample, y = psi, colour = role)) + geom_point(size = 2.5) +
    facet_wrap(~event, ncol = 3) +
    scale_colour_manual(values = c("#D55E00", "#0072B2"), name = "condition") +
    labs(title = paste0("Per-sample PSI of the top rMATS events, ", contrast), x = "sample", y = "PSI (IncLevel)") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# volcano: x = dPSI (treatment minus control), y = -log10 adjusted p-value, thresholds as dashed lines, the top 10 labelled
volcano_plot <- function(tab, type, fdr, dpsi, contrast, treatment, control) {
  d <- tab[!is.na(tab$padj) & !is.na(tab$effect), , drop = FALSE]
  d$nlp <- -log10(pmax(d$padj, 1e-300))
  d$status <- ifelse(d$low_support, "low support", ifelse(sig_flag(d, fdr, dpsi), "significant", "not significant"))
  lab <- top_events(d, fdr, dpsi, 10)
  lab$label <- ifelse(is.na(lab$gene_symbol) | !nzchar(lab$gene_symbol), lab$feature_id, lab$gene_symbol)
  lab$nlp <- -log10(pmax(lab$padj, 1e-300))
  p <- ggplot(d, aes(x = effect, y = nlp, colour = status)) + geom_point(alpha = 0.7, size = 1.6) +
    geom_vline(xintercept = c(-dpsi, dpsi), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(fdr), linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = c("significant" = "#D55E00", "not significant" = "grey60", "low support" = "#56B4E9")) +
    labs(title = paste0(type, ", ", contrast), x = paste0("dPSI (", treatment, " minus ", control, ")"), y = "-log10 adjusted p-value") +
    theme_bw()
  if (nrow(lab) > 0) {
    p <- if (requireNamespace("ggrepel", quietly = TRUE)) {
      p + ggrepel::geom_text_repel(data = lab, aes(label = label), colour = "black", size = 3, max.overlaps = Inf)
    } else {
      p + geom_text(data = lab, aes(label = label), colour = "black", size = 3, check_overlap = TRUE, vjust = -0.6)
    }
  }
  p
}
effect_hist <- function(tab, contrast, treatment, control) {
  d <- tab[!is.na(tab$effect), , drop = FALSE]
  ggplot(d, aes(x = effect)) + geom_histogram(bins = 40, fill = "grey55", colour = "white") +
    geom_vline(xintercept = 0, linetype = "dashed") + facet_wrap(~feature_type, scales = "free_y") +
    labs(title = paste0("dPSI of all rMATS events, ", contrast), x = paste0("dPSI (", treatment, " minus ", control, ")"), y = "events") +
    theme_bw()
}
save_plot <- function(p, file, width = 7, height = 5) {
  ggplot2::ggsave(filename = file, plot = p, width = width, height = height, dpi = 150)
  invisible(file)
}
show_image <- function(file, alt) cat("![", alt, "](", file, ")\n\n", sep = "")
```

```{r report, results="asis"}
modules <- strsplit(trimws(params$modules), "[[:space:]]+")[[1]]
fig_dir <- file.path(params$report_dir, "figures"); tab_dir <- file.path(params$report_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE); dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
contrasts <- read_contrasts(params$contrasts)
ss <- read_delim_chr(params$samplesheet, ",")
need_cols(ss, c("sample", "condition"), params$samplesheet)
if (!"rmats" %in% modules) cat("rMATS is not among the modules of this report (", params$modules, "): there is nothing to show.\n\n", sep = "")
for (i in seq_len(nrow(contrasts))) {
  if (!"rmats" %in% modules) break
  C <- contrasts$contrast[i]; Tr <- contrasts$treatment[i]; Co <- contrasts$control[i]
  cat("\n\n## Contrast ", C, "\n\n", sep = "")
  cat("Treatment `", Tr, "`, control `", Co, "`: every effect below is `", Tr, "` minus `", Co, "`. A positive dPSI: ",
      direction_wording("treatment_up"), "; a negative dPSI: ", direction_wording("control_up"), ".\n\n", sep = "")
  cat("Thresholds: FDR ", format(params$fdr), ", abs(dPSI) ", format(params$dpsi), "; rMATS minimum read support ", params$min_reads,
      "; rMATS counts ", params$rmats_type, ". rMATS statistics: ", if (isTRUE(params$paired_rmats)) "paired statistics" else "unpaired statistics",
      " (read from the run's parameters).\n\n", sep = "")
  f <- file.path(params$report_dir, "rds", paste0(C, ".rmats.rds"))
  if (!file.exists(f)) {
    cat("The file `", f, "` was not found: the import found no rMATS table for this contrast (the files it looked for are listed in `files_missing` of `rds/",
        C, ".qc.rds`). No rMATS result is shown for this contrast.\n\n", sep = "")
    qf <- file.path(params$report_dir, "rds", paste0(C, ".qc.rds"))
    if (file.exists(qf)) { qc <- readRDS(qf); if (length(qc$files_missing)) cat("Files not found:\n\n", paste0("- `", qc$files_missing, "`\n"), "\n", sep = "") }
    next
  }
  tab <- readRDS(f)
  need_cols(tab, c("gene_id", "gene_symbol", "feature_type", "feature_id", "effect", "padj", "direction", "n_support", "low_support"), f)
  lr <- low_support_report(tab, params$fdr, params$dpsi)
  cat("Events flagged low_support (read support below ", params$min_reads, "): ", lr[["flagged"]], "; of these ", lr[["would_be_significant"]],
      " would be significant at FDR ", format(params$fdr), ", |dPSI| ", format(params$dpsi), " if they were not flagged. Flagged events stay in the tables and plots, ",
      "are left out of the significant counts below, and keep the FDR that rMATS wrote.\n\n", sep = "")
  cn <- split_counts(tab, params$fdr, params$dpsi, RMATS_TYPES)
  utils::write.csv(cn, file.path(tab_dir, paste0(C, ".rmats.counts.csv")), row.names = FALSE)
  cat("Events per type and significant events split by direction (counts file `tables/", C, ".rmats.counts.csv`):\n\n", sep = "")
  cat(knitr::kable(cn, format = "markdown"), sep = "\n")
  cat("\n\n")
  if (!isTRUE(params$events)) {
    cat("The event tables and volcano plots are switched off for this report (events = FALSE).\n\n")
  } else {
    for (ty in RMATS_TYPES) {
      cat("\n\n### ", ty, ", ", C, "\n\n", sep = "")
      st <- tab[tab$feature_type == ty, , drop = FALSE]
      sig <- st[sig_flag(st, params$fdr, params$dpsi), , drop = FALSE]
      sig <- sig[order(sig$padj, 1 - abs(sig$effect), sig$feature_id), , drop = FALSE]
      sig$direction_wording <- direction_wording(sig$direction)
      utils::write.csv(sig[, c("contrast", "gene_id", "gene_symbol", "feature_type", "feature_id", "effect", "padj", "direction", "direction_wording", "n_support")],
                       file.path(tab_dir, sprintf("%s.rmats.%s.significant.csv", C, ty)), row.names = FALSE)
      if (nrow(st) == 0) { cat("No ", ty, " event in the table.\n\n", sep = ""); next }
      cat(nrow(sig), " significant ", ty, " events: ", sum(sig$direction == "treatment_up"), " ", direction_wording("treatment_up"), ", ",
          sum(sig$direction == "control_up"), " ", direction_wording("control_up"), " (table `tables/", C, ".rmats.", ty, ".significant.csv`).\n\n", sep = "")
      if (nrow(sig) == 0) {
        cat("**", none_message(params$fdr, params$dpsi), "**\n\n", sep = "")
        nr <- nearest_events(st, 3)
        cat("The three nearest events by adjusted p-value (events flagged low_support are not listed):\n\n")
        if (nrow(nr) > 0) cat(knitr::kable(show_events(nr), format = "markdown"), sep = "\n") else cat("No unflagged event of this type has an adjusted p-value.\n")
        cat("\n\n")
      } else {
        cat(knitr::kable(show_events(head(sig, 20)), format = "markdown"), sep = "\n")
        cat("\n\n")
      }
      if (any(!is.na(st$padj) & !is.na(st$effect))) {
        vf <- file.path(fig_dir, sprintf("%s.rmats.volcano.%s.png", C, ty))
        save_plot(volcano_plot(st, ty, params$fdr, params$dpsi, C, Tr, Co), vf)
        show_image(vf, paste0("Volcano plot of the ", ty, " events, ", C))
      } else cat("No ", ty, " event has both an effect and an adjusted p-value: no volcano plot.\n\n", sep = "")
    }
    cat("\n\n### Distribution of the effects, ", C, "\n\n", sep = "")
    hf <- file.path(fig_dir, paste0(C, ".rmats.effect_hist.png"))
    save_plot(effect_hist(tab, C, Tr, Co), hf)
    show_image(hf, paste0("Histogram of the dPSI of all events, ", C))
  }
  if (!isTRUE(params$topgene)) {
    cat("The top-gene plots are switched off for this report (topgene = FALSE).\n\n")
  } else {
    cat("\n\n### Per-sample PSI of the top events, ", C, "\n\n", sep = "")
    top <- top_events(tab, params$fdr, params$dpsi, 6)
    rdir <- find_rmats_dir(params$outdir, C, params$paired_rmats, params$rmats_type)
    if (nrow(top) == 0) {
      cat("No unflagged event with an adjusted p-value: no PSI plot.\n\n")
    } else if (is.na(rdir)) {
      cat("No rMATS table of type ", params$rmats_type, " was found for this contrast: no PSI plot.\n\n", sep = "")
    } else {
      inc <- read_inclevels(rdir, top$feature_id, params$rmats_type)
      bam_t <- read_bamlist(params$outdir, Tr); bam_c <- read_bamlist(params$outdir, Co)
      psi <- psi_long(inc, ss, Tr, Co, bam_t, bam_c)
      cat("The six events that are significant, or nearest to significance when fewer than six are (events flagged low_support are left out). ",
          "IncLevel1 holds the samples of the treatment `", Tr, "` and IncLevel2 those of the control `", Co,
          "`, read again from the tables in `", rdir, "`. The samples are in the order of the sample lists rMATS was started with (`star/rmats/bamlist/",
          Tr, "_bamlist.txt` and `", Co, "_bamlist.txt`).\n\n", sep = "")
      for (w in attr(psi, "order_warnings")) cat("**", w, "**\n\n", sep = "")
      pf <- file.path(fig_dir, paste0(C, ".rmats.psi_top.png"))
      save_plot(psi_plot(psi, top[, c("feature_id", "gene_symbol")], Tr, Co, C), pf, width = 9, height = 6)
      show_image(pf, paste0("Per-sample PSI of the top events, ", C))
    }
  }
}
```
````

The SUPPA2 report `03_suppa.Rmd` reads `rds/{CONTRAST}.suppa.rds` for every contrast. The RDS holds the local events and, when the isoform file exists, the isoforms (`feature_type` is `isoform`); the event type of a local event is read from its event identifier `GENE;TYPE:chr:...`, split at the first `;` only (`SE`, `A5`, `A3`, `MX`, `RI`, `AF`, `AL`). It writes, per contrast, `tables/{CONTRAST}.suppa.counts.csv`, `tables/{CONTRAST}.suppa.{type}.significant.csv`, `figures/{CONTRAST}.suppa.volcano.{type}.png`, `figures/{CONTRAST}.suppa.effect_hist.png` and `figures/{CONTRAST}.suppa.top.png` (`type` is an event type or `isoform`). The adjusted p-value is the Benjamini-Hochberg value computed by `01_import.Rmd`. The SUPPA2 difference files hold dPSI and a p-value for each event, not the PSI of each sample, so the top-event plot shows the dPSI of the top events. A contrast whose RDS is missing is reported as missing and the report still renders.

**File `03_suppa.Rmd`.**
````rmd
---
title: "SUPPA2 differential splicing"
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
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE, dev = "png", dpi = 150, fig.width = 7, fig.height = 5)
```

This report shows the SUPPA2 local events (skipping exon SE, alternative 5' and 3' splice sites A5 and A3, mutually exclusive exons MX, retained intron RI, alternative first and last exons AF and AL) and the transcript isoforms of each contrast, read from the tables that `01_import.Rmd` wrote to `r file.path(params$report_dir, "rds")`. Every effect is stated as **treatment minus control**: a positive effect is higher in the treatment group (more inclusion for rMATS and SUPPA2, higher usage for the exon and transcript usage tools). For an event, a positive dPSI means it is more included in treatment, a negative dPSI means it is more skipped in treatment. The effect is the harmonised column of the RDS and is used as it is.

SUPPA2 reports a p-value and no adjusted p-value. The adjusted p-value used here is the Benjamini-Hochberg correction that `01_import.Rmd` applied, separately for each contrast and for each file (local events, isoforms). An event whose dPSI is `nan` was not tested: it is kept in the table with a missing effect, counted as untested, left out of the correction and not drawn. An event is significant when its adjusted p-value is below `r params$fdr` and its absolute dPSI is at least `r params$dpsi`.

```{r functions}
suppressPackageStartupMessages(library(ggplot2))
SUPPA_TYPES <- c("SE", "A5", "A3", "MX", "RI", "AF", "AL", "isoform")
DIRECTION_WORDING <- c(treatment_up = "more included in treatment", control_up = "more skipped in treatment")
direction_wording <- function(direction) unname(DIRECTION_WORDING[direction])
none_message <- function(fdr, dpsi) paste0("none at FDR ", format(fdr), ", |dPSI| ", format(dpsi))

# the event identifier of SUPPA2 is GENE;TYPE:chr:...: the type is the text before the first : after the FIRST ; (a later ; stays in the event)
suppa_event_type <- function(id) {
  ev <- sub("^[^;]*;", "", id)
  sub(":.*$", "", ev)
}
# the group of a row: isoform for a row of the isoform file, else the event type read from the identifier
suppa_group <- function(tab) as.character(ifelse(tab$feature_type == "isoform", "isoform", suppa_event_type(tab$feature_id)))

# significant: padj below the FDR and absolute effect at least the dPSI threshold (untested rows have no padj and no effect)
sig_flag <- function(tab, fdr, dpsi) {
  !is.na(tab$padj) & tab$padj < fdr & !is.na(tab$effect) & abs(tab$effect) >= dpsi
}

# one row per group: events, untested (no effect), tested (padj present), significant, and the significant ones split by direction
split_counts <- function(tab, fdr, dpsi, types) {
  sig <- sig_flag(tab, fdr, dpsi)
  g <- suppa_group(tab)
  rows <- lapply(types, function(ty) {
    i <- g == ty
    r <- data.frame(type = ty, events = sum(i), untested = sum(i & is.na(tab$effect)), tested = sum(i & !is.na(tab$padj)),
                    significant = sum(i & sig), stringsAsFactors = FALSE)
    r[[direction_wording("treatment_up")]] <- sum(i & sig & tab$direction == "treatment_up")
    r[[direction_wording("control_up")]] <- sum(i & sig & tab$direction == "control_up")
    r
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# tested events (with an adjusted p-value), nearest first (smallest padj, then the larger change)
nearest_events <- function(tab, n = 3) {
  ok <- tab[!is.na(tab$padj), , drop = FALSE]
  ok <- ok[order(ok$padj, 1 - abs(ok$effect), ok$feature_id), , drop = FALSE]
  head(ok, n)
}
# the significant events first, then the nearest ones
top_events <- function(tab, fdr, dpsi, n) {
  ok <- tab[!is.na(tab$padj), , drop = FALSE]
  ok <- ok[order(!sig_flag(ok, fdr, dpsi), ok$padj, 1 - abs(ok$effect), ok$feature_id), , drop = FALSE]
  head(ok, n)
}
show_events <- function(d) {
  w <- direction_wording(d$direction)
  data.frame(gene = d$gene_id, event = d$feature_id, effect = round(d$effect, 3), padj = signif(d$padj, 3),
             direction = ifelse(is.na(w), "no change", w), stringsAsFactors = FALSE)
}

strip_quotes <- function(x) gsub("^\"|\"$", "", x)
read_delim_chr <- function(path, sep) {
  if (!file.exists(path)) stop("file ", path, " not found: none of its columns can be read")
  l <- sub("\r$", "", readLines(path, warn = FALSE))
  l <- l[nzchar(l)]
  if (length(l) == 0) stop("file ", path, " is empty: no header line, so no column can be found")
  nm <- strsplit(l[1], sep, fixed = TRUE)[[1]]
  if (length(l) == 1) {
    return(as.data.frame(stats::setNames(replicate(length(nm), character(0), simplify = FALSE), nm),
                         check.names = FALSE, stringsAsFactors = FALSE))
  }
  tab <- utils::read.table(text = l[-1], sep = sep, quote = "", header = FALSE, colClasses = "character", na.strings = character(0),
                           comment.char = "", col.names = paste0("V", seq_along(nm)), fill = FALSE, stringsAsFactors = FALSE)
  names(tab) <- nm
  tab
}
need_cols <- function(tab, cols, path) {
  miss <- setdiff(cols, names(tab))
  if (length(miss)) stop("file ", path, ": column not found: ", paste(miss, collapse = ", "),
                         " (columns present: ", paste(names(tab), collapse = ", "), ")")
}
read_contrasts <- function(path) {
  tab <- read_delim_chr(path, ",")
  need_cols(tab, c("contrast", "treatment", "control"), path)
  tab <- tab[nzchar(tab$contrast), c("contrast", "treatment", "control")]
  rownames(tab) <- NULL
  tab
}

# volcano: x = dPSI (treatment minus control), y = -log10 adjusted p-value, thresholds as dashed lines, the top 10 labelled
volcano_plot <- function(tab, type, fdr, dpsi, contrast, treatment, control) {
  d <- tab[!is.na(tab$padj) & !is.na(tab$effect), , drop = FALSE]
  d$nlp <- -log10(pmax(d$padj, 1e-300))
  d$status <- ifelse(sig_flag(d, fdr, dpsi), "significant", "not significant")
  lab <- top_events(d, fdr, dpsi, 10)
  lab$label <- lab$gene_id
  lab$nlp <- -log10(pmax(lab$padj, 1e-300))
  p <- ggplot(d, aes(x = effect, y = nlp, colour = status)) + geom_point(alpha = 0.7, size = 1.6) +
    geom_vline(xintercept = c(-dpsi, dpsi), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(fdr), linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = c("significant" = "#D55E00", "not significant" = "grey60")) +
    labs(title = paste0("SUPPA2 ", type, ", ", contrast), x = paste0("dPSI (", treatment, " minus ", control, ")"), y = "-log10 adjusted p-value") +
    theme_bw()
  if (nrow(lab) > 0) {
    p <- if (requireNamespace("ggrepel", quietly = TRUE)) {
      p + ggrepel::geom_text_repel(data = lab, aes(label = label), colour = "black", size = 3, max.overlaps = Inf)
    } else {
      p + geom_text(data = lab, aes(label = label), colour = "black", size = 3, check_overlap = TRUE, vjust = -0.6)
    }
  }
  p
}
effect_hist <- function(tab, contrast, treatment, control) {
  d <- tab[!is.na(tab$effect), , drop = FALSE]
  d$group <- suppa_group(d)
  ggplot(d, aes(x = effect)) + geom_histogram(bins = 40, fill = "grey55", colour = "white") +
    geom_vline(xintercept = 0, linetype = "dashed") + facet_wrap(~group, scales = "free_y") +
    labs(title = paste0("dPSI of all tested SUPPA2 events, ", contrast), x = paste0("dPSI (", treatment, " minus ", control, ")"), y = "events") +
    theme_bw()
}
# the dPSI of the top events (the SUPPA2 difference files hold no PSI per sample)
top_plot <- function(top, treatment, control, contrast) {
  top$label <- make.unique(paste0(top$gene_id, " ", suppa_group(top)))
  top$label <- factor(top$label, levels = rev(top$label))
  w <- direction_wording(top$direction)
  top$wording <- factor(ifelse(is.na(w), "no change", w), levels = c(unname(DIRECTION_WORDING), "no change"))
  ggplot(top, aes(x = effect, y = label, fill = wording)) + geom_col() +
    scale_fill_manual(values = c("#D55E00", "#0072B2", "grey60"), drop = FALSE, name = "direction") +
    labs(title = paste0("dPSI of the top SUPPA2 events, ", contrast), x = paste0("dPSI (", treatment, " minus ", control, ")"), y = NULL) +
    theme_bw()
}
save_plot <- function(p, file, width = 7, height = 5) {
  ggplot2::ggsave(filename = file, plot = p, width = width, height = height, dpi = 150)
  invisible(file)
}
show_image <- function(file, alt) cat("![", alt, "](", file, ")\n\n", sep = "")
```

```{r report, results="asis"}
modules <- strsplit(trimws(params$modules), "[[:space:]]+")[[1]]
fig_dir <- file.path(params$report_dir, "figures"); tab_dir <- file.path(params$report_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE); dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
contrasts <- read_contrasts(params$contrasts)
if (!"suppa" %in% modules) cat("SUPPA2 is not among the modules of this report (", params$modules, "): there is nothing to show.\n\n", sep = "")
for (i in seq_len(nrow(contrasts))) {
  if (!"suppa" %in% modules) break
  C <- contrasts$contrast[i]; Tr <- contrasts$treatment[i]; Co <- contrasts$control[i]
  cat("\n\n## Contrast ", C, "\n\n", sep = "")
  cat("Treatment `", Tr, "`, control `", Co, "`: every effect below is `", Tr, "` minus `", Co, "`. A positive dPSI: ",
      direction_wording("treatment_up"), "; a negative dPSI: ", direction_wording("control_up"), ".\n\n", sep = "")
  cat("Thresholds: FDR ", format(params$fdr), ", abs(dPSI) ", format(params$dpsi), ". SUPPA2 statistics: ",
      if (isTRUE(params$paired_suppa)) "paired statistics" else "unpaired statistics",
      " (read from the run's parameters). The adjusted p-value is Benjamini-Hochberg, computed by the import for each file.\n\n", sep = "")
  f <- file.path(params$report_dir, "rds", paste0(C, ".suppa.rds"))
  if (!file.exists(f)) {
    cat("The file `", f, "` was not found: the import found no SUPPA2 table for this contrast (the files it looked for are listed in `files_missing` of `rds/",
        C, ".qc.rds`). No SUPPA2 result is shown for this contrast.\n\n", sep = "")
    qf <- file.path(params$report_dir, "rds", paste0(C, ".qc.rds"))
    if (file.exists(qf)) { qc <- readRDS(qf); if (length(qc$files_missing)) cat("Files not found:\n\n", paste0("- `", qc$files_missing, "`\n"), "\n", sep = "") }
    next
  }
  tab <- readRDS(f)
  need_cols(tab, c("gene_id", "feature_type", "feature_id", "effect", "padj", "direction"), f)
  tab$group <- suppa_group(tab)
  types <- c(SUPPA_TYPES, setdiff(unique(tab$group), SUPPA_TYPES))
  if (!any(tab$group == "isoform")) cat("No isoform table was read for this contrast (`feature_type` isoform is absent from the RDS): the isoform section is empty.\n\n")
  cat("Rows of the table: ", nrow(tab), "; not tested (dPSI nan, no effect and no adjusted p-value): ", sum(is.na(tab$effect)), ".\n\n", sep = "")
  cn <- split_counts(tab, params$fdr, params$dpsi, types)
  utils::write.csv(cn, file.path(tab_dir, paste0(C, ".suppa.counts.csv")), row.names = FALSE)
  cat("Events per type and significant events split by direction (counts file `tables/", C, ".suppa.counts.csv`):\n\n", sep = "")
  cat(knitr::kable(cn, format = "markdown"), sep = "\n")
  cat("\n\n")
  if (!isTRUE(params$events)) {
    cat("The event tables and volcano plots are switched off for this report (events = FALSE).\n\n")
  } else {
    for (ty in types) {
      cat("\n\n### ", ty, ", ", C, "\n\n", sep = "")
      st <- tab[tab$group == ty, , drop = FALSE]
      sig <- st[sig_flag(st, params$fdr, params$dpsi), , drop = FALSE]
      sig <- sig[order(sig$padj, 1 - abs(sig$effect), sig$feature_id), , drop = FALSE]
      sig$direction_wording <- direction_wording(sig$direction)
      utils::write.csv(sig[, c("contrast", "gene_id", "feature_type", "feature_id", "effect", "padj", "direction", "direction_wording")],
                       file.path(tab_dir, sprintf("%s.suppa.%s.significant.csv", C, ty)), row.names = FALSE)
      if (nrow(st) == 0) { cat("No ", ty, " row in the table.\n\n", sep = ""); next }
      cat(nrow(sig), " significant ", ty, " rows: ", sum(sig$direction == "treatment_up"), " ", direction_wording("treatment_up"), ", ",
          sum(sig$direction == "control_up"), " ", direction_wording("control_up"), " (table `tables/", C, ".suppa.", ty, ".significant.csv`).\n\n", sep = "")
      if (nrow(sig) == 0) {
        cat("**", none_message(params$fdr, params$dpsi), "**\n\n", sep = "")
        nr <- nearest_events(st, 3)
        cat("The three nearest events by adjusted p-value:\n\n")
        if (nrow(nr) > 0) cat(knitr::kable(show_events(nr), format = "markdown"), sep = "\n") else cat("No row of this group was tested (every dPSI is nan).\n")
        cat("\n\n")
      } else {
        cat(knitr::kable(show_events(head(sig, 20)), format = "markdown"), sep = "\n")
        cat("\n\n")
      }
      if (any(!is.na(st$padj) & !is.na(st$effect))) {
        vf <- file.path(fig_dir, sprintf("%s.suppa.volcano.%s.png", C, ty))
        save_plot(volcano_plot(st, ty, params$fdr, params$dpsi, C, Tr, Co), vf)
        show_image(vf, paste0("Volcano plot of the SUPPA2 ", ty, " rows, ", C))
      } else cat("No ", ty, " row has both an effect and an adjusted p-value: no volcano plot.\n\n", sep = "")
    }
    cat("\n\n### Distribution of the effects, ", C, "\n\n", sep = "")
    if (any(!is.na(tab$effect))) {
      hf <- file.path(fig_dir, paste0(C, ".suppa.effect_hist.png"))
      save_plot(effect_hist(tab, C, Tr, Co), hf)
      show_image(hf, paste0("Histogram of the dPSI of all tested events, ", C))
    } else cat("No event was tested: no histogram.\n\n")
  }
  if (!isTRUE(params$topgene)) {
    cat("The top-gene plots are switched off for this report (topgene = FALSE).\n\n")
  } else {
    cat("\n\n### dPSI of the top events, ", C, "\n\n", sep = "")
    top <- top_events(tab, params$fdr, params$dpsi, 6)
    if (nrow(top) == 0) {
      cat("No event was tested: no plot.\n\n")
    } else {
      cat("The six events that are significant, or nearest to significance when fewer than six are. The SUPPA2 difference files hold the dPSI and the p-value of each event, not the PSI of each sample.\n\n")
      pf <- file.path(fig_dir, paste0(C, ".suppa.top.png"))
      save_plot(top_plot(top, Tr, Co, C), pf, width = 9, height = 5)
      show_image(pf, paste0("dPSI of the top SUPPA2 events, ", C))
    }
  }
}
```
````

The exon and transcript usage report `04_usage.Rmd` reads `rds/{CONTRAST}.exon_dexseq.rds`, `rds/{CONTRAST}.exon_edger.rds` and `rds/{CONTRAST}.dtu.rds` (written by `01_import.Rmd`) for every contrast and the gene-level files of the run, which the RDS does not hold and which the report reads again (read only): `star/dexseq_exon/results/perGeneQValue.{CONTRAST}.csv` (columns `groupID`, `padj`), `star/edger/contrast_{CONTRAST}.usage.simes.csv` and `.usage.gene.csv` (columns `Geneid`, `FDR`), `salmon/dexseq_dtu/results/stager/getAdjustedPValues.{T}-{C}.tsv` (columns `geneID`, `txID`, `gene`, `transcript`), the stageR object `stageRTx.{T}-{C}.rds` next to it and `salmon/dexseq_dtu/results/dexseq/DEXSeqResults.{T}-{C}.tsv` (for the two group columns, named by the treatment and the control). It writes, per contrast, `tables/{CONTRAST}.usage.counts.csv`, `tables/{CONTRAST}.exon_dexseq.genes.csv` and `.exon_dexseq.bins.csv`, `tables/{CONTRAST}.exon_edger.genes.csv` and `.exon_edger.exons.csv`, `tables/{CONTRAST}.dtu.genes.csv` and `.dtu.transcripts.csv`, `figures/{CONTRAST}.usage.volcano.{exon_dexseq,exon_edger,dtu}.png`, `figures/{CONTRAST}.usage.exon_top.png` and `figures/{CONTRAST}.usage.dtu_top.png`. Exon bins, exons, genes and transcripts are significant when their adjusted p-value is below `{FDR}`; `{DPSI}` is not applied. For DEXSeq exon usage a gene is significant when its `perGeneQValue` is below `{FDR}`; for edgeR exon usage when its Simes FDR is below `{FDR}`; for DTU, stageR's `gene` column (the screened adjusted p-value at the gene level) is compared with `{FDR}` for the genes and its `transcript` column (the screened adjusted p-value at the transcript level, `NA` where the gene was not screened in) for the transcripts. stageR's own overall alpha is the run's setting: the report prints it when the `stageRTx` object is readable with the stageR package of the image and otherwise says that it is not recorded in the table. stageR's adjusted p-values are valid only at the alpha they were computed for: when the report's FDR differs from the alpha of the `stageRTx` object (read with `getAlpha`), the stage-wise adjustment is re-run at the report's FDR from that object with the stageR package of the image, and when the object cannot be read the values of the run are used with a printed warning. The DEXSeq files of the DTU run hold the transcript counts (`countData.N`); with the sample annotation of `salmon/dexseq_dtu/results/dexseq/DEXSeqDataSet.{T}-{C}.rds` (read only, the sample and condition of each column) the isoform plot shows the proportion of each transcript in its gene for each sample (`tables/{CONTRAST}.usage.dtu_proportions.csv`), and without that object it shows the two group columns, labelled as DEXSeq usage coefficients that are not proportions. `params$events` switches the tables and the volcano plots, `params$topgene` the top-gene plots. A tool whose RDS is missing is reported as missing and the report still renders.

**File `04_usage.Rmd`.**
````rmd
---
title: "Exon and transcript usage"
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
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE, dev = "png", dpi = 150, fig.width = 7, fig.height = 5)
```

This report shows the exon usage results of DEXSeq and edgeR and the transcript usage (DTU) results of DEXSeq with stageR screening for each contrast, read from the tables that `01_import.Rmd` wrote to `r file.path(params$report_dir, "rds")` and from the gene-level files of the run (read only). Every effect is stated as **treatment minus control**: a positive effect is higher in the treatment group (more inclusion for rMATS and SUPPA2, higher usage for the exon and transcript usage tools). For an exon bin, an exon or a transcript, a positive log2 fold change means it is more used in treatment and a negative one that it is more used in control. The effect is the harmonised column of the RDS and is used as it is.

Exon usage and transcript usage are cut by the adjusted p-value only: a gene, an exon bin or a transcript is significant when its adjusted p-value is below `r params$fdr` (the dPSI threshold of the rMATS and SUPPA2 reports is not applied here). DEXSeq exon usage: a gene is significant when its per-gene q-value (`perGeneQValue`, read from the run's file) is below the FDR, an exon bin when its adjusted p-value is. edgeR exon usage: a gene is significant when its Simes FDR is below the FDR, an exon when its exon-level FDR is; the F-test gene table is counted as well. DTU: stageR screening. In the stageR table the column `gene` holds the screened adjusted p-value at the gene level and the column `transcript` the screened adjusted p-value at the transcript level (`NA` where the gene was not screened in); a gene is significant when `gene` is below the FDR and a transcript when `transcript` is. stageR's adjusted p-values are valid only at the alpha they were computed for: when the FDR of this report differs from the alpha of the run (read from the `stageRTx` object with `getAlpha`), the stage-wise adjustment is re-run at the FDR of this report from that object, and when the object cannot be read the values of the run are used with a warning. The direction of a gene is the sign of the effect of its most significant feature (smallest adjusted p-value among the significant bins, exons or transcripts of the gene); it is empty when the gene is significant but none of its features is.

```{r functions}
suppressPackageStartupMessages(library(ggplot2))
USAGE_WORDING <- c(treatment_up = "more used in treatment", control_up = "more used in control")
usage_wording <- function(direction) unname(USAGE_WORDING[direction])
none_message <- function(fdr) paste0("none at FDR ", format(fdr))

# significant: adjusted p-value present and below the FDR
sig_usage <- function(tab, fdr) !is.na(tab$padj) & tab$padj < fdr

# per gene: the significant features, the one with the smallest adjusted p-value (then the larger change, then the identifier) and its direction
top_by_gene <- function(feat, sig) {
  s <- feat[sig & !is.na(feat$padj), c("gene_id", "feature_id", "padj", "effect", "direction"), drop = FALSE]
  s <- s[order(s$padj, 1 - abs(s$effect), s$feature_id), , drop = FALSE]
  n <- table(s$gene_id)
  top <- s[!duplicated(s$gene_id), , drop = FALSE]
  data.frame(gene_id = top$gene_id, n_sig_features = as.integer(n[top$gene_id]), top_feature = top$feature_id, top_effect = top$effect,
             top_padj = top$padj, direction = ifelse(top$direction %in% c("treatment_up", "control_up"), top$direction, NA_character_),
             stringsAsFactors = FALSE)
}
# the significant genes of a gene-level table (gene_id, padj), with their top feature and direction
gene_summary <- function(gene_tab, feat, sig_feat, fdr) {
  g <- gene_tab[sig_usage(gene_tab, fdr), , drop = FALSE]
  g <- g[order(g$padj, g$gene_id), , drop = FALSE]
  tp <- top_by_gene(feat, sig_feat)
  k <- match(g$gene_id, tp$gene_id)
  out <- data.frame(gene_id = g$gene_id, gene_padj = g$padj, n_sig_features = tp$n_sig_features[k], top_feature = tp$top_feature[k],
                    top_effect = tp$top_effect[k], top_padj = tp$top_padj[k], direction = tp$direction[k], stringsAsFactors = FALSE)
  out$n_sig_features[is.na(out$n_sig_features)] <- 0L
  out$direction_wording <- usage_wording(out$direction)
  out
}
# the three rows with the smallest adjusted p-value
nearest_rows <- function(tab, n = 3) {
  ok <- tab[!is.na(tab$padj), , drop = FALSE]
  head(ok[order(ok$padj, ok$gene_id), , drop = FALSE], n)
}
# one row of the counts table: tested (adjusted p-value present), significant, and the significant ones by direction (NA at gene level)
usage_counts <- function(tool, level, p, sig, direction = NULL) {
  r <- data.frame(tool = tool, level = level, tested = sum(!is.na(p)), significant = sum(sig), stringsAsFactors = FALSE)
  r[[usage_wording("treatment_up")]] <- if (is.null(direction)) NA_integer_ else sum(sig & direction %in% "treatment_up")
  r[[usage_wording("control_up")]] <- if (is.null(direction)) NA_integer_ else sum(sig & direction %in% "control_up")
  r
}
# the exon number of a DEXSeq bin identifier (GENE:E001); a merged group GENE1+GENE2 has the same form
bin_number <- function(feature_id) as.integer(sub("^.*:E", "", feature_id))

strip_quotes <- function(x) gsub("^\"|\"$", "", x)
read_delim_chr <- function(path, sep) {
  if (!file.exists(path)) stop("file ", path, " not found: none of its columns can be read")
  l <- sub("\r$", "", readLines(path, warn = FALSE))
  l <- l[nzchar(l)]
  if (length(l) == 0) stop("file ", path, " is empty: no header line, so no column can be found")
  nm <- strsplit(l[1], sep, fixed = TRUE)[[1]]
  if (length(l) == 1) {
    return(as.data.frame(stats::setNames(replicate(length(nm), character(0), simplify = FALSE), nm),
                         check.names = FALSE, stringsAsFactors = FALSE))
  }
  tab <- utils::read.table(text = l[-1], sep = sep, quote = "", header = FALSE, colClasses = "character", na.strings = character(0),
                           comment.char = "", col.names = paste0("V", seq_along(nm)), fill = FALSE, stringsAsFactors = FALSE)
  names(tab) <- nm
  tab
}
need_cols <- function(tab, cols, path) {
  miss <- setdiff(cols, names(tab))
  if (length(miss)) stop("file ", path, ": column not found: ", paste(miss, collapse = ", "),
                         " (columns present: ", paste(names(tab), collapse = ", "), ")")
}
to_num <- function(x, path, col) {
  x <- trimws(strip_quotes(x))
  miss <- x %in% c("", "NA", "nan", "NaN", "-nan")
  v <- suppressWarnings(as.numeric(ifelse(miss, NA, x)))
  bad <- is.na(v) & !miss
  if (any(bad)) stop("file ", path, ": column ", col, " has a value that is not a number: '", x[bad][1], "'")
  v
}
read_contrasts <- function(path) {
  tab <- read_delim_chr(path, ",")
  need_cols(tab, c("contrast", "treatment", "control"), path)
  tab <- tab[nzchar(tab$contrast), c("contrast", "treatment", "control")]
  rownames(tab) <- NULL
  tab
}

# a gene-level table: the gene identifier and the adjusted p-value (DEXSeq perGeneQValue: groupID, padj; edgeR: Geneid, FDR)
read_gene_table <- function(path, sep, id_col, p_col) {
  tab <- read_delim_chr(path, sep)
  need_cols(tab, c(id_col, p_col), path)
  data.frame(gene_id = strip_quotes(tab[[id_col]]), padj = to_num(tab[[p_col]], path, p_col), stringsAsFactors = FALSE)
}
# stageR: gene and transcript are the screened adjusted p-values at the two levels, NA where not screened
read_stager <- function(path) {
  tab <- read_delim_chr(path, "\t")
  need_cols(tab, c("geneID", "txID", "gene", "transcript"), path)
  data.frame(gene_id = strip_quotes(tab$geneID), tx_id = strip_quotes(tab$txID), gene = to_num(tab$gene, path, "gene"),
             transcript = to_num(tab$transcript, path, "transcript"), stringsAsFactors = FALSE)
}
# one row per gene of the stageR table (the gene-level value is repeated on every transcript row)
stager_genes <- function(st) {
  g <- st[!duplicated(st$gene_id), , drop = FALSE]
  data.frame(gene_id = g$gene_id, padj = g$gene, stringsAsFactors = FALSE)
}
# an Ensembl gene identifier loses its version suffix (.5) only when it has one; a transcript identifier loses a trailing version in the same way
# (the pipeline strips the versions of the stageR identifiers, not those of the DEXSeq table: both sides of a join are stripped here)
strip_version <- function(x) sub("^(ENS[A-Z]*G[0-9]+)\\.[0-9]+$", "\\1", x)
strip_tx_version <- function(x) sub("\\.[0-9]+$", "", x)
# the overall alpha stageR was run with, read from the stageRTx object with getAlpha; NA when the file is absent or cannot be read
stager_alpha <- function(path) {
  if (!file.exists(path) || !requireNamespace("stageR", quietly = TRUE)) return(NA_real_)
  tryCatch(as.numeric(stageR::getAlpha(readRDS(path)))[1], error = function(e) NA_real_)
}
# stageR's adjusted p-values are valid only at the alpha they were computed for. At another FDR the stage-wise adjustment is re-run
# from the stageRTx object of the run (stageR is in the image); when the object cannot be read the values of the run are used as they stand
stager_for_fdr <- function(st, path, fdr) {
  alpha <- stager_alpha(path)
  if (is.na(alpha)) return(list(st = st, alpha = NA_real_, mode = "unreadable"))
  if (isTRUE(all.equal(alpha, fdr))) return(list(st = st, alpha = alpha, mode = "run"))
  re <- tryCatch({
    obj <- stageR::stageWiseAdjustment(readRDS(path), method = "dtu", alpha = fdr, allowNA = TRUE)
    r <- stageR::getAdjustedPValues(obj, onlySignificantGenes = FALSE, order = FALSE)
    data.frame(gene_id = as.character(r$geneID), tx_id = as.character(r$txID), gene = as.numeric(r$gene), transcript = as.numeric(r$transcript),
               stringsAsFactors = FALSE)
  }, error = function(e) NULL)
  if (is.null(re)) return(list(st = st, alpha = alpha, mode = "unreadable"))
  list(st = re, alpha = alpha, mode = "rerun")
}
# the sentence that tells which stageR values are used
stager_note <- function(res, fdr) {
  if (res$mode == "run") return(paste0("stageR was run with an overall alpha of ", format(res$alpha), " (read from the `stageRTx` object of the run), equal to the FDR of this report: the stageR table of the run is used as it is."))
  if (res$mode == "rerun") return(paste0("stageR's adjusted p-values are valid only at the alpha they were computed for. The run used an overall alpha of ", format(res$alpha),
    " (read from the `stageRTx` object of the run) and the FDR of this report is ", format(fdr), ", so the stage-wise adjustment was re-run at ", format(fdr),
    " from that object with the stageR package of the image; the values of the run are not used."))
  paste0("**WARNING: the stageR values of the run are valid only at the alpha the run used, which is the run's setting and is not recorded in the table.** stageR's overall alpha is the run's setting; the `stageRTx` object was not found or could not be read, or the adjustment could not be re-run, so the values of the run are used as they stand for the FDR of this report, ",
         format(fdr), ", and the transcript-level values may not be adjusted p-values at that level.")
}
# the two group columns of the DEXSeq transcript table are named by the conditions of the contrast
read_dtu_groups <- function(path, treatment, control) {
  tab <- read_delim_chr(path, "\t")
  need_cols(tab, c("groupID", "featureID", treatment, control), path)
  data.frame(gene_id = strip_version(strip_quotes(tab$groupID)), feature_id = strip_tx_version(strip_quotes(tab$featureID)), treatment = to_num(tab[[treatment]], path, treatment),
             control = to_num(tab[[control]], path, control), stringsAsFactors = FALSE)
}
# the sample and the condition of each column of the DEXSeq transcript table (countData.1, countData.2, ...): the colData of the DEXSeqDataSet of the run
read_dtu_samples <- function(path) {
  if (!file.exists(path) || !requireNamespace("DEXSeq", quietly = TRUE)) return(NULL)
  tryCatch({
    sa <- as.data.frame(DEXSeq::sampleAnnotation(readRDS(path)))
    nm <- if ("sample.1" %in% names(sa)) sa[["sample.1"]] else sa[["sample"]]
    data.frame(sample = as.character(nm), condition = as.character(sa$condition), stringsAsFactors = FALSE)
  }, error = function(e) NULL)
}
# the raw counts of the transcripts: the columns countData.1 ... countData.N of the table, named by the samples; NULL when the columns do not match the samples
read_dtu_counts <- function(path, samples) {
  tab <- read_delim_chr(path, "\t")
  need_cols(tab, c("groupID", "featureID"), path)
  cc <- grep("^countData\\.[0-9]+$", names(tab), value = TRUE)
  cc <- cc[order(as.integer(sub("^countData\\.", "", cc)))]
  if (is.null(samples) || length(cc) != nrow(samples)) return(NULL)
  out <- data.frame(gene_id = strip_version(strip_quotes(tab$groupID)), feature_id = strip_tx_version(strip_quotes(tab$featureID)), stringsAsFactors = FALSE)
  for (i in seq_along(cc)) out[[paste0("count_", i)]] <- to_num(tab[[cc[i]]], path, cc[i])
  out
}
# the proportion of each transcript in its gene, for each sample: the count of the transcript divided by the sum of the counts of all the
# transcripts of the gene in that sample (NA when the gene has no count in the sample); the condition of the sample decides treatment or control
dtu_proportions <- function(cnt, samples, genes, treatment, control) {
  d <- cnt[cnt$gene_id %in% genes, , drop = FALSE]
  rows <- lapply(seq_len(nrow(samples)), function(i) {
    v <- d[[paste0("count_", i)]]
    tot <- stats::ave(v, d$gene_id, FUN = function(x) sum(x, na.rm = TRUE))
    role <- if (samples$condition[i] == treatment) "treatment" else if (samples$condition[i] == control) "control" else NA_character_
    data.frame(gene_id = d$gene_id, feature_id = d$feature_id, sample = samples$sample[i], condition = samples$condition[i], role = role,
               proportion = ifelse(tot > 0, v / tot, NA_real_), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out[!is.na(out$role), , drop = FALSE]
}

# volcano: x = log2 fold change (treatment minus control), y = -log10 adjusted p-value, the FDR as a dashed line, the top 10 significant labelled
volcano_usage <- function(tab, sig, title, fdr, treatment, control) {
  keep <- !is.na(tab$padj) & !is.na(tab$effect)
  d <- tab[keep, , drop = FALSE]
  d$sig <- sig[keep]
  d$nlp <- -log10(pmax(d$padj, 1e-300))
  d$status <- ifelse(d$sig, "significant", "not significant")
  lab <- d[d$sig, , drop = FALSE]
  lab <- head(lab[order(lab$padj, lab$gene_id), , drop = FALSE], 10)
  p <- ggplot(d, aes(x = effect, y = nlp, colour = status)) + geom_point(alpha = 0.7, size = 1.6) +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey40") +
    geom_hline(yintercept = -log10(fdr), linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = c("significant" = "#D55E00", "not significant" = "grey60")) +
    labs(title = title, x = paste0("log2 fold change (", treatment, " minus ", control, ")"), y = "-log10 adjusted p-value") +
    theme_bw()
  if (nrow(lab) > 0) {
    p <- if (requireNamespace("ggrepel", quietly = TRUE)) {
      p + ggrepel::geom_text_repel(data = lab, aes(label = gene_id), colour = "black", size = 3, max.overlaps = Inf)
    } else {
      p + geom_text(data = lab, aes(label = gene_id), colour = "black", size = 3, check_overlap = TRUE, vjust = -0.6)
    }
  }
  p
}
# exon usage along the exon number for some genes: one panel per gene, the bars are the log2 fold change of each bin
exon_usage_plot <- function(feat, genes, fdr, treatment, control, contrast) {
  d <- feat[feat$gene_id %in% genes & !is.na(feat$effect), , drop = FALSE]
  d$bin <- bin_number(d$feature_id)
  d$status <- ifelse(sig_usage(d, fdr), "significant", "not significant")
  d$gene_id <- factor(d$gene_id, levels = genes)
  ggplot(d, aes(x = bin, y = effect, fill = status)) + geom_col(width = 0.8) + geom_hline(yintercept = 0) +
    facet_wrap(~gene_id, scales = "free_x", ncol = 2) +
    scale_fill_manual(values = c("significant" = "#D55E00", "not significant" = "grey60")) +
    labs(title = paste0("Exon usage of the top DEXSeq genes, ", contrast), x = "exon bin number",
         y = paste0("log2 fold change in exon usage (", treatment, " minus ", control, ")")) +
    theme_bw()
}
# the proportion of each transcript in its gene for each sample (points) and the mean of the group (bar), treatment next to control; confirmed = keys "gene transcript" kept by stageR
dtu_proportion_plot <- function(prop, genes, confirmed, treatment, control, contrast) {
  d <- prop[prop$gene_id %in% genes & !is.na(prop$proportion), , drop = FALSE]
  d$label <- paste0(d$feature_id, ifelse(paste(d$gene_id, d$feature_id) %in% confirmed, " *", ""))
  d$role <- factor(d$role, levels = c("treatment", "control"), labels = c(paste0("treatment (", treatment, ")"), paste0("control (", control, ")")))
  d$gene_id <- factor(d$gene_id, levels = genes)
  ggplot(d, aes(x = label, y = proportion, colour = role, fill = role)) +
    stat_summary(fun = mean, geom = "col", position = position_dodge(width = 0.8), width = 0.7, alpha = 0.35) +
    geom_point(position = position_dodge(width = 0.8), size = 1.8) +
    facet_wrap(~gene_id, scales = "free_x", ncol = 2) +
    scale_colour_manual(values = c("#D55E00", "#0072B2"), name = "group") + scale_fill_manual(values = c("#D55E00", "#0072B2"), name = "group") +
    labs(title = paste0("Transcript proportions of the top DTU genes, ", contrast), x = "transcript (* confirmed by stageR)",
         y = "proportion of the gene's transcript counts") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
# when the counts and the samples are not available: the two group columns of the DEXSeq transcript table (usage coefficients, not proportions); confirmed = keys "gene transcript" kept by stageR
dtu_group_plot <- function(grp, genes, confirmed, treatment, control, contrast) {
  d <- grp[grp$gene_id %in% genes, , drop = FALSE]
  d$label <- paste0(d$feature_id, ifelse(paste(d$gene_id, d$feature_id) %in% confirmed, " *", ""))
  long <- rbind(data.frame(gene_id = d$gene_id, label = d$label, role = "treatment", value = d$treatment, stringsAsFactors = FALSE),
                data.frame(gene_id = d$gene_id, label = d$label, role = "control", value = d$control, stringsAsFactors = FALSE))
  long <- long[!is.na(long$value), , drop = FALSE]
  long$role <- factor(long$role, levels = c("treatment", "control"), labels = c(paste0("treatment (", treatment, ")"), paste0("control (", control, ")")))
  long$gene_id <- factor(long$gene_id, levels = genes)
  ggplot(long, aes(x = label, y = value, fill = role)) + geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    facet_wrap(~gene_id, scales = "free_x", ncol = 2) +
    scale_fill_manual(values = c("#D55E00", "#0072B2"), name = "group") +
    labs(title = paste0("Transcripts of the top DTU genes, ", contrast), x = "transcript (* confirmed by stageR)",
         y = "DEXSeq usage coefficient (fitted, variance-stabilised; not a proportion)") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
save_plot <- function(p, file, width = 7, height = 5) {
  ggplot2::ggsave(filename = file, plot = p, width = width, height = height, dpi = 150)
  invisible(file)
}
show_image <- function(file, alt) cat("![", alt, "](", file, ")\n\n", sep = "")
```

```{r report, results="asis"}
modules <- strsplit(trimws(params$modules), "[[:space:]]+")[[1]]
usage_modules <- c("dexseq_exon", "edger_exon", "dtu")
fig_dir <- file.path(params$report_dir, "figures"); tab_dir <- file.path(params$report_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE); dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
contrasts <- read_contrasts(params$contrasts)
path_of <- function(...) file.path(params$outdir, ...)
write_tab <- function(d, name) utils::write.csv(d, file.path(tab_dir, paste0(C, ".", name, ".csv")), row.names = FALSE)
show_tab <- function(d, n = 20) { if (nrow(d) > 0) cat(knitr::kable(head(d, n), format = "markdown"), sep = "\n"); cat("\n\n") }
load_rds <- function(tool, label) {
  f <- file.path(params$report_dir, "rds", paste0(C, ".", tool, ".rds"))
  if (!file.exists(f)) {
    cat("The file `", f, "` was not found: the import found no ", label, " table for this contrast (the files it looked for are listed in `files_missing` of `rds/",
        C, ".qc.rds`). No ", label, " result is shown for this contrast.\n\n", sep = "")
    return(NULL)
  }
  tab <- readRDS(f)
  need_cols(tab, c("gene_id", "feature_id", "effect", "padj", "direction"), f)
  tab
}
missing_file <- function(f, what) cat("The file `", f, "` was not found: ", what, "\n\n", sep = "")
# the empty-set message with the three nearest rows by adjusted p-value
none_block <- function(d, what) {
  cat("**", none_message(params$fdr), "** (", what, ")\n\n", sep = "")
  nr <- nearest_rows(d, 3)
  cat("The three nearest by adjusted p-value:\n\n")
  if (nrow(nr) > 0) show_tab(data.frame(id = nr$gene_id, padj = signif(nr$padj, 3), stringsAsFactors = FALSE)) else cat("Nothing was tested.\n\n")
}
if (!any(usage_modules %in% modules)) cat("None of the modules dexseq_exon, edger_exon and dtu is among the modules of this report (", params$modules, "): there is nothing to show.\n\n", sep = "")
for (i in seq_len(nrow(contrasts))) {
  if (!any(usage_modules %in% modules)) break
  C <- contrasts$contrast[i]; Tr <- contrasts$treatment[i]; Co <- contrasts$control[i]
  cat("\n\n## Contrast ", C, "\n\n", sep = "")
  cat("Treatment `", Tr, "`, control `", Co, "`: every effect below is `", Tr, "` minus `", Co, "`. A positive log2 fold change: ",
      usage_wording("treatment_up"), "; a negative log2 fold change: ", usage_wording("control_up"), ".\n\n", sep = "")
  cat("Threshold: FDR ", format(params$fdr), " (exon and transcript usage are cut by the adjusted p-value only). The exon and transcript usage tools use the design of the sample sheet as it is; this report does not change it.\n\n", sep = "")
  cnt <- list(); lost <- FALSE

  if ("dexseq_exon" %in% modules) {
    cat("\n\n### DEXSeq exon usage, ", C, "\n\n", sep = "")
    tab <- load_rds("exon_dexseq", "DEXSeq exon usage")
    if (is.null(tab)) lost <- TRUE else {
      sigb <- sig_usage(tab, params$fdr)
      cnt$dex_bin <- usage_counts("exon_dexseq", "exon_bin", tab$padj, sigb, tab$direction)
      pgf <- path_of("star", "dexseq_exon", "results", sprintf("perGeneQValue.%s.csv", C))
      pg <- NULL
      if (file.exists(pgf)) pg <- read_gene_table(pgf, ",", "groupID", "padj") else missing_file(pgf, "no gene-level result is shown for DEXSeq exon usage.")
      cat(sum(sigb), " exon bins with an adjusted p-value below ", format(params$fdr), " (of ", sum(!is.na(tab$padj)), " tested).\n\n", sep = "")
      if (!is.null(pg)) {
        cnt$dex_gene <- usage_counts("exon_dexseq", "gene", pg$padj, sig_usage(pg, params$fdr))
        gt <- gene_summary(pg, tab, sigb, params$fdr)
        cat(nrow(gt), " genes with a per-gene q-value (`perGeneQValue`) below ", format(params$fdr), " (of ", sum(!is.na(pg$padj)), " tested; table `tables/", C, ".exon_dexseq.genes.csv`).\n\n", sep = "")
        if (isTRUE(params$events)) {
          write_tab(gt, "exon_dexseq.genes")
          if (nrow(gt) > 0) show_tab(gt[, c("gene_id", "gene_padj", "n_sig_features", "top_feature", "top_effect", "direction_wording")]) else none_block(pg, "genes by per-gene q-value")
        }
      }
      if (isTRUE(params$events)) {
        bins <- tab[sigb, c("gene_id", "feature_id", "padj", "effect", "direction"), drop = FALSE]
        bins <- bins[order(bins$padj, bins$feature_id), , drop = FALSE]
        bins$direction_wording <- usage_wording(bins$direction)
        write_tab(bins, "exon_dexseq.bins")
        cat("Exon bins below the FDR (table `tables/", C, ".exon_dexseq.bins.csv`):\n\n", sep = "")
        if (nrow(bins) > 0) show_tab(transform(bins, padj = signif(padj, 3), effect = round(effect, 3))[, c("feature_id", "padj", "effect", "direction_wording")]) else
          none_block(data.frame(gene_id = tab$feature_id, padj = tab$padj, stringsAsFactors = FALSE), "exon bins")
        if (any(!is.na(tab$padj) & !is.na(tab$effect))) {
          vf <- file.path(fig_dir, sprintf("%s.usage.volcano.exon_dexseq.png", C))
          save_plot(volcano_usage(tab, sigb, paste0("DEXSeq exon usage, ", C), params$fdr, Tr, Co), vf)
          show_image(vf, paste0("Volcano plot of the DEXSeq exon bins, ", C))
        } else cat("No exon bin has both an effect and an adjusted p-value: no volcano plot.\n\n")
      } else cat("The tables and the volcano plot are switched off for this report (events = FALSE).\n\n")
      if (isTRUE(params$topgene)) {
        if (is.null(pg)) cat("No gene-level table: no exon usage plot.\n\n") else {
          cand <- pg[pg$gene_id %in% tab$gene_id[!is.na(tab$effect)] & !is.na(pg$padj), , drop = FALSE]
          top <- head(cand$gene_id[order(cand$padj, cand$gene_id)], 6)
          if (length(top) == 0) cat("No gene with a tested exon bin: no exon usage plot.\n\n") else {
            cat("Exon usage of the six genes with the smallest per-gene q-value (the significant genes first); the bars are the log2 fold change of each tested exon bin, coloured when the bin is below the FDR.\n\n")
            ef <- file.path(fig_dir, paste0(C, ".usage.exon_top.png"))
            save_plot(exon_usage_plot(tab, top, params$fdr, Tr, Co, C), ef, width = 9, height = 2 + 2 * ceiling(length(top) / 2))
            show_image(ef, paste0("Exon usage of the top DEXSeq genes, ", C))
          }
        }
      } else cat("The top-gene plots are switched off for this report (topgene = FALSE).\n\n")
    }
  }

  if ("edger_exon" %in% modules) {
    cat("\n\n### edgeR exon usage, ", C, "\n\n", sep = "")
    tab <- load_rds("exon_edger", "edgeR exon usage")
    if (is.null(tab)) lost <- TRUE else {
      sige <- sig_usage(tab, params$fdr)
      cnt$edger_exon <- usage_counts("exon_edger", "exon", tab$padj, sige, tab$direction)
      sf <- path_of("star", "edger", sprintf("contrast_%s.usage.simes.csv", C)); gf <- path_of("star", "edger", sprintf("contrast_%s.usage.gene.csv", C))
      simes <- NULL
      if (file.exists(sf)) simes <- read_gene_table(sf, ",", "Geneid", "FDR") else missing_file(sf, "no Simes gene-level result is shown for edgeR exon usage.")
      if (file.exists(gf)) {
        ft <- read_gene_table(gf, ",", "Geneid", "FDR")
        cnt$edger_ftest <- usage_counts("exon_edger", "gene_ftest", ft$padj, sig_usage(ft, params$fdr))
        cat(sum(sig_usage(ft, params$fdr)), " genes with an F-test FDR below ", format(params$fdr), " (gene-level F-test table, of ", sum(!is.na(ft$padj)), " tested).\n\n", sep = "")
      } else missing_file(gf, "the gene-level F-test table is not counted.")
      cat(sum(sige), " exons with an exon-level FDR below ", format(params$fdr), " (of ", sum(!is.na(tab$padj)), " tested).\n\n", sep = "")
      if (!is.null(simes)) {
        cnt$edger_simes <- usage_counts("exon_edger", "gene_simes", simes$padj, sig_usage(simes, params$fdr))
        gt <- gene_summary(simes, tab, sige, params$fdr)
        cat(nrow(gt), " genes with a Simes FDR below ", format(params$fdr), " (of ", sum(!is.na(simes$padj)), " tested; table `tables/", C, ".exon_edger.genes.csv`).\n\n", sep = "")
        if (isTRUE(params$events)) {
          write_tab(gt, "exon_edger.genes")
          if (nrow(gt) > 0) show_tab(gt[, c("gene_id", "gene_padj", "n_sig_features", "top_feature", "top_effect", "direction_wording")]) else none_block(simes, "genes by Simes FDR")
        }
      }
      if (isTRUE(params$events)) {
        ex <- tab[sige, c("gene_id", "feature_id", "padj", "effect", "direction"), drop = FALSE]
        ex <- ex[order(ex$padj, ex$feature_id), , drop = FALSE]
        ex$direction_wording <- usage_wording(ex$direction)
        write_tab(ex, "exon_edger.exons")
        cat("Exons below the FDR (table `tables/", C, ".exon_edger.exons.csv`):\n\n", sep = "")
        if (nrow(ex) > 0) show_tab(transform(ex, padj = signif(padj, 3), effect = round(effect, 3))[, c("feature_id", "padj", "effect", "direction_wording")]) else
          none_block(data.frame(gene_id = tab$feature_id, padj = tab$padj, stringsAsFactors = FALSE), "exons")
        if (any(!is.na(tab$padj) & !is.na(tab$effect))) {
          vf <- file.path(fig_dir, sprintf("%s.usage.volcano.exon_edger.png", C))
          save_plot(volcano_usage(tab, sige, paste0("edgeR exon usage, ", C), params$fdr, Tr, Co), vf)
          show_image(vf, paste0("Volcano plot of the edgeR exons, ", C))
        } else cat("No exon has both an effect and an adjusted p-value: no volcano plot.\n\n")
      } else cat("The tables and the volcano plot are switched off for this report (events = FALSE).\n\n")
    }
  }

  if ("dtu" %in% modules) {
    cat("\n\n### Transcript usage (DEXSeq DTU with stageR screening), ", C, "\n\n", sep = "")
    tab <- load_rds("dtu", "DEXSeq transcript usage")
    if (is.null(tab)) lost <- TRUE else {
      sf <- path_of("salmon", "dexseq_dtu", "results", "stager", sprintf("getAdjustedPValues.%s-%s.tsv", Tr, Co))
      if (!file.exists(sf)) {
        missing_file(sf, "no stageR result is shown. The DEXSeq transcript adjusted p-values alone are counted below.")
        cnt$dtu_dexseq <- usage_counts("dtu", "transcript_dexseq", tab$padj, sig_usage(tab, params$fdr), tab$direction)
      } else {
        st <- read_stager(sf)
        sres <- stager_for_fdr(st, path_of("salmon", "dexseq_dtu", "results", "stager", sprintf("stageRTx.%s-%s.rds", Tr, Co)), params$fdr)
        st <- sres$st
        st$gene_id <- strip_version(st$gene_id); st$tx_id <- strip_tx_version(st$tx_id)
        tab$gene_id <- strip_version(tab$gene_id); tab$feature_id <- strip_tx_version(tab$feature_id)
        cat(stager_note(sres, params$fdr), "\n\n")
        cat("The column `gene` of the stageR table is the screened adjusted p-value at the gene level and the column `transcript` the screened adjusted p-value at the transcript level (`NA` where the gene was not screened in).\n\n")
        sg <- stager_genes(st)
        txsig <- !is.na(st$transcript) & st$transcript < params$fdr
        confirmed <- paste(st$gene_id, st$tx_id)[txsig]
        sigt <- paste(tab$gene_id, tab$feature_id) %in% confirmed
        k <- match(paste(st$gene_id, st$tx_id), paste(tab$gene_id, tab$feature_id))
        cnt$dtu_gene <- usage_counts("dtu", "gene_stager", sg$padj, sig_usage(sg, params$fdr))
        cnt$dtu_tx <- usage_counts("dtu", "transcript_stager", st$transcript, txsig, tab$direction[k])
        cnt$dtu_dexseq <- usage_counts("dtu", "transcript_dexseq", tab$padj, sig_usage(tab, params$fdr), tab$direction)
        gt <- gene_summary(sg, tab, sigt, params$fdr)
        cat(nrow(gt), " genes with the stageR gene-level adjusted p-value (`gene`) below ", format(params$fdr), " (of ", sum(!is.na(sg$padj)), " screened; table `tables/", C, ".dtu.genes.csv`). ",
            sum(txsig), " transcripts with the stageR transcript-level adjusted p-value (`transcript`) below ", format(params$fdr), " (table `tables/", C, ".dtu.transcripts.csv`). ",
            "For comparison, ", sum(sig_usage(tab, params$fdr)), " transcripts have a DEXSeq adjusted p-value below ", format(params$fdr), " before the screening.\n\n", sep = "")
        if (isTRUE(params$events)) {
          tx <- data.frame(gene_id = st$gene_id[txsig], transcript_id = st$tx_id[txsig], stager_transcript_padj = st$transcript[txsig], dexseq_padj = tab$padj[k][txsig],
                           effect = tab$effect[k][txsig], direction = tab$direction[k][txsig], stringsAsFactors = FALSE)
          tx <- tx[order(tx$stager_transcript_padj, tx$dexseq_padj, tx$transcript_id), , drop = FALSE]
          tx$direction_wording <- usage_wording(tx$direction)
          write_tab(gt, "dtu.genes"); write_tab(tx, "dtu.transcripts")
          if (nrow(gt) > 0) show_tab(gt[, c("gene_id", "gene_padj", "n_sig_features", "top_feature", "top_effect", "direction_wording")]) else none_block(sg, "genes by the stageR gene-level adjusted p-value")
          if (nrow(tx) > 0) show_tab(transform(tx, effect = round(effect, 3), stager_transcript_padj = signif(stager_transcript_padj, 3), dexseq_padj = signif(dexseq_padj, 3))[, c("gene_id", "transcript_id", "stager_transcript_padj", "effect", "direction_wording")]) else
            none_block(data.frame(gene_id = tab$feature_id, padj = tab$padj, stringsAsFactors = FALSE), "transcripts by the stageR transcript-level adjusted p-value")
          if (any(!is.na(tab$padj) & !is.na(tab$effect))) {
            vf <- file.path(fig_dir, sprintf("%s.usage.volcano.dtu.png", C))
            save_plot(volcano_usage(tab, sigt, paste0("DEXSeq transcript usage, ", C), params$fdr, Tr, Co), vf)
            cat("Volcano plot of the transcripts: the y axis is the DEXSeq adjusted p-value of the transcript; the transcripts confirmed by stageR (`transcript` below the FDR) are coloured.\n\n")
            show_image(vf, paste0("Volcano plot of the DTU transcripts, ", C))
          } else cat("No transcript has both an effect and an adjusted p-value: no volcano plot.\n\n")
        } else cat("The tables and the volcano plot are switched off for this report (events = FALSE).\n\n")
        if (isTRUE(params$topgene)) {
          dxf <- path_of("salmon", "dexseq_dtu", "results", "dexseq", sprintf("DEXSeqResults.%s-%s.tsv", Tr, Co))
          cand <- sg[!is.na(sg$padj), , drop = FALSE]
          top <- head(cand$gene_id[order(cand$padj, cand$gene_id)], 6)
          if (!file.exists(dxf)) missing_file(dxf, "no transcript plot of the top genes.") else if (length(top) == 0) cat("No gene was screened: no transcript plot.\n\n") else {
            grp <- read_dtu_groups(dxf, Tr, Co)
            pf <- file.path(fig_dir, paste0(C, ".usage.dtu_top.png"))
            samples <- read_dtu_samples(path_of("salmon", "dexseq_dtu", "results", "dexseq", sprintf("DEXSeqDataSet.%s-%s.rds", Tr, Co)))
            cnts <- read_dtu_counts(dxf, samples)
            if (!is.null(cnts) && all(samples$condition %in% c(Tr, Co))) {
              prop <- dtu_proportions(cnts, samples, top, Tr, Co)
              write_tab(prop, "usage.dtu_proportions")
              cat("The six genes with the smallest stageR gene-level adjusted p-value. For each sample the proportion of a transcript is its count (columns `countData.N` of the DEXSeq table) divided by the sum of the counts of all the transcripts of the gene in that sample; the group of each sample is read from the sample annotation of the `DEXSeqDataSet` of the run (`", Tr,
                  "` treatment, `", Co, "` control). Points are samples, bars the mean of the group; transcripts confirmed by stageR are marked with * (table `tables/", C, ".usage.dtu_proportions.csv`).\n\n", sep = "")
              save_plot(dtu_proportion_plot(prop, top, confirmed, Tr, Co, C), pf, width = 9, height = 2 + 2 * ceiling(length(top) / 2))
            } else {
              cat("The `DEXSeqDataSet` of the run (`DEXSeqDataSet.", Tr, "-", Co, ".rds`) was not found or does not match the counts of the DEXSeq table, so the transcript proportions of the samples cannot be drawn. ",
                  "The plot below shows the two group columns of the DEXSeq table, `", Tr, "` (treatment) and `", Co, "` (control): DEXSeq usage coefficient (fitted, variance-stabilised; not a proportion). ",
                  "Transcripts confirmed by stageR are marked with *.\n\n", sep = "")
              save_plot(dtu_group_plot(grp, top, confirmed, Tr, Co, C), pf, width = 9, height = 2 + 2 * ceiling(length(top) / 2))
            }
            show_image(pf, paste0("Transcripts of the top DTU genes, ", C))
          }
        } else cat("The top-gene plots are switched off for this report (topgene = FALSE).\n\n")
      }
    }
  }

  cat("\n\n### Counts, ", C, "\n\n", sep = "")
  cn <- if (length(cnt)) do.call(rbind, unname(cnt)) else usage_counts("none", "none", numeric(0), logical(0))[0, ]
  rownames(cn) <- NULL
  write_tab(cn, "usage.counts")
  cat("Tested and significant counts of every level (counts file `tables/", C, ".usage.counts.csv`); the direction columns are empty for the gene-level rows:\n\n", sep = "")
  cat(knitr::kable(cn, format = "markdown"), sep = "\n")
  cat("\n\n")
  if (lost) {
    qf <- file.path(params$report_dir, "rds", paste0(C, ".qc.rds"))
    if (file.exists(qf)) { qc <- readRDS(qf); if (length(qc$files_missing)) cat("Files the import did not find:\n\n", paste0("- `", qc$files_missing, "`\n"), "\n", sep = "") }
  }
}
```
````

The cross-tool report `05_cross_tool.Rmd` is generated only when two or more tools ran (Step 10). It compares, for every contrast, the genes called by each tool that has a result. It reads `rds/{CONTRAST}.{tool}.rds` (written by `01_import.Rmd`) and, for the gene-level calls of DEXSeq exon usage, edgeR exon usage and DTU, the same gene-level files of the run as `04_usage.Rmd` (`perGeneQValue.{CONTRAST}.csv`, `usage.simes.csv`, stageR's `getAdjustedPValues.{T}-{C}.tsv`), read only. A gene is called by rMATS or SUPPA2 when at least one of its events is significant (adjusted p-value below `{FDR}`, absolute dPSI at least `{DPSI}`, and for rMATS not flagged `low_support`); by DEXSeq exon usage when its `perGeneQValue` is below `{FDR}`; by edgeR exon usage when its Simes FDR is below `{FDR}`; by DTU when stageR's `gene` column is below `{FDR}` (stageR's adjusted p-values are valid only at the alpha they were computed for: at another FDR the stage-wise adjustment is re-run from the `stageRTx` object of the run, or the values of the run are used with a warning). Gene identifiers are Ensembl identifiers; a version suffix is removed only where one is present, and a group of overlapping genes (a DEXSeq group `GENE1+GENE2` or a SUPPA2 locus `GENE1_and_GENE2`) counts for each of its genes. The direction of a gene in a tool is the sign of the effect of the feature with the smallest adjusted p-value among that tool's significant features in the gene (rMATS and SUPPA2 on dPSI, the other tools on log2 fold change): `treatment_up`, `control_up`, or empty when the gene is called at gene level and none of its features is significant. The concordance table counts, for each pair of tools, the genes called by both with the same direction and with the opposite direction; directions compared only, never effect sizes across scales. The background is stated in the report: a gene outside a tool's set was either not significant or not tested by that tool. It writes, per contrast, `tables/{CONTRAST}.cross_tool.sets.csv`, `.cross_tool.membership.csv`, `.cross_tool.intersections.csv`, `.cross_tool.genes.csv` (the genes called by two or more tools, with the direction and the harmonised effect of each tool: `direction_{tool}`, `effect_{tool}`), `.cross_tool.concordance.csv` and `figures/{CONTRAST}.upset.png` (ComplexUpset when the image has it, otherwise UpSetR). A tool whose RDS or gene-level file is missing is left out and named; fewer than two tools left means no comparison, and the report still renders.

**File `05_cross_tool.Rmd`.**
````rmd
---
title: "Cross-tool comparison"
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
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE, dev = "png", dpi = 150, fig.width = 8, fig.height = 5)
```

This report compares the genes called by each tool that ran, for each contrast, from the tables that `01_import.Rmd` wrote to `r file.path(params$report_dir, "rds")` and from the gene-level files of the run (read only). Every effect is stated as **treatment minus control**: a positive effect is higher in the treatment group (more inclusion for rMATS and SUPPA2, higher usage for the exon and transcript usage tools). The direction of a gene is `treatment_up` for a positive effect and `control_up` for a negative one, as in the harmonised column `direction` of the RDS, which is used as it is.

A gene is called by a tool when: rMATS and SUPPA2, at least one event of the gene has an adjusted p-value below `r params$fdr`, an absolute dPSI of at least `r params$dpsi`, and (rMATS) is not flagged `low_support`; DEXSeq exon usage, the per-gene q-value (`perGeneQValue`) is below the FDR; edgeR exon usage, the Simes FDR is below the FDR; DTU, the stageR column `gene` (the screened adjusted p-value at the gene level) is below the FDR (stageR's adjusted p-values are valid only at the alpha they were computed for: at another FDR the stage-wise adjustment is re-run from the `stageRTx` object of the run, or the values of the run are used with a warning). Gene identifiers are Ensembl identifiers; a version suffix is removed only where one is present, and a group of overlapping genes (a DEXSeq group `GENE1+GENE2` or a SUPPA2 locus `GENE1_and_GENE2`) counts for each of its genes.

The direction of a gene in a tool is the sign of the effect of the feature with the smallest adjusted p-value among that tool's significant features in the gene (dPSI for rMATS and SUPPA2, log2 fold change for the other tools); it is empty when the gene is called at gene level but none of its features is significant. The cross-tool table has, for each tool, `direction_<tool>` and `effect_<tool>`, the effect of that feature (dPSI for rMATS and SUPPA2, log2 fold change for the others; the two scales are not comparable). The concordance table counts, for each pair of tools, the genes called by both with the same direction and with the opposite direction (directions compared only, never effect sizes across scales: dPSI and log2 fold change are not on one scale); genes without a direction in one of the two tools are counted separately.

Background: the tools do not test the same genes. A gene outside the set of a tool was either not significant or not tested by that tool, so the overlaps describe the genes called and say nothing about genes a tool could not test; the number of genes each tool tested is in the sets table.

```{r functions}
suppressPackageStartupMessages(library(ggplot2))
TOOL_LABELS <- c(rmats = "rMATS", suppa = "SUPPA2", exon_dexseq = "DEXSeq exon", exon_edger = "edgeR exon", dtu = "DEXSeq DTU")
MODULE_OF_TOOL <- c(rmats = "rmats", suppa = "suppa", exon_dexseq = "dexseq_exon", exon_edger = "edger_exon", dtu = "dtu")

# an Ensembl gene identifier loses its version suffix (.5) only when it has one
strip_version <- function(x) sub("^(ENS[A-Z]*G[0-9]+)\\.[0-9]+$", "\\1", x)
# a transcript identifier loses a trailing version in the same way (the pipeline strips the versions of the stageR identifiers, not those of the DEXSeq table: both sides of a join are stripped)
strip_tx_version <- function(x) sub("\\.[0-9]+$", "", x)
# a DEXSeq group of overlapping genes (GENE1+GENE2) or a SUPPA2 locus of overlapping genes (GENE1_and_GENE2) becomes one row per gene
expand_genes <- function(df) {
  parts <- strsplit(df$gene_id, "\\+|_and_")
  out <- df[rep(seq_len(nrow(df)), lengths(parts)), , drop = FALSE]
  out$gene_id <- as.character(unlist(parts))
  rownames(out) <- NULL
  out
}
tidy_ids <- function(df) { df$gene_id <- strip_version(df$gene_id); df }

# significant: adjusted p-value below the FDR; for events also the absolute dPSI threshold and not low_support
sig_usage <- function(tab, fdr) !is.na(tab$padj) & tab$padj < fdr
sig_events <- function(tab, fdr, dpsi) {
  !is.na(tab$padj) & tab$padj < fdr & !is.na(tab$effect) & abs(tab$effect) >= dpsi & !tab$low_support
}

# per gene: the significant features, the one with the smallest adjusted p-value (then the larger change, then the identifier) and its direction
top_by_gene <- function(feat, sig) {
  s <- feat[sig & !is.na(feat$padj), c("gene_id", "feature_id", "padj", "effect", "direction"), drop = FALSE]
  s <- s[order(s$padj, 1 - abs(s$effect), s$feature_id), , drop = FALSE]
  n <- table(s$gene_id)
  top <- s[!duplicated(s$gene_id), , drop = FALSE]
  data.frame(gene_id = top$gene_id, n_sig_features = as.integer(n[top$gene_id]), top_feature = top$feature_id, top_effect = top$effect,
             top_padj = top$padj, direction = ifelse(top$direction %in% c("treatment_up", "control_up"), top$direction, NA_character_),
             stringsAsFactors = FALSE)
}
# the calls of one tool: the genes called, the genes tested, and the direction of every gene that has a significant feature
make_calls <- function(genes, tested, feat, sig_feat) {
  list(genes = sort(unique(as.character(genes))), tested = sort(unique(as.character(tested))), dir = top_by_gene(feat, sig_feat))
}
# rMATS and SUPPA2: a gene is called when one of its events is significant
calls_from_events <- function(tab, fdr, dpsi) {
  tab <- tidy_ids(expand_genes(tab))
  sig <- sig_events(tab, fdr, dpsi)
  make_calls(tab$gene_id[sig], tab$gene_id[!is.na(tab$padj)], tab, sig)
}
# DEXSeq and edgeR exon usage: a gene is called by the gene-level table (gene_id, padj); the features are the exon bins or exons of the RDS
calls_from_gene_table <- function(gene_tab, feat, fdr) {
  g <- tidy_ids(expand_genes(gene_tab))
  f <- tidy_ids(expand_genes(feat))
  make_calls(g$gene_id[sig_usage(g, fdr)], g$gene_id[!is.na(g$padj)], f, sig_usage(f, fdr))
}
# DTU: a gene is called when stageR's gene is below the FDR; its significant features are the transcripts with stageR's transcript below the FDR
calls_from_stager <- function(st, feat, fdr) {
  s <- tidy_ids(st)
  f <- tidy_ids(feat)
  s$tx_id <- strip_tx_version(s$tx_id)
  f$feature_id <- strip_tx_version(f$feature_id)
  txsig <- !is.na(s$transcript) & s$transcript < fdr
  confirmed <- paste(s$gene_id, s$tx_id)[txsig]
  sigf <- paste(f$gene_id, f$feature_id) %in% confirmed
  make_calls(s$gene_id[!is.na(s$gene) & s$gene < fdr], s$gene_id[!is.na(s$gene)], f, sigf)
}

# one row per gene of the union, one logical column per tool
membership_table <- function(calls) {
  genes <- sort(unique(as.character(unlist(lapply(calls, function(x) x$genes)))))
  out <- data.frame(gene_id = genes, stringsAsFactors = FALSE)
  for (t in names(calls)) out[[t]] <- out$gene_id %in% calls[[t]]$genes
  out
}
# the UpSet intersections: the number of genes for each exact combination of tools
intersection_counts <- function(calls) {
  m <- membership_table(calls)
  tools <- names(calls)
  if (nrow(m) == 0) return(data.frame(combination = character(0), n_genes = integer(0), stringsAsFactors = FALSE))
  mat <- as.matrix(m[tools])
  combo <- vapply(seq_len(nrow(mat)), function(i) paste(tools[mat[i, ]], collapse = "&"), character(1))
  tb <- table(combo)
  out <- data.frame(combination = names(tb), n_genes = as.integer(tb), stringsAsFactors = FALSE)
  out <- out[order(out$n_genes, out$combination, decreasing = c(TRUE, FALSE), method = "radix"), , drop = FALSE]
  rownames(out) <- NULL
  out
}
# the genes called by min_tools or more tools, with the direction of each tool (empty when the tool did not call the gene or has no direction)
cross_table <- function(calls, min_tools = 2) {
  m <- membership_table(calls)
  tools <- names(calls)
  mat <- as.matrix(m[tools])
  n <- rowSums(mat)
  keep <- n >= min_tools
  out <- data.frame(gene_id = m$gene_id[keep], n_tools = as.integer(n[keep]),
                    tools = vapply(which(keep), function(i) paste(tools[mat[i, ]], collapse = "; "), character(1)), stringsAsFactors = FALSE)
  for (t in tools) {
    d <- calls[[t]]$dir
    out[[paste0("direction_", t)]] <- ifelse(out$gene_id %in% calls[[t]]$genes, d$direction[match(out$gene_id, d$gene_id)], NA_character_)
    out[[paste0("effect_", t)]] <- ifelse(out$gene_id %in% calls[[t]]$genes, d$top_effect[match(out$gene_id, d$gene_id)], NA_real_)
  }
  out <- out[order(out$n_tools, out$gene_id, decreasing = c(TRUE, FALSE), method = "radix"), , drop = FALSE]
  rownames(out) <- NULL
  out
}
# for each pair of tools: the genes called by both, with the same direction, with the opposite direction, and without a direction in one of them
# (directions are compared, never the sizes of the effects: dPSI and log2 fold change are not on one scale)
concordance_table <- function(calls) {
  tools <- names(calls)
  rows <- list()
  if (length(tools) >= 2) for (i in seq_len(length(tools) - 1)) for (j in (i + 1):length(tools)) {
    a <- tools[i]; b <- tools[j]
    both <- intersect(calls[[a]]$genes, calls[[b]]$genes)
    da <- calls[[a]]$dir$direction[match(both, calls[[a]]$dir$gene_id)]
    db <- calls[[b]]$dir$direction[match(both, calls[[b]]$dir$gene_id)]
    known <- !is.na(da) & !is.na(db)
    same <- sum(known & da == db)
    opposite <- sum(known & da != db)
    rows[[length(rows) + 1]] <- data.frame(tool_a = a, tool_b = b, genes_in_both = length(both), same_direction = same, opposite_direction = opposite,
                                           no_direction = length(both) - same - opposite, stringsAsFactors = FALSE)
  }
  if (length(rows)) do.call(rbind, rows) else
    data.frame(tool_a = character(0), tool_b = character(0), genes_in_both = integer(0), same_direction = integer(0), opposite_direction = integer(0),
               no_direction = integer(0), stringsAsFactors = FALSE)
}

# the UpSet plot of the sets that are not empty: ComplexUpset when the image has it, otherwise UpSetR; FALSE when fewer than two sets have a gene
upset_engine <- function() {
  if (requireNamespace("ComplexUpset", quietly = TRUE)) "ComplexUpset" else if (requireNamespace("UpSetR", quietly = TRUE)) "UpSetR" else "none"
}
make_upset <- function(calls, file, engine = upset_engine()) {
  sets <- lapply(calls, function(x) x$genes)
  names(sets) <- unname(TOOL_LABELS[names(calls)])
  sets <- sets[vapply(sets, length, integer(1)) > 0]
  if (length(sets) < 2 || engine == "none") return(FALSE)
  if (engine == "ComplexUpset") {
    genes <- sort(unique(unlist(sets)))
    m <- data.frame(gene_id = genes, stringsAsFactors = FALSE)
    for (s in names(sets)) m[[s]] <- genes %in% sets[[s]]
    ggplot2::ggsave(filename = file, plot = ComplexUpset::upset(m, names(sets), name = "genes", width_ratio = 0.2), width = 8, height = 5, dpi = 150)
  } else {
    grDevices::png(file, width = 8, height = 5, units = "in", res = 150)
    print(UpSetR::upset(UpSetR::fromList(sets), nsets = length(sets), order.by = "freq"))
    grDevices::dev.off()
  }
  TRUE
}

strip_quotes <- function(x) gsub("^\"|\"$", "", x)
read_delim_chr <- function(path, sep) {
  if (!file.exists(path)) stop("file ", path, " not found: none of its columns can be read")
  l <- sub("\r$", "", readLines(path, warn = FALSE))
  l <- l[nzchar(l)]
  if (length(l) == 0) stop("file ", path, " is empty: no header line, so no column can be found")
  nm <- strsplit(l[1], sep, fixed = TRUE)[[1]]
  if (length(l) == 1) {
    return(as.data.frame(stats::setNames(replicate(length(nm), character(0), simplify = FALSE), nm),
                         check.names = FALSE, stringsAsFactors = FALSE))
  }
  tab <- utils::read.table(text = l[-1], sep = sep, quote = "", header = FALSE, colClasses = "character", na.strings = character(0),
                           comment.char = "", col.names = paste0("V", seq_along(nm)), fill = FALSE, stringsAsFactors = FALSE)
  names(tab) <- nm
  tab
}
need_cols <- function(tab, cols, path) {
  miss <- setdiff(cols, names(tab))
  if (length(miss)) stop("file ", path, ": column not found: ", paste(miss, collapse = ", "),
                         " (columns present: ", paste(names(tab), collapse = ", "), ")")
}
to_num <- function(x, path, col) {
  x <- trimws(strip_quotes(x))
  miss <- x %in% c("", "NA", "nan", "NaN", "-nan")
  v <- suppressWarnings(as.numeric(ifelse(miss, NA, x)))
  bad <- is.na(v) & !miss
  if (any(bad)) stop("file ", path, ": column ", col, " has a value that is not a number: '", x[bad][1], "'")
  v
}
read_contrasts <- function(path) {
  tab <- read_delim_chr(path, ",")
  need_cols(tab, c("contrast", "treatment", "control"), path)
  tab <- tab[nzchar(tab$contrast), c("contrast", "treatment", "control")]
  rownames(tab) <- NULL
  tab
}
# a gene-level table: the gene identifier and the adjusted p-value (DEXSeq perGeneQValue: groupID, padj; edgeR Simes: Geneid, FDR)
read_gene_table <- function(path, sep, id_col, p_col) {
  tab <- read_delim_chr(path, sep)
  need_cols(tab, c(id_col, p_col), path)
  data.frame(gene_id = strip_quotes(tab[[id_col]]), padj = to_num(tab[[p_col]], path, p_col), stringsAsFactors = FALSE)
}
# stageR: gene and transcript are the screened adjusted p-values at the two levels, NA where not screened
read_stager <- function(path) {
  tab <- read_delim_chr(path, "\t")
  need_cols(tab, c("geneID", "txID", "gene", "transcript"), path)
  data.frame(gene_id = strip_quotes(tab$geneID), tx_id = strip_quotes(tab$txID), gene = to_num(tab$gene, path, "gene"),
             transcript = to_num(tab$transcript, path, "transcript"), stringsAsFactors = FALSE)
}
# the overall alpha stageR was run with, read from the stageRTx object with getAlpha; NA when the file is absent or cannot be read
stager_alpha <- function(path) {
  if (!file.exists(path) || !requireNamespace("stageR", quietly = TRUE)) return(NA_real_)
  tryCatch(as.numeric(stageR::getAlpha(readRDS(path)))[1], error = function(e) NA_real_)
}
# stageR's adjusted p-values are valid only at the alpha they were computed for. At another FDR the stage-wise adjustment is re-run
# from the stageRTx object of the run (stageR is in the image); when the object cannot be read the values of the run are used as they stand
stager_for_fdr <- function(st, path, fdr) {
  alpha <- stager_alpha(path)
  if (is.na(alpha)) return(list(st = st, alpha = NA_real_, mode = "unreadable"))
  if (isTRUE(all.equal(alpha, fdr))) return(list(st = st, alpha = alpha, mode = "run"))
  re <- tryCatch({
    obj <- stageR::stageWiseAdjustment(readRDS(path), method = "dtu", alpha = fdr, allowNA = TRUE)
    r <- stageR::getAdjustedPValues(obj, onlySignificantGenes = FALSE, order = FALSE)
    data.frame(gene_id = as.character(r$geneID), tx_id = as.character(r$txID), gene = as.numeric(r$gene), transcript = as.numeric(r$transcript),
               stringsAsFactors = FALSE)
  }, error = function(e) NULL)
  if (is.null(re)) return(list(st = st, alpha = alpha, mode = "unreadable"))
  list(st = re, alpha = alpha, mode = "rerun")
}
# the sentence that tells which stageR values are used
stager_note <- function(res, fdr) {
  if (res$mode == "run") return(paste0("stageR was run with an overall alpha of ", format(res$alpha), " (read from the `stageRTx` object of the run), equal to the FDR of this report: the stageR table of the run is used as it is."))
  if (res$mode == "rerun") return(paste0("stageR's adjusted p-values are valid only at the alpha they were computed for. The run used an overall alpha of ", format(res$alpha),
    " (read from the `stageRTx` object of the run) and the FDR of this report is ", format(fdr), ", so the stage-wise adjustment was re-run at ", format(fdr),
    " from that object with the stageR package of the image; the values of the run are not used."))
  paste0("**WARNING: the stageR values of the run are valid only at the alpha the run used, which is the run's setting and is not recorded in the table.** stageR's overall alpha is the run's setting; the `stageRTx` object was not found or could not be read, or the adjustment could not be re-run, so the values of the run are used as they stand for the FDR of this report, ",
         format(fdr), ", and the transcript-level values may not be adjusted p-values at that level.")
}
show_image <- function(file, alt) cat("![", alt, "](", file, ")\n\n", sep = "")
```

```{r report, results="asis"}
modules <- strsplit(trimws(params$modules), "[[:space:]]+")[[1]]
fig_dir <- file.path(params$report_dir, "figures"); tab_dir <- file.path(params$report_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE); dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
contrasts <- read_contrasts(params$contrasts)
path_of <- function(...) file.path(params$outdir, ...)
write_tab <- function(d, name) utils::write.csv(d, file.path(tab_dir, paste0(C, ".cross_tool.", name, ".csv")), row.names = FALSE)
if (!isTRUE(params$cross)) cat("The cross-tool comparison is switched off for this report (cross = FALSE).\n\n")
for (i in seq_len(nrow(contrasts))) {
  if (!isTRUE(params$cross)) break
  C <- contrasts$contrast[i]; Tr <- contrasts$treatment[i]; Co <- contrasts$control[i]
  cat("\n\n## Contrast ", C, "\n\n", sep = "")
  cat("Treatment `", Tr, "`, control `", Co, "`: every effect and direction below is `", Tr, "` minus `", Co, "`; `treatment_up` is higher in `", Tr,
      "`, `control_up` is higher in `", Co, "`.\n\n", sep = "")
  cat("Thresholds: FDR ", format(params$fdr), ", abs(dPSI) ", format(params$dpsi), " (rMATS and SUPPA2 events only); rMATS minimum read support ", params$min_reads, ".\n\n", sep = "")
  calls <- list(); notfound <- character(0); stnote <- NULL
  for (tool in names(TOOL_LABELS)) {
    if (!MODULE_OF_TOOL[[tool]] %in% modules) next
    f <- file.path(params$report_dir, "rds", paste0(C, ".", tool, ".rds"))
    if (!file.exists(f)) { notfound <- c(notfound, f); next }
    tab <- readRDS(f)
    need_cols(tab, c("gene_id", "feature_id", "effect", "padj", "direction", "low_support"), f)
    cl <- NULL
    if (tool %in% c("rmats", "suppa")) {
      cl <- calls_from_events(tab, params$fdr, params$dpsi)
    } else if (tool == "exon_dexseq") {
      pgf <- path_of("star", "dexseq_exon", "results", sprintf("perGeneQValue.%s.csv", C))
      if (file.exists(pgf)) { pg <- read_gene_table(pgf, ",", "groupID", "padj"); cl <- calls_from_gene_table(pg, tab, params$fdr) } else notfound <- c(notfound, pgf)
    } else if (tool == "exon_edger") {
      sf <- path_of("star", "edger", sprintf("contrast_%s.usage.simes.csv", C))
      if (file.exists(sf)) { sm <- read_gene_table(sf, ",", "Geneid", "FDR"); cl <- calls_from_gene_table(sm, tab, params$fdr) } else notfound <- c(notfound, sf)
    } else if (tool == "dtu") {
      gf <- path_of("salmon", "dexseq_dtu", "results", "stager", sprintf("getAdjustedPValues.%s-%s.tsv", Tr, Co))
      if (file.exists(gf)) {
        st <- read_stager(gf)
        sres <- stager_for_fdr(st, path_of("salmon", "dexseq_dtu", "results", "stager", sprintf("stageRTx.%s-%s.rds", Tr, Co)), params$fdr)
        stnote <- stager_note(sres, params$fdr)
        cl <- calls_from_stager(sres$st, tab, params$fdr)
      } else notfound <- c(notfound, gf)
    }
    if (!is.null(cl)) calls[[tool]] <- cl
  }
  if (length(notfound)) cat("Files not found (the tool is left out of the comparison):\n\n", paste0("- `", notfound, "`\n"), "\n", sep = "")
  if (!is.null(stnote)) cat("DTU, stageR: ", stnote, "\n\n", sep = "")
  if (length(calls) < 2) {
    cat("There are fewer than two tools with a result for this contrast", if (length(calls) == 1) paste0(" (only ", TOOL_LABELS[[names(calls)]], ")") else "",
        ": no comparison is made.\n\n", sep = "")
    next
  }
  sets <- data.frame(tool = names(calls), label = unname(TOOL_LABELS[names(calls)]), tested_genes = vapply(calls, function(x) length(x$tested), integer(1)),
                     significant_genes = vapply(calls, function(x) length(x$genes), integer(1)), stringsAsFactors = FALSE)
  rownames(sets) <- NULL
  write_tab(sets, "sets")
  write_tab(membership_table(calls), "membership")
  ic <- intersection_counts(calls)
  write_tab(ic, "intersections")
  cat("### Genes called by each tool, ", C, "\n\n", sep = "")
  cat("Tested genes and significant genes per tool (sets file `tables/", C, ".cross_tool.sets.csv`):\n\n", sep = "")
  cat(knitr::kable(sets, format = "markdown"), sep = "\n")
  cat("\n\n")
  uf <- file.path(fig_dir, paste0(C, ".upset.png"))
  engine <- upset_engine()
  cat("### UpSet plot, ", C, "\n\n", sep = "")
  if (isTRUE(make_upset(calls, uf, engine))) {
    cat("Intersections of the gene sets of the tools with at least one gene, drawn with ", engine, " (figure `figures/", C, ".upset.png`).\n\n", sep = "")
    show_image(uf, paste0("UpSet plot of the significant genes of each tool, ", C))
  } else cat("Fewer than two tools have a significant gene (or neither ComplexUpset nor UpSetR is in the image): no UpSet plot.\n\n")
  cat("Exact intersections (genes in exactly these tools; table `tables/", C, ".cross_tool.intersections.csv`):\n\n", sep = "")
  if (nrow(ic) > 0) cat(knitr::kable(ic, format = "markdown"), sep = "\n") else cat("No gene is called by any tool.\n")
  cat("\n\n### Genes called by two or more tools, ", C, "\n\n", sep = "")
  xt <- cross_table(calls)
  write_tab(xt, "genes")
  if (nrow(xt) == 0) {
    cat("**none**: no gene is called by two or more tools at FDR ", format(params$fdr), ", |dPSI| ", format(params$dpsi), ".\n\n", sep = "")
  } else {
    cat(nrow(xt), " genes are called by two or more tools, with the direction of each tool (table `tables/", C, ".cross_tool.genes.csv`):\n\n", sep = "")
    cat(knitr::kable(head(xt, 30), format = "markdown"), sep = "\n")
    cat("\n\n")
  }
  cat("### Direction concordance, ", C, "\n\n", sep = "")
  cc <- concordance_table(calls)
  write_tab(cc, "concordance")
  cat("For each pair of tools, the genes called by both with the same direction and with the opposite direction (directions compared only, never effect sizes across scales; counts file `tables/",
      C, ".cross_tool.concordance.csv`). `no_direction` counts the genes called by both for which at least one tool has no direction.\n\n", sep = "")
  cat(knitr::kable(cc, format = "markdown"), sep = "\n")
  cat("\n\n")
}
```
````

The enrichment report `06_enrichment.Rmd` runs a GO enrichment (biological process, molecular function, cellular component) with clusterProfiler on the significant genes of each tool, for each contrast. It reads `rds/{CONTRAST}.{tool}.rds` (written by `01_import.Rmd`) and, for DEXSeq exon usage, edgeR exon usage and DTU, the same gene-level files of the run as `05_cross_tool.Rmd` (read only), and builds the same gene sets: rMATS and SUPPA2 genes with a significant event (adjusted p-value below `{FDR}`, absolute dPSI at least `{DPSI}`, and for rMATS not flagged `low_support`), DEXSeq exon usage genes with a `perGeneQValue` below `{FDR}`, edgeR exon usage genes with a Simes FDR below `{FDR}`, DTU genes with stageR's `gene` below `{FDR}` (re-run at the report's FDR when the run's alpha differs, as in `04_usage.Rmd`). Gene identifiers are Ensembl identifiers (versions removed, groups of overlapping genes split as in `05_cross_tool.Rmd`) and are mapped to Entrez identifiers with `org.Hs.eg.db` (`{ORGANISM}` human) or `org.Mm.eg.db` (mouse); the identifiers that cannot be mapped are counted and reported for every gene list and every universe. The universe of an enrichment is the set of genes that were testable in that tool (an adjusted p-value present in the harmonised table, or in the gene-level table), never the whole genome. There is one enrichment per tool and per group: the genes whose direction in the tool is `treatment_up`, those whose direction is `control_up` (the direction of a gene as defined in `05_cross_tool.Rmd`), and all significant genes; a group is tested only when it has at least 10 genes, and otherwise the report prints `skipped: fewer than 10 genes (n=...)`. Terms are kept at a Benjamini-Hochberg adjusted p-value and q-value of 0.05. It writes, per contrast, `tables/{CONTRAST}.enrich.{tool}.{group}.{go_bp,go_mf,go_cc}.csv` (`tool` as in the RDS names, `group` one of `all`, `treatment_up`, `control_up`), `tables/{CONTRAST}.enrich.summary.csv` (genes, mapped and unmapped counts, universe size, status) and `figures/{CONTRAST}.enrich.dotplot.png` (the top terms of the groups that ran; no figure when none ran). The report needs no network: it reads only the annotation packages of the image and never downloads or installs anything. `{ENRICH_DBS}` is `GO` unless the user asked for more; a request for KEGG or Reactome is answered in the report with `KEGG/Reactome require an image rebuild` (the packages ReactomePA, KEGG.db and reactome.db are not in the image), an organism other than human or mouse with `organism not supported`, and in both cases the report still renders and the job ends with exit code 0. A tool whose RDS or gene-level file is missing is named and left out.

**File `06_enrichment.Rmd`.**
````rmd
---
title: "GO enrichment of the significant genes"
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
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE, dev = "png", dpi = 150, fig.width = 8, fig.height = 5)
```

This report tests the significant genes of each tool for enriched Gene Ontology terms, for each contrast, from the tables that `01_import.Rmd` wrote to `r file.path(params$report_dir, "rds")` and from the gene-level files of the run (read only). Every effect is stated as **treatment minus control**: a positive effect is higher in the treatment group (more inclusion for rMATS and SUPPA2, higher usage for the exon and transcript usage tools). The direction of a gene in a tool is the sign of the effect of its most significant feature, as in `05_cross_tool.Rmd`: `treatment_up` for a positive effect, `control_up` for a negative one.

The gene lists are the genes called by each tool at FDR `r params$fdr` (rMATS and SUPPA2: a significant event with an absolute dPSI of at least `r params$dpsi`, and for rMATS not flagged `low_support`; DEXSeq exon usage: `perGeneQValue`; edgeR exon usage: Simes FDR; DTU: stageR's `gene`, re-run at the FDR of this report when the alpha of the run differs). Each tool is tested on its own, in three groups: the genes with direction `treatment_up`, the genes with direction `control_up`, and all significant genes. A group is tested only when it has at least 10 genes. The universe (the background) of each enrichment is the set of genes that were testable in that tool, the genes with an adjusted p-value in the harmonised table or in the gene-level table, and not the whole genome: a gene that the tool could not test cannot be a hit, so it is not in the background. Ensembl identifiers are mapped to Entrez identifiers with the annotation package of the organism (`r params$organism`); the identifiers that cannot be mapped are counted and left out. No network access is used: the report reads only the packages of the image.

```{r functions}
suppressPackageStartupMessages(library(ggplot2))
TOOL_LABELS <- c(rmats = "rMATS", suppa = "SUPPA2", exon_dexseq = "DEXSeq exon", exon_edger = "edgeR exon", dtu = "DEXSeq DTU")
MODULE_OF_TOOL <- c(rmats = "rmats", suppa = "suppa", exon_dexseq = "dexseq_exon", exon_edger = "edger_exon", dtu = "dtu")
GROUP_LABELS <- c(all = "all significant genes", treatment_up = "treatment_up (higher in treatment)", control_up = "control_up (higher in control)")
ONTOLOGIES <- c(go_bp = "BP", go_mf = "MF", go_cc = "CC")
ORGDB <- c(human = "org.Hs.eg.db", mouse = "org.Mm.eg.db")
MIN_GENES <- 10
ENRICH_CUTOFF <- 0.05

# the databases asked for in enrich_dbs (GO by default), lower case, separated by commas, semicolons or spaces
requested_dbs <- function(x) { v <- tolower(trimws(unlist(strsplit(as.character(x), "[,; ]+")))); v[nzchar(v)] }

# an Ensembl gene identifier loses its version suffix (.5) only when it has one
strip_version <- function(x) sub("^(ENS[A-Z]*G[0-9]+)\\.[0-9]+$", "\\1", x)
# a transcript identifier loses a trailing version in the same way (the pipeline strips the versions of the stageR identifiers, not those of the DEXSeq table: both sides of a join are stripped)
strip_tx_version <- function(x) sub("\\.[0-9]+$", "", x)
# a DEXSeq group of overlapping genes (GENE1+GENE2) or a SUPPA2 locus of overlapping genes (GENE1_and_GENE2) becomes one row per gene
expand_genes <- function(df) {
  parts <- strsplit(df$gene_id, "\\+|_and_")
  out <- df[rep(seq_len(nrow(df)), lengths(parts)), , drop = FALSE]
  out$gene_id <- as.character(unlist(parts))
  rownames(out) <- NULL
  out
}
tidy_ids <- function(df) { df$gene_id <- strip_version(df$gene_id); df }

# significant: adjusted p-value below the FDR; for events also the absolute dPSI threshold and not low_support
sig_usage <- function(tab, fdr) !is.na(tab$padj) & tab$padj < fdr
sig_events <- function(tab, fdr, dpsi) {
  !is.na(tab$padj) & tab$padj < fdr & !is.na(tab$effect) & abs(tab$effect) >= dpsi & !tab$low_support
}
# per gene: the significant features, the one with the smallest adjusted p-value (then the larger change, then the identifier) and its direction
top_by_gene <- function(feat, sig) {
  s <- feat[sig & !is.na(feat$padj), c("gene_id", "feature_id", "padj", "effect", "direction"), drop = FALSE]
  s <- s[order(s$padj, 1 - abs(s$effect), s$feature_id), , drop = FALSE]
  n <- table(s$gene_id)
  top <- s[!duplicated(s$gene_id), , drop = FALSE]
  data.frame(gene_id = top$gene_id, n_sig_features = as.integer(n[top$gene_id]), top_feature = top$feature_id, top_effect = top$effect,
             top_padj = top$padj, direction = ifelse(top$direction %in% c("treatment_up", "control_up"), top$direction, NA_character_),
             stringsAsFactors = FALSE)
}
# the calls of one tool: the genes called, the genes tested (the universe), and the direction of every gene that has a significant feature
make_calls <- function(genes, tested, feat, sig_feat) {
  list(genes = sort(unique(as.character(genes))), tested = sort(unique(as.character(tested))), all = sort(unique(as.character(feat$gene_id))),
       dir = top_by_gene(feat, sig_feat))
}
# rMATS and SUPPA2: a gene is called when one of its events is significant
calls_from_events <- function(tab, fdr, dpsi) {
  tab <- tidy_ids(expand_genes(tab))
  sig <- sig_events(tab, fdr, dpsi)
  make_calls(tab$gene_id[sig], tab$gene_id[!is.na(tab$padj)], tab, sig)
}
# DEXSeq and edgeR exon usage: a gene is called by the gene-level table (gene_id, padj); the features are the exon bins or exons of the RDS
calls_from_gene_table <- function(gene_tab, feat, fdr) {
  g <- tidy_ids(expand_genes(gene_tab))
  f <- tidy_ids(expand_genes(feat))
  make_calls(g$gene_id[sig_usage(g, fdr)], g$gene_id[!is.na(g$padj)], f, sig_usage(f, fdr))
}
# DTU: a gene is called when stageR's gene is below the FDR; its significant features are the transcripts with stageR's transcript below the FDR
calls_from_stager <- function(st, feat, fdr) {
  s <- tidy_ids(st)
  f <- tidy_ids(feat)
  s$tx_id <- strip_tx_version(s$tx_id)
  f$feature_id <- strip_tx_version(f$feature_id)
  txsig <- !is.na(s$transcript) & s$transcript < fdr
  confirmed <- paste(s$gene_id, s$tx_id)[txsig]
  sigf <- paste(f$gene_id, f$feature_id) %in% confirmed
  make_calls(s$gene_id[!is.na(s$gene) & s$gene < fdr], s$gene_id[!is.na(s$gene)], f, sigf)
}
# the three gene groups of a tool: the genes whose direction is treatment_up, those whose direction is control_up, and all called genes
group_genes <- function(cl) {
  d <- cl$dir
  list(all = cl$genes,
       treatment_up = intersect(cl$genes, d$gene_id[d$direction %in% "treatment_up"]),
       control_up = intersect(cl$genes, d$gene_id[d$direction %in% "control_up"]))
}

strip_quotes <- function(x) gsub("^\"|\"$", "", x)
read_delim_chr <- function(path, sep) {
  if (!file.exists(path)) stop("file ", path, " not found: none of its columns can be read")
  l <- sub("\r$", "", readLines(path, warn = FALSE))
  l <- l[nzchar(l)]
  if (length(l) == 0) stop("file ", path, " is empty: no header line, so no column can be found")
  nm <- strsplit(l[1], sep, fixed = TRUE)[[1]]
  if (length(l) == 1) {
    return(as.data.frame(stats::setNames(replicate(length(nm), character(0), simplify = FALSE), nm),
                         check.names = FALSE, stringsAsFactors = FALSE))
  }
  tab <- utils::read.table(text = l[-1], sep = sep, quote = "", header = FALSE, colClasses = "character", na.strings = character(0),
                           comment.char = "", col.names = paste0("V", seq_along(nm)), fill = FALSE, stringsAsFactors = FALSE)
  names(tab) <- nm
  tab
}
need_cols <- function(tab, cols, path) {
  miss <- setdiff(cols, names(tab))
  if (length(miss)) stop("file ", path, ": column not found: ", paste(miss, collapse = ", "),
                         " (columns present: ", paste(names(tab), collapse = ", "), ")")
}
to_num <- function(x, path, col) {
  x <- trimws(strip_quotes(x))
  miss <- x %in% c("", "NA", "nan", "NaN", "-nan")
  v <- suppressWarnings(as.numeric(ifelse(miss, NA, x)))
  bad <- is.na(v) & !miss
  if (any(bad)) stop("file ", path, ": column ", col, " has a value that is not a number: '", x[bad][1], "'")
  v
}
read_contrasts <- function(path) {
  tab <- read_delim_chr(path, ",")
  need_cols(tab, c("contrast", "treatment", "control"), path)
  tab <- tab[nzchar(tab$contrast), c("contrast", "treatment", "control")]
  rownames(tab) <- NULL
  tab
}
# a gene-level table: the gene identifier and the adjusted p-value (DEXSeq perGeneQValue: groupID, padj; edgeR Simes: Geneid, FDR)
read_gene_table <- function(path, sep, id_col, p_col) {
  tab <- read_delim_chr(path, sep)
  need_cols(tab, c(id_col, p_col), path)
  data.frame(gene_id = strip_quotes(tab[[id_col]]), padj = to_num(tab[[p_col]], path, p_col), stringsAsFactors = FALSE)
}
# stageR: gene and transcript are the screened adjusted p-values at the two levels, NA where not screened
read_stager <- function(path) {
  tab <- read_delim_chr(path, "\t")
  need_cols(tab, c("geneID", "txID", "gene", "transcript"), path)
  data.frame(gene_id = strip_quotes(tab$geneID), tx_id = strip_quotes(tab$txID), gene = to_num(tab$gene, path, "gene"),
             transcript = to_num(tab$transcript, path, "transcript"), stringsAsFactors = FALSE)
}
# the overall alpha stageR was run with, read from the stageRTx object with getAlpha; NA when the file is absent or cannot be read
stager_alpha <- function(path) {
  if (!file.exists(path) || !requireNamespace("stageR", quietly = TRUE)) return(NA_real_)
  tryCatch(as.numeric(stageR::getAlpha(readRDS(path)))[1], error = function(e) NA_real_)
}
# stageR's adjusted p-values are valid only at the alpha they were computed for. At another FDR the stage-wise adjustment is re-run
# from the stageRTx object of the run (stageR is in the image); when the object cannot be read the values of the run are used as they stand
stager_for_fdr <- function(st, path, fdr) {
  alpha <- stager_alpha(path)
  if (is.na(alpha)) return(list(st = st, alpha = NA_real_, mode = "unreadable"))
  if (isTRUE(all.equal(alpha, fdr))) return(list(st = st, alpha = alpha, mode = "run"))
  re <- tryCatch({
    obj <- stageR::stageWiseAdjustment(readRDS(path), method = "dtu", alpha = fdr, allowNA = TRUE)
    r <- stageR::getAdjustedPValues(obj, onlySignificantGenes = FALSE, order = FALSE)
    data.frame(gene_id = as.character(r$geneID), tx_id = as.character(r$txID), gene = as.numeric(r$gene), transcript = as.numeric(r$transcript),
               stringsAsFactors = FALSE)
  }, error = function(e) NULL)
  if (is.null(re)) return(list(st = st, alpha = alpha, mode = "unreadable"))
  list(st = re, alpha = alpha, mode = "rerun")
}
# the sentence that tells which stageR values are used
stager_note <- function(res, fdr) {
  if (res$mode == "run") return(paste0("stageR was run with an overall alpha of ", format(res$alpha), " (read from the `stageRTx` object of the run), equal to the FDR of this report: the stageR table of the run is used as it is."))
  if (res$mode == "rerun") return(paste0("stageR's adjusted p-values are valid only at the alpha they were computed for. The run used an overall alpha of ", format(res$alpha),
    " (read from the `stageRTx` object of the run) and the FDR of this report is ", format(fdr), ", so the stage-wise adjustment was re-run at ", format(fdr),
    " from that object with the stageR package of the image; the values of the run are not used."))
  paste0("**WARNING: the stageR values of the run are valid only at the alpha the run used, which is the run's setting and is not recorded in the table.** stageR's overall alpha is the run's setting; the `stageRTx` object was not found or could not be read, or the adjustment could not be re-run, so the values of the run are used as they stand for the FDR of this report, ",
         format(fdr), ", and the transcript-level values may not be adjusted p-values at that level.")
}
show_image <- function(file, alt) cat("![", alt, "](", file, ")\n\n", sep = "")

# the gene sets of every tool of one contrast, as in 05_cross_tool.Rmd: the calls of each tool that has a result, the files that were not found,
# and the sentence about the stageR values
collect_calls <- function(C, Tr, Co, modules, report_dir, outdir, fdr, dpsi) {
  path_of <- function(...) file.path(outdir, ...)
  calls <- list(); notfound <- character(0); stnote <- NULL
  for (tool in names(TOOL_LABELS)) {
    if (!MODULE_OF_TOOL[[tool]] %in% modules) next
    f <- file.path(report_dir, "rds", paste0(C, ".", tool, ".rds"))
    if (!file.exists(f)) { notfound <- c(notfound, f); next }
    tab <- readRDS(f)
    need_cols(tab, c("gene_id", "feature_id", "effect", "padj", "direction", "low_support"), f)
    cl <- NULL
    if (tool %in% c("rmats", "suppa")) {
      cl <- calls_from_events(tab, fdr, dpsi)
    } else if (tool == "exon_dexseq") {
      pgf <- path_of("star", "dexseq_exon", "results", sprintf("perGeneQValue.%s.csv", C))
      if (file.exists(pgf)) { pg <- read_gene_table(pgf, ",", "groupID", "padj"); cl <- calls_from_gene_table(pg, tab, fdr) } else notfound <- c(notfound, pgf)
    } else if (tool == "exon_edger") {
      sf <- path_of("star", "edger", sprintf("contrast_%s.usage.simes.csv", C))
      if (file.exists(sf)) { sm <- read_gene_table(sf, ",", "Geneid", "FDR"); cl <- calls_from_gene_table(sm, tab, fdr) } else notfound <- c(notfound, sf)
    } else if (tool == "dtu") {
      gf <- path_of("salmon", "dexseq_dtu", "results", "stager", sprintf("getAdjustedPValues.%s-%s.tsv", Tr, Co))
      if (file.exists(gf)) {
        st <- read_stager(gf)
        sres <- stager_for_fdr(st, path_of("salmon", "dexseq_dtu", "results", "stager", sprintf("stageRTx.%s-%s.rds", Tr, Co)), fdr)
        stnote <- stager_note(sres, fdr)
        cl <- calls_from_stager(sres$st, tab, fdr)
      } else notfound <- c(notfound, gf)
    }
    if (!is.null(cl)) calls[[tool]] <- cl
  }
  list(calls = calls, notfound = notfound, stnote = stnote)
}

# the annotation object of an organism package (org.Hs.eg.db holds the object org.Hs.eg.db)
orgdb_object <- function(pkg) getExportedValue(pkg, pkg)
# Ensembl identifiers to Entrez identifiers; the identifiers without an Entrez identifier are counted as unmapped
to_entrez <- function(ids, pkg) {
  ids <- unique(as.character(ids)); ids <- ids[!is.na(ids) & nzchar(ids)]
  db <- orgdb_object(pkg)
  known <- ids[ids %in% AnnotationDbi::keys(db, keytype = "ENSEMBL")]
  m <- if (length(known)) suppressMessages(AnnotationDbi::select(db, keys = known, keytype = "ENSEMBL", columns = "ENTREZID")) else
    data.frame(ENSEMBL = character(0), ENTREZID = character(0), stringsAsFactors = FALSE)
  m <- m[!is.na(m$ENTREZID), , drop = FALSE]
  mapped <- length(unique(m$ENSEMBL))
  list(entrez = sort(unique(as.character(m$ENTREZID))), n_ids = length(ids), n_mapped = mapped, n_unmapped = length(ids) - mapped)
}
GO_COLUMNS <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
# one GO enrichment (ont is BP, MF or CC) of the Entrez genes in the universe; terms at the cutoff, Benjamini-Hochberg; a failure gives no term and a message
run_go <- function(entrez, universe, pkg, ont) {
  err <- NULL
  res <- tryCatch(clusterProfiler::enrichGO(gene = entrez, universe = universe, OrgDb = orgdb_object(pkg), keyType = "ENTREZID", ont = ont,
                                            pAdjustMethod = "BH", pvalueCutoff = ENRICH_CUTOFF, qvalueCutoff = ENRICH_CUTOFF, readable = FALSE),
                  error = function(e) { err <<- conditionMessage(e); NULL })
  d <- if (is.null(res)) NULL else as.data.frame(res)
  if (is.null(d) || nrow(d) == 0) d <- as.data.frame(stats::setNames(replicate(length(GO_COLUMNS), character(0), simplify = FALSE), GO_COLUMNS), stringsAsFactors = FALSE)
  d <- d[, GO_COLUMNS, drop = FALSE]
  rownames(d) <- NULL
  list(table = d, error = err)
}
ratio_of <- function(x) as.numeric(sub("/.*$", "", x)) / as.numeric(sub("^.*/", "", x))
# the top terms of every group that ran (n per group and ontology), as one table for the dotplot
dotplot_data <- function(res, n = 5) {
  if (nrow(res) == 0) return(res)
  parts <- split(res, paste(res$run, res$ont))
  out <- do.call(rbind, lapply(parts, function(d) head(d[order(d$p.adjust, d$ID), , drop = FALSE], n)))
  out$ratio <- ratio_of(out$GeneRatio)
  out$label <- ifelse(nchar(out$Description) > 60, paste0(substr(out$Description, 1, 57), "..."), out$Description)
  rownames(out) <- NULL
  out
}
enrich_dotplot <- function(res, file) {
  d <- dotplot_data(res)
  d$label <- factor(d$label, levels = rev(unique(d$label[order(d$ont, d$p.adjust)])))
  p <- ggplot(d, aes(x = run, y = label, size = Count, colour = p.adjust)) + geom_point() +
    facet_grid(ont ~ ., scales = "free_y", space = "free_y") +
    scale_colour_gradient(low = "#b2182b", high = "#2166ab") +
    labs(x = NULL, y = NULL, colour = "adjusted p-value", size = "genes") + theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggplot2::ggsave(filename = file, plot = p, width = 9, height = max(4, 0.25 * nlevels(d$label) + 2), dpi = 150, limitsize = FALSE)
  invisible(file)
}
```

The enrichment is run only for human and mouse (the image holds the annotation packages of these two organisms), with the Gene Ontology through clusterProfiler. A term is kept at a Benjamini-Hochberg adjusted p-value and q-value of `r ENRICH_CUTOFF`, whatever the FDR used to call the genes, and a group needs at least `r MIN_GENES` genes.

```{r report, results="asis"}
modules <- strsplit(trimws(params$modules), "[[:space:]]+")[[1]]
fig_dir <- file.path(params$report_dir, "figures"); tab_dir <- file.path(params$report_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE); dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
contrasts <- read_contrasts(params$contrasts)
dbs <- requested_dbs(params$enrich_dbs)
pkg <- unname(ORGDB[params$organism])
can_run <- isTRUE(params$enrich)
if (!can_run) {
  cat("The enrichment is switched off for this report (enrich = FALSE).\n\n")
} else {
  if (any(c("kegg", "reactome") %in% dbs)) {
    if (all(vapply(c("ReactomePA", "KEGG.db", "reactome.db"), requireNamespace, logical(1), quietly = TRUE))) {
      cat("KEGG or Reactome was requested in `enrich_dbs`: this report runs GO only.\n\n")
    } else {
      cat("KEGG or Reactome was requested in `enrich_dbs`: KEGG/Reactome require an image rebuild (ReactomePA, KEGG.db and reactome.db are not in this image). GO only.\n\n")
    }
  }
  if (is.na(pkg)) {
    cat("**organism not supported**: the organism is `", params$organism, "`, and the image holds annotation packages for human and mouse only. No enrichment is run.\n\n", sep = "")
    can_run <- FALSE
  } else if (!all(vapply(c(pkg, "clusterProfiler", "GO.db", "AnnotationDbi"), requireNamespace, logical(1), quietly = TRUE))) {
    cat("The packages `", pkg, "`, `clusterProfiler`, `GO.db` and `AnnotationDbi` are not all in this image: no enrichment is run (an image rebuild is needed).\n\n", sep = "")
    can_run <- FALSE
  } else if (!"go" %in% dbs) {
    cat("GO is not in `enrich_dbs` (`", params$enrich_dbs, "`): this report has nothing to run.\n\n", sep = "")
    can_run <- FALSE
  }
}
for (i in seq_len(nrow(contrasts))) {
  if (!can_run) break
  C <- contrasts$contrast[i]; Tr <- contrasts$treatment[i]; Co <- contrasts$control[i]
  cat("\n\n## Contrast ", C, "\n\n", sep = "")
  cat("Treatment `", Tr, "`, control `", Co, "`: `treatment_up` genes have a higher effect in `", Tr, "`, `control_up` genes in `", Co, "`.\n\n", sep = "")
  cat("Thresholds: FDR ", format(params$fdr), ", abs(dPSI) ", format(params$dpsi), " (rMATS and SUPPA2 events only); rMATS minimum read support ", params$min_reads,
      ". Enrichment: GO (BP, MF, CC), Benjamini-Hochberg cutoff ", format(ENRICH_CUTOFF), ", groups of at least ", MIN_GENES, " genes.\n\n", sep = "")
  cc <- collect_calls(C, Tr, Co, modules, params$report_dir, params$outdir, params$fdr, params$dpsi)
  if (length(cc$notfound)) cat("Files not found (the tool is left out):\n\n", paste0("- `", cc$notfound, "`\n"), "\n", sep = "")
  if (!is.null(cc$stnote)) cat("DTU, stageR: ", cc$stnote, "\n\n", sep = "")
  if (length(cc$calls) == 0) { cat("No tool has a result for this contrast: there is no gene list to test.\n\n"); next }
  summ <- list(); allres <- list()
  for (tool in names(cc$calls)) {
    cl <- cc$calls[[tool]]
    cat("### ", TOOL_LABELS[[tool]], ", ", C, "\n\n", sep = "")
    uni <- to_entrez(cl$tested, pkg)
    cat("Universe: the ", uni$n_ids, " genes testable in ", TOOL_LABELS[[tool]], " (an adjusted p-value is present); ", uni$n_mapped, " mapped to Entrez, ",
        uni$n_unmapped, " unmapped.\n\n", sep = "")
    gs <- group_genes(cl)
    for (g in names(gs)) {
      ids <- gs[[g]]
      row <- data.frame(tool = tool, group = g, genes = length(ids), genes_mapped = NA_integer_, genes_unmapped = NA_integer_, universe_genes = uni$n_ids,
                        universe_mapped = uni$n_mapped, universe_unmapped = uni$n_unmapped, status = "skipped", terms_go_bp = NA_integer_,
                        terms_go_mf = NA_integer_, terms_go_cc = NA_integer_, stringsAsFactors = FALSE)
      cat("#### ", GROUP_LABELS[[g]], "\n\n", sep = "")
      if (length(ids) < MIN_GENES) {
        cat("skipped: fewer than ", MIN_GENES, " genes (n=", length(ids), ")\n\n", sep = "")
        summ[[length(summ) + 1]] <- row
        next
      }
      m <- to_entrez(ids, pkg)
      row$genes_mapped <- m$n_mapped; row$genes_unmapped <- m$n_unmapped
      cat(length(ids), " genes, ", m$n_mapped, " mapped to Entrez, ", m$n_unmapped, " unmapped.\n\n", sep = "")
      if (length(m$entrez) < MIN_GENES) {
        cat("skipped: fewer than ", MIN_GENES, " genes (n=", length(m$entrez), " with an Entrez identifier)\n\n", sep = "")
        summ[[length(summ) + 1]] <- row
        next
      }
      row$status <- "ran"
      for (k in names(ONTOLOGIES)) {
        r <- run_go(m$entrez, uni$entrez, pkg, ONTOLOGIES[[k]])
        utils::write.csv(r$table, file.path(tab_dir, paste0(C, ".enrich.", tool, ".", g, ".", k, ".csv")), row.names = FALSE)
        row[[paste0("terms_", k)]] <- nrow(r$table)
        cat("GO ", ONTOLOGIES[[k]], ": ", nrow(r$table), " terms at the cutoff (table `tables/", C, ".enrich.", tool, ".", g, ".", k, ".csv`)", sep = "")
        if (!is.null(r$error)) cat(" (the enrichment failed: ", r$error, ")", sep = "")
        cat("\n\n")
        if (nrow(r$table) > 0) {
          show <- head(r$table[, c("ID", "Description", "GeneRatio", "BgRatio", "p.adjust", "Count")], 10)
          show$p.adjust <- signif(show$p.adjust, 3)
          cat(knitr::kable(show, format = "markdown", row.names = FALSE), sep = "\n")
          cat("\n\n")
          allres[[length(allres) + 1]] <- data.frame(run = paste0(TOOL_LABELS[[tool]], "\n", g), ont = ONTOLOGIES[[k]], r$table[, c("ID", "Description", "GeneRatio", "p.adjust", "Count")], stringsAsFactors = FALSE)
        }
      }
      summ[[length(summ) + 1]] <- row
    }
  }
  sm <- do.call(rbind, summ); rownames(sm) <- NULL
  utils::write.csv(sm, file.path(tab_dir, paste0(C, ".enrich.summary.csv")), row.names = FALSE)
  cat("### Summary, ", C, "\n\n", sep = "")
  cat("Genes per group, genes mapped and unmapped, universe (testable genes) and status (table `tables/", C, ".enrich.summary.csv`):\n\n", sep = "")
  cat(knitr::kable(sm, format = "markdown"), sep = "\n")
  cat("\n\n### Dotplot, ", C, "\n\n", sep = "")
  res <- if (length(allres)) do.call(rbind, allres) else data.frame()
  if (nrow(res) > 0) {
    df <- file.path(fig_dir, paste0(C, ".enrich.dotplot.png"))
    enrich_dotplot(res, df)
    cat("The top terms of the groups that ran (up to 5 per group and ontology; figure `figures/", C, ".enrich.dotplot.png`).\n\n", sep = "")
    show_image(df, paste0("GO enrichment, top terms, ", C))
  } else {
    cat("No group ran, or no term passed the cutoff: no dotplot for this contrast.\n\n")
  }
}
```
````

The summary page `07_summary.Rmd` is the index of the report directory. It runs last and starts when the other jobs have ended, whatever their outcome, so some reports may be missing. It scans `{REPORT_DIR}` and lists, with links, the reports `01_import.html` to `06_enrichment.html` that exist (a report that does not exist is listed as `not produced` and is not linked), every figure in `figures/` and every table in `tables/` that exists (nothing is listed that is not in the directory), the thresholds (`{FDR}`, `{DPSI}`, `{RMATS_TYPE}`, `{MIN_READS}`), the contrasts with their treatment and control, the modules found and the sign rules. For each contrast it shows the counts tables that the other reports wrote (`tables/{CONTRAST}.rmats.counts.csv`, `.suppa.counts.csv`, `.usage.counts.csv`), the cross-tool tables (`.cross_tool.sets.csv`, `.cross_tool.concordance.csv`) and the enrichment summary (`.enrich.summary.csv`); a table that does not exist is reported as not found and the page still renders. A Warnings block collects the lines that begin with `WARNING` that the other reports print into their text, recovered from the HTML of those reports; warnings that R raised and that knitr suppressed are not in the HTML and are not listed. The page is written to `report_summary.html`, the entry point of the report directory.

**File `07_summary.Rmd`.**
````rmd
---
title: "Summary of the nf-core/rnasplice reports"
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

This page lists the reports, figures and tables of the report directory `r params$report_dir` for the nf-core/rnasplice run in `r params$outdir`, the thresholds, the contrasts, the sign rules, the counts of every contrast and the warnings the reports printed. Only files that exist are listed. Every effect is stated as **treatment minus control**: a positive effect is higher in the treatment group (more inclusion for rMATS and SUPPA2, higher usage for the exon and transcript usage tools).

```{r functions}
REPORTS <- data.frame(
  file = c("01_import.html", "02_rmats.html", "03_suppa.html", "04_usage.html", "05_cross_tool.html", "06_enrichment.html"),
  title = c("Import and harmonise", "rMATS events", "SUPPA2 events", "Exon and transcript usage", "Cross-tool comparison", "GO enrichment"),
  stringsAsFactors = FALSE)
# the counts tables that the other reports write for each contrast: file part, title, and the module it belongs to (NA: not tied to one module)
COUNT_TABLES <- data.frame(
  part = c("rmats.counts", "suppa.counts", "usage.counts", "cross_tool.sets", "cross_tool.concordance", "enrich.summary"),
  title = c("rMATS events and significant events", "SUPPA2 events and significant events", "Exon and transcript usage: tested and significant",
            "Cross-tool: genes tested and called by each tool", "Cross-tool: direction concordance for each pair of tools", "GO enrichment: groups, genes and universe"),
  module = c("rmats", "suppa", "usage", NA, NA, NA), stringsAsFactors = FALSE)
USAGE_MODULES <- c("dexseq_exon", "edger_exon", "dtu")

html_escape <- function(x) { x <- gsub("&", "&amp;", x, fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE); x <- gsub(">", "&gt;", x, fixed = TRUE); gsub("\"", "&quot;", x, fixed = TRUE) }
# the files of a directory that match a pattern, sorted; nothing when the directory does not exist
list_files <- function(dir, pattern) if (dir.exists(dir)) sort(list.files(dir, pattern = pattern)) else character(0)
# a table of the report directory with every field kept as text (what the other report wrote is shown as it stands); NULL when it is missing or unreadable
read_counts <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(utils::read.csv(path, colClasses = "character", check.names = FALSE, na.strings = character(0), stringsAsFactors = FALSE), error = function(e) NULL)
}
# the lines that begin with WARNING in the text of an HTML report: tags removed, entities decoded, each line once
html_warnings <- function(path) {
  if (!file.exists(path)) return(character(0))
  h <- paste(readLines(path, warn = FALSE), collapse = "\n")
  h <- gsub("(?s)<(script|style)[^>]*>.*?</\\1>", " ", h, perl = TRUE)
  h <- gsub("</(p|li|tr|h[1-6]|div|pre)>|<br ?/?>", "\n", h, perl = TRUE)
  h <- gsub("<[^>]+>", "", h)
  h <- gsub("&amp;", "&", gsub("&lt;", "<", gsub("&gt;", ">", gsub("&quot;", "\"", h, fixed = TRUE), fixed = TRUE), fixed = TRUE), fixed = TRUE)
  l <- trimws(gsub("[[:space:]]+", " ", strsplit(h, "\n", fixed = TRUE)[[1]]))
  l <- sub("^.*?(WARNING)", "\\1", l[grepl("WARNING", l, fixed = TRUE)], perl = TRUE)
  unique(l)
}
strip_quotes <- function(x) gsub("^\"|\"$", "", x)
read_delim_chr <- function(path, sep) {
  if (!file.exists(path)) stop("file ", path, " not found: none of its columns can be read")
  l <- sub("\r$", "", readLines(path, warn = FALSE))
  l <- l[nzchar(l)]
  if (length(l) == 0) stop("file ", path, " is empty: no header line, so no column can be found")
  nm <- strsplit(l[1], sep, fixed = TRUE)[[1]]
  if (length(l) == 1) {
    return(as.data.frame(stats::setNames(replicate(length(nm), character(0), simplify = FALSE), nm),
                         check.names = FALSE, stringsAsFactors = FALSE))
  }
  tab <- utils::read.table(text = l[-1], sep = sep, quote = "", header = FALSE, colClasses = "character", na.strings = character(0),
                           comment.char = "", col.names = paste0("V", seq_along(nm)), fill = FALSE, stringsAsFactors = FALSE)
  names(tab) <- nm
  tab
}
need_cols <- function(tab, cols, path) {
  miss <- setdiff(cols, names(tab))
  if (length(miss)) stop("file ", path, ": column not found: ", paste(miss, collapse = ", "),
                         " (columns present: ", paste(names(tab), collapse = ", "), ")")
}
read_contrasts <- function(path) {
  tab <- read_delim_chr(path, ",")
  need_cols(tab, c("contrast", "treatment", "control"), path)
  tab <- tab[nzchar(tab$contrast), c("contrast", "treatment", "control")]
  rownames(tab) <- NULL
  tab
}
# a list of links to files of the report directory, as HTML (the paths are relative to the directory)
link_list <- function(rel, label = rel) {
  if (length(rel) == 0) return("none\n")
  paste0("<ul>\n", paste0("<li><a href=\"", html_escape(rel), "\">", html_escape(label), "</a></li>", collapse = "\n"), "\n</ul>\n")
}
```

```{r report, results="asis"}
modules <- strsplit(trimws(params$modules), "[[:space:]]+")[[1]]
modules <- modules[nzchar(modules)]
rd <- params$report_dir
cat("## Reports\n\n")
cat("The reports of `", rd, "`. A report missing from the directory is marked in the list (its job failed or did not run, or the section was not requested).\n\n", sep = "")
for (i in seq_len(nrow(REPORTS))) {
  f <- REPORTS$file[i]
  if (file.exists(file.path(rd, f))) cat("- <a href=\"", html_escape(f), "\">", html_escape(f), "</a>: ", REPORTS$title[i], "\n", sep = "")
  else cat("- ", REPORTS$title[i], ": not produced\n", sep = "")
}
cat("\n\n## Thresholds and run\n\n")
cat("- FDR (adjusted p-value cutoff): ", format(params$fdr), "\n", sep = "")
cat("- dPSI (minimum absolute dPSI of an rMATS or SUPPA2 event; not applied to exon or transcript usage): ", format(params$dpsi), "\n", sep = "")
cat("- rMATS count type: ", params$rmats_type, "\n", sep = "")
cat("- rMATS minimum read support (mean of IJC plus SJC per group, smaller of the two groups): ", params$min_reads, "\n", sep = "")
cat("- Pairing: rMATS ", if (isTRUE(params$paired_rmats)) "paired" else "unpaired", " statistics, SUPPA2 ", if (isTRUE(params$paired_suppa)) "paired" else "unpaired", " statistics (read from the run's parameters)\n", sep = "")
cat("- Organism: ", params$organism, "; enrichment databases requested: ", params$enrich_dbs, "\n", sep = "")
cat("- Modules found: ", if (length(modules)) paste(modules, collapse = ", ") else "none", "\n", sep = "")
cat("- Run output (read only): `", params$outdir, "`\n\n", sep = "")
cat("## Contrasts\n\n")
contrasts <- read_contrasts(params$contrasts)
cat("Each contrast is treatment against control; every sign below is fixed by this choice.\n\n")
cat(knitr::kable(contrasts, format = "markdown"), sep = "\n")
cat("\n\n## Sign rules\n\n")
cat("| Tool | Column read | Convention of the tool | Conversion to treatment minus control |\n|---|---|---|---|\n")
cat("| rMATS | `IncLevelDifference` | IncLevel1 minus IncLevel2, group 1 (b1) is the treatment: treatment minus control | kept |\n")
cat("| SUPPA2 | `{kind}_{T}-{kind}_{C}_dPSI` | control minus treatment | sign flipped |\n")
cat("| DEXSeq exon usage | `log2fold_{C}_{T}` | control versus treatment | sign flipped |\n")
cat("| DEXSeq transcript usage (DTU) | `log2fold_{C}_{T}` | control versus treatment | sign flipped |\n")
cat("| edgeR exon usage | `logFC` | treatment minus control | kept |\n\n")
cat("`direction` is `treatment_up` for a positive effect, `control_up` for a negative one and `none` for zero or a missing effect.\n\n")
cat("## Counts per contrast\n\n")
for (i in seq_len(nrow(contrasts))) {
  C <- contrasts$contrast[i]
  cat("### ", C, "\n\n", sep = "")
  cat("Treatment `", contrasts$treatment[i], "`, control `", contrasts$control[i], "`.\n\n", sep = "")
  for (j in seq_len(nrow(COUNT_TABLES))) {
    rel <- file.path("tables", paste0(C, ".", COUNT_TABLES$part[j], ".csv"))
    d <- read_counts(file.path(rd, rel))
    cat("**", COUNT_TABLES$title[j], "** (`", rel, "`)\n\n", sep = "")
    if (is.null(d)) {
      mod <- COUNT_TABLES$module[j]
      absent <- !is.na(mod) && !any((if (mod == "usage") USAGE_MODULES else mod) %in% modules)
      cat(if (absent) "not found: the module is not in the run.\n\n" else "not found: the report that writes it was not produced, or the table is empty or unreadable.\n\n")
    } else if (nrow(d) == 0) {
      cat("The table has no rows.\n\n")
    } else {
      cat(knitr::kable(d, format = "markdown", row.names = FALSE), sep = "\n")
      cat("\n\n")
    }
  }
}
cat("## Figures\n\n")
figs <- list_files(file.path(rd, "figures"), "\\.(png|pdf|svg|jpg)$")
cat(link_list(file.path("figures", figs), figs), "\n")
cat("## Tables\n\n")
tabs <- list_files(file.path(rd, "tables"), "\\.(csv|tsv|txt)$")
cat(link_list(file.path("tables", tabs), tabs), "\n")
cat("## Warnings\n\n")
cat("Lines that begin with WARNING in the text of the other reports, recovered from the HTML of those reports (a warning that R raised and that knitr suppressed is not in the HTML and is not listed here).\n\n")
nw <- 0
for (i in seq_len(nrow(REPORTS))) {
  w <- html_warnings(file.path(rd, REPORTS$file[i]))
  if (length(w)) { nw <- nw + length(w); cat("**", REPORTS$file[i], "**\n\n", paste0("- ", w, "\n"), "\n", sep = "") }
}
if (nw == 0) cat("No WARNING line was found in the reports that exist.\n\n")
```
````
