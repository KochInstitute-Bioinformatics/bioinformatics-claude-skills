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
- `{MIN_READS}`: minimum read support of an rMATS event: for each group, the mean over its samples of `IJC` plus `SJC`, required in both groups (the smaller of the two group means must reach it). Default: 10. Events below it are flagged (`low_support`), kept in the tables and excluded from the significance counts of the reports, never deleted; their number is reported, and the rMATS FDR is not recomputed (it stays as rMATS wrote it). On very shallow data (like the nf-core test data) most events fall below 10; the user may lower it. Asked only when `rmats` is in `{MODULES}`.

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

The rMATS report `02_rmats.Rmd` reads `rds/{CONTRAST}.rmats.rds` (written by `01_import.Rmd`) for every contrast and writes, per contrast, `tables/{CONTRAST}.rmats.counts.csv`, `tables/{CONTRAST}.rmats.{type}.significant.csv`, `figures/{CONTRAST}.rmats.volcano.{type}.png`, `figures/{CONTRAST}.rmats.effect_hist.png` and `figures/{CONTRAST}.rmats.psi_top.png` (`type` is one of `SE`, `A5SS`, `A3SS`, `MXE`, `RI`). An event is significant when its adjusted p-value is below `{FDR}`, its absolute dPSI is at least `{DPSI}`, and it is not flagged `low_support`. The per-sample PSI of the top events is not in the RDS; the report reads `IncLevel1` and `IncLevel2` again from the rMATS table of the same directory the import used (the contrast's own, or the one with the suffix `_paired`). A contrast whose RDS is missing is reported as missing and the report still renders.

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
# one row per event and sample. IncLevel1 lists the samples of group 1 (the treatment), IncLevel2 those of group 2 (the control),
# each in the order of the sample sheet
psi_long <- function(inc, ss, treatment, control) {
  s1 <- ss$sample[ss$condition == treatment]
  s2 <- ss$sample[ss$condition == control]
  rows <- lapply(seq_len(nrow(inc)), function(i) {
    v1 <- psi_values(inc$IncLevel1[i]); v2 <- psi_values(inc$IncLevel2[i])
    if (length(v1) != length(s1) || length(v2) != length(s2))
      stop("event ", inc$feature_id[i], ": IncLevel1 and IncLevel2 list ", length(v1), " and ", length(v2), " samples, the sample sheet has ",
           length(s1), " samples of ", treatment, " and ", length(s2), " samples of ", control)
    data.frame(feature_id = inc$feature_id[i], sample = c(s1, s2), condition = rep(c(treatment, control), c(length(s1), length(s2))),
               role = rep(c("treatment", "control"), c(length(s1), length(s2))), psi = c(v1, v2), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
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
    labs(title = paste0("Per-sample PSI of the top rMATS events, ", contrast), x = "sample (order of the sample sheet)", y = "PSI (IncLevel)") +
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
      cat("The six events that are significant, or nearest to significance when fewer than six are (events flagged low_support are left out). ",
          "IncLevel1 holds the samples of the treatment `", Tr, "` and IncLevel2 those of the control `", Co,
          "`, read again from the tables in `", rdir, "`; the samples are in the order of the sample sheet.\n\n", sep = "")
      pf <- file.path(fig_dir, paste0(C, ".rmats.psi_top.png"))
      save_plot(psi_plot(psi_long(inc, ss, Tr, Co), top[, c("feature_id", "gene_symbol")], Tr, Co, C), pf, width = 9, height = 6)
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
