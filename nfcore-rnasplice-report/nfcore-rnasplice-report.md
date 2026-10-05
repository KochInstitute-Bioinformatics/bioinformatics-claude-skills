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

**Module detection.** A module counts as found only when its main result file exists and has more than a header line; an empty file or a file with only a header is not found. The files are named from the contrasts file (treatment `{T}`, control `{C}`, contrast name `{CONTRAST}`), never by splitting file names: rMATS `star/rmats/{CONTRAST}/{SE,A5SS,A3SS,MXE,RI}.MATS.{JC,JCEC}.txt`; SUPPA2 `salmon/suppa/diffsplice/per_local_event/local_{T}-{C}.dpsi`; DEXSeq exon `star/dexseq_exon/results/DEXSeqResults.{CONTRAST}.csv`; edgeR exon `star/edger/contrast_{CONTRAST}.usage.exon.csv`; DTU `salmon/dexseq_dtu/results/dexseq/DEXSeqResults.{T}-{C}.tsv`. SUPPA2 and DTU file names use `{T}-{C}` (treatment, hyphen, control); the other three use the contrast name. A BAM-input run has no `salmon/` directory, so SUPPA2 and DTU are simply not found. The rMATS tables have the column `ID` twice (first column, and again after the coordinate columns); detection only counts lines and never reads columns, so this does not matter here, but every reader of these tables must account for it. Run, in one Bash call (replace `{OUTDIR_RNASPLICE}` and `{CONTRASTS_CSV}` with the literal paths):
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
  for ty in SE A5SS A3SS MXE RI; do for m in JC JCEC; do has_rows "$OUT/star/rmats/$c/$ty.MATS.$m.txt" && hit="$hit $ty.$m"; done; done
  if [ -n "$hit" ]; then echo "found rmats $c$hit"; add rmats; else echo "missing rmats $c $OUT/star/rmats/$c/ (no MATS table with data)"; fi
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
