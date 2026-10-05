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
- `{MIN_READS}`: minimum read support of an rMATS event, the sum of `IJC` and `SJC` over the samples of a group, required in both groups. Default: 10. Events below it are dropped from the tests, and their number is reported. Asked only when `rmats` is in `{MODULES}`.

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
