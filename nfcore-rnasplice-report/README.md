# `/nfcore-rnasplice-report` — nf-core/rnasplice Alternative-Splicing Report Skill

An interactive Claude Code skill that turns the output of an [nf-core/rnasplice](https://nf-co.re/rnasplice) run (set up with `/nfcore-rnasplice-setup`) into R Markdown reports and SLURM scripts: rMATS splicing events, SUPPA2 event and isoform PSI, DEXSeq and edgeR differential exon usage, and DEXSeq differential transcript usage (stageR), a cross-tool comparison, GO enrichment and one summary page. The report adapts to the modules that produced output. Every effect direction in every table and plot is stated as treatment minus control. The skill reads the run's tables and writes only to its own report directory; it never runs the pipeline.

---

## Installation

```bash
cp nfcore-rnasplice-report.md ~/.claude/commands/
```

Then invoke it in Claude Code with `/nfcore-rnasplice-report`, from the project directory that holds the run's sample sheet, contrasts sheet and params file.

---

## Prerequisites

| Requirement | Notes |
|-------------|-------|
| SLURM scheduler | Scripts target the `bcc` partition; R runs only on compute nodes |
| Singularity | Loaded with `module add singularity/3.10.4` |
| R image | `/net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif`. Mandatory packages: DEXSeq, edgeR, stageR, ggplot2, dplyr, tidyr, readr, knitr, rmarkdown. Optional, used when present: UpSetR or ComplexUpset (cross-tool plot), clusterProfiler with org.Hs.eg.db or org.Mm.eg.db (GO enrichment), enrichplot, ggrepel, scales. The wizard checks the packages in a compute-node job and stops when a mandatory one is missing; packages are never installed at run time |
| Internet | Not needed for the report: the Rmd files make no network call and install nothing. Only pulling a missing image (a compute-node job) needs it |
| A finished nf-core/rnasplice run | At least one of rMATS, SUPPA2, DEXSeq exon, edgeR exon or DEXSeq DTU output with more than a header line |

---

## What the skill does

The wizard asks one question at a time and states the default with each:

| Step | Content |
|------|---------|
| 0 | Working directory, report directory `results/{TODAY_ISO}_{WD_NAME}_report`; an existing directory is never overwritten silently |
| 1 | Email for SLURM notifications (default: none) |
| 2 | R image and its package check (compute node) |
| 3 | Locate the run and detect the modules: a module counts as found only when its main result file exists and has more than a header |
| 4 | Read the sample sheet and contrasts sheet; the user confirms the treatment and the control of each contrast, which fixes every sign |
| 5 | Annotation (GTF) and organism, default from the run's params |
| 6 | Pairing, read from the run's params (`rmats_paired_stats`, `diffsplice_paired`) and printed; the report never claims unpaired statistics for a paired run |
| 7 | Thresholds: FDR, dPSI, rMATS count type, minimum read support |
| 8 | Analyses menu: `events`, `cross`, `topgene`, `enrich` (only what the run supports is offered) |
| 9 | Memory (default 32G), time (default 4:00:00), project name |
| 10 | Writes the Rmd files and SLURM scripts and prints one submit command |

Several contrasts are handled one by one: file names, tables and figures carry the contrast name, and contrasts are never mixed. File and column names are built from the contrasts sheet, never recovered by splitting on `-` or `_`.

### Inputs and outputs

The report directory is `results/{TODAY_ISO}_{WD_NAME}_report`, next to (never inside) the pipeline's own output directory `results/{TODAY_ISO}_{WD_NAME}`. The output of the run is read-only.

| Rmd | Job | Content | Generated when |
|-----|-----|---------|----------------|
| `01_import.Rmd` | `run_01.sbatch` | Reads every module's tables, harmonises them (one table per tool and contrast, saved as RDS), counts, quality notes | always |
| `02_rmats.Rmd` | `run_02.sbatch` | rMATS events by type, volcano plots, PSI of the top events per sample | `rmats` found |
| `03_suppa.Rmd` | `run_03.sbatch` | SUPPA2 local events and isoforms, volcano plots | `suppa` found |
| `04_usage.Rmd` | `run_04.sbatch` | DEXSeq and edgeR exon usage, DEXSeq DTU with stageR, top-gene plots | any of `dexseq_exon`, `edger_exon`, `dtu` found |
| `05_cross_tool.Rmd` | `run_05.sbatch` | UpSet plot of the significant genes per tool, direction concordance | `cross` chosen and two or more tools found |
| `06_enrichment.Rmd` | `run_06.sbatch` | GO enrichment (BP, MF, CC) per tool and direction | `enrich` chosen |
| `07_summary.Rmd` | `run_07.sbatch` | Summary page that indexes every report, figure and table | always, last |

`run_all.sh` only calls sbatch: jobs 02 to 06 start when job 01 has finished successfully, job 07 when the others have ended, whatever their outcome. Start with **`report_summary.html`** in the report directory: it is the entry page, written by job 07. Tables go to `tables/`, figures to `figures/`, the harmonised tables to `rds/`.

### How the scripts are made

Step 10 cuts the seven Rmd files out of the installed skill file into the report directory, writes two small template files to `templates/` there (`run_job.sbatch.tpl` and `run_all.tpl`), and renders the job scripts and `run_all.sh` from them with three short Bash blocks (checks and job list, settings, rendering), so that no block is long or nests a heredoc. The wizard refuses typed values (e-mail address, paths, project name) that contain a single quote, a double quote, a backslash or a line break; a `$`, a backtick or `$(...)` in a value is kept literally and no command in it is run.

### Harmonised-table contract

Rmd 01 turns each tool's table into one common table per contrast and tool, `rds/{CONTRAST}.{tool}.rds` (tool: `rmats`, `suppa`, `exon_dexseq`, `exon_edger`, `dtu`), with the columns `contrast, tool, gene_id, gene_symbol, feature_type, feature_id, effect, effect_scale, padj, direction, n_support, low_support`. `direction` is `treatment_up`, `control_up` or `none`. `low_support` flags rMATS events with too little read support (see Thresholds); it is `FALSE` for the other tools. A file that is missing is listed in the quality notes, and the later reports skip that tool and say so.

### Sign table

Every effect in the report is treatment minus control. The tools write their own signs; the report converts them in one place (`harmonise_effect` in Rmd 01):

| Tool | Value read | Conversion | Result |
|------|-----------|------------|--------|
| rMATS | `IncLevelDifference` (group 1 = treatment) | none | treatment minus control |
| SUPPA2 | `dPSI` (control minus treatment) | flipped | treatment minus control |
| DEXSeq exon usage | `log2fold_{CONTROL}_{TREATMENT}` | flipped | treatment minus control |
| DEXSeq DTU | `log2fold_{CONTROL}_{TREATMENT}` | flipped | treatment minus control |
| edgeR exon usage | `logFC` | none | treatment minus control |

How the signs were checked, on the output of the nf-core test-data run (test-data facts, not statements about other data):

- rMATS: `IncLevelDifference` equals mean(`IncLevel1`) minus mean(`IncLevel2`) on 202 SE rows, with group 1 the treatment group.
- SUPPA2: the sign of `dPSI` equals the sign of YRI minus GBR in the mean PSI on 582 of the 1096 tested local rows (the rest have a dPSI of 0) and on 860 isoform rows (of 1854 tested), with 0 opposite.
- DEXSeq exon usage and DTU: `log2fold_YRI_GBR` matches YRI minus GBR on 7420 and 218 rows.
- edgeR exon usage agrees in sign with DEXSeq (after the DEXSeq flip) on 3 exon bins with identical coordinates.

Effect sizes are never compared across tools (dPSI and log2 fold change are not on one scale); the cross-tool report compares directions only.

### Thresholds and defaults

| Setting | Default | Notes |
|---------|---------|-------|
| FDR | 0.05 | Adjusted p-value cutoff for an event, gene or transcript |
| dPSI | 0.1 | Minimum absolute dPSI of an rMATS or SUPPA2 event; not applied to exon or transcript usage |
| rMATS count type | `JC` | `JCEC` is offered when the run has those files |
| Minimum read support (`min_reads`) | 10 | For each group, the mean over its samples of IJC plus SJC; the smaller of the two group means must reach it |

- Events below `min_reads` are flagged `low_support`, kept in the tables and excluded from the significance counts; they are never deleted, and the rMATS FDR is not recomputed. On very shallow data (like the nf-core test data) most events fall below 10, and the value may be lowered.
- SUPPA2 gives a p-value but no adjusted p-value. The report applies Benjamini-Hochberg per file over the tested rows (local events and isoforms separately) and says so.
- stageR adjusted p-values are valid at the alpha the run used. When the chosen FDR differs, the report re-runs stageR's stage-wise adjustment from the run's `stageRTx` object (when it is readable) and says which values it used.
- Enrichment is GO only (clusterProfiler; human and mouse). The image has no ReactomePA, KEGG.db or reactome.db, so KEGG/Reactome require an image rebuild; the wizard and the report both say so. The background of each enrichment is the set of genes testable in that tool, and a group needs at least 10 genes.

### Paired runs

rMATS writes the results of a paired run (`rmats_paired_stats: true`) to `star/rmats/{CONTRAST}_paired`. The report looks in `star/rmats/{CONTRAST}` and in `star/rmats/{CONTRAST}_paired`, and when both hold tables, takes the one that matches the run's pairing switch. The pairing of rMATS and SUPPA2 is printed in Step 6.

---

## Tests

```bash
bash nfcore-rnasplice-report/tests/run_all_tests.sh           # default: static suites, no R
bash nfcore-rnasplice-report/tests/run_all_tests.sh --full    # adds the R tests and the proof-tool self-test
```

- The default suite uses bash, awk and sed only (R is never started) and runs each test under `env -i` with an explicit PATH. It checks the skill text, this README and the root `README.md` row, and runs the mutation rows (each mutates a copy of the skill, the README or the runner and must turn the checker red). It takes about 4 minutes (191 to 231 s measured).
- `--full` also runs the R tests through sbatch on a compute node (the Rmd functions on fixtures from the nf-core test-data output, and a render of Rmds 01 to 07), the stub test of the R test runner, the proof-tool self-test and the R mutation rows. It takes about 40 minutes. It never runs R on a login node.
- The proof-tool self-test needs git 1.8.5 or newer; `/usr/bin/git` on the cluster is 1.8.3. Set `RNASPLICE_REPORT_TEST_GIT_DIR` to a directory with a newer git, which is put first on the PATH of every test:

```bash
RNASPLICE_REPORT_TEST_GIT_DIR=$HOME/.conda/envs/git-new/bin bash nfcore-rnasplice-report/tests/run_all_tests.sh --full
```

---

## Validation status

What has run:

- The static suite (`tests/run_all_tests.sh`, default run), including the mutation rows: ALL TESTS PASS.
- The R tests of the report's functions (four scripts: 101, 105, 147 and 100 checks, all pass) on the real output of the nf-core test-data run of nf-core/rnasplice (4 samples, chrX, 50,000 read pairs), and the rendering of Rmds 01 to 07 on that output.
- The wizard end to end and the generated SLURM scripts on the cluster (2026-10-06), on the output of the nf-core test-data run (4 samples, chrX, 50,000 read pairs), in four runs: all modules; rMATS only; BAM input; all modules with `min_reads` 0. The wizard completed in every run, and every Bash block was accepted by the Claude Code safety check (skill at commit 0abe782).
- The 24 report jobs of the four runs all completed, and the chain dependencies shown by `scontrol` were as designed: jobs 02 to 06 `afterok` on job 01 with cancellation when job 01 fails, job 07 `afterany` on every earlier job.
- With `min_reads` 0 the rMATS counts per type (SE 5, A5SS 3, A3SS 2, MXE 3, RI 1) equal an independent recomputation with awk from the raw tables; the pipeline output directories were not modified.
- Chain wall times were 23 s, 17 s, 18 s and 71 s (all modules, rMATS only, BAM input, `min_reads` 0), and the four reports take 26 MB of output in total.

What has not run:

- The live run above used the skill at commit 0abe782, in which the Rmd files were cut out with a separate awk command. The block that extracts them from the installed skill file, the lines of Step 2 that drop the mail lines when no mail is wanted, and the exact-marker check of the rendering block did not exist then: they have static and render tests only, and have not run on the cluster.
- A later version of the block that extracts the Rmd files was refused by the Claude Code safety check (a removal through a shell variable in a cleanup step), and was replaced by a block that validates the skill file first and then writes the files, with no cleanup. The live run of the final blocks is part of the real-data run, which has not completed.
- Real data: not yet run.
- Peak memory of the report jobs: not measured.
- The behaviour when a report job fails: not yet run.

The sign checks above come from the test-data run only.
