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
| Internet | Not needed for the report: the Rmd files make no network call and install nothing, and the HTML pages load nothing from the internet (MathJax is switched off in every Rmd header). Only pulling a missing image (a compute-node job) needs it |
| A finished nf-core/rnasplice run | At least one of rMATS, SUPPA2, DEXSeq exon, edgeR exon or DEXSeq DTU output with more than a header line |

---

## What the skill does

The wizard asks one question at a time and states the default with each:

| Step | Content |
|------|---------|
| 0 | Working directory, report directory `results/{TODAY_ISO}_{WD_NAME}_report`; an existing directory is never overwritten silently: a new one takes the first free numbered suffix (`_2`, `_3`, ...), and when the old one is reused the summary page marks the files an earlier run left |
| 1 | Email for SLURM notifications (default: none) |
| 2 | R image and its package check (compute node) |
| 3 | Locate the run and detect the modules: a module counts as found only when its main result file exists and has more than a header. A relative `outdir`, `input` or `contrasts` in the params file is resolved from the directory that holds the params file (the run directory), not from the directory the wizard runs in |
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

`run_all.sh` only calls sbatch: jobs 02 to 06 start when job 01 has finished successfully, job 07 when the others have ended, whatever their outcome. Start with **`report_summary.html`** in the report directory: it is the entry page, written by job 07. Tables go to `tables/`, figures to `figures/`, the harmonised tables to `rds/`. Job 01 writes a run stamp, `rds/run_stamp.rds` (time and parameters of the run), when it starts; the summary page lists a report, table or figure as written by this run only when its file is not older than the stamp, and lists older files (left by an earlier run in the same directory, or by a job of this run that failed or did not run) under "From an earlier run (not updated by this run)", with this run's thresholds at the top.

### How the scripts are made

Step 10 cuts the seven Rmd files out of the installed skill file into the report directory, writes two small template files to `templates/` there (`run_job.sbatch.tpl` and `run_all.tpl`), and renders the job scripts and `run_all.sh` from them with three short Bash blocks (checks and job list, settings, rendering), so that no block is long or nests a heredoc. The wizard refuses typed values (e-mail address, paths, project name) that contain a single quote, a double quote, a backslash or a line break; a `$`, a backtick or `$(...)` in a value is kept literally and no command in it is run. Typed numbers (memory, time, FDR, dPSI, minimum read support) are accepted only in their form, and the blocks that derive a value (the working directory, the skill file path, the paths of the params file, the absolute output directory) stop when it holds a quote, a backslash or a line break. The report directory is checked (under `results/`, not inside the pipeline output) in Step 3 before anything is written into it, and again in the blocks of Step 10.

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

How the signs were checked on real data (the human hiPSC-CM run of 2026-10-06, KO versus WT, treatment KO; see Validation status), with 0 mismatches of the conversion:

- rMATS: the run's rMATS step was started with the KO sample list as group 1, and `IncLevelDifference` equals mean(`IncLevel1`) minus mean(`IncLevel2`) on 1998 of 1998 SE rows checked; 10 significant SE events (among them RBPMS2 -0.89 and -0.854, FGFR2 -0.70 and +0.58, TFPI +0.654) have the report effect equal to the raw value.
- SUPPA2: on 10 events the raw `local_KO-local_WT_dPSI` has the sign of WT minus KO in the mean PSI per group (for example +0.970 with KO 0.030 and WT 1.000), and the report effect is its negative.
- DEXSeq exon usage: `log2fold_WT_KO` has the sign of the WT column minus the KO column on all 105501 significant bins with a fold change (0 opposite), and 5 bins were checked by hand.
- DTU: the same holds on the 6169 transcripts confirmed by stageR (stageR `gene` and `transcript` below 0.05) and on the 6654 transcripts with a DEXSeq adjusted p-value below 0.05 (0 opposite in both sets, recomputed with awk from the raw tables), and 5 transcripts were checked by hand.
- edgeR exon usage: on the 23930 exon bins whose coordinates are the same in edgeR and DEXSeq, edgeR's `logFC` has the sign of the flipped DEXSeq fold change on 22765 of the 22766 bins significant in both (1 disagrees).

`tests/fixtures/sign_rows.tsv` holds 34 of these real rows (at least 6 per tool, each with how its sign was checked), and the R tests check them against the trimmed real tables.

How the direction of a gene is chosen (cross-tool comparison, enrichment groups, gene tables of the usage report): the sign of the feature with the smallest adjusted p-value among the gene's significant features in that tool. On real data many adjusted p-values are exactly 0, so ties are common: among features with equal smallest adjusted p-value the one with the largest absolute effect is taken, and if that is tied too, the first by feature identifier. The result does not depend on the order of the rows. This tie-break is a convention, not a test: a gene whose tied features point in opposite directions (TEAD3 on the real data: rMATS SE -0.467 and MXE +0.4, both FDR 0) takes the direction of the larger change and is not marked separately.

Effect sizes are never compared across tools (dPSI and log2 fold change are not on one scale); the cross-tool report compares directions only.

### Thresholds and defaults

| Setting | Default | Notes |
|---------|---------|-------|
| FDR | 0.05 | Adjusted p-value cutoff for an event, gene or transcript |
| dPSI | 0.1 | Minimum absolute dPSI of an rMATS or SUPPA2 event; not applied to exon or transcript usage |
| rMATS count type | `JC` | `JCEC` is offered when the run has those files |
| Minimum read support (`min_reads`) | 10 | For each group, the mean over its samples of IJC plus SJC; the smaller of the two group means must reach it |

- Events below `min_reads` are flagged `low_support`, kept in the tables and excluded from the significance counts; they are never deleted, and the rMATS FDR is not recomputed. On very shallow data (like the nf-core test data) most events fall below 10, and the value may be lowered.
- DTU log2 fold changes beyond 10 in absolute value come from transcripts with (almost) no counts in one group: DEXSeq then reports values of 20 or more (on the real data up to about 22) that are not estimates of a change in usage. The report keeps the effect as it is (no capping) and marks it `zero-count artefact` in a column `effect_flag` of the DTU tables of the usage report and of `cross_tool.genes.csv` (for `effect_dtu`). The threshold of 10 is a heuristic, not a statistical test.
- The DEXSeq exon bins table lists the bins with an effect first: a bin of a group of overlapping genes (`GENE1+GENE2+GENE3`) can have no fold change, and such bins would otherwise fill the top of the table.
- A DEXSeq group of overlapping genes or a SUPPA2 locus `GENE1_and_GENE2` counts for each of its genes in the cross-tool comparison and the enrichment. The enrichment report prints both numbers for each tool: the identifiers as the tool writes them (for DEXSeq, the pipeline's own count of `perGeneQValue` groups) and the genes after the split (on the real data 9567 groups and 12867 genes); the cross-tool sets table has both columns.
- The minimum read support (`min_reads`) is a heuristic too: on the real data it changes the rMATS count from 7079 events (0) to 4497 (10).
- SUPPA2 gives a p-value but no adjusted p-value. The report applies Benjamini-Hochberg per file over the tested rows (local events and isoforms separately) and says so.
- stageR adjusted p-values are valid at the alpha the run used. When the chosen FDR differs, the report re-runs stageR's stage-wise adjustment from the run's `stageRTx` object (when it is readable) and says which values it used.
- Enrichment is GO only (clusterProfiler; human and mouse). The report never runs KEGG or Reactome: when they are requested it says that they are not offered by this report and, when the image lacks ReactomePA, KEGG.db or reactome.db (the shared image lacks all three), that KEGG/Reactome require an image rebuild. The background of each enrichment is the set of genes testable in that tool, and a group needs at least 10 genes.

### Paired runs

rMATS writes the results of a paired run (`rmats_paired_stats: true`) to `star/rmats/{CONTRAST}_paired`. The report looks in `star/rmats/{CONTRAST}` and in `star/rmats/{CONTRAST}_paired`, and when both hold tables, takes the one that matches the run's pairing switch. The pairing of rMATS and SUPPA2 is printed in Step 6.

---

## Tests

```bash
bash nfcore-rnasplice-report/tests/run_all_tests.sh           # default: static suites, no R
bash nfcore-rnasplice-report/tests/run_all_tests.sh --full    # adds the R tests and the proof-tool self-test
```

- The default suite uses bash, awk and sed only (R is never started) and runs each test under `env -i` with an explicit PATH. It checks the skill text, this README and the root `README.md` row, and runs the mutation rows (each mutates a copy of the skill, the README or the runner and must turn the checker, the runner test or the render test red). With the current mutation rows it took 17 min 14 s on a compute node (2026-10-06).
- `--full` also runs the R tests through sbatch on a compute node (the Rmd functions on fixtures from the nf-core test-data output and on trimmed tables of the real hiPSC-CM run, and renders of Rmds 01 to 07 on both), the stub test of the R test runner, the proof-tool self-test and the R mutation rows. It took about 40 minutes before the test on the trimmed real tables was added; that test adds about 2 minutes to every run of the R tests, and every R mutation row runs them again. It never runs R on a login node.
- The trimmed real tables are in `tests/fixtures/real_hipscm/` (2.8 MB; how they were cut is in `tests/fixtures/README.md`). In that test the enrichment report runs on the DEXSeq exon usage genes only, because each GO enrichment takes 6 to 25 s in the image.
- The proof-tool self-test needs git 1.8.5 or newer; `/usr/bin/git` on the cluster is 1.8.3. Set `RNASPLICE_REPORT_TEST_GIT_DIR` to a directory with a newer git, which is put first on the PATH of every test:

```bash
RNASPLICE_REPORT_TEST_GIT_DIR=$HOME/.conda/envs/git-new/bin bash nfcore-rnasplice-report/tests/run_all_tests.sh --full
```

---

## Validation status

What has run:

- The static suite (`tests/run_all_tests.sh`, default run), including the mutation rows: ALL TESTS PASS.
- The R tests of the report's functions (five scripts: 135, 106, 153, 118 and 93 checks, all pass) on the real output of the nf-core test-data run of nf-core/rnasplice (4 samples, chrX, 50,000 read pairs) and on trimmed tables of the real hiPSC-CM run below, and the rendering of Rmds 01 to 07 on both.
- The wizard end to end and the generated SLURM scripts on the cluster (2026-10-06), on the output of the nf-core test-data run (4 samples, chrX, 50,000 read pairs), in four runs: all modules; rMATS only; BAM input; all modules with `min_reads` 0. The wizard completed in every run, and every Bash block was accepted by the Claude Code safety check (skill at commit 0abe782).
- The 24 report jobs of the four runs all completed, and the chain dependencies shown by `scontrol` were as designed: jobs 02 to 06 `afterok` on job 01 with cancellation when job 01 fails, job 07 `afterany` on every earlier job.
- With `min_reads` 0 the rMATS counts per type (SE 5, A5SS 3, A3SS 2, MXE 3, RI 1) equal an independent recomputation with awk from the raw tables; the pipeline output directories were not modified.
- Chain wall times were 23 s, 17 s, 18 s and 71 s (all modules, rMATS only, BAM input, `min_reads` 0), and the four reports take 26 MB of output in total.

### Real data (2026-10-06)

The skill (at commit 6e0350b) ran on the output of a real nf-core/rnasplice run: human hiPSC-derived cardiomyocytes, KO versus WT (treatment KO, control WT), 3 versus 3 samples, reverse-stranded (the Salmon strandedness helper of `/nfcore-rnasplice-setup`: ISR, reverse fraction 0.9999), nf-core/rnasplice at the development revision 1b447239488097651d8eac44bca2c1556865eb0f, STAR and Salmon, all five modules (rMATS, SUPPA2, DEXSeq exon usage, edgeR exon usage, DEXSeq DTU with stageR). The pipeline's execution trace shows 103 tasks completed; one attempt of the DEXSeq exon-usage task was killed (exit 143) and its retry completed. The report ran twice, with `min_reads` 10 (the default) and 0, from a working directory separate from the run directory. The second render (`min_reads` 0) did not run the checks block of Step 10 and the package check again (same image and modules): its job list was written with the same seven lines and its two template files were copied from the first render (identical); every other block ran as written.

- The wizard completed; every Bash block was accepted by the Claude Code safety check, including the block that extracts the Rmd files (validate first, then write), which wrote the seven Rmd files identical to the blocks of the skill. Step 3 found the five modules and Step 4 the columns `log2fold_WT_KO` and `local_KO-local_WT_dPSI` in the raw files.
- Counts at the defaults (FDR 0.05, |dPSI| 0.1, `min_reads` 10, rMATS JC), each equal to an independent recomputation with awk from the raw tables: rMATS SE 2049, A5SS 257, A3SS 429, MXE 1394, RI 368 (4497 events, 2453 genes); SUPPA2 536 local events and 576 isoforms; DEXSeq exon usage 105669 bins and 9567 `perGeneQValue` groups (12867 genes after splitting the groups of overlapping genes); edgeR exon usage 5944 genes (Simes) of 16804, 32253 exons, 8818 genes by the F-test; DTU 3628 genes (stageR) of 11452 and 6169 transcripts (stageR `transcript` below 0.05; their `gene` is below 0.05 too), 6654 transcripts by the DEXSeq adjusted p-value alone (two different sets).
- With `min_reads` 0 the rMATS counts are SE 3206, A5SS 537, A3SS 849, MXE 1867, RI 620 (7079 events), also equal to the recomputation; with 10 they are the 4497 above.
- Signs: see "How the signs were checked on real data" above (0 mismatches; edgeR agrees with DEXSeq on 22765 of 22766 shared significant bins). RBPMS2 (the KO gene's skipped exon chr15:64,750,343-64,750,381): rMATS SE -0.854 (FDR 0), SUPPA2 raw +0.454 reported as -0.454, edgeR -1.588 on the same exon, DEXSeq bins more used in WT: every tool says more included or more used in WT, as in the real-data record of `/nfcore-rnasplice-setup`.
- The published analysis of this dataset reports 2,679 events; this report finds 4497 rMATS events with `min_reads` 10 (1.68 times) and 7079 with `min_reads` 0 (2.64 times). An earlier run of the same data found 7,031 events. The reason for the difference to the published number has not been determined.
- Report jobs (32 GB, 4 hours, 4 CPUs each): 01 up to 2 min 21 s, 02 about 15 s, 03, 04 and 05 about 37 s each, 06 (enrichment) 15 min 17 s, 07 4 s; the whole chain 17.7 min. No limit was reached (no out-of-memory or kill message in any log), and each report directory takes 79 MB.
- The enrichment ran for all 15 groups (five tools, three groups each); the terms found are recorded as observations only: translation and ribosome, RNA splicing and focal adhesion dominate, with some muscle terms.

The tie-break text, the DTU zero-count flag, the DEXSeq bins order, Rmd 06's two gene counts (with the identifiers before the split in the cross-tool sets table) and the Step 3 outdir rule are tested on the trimmed real tables (`tests/fixtures/real_hipscm/`) and by the static and render tests; they have not been run on the whole real output. The blocks of the skill at commit 6e0350b ran in the real-data run; the Step 3 block that resolves a relative `outdir` changed after that run (commit ca78ad5) and has static and render tests only.

What has not run:

- The live run on the test data used the skill at commit 0abe782, in which the Rmd files were cut out with a separate awk command; the real-data run used the extraction block and the rendering block with its exact-marker check. The lines of Step 2 that drop the mail lines when no mail is wanted have static and render tests only; the real-data record does not say whether they ran.
- The block that extracts the Rmd files validates the skill file first and then writes the files; it deletes nothing. That block ran in the real-data run.
- Not run, neither live nor on real data: several contrasts in one run (Rmds 01 to 07 were rendered with two contrasts only on a copy of the test fixtures), paired rMATS statistics (`rmats_paired_stats: true`, the `{CONTRAST}_paired` directory; only the choice of the directory has R tests), a mouse run (Rmd 06 has an R test with mouse identifiers only), the rMATS JCEC tables, the image pull job (`image_check/pull_image.sbatch` was never submitted), and a report directory with a numbered suffix (`_2`).
- Changed after the real-data run, with static, render and R tests only: the run stamp and the section From an earlier run of the summary page, the first free numbered suffix of Step 0, the report directory check of Step 3 and its repetition in the extraction block, the checks of typed numbers and derived values, the check that the opening line of an Rmd block is the line right after its heading, the stageR identifier rule (identifiers cut at the first dot on both sides of a join, with a WARNING for unmatched transcripts), the KEGG/Reactome wording of Rmd 06, and MathJax switched off in the HTML pages.
- Peak memory of the report jobs: not measured. On the real data the job accounting database was unavailable, so peak memory and the exit codes of jobs 01 to 05 were not read. Their outcome was read from the output files and the chain instead: the `afterok` dependencies released jobs 02 to 06, every log ends with the written HTML file and no error, and jobs 06 and 07 completed with exit code 0 (read with `scontrol`).
- The figures of the real-data reports were not viewed; only the tables and the text of the HTML pages were checked.
- The behaviour when a report job fails: not yet run.

The sign checks above come from the test-data run and the real-data run.
