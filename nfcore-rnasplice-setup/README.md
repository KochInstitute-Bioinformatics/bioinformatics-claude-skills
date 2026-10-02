# `/nfcore-rnasplice-setup` — nf-core/rnasplice Differential Splicing Setup Skill

An interactive Claude Code skill that sets up and submits [nf-core/rnasplice](https://nf-co.re/rnasplice) for differential alternative splicing of bulk RNA-seq on an HPC cluster with SLURM and Singularity: rMATS (splicing events), SUPPA2 (event and isoform PSI), DEXSeq and edgeR differential exon usage, and DEXSeq differential transcript usage (DRIMSeq filter, DEXSeq, stageR). It follows the style of `/nfcore-rnavar-setup`: every pipeline parameter goes into one params file, and every analysis switch is written explicitly. It sets up the run; it does not analyse the results.

The skill is pinned to nf-core/rnasplice 1b447239488097651d8eac44bca2c1556865eb0f (`dev-1b44723`), an unreleased development commit of the pipeline, not a release (see Known limitations).

---

## Installation

```bash
cp nfcore-rnasplice-setup.md ~/.claude/commands/
```

Then invoke it in Claude Code with `/nfcore-rnasplice-setup`.

---

## Prerequisites

| Requirement | Notes |
|-------------|-------|
| SLURM scheduler | Scripts target the `bcc` partition |
| Singularity | Loaded with `module add singularity/3.10.4` |
| Nextflow 26.04.0 or newer | In a conda environment; verified with Nextflow 26.04.6 (env `nf-env`). The submission script checks the version on the compute node |
| Internet on compute nodes | The job downloads the pinned revision of the pipeline from GitHub and the containers (depot.galaxyproject.org and community-cr-prod.seqera.io) |
| Internet on the login node | The wizard uses `WebFetch` (GitHub API, Ensembl) |

---

## What the skill does

| Step | Topic | Asked or detected |
|------|-------|-------------------|
| 0 | Working directory | Detected (`pwd`); every file is written there |
| 1 | Email for SLURM notifications | Asked |
| 2 | Conda environment with Nextflow | Asked (the version is checked later, by the submission script on the compute node) |
| 3 | Pipeline version | Pinned: nf-core/rnasplice 1b447239488097651d8eac44bca2c1556865eb0f. The latest release is looked up (GitHub API); a release newer than 1.0.4 is offered only when its schema has every key the skill writes, and is not verified by this skill |
| 4 | Input and samplesheet rows | FASTQ (files found, paired-end detected, names sanitised to valid R identifiers, name collisions flagged) or genome BAM; strandedness asked (or read from the RSeQC results of an nf-core/rnaseq run of the same samples) |
| 5 | Conditions, contrasts, paired design | Asked; both sheets validated in a scratch directory before they are written |
| 6 | Read length for rMATS | Detected from the first 1000 reads of every FASTQ file (asked for BAM input); warns when read lengths differ, and stops by default when they differ between conditions |
| 7 | Organism and genome | Ensembl release in the shared genome folder (missing FASTA/GTF: a download helper) or a custom FASTA/GTF. An existing STAR index is reused only when compatible (`versionGenome 2.7.4a`); the Salmon index is always built by the pipeline, never reused |
| 8 | Analyses and their settings | Asked: rMATS (default), SUPPA2, DEXSeq exon usage, edgeR exon usage, DEXSeq DTU; every module switch is written |
| 9 | Trimming and QC | Pipeline defaults (nothing asked) |
| 10 | MultiQC title, output directory, `nextflow.config` | Asked; config written only if absent |
| 11 | Params file, submission script, download helper | Generated |
| 12 | Summary and hand-off | Where each result is, which question it answers, and the sign of each difference column |

### Key design points

- **Every analysis switch is explicit.** In the pipeline's own configuration every analysis module except LeafCutter is on, so a switch left out of the params file runs that module anyway. The params file always contains all of them, and every option key that has a non-null pipeline default, at that default (except the settings below).
- **MISO is not used** (`sashimi_plot: false`): in this pipeline MISO only draws sashimi plots for a short gene list (by default three human Ensembl IDs), it is no genome-wide splicing test, and MISO itself is unmaintained Python 2 software. For sashimi plots, use rmats2sashimiplot or ggsashimi on the BAM files afterwards. IsoformSwitchAnalyzeR and LeafCutter are written `false` and not offered.
- **Strandedness is asked**, never guessed (rnasplice has no `auto`), one value for all samples (rMATS needs it); the wizard can read the RSeQC `infer_experiment.txt` files of an nf-core/rnaseq run of the same samples, with nf-core/rnaseq's rule.
- **`rmats_read_len` is always set** from the reads (the pipeline default of 40 is wrong for almost all data).
- **Paired statistics only when confirmed.** At this revision `rmats_paired_stats` defaults to `false` and SUPPA2's `diffsplice_paired` to `true`. The skill writes both from the confirmed pairing: `false` unless the user confirms a paired design, which is asked only with exactly two conditions of equal size. In a paired design the rows are sorted by subject within each condition, so that the i-th sample of one condition is paired with the i-th sample of the other (rMATS pairs them by samplesheet position, verified for this revision); the validator checks it with a pairs file.
- **DTU filter values are explicit:** the pipeline defaults fit one design size only, so the skill uses the rule of Love et al. 2018 (as `/bulk-rnaseq-pipeline` does): `min_samps_gene_expr` = the number of samples, `min_samps_feature_expr` and `min_samps_feature_prop` = the size of the smallest compared condition.
- **Genome BAM input:** strandedness and read type are written to the BAM samplesheet (`sample,condition,genome_bam,strandedness,single_end`); this revision passes them to rMATS, DEXSeq and featureCounts. BAMs from a splice-aware aligner (for example the STAR BAMs of an nf-core/rnaseq run), aligned to the same FASTA and GTF; no `.bai` needed. With BAM input only rMATS, DEXSeq exon usage and edgeR exon usage run (DTU and SUPPA2 need Salmon quantification from reads).
- **Salmon route:** `aligner: "star"` and `pseudo_aligner: "salmon"`: the STAR alignments feed rMATS, DEXSeq and edgeR, and Salmon on the reads feeds DTU and SUPPA2 once (the pipeline default quantifies twice). With FASTQ input, Salmon (index build and quantification) runs even in an rMATS-only run, because `pseudo_aligner` has no off value at this revision.
- **Indexes:** only an existing STAR index is reused, when its `versionGenome` is `2.7.4a` (a note is shown when its `sjdbOverhang` is not 100, the pipeline's own value; at this revision the pipeline copies a given index into its `work/` directory, a path this skill's verification did not exercise: unverified). The Salmon index is always built by the pipeline from the transcripts it extracts from the GTF: an index from `/nfcore-rnaseq-setup` is built from Ensembl cDNA and lacks the non-coding transcripts. A decoy-aware Salmon index of a human or mouse genome likely takes more than an hour to build (unverified). Built indexes are not kept (`save_reference: false`).

### Resources

The `nextflow.config` that the skill writes (only if none exists) sets a process default of 2 CPUs, 8 GB and 4 h, raises the processes below with `withName` selectors, and caps every task at 16 CPUs, 64 GB and 24 h (`resourceLimits`). Every other process keeps the resources of the pipeline's own process labels, which take precedence over that default: in the verification run they requested up to 12 CPUs, 64 GB and 16 h (the SUPPA2 diffSplice and cluster steps and `MAKE_TRANSCRIPTS_FASTA`).

| Selector | CPUs | Memory | Time |
|----------|------|--------|------|
| `STAR_GENOMEGENERATE` | 8 | 64 GB | 8h |
| `STAR_ALIGN` | 8 | 48 GB | 8h |
| `RMATS_PREP` | 4 | 16 GB | 8h |
| `RMATS_POST` | 8 | 32 GB | 16h |
| `DEXSEQ_COUNT` | 2 | 8 GB | 8h |
| `DEXSEQ_EXON` | 8 | 32 GB | 8h |
| `DEXSEQ_DTU` | 8 | 32 GB | 8h |
| `SALMON_QUANT.*` | 8 | 36 GB | 8h |

These values are judgement and, for `SALMON_QUANT`, the memory and time of the pipeline's own tested process label; they are NOT measured on a real genome. The first seven selectors were applied in the verification run on the nf-core test data; the `SALMON_QUANT.*` selector was applied in the cluster acceptance run on the same data: its SALMON_QUANT_SALMON tasks requested 8 CPUs, 36 GB and 8 h (in the verification run, whose config's `.*:SALMON_QUANT` matched no process, Salmon ran with its label's 6 CPUs, 36 GB and 8 h). With this skill's route no SALMON_QUANT_STAR task runs; the pattern covers both names. Compare them with `pipeline_info/execution_trace.txt` after the first real run. The head job requests 2 CPUs, 8 GB and 48 h (`-n 2 --mem=8G -t 48:00:00`): it only coordinates the pipeline but must outlive every task, so it asks for more than the 4 h that this repository's skills use by default (a deliberate exception). The genome download helper requests 2 CPUs, 8 GB and 4 h.

---

## Output files

All files are written in the working directory; `{SHEET_PREFIX}` is the file prefix chosen in Step 5 (for example `{SEQ_DATE}_{WD_NAME}`: the sequencing date and the name of the working directory).

| File | Description |
|------|-------------|
| `{SHEET_PREFIX}_samplesheet.csv` | Samplesheet (`sample,fastq_1,fastq_2,strandedness,condition`, or `sample,condition,genome_bam,strandedness,single_end` for genome BAM) |
| `{SHEET_PREFIX}_contrasts.csv` | Contrasts (`contrast,treatment,control`) |
| `{SHEET_PREFIX}_pairs.csv` | Paired design only: `sample,subject`, a record of the pairing (the pipeline does not read it) |
| `{SHEET_PREFIX}_params.yaml` | Every pipeline parameter of the run, passed with `-params-file` |
| `nextflow.config` | SLURM + Singularity profiles with process selectors (written only if absent) |
| `nf-core_rnasplice_dev-1b44723.sh` | Pipeline submission script (changes to the project directory, checks its files and the Nextflow version, runs the pipeline) |
| `download_genome_{REF_TAG}.sh` | Ensembl FASTA/GTF download (only if missing; URLs verified in the session) |

## Where the results are

Relative to the output directory, as recorded by the verification run (`{CONTRAST}` = the contrast name, `{TREATMENT}` and `{CONTROL}` = its two conditions); the hand-off note of the skill lists only the analyses that ran:

| Path | Content |
|------|---------|
| `star/rmats/{CONTRAST}/` | rMATS, one directory per contrast: for each event type (SE, A5SS, A3SS, MXE, RI) a .MATS.JC.txt and a .MATS.JCEC.txt table |
| `star/rmats/{CONTRAST}/SE.MATS.JC.txt` | example: skipped exons, junction reads only |
| `star/dexseq_exon/results/` | DEXSeq differential exon usage: DEXSeqResults.{CONTRAST}.csv and perGeneQValue.{CONTRAST}.csv |
| `star/edger/` | edgeR differential exon usage: contrast_{CONTRAST}.usage.exon.csv, .usage.gene.csv and .usage.simes.csv |
| `salmon/dexseq_dtu/results/` | DEXSeq differential transcript usage with stageR: dexseq/DEXSeqResults.{TREATMENT}-{CONTROL}.tsv and stager/getAdjustedPValues.{TREATMENT}-{CONTROL}.tsv |
| `salmon/suppa/` | SUPPA2 event and isoform PSI; differential splicing in diffsplice/per_local_event/local_{TREATMENT}-{CONTROL}.dpsi and diffsplice/per_isoform/transcript_{TREATMENT}-{CONTROL}.dpsi |
| `multiqc/` | MultiQC report: {MULTIQC_TITLE}_multiqc_report.html |

With a paired design the rMATS directory is `star/rmats/{CONTRAST}_paired/`. DTU and SUPPA2 files are named `{TREATMENT}-{CONTROL}`, not after the contrast. `pipeline_info/` holds the execution report, timeline and trace. `star/` also keeps a sorted copy of every BAM file (disk use).

**Direction of the difference columns** (as the skill's Step 12 states them):

| Output | Column | Positive value means | Source |
|--------|--------|----------------------|--------|
| rMATS | `IncLevelDifference` = mean(`IncLevel1`) − mean(`IncLevel2`); `b1` = the treatment samples | more inclusion in the treatment | verification run (bamlists and values) |
| SUPPA2 | dPSI = mean PSI of the control − mean PSI of the treatment, although the header reads `local_{TREATMENT}-local_{CONTROL}_dPSI` | more inclusion in the control | verification run (every event, values) and pinned code |
| DEXSeq DTU | `log2fold_{CONTROL}_{TREATMENT}` | a larger share of the transcript in the control | verification run (header plus checked rows: all 217 non-NA rows) and pinned code |
| DEXSeq exon usage | `log2fold_{CONTROL}_{TREATMENT}` | more usage of the exon bin in the control | verification run (header plus checked rows: all 7420 non-NA rows) and pinned code |
| edgeR exon usage | `logFC` = treatment − control | more usage of the exon in the treatment | pinned code only, not checked numerically |

---

## Typical workflow

```bash
# Submit from the project directory. With a genome download helper, the pipeline waits for it:
jid=$(sbatch --parsable download_genome_{REF_TAG}.sh) && sbatch --dependency=afterok:$jid nf-core_rnasplice_dev-1b44723.sh
# Without a helper:
sbatch nf-core_rnasplice_dev-1b44723.sh
```

The job downloads the pinned pipeline revision the first time it runs it (also when `~/.bashrc` sets `NXF_OFFLINE=TRUE`, verified for this revision). If the head job reaches its 48 h limit, or after raising a selector in `nextflow.config`, add ` -resume` at the end of the `nextflow run` line and submit again.

---

## Validation status

- Static tests: `bash tests/run_all_tests.sh` runs, each under `env -i` with an explicit PATH: the checker (every parameter against the recorded schema of nf-core/rnasplice 1b447239488097651d8eac44bca2c1556865eb0f, module switches, defaults against the recorded pipeline configuration, process selectors against the processes of the verification runs, the launch line, guarded module loads, this README), the tests that run the skill's own code blocks (sample names, sheet validation, strandedness from RSeQC, BAM rule, read length, params rendering), the stub dry runs of the submission script and the download helper, a test of the test runner itself on a fake tree (`test_run_all.sh`: failures propagate, every listed test runs, the git version check), and the mutation proofs of every structural check. The default run takes 13 to 15 minutes (bash, awk and sed only) and skips the proof-tool self-test (`test_proof_tools.sh`, about 30 minutes), and says so; `bash tests/run_all_tests.sh --full` includes it. The self-test needs git 1.8.5 or newer (`/usr/bin/git` on the cluster is 1.8.3): set `RNASPLICE_TEST_GIT_DIR` to the directory of a newer git (for example `RNASPLICE_TEST_GIT_DIR=$HOME/.conda/envs/git-new/bin bash tests/run_all_tests.sh --full`).
- Verification gate (2026-10-01): nf-core/rnasplice 1b447239488097651d8eac44bca2c1556865eb0f with Nextflow 26.04.6 on the cluster, on nf-core's tiny test dataset: 4 paired-end human chrX samples (about 30 MB of FASTQ, 75 bp reads, two conditions). Runs, in `/net/bmc-lab3/data/bcc/yannvrb/rnasplice_gate/results/` (on the author's cluster account; these run directories are not distributed, and the recorded evidence is in `tests/fixtures/`):
  - `2026-10-01_g2a_dev`: the pipeline's own test profile;
  - `2026-10-01_g3_typed`: a fully typed params file with all modules (the reference run);
  - `2026-10-01_g4_rmats_reverse`: rMATS only, `reverse` strandedness;
  - `2026-10-01_g4o_order_paired`: rows reversed within each condition, paired rMATS statistics (bamlist order = samplesheet order);
  - `2026-10-01_g5_bam`: genome BAM with the `strandedness` and `single_end` columns;
  - `2026-10-01_g5b_bam_nocols`: genome BAM without them;
  - `2026-10-01_g3s_save_reference`: the reference run resumed with `save_reference: true`.

  Details: `tests/fixtures/gate_report.md`.

Cluster acceptance (nf-core test data): DONE (2026-10-02)

The wizard was followed step by step on the cluster with the same nf-core test data, and three runs were submitted with the generated scripts:

- Static gate: the default `tests/run_all_tests.sh` under `env -i` (git 2.49.0): ALL TESTS PASS (MUTATIONS PASS (524), RUN ALL PASS (24)), 898 s.
- Run (a), all five analyses: job 11380326, wall time 11 min 49 s, 80/80 tasks COMPLETED. Checked: the params file is identical to the output of the tested renderer (`tests/render_params.sh`); the launch line carries no `--` flag; the version banner reads `Nextflow 26.04.6 (verified with 26.04.6)`; every task COMPLETED; the eight selectors against the execution trace, including `SALMON_QUANT.*` (its 4 SALMON_QUANT_SALMON tasks requested 8 CPUs, 36 GB and 8 h); `--readLength 75` in every rMATS prep command; every Step 12 hand-off path exists; the rMATS sign on 202/202 rows of SE.MATS.JC.txt and the SUPPA2 sign on 1098/1098 local and 1854/1854 isoform events, as in the direction table above.
- Run (b), rMATS only: job 11380328, 8 min 11 s (from the Nextflow log), 51/51 tasks COMPLETED; zero DEXSeq, edgeR, DRIMSeq, stageR, SUPPA2 or MISO tasks (Salmon still ran: `pseudo_aligner` has no off value).
- Run (c), genome-BAM input (the BAM files of run (a)): job 11380459, wall time 6 min 16 s, 43/43 tasks COMPLETED; only rMATS, DEXSeq exon usage and edgeR exon usage were offered; no Salmon task and no STAR alignment or index task.
- Not exercised: real data; the Ensembl download helper; a paired design; strandedness other than unstranded, and single-end input; a human- or mouse-size genome; the time and memory of any process on real data.

- No real biological data was run. The resources in `nextflow.config` are judgement, not measured on a real genome, and nothing in this skill is validated on real data.

---

## Known limitations

- **Unreleased pin.** nf-core/rnasplice 1b447239488097651d8eac44bca2c1556865eb0f is a development commit (`dev-1b44723`), not a release. Release 1.0.4 does not launch under Nextflow 26.04.6 (config parse error on `def check_max`, seen at the verification gate), so it was not used. The pinned commit is the merge of upstream PR #291 (2026-09-24), which ported the rMATS subworkflow to the nf-core structure (a rewrite of the subworkflow; evidence in `tests/fixtures/gate_report.md`): young code, a risk. (Gate observation, for a later choice of pin: release 1.0.4 on Nextflow 24.04.4 completed the test profile, but one rMATS 4.1.2 task hung for 43 minutes until it was killed by hand (the retry then finished in 2 s); cause not verified.)
- **Bumping the pin** (a newer commit or a release) means running the verification gate again before relying on the skill: the test profile and the gate runs above on the cluster; record the schema (`rnasplice_schema.json`, `schema_input*.json`), the pipeline's `params {}` block (`config_params.txt`), the process names and requested resources (`trace_process_names.txt`, `trace_resources_g3.tsv`), the output tree (`output_tree.txt`), the task command lines (`command_lines.txt`: strandedness, read type, bamlist order) and the container URLs (`verified_urls.txt`); set every key of `gate_values.tsv` (the checker validates each against its allowed values) and update `gate_report.md`; check every sign on the values or the code, never from a header alone; then update the skill until `tests/run_all_tests.sh --full` passes.
- **Option keys first run in the cluster acceptance run:** `ignore_tx_version`, `miso_genes`, `miso_read_len`, `fig_height`, `fig_width`, `isoformswitchanalyzer_alpha` and `isoformswitchanalyzer_dIF` are in the params file (at the pipeline defaults) but were not in the gate's typed run; they were first used in the cluster acceptance run (run (a) below). `star_index` (reuse of an existing STAR index) was not run either.
- **Single-end input and `forward` strandedness:** single-end BAM input and `forward` strandedness were verified from the pipeline code only, not by a run (the gate ran paired-end samples only: `unstranded` and `reverse`, from FASTQ and from genome BAM); single-end FASTQ input was not run either.
- One contrast was run with the skill's settings; several contrasts only in the pipeline's own test profile.
- Strandedness is asked, not detected (except from the RSeQC results of an existing nf-core/rnaseq run); rMATS needs one strandedness and one read type for all samples.
- Paired designs need exactly two conditions of equal size.
- Inputs: FASTQ or genome BAM only. No Salmon-results or transcriptome-BAM input, and no iGenomes `genome` key (FASTA and GTF are always given). Organisms other than human and mouse need a custom FASTA/GTF.
- With FASTQ input, Salmon (index and quantification) always runs, also when neither DTU nor SUPPA2 is chosen.
- Trimming, QC and SUPPA2 clustering use the pipeline defaults; edgeR exon usage uses diffSpliceDGE glmQLFit glmQLFTest in this revision.
- MISO (dropped: see Key design points), IsoformSwitchAnalyzeR and LeafCutter are off and not offered.
- The skill does not analyse the results; there is no downstream report skill for the splicing tables yet (planned later).

---

## rnasplice or `/bulk-rnaseq-pipeline`

`/bulk-rnaseq-pipeline` does gene-level differential expression, GSEA and, optionally, differential transcript usage (DRIMSeq → DEXSeq → stageR) on the Salmon output of an nf-core/rnaseq run. rnasplice answers which exons and splice events change (rMATS, SUPPA2, DEXSeq/edgeR exon usage). Its DEXSeq DTU runs the same steps as the bulk skill's DTU module (on another Salmon quantification): run it in one place only.
