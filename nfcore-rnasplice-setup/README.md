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
| bash 4 or newer | The wizard's helper functions use associative arrays (the cluster has bash 4.2.46) |
| Disk space | About 70 GB plus 18 GB per sample for `work/` and the results together (6 samples: about 175 GB; FASTQ files not counted), estimated from the disk use of the real-data test on one human dataset of about 38 M read pairs of 78 bp per sample, with the STAR_ALIGN override of the `nextflow.config` template (about 38 GB or more per sample without it); it grows with the reads. The fixed part is mostly the STAR index copy (28.5 GB), the Salmon index (about 22 GB), the transcript FASTA (about 7.4 GB) and the gene-filter step (about 4.4 GB). Nextflow's work directory is `work/` in the project directory: resuming a run needs it, and after a successful run it can be deleted; the skill never deletes it |

---

## What the skill does

| Step | Topic | Asked or detected |
|------|-------|-------------------|
| 0 | Working directory | Detected (`pwd`); every file is written there |
| 1 | Email for SLURM notifications | Asked |
| 2 | Conda environment with Nextflow | Asked (the version is checked later, by the submission script on the compute node) |
| 3 | Pipeline version | Pinned: nf-core/rnasplice 1b447239488097651d8eac44bca2c1556865eb0f. The latest release is looked up (GitHub API); a release newer than 1.0.4 is offered only when its schema has every key the skill writes, and is not verified by this skill |
| 4 | Input and samplesheet rows | FASTQ (files found, paired-end detected, names sanitised to valid R identifiers, name collisions flagged) or genome BAM; strandedness asked (or read from the RSeQC results of an nf-core/rnaseq run of the same samples, or, for FASTQ input, from a Salmon helper that detects the library type from up to 1,000,000 reads of a few samples) |
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
- **Strandedness is asked**, never guessed (rnasplice has no `auto`), one value for all samples (rMATS needs it); the wizard can read the RSeQC `infer_experiment.txt` files of an nf-core/rnaseq run of the same samples, with nf-core/rnaseq's rule, or (FASTQ input) generate an optional sbatch helper, `infer_strandedness_salmon_{WD_NAME}.sh`, which runs Salmon's library-type detection (`-l A`) on at most the first 1,000,000 reads (read pairs) of a few samples against a Salmon index of the organism (found in the shared genome folder, or named by the user), and reads the result with the same thresholds (the pipeline never receives that index). A kit name alone does not determine the direction: in the real-data test the GEO metadata said "NEB Ultra II protocol" (without "Directional"), yet the reads were reverse-stranded (ISR); the hand-written helper found it in about 6 minutes; the generated helper took 11 s on two samples of the 50,000-pair nf-core test data and found them unstranded (ISF about equal to ISR, IU 0).
- **`rmats_read_len` is always set** from the reads (the pipeline default of 40 is wrong for almost all data).
- **Paired statistics only when confirmed.** At this revision `rmats_paired_stats` defaults to `false` and SUPPA2's `diffsplice_paired` to `true`. The skill writes both from the confirmed pairing: `false` unless the user confirms a paired design, which is asked only with exactly two conditions of equal size. In a paired design the rows are sorted by subject within each condition, so that the i-th sample of one condition is paired with the i-th sample of the other (rMATS pairs them by samplesheet position, verified for this revision); the validator checks it with a pairs file.
- **DTU filter values are explicit:** the pipeline defaults fit one design size only, so the skill uses the rule of Love et al. 2018 (as `/bulk-rnaseq-pipeline` does): `min_samps_gene_expr` = the number of samples, `min_samps_feature_expr` and `min_samps_feature_prop` = the size of the smallest compared condition.
- **Genome BAM input:** strandedness and read type are written to the BAM samplesheet (`sample,condition,genome_bam,strandedness,single_end`); this revision passes them to rMATS, DEXSeq and featureCounts. BAMs from a splice-aware aligner (for example the STAR BAMs of an nf-core/rnaseq run), aligned to the same FASTA and GTF; no `.bai` needed. With BAM input only rMATS, DEXSeq exon usage and edgeR exon usage run (DTU and SUPPA2 need Salmon quantification from reads).
- **Salmon route:** `aligner: "star"` and `pseudo_aligner: "salmon"`: the STAR alignments feed rMATS, DEXSeq and edgeR, and Salmon on the reads feeds DTU and SUPPA2 once (the pipeline default quantifies twice). With FASTQ input, Salmon (index build and quantification) runs even in an rMATS-only run, because `pseudo_aligner` has no off value at this revision.
- **Indexes:** only an existing STAR index is reused, when its `versionGenome` is `2.7.4a` (a note is shown when its `sjdbOverhang` is not 100, the pipeline's own value; at this revision the pipeline copies a given index into its `work/` directory: in the real-data test, 28.5 GB in about 1 minute, with 88.5-92.7% of the reads of each sample mapped uniquely). The Salmon index is always built by the pipeline from the transcripts it extracts from the GTF: an index from `/nfcore-rnaseq-setup` is built from Ensembl cDNA and lacks the non-coding transcripts. The pipeline's decoy-aware Salmon index of the human genome took 33 min (peak 19.8 GB) in the real-data test; a mouse genome was not measured, and STAR_GENOMEGENERATE was not measured at all (the real-data test reused a STAR index). Built indexes are not kept (`save_reference: false`).
- **No STAR transcriptome BAM.** The pinned `conf/modules.config` always makes STAR write a transcriptome BAM (an estimated 20 GB or more per sample of 38 M read pairs, projected from a partial file in the real-data test) that this skill's route (`aligner: "star"`) never uses. The generated `nextflow.config` therefore repeats the pinned STAR_ALIGN `ext.args` without `--quantMode TranscriptomeSAM` and `--quantTranscriptomeSAMoutput BanSingleEnd`; the launch line is unchanged. Verified on the real-data test with the same selector and arguments in a separate config file (resumed run: 22 tasks cached, 0 failed); the form written into `nextflow.config` was applied in a run on the nf-core test data (2026-10-03: no STAR_ALIGN task had `--quantMode`, no transcriptome BAM was written; see Validation status). Not verified: an unmodified STAR_ALIGN run to completion at full size. BAM input is unaffected (no STAR_ALIGN task).

### Resources

The `nextflow.config` that the skill writes (only if none exists) sets a process default of 2 CPUs, 8 GB and 4 h, raises the processes below with `withName` selectors, and caps every task at 16 CPUs, 64 GB and 24 h (`resourceLimits`). Every other process keeps the resources of the pipeline's own process labels, which take precedence over that default: in the verification run they requested up to 12 CPUs, 64 GB and 16 h (the SUPPA2 cluster steps; the SUPPA2 diffSplice steps and `MAKE_TRANSCRIPTS_FASTA` had the same label and now have selectors).

| Selector | CPUs | Memory | Time | Evidence |
|----------|------|--------|------|----------|
| `STAR_GENOMEGENERATE` | 8 | 64 GB | 8h | verification gate (nf-core test data); not run in the real-data test (index reused) |
| `STAR_ALIGN` | 8 | 48 GB | 4h | real-data test: peak 37.7 GB, at most 15 min |
| `RMATS_PREP` | 2 | 8 GB | 4h | real-data test: single-threaded, 1.3 GB, 8 min |
| `RMATS_POST` | 8 | 16 GB | 8h | real-data test: 2.0 GB, 3 min (one contrast) |
| `DEXSEQ_COUNT` | 1 | 4 GB | 8h | real-data test: single-threaded, 1.15 GB, 55 min |
| `DEXSEQ_EXON` | 8 | 32 GB | 8h | verification gate; real-data test: 83 min, true memory peak not measured |
| `DEXSEQ_DTU` | 8 | 16 GB | 4h | real-data test: 11.7 GB, 2 min |
| `SALMON_QUANT.*` | 8 | 32 GB | 4h | real-data test: 20.2 GB, 13.5 min |
| `DIFFSPLICE_IO[EI]` | 1 | 8 GB | 16h | real-data test: single-threaded, 3.2-4.2 GB, up to 3 h 03 min |
| `MAKE_TRANSCRIPTS_FASTA` | 2 | 8 GB | 2h | real-data test: 4.2 GB, 2 min |

The values marked "real-data test" were set from what each process used in the real-data test (below): ONE human dataset, 6 samples of about 38 M read pairs of 78 bp, one contrast, all five analyses. They rest on that one dataset (the test ran with the earlier values); they were then applied in a run on the nf-core test data (2026-10-03), which shows that they are applied, not that they are enough for real data. `STAR_GENOMEGENERATE` and `DEXSEQ_EXON` keep the values of the verification gate. Without a selector, SUPPA2's `DIFFSPLICE_IOE`/`DIFFSPLICE_IOI` and `MAKE_TRANSCRIPTS_FASTA` would request their process label's 12 CPUs, 64 GB and 16 h; the SUPPA2 steps run single-threaded (about 3 h each in the real-data test, the longest tasks of the run). With this skill's route no SALMON_QUANT_STAR task runs; the pattern covers both names. The executor's `queueSize` is 20 (it was 10: in the real-data test, RMATS_POST and two Salmon tasks waited 11 to 13 min for a free slot while six 55-min DEXSEQ_COUNT tasks held six of the ten). The trace records `attempt` and `%cpu` too. Compare the selectors with `pipeline_info/execution_trace.txt` after the first real run; for `DEXSEQ_EXON`, the trace's `peak_rss` is misleading (187 GB reported in the real-data test for a 32 GB task that was not killed; `sstat` showed 8.4 GB; probably forked workers counted separately): compare with SLURM's MaxRSS instead. The head job requests 2 CPUs, 8 GB and 48 h (`-n 2 --mem=8G -t 48:00:00`): it only coordinates the pipeline but must outlive every task, so it asks for more than the 4 h that this repository's skills use by default (a deliberate exception). The genome download helper requests 2 CPUs, 8 GB and 4 h.

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
| `infer_strandedness_salmon_{WD_NAME}.sh` | Optional (FASTQ input): Salmon library-type detection on a subsample of the reads, on a compute node (Step 4) |

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
| edgeR exon usage | `logFC` = treatment − control | more usage of the exon in the treatment | pinned code; checked on one exon of the real-data test (logFC −1.59 on the RBPMS2 exon skipped in the KO) |

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

- Static tests: `bash tests/run_all_tests.sh` runs, each under `env -i` with an explicit PATH: the checker (every parameter against the recorded schema of nf-core/rnasplice 1b447239488097651d8eac44bca2c1556865eb0f, module switches, defaults against the recorded pipeline configuration, process selectors against the processes of the verification runs, the launch line, guarded module loads, this README), the tests that run the skill's own code blocks (sample names, sheet validation, strandedness from RSeQC and from Salmon, including the real Salmon output of the real-data test, BAM rule, read length, params rendering), the stub dry runs of the submission script, the download helper and the Salmon strandedness helper, a test of the test runner itself on a fake tree (`test_run_all.sh`: failures propagate, every listed test runs, the git version check), and the mutation proofs of every structural check. The default run takes about 50 minutes (bash, awk and sed only; measured on 2026-10-02, it varies with the load of the login node) and skips the proof-tool self-test (`test_proof_tools.sh`, about 100 minutes), and says so; `bash tests/run_all_tests.sh --full` includes it (about 2 h 30 min in all). The self-test needs git 1.8.5 or newer (`/usr/bin/git` on the cluster is 1.8.3): set `RNASPLICE_TEST_GIT_DIR` to the directory of a newer git (for example `RNASPLICE_TEST_GIT_DIR=$HOME/.conda/envs/git-new/bin bash tests/run_all_tests.sh --full`).
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

The wizard was followed step by step on the cluster with the same nf-core test data, and three runs were submitted with the generated scripts. The runs and the static gate below were made on the skill text of commit bcb63e9. The later changes (commit dd69f63 and the final review fixes) touched only the Step 4 helper functions (R2 mates, BAM sample names, BAM path style), the Step 5 file-name rules (one option when the dates are equal, the custom-name characters, the file-name check of `install_rnasplice_sheets`), wording, the tests and the checker; they were tested statically only. The generated files of the three runs match the skill as merged (commit 77bfb5d): in `a/`, `b/` and `c/`, `nextflow.config` is byte-identical to the Step 10 template of that commit, and `nf-core_rnasplice_dev-1b44723.sh` and `261002_{a,b,c}_params.yaml` match the Step 11 templates (with the Step 8 block) line for line, with only the placeholders filled; the module block, the params template, the submission script, the download helper and the read-length block are byte-identical to those of bcb63e9, and the sheet-validation block differs only in that file-name check. The `nextflow.config` template was changed after the real-data test (see **Real-data test** below); these three runs did not use its new form.

- Static gate (commit bcb63e9): the default `tests/run_all_tests.sh` under `env -i` (git 2.49.0): ALL TESTS PASS (MUTATIONS PASS (524), RUN ALL PASS (24)), 898 s.
- Run (a), all five analyses: job 11380326, wall time 11 min 49 s, 80/80 tasks COMPLETED. Checked: the params file is identical to the output of the tested renderer (`tests/render_params.sh`); the launch line carries no `--` flag; the version banner reads `Nextflow 26.04.6 (verified with 26.04.6)`; every task COMPLETED; the eight selectors against the execution trace, including `SALMON_QUANT.*` (its 4 SALMON_QUANT_SALMON tasks requested 8 CPUs, 36 GB and 8 h); `--readLength 75` in every rMATS prep command; every Step 12 hand-off path exists; the rMATS sign on 202/202 rows of SE.MATS.JC.txt and the SUPPA2 sign on 1098/1098 local and 1854/1854 isoform events, as in the direction table above.
- Run (b), rMATS only: job 11380328, 8 min 11 s (from the Nextflow log), 51/51 tasks COMPLETED; zero DEXSeq, edgeR, DRIMSeq, stageR, SUPPA2 or MISO tasks (Salmon still ran: `pseudo_aligner` has no off value).
- Run (c), genome-BAM input (the BAM files of run (a)): job 11380459, wall time 6 min 16 s, 43/43 tasks COMPLETED; its automatic sample names still ended in `_sorted` (made before Step 4b strips `_sorted.bam`; that fix is unit-tested, run (c) was not repeated); only rMATS, DEXSeq exon usage and edgeR exon usage were offered; no Salmon task and no STAR alignment or index task.
- Not exercised by these runs (see **Real-data test** below for the real-data run): real data; the Ensembl download helper; a paired design; strandedness other than unstranded, and single-end input; a human- or mouse-size genome; the time and memory of any process on real data.

Live re-run of the changed template (nf-core test data): DONE (2026-10-03)

After the changes that followed the real-data test, the wizard was followed again on the same nf-core test data with the changed skill (commit 30ee087) and Nextflow 26.04.6.
- `nextflow.config` was byte-identical to the Step 10 template and was accepted by Nextflow 26.04.6; the launch line was unchanged.
- Run (a), all five analyses: 80/80 tasks COMPLETED; run (b), rMATS only: 51/51; run (c), genome-BAM input from the BAM files of run (a): 43/43; the Salmon strandedness helper: 11 s for two samples of 50,000 read pairs, both found unstranded (the known library type).
- Every selector was applied: each task requested exactly the values of the Resources table; the trace has the `attempt` and `%cpu` columns; no task was retried.
- STAR_ALIGN: no task had `--quantMode`, and no transcriptome BAM was written. The rMATS JC and JCEC row counts and IncLevelDifference values were identical to the earlier acceptance run (other tables differ slightly, because Salmon's quantification is not bit-identical between runs).
- This was applied in a run on the test data: it shows that the template and its values are applied, not that they are enough for real data. Not exercised: `queueSize` 20 at scale (at most 6 tasks ran at once), the stop-a-run guidance, the disk estimate at real size.

Real-data test (one human dataset): DONE (2026-10-02), with the deviations listed below

The wizard was followed step by step on one real dataset with the skill as merged (commit 77bfb5d), the generated scripts were submitted, and the results were compared with the paper. This is one dataset, not a validation on real data in general.

- Dataset: Akerberg et al. 2022 (GEO GSE207681), RBPMS2-null versus wild-type human iPSC-derived cardiomyocytes, a 3 vs 3 subset: WT_1, WT_2, WT_3 = SRR20021261, SRR20021260, SRR20021259 (GSM6307608-GSM6307610) and KO_1, KO_2, KO_3 = SRR20021257, SRR20021256, SRR20021255 (GSM6307612-GSM6307614). Paired-end, 78 bp reads (the metadata said 75 bp; Step 6 measured 78), about 38 M read pairs per sample, reverse-stranded (found with a Salmon `-l A` subsample, ISR); Ensembl 116 genome, with the existing STAR index reused.
- What ran: all five analyses (rMATS, SUPPA2, DEXSeq and edgeR exon usage, DEXSeq DTU), one contrast (`KO_vs_WT`), no paired design.
- Completion: "Pipeline completed successfully", 81 tasks completed and 22 cached, 0 failed, 0 retried; the resumed run took 4 h 11 min (4 h 34 min from the first submission).
- Biology: RBPMS2 (ENSG00000166831) expression (Salmon gene TPM) was 67% lower in KO (mean 37.5 vs 115.4). The RBPMS2 exon chr15:64,750,343-64,750,381 skipped in the KO was found by all four splicing analyses (rMATS SE IncLevelDifference −0.854, FDR 0; DEXSeq exon bins padj 1.8e-56; edgeR exon logFC −1.59, FDR 1.7e-12; SUPPA2 SE dPSI +0.454). On the 1,000 most variable genes, PC1 (65.9% of the variance) separated WT from KO, and hierarchical clustering gave exactly the two groups; on rMATS SE PSI the separation was partial (one KO sample was an outlier).
- Comparison with the paper (the paper's cut-offs: FDR < 0.05, |IncLevelDifference| > 0.1, 0 uncalled replicates): rMATS found 7,031 significant events (junction-count tables), about 2.6 times the 2,679 events the paper reports. The reason was not determined (possible causes: 3 vs 3 samples here against 4 WT and 5 KO in the paper, Ensembl 116 against 101, rMATS version and JC or JCEC choice of the paper not stated in its accessible text). As in the paper, skipped exons were the largest class (45% here, 52% in the paper) and exon exclusion in the KO was more frequent than inclusion (1,816 vs 1,372 skipped exons).
- Directions: every sign statement of Step 12 held on real rows: rMATS on the RBPMS2 row (with the KO samples as `--b1`), SUPPA2, DTU and DEXSeq exon usage on every row, edgeR (logFC = KO − WT) on one exon; every hand-off path existed.
- Deviations of the test: the first head job was cancelled after 22 min, because the STAR transcriptome BAM files (an estimated 20 GB or more per sample) were filling the disk budget, and seven of its task jobs, still running, were cancelled by hand. The run was resumed with a separate config file (`-c`) holding the STAR_ALIGN `ext.args` override that is now part of the `nextflow.config` template, and dead intermediate files (the killed task directories, the unsorted BAM files, the STAR index copy) were deleted by hand during the run.
- Resource findings: the selector values in **Resources** above, `queueSize` 20, the trace fields `attempt` and `%cpu`, the misleading `peak_rss` of DEXSEQ_EXON, and the disk estimate in **Prerequisites** come from this run.
- Not verified by this test: STAR_GENOMEGENERATE (the index was reused); an unmodified STAR_ALIGN run to completion; the true memory peak of DEXSEQ_EXON; the new selector values and the new `nextflow.config` template (the run used the earlier values and a separate file); a paired design; single-end input; `forward` strandedness; genome-BAM input; the Ensembl download helper; the generated Salmon strandedness helper (the test used a hand-written job with the same Salmon command).

---

## Known limitations

- **Unreleased pin.** nf-core/rnasplice 1b447239488097651d8eac44bca2c1556865eb0f is a development commit (`dev-1b44723`), not a release. Release 1.0.4 does not launch under Nextflow 26.04.6 (config parse error on `def check_max`, seen at the verification gate), so it was not used. The pinned commit is the merge of upstream PR #291 (2026-09-24), which ported the rMATS subworkflow to the nf-core structure (a rewrite of the subworkflow; evidence in `tests/fixtures/gate_report.md`): young code, a risk. (Gate observation, for a later choice of pin: release 1.0.4 on Nextflow 24.04.4 completed the test profile, but one rMATS 4.1.2 task hung for 43 minutes until it was killed by hand (the retry then finished in 2 s); cause not verified.)
- **Bumping the pin** (a newer commit or a release) means running the verification gate again before relying on the skill: the test profile and the gate runs above on the cluster; record the schema (`rnasplice_schema.json`, `schema_input*.json`), the pipeline's `params {}` block (`config_params.txt`), the process names and requested resources (`trace_process_names.txt`, `trace_resources_g3.tsv`), the output tree (`output_tree.txt`), the task command lines (`command_lines.txt`: strandedness, read type, bamlist order) and the container URLs (`verified_urls.txt`); set every key of `gate_values.tsv` (the checker validates each against its allowed values) and update `gate_report.md`; check every sign on the values or the code, never from a header alone; then update the skill until `tests/run_all_tests.sh --full` passes.
- **Option keys first run in the cluster acceptance run:** `ignore_tx_version`, `miso_genes`, `miso_read_len`, `fig_height`, `fig_width`, `isoformswitchanalyzer_alpha` and `isoformswitchanalyzer_dIF` are in the params file (at the pipeline defaults) but were not in the gate's typed run; they were first used in the cluster acceptance run (run (a) below). `star_index` (reuse of an existing STAR index) was first run in the real-data test.
- **Single-end input and `forward` strandedness:** single-end BAM input and `forward` strandedness were verified from the pipeline code only, not by a run (the gate ran paired-end samples only: `unstranded` and `reverse`, from FASTQ and from genome BAM); single-end FASTQ input was not run either.
- One contrast was run with the skill's settings; several contrasts only in the pipeline's own test profile.
- Strandedness is asked, not detected automatically: the wizard proposes a value only from the RSeQC results of an existing nf-core/rnaseq run or from the optional Salmon helper (FASTQ input, needs a Salmon index of the organism in the shared genome folder), and the user confirms it; rMATS needs one strandedness and one read type for all samples. The Salmon helper was run by hand in the real-data test (reverse-stranded paired-end data); its generated form was tested with stubs only.
- Paired designs need exactly two conditions of equal size.
- FASTA, GTF and BAM contig names must match (for example all Ensembl `1` or all UCSC `chr1`); the skill does not check it. With FASTQ input a mismatch stops the pipeline when it builds the STAR index; with genome-BAM input it gives no error, only zero counts and empty rMATS tables (compare the `@SQ` names of `samtools view -H` of one BAM with the FASTA headers and the first GTF column, on a compute node).
- Inputs: FASTQ or genome BAM only. No Salmon-results or transcriptome-BAM input, and no iGenomes `genome` key (FASTA and GTF are always given). Organisms other than human and mouse need a custom FASTA/GTF.
- With FASTQ input, Salmon (index and quantification) always runs, also when neither DTU nor SUPPA2 is chosen.
- Trimming, QC and SUPPA2 clustering use the pipeline defaults; edgeR exon usage uses diffSpliceDGE glmQLFit glmQLFTest in this revision.
- MISO (dropped: see Key design points), IsoformSwitchAnalyzeR and LeafCutter are off and not offered.
- The skill does not analyse the results; there is no downstream report skill for the splicing tables yet (planned later).

---

## rnasplice or `/bulk-rnaseq-pipeline`

`/bulk-rnaseq-pipeline` does gene-level differential expression, GSEA and, optionally, differential transcript usage (DRIMSeq → DEXSeq → stageR) on the Salmon output of an nf-core/rnaseq run. rnasplice answers which exons and splice events change (rMATS, SUPPA2, DEXSeq/edgeR exon usage). Its DEXSeq DTU runs the same steps as the bulk skill's DTU module (on another Salmon quantification): run it in one place only.
