# Design: `nfcore-rnasplice-setup` skill (nf-core/rnasplice on SLURM + Singularity)

Status: design agreed with the user section by section on 2026-10-01; this file is the written spec for review. Facts about the pipeline come from `.superpowers/sdd/rnasplice-research.md` (fetched 2026-10-01, git-ignored) and are re-verified in Task 0.

## Purpose

An interactive Claude Code wizard, in the style of `nfcore-rnavar-setup`, that configures and submits **nf-core/rnasplice** for differential alternative-splicing analysis of bulk RNA-seq on the Koch cluster (SLURM, Singularity, Nextflow from a conda environment). It replaces the user's earlier MISO + rMATS workflows (example: Akerberg et al. 2022, Circ Res, PMC9770155: rMATS for all five event types, FDR < 0.05 or 0.1, |IncLevelDifference| > 0.1) with the rnasplice toolset. It does not analyse the results: a downstream report skill is a separate, later project.

## Decisions already made (with the user, 2026-10-01)

| Decision | Choice |
|---|---|
| Skill shape | Standalone `nfcore-rnasplice-setup` (folder with `.md`, `README.md`, `tests/`), same pattern as `nfcore-rnavar-setup`: all pipeline parameters in a **params file**, explicit true/false for every module, a launch line that carries no `--flag` |
| Modules offered | **rMATS** (default on) plus opt-in **SUPPA2**, **DEXSeq / edgeR differential exon usage**, **DEXSeq DTU** (DRIMSeq + DEXSeq + stageR). The DTU option warns that it duplicates the optional DTU module of `bulk-rnaseq-pipeline` |
| MISO | **Dropped.** MISO is archived (no commits since 2019, Python 2) and rnasplice uses it only for sashimi plots, not PSI. `sashimi_plot` is always written `false`; the wizard says why. A standalone MISO (Bayes factor) workflow is out of scope |
| Input | **FASTQ and BAM** (`source` = `fastq` or `genome_bam`), with the release-1.0.4 limitations of BAM input stated and handled (see BAM input) |
| Downstream report skill | Later, after the setup skill, once a real run shows the output layout |
| Relationship to `bulk-rnaseq-pipeline` | Kept as is. The wizard states when to use which: DE and DTU on Salmon output from `nfcore-rnaseq-setup` stay in `bulk-rnaseq-pipeline`; event-level splicing (cassette exons, A5/A3, RI, MXE) goes to rnasplice |

## Architecture

New folder `nfcore-rnasplice-setup/`: `nfcore-rnasplice-setup.md` (the skill), `README.md` (steps, outputs, validation status, known limitations), `tests/check_skill.sh` (static checker), `tests/fixtures/` (recorded `nextflow_schema.json` and the effective `nextflow.config` defaults of the pinned release, the nf-core test sheets, verification notes). One row is added to the root `README.md`. The skill is installed by copying the `.md` to `~/.claude/commands/`, like the others.

The skill reuses, step by step, the structure of `nfcore-rnavar-setup` (Steps 0-4 working directory, email, conda environment, version, input scan and samplesheet building with sample-name sanitisation, naming, overwrite checks; the shared genome folder convention `{genome_base}/{organism}/{assembly}_ens{version}/`; the params-file typing rules; the `nextflow.config` skeleton with SLURM selectors; the submission script; helper scripts; the Lmod hazard rule: a failed `module add` silently breaks later ones, so every helper guards each one with `|| exit 1`).

## Task 0: verification gate (before any skill text)

Run the nf-core `test` profile of the candidate release on a compute node with the cluster's Nextflow (26.04.6): `-profile test,slurm,singularity`. The test data is 4 paired-end human chrX samples (about 30 MB FASTQ) with a 46 MB FASTA and a 26 MB GTF; the test profile turns every module on, including MISO with human gene IDs. The gate decides and records:

1. **Does release 1.0.4 launch under Nextflow 26.04.6?** It uses `nf-validation` 1.1.3 and `check_max`; the `dev` branch (1.1.0dev) was raised to Nextflow >= 26.04 and `nf-schema` 2.7.2. If 1.0.4 fails, the options are `-r dev` (unreleased, a moving target: pin a commit) or an older Nextflow in the conda environment; the user decides before the plan proceeds.
2. **Real process names, resources and the output tree** (for the `withName` selectors and for the README's output table). The research could not verify the STAR index `sjdbOverhang`, the label of the STAR alignment step, the edgeR DEU function, or exact output file names.
3. **Whether the pipeline accepts a params file with every parameter typed** (the Nextflow 26.04 string-typing problem that rejected numeric CLI parameters in rnavar), and which parameters the validator rejects.
4. **How strandedness reaches rMATS for BAM input**: release 1.0.4 BAM sheets (`sample,condition,genome_bam`) have no strandedness or `single_end` column (added only in dev). The gate finds out what rMATS is run with, and whether a BAM run on the test data succeeds.
5. Container availability for rMATS, DEXSeq, stageR, SUPPA (the Singularity pulls) and whether the default `miso_genes` exist in the test GTF.

If the gate fails in a way that makes the skill impossible, the project stops and says so, as for phASER.

## Wizard (planned steps)

0. Working directory, `{WD_NAME}`, `{TODAY}`. 1. Email. 2. Conda environment (asked separately). 3. Latest release by `WebFetch` of the GitHub API (never `gh`), schema re-validation for newer releases: the version defaults to the gated release.
4. **Input source**: FASTQ or BAM.
   - FASTQ: scan, detect paired-end, build `sample,fastq_1,fastq_2,strandedness,condition` with the rnavar sanitisation and uniqueness rules (rows with the same sample name are technical replicates and are merged). **Strandedness is asked** (`unstranded` | `forward` | `reverse`; the pipeline has no `auto`); all samples must share strandedness and read type, because rMATS requires it.
   - BAM: `sample,condition,genome_bam`, from a splice-aware aligner (for example the STAR BAMs of `nfcore-rnaseq-setup`); the 1.0.4 limitations from Task 0 are shown and enforced (what the wizard does about strandedness and single-end).
5. **Design**: assign `condition` (asked, with a plain-language description of the groups); build the **contrasts sheet** `contrast,treatment,control` (header required; every treatment/control label must be a value of `condition`): all pairwise or user-defined, with validation. Ask whether the design is **paired** (needed for `rmats_paired_stats`; the 1.0.4 default `true` is a trap: the wizard writes `false` unless the user confirms pairing and two conditions).
6. **Read length**: `rmats_read_len` is always detected from the FASTQ (`zcat | awk` on several samples; for BAM, asked); the 1.0.4 default 40 is wrong for most data. There is no `read_length` or `sjdbOverhang` parameter, so unlike rnavar there is no per-read-length STAR index.
7. **Organism and genome**: reuse the shared folder; FASTA and GTF are required; Ensembl/GENCODE detection (`gencode: true` for GENCODE transcript FASTA); existing `star_index` and `salmon_index` may be passed if the user has them (asked; otherwise the pipeline builds its own, which needs the heavy STAR resources below).
8. **Modules**, every switch explicit: `rmats`, `dexseq_exon`, `edger_exon`, `dexseq_dtu`, `suppa`, `sashimi_plot` (false). In 1.0.4 `nextflow.config` all of them default to **true** while the schema shows no default, so an omitted switch silently runs everything. rMATS options asked: `rmats_novel_splice_site`, `rmats_splice_diff_cutoff` (default 0.0001), `rmats_paired_stats`. DTU filters (`min_samps_gene_expr`, `min_samps_feature_expr`, `min_samps_feature_prop`, `min_gene_expr`, `min_feature_expr`, `min_feature_prop`, `dtu_txi`) are written explicitly because the schema (6/0/0) and the config (4/2/2) disagree; the wizard states which it uses and why. SUPPA2 and the DEU modules: defaults, written explicitly.
9. Trimming and QC: pipeline defaults, not tuned in v1 (written only if the user changes them).
10. MultiQC title, output directory, `nextflow.config` (SLURM selectors; heavy steps `RMATS_PREP`/`RMATS_POST`, `DEXSEQ_*`, STAR alignment and index generation; `max_cpus`/`max_memory`/`max_time` go in the params file because 1.0.4 still has them), params file, submission script. Helper scripts only for missing resources (genome download).
11. Final summary: what was written, how to submit, which output answers which question (events: rMATS and SUPPA2; exon usage: DEXSeq and edgeR; transcript usage: DEXSeq DTU), and where the DTU overlap with `bulk-rnaseq-pipeline` is.

## Parameter policy

All pipeline parameters are written to `{date}_{project}_params.yaml` (typing rules of `nfcore-rnavar-setup`: paths and strings quoted, numbers and booleans bare, `multiqc_title` always quoted, empty placeholder lines deleted). Never emitted, because they are not parameters of 1.0.4: `isoformswitchanalyzer`, `leafcutter`, `ignore_tx_version`, `rmats_variable_read_len` (rMATS always runs with `--variable-read-length --allow-clipping`); the doc typo `--local_events` is not a parameter (the real one is `generateevents_event_type`). The checker confirms every emitted parameter exists in the recorded schema of the gated release.

## Testing and acceptance

- **Static checker** `tests/check_skill.sh`, in the style of `nfcore-rnavar-setup/tests/check_skill.sh` (no python on the host): every parameter name used in the skill is in the recorded schema or an allowlist; need/forbid strings per task; no `--flag` on the launch line; mutation proofs for every structural pin.
- **Cluster acceptance** on the nf-core test data with the skill's own params file: (a) all modules on; (b) rMATS only; (c) a BAM-input run, if the gate shows it is feasible. Each must complete; the process selectors are compared with the real task names; the README records the run honestly (Nextflow version, release, what was verified).
- No real biological data is part of the acceptance (the test dataset is tiny); the README says so.

## Risks and known limitations

Release 1.0.4 is about 2.5 years old and has had no release since (Task 0 decides pinning). BAM input in 1.0.4 has no strandedness column. rMATS needs identical strandedness and read type for all samples and runs one prep/post pair per contrast; `rmats_read_len` must be right. rMATS and DEXSeq steps are `process_high` (12 CPUs, 60 GB, 16 h by default), and the STAR index step is heavy for a human genome. Output file names and columns are not verified until Task 0.

## Out of scope for this project

MISO as an analysis; a downstream report skill; `salmon_results` and `transcriptome_bam` input sources; iGenomes (`--genome`); IsoformSwitchAnalyzeR and LeafCutter (dev only); trimming and SUPPA clustering tuning; cross-species orthology.

## Staging

A single stage, with the Task 0 gate first. The process is the same as for the earlier skills: plan, subagent-driven implementation with per-task review, a final whole-branch review, a live cluster acceptance run, and the user's approval before any merge or push.
