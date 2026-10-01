# nfcore-rnasplice-setup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an interactive `nfcore-rnasplice-setup` Claude Code skill that configures and submits nf-core/rnasplice (differential alternative splicing: rMATS, SUPPA2, DEXSeq/edgeR exon usage, DEXSeq DTU) on the Koch cluster (SLURM + Singularity, Nextflow from conda), with every pipeline parameter in a params file and every module switch explicit.

**Architecture:** The deliverable is one self-contained Markdown skill (`nfcore-rnasplice-setup/nfcore-rnasplice-setup.md`) plus a README, in the pattern of `nfcore-rnavar-setup`. A verification gate (Task 0) runs the pipeline on the cluster first and records its real behaviour as fixtures (`tests/fixtures/`: schema, config defaults, process names, output tree, gate values). The skill's logic that can go wrong silently (sheet validation, strandedness inference, BAM policy, read-length detection, params rendering, submission script, download helper) lives in literal code blocks that the tests cut out of the skill and run; a static checker (`tests/check_skill.sh`, bash + awk, no python, no network) pins structure against the fixtures, and every structural pin is proven by a committed mutation list. A controller-run acceptance on the nf-core test data closes the loop.

**Tech Stack:** Markdown skill file; bash, awk (GNU Awk 4.0.2), grep, sed, gzip on the host (no python, no jq); nf-core/rnasplice 1.0.4 (or the gate's fallback), Nextflow 26.04.6 in conda env `nf-env`, Singularity 3.10.4, SLURM partition `bcc`.

**Spec:** `docs/superpowers/specs/2026-10-01-nfcore-rnasplice-setup-design.md` (pipeline facts: `.superpowers/sdd/rnasplice-research.md`, git-ignored; template: `nfcore-rnavar-setup/`).

All repo paths are relative to `/net/bmc-lab3/data/bcc/yannvrb/Github_repos/bioinformatics-claude-skills`. `D` below means `nfcore-rnasplice-setup`.

## Global Constraints

- User rules (from `~/.claude/CLAUDE.md`, binding on every task and every generated script): raw data is read-only (FASTQ, BAM, `*.h5ad`); outputs go to `results/YYYY-MM-DD_*`; never run heavy work on login nodes: Nextflow, Java, conda, Python and R run only on compute nodes through `sbatch` (debug with `salloc`/`srun`); default SLURM requests at most `--mem=64G` and `-t 4:00:00` unless the user asked otherwise; show `git status` and `git diff` before every commit; never `git push` without the user's explicit approval.
- On the login node only: bash, awk, grep, sed, `zcat | head`, `gzip -t`, git, `curl`/`WebFetch` of small text files (schemas, configs, GitHub API JSON). No python on the host (`/usr/bin/python` exists; never use it). No `gh` CLI in the skill (use `WebFetch`).
- Commit messages end with the attribution line given in the session's system-reminder.
- Pipeline: `nf-core/rnasplice`, the revision recorded by Task 0 in `D/tests/fixtures/gate_values.tsv` (`PIPELINE_REVISION`), default release 1.0.4.
- Every pipeline parameter goes in `{PARAMS_YAML}`; the launch line is exactly `nextflow run nf-core/rnasplice -r {VERSION} -c nextflow.config -profile slurm,singularity -params-file {PARAMS_YAML}` (no `--flag`). Params-file typing: paths and strings double-quoted, numbers and booleans bare, `multiqc_title` always quoted; only the `star_index` and `salmon_index` lines are deleted when empty.
- Every module switch is written explicitly: `rmats`, `dexseq_exon`, `edger_exon`, `dexseq_dtu`, `suppa`, and `sashimi_plot: false` (MISO dropped). Never emitted (not parameters of 1.0.4): `isoformswitchanalyzer`, `leafcutter`, `ignore_tx_version`, `rmats_variable_read_len`, `local_events`; no iGenomes `genome`. Every emitted key must exist in the recorded schema of the gated revision.
- `rmats_read_len` is always detected (FASTQ) or asked (BAM) and written; strandedness is always asked (`unstranded`|`forward`|`reverse`), one value for all samples; `rmats_paired_stats` is written `false` unless the user confirms a paired design.
- `nextflow.config` is written only if absent, never overwritten. Existing samplesheet, contrasts sheet, params file, submission script or helper: ask `1. overwrite · 2. choose another filename` before writing.
- Every `module add` in the skill is guarded (`|| { echo "ERROR: ..." >&2; exit 1; }`): under Lmod a failed `module add` silently breaks every later one.
- Generated scripts and tests delete only files they created themselves, by explicit name, or a `mktemp -d` directory whose path was checked first; never `rm -rf` of a variable that could be empty or user-supplied.
- Tests that execute extracted shell code run under `env -i` with `PATH=<stub dir>:/usr/bin:/bin`, type-check every stub before the first scenario, and stub every command that must never run on the login node (nextflow, conda, java, python, R, singularity, sbatch).
- An implementer never weakens a gate, a test or a checker line to make it pass: report `BLOCKED` with the evidence; the controller rules and records the ruling in the ledger.
- `cut_block.sh` and `eval` are used only on the skill file, never on this plan: the plan mentions every anchor in prose before its block, so a cut from the plan returns the wrong block (the planner once ran the G0 fetch script that way by accident; it was cleaned up).
- The skill is self-contained: it may name `/bulk-rnaseq-pipeline` and `/nfcore-rnaseq-setup` (spec requirement) but never refers to their step numbers.

## Review Focus

1. **A module switch left out of the params file** → the 1.0.4 `nextflow.config` defaults every module to `true`, so an omitted switch silently runs everything (hours of DEXSeq/MISO on a human genome). Expected: the params file always has all six switch lines with explicit values. Pinned in Task 5 (checker rule "module switch ... exactly once", mutations T5-switch, T5-dup, T5-sashimi) and Task 9 (rMATS-only run: the trace has no DEXSeq, edgeR, DRIMSeq, stageR, SUPPA or MISO task).
2. **Wrong strandedness** (for example `unstranded` declared for a dUTP library, or mixed values) → rMATS runs with the wrong `--libType` and loses or misassigns junction reads, with no error. Expected: strandedness is asked, never defaulted, one value for all samples, checked against RSeQC output when an nf-core/rnaseq run exists. Pinned in Task 2 (`test_strandedness.sh`, mutations T2-thr, T2-swap) and Task 3 (`test_validate_sheets.sh` cases `mixed_strandedness`, `strandedness_auto`); Task 0 G4 records how a declared `reverse` reaches rMATS.
3. **Wrong `rmats_read_len`** (the default 40, or the length of one odd file) → rMATS counts reads against the wrong length. Expected: the most common read length over several samples, written as `rmats_read_len: {READ_LENGTH}`. Pinned in Task 4 (`test_read_length.sh` cases `mixed`, `tie`; mutations T4-tie, T4-mode), Task 5 (checker line, mutation T5-readlen) and Task 9 (the RMATS_PREP command line carries the detected length).
4. **A contrasts sheet whose treatment/control labels are not values of the `condition` column** (typo, case difference, `-` vs `_`) → the pipeline fails late or skips the contrast. Expected: the wizard validates both sheets before writing them and names the bad label. Pinned in Task 3 (`test_validate_sheets.sh` cases `label_typo_treatment`, `label_typo_control`, `contrasts_no_header`, `self_contrast`, `duplicate_contrast`; mutations T3-label, T3-header).
5. **BAM input without strandedness** (release 1.0.4 BAM sheets have no strandedness or single-end column) → rMATS, DEXSeq and featureCounts run with whatever the revision hard-codes, possibly the wrong strand for a stranded library. Expected: BAM input is allowed only when the library matches what the pinned revision does with BAM input (recorded by Task 0 G5), otherwise refused with the advice to start from FASTQ. Pinned in Task 2 (`test_bam_policy.sh` 3 x 2 matrix against the gate values; checker lines tying the function's literals to the gate values; mutation T2-bam-lt) and Task 9 (BAM run).

## Planner decisions (need user review)

1. Strandedness is asked (`unstranded`/`forward`/`reverse`), never detected from the reads; the wizard can propose it from an existing nf-core/rnaseq run's RSeQC `infer_experiment` files. Reason: rnasplice has no `auto`, and detection would need an extra alignment job.
2. `rmats_paired_stats` and SUPPA's `diffsplice_paired` default to `false`; `true` only for a confirmed subject pairing with exactly two conditions and equal group sizes, and only if Task 0 shows a fixed rMATS sample order. Reason: the upstream `true` silently assumes pairing.
3. Existing indexes in the shared genome folder (`index/star/`, `index/salmon/`, as built by `/nfcore-rnaseq-setup`) are offered for reuse only when their index-format version equals the one recorded by Task 0; otherwise the pipeline builds its own and does not save it (`save_reference: false`). Reason: an incompatible index fails late or silently.
4. One helper script only, `download_genome_{REF_TAG}.sh`, for a missing Ensembl FASTA/GTF; no index-build helper. Reason: the pipeline builds STAR/Salmon indexes itself.
5. Module menu defaults: rMATS on; SUPPA2, DEXSeq DEU, edgeR DEU and DEXSeq DTU off unless chosen; `sashimi_plot: false` always. Reason: spec; avoids duplicating the bulk-rnaseq-pipeline DTU by default.
6. Every option key of every module is written to the params file (also for modules switched off), at the recorded config default except named deviations. Reason: a complete record that the checker can compare with the recorded config.
7. DTU filters: `min_samps_gene_expr` = number of samples (per contrast if Task 0 shows per-contrast filtering), `min_samps_feature_expr` = `min_samps_feature_prop` = smallest compared condition, `min_gene_expr` 10, `min_feature_expr` 10, `min_feature_prop` 0.1 (Love et al. 2018, as in bulk-rnaseq-pipeline). Reason: neither the schema (6/0/0) nor the config (4/2/2) fits a given design.
8. Every compared condition needs at least 2 samples (refused otherwise); sample and condition names use letters, digits and `_`, start with a letter (a leading digit gets an `S` prefix). Reason: DEXSeq/edgeR need replicates; R mangles names that start with a digit.
9. All-pairwise contrasts follow the user's reference-first condition order: control = the earlier condition, name `{treatment}_vs_{control}`.
10. Salmon route chosen by Task 0, preferring `aligner: "star"` + `pseudo_aligner: "salmon"` so DTU and SUPPA2 run once (from Salmon) instead of twice. Reason: the 1.0.4 defaults (`star_salmon` + `salmon`) can run both Salmon branches.
11. File names: `{prefix}_samplesheet.csv`, `{prefix}_contrasts.csv`, `{prefix}_params.yaml`, `nf-core_rnasplice_{VERSION_TAG}.sh`; default output directory `results/{TODAY_ISO}_{WD_NAME}` (user rule).
12. Version check: the submission script reads `nextflow -version` on the compute node and exits 1 outside the range verified by Task 0; the wizard never runs Nextflow on the login node.
13. Pipeline head job `-n 2 --mem=8G -t 48:00:00`, beyond the 4 h default, because the head job must outlive every task (tasks are capped at 24 h).
14. The submission script runs `NXF_OFFLINE=false nextflow pull nf-core/rnasplice -r {VERSION}` before the run (the user's `~/.bashrc` exports `NXF_OFFLINE=TRUE`, which blocks a first download) and continues with the cached copy if the pull fails.
15. BAM input of a library that the pinned revision cannot describe (strandedness or read type) is refused, not allowed with a warning.
16. `gencode: false` always: the pipeline builds the transcript FASTA from the GTF (gffread), so no GENCODE-format transcript FASTA is ever given.
17. Selector resources: STAR_GENOMEGENERATE 8 CPU/64 GB/8 h, STAR_ALIGN 8/48 GB/8 h, RMATS_PREP 4/16 GB/8 h, RMATS_POST 8/32 GB/16 h, DEXSEQ_COUNT 2/8 GB/8 h, DEXSEQ_EXON and DEXSEQ_DTU 8/32 GB/8 h, SALMON_QUANT 8/16 GB/4 h; caps 16 CPU/64 GB/24 h. Reason: judgement, not measured on real data (README says so).

## Proven facts (from the research report; source per line)

Verified from raw files read with curl on 2026-10-01 (`.superpowers/sdd/rnasplice-research.md` §0 lists the URLs):
- Latest release 1.0.4 (2024-05-09; GitHub API `releases/latest`); `nextflowVersion = '!>=23.04.0'`, plugin `nf-validation@1.1.3` (1.0.4 `nextflow.config`). `dev` = 1.1.0dev, `!>=26.04.0`, `nf-schema@2.7.2` (dev `nextflow.config`).
- In 1.0.4 `nextflow.config`: `rmats`, `dexseq_exon`, `edger_exon`, `dexseq_dtu`, `suppa`, `sashimi_plot` all `true`; `aligner = 'star_salmon'`, `pseudo_aligner = 'salmon'`; the schema shows no default for the switches and `aligner` default `"star"` (research §2.2).
- DTU filters: config 4/2/2 (`min_samps_gene_expr`/`min_samps_feature_expr`/`min_samps_feature_prop`), schema 6/0/0, usage.md "all set to 0" (research §2.2).
- `rmats_paired_stats` default `true`; `rmats_read_len` default 40; `rmats_splice_diff_cutoff` 0.0001; `rmats_novel_splice_site` false; `rmats_min_intron_len` 50, `rmats_max_exon_len` 500, passed only with novel splice sites (rmats_prep.nf lines 48-56) (research §2.4).
- rMATS always runs with `--variable-read-length --allow-clipping` (rmats_prep/post); strandedness mapping in rmats_prep.nf: `forward` -> `fr-secondstrand`, `reverse` -> `fr-firststrand`, `unstranded` -> `fr-unstranded` (research §6 item 8, verified).
- Required: `input`, `source`, `outdir`, `contrasts` (schema group `input_output_options`); plus `fasta` and `gtf` or `gff` (`WorkflowRnasplice.initialise`) (research §2.3).
- Sheets: FASTQ `sample,fastq_1,fastq_2,strandedness,condition` (strandedness enum unstranded|forward|reverse, no `auto`; same `sample` = technical replicates merged by CAT_FASTQ; `.fastq.gz`/`.fq.gz`); genome BAM `sample,condition,genome_bam`; contrasts `contrast,treatment,control`, header required, labels must be `condition` values; test contrasts `GBR-YRI,GBR,YRI` and `YRI-GBR,YRI,GBR` (research §2.5).
- rMATS from FASTQ needs one strandedness and one read type ("Cannot run rMats with mixed single and paired end samples" / "mixed stranded samples"); `clusterevents_local_event` needs `diffsplice_local_event`, `clusterevents_isoform` needs `diffsplice_isoform` (research §2.4).
- 1.0.4 still has `max_cpus` (16), `max_memory` ("128.GB", pattern `^\d+(\.\d+)?\.?\s*(K|M|G|T)?B$`), `max_time` ("240.h") and `check_max` in `base.config` (research §2.4).
- Labels in 1.0.4: RMATS_PREP/RMATS_POST, DEXSEQ_EXON, DEXSEQ_DTU = `process_high` (12 CPU/60 GB/16 h); STAGER, MISO_RUN = `process_medium`; SUPPA_GENERATEEVENTS = `process_low`; EDGER_EXON, MISO_SASHIMI = `process_single` (research §2.7).
- Test profile: 4 paired-end human chrX samples (ERR188383, ERR188428 = GBR; ERR188454, ERR204916 = YRI), `unstranded`; `X.fa.gz` 45.8 MB, `genes_chrX.gtf` 25.7 MB, FASTQ 3.6-3.8 MB each; `max_cpus = 2, max_memory = '6.GB', max_time = '6.h'` (research §5).
- Cluster: Nextflow 26.04.6 in conda env `nf-env`; the rnavar run used `module add miniconda3/v4; source /home/software/conda/miniconda3/bin/condainit; conda activate nf-env; module add singularity/3.10.4` and `-c nextflow.config` (`/net/bmc-lab3/data/bcc/yannvrb/rnavar_test2/nf-core_rnavar_1.3.0.sh`, job 11375465). The user's `~/.bashrc` exports `NXF_SINGULARITY_CACHEDIR=/home/yannvrb/.singularity/cache` and `NXF_OFFLINE='TRUE'` (read 2026-10-01). Plugins `nf-validation-1.1.3` and `nf-schema-2.7.2` are already in `~/.nextflow/plugins/`; nf-core/rnasplice is not yet in `~/.nextflow/assets/`.
- Under Nextflow 26.04 numeric CLI parameters reached the nf-schema validator as strings (`--read_length 151` rejected, rnavar run 2026-09-28); a params file avoided it.

NOT VERIFIED until Task 0 (each is a gate item):
- Whether 1.0.4 launches under Nextflow 26.04.6 (G1); whether every typed parameter of the params file is accepted (G3).
- Real process names, labels, requested resources; STAR_ALIGN label; STAR `sjdbOverhang` of STAR_GENOMEGENERATE; STAR and Salmon index-format versions (G3, G8).
- Output tree and file names (rMATS JC/JCEC tables, DEXSeq, edgeR, SUPPA, DTU) (G3).
- Which bamlist is rMATS `b1` (sign of IncLevelDifference) and whether the bamlist order follows the samplesheet (G3).
- Whether DTU/SUPPA run on both Salmon branches with the defaults, and whether `aligner: "star"` + `pseudo_aligner: "salmon"` runs them once (G1, G3).
- DRIMSeq filter scope (all samples or per contrast); edgeR DEU function (G3, from the pipeline's `bin/` R scripts).
- How strandedness and read type reach rMATS, DEXSeq count and featureCounts for genome BAM input in 1.0.4; whether `.bai` files are needed (G5).
- Container availability (rMATS, DEXSeq, stageR, SUPPA, edgeR, MISO images) (G6); whether the default `miso_genes` exist in the test GTF (G7).
- Whether `NXF_OFFLINE=TRUE` blocks the first run, and whether compute nodes reach GitHub and depot.galaxyproject.org (G0b, G1).

---

## Conventions used in this plan

- `«KEY»` inside skill text to be written means: substitute the literal value of `KEY` from `D/tests/fixtures/gate_values.tsv` (written by Task 0) when writing the skill. The skill itself never contains `«`; the checker verifies the substituted literals (it reads the same file). Example: `«PIPELINE_REVISION»` becomes `1.0.4` if the gate pinned release 1.0.4.
- `{NAME}` is a wizard placeholder that stays in the skill text (bound by an earlier wizard step, substituted by the assistant when it writes files for the user).
- Paragraphs marked **[A/C]** are written when `GATE_OUTCOME` is `A` (1.0.4 works under Nextflow 26.04.6) or `C` (1.0.4 with an older Nextflow); **[B]** when it is `B` (pinned dev commit). Only one variant goes into the skill.
- Anchors: a bold line such as `**Sheet validation.**` that precedes a fenced code block is an anchor; `tests/cut_block.sh <skill> "<anchor>"` prints the first fenced block after it. Every anchor appears exactly once in the skill (checker rule `anchor_once`).
- `BASE` in a task = the commit at which the task started (`git rev-parse HEAD` before its first edit).

## File Structure

| File | Responsibility | Task |
|---|---|---|
| `D/nfcore-rnasplice-setup.md` | The skill (Steps 0-12, Notes) | 1-7 |
| `D/README.md` | Steps, outputs, validation status, limitations | 8, 9 |
| `D/tests/check_skill.sh` | Static checker against the fixtures | 1, extended 2-8 |
| `D/tests/cut_block.sh` | Prints the fenced block after an anchor (copy of `ase-pipeline/tests/synthetic/cut_block.sh`) | 1 |
| `D/tests/prove_red.sh` | Proves every new `need` fails on the skill at `BASE` | 1 |
| `D/tests/prove_mutations.sh`, `D/tests/mutations.tsv` | Proves every structural pin catches its mutation | 1, rows added 2-8 |
| `D/tests/test_strandedness.sh` | Runs `infer_strandedness` from the skill on RSeQC fixtures | 2 |
| `D/tests/test_bam_policy.sh` | Runs `bam_input_allowed` from the skill against the gate values | 2 |
| `D/tests/test_validate_sheets.sh` | Runs `validate_rnasplice_sheets` from the skill on good/bad sheets | 3 |
| `D/tests/test_read_length.sh` | Runs `detect_read_length` from the skill on generated FASTQ | 4 |
| `D/tests/render_params.sh`, `D/tests/test_render_params.sh` | Renders the params file exactly as Step 11 says; tests it | 6 |
| `D/tests/dry_run_submit.sh` | Runs the submission script from the skill against stubs (`env -i`) | 6 |
| `D/tests/dry_run_helpers.sh` | Runs the download helper (and, branch C, the env helper) against stubs | 6 |
| `D/tests/run_all_tests.sh` | Runs everything; prints `ALL TESTS PASS` | 8 |
| `D/tests/fixtures/*` | Recorded gate evidence (list in Task 0) | 0, 6 |
| `README.md` (root) | Skills table row and intro bullet | 1 (stub), 8 |

## Ledger instructions for the controller

- Ledger: `.superpowers/sdd/2026-10-01-nfcore-rnasplice-setup/progress.md` (create it; first lines: plan path, spec path, branch `feat/nfcore-rnasplice-setup`, base commit, "No push without user approval", standing rules: Nextflow/conda/Java/R/python only via `sbatch` on compute nodes; the controller never fixes skill text itself).
- Before Task 1: a pre-flight table (pairs of tasks sharing a file or interface, as in the earlier ledgers) and the gate outcome line `GATE: <A|B|C|STOP> revision <PIPELINE_REVISION> Nextflow <NEXTFLOW_TESTED>; user decisions: <list>`.
- Per task: `Task N: implemented (commit X, base Y) <status>`; review result; each Important/Critical finding with a `Ruling:` line that states the fix and `Cost if wrong: ...`; deferred minors as `minor (deferred): ...`; `Task N: complete (commits Y..X, review clean after k fix rounds)`.
- `BLOCKED` from an implementer: record the evidence, rule (never by weakening a test or gate), record the ruling, re-dispatch.
- Tasks 0 and 9 are controller-run (they need cluster jobs and the user's email/conda env: `yannvrb@mit.edu`, `nf-env`, as used in the rnavar acceptance); the controller may dispatch one agent to run the gate jobs but checks every result file itself. Any incident (a heavy command on the login node, a deletion outside a scratch path) is recorded verbatim with its cleanup.
- After Task 10: record the final review verdict, the fix wave, the re-review, and the user's install/merge decision. Nothing is pushed.

---

### Task 0: Verification gate (controller-run; ends with a CONTROLLER CHECKPOINT)

Nothing of the skill is written before this task is finished and the user has answered the checkpoint. Every Nextflow, conda and container command runs inside an `sbatch` job; the login node only fetches small text files with `curl`, reads result files with `grep`/`awk`, and writes the fixtures.

**Files:**
- Create (scratch, outside the repo): `/net/bmc-lab3/data/bcc/yannvrb/rnasplice_gate/` with `src/`, `data/`, `jobs/`, `results/` (results named `results/YYYY-MM-DD_<run>`).
- Create: `D/tests/fixtures/rnasplice_schema.json`, `config_params.txt`, `gate_values.tsv`, `gate_report.md`, `trace_process_names.txt`, `output_tree.txt`, `test_samplesheet.csv`, `test_contrastsheet.csv`, `command_lines.txt`, `verified_urls.txt`, `typing_test.txt`.

**Interfaces:**
- Consumes: nothing from the repo except this plan.
- Produces: `D/tests/fixtures/gate_values.tsv`, one `KEY<TAB>VALUE` line per key, read by the checker (Task 1 on) and by every `«KEY»` substitution. Keys and allowed values (all required unless marked informational):

| Key | Allowed values | Used by |
|---|---|---|
| `GATE_OUTCOME` | `A` (1.0.4 + Nextflow 26.04.6), `B` (pinned dev commit + 26.04.6), `C` (1.0.4 + older Nextflow), `STOP` | all |
| `PIPELINE_REVISION` | `1.0.4` or a 40-character commit | T1, T6 |
| `VERSION_TAG` | `1.0.4`, or `dev-` + first 7 characters of the commit | T1, T6 |
| `NEXTFLOW_TESTED` | the Nextflow version of the passing run, e.g. `26.04.6` | T1, T6, T8 |
| `NEXTFLOW_MIN` | the revision's manifest minimum without `!>=`, e.g. `23.04.0` | T1, T6 |
| `NEXTFLOW_MAX_EXCL` | `none` (A, B) or the first version family known to fail (C: e.g. `26.04.0`) | T1, T6 |
| `CONDA_ENV_TESTED` | `nf-env` (A, B) or the env created in G2b (C) | T1, T8 |
| `HAS_MAX_PARAMS` | `yes` if `max_memory` is a schema parameter, else `no` | T6 |
| `HAS_RESOURCE_LIMITS` | `yes` if `NEXTFLOW_TESTED` >= 24.04.0, else `no` | T6 |
| `SALMON_ROUTE` | `pseudo_only`, `star_salmon_only`, `star_salmon_both` | T5, T7 |
| `PSEUDO_OFF_LINE` | the params line that disables the separate Salmon pseudo-alignment (route `star_salmon_only`), else `none` | T5 |
| `BAM_SHEET_HEADER` | exact genome-BAM samplesheet header of the revision | T2, T3 |
| `BAM_RMATS_LIBTYPE` | `fr-unstranded`, `fr-firststrand`, `fr-secondstrand`, `from_sheet`, `none` | T2 |
| `BAM_RMATS_READTYPE` | `paired`, `single`, `from_sheet`, `none` | T2 |
| `BAM_DEXSEQ_STRAND` | `no`, `yes`, `reverse`, `from_sheet`, `none` | T2 |
| `BAM_FC_STRAND` | `0`, `1`, `2`, `from_sheet`, `none` | T2 |
| `BAM_NEEDS_BAI` | `yes`, `no`, `n/a` | T2 |
| `RMATS_BAMLIST_ORDER` | `sheet`, `sorted_by_name`, `unordered` | T3 |
| `RMATS_B1_GROUP` | `treatment` or `control` | T7 |
| `STAR_VERSION_GENOME` | `versionGenome` of the index STAR_GENOMEGENERATE built, e.g. `2.7.4a` | T4 |
| `SALMON_INDEX_VERSION` | `indexVersion` of the index SALMON_INDEX built | T4 |
| `DTU_FILTER_SCOPE` | `all_samples` or `per_contrast` | T5 |
| `TEST_CONTRAST` | `GBR_vs_YRI` (the contrast name used in G3) | T6, T7 |
| `TEST_READ_LENGTH` | the read length of the test FASTQ from the G3 step 1 one-liner | T4, T6, T9 |
| `EDGER_DEU_FUNCTION` | edgeR function(s) found in the pipeline's edgeR exon script | T7, T8 |
| informational (free text; `GATE_DATE` and `RUN_DIRS` are required because the README quotes them): `GATE_DATE`, `STAR_VERSION`, `SALMON_VERSION`, `RMATS_VERSION`, `MISO_GENES_IN_TEST_GTF`, `OFFLINE_PULL_NEEDED`, `COMPUTE_INTERNET`, `VALIDATOR_REJECTED`, `RUN_DIRS` | free text without tabs | T8 |

- [ ] **Step G0: Fetch the revision's source files (login node, `curl` of small text files only)**

```bash
GATE=/net/bmc-lab3/data/bcc/yannvrb/rnasplice_gate
mkdir -p "$GATE"/src/1.0.4 "$GATE"/data "$GATE"/jobs "$GATE"/results
cd "$GATE"/src/1.0.4 || exit 1
curl -sf https://api.github.com/repos/nf-core/rnasplice/git/trees/1.0.4?recursive=1 \
  | grep -oE '"path": *"[^"]*"' | sed -E 's/"path": *"//; s/"$//' > tree.txt
curl -sf https://api.github.com/repos/nf-core/rnasplice/git/ref/tags/1.0.4 > tag_ref.json
for f in nextflow_schema.json nextflow.config conf/base.config conf/modules.config conf/test.config \
         workflows/rnasplice.nf lib/WorkflowRnasplice.groovy assets/schema_input.json \
         $(grep -E '^(subworkflows/local/|modules/local/|bin/).*(rmats|bamlist|drimseq|dexseq|edger|stager|suppa|samplesheet|check_samplesheet)' tree.txt); do
  mkdir -p "$(dirname "$f")"
  curl -sfL "https://raw.githubusercontent.com/nf-core/rnasplice/1.0.4/$f" -o "$f" || echo "MISSING $f"
done
awk '/^params[ \t]*\{/ {f = 1} f {print} f && /^\}/ {exit}' nextflow.config > config_params.txt
grep -n "nextflowVersion\|nf-validation\|nf-schema" nextflow.config
ls -R | head -80
```
Expected: `tree.txt` lists the repository; `nextflow_schema.json` about 44 kB; `config_params.txt` starts with `params {` and ends with `}`; no `MISSING` line for the first eight files. Record the commit of the tag from `tag_ref.json` in the gate report.

- [ ] **Step G0b: Fetch the test data on a compute node (also proves compute-node internet)**

Write `$GATE/jobs/g0_fetch.sh`:
```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=2G
#SBATCH -t 0:30:00
#SBATCH -p bcc
#SBATCH -o /net/bmc-lab3/data/bcc/yannvrb/rnasplice_gate/jobs/g0_fetch.log
B=https://raw.githubusercontent.com/nf-core/test-datasets/rnasplice
cd /net/bmc-lab3/data/bcc/yannvrb/rnasplice_gate/data || exit 1
for f in samplesheet/samplesheet.csv samplesheet/contrastsheet.csv reference/X.fa.gz reference/genes_chrX.gtf; do
  wget -q -c -O "$(basename "$f")" "$B/$f" || { echo "FAIL $f"; exit 1; }
done
head -1 samplesheet.csv
for u in $(awk -F',' 'NR > 1 {print $2; if ($3 != "") print $3}' samplesheet.csv); do
  wget -q -c -O "$(basename "$u")" "$u" || { echo "FAIL $u"; exit 1; }
done
ls -la
echo "G0B DONE"
```
Run: `sbatch $GATE/jobs/g0_fetch.sh`. Expected: `G0B DONE` in the log; 8 FASTQ files of 3.6-3.8 MB, `X.fa.gz` about 46 MB, `genes_chrX.gtf` about 26 MB. Record `COMPUTE_INTERNET` = `yes` (or `no`, then the gate stages everything from the login node and records it). If the header printed is not `sample,fastq_1,fastq_2,strandedness,condition`, record the real header; it becomes the fixture and Task 2 uses the real one.

- [ ] **Step G1: Run release 1.0.4 with its own test profile under the cluster's Nextflow**

Write `$GATE/nextflow.config` (this is also the planned skill config; Task 6 keeps the selectors that G8 confirms):
```nextflow
// nextflow.config — nf-core/rnasplice on SLURM + Singularity
profiles {
    slurm {
        process {
            executor = 'slurm'
            queue = 'bcc'
            cpus = 2
            memory = '8 GB'
            time = '4h'

            withName: '.*:STAR_GENOMEGENERATE' {
                cpus = 8
                memory = '64 GB'
                time = '8h'
            }
            withName: '.*:STAR_ALIGN' {
                cpus = 8
                memory = '48 GB'
                time = '8h'
            }
            withName: '.*:RMATS_PREP' {
                cpus = 4
                memory = '16 GB'
                time = '8h'
            }
            withName: '.*:RMATS_POST' {
                cpus = 8
                memory = '32 GB'
                time = '16h'
            }
            withName: '.*:DEXSEQ_COUNT' {
                cpus = 2
                memory = '8 GB'
                time = '8h'
            }
            withName: '.*:DEXSEQ_EXON' {
                cpus = 8
                memory = '32 GB'
                time = '8h'
            }
            withName: '.*:DEXSEQ_DTU' {
                cpus = 8
                memory = '32 GB'
                time = '8h'
            }
            withName: '.*:SALMON_QUANT' {
                cpus = 8
                memory = '16 GB'
                time = '4h'
            }
        }
        executor {
            queueSize = 10
            submitRateLimit = '10/1min'
            pollInterval = '30s'
        }
    }
    singularity {
        singularity {
            enabled = true
            autoMounts = true
        }
    }
}

process {
    resourceLimits = [ cpus: 16, memory: '64 GB', time: '24h' ]
}

timeline { enabled = true; overwrite = true; file = "${params.outdir}/pipeline_info/execution_timeline.html" }
report   { enabled = true; overwrite = true; file = "${params.outdir}/pipeline_info/execution_report.html"   }
trace    { enabled = true; overwrite = true; file = "${params.outdir}/pipeline_info/execution_trace.txt"; fields = 'task_id,hash,native_id,name,status,exit,cpus,memory,time,realtime,peak_rss' }
dag      { enabled = true; overwrite = true; file = "${params.outdir}/pipeline_info/pipeline_dag.svg"        }
```

Write `$GATE/jobs/g1.sh` (submitted from a normal login shell, so the job inherits `NXF_OFFLINE=TRUE` exactly as a user's job would):
```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=8G
#SBATCH -t 4:00:00
#SBATCH -p bcc
#SBATCH -o /net/bmc-lab3/data/bcc/yannvrb/rnasplice_gate/jobs/g1.log
cd /net/bmc-lab3/data/bcc/yannvrb/rnasplice_gate || exit 1
module add miniconda3/v4 || { echo "ERROR: miniconda3"; exit 1; }
source /home/software/conda/miniconda3/bin/condainit || { echo "ERROR: condainit"; exit 1; }
conda activate nf-env || { echo "ERROR: conda activate"; exit 1; }
module add singularity/3.10.4 || { echo "ERROR: singularity"; exit 1; }
echo "NXF_OFFLINE=${NXF_OFFLINE:-unset} NXF_SINGULARITY_CACHEDIR=${NXF_SINGULARITY_CACHEDIR:-unset}"
nextflow -version
printf 'outdir: "results/%s_g1_test_profile"\n' "$(date +%F)" > g1.yaml
nextflow run nf-core/rnasplice -r 1.0.4 -profile test,slurm,singularity -c nextflow.config -params-file g1.yaml
rc=$?; echo "G1A EXIT $rc"
if [ $rc -ne 0 ]; then
  NXF_OFFLINE=false nextflow pull nf-core/rnasplice -r 1.0.4; echo "PULL EXIT $?"
  nextflow run nf-core/rnasplice -r 1.0.4 -profile test,slurm,singularity -c nextflow.config -params-file g1.yaml
  echo "G1B EXIT $?"
fi
```
Run: `sbatch $GATE/jobs/g1.sh`; wait with a Monitor until-loop on `squeue -j <id>` (never a foreground sleep). Then evaluate:
```bash
cd $GATE
grep -nE "G1A EXIT|PULL EXIT|G1B EXIT|Pipeline completed successfully|ERROR|offline|Unknown|Invalid|not valid|nf-validation" jobs/g1.log | head -60
TR=$(ls -d results/*_g1_test_profile)/pipeline_info/execution_trace.txt
awk -F'\t' 'NR > 1 {print $5}' "$TR" | sort | uniq -c
awk -F'\t' 'NR > 1 && $5 != "COMPLETED" && $5 != "CACHED" {print $4, $5, $6}' "$TR"
```
Outcome rules (write the evidence lines into the gate report):
- `OFFLINE_PULL_NEEDED` = `yes` if `G1A EXIT` is non-zero with an offline/"cannot find" message and the run after the pull proceeds; else `no`.
- **G1 PASS**: `Pipeline completed successfully` and every task `COMPLETED` (a `FAILED` row is allowed only when a retry of the same task name completed). Set `GATE_OUTCOME` = `A`, `PIPELINE_REVISION` = `VERSION_TAG` = `1.0.4`, `NEXTFLOW_TESTED` = the printed version, `CONDA_ENV_TESTED` = `nf-env`, `NEXTFLOW_MAX_EXCL` = `none`. Continue at G3.
- **Only `MISO_*` tasks failed**: write `g1m.yaml` = `g1.yaml` plus `sashimi_plot: false`, rerun the job with `-params-file g1m.yaml` and outdir `..._g1m_no_miso`; if that completes, G1 PASS as above, and record the MISO failure (relevant to G7).
- **Launch failure** (an error before the first task: config/script compile error, plugin load failure, unsupported syntax under 26.04.6): 1.0.4 is incompatible; go to G2 and then the mid-gate checkpoint.
- **Any other task failure**: use superpowers:systematic-debugging on the task's `.command.err`/`.command.log`. Infrastructure causes (container pull, network, queue) are fixed and the run repeated; a pipeline bug in a module the skill offers is recorded as a checkpoint decision (possibly `STOP`).

- [ ] **Step G2 (only if G1 is a launch failure): run both fallback candidates, then a mid-gate checkpoint**

G2a, pinned dev commit:
```bash
cd $GATE
SHA=$(curl -sf https://api.github.com/repos/nf-core/rnasplice/commits/dev | grep -m1 '"sha"' | sed -E 's/.*"sha": *"([0-9a-f]{40})".*/\1/')
echo "$SHA"; [ ${#SHA} -eq 40 ] || echo "BAD SHA"
mkdir -p src/dev-${SHA:0:7}
```
Fetch the same file list into `src/dev-${SHA:0:7}/` (Step G0 loop with `1.0.4` replaced by `$SHA`), then a copy of `g1.sh` named `g2a.sh` with `-r 1.0.4` replaced by `-r $SHA` (literal value), outdir `..._g2a_dev`, log `jobs/g2a.log`.

G2b, release 1.0.4 with an older Nextflow in a new conda env (conda only on a compute node):
```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 1:00:00
#SBATCH -p bcc
#SBATCH -o /net/bmc-lab3/data/bcc/yannvrb/rnasplice_gate/jobs/g2b_env.log
module add miniconda3/v4 || { echo "ERROR: miniconda3"; exit 1; }
source /home/software/conda/miniconda3/bin/condainit || { echo "ERROR: condainit"; exit 1; }
conda create -y -n nf-24.04.4 -c conda-forge -c bioconda nextflow=24.04.4 || { echo "ERROR: conda create"; exit 1; }
conda activate nf-24.04.4 && nextflow -version
```
Then `g2b.sh` = `g1.sh` with `conda activate nf-24.04.4`, outdir `..._g2b_nf24`, log `jobs/g2b.log`. If it also fails at launch, repeat once with `nextflow=23.10.1` (env `nf-23.10.1`).

Mid-gate checkpoint (controller, before G3): write the G1/G2a/G2b evidence into `gate_report.md` section "Pin candidates" and ask the user (numbered): 1. pinned dev commit `<SHA>` with Nextflow 26.04.6 (`GATE_OUTCOME` B) · 2. release 1.0.4 with Nextflow `<version>` in env `<env>` (`GATE_OUTCOME` C) · 3. stop the project. Continue G3-G8 with the chosen candidate only (its `-r` value, its conda env, its source files). B: `NEXTFLOW_MIN` from the dev manifest, `NEXTFLOW_MAX_EXCL` = `none`. C: `NEXTFLOW_MAX_EXCL` = `26.04.0`.

- [ ] **Step G3: Typed params-file run with the settings the skill will write (all modules, local test data)**

Login node, light:
```bash
cd $GATE
RL=$(zcat data/ERR188383_1.fastq.gz | head -n 4000 | awk 'NR % 4 == 2 {print length($0)}' | sort -n | uniq -c | sort -k1,1nr -k2,2nr | head -1 | awk '{print $2}')
echo "read length $RL"
awk -F',' -v d="$GATE/data" 'NR == 1 {print; next} {n = split($2, a, "/"); m = split($3, b, "/"); print $1 "," d "/" a[n] "," (m ? d "/" b[m] : "") "," $4 "," $5}' data/samplesheet.csv > g3_samplesheet.csv
printf 'contrast,treatment,control\nGBR_vs_YRI,GBR,YRI\n' > g3_contrasts.csv
cat g3_samplesheet.csv
```
(The awk assumes the header `sample,fastq_1,fastq_2,strandedness,condition`; G0b recorded it. If it differs, rebuild the columns by name and record that.)

Write `g3.yaml` with the values below; the DEXSeq, edgeR and SUPPA values are those of the research report and must equal `src/<rev>/config_params.txt`: compare every line first, and where `config_params.txt` differs, use its value and record the difference in the gate report. Replace `<RL>` by `$RL`, `<DATE>` by `$(date +%F)`, `<GATE>` by the path:
```yaml
input: "<GATE>/g3_samplesheet.csv"
contrasts: "<GATE>/g3_contrasts.csv"
source: "fastq"
outdir: "results/<DATE>_g3_typed"
multiqc_title: "261001"
fasta: "<GATE>/data/X.fa.gz"
gtf: "<GATE>/data/genes_chrX.gtf"
gencode: false
save_reference: false
max_cpus: 4
max_memory: "16.GB"
max_time: "4.h"
aligner: "star"
pseudo_aligner: "salmon"
rmats: true
dexseq_exon: true
edger_exon: true
dexseq_dtu: true
suppa: true
sashimi_plot: false
rmats_read_len: <RL>
rmats_paired_stats: false
rmats_splice_diff_cutoff: 0.0001
rmats_novel_splice_site: false
rmats_min_intron_len: 50
rmats_max_exon_len: 500
alignment_quality: 10
aggregation: true
save_dexseq_annotation: false
save_dexseq_plot: true
n_dexseq_plot: 10
save_edger_plot: true
n_edger_plot: 10
dtu_txi: "dtuScaledTPM"
min_samps_gene_expr: 4
min_samps_feature_expr: 2
min_samps_feature_prop: 2
min_gene_expr: 10
min_feature_expr: 10
min_feature_prop: 0.1
suppa_per_local_event: true
suppa_per_isoform: true
generateevents_pool_genes: true
generateevents_event_type: "SE SS MX RI FL"
generateevents_boundary: "S"
generateevents_threshold: 10
generateevents_exon_length: 100
psiperevent_total_filter: 0
diffsplice_local_event: true
diffsplice_isoform: true
diffsplice_method: "empirical"
diffsplice_area: 1000
diffsplice_lower_bound: 0
diffsplice_gene_correction: true
diffsplice_paired: false
diffsplice_alpha: 0.05
diffsplice_median: false
diffsplice_tpm_threshold: 0
diffsplice_nan_threshold: 0
clusterevents_local_event: true
clusterevents_isoform: true
clusterevents_dpsithreshold: 0.05
clusterevents_eps: 0.05
clusterevents_metric: "euclidean"
clusterevents_min_pts: 20
clusterevents_method: "DBSCAN"
```
For candidate B: drop `max_cpus`/`max_memory`/`max_time` if they are not in the dev schema (record), and add `isoformswitchanalyzer: false` and `leafcutter: false` if those are dev schema parameters.

`jobs/g3.sh` = `g1.sh` without the test profile: `nextflow run nf-core/rnasplice -r <rev> -profile slurm,singularity -c nextflow.config -params-file g3.yaml` (no G1A/G1B branch; the pipeline is pulled by now), log `jobs/g3.log`. Evaluate as in G1, plus the typing test: `grep -nE "not valid|Invalid|should be|is not|Unrecognised|unrecognised|WARN" jobs/g3.log | head -40` and the parameter summary printed at the top of the log (numbers shown as numbers, `multiqc_title` as `261001`). Record all of it in `typing_test.txt` (the yaml, the validator lines, the verdict). A rejected key is fixed only in its YAML representation (quoting) per the schema type and rerun; a key rejected as unknown is removed and listed in `VALIDATOR_REJECTED`; a rejected module switch is a checkpoint decision (never silently dropped).

Salmon route, from the G3 trace:
```bash
TR=$(ls -d results/*_g3_typed)/pipeline_info/execution_trace.txt
awk -F'\t' 'NR > 1 {n = $4; sub(/ \(.*$/, "", n); print n}' "$TR" | grep -E 'SALMON|DRIMSEQ|DEXSEQ_DTU|STAGER|SUPPA' | sort | uniq -c
```
- `SALMON_ROUTE` = `pseudo_only` if DTU (`DEXSEQ_DTU`) and SUPPA tasks ran, each once per contrast, only in subworkflows without `STAR_SALMON` in their name. `PSEUDO_OFF_LINE` = `none`. G3 is then the reference run.
- Otherwise read `grep -n "pseudo_aligner" src/<rev>/workflows/rnasplice.nf src/<rev>/lib/WorkflowRnasplice.groovy`, and run `g3c` = `g3.yaml` with `aligner: "star_salmon"` and, as the pseudo-aligner line, the first of `pseudo_aligner: null` / `pseudo_aligner: false` that passes validation and leaves no task of the pseudo-alignment Salmon branch. If one works: `SALMON_ROUTE` = `star_salmon_only`, `PSEUDO_OFF_LINE` = that line; else `star_salmon_both` (both branches run; `PSEUDO_OFF_LINE` = `none`; the run with `aligner: "star_salmon"`, `pseudo_aligner: "salmon"` becomes the reference run). The fixture tree and process names come from the reference run.

- [ ] **Step G4: rMATS-only run with `reverse` declared (switch plumbing and strand mapping)**

`g4_samplesheet.csv` = `g3_samplesheet.csv` with column 4 set to `reverse`; `g4.yaml` = the reference yaml with `input` -> g4 sheet, outdir `..._g4_rmats_reverse`, and `dexseq_exon`, `edger_exon`, `dexseq_dtu`, `suppa` set to `false` (`rmats: true`, `sashimi_plot: false`). Run as G3. Expected and recorded:
```bash
TR=$(ls -d results/*_g4_rmats_reverse)/pipeline_info/execution_trace.txt
awk -F'\t' 'NR > 1 {print $4}' "$TR" | grep -cE 'DEXSEQ|EDGER|DRIMSEQ|STAGER|SUPPA|MISO'     # expected 0
```
and the `--libType` of every RMATS_PREP task (command recipe in G8) is `fr-firststrand`. A non-zero count means an explicit `false` does not switch a module off: that is a checkpoint decision.

- [ ] **Step G5: Genome-BAM input run**

BAMs: `find $(ls -d results/*_g3_typed) -name "*.bam" ! -name "*toTranscriptome*" | sort`. Use the coordinate-sorted per-sample genome BAMs published by the reference run (one per sample). If none are published, list `curl -sf "https://api.github.com/repos/nf-core/test-datasets/contents/testdata?ref=rnasplice" | grep -oE '"download_url": *"[^"]*Aligned\.out\.bam"'` and fetch those four with a `g0_fetch.sh`-style job. If the pipeline then needs sorted/indexed BAMs, sort and index them in an sbatch job with the cached samtools container (`singularity exec --bind $GATE /home/yannvrb/.singularity/cache/depot.galaxyproject.org-singularity-samtools-1.21--h50ea8bc_0.img samtools sort -o OUT IN` and `samtools index OUT`), and record that.

`g5_samplesheet.csv`: header = the revision's genome-BAM header (1.0.4: `sample,condition,genome_bam`; B: the header in `src/<rev>/assets/schema_input.json` or the usage docs, with `strandedness` = `reverse` written if that column exists); rows = sample, its condition from the G3 sheet, BAM path. `g5.yaml` = reference yaml with `source: "genome_bam"`, `input` -> g5 sheet, outdir `..._g5_bam`, `rmats`, `dexseq_exon`, `edger_exon` `true`, `dexseq_dtu`, `suppa` `false`. Run as G3. Record (G8 recipe) the `--libType` and `-t` of RMATS_PREP, the `-s` of the DEXSeq count task, the `-s` and `-p` of the featureCounts task, and whether `.bai` files were required (`BAM_NEEDS_BAI`). Values:
- run completes, columns absent: the four `BAM_*` values = what the command lines show (`fr-unstranded`/`paired`/`no`/`0` expected for 1.0.4 if it defaults to unstranded paired-end);
- B with `strandedness` column and `reverse` declared, command lines show `fr-firststrand` / `reverse` / `2`: all four = `from_sheet`;
- run fails for a pipeline reason (not input preparation): all four = `none` (BAM input becomes a checkpoint decision);
- `BAM_SHEET_HEADER` = the header used.

- [ ] **Step G6: Container availability**

```bash
CL=$(ls -d ~/.nextflow/assets/.repos/nf-core/rnasplice/clones/* 2>/dev/null | head -1); [ -n "$CL" ] || CL=$(ls -d ~/.nextflow/assets/nf-core/rnasplice)
grep -rhoE "https://depot.galaxyproject.org/singularity/[^ '\"]+" "$CL/modules" | sort -u > $GATE/container_urls.txt
while read -r u; do printf '%s %s %s\n' "$u" "$(curl -sI -o /dev/null -w '%{http_code}' "$u")" "$(date +%F)"; done < $GATE/container_urls.txt > $GATE/verified_urls.txt
grep -ciE 'rmats|dexseq|stager|suppa|edger|drimseq|misopy' $GATE/verified_urls.txt
ls -la "$NXF_SINGULARITY_CACHEDIR" | grep -iE 'rmats|dexseq|stager|suppa|misopy|drimseq|edger|subread|htseq'
```
Expected: HTTP 200 for every URL of the rMATS, DEXSeq, stageR, SUPPA and edgeR images (and the images present in the cache after G3). A non-200 for an offered module is a checkpoint decision.

- [ ] **Step G7: Default `miso_genes` in the test GTF**

```bash
for id in ENSG00000004961 ENSG00000005302 ENSG00000147403; do printf '%s\t%s\n' "$id" "$(grep -c "gene_id \"$id\"" $GATE/data/genes_chrX.gtf)"; done
awk -F'\t' 'NR > 1 && $4 ~ /MISO/ {print $4, $5}' $(ls -d $GATE/results/*_g1_*)/pipeline_info/execution_trace.txt
```
Record `MISO_GENES_IN_TEST_GTF` = `<n found>/3` and the MISO task statuses (informational: the skill always writes `sashimi_plot: false`).

- [ ] **Step G8: Extract the evidence (login node: awk/grep on result files)**

Command lines of named tasks, for each run (`G` = run dir name):
```bash
cd $GATE
for G in $(ls -d results/*_g3_typed results/*_g4_rmats_reverse results/*_g5_bam 2>/dev/null); do
  TR=$G/pipeline_info/execution_trace.txt
  awk -F'\t' 'NR > 1 && $4 ~ /RMATS_PREP|RMATS_POST|DEXSEQ_COUNT|FEATURECOUNTS|STAR_GENOMEGENERATE|STAR_ALIGN|SALMON_INDEX|CREATE_BAMLIST/ {print $2 "\t" $4}' "$TR" |
  while IFS=$'\t' read -r hash name; do
    d=$(ls -d work/${hash}* | head -1)
    echo "== $G | $name | $d"
    grep -E 'rmats|libType|--b1|--b2| -t |readLength|dexseq_count|featureCounts|sjdbOverhang|--genomeDir|salmon index' "$d/.command.sh"
  done
done > command_lines.txt
```
Then:
- process names: `for G in <every completed run>; do awk -F'\t' 'NR > 1 {n = $4; sub(/ \(.*$/, "", n); print n}' $G/pipeline_info/execution_trace.txt; done | sort -u > trace_process_names.txt`;
- requested resources for the gate report: `awk -F'\t' 'NR > 1 {n = $4; sub(/ \(.*$/, "", n); print n "\t" $7 "\t" $8 "\t" $9}' <reference trace> | sort -u` (confirms which `withName` selectors applied: STAR_ALIGN must show 8 CPUs and 48 GB);
- output tree of the reference run: `(cd <reference run dir> && find . -not -path './pipeline_info*' | sort) > output_tree.txt`;
- versions: `grep -iE 'star|salmon|rmats|dexseq|edger|suppa' <reference run>/pipeline_info/*versions*.yml`;
- `STAR_VERSION_GENOME`: `grep -h '^versionGenome' $(find $(ls -d work/<STAR_GENOMEGENERATE hash>*) -name genomeParameters.txt)` and its `sjdbOverhang` line (gate report);
- `SALMON_INDEX_VERSION`: `grep -h '"indexVersion"' $(find $(ls -d work/<SALMON_INDEX hash>*) -name versionInfo.json)`;
- `RMATS_B1_GROUP`: in the RMATS_POST (or RMATS_PREP) work dir of contrast `GBR_vs_YRI`, read the file passed as `--b1` and map its BAM names to samples and conditions: `treatment` if they are the GBR samples, `control` if YRI;
- `RMATS_BAMLIST_ORDER`: `grep -nE 'sort|toSortedList|collect|toList' src/<rev>/subworkflows/local/*rmats* src/<rev>/modules/local/*bamlist*`: `sorted_by_name` if the BAM list is sorted by sample name, `sheet` only if the code orders explicitly by samplesheet position, otherwise `unordered` (channel arrival order is not guaranteed);
- `DTU_FILTER_SCOPE`: read the DRIMSeq filter script and its subworkflow (`src/<rev>/bin/*drimseq*`, the DTU subworkflow): `per_contrast` if the filter runs once per contrast on that contrast's samples, else `all_samples`;
- `EDGER_DEU_FUNCTION`: `grep -oE 'diffSpliceDGE|diffSplice|exactTest|glmQLFTest|glmLRT|glmFit' src/<rev>/bin/*edger* | sort -u | tr '\n' ' '`.

- [ ] **Step G9: Write the fixtures and the gate report**

```bash
F=/net/bmc-lab3/data/bcc/yannvrb/Github_repos/bioinformatics-claude-skills/nfcore-rnasplice-setup/tests/fixtures
mkdir -p "$F"
cp $GATE/src/<rev>/nextflow_schema.json "$F/rnasplice_schema.json"
cp $GATE/src/<rev>/config_params.txt "$F/config_params.txt"
cp $GATE/data/samplesheet.csv "$F/test_samplesheet.csv"; cp $GATE/data/contrastsheet.csv "$F/test_contrastsheet.csv"
cp $GATE/trace_process_names.txt $GATE/output_tree.txt $GATE/command_lines.txt $GATE/verified_urls.txt $GATE/typing_test.txt "$F/"
```
Write `$F/gate_values.tsv` (tab-separated, one line per key of the table above, no header, no empty values; `none` where the table allows it) and `$F/gate_report.md` with these sections: `## GATE VERDICT` (outcome A/B/C/STOP in one line, then one paragraph); `## Environment` (Nextflow version, conda env, module lines, `NXF_*` values seen in the job, compute-node internet, pull behaviour); `## Runs` (one line per job: id, command, log path, result dir, verdict); `## Typing test`; `## Process names, labels and requested resources` (table from G8); `## Output tree` (summary, rMATS JC/JCEC file names, DEXSeq/edgeR/SUPPA/DTU result files); `## Strandedness and read type` (G4 and G5 command lines); `## Salmon route`; `## Containers`; `## miso_genes`; `## Facts that change the plan` (each with the task it affects); `## Decisions for the user`.

Verify the fixtures:
```bash
cd $F && for f in rnasplice_schema.json config_params.txt gate_values.tsv gate_report.md trace_process_names.txt output_tree.txt test_samplesheet.csv test_contrastsheet.csv command_lines.txt verified_urls.txt typing_test.txt; do [ -s "$f" ] || echo "EMPTY $f"; done
awk -F'\t' 'NF != 2 || $2 == "" {print "BAD LINE " NR ": " $0}' gate_values.tsv
grep -c "PIPELINE_REVISION" gate_values.tsv
```
Expected: no `EMPTY`, no `BAD LINE`, count 1.

- [ ] **Step G10: Commit the fixtures**

```bash
git status && git diff --stat && git diff
git add nfcore-rnasplice-setup/tests/fixtures
git commit -m "nfcore-rnasplice-setup: Task 0 verification gate fixtures"
```

- [ ] **Step G11: CONTROLLER CHECKPOINT (stop until the user answers)**

Show the user `gate_report.md` (verdict, environment, runs) and ask, numbered, for each decision that applies:
1. **Pin**: A — confirm release 1.0.4 with Nextflow `«NEXTFLOW_TESTED»` (env `nf-env`); B/C — confirm the candidate chosen at the mid-gate checkpoint, or stop.
2. **BAM input** (from G5): `fixed` values — "BAM input will be accepted only for `<libtype>`, `<readtype>`-end libraries; others are refused and pointed to FASTQ" (confirm, or drop BAM input); `none` — "BAM input does not work in this revision: drop it from the skill (spec change)" (confirm or stop); `from_sheet` — no decision.
3. **Salmon route**: `star_salmon_both` — confirm that DTU/SUPPA2 results will appear twice and the skill points the user to one of them (`star_salmon/`), or drop the two Salmon modules.
4. **Paired rMATS**: `unordered` — confirm that paired designs are not offered (`rmats_paired_stats` always `false`).
5. Any `VALIDATOR_REJECTED` key, a module switch that did not switch a module off (G4), a container not available (G6), or a pipeline bug in an offered module: one decision each, as written in the report.
6. The "Planner decisions (need user review)" list of this plan.
If the outcome is `STOP` (no candidate runs, or an offered core module, rMATS, is broken), the project stops: write `## GATE VERDICT: STOP` with the evidence and do not start Task 1. Record the answers in the ledger (`GATE: ...; user decisions: ...`) and, where a decision changes plan text (for example dropping BAM input), edit this plan's affected task text, commit it, and tell the user before Task 1.

---

### Task 1: Skeleton, Steps 0-3, checker scaffold, proof tools, root README row stub

**Files:**
- Create: `D/tests/cut_block.sh` (byte copy of `ase-pipeline/tests/synthetic/cut_block.sh`)
- Create: `D/tests/check_skill.sh`, `D/tests/prove_red.sh`, `D/tests/prove_mutations.sh`, `D/tests/mutations.tsv`
- Create: `D/nfcore-rnasplice-setup.md` (header, Steps 0-3)
- Modify: `README.md` (root; one table row)

**Interfaces:**
- Consumes: `D/tests/fixtures/*` from Task 0 (the checker exits 2 if one is missing or a gate key is absent).
- Produces: in `check_skill.sh`: `need "<text>"`, `forbid "<text>"`, `needr "<text>"` (README), `anchor_once "<anchor>"`, `gv KEY` (gate value), `$schema_names` (one parameter per line), `$config_kv` (`key<TAB>value` of the recorded config params block), `$CUT` (path of cut_block.sh), and the insertion point "before the final `[ $fail -eq 0 ]` line" for later task blocks. `CHECK_LIST_NEEDS=1` makes need/needr also print `NEED <checker line> <text>`. `prove_red.sh <BASE>` and `prove_mutations.sh <skill>`; `mutations.tsv` columns `id<TAB>runner<TAB>target<TAB>sed program<TAB>expected text`, runners `check strand bam validate readlen render submit helper`, targets `skill readme`, `@KEY@` = gate value. Skill placeholders `{CWD}`, `{WD_NAME}`, `{TODAY_YYMMDD}`, `{TODAY_ISO}`, `{USER_EMAIL}`, `{CONDA_ENV}`, `{VERSION}`, `{VERSION_TAG}`.

- [ ] **Step 1: Create the test tools**

`cp ase-pipeline/tests/synthetic/cut_block.sh nfcore-rnasplice-setup/tests/cut_block.sh` (byte copy; `cmp` must report no difference).

Create `D/tests/check_skill.sh`:
```bash
#!/bin/bash
# Static checks for nfcore-rnasplice-setup.md. Usage: check_skill.sh <skill.md>
# Uses only the recorded fixtures of the gated pipeline revision (tests/fixtures/); no network, no python.
# The README checked is the README.md next to the skill file. CHECK_LIST_NEEDS=1 also prints "NEED <line> <text>" per need.
set -u
SKILL=${1:?usage: check_skill.sh <skill.md>}
HERE=$(cd "$(dirname "$0")" && pwd); FIX="$HERE/fixtures"; CUT="$HERE/cut_block.sh"
README="$(dirname "$SKILL")/README.md"
fail=0
[ -s "$SKILL" ] || { echo "FAIL: skill file missing or empty: $SKILL"; exit 1; }
for f in rnasplice_schema.json config_params.txt gate_values.tsv trace_process_names.txt output_tree.txt test_samplesheet.csv test_contrastsheet.csv; do
  [ -s "$FIX/$f" ] || { echo "cannot read fixture $FIX/$f"; exit 2; }
done
GATE_KEYS="GATE_OUTCOME PIPELINE_REVISION VERSION_TAG NEXTFLOW_TESTED NEXTFLOW_MIN NEXTFLOW_MAX_EXCL CONDA_ENV_TESTED HAS_MAX_PARAMS HAS_RESOURCE_LIMITS SALMON_ROUTE PSEUDO_OFF_LINE BAM_SHEET_HEADER BAM_RMATS_LIBTYPE BAM_RMATS_READTYPE BAM_DEXSEQ_STRAND BAM_FC_STRAND BAM_NEEDS_BAI RMATS_BAMLIST_ORDER RMATS_B1_GROUP STAR_VERSION_GENOME SALMON_INDEX_VERSION DTU_FILTER_SCOPE TEST_CONTRAST TEST_READ_LENGTH EDGER_DEU_FUNCTION"
for k in $GATE_KEYS; do
  [ "$(awk -F'\t' -v k="$k" '$1 == k && $2 != ""' "$FIX/gate_values.tsv" | wc -l)" -eq 1 ] \
    || { echo "gate value $k must appear exactly once, non-empty, in fixtures/gate_values.tsv"; exit 2; }
done
gv() { awk -F'\t' -v k="$1" '$1 == k {print $2; exit}' "$FIX/gate_values.tsv"; }

need()   { [ -n "${CHECK_LIST_NEEDS:-}" ] && echo "NEED ${BASH_LINENO[0]} $1"; grep -qF -- "$1" "$SKILL" || { echo "FAIL: missing required text: $1"; fail=1; }; }
needr()  { [ -n "${CHECK_LIST_NEEDS:-}" ] && echo "NEED ${BASH_LINENO[0]} $1"; grep -qF -- "$1" "$README" 2>/dev/null || { echo "FAIL: README missing required text: $1"; fail=1; }; }
forbid() { ! grep -qF -- "$1" "$SKILL" || { echo "FAIL: forbidden text present: $1"; fail=1; }; }
anchor_once() { [ "$(grep -cF -- "$1" "$SKILL")" -eq 1 ] || { echo "FAIL: anchor must appear exactly once: $1"; fail=1; }; }

# Schema parameters = keys whose parent object key is "properties" (group names and structural keys are excluded).
schema_names=$(awk '
  { s = $0
    while (s != "") {
      if (match(s, /^[ \t]*"[^"]*"[ \t]*:[ \t]*\{/)) {
        k = s; sub(/^[ \t]*"/, "", k); sub(/".*$/, "", k)
        if (d > 0 && st[d] == "properties") print k
        st[++d] = k; s = substr(s, RSTART + RLENGTH); continue
      }
      c = substr(s, 1, 1); s = substr(s, 2)
      if (c == "\"") { while (s != "") { c2 = substr(s, 1, 1); s = substr(s, 2); if (c2 == "\\") { s = substr(s, 2); continue }; if (c2 == "\"") break } }
      else if (c == "{") st[++d] = ""
      else if (c == "}") d--
    } }' "$FIX/rnasplice_schema.json" | sort -u)
for p in input contrasts source outdir rmats rmats_read_len; do
  echo "$schema_names" | grep -qx "$p" || { echo "schema extraction broken: $p missing"; exit 2; }
done
for p in properties definitions input_output_options; do
  ! echo "$schema_names" | grep -qx "$p" || { echo "schema extraction broken: structural key $p listed"; exit 2; }
done
# Recorded pipeline config defaults: key<TAB>value with surrounding quotes and trailing // comments removed.
config_kv=$(awk -v q="'" '/^params[ \t]*\{/ {f = 1; next} f && /^\}/ {exit}
  f && /^[ \t]*[a-z_0-9]+[ \t]*=/ { k = $0; sub(/^[ \t]*/, "", k); sub(/[ \t]*=.*$/, "", k)
    v = $0; sub(/^[^=]*=[ \t]*/, "", v); sub(/[ \t]+\/\/.*$/, "", v); sub(/[ \t]+$/, "", v)
    gsub("^[\"" q "]|[\"" q "]$", "", v); print k "\t" v }' "$FIX/config_params.txt")
echo "$config_kv" | grep -q "^rmats	" || { echo "config extraction broken: rmats missing"; exit 2; }

# (a) no pipeline parameter as a --flag anywhere; only allowlisted tool flags (never a schema name)
allow="mail-type mail-user mem dependency parsable variable-read-length allow-clipping"
for a in $allow; do ! echo "$schema_names" | grep -qx -- "$a" || { echo "allowlist contains a pipeline parameter: $a"; exit 2; }; done
for flag in $(grep -oE '(^|[^A-Za-z0-9_-])--[A-Za-z_][A-Za-z_0-9-]*' "$SKILL" | sed -E 's/^[^-]*--//' | sort -u); do
  echo "$allow" | tr ' ' '\n' | grep -qx -- "$flag" || { echo "FAIL: --$flag is not an allowlisted tool flag (pipeline parameters go in the params file, never as --flags)"; fail=1; }
done
# (b) every key in a fenced yaml block is a parameter of the recorded schema
while IFS= read -r k; do
  echo "$schema_names" | grep -qx -- "$k" || { echo "FAIL: params key not in the recorded schema: $k"; fail=1; }
done < <(awk '/^```yaml$/ {f = 1; next} /^```$/ {f = 0} f' "$SKILL" | grep -oE '^[A-Za-z_0-9-]+:' | tr -d ':' | sort -u)
forbid "«"
forbid "»"
forbid "python3"
forbid "python -c"

# (c) content checks, appended per task -------------------------------
# --- Task 1 (Steps 0-3)
need "# nf-core/rnasplice Pipeline Setup Skill"
for n in 0 1 2 3; do need "## Step $n — "; done
need "{TODAY_ISO}"
need "api.github.com/repos/nf-core/rnasplice/releases/latest"
need "Do NOT use the \`gh\` CLI"
forbid "gh api"
need "verified against nf-core/rnasplice $(gv PIPELINE_REVISION) with Nextflow $(gv NEXTFLOW_TESTED)"
need "Nextflow $(gv NEXTFLOW_MIN) or newer"
[ "$(gv NEXTFLOW_MAX_EXCL)" = none ] || need "and older than $(gv NEXTFLOW_MAX_EXCL)"
need "Set \`{VERSION_TAG}\`"
need "Do not run \`nextflow\` here"
case "$(gv GATE_OUTCOME)" in
  B) need "that is a commit of the development branch (\`$(gv VERSION_TAG)\`)" ;;
  C) need "create_nextflow_env_$(gv NEXTFLOW_TESTED).sh" ;;
esac
# --- end Task 1

[ $fail -eq 0 ] && echo "PASS" || exit 1
```

Create `D/tests/prove_red.sh`:
```bash
#!/bin/bash
# Usage: prove_red.sh <base commit>   (run from anywhere inside the repository, before committing)
# Every need/needr line added to tests/check_skill.sh since <base> (working tree included) must fail on the skill and the
# README as they were at <base>. Prints "RED PASS (...)" or exits 1.
set -u
BASE=${1:?usage: prove_red.sh <base commit>}
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(git -C "$HERE" rev-parse --show-toplevel) || exit 1
D=nfcore-rnasplice-setup
T=$(mktemp -d) || exit 1
case "$T" in /tmp/tmp.*) ;; *) echo "unexpected temp dir: $T"; exit 1 ;; esac
trap 'rm -f "$T/nfcore-rnasplice-setup.md" "$T/README.md" "$T/out" "$T/list"; rmdir "$T"' EXIT
git -C "$ROOT" show "$BASE:$D/nfcore-rnasplice-setup.md" > "$T/nfcore-rnasplice-setup.md" 2>/dev/null || : > "$T/nfcore-rnasplice-setup.md"
git -C "$ROOT" show "$BASE:$D/README.md" > "$T/README.md" 2>/dev/null || : > "$T/README.md"
if [ ! -s "$T/nfcore-rnasplice-setup.md" ]; then echo "RED PASS (no skill file at $BASE: every need fails)"; exit 0; fi
bash "$HERE/check_skill.sh" "$T/nfcore-rnasplice-setup.md" > "$T/out" 2>&1
CHECK_LIST_NEEDS=1 bash "$HERE/check_skill.sh" "$ROOT/$D/nfcore-rnasplice-setup.md" 2>/dev/null | grep '^NEED ' > "$T/list"
added=$(git -C "$ROOT" diff -U0 "$BASE" -- "$D/tests/check_skill.sh" \
  | awk '/^@@/ { split($3, p, ","); s = substr(p[1], 2); n = (p[2] == "") ? 1 : p[2]; for (i = 0; i < n; i++) print s + i }')
n=0; bad=0
while IFS= read -r line; do
  ln=${line#NEED }; ln=${ln%% *}; text=${line#NEED $ln }
  echo "$added" | grep -qx -- "$ln" || continue
  n=$((n + 1))
  grep -qxF -- "FAIL: missing required text: $text" "$T/out" || grep -qxF -- "FAIL: README missing required text: $text" "$T/out" \
    || { echo "NOT RED on $BASE (checker line $ln): $text"; bad=1; }
done < "$T/list"
[ $bad -eq 0 ] && echo "RED PASS ($n new needs fail on $BASE)" || exit 1
```

Create `D/tests/prove_mutations.sh`:
```bash
#!/bin/bash
# Usage: prove_mutations.sh <skill.md>
# Each row of mutations.tsv (id<TAB>runner<TAB>target<TAB>sed program<TAB>expected text) is applied to a copy of the skill and
# its README (target: skill or readme); @KEY@ in the program and the text is replaced by the gate value KEY. A row passes when
# the sed program changed the target, the runner exits non-zero, and the runner printed the expected text.
# Runners: check strand bam validate readlen render submit helper (see runner() below).
# Prints "MUTATIONS PASS (<n>)" or exits 1.
set -u
SKILL=${1:?usage: prove_mutations.sh <skill.md>}
case "$SKILL" in /*) ;; *) SKILL="$PWD/$SKILL" ;; esac
HERE=$(cd "$(dirname "$0")" && pwd); README="$(dirname "$SKILL")/README.md"
T=$(mktemp -d) || exit 1
case "$T" in /tmp/tmp.*) ;; *) echo "unexpected temp dir: $T"; exit 1 ;; esac
trap 'rm -f "$T/nfcore-rnasplice-setup.md" "$T/README.md" "$T/out"; rmdir "$T"' EXIT
gv() { awk -F'\t' -v k="$1" '$1 == k {print $2; exit}' "$HERE/fixtures/gate_values.tsv"; }
subst() { local s=$1 k; while [[ $s =~ @([A-Z_0-9]+)@ ]]; do k=${BASH_REMATCH[1]}; s=${s//@$k@/$(gv "$k")}; done; printf '%s' "$s"; }
runner() { case "$1" in
  check) echo check_skill.sh ;; strand) echo test_strandedness.sh ;; bam) echo test_bam_policy.sh ;;
  validate) echo test_validate_sheets.sh ;; readlen) echo test_read_length.sh ;; render) echo test_render_params.sh ;;
  submit) echo dry_run_submit.sh ;; helper) echo dry_run_helpers.sh ;; *) echo "" ;; esac; }
for r in $(awk -F'\t' '$1 !~ /^#/ && NF >= 5 {print $2}' "$HERE/mutations.tsv" | sort -u); do
  s=$(runner "$r"); [ -n "$s" ] || { echo "unknown runner: $r"; exit 1; }
  bash "$HERE/$s" "$SKILL" > "$T/out" 2>&1 || { echo "$s fails on the unmutated skill:"; tail -20 "$T/out"; exit 1; }
done
n=0; bad=0
while IFS=$'\t' read -r id run target prog expect; do
  case "$id" in ''|'#'*) continue ;; esac
  n=$((n + 1)); prog=$(subst "$prog"); expect=$(subst "$expect")
  cp "$SKILL" "$T/nfcore-rnasplice-setup.md"; cp "$README" "$T/README.md" 2>/dev/null || : > "$T/README.md"
  case "$target" in
    skill) f="$T/nfcore-rnasplice-setup.md"; orig=$SKILL ;;
    readme) f="$T/README.md"; orig=$README ;;
    *) echo "MUTATION $id: bad target $target"; bad=1; continue ;;
  esac
  sed -i -e "$prog" "$f" || { echo "MUTATION $id: sed error"; bad=1; continue; }
  if cmp -s "$orig" "$f"; then echo "MUTATION $id: the sed program changed nothing"; bad=1; continue; fi
  s=$(runner "$run")
  if bash "$HERE/$s" "$T/nfcore-rnasplice-setup.md" > "$T/out" 2>&1; then echo "MUTATION $id: $s passed (it must fail)"; bad=1; continue; fi
  grep -qF -- "$expect" "$T/out" || { echo "MUTATION $id: $s failed without printing '$expect':"; tail -20 "$T/out"; bad=1; }
done < "$HERE/mutations.tsv"
[ $bad -eq 0 ] && echo "MUTATIONS PASS ($n)" || exit 1
```

Create `D/tests/mutations.tsv` with the header comment line only (tab-separated):
```
# id	runner	target	sed program	expected text
```

- [ ] **Step 2: Run the checker to verify it fails**

Run: `bash nfcore-rnasplice-setup/tests/check_skill.sh nfcore-rnasplice-setup/nfcore-rnasplice-setup.md`
Expected: `FAIL: skill file missing or empty: ...`, exit 1. Also run it against a scratch copy of `nfcore-rnavar-setup/nfcore-rnavar-setup.md` placed in a `mktemp -d` directory: expected exit 1 with `FAIL: missing required text: # nf-core/rnasplice Pipeline Setup Skill` and several `--flag ... not an allowlisted tool flag` lines (the rnavar skill names STAR flags), which proves the flag rule and the needs are live.

- [ ] **Step 3: Write the skill header and Steps 0-3**

Create `D/nfcore-rnasplice-setup.md` with this text, choosing the [A/C] or [B] variants and substituting every `«KEY»` from `gate_values.tsv`:

````markdown
# nf-core/rnasplice Pipeline Setup Skill

You are helping the user set up and submit nf-core/rnasplice, a pipeline for differential alternative splicing of bulk RNA-seq, on an HPC cluster using SLURM and Singularity. rnasplice aligns the reads with STAR, quantifies transcripts with Salmon, and runs the analyses the user chooses: rMATS (differential splicing events: skipped exons, alternative 5' and 3' splice sites, mutually exclusive exons, retained introns), SUPPA2 (event and isoform PSI from transcript abundance), DEXSeq and edgeR differential exon usage, and DEXSeq differential transcript usage (DRIMSeq filter, DEXSeq, stageR). Walk through each step below in order, ask only what you need, and perform automated steps silently. This skill configures and submits the run; it does not analyse the results.

---

## Step 0 — Working directory

Before anything else, run `pwd` to record the current working directory (`{CWD}`). The samplesheet, the contrasts sheet, the params file, `nextflow.config`, the scripts and the output directory are written here. Tell the user: "Working directory: `{CWD}`"

Also derive:
- `{WD_NAME}` = basename of `{CWD}`
- `{TODAY_YYMMDD}` = today's date as `YYMMDD`
- `{TODAY_ISO}` = today's date as `YYYY-MM-DD` (used for the output directory `results/{TODAY_ISO}_{WD_NAME}`)

---

## Step 1 — Email address

Ask the user: "What email address should SLURM use for pipeline notifications (completion/failure)?"

Do **not** pre-suggest or pre-fill any email address. Wait for the answer before Step 2. Store `{USER_EMAIL}`.

---

## Step 2 — Conda environment

Ask the user: "Which conda environment should be activated for Nextflow?"

Do **not** pre-suggest any environment name. Ask this as a separate question after Step 1 — never combine Steps 1 and 2 in one message. Store `{CONDA_ENV}`.

[A/B] Tell the user that this pipeline revision needs Nextflow «NEXTFLOW_MIN» or newer (verified with Nextflow «NEXTFLOW_TESTED»).
[C] Tell the user that this pipeline revision needs Nextflow «NEXTFLOW_MIN» or newer and older than «NEXTFLOW_MAX_EXCL» (verified with Nextflow «NEXTFLOW_TESTED»; newer Nextflow versions cannot launch it). If the user has no environment with such a Nextflow, offer to generate `create_nextflow_env_«NEXTFLOW_TESTED».sh` (Step 11), which creates the environment `nf-«NEXTFLOW_TESTED»` on a compute node; then `{CONDA_ENV}` = `nf-«NEXTFLOW_TESTED»`.

Do not run `nextflow` here to check the version: Nextflow runs only on compute nodes. The submission script (Step 11) checks the version when the job starts and stops with a clear message if it is outside the verified range.

---

## Step 3 — nf-core/rnasplice version

[A/C] This skill was verified against nf-core/rnasplice «PIPELINE_REVISION» with Nextflow «NEXTFLOW_TESTED» (see the README, Validation status).
[B] This skill was verified against nf-core/rnasplice «PIPELINE_REVISION» with Nextflow «NEXTFLOW_TESTED»; that is a commit of the development branch (`«VERSION_TAG»`), pinned because release 1.0.4 does not run under that Nextflow (see the README).

Fetch the latest release with the `WebFetch` tool against:
```
https://api.github.com/repos/nf-core/rnasplice/releases/latest
```
Extract `tag_name` (strip a leading `v`) as `{TAG}`.
[A/C]
- If `{TAG}` is «PIPELINE_REVISION»: set `{VERSION}` = «PIPELINE_REVISION» and tell the user.
- If `{TAG}` is newer: tell the user, fetch `https://raw.githubusercontent.com/nf-core/rnasplice/{TAG}/nextflow_schema.json` with `WebFetch`, and list every key this skill writes (the Step 8 module block and the Step 11 template) that is not a parameter there. Ask (numbered): 1. Use «PIPELINE_REVISION» (verified; default) · 2. Use {TAG} (not verified by this skill; offer it only when no key is missing). Set `{VERSION}` to the choice.
[B]
- Set `{VERSION}` = «PIPELINE_REVISION» (`nextflow run -r` accepts a commit).
- If `{TAG}` is newer than 1.0.4: tell the user, fetch its schema as above, list the keys this skill writes that are missing from it, and ask (numbered): 1. Use the pinned commit (verified; default) · 2. Use release {TAG} (not verified by this skill; offer it only when no key is missing). If the user picks the release, `{VERSION}` = {TAG}.

Set `{VERSION_TAG}` = `{VERSION}` when it is a release tag, or `dev-` followed by the first 7 characters when it is a commit; it is used in file names.

**Important:** Do NOT use the `gh` CLI — it is not installed on this HPC cluster. Always use `WebFetch`.
````

- [ ] **Step 4: Run the checker to verify it passes**

Run: `bash nfcore-rnasplice-setup/tests/check_skill.sh nfcore-rnasplice-setup/nfcore-rnasplice-setup.md`
Expected: `PASS`.

- [ ] **Step 5: Add the Task 1 mutation rows and prove them**

Append to `D/tests/mutations.tsv` (fields separated by one tab):
```
T1-flag	check	skill	s/^Extract `tag_name`/Extract --rmats_read_len `tag_name`/	--rmats_read_len is not an allowlisted tool flag
T1-gh	check	skill	s/Always use `WebFetch`\./Always use `gh api`./	forbidden text present: gh api
T1-yaml	check	skill	$s/$/\n```yaml\nbogus_param: 1\n```/	params key not in the recorded schema: bogus_param
T1-group	check	skill	$s/$/\n```yaml\ninput_output_options: 1\n```/	params key not in the recorded schema: input_output_options
T1-rev	check	skill	s/nf-core\/rnasplice @PIPELINE_REVISION@ with Nextflow/nf-core\/rnasplice 0.0.0 with Nextflow/	missing required text: verified against nf-core/rnasplice @PIPELINE_REVISION@ with Nextflow @NEXTFLOW_TESTED@
T1-nfmin	check	skill	s/Nextflow @NEXTFLOW_MIN@ or newer/Nextflow 0.0.1 or newer/g	missing required text: Nextflow @NEXTFLOW_MIN@ or newer
T1-guil	check	skill	s/^# nf-core\/rnasplice Pipeline Setup Skill$/&\n«X»/	forbidden text present: «
```
Run: `bash nfcore-rnasplice-setup/tests/prove_mutations.sh nfcore-rnasplice-setup/nfcore-rnasplice-setup.md`
Expected: `MUTATIONS PASS (7)`.

- [ ] **Step 6: Root README row stub**

In `README.md` (root), insert after the `nf-core/rnavar setup` row:
```
| nf-core/rnasplice setup | `/nfcore-rnasplice-setup` | (in development, not yet validated) Interactive setup wizard for nf-core/rnasplice differential alternative splicing (rMATS, SUPPA2, DEXSeq/edgeR exon usage, DEXSeq DTU) on SLURM + Singularity |
```

- [ ] **Step 7: RED proof and commit**

```bash
bash nfcore-rnasplice-setup/tests/prove_red.sh <BASE>
bash nfcore-rnasplice-setup/tests/check_skill.sh nfcore-rnasplice-setup/nfcore-rnasplice-setup.md
git status && git diff --stat && git diff
git add nfcore-rnasplice-setup README.md
git commit -m "nfcore-rnasplice-setup: skeleton (Steps 0-3), fixture-based checker, proof tools"
```
Expected: `RED PASS (no skill file at <BASE>: every need fails)`, `PASS`.

---

### Task 2: Input source, FASTQ rows, strandedness, BAM rows (Step 4)

**Files:**
- Create: `D/tests/test_env.sh`, `D/tests/test_strandedness.sh`, `D/tests/test_bam_policy.sh`
- Modify: `D/tests/check_skill.sh` (Task 2 block), `D/tests/mutations.tsv` (T2 rows)
- Modify: `D/nfcore-rnasplice-setup.md` (append Step 4)

**Interfaces:**
- Consumes: `need`, `forbid`, `anchor_once`, `gv`, `$FIX` (Task 1); `{CWD}`, `{TODAY_YYMMDD}`.
- Produces: skill placeholders `{SOURCE}` (`fastq`|`genome_bam`), `{FASTQ_DIR}`, `{LAYOUT}` (`paired`|`single`), `{SEQ_DATE}`, `{STRANDEDNESS}` (`unstranded`|`forward`|`reverse`), `{RNASEQ_OUTDIR}`; in-memory rows (sample, fastq_1, fastq_2 or genome_bam) used by Step 5; skill functions `infer_strandedness <file>` (prints `forward`|`reverse`|`unstranded`|`unclear`) and `bam_input_allowed <strandedness> <paired|single>` (prints `ALLOWED: ...` rc 0 or `REFUSED: ...` rc 1); anchors `**Strandedness from an existing nf-core/rnaseq run.**`, `**BAM input rule.**`. `D/tests/test_env.sh` (sourced by every test that runs extracted code): `$TEST_TMP`, `$STUBS`, `forbidden_ran`.

- [ ] **Step 1: Write the test environment and the two tests (they fail: no blocks yet)**

Create `D/tests/test_env.sh`:
```bash
# Sourced by the test scripts (source "$HERE/test_env.sh" "$@"): re-runs the calling script under env -i with a clean PATH,
# creates a stub directory in which every command that must never run during a test (Nextflow, conda, Java, Python, R,
# Singularity, SLURM, Lmod, STAR, Salmon) is a stub that fails and leaves a trace, and type-checks the stubs.
# Provides: $TEST_TMP (a mktemp -d directory, removed at exit), $STUBS, forbidden_ran (true when a forbidden stub ran).
if [ "${TEST_CLEAN_ENV:-}" != 1 ]; then
  exec env -i TEST_CLEAN_ENV=1 HOME="${HOME:-/nonexistent}" PATH=/usr/bin:/bin /bin/bash --noprofile --norc "$0" "$@"
fi
TEST_TMP=$(mktemp -d) || exit 1
case "$TEST_TMP" in /tmp/tmp.?*) ;; *) echo "unexpected temp dir: $TEST_TMP"; exit 1 ;; esac
trap 'rm -rf -- "$TEST_TMP"' EXIT
STUBS="$TEST_TMP/stubs"; mkdir -p "$STUBS" || exit 1
FORBIDDEN_CMDS="nextflow conda java python python3 R Rscript singularity sbatch srun squeue module STAR salmon"
for c in $FORBIDDEN_CMDS; do
  printf '#!/bin/bash\necho "FORBIDDEN %s $*" >> "%s/forbidden"\nexit 99\n' "$c" "$TEST_TMP" > "$STUBS/$c"; chmod +x "$STUBS/$c"
done
export PATH="$STUBS:/usr/bin:/bin"
for c in $FORBIDDEN_CMDS; do [ "$(type -P "$c")" = "$STUBS/$c" ] || { echo "stub $c does not resolve to $STUBS/$c"; exit 1; }; done
forbidden_ran() { if [ -s "$TEST_TMP/forbidden" ]; then echo "FAIL: a forbidden command ran:"; cat "$TEST_TMP/forbidden"; return 0; fi; return 1; }
```

Create `D/tests/test_strandedness.sh`:
```bash
#!/bin/bash
# Usage: test_strandedness.sh <skill.md>
# Runs infer_strandedness from the skill (block after "**Strandedness from an existing nf-core/rnaseq run.**") on RSeQC
# infer_experiment.txt files written here. Prints "STRANDEDNESS PASS" or exits 1.
set -u
source "$(cd "$(dirname "$0")" && pwd)/test_env.sh" "$@"
SKILL=${1:?usage: test_strandedness.sh <skill.md>}
HERE=$(cd "$(dirname "$0")" && pwd)
CODE=$(bash "$HERE/cut_block.sh" "$SKILL" "**Strandedness from an existing nf-core/rnaseq run.**") || { echo "FAIL: strandedness block not found"; exit 1; }
eval "$CODE"
declare -F infer_strandedness >/dev/null || { echo "FAIL: infer_strandedness is not defined by the block"; exit 1; }
mk() { # mk <name> <PE|SE> <fraction explained by sense pattern> <fraction explained by antisense pattern>
  if [ "$2" = PE ]; then
    printf '\n\nThis is PairEnd Data\nFraction of reads failed to determine: 0.0100\nFraction of reads explained by "1++,1--,2+-,2-+": %s\nFraction of reads explained by "1+-,1-+,2++,2--": %s\n' "$3" "$4" > "$TEST_TMP/$1"
  else
    printf '\n\nThis is SingleEnd Data\nFraction of reads failed to determine: 0.0100\nFraction of reads explained by "++,--": %s\nFraction of reads explained by "+-,-+": %s\n' "$3" "$4" > "$TEST_TMP/$1"
  fi; }
fail=0
check() { local got; got=$(infer_strandedness "$TEST_TMP/$1"); [ "$got" = "$2" ] || { echo "FAIL: case $1: expected $2, got '$got'"; fail=1; }; }
mk pe_forward PE 0.9512 0.0388;    check pe_forward forward
mk pe_reverse PE 0.0210 0.9690;    check pe_reverse reverse
mk pe_unstranded PE 0.4950 0.4950; check pe_unstranded unstranded
mk se_reverse SE 0.0500 0.9400;    check se_reverse reverse
mk se_forward SE 0.8800 0.1100;    check se_forward forward
mk pe_ambiguous PE 0.7000 0.2900;  check pe_ambiguous unclear
: > "$TEST_TMP/empty";             check empty unclear
printf 'Fraction of reads explained by "1++,1--,2+-,2-+": 0.95\n' > "$TEST_TMP/one_line"; check one_line unclear
forbidden_ran && fail=1
[ $fail -eq 0 ] && echo "STRANDEDNESS PASS" || exit 1
```

Create `D/tests/test_bam_policy.sh`:
```bash
#!/bin/bash
# Usage: test_bam_policy.sh <skill.md>
# Runs bam_input_allowed from the skill (block after "**BAM input rule.**") for every strandedness x read type and compares
# each answer with the expectation computed here from fixtures/gate_values.tsv. Prints "BAM POLICY PASS" or exits 1.
set -u
source "$(cd "$(dirname "$0")" && pwd)/test_env.sh" "$@"
SKILL=${1:?usage: test_bam_policy.sh <skill.md>}
HERE=$(cd "$(dirname "$0")" && pwd)
gv() { awk -F'\t' -v k="$1" '$1 == k {print $2; exit}' "$HERE/fixtures/gate_values.tsv"; }
CODE=$(bash "$HERE/cut_block.sh" "$SKILL" "**BAM input rule.**") || { echo "FAIL: BAM input rule block not found"; exit 1; }
eval "$CODE"
declare -F bam_input_allowed >/dev/null || { echo "FAIL: bam_input_allowed is not defined by the block"; exit 1; }
G_LT=$(gv BAM_RMATS_LIBTYPE); G_RT=$(gv BAM_RMATS_READTYPE); G_DS=$(gv BAM_DEXSEQ_STRAND); G_FC=$(gv BAM_FC_STRAND)
fail=0
for lib in unstranded forward reverse; do
  case $lib in unstranded) lt=fr-unstranded ds=no fc=0 ;; forward) lt=fr-secondstrand ds=yes fc=1 ;; reverse) lt=fr-firststrand ds=reverse fc=2 ;; esac
  for lay in paired single; do
    if [ "$G_LT" = none ]; then exp=REFUSED
    elif [ "$G_LT" = from_sheet ]; then exp=ALLOWED
    elif [ "$lt" = "$G_LT" ] && [ "$lay" = "$G_RT" ] && { [ "$G_DS" = n/a ] || [ "$ds" = "$G_DS" ]; } && { [ "$G_FC" = n/a ] || [ "$fc" = "$G_FC" ]; }; then exp=ALLOWED
    else exp=REFUSED; fi
    out=$(bam_input_allowed "$lib" "$lay"); rc=$?
    if [ "$exp" = ALLOWED ]; then
      [ $rc -eq 0 ] && [[ $out == ALLOWED* ]] || { echo "FAIL: case $lib/$lay: expected ALLOWED, got rc=$rc '$out'"; fail=1; }
    else
      [ $rc -ne 0 ] && [[ $out == REFUSED* ]] && [[ $out == *"Start from FASTQ"* ]] || { echo "FAIL: case $lib/$lay: expected a REFUSED answer pointing to FASTQ, got rc=$rc '$out'"; fail=1; }
    fi
  done
done
out=$(bam_input_allowed auto paired); rc=$?; [ $rc -ne 0 ] && [[ $out == REFUSED* ]] || { echo "FAIL: case invalid_strandedness"; fail=1; }
out=$(bam_input_allowed unstranded both); rc=$?; [ $rc -ne 0 ] && [[ $out == REFUSED* ]] || { echo "FAIL: case invalid_layout"; fail=1; }
forbidden_ran && fail=1
[ $fail -eq 0 ] && echo "BAM POLICY PASS" || exit 1
```

- [ ] **Step 2: Append the Task 2 checker block and run everything to see it fail**

Insert before the final `[ $fail -eq 0 ]` line of `check_skill.sh`:
```bash
# --- Task 2 (Step 4)
need "## Step 4 — Input files and samplesheet rows"
need "sample,fastq_1,fastq_2,strandedness,condition"
need "$(head -n1 "$FIX/test_samplesheet.csv" | tr -d '\r')"
need "Strandedness is asked, never guessed."
need "1. unstranded · 2. forward · 3. reverse · 4. I don't know"
need "rnasplice has no \`auto\` strandedness"
forbid "strandedness to \`auto\`"
forbid "strandedness: auto"
need "Never continue with \"I don't know\""
need "rMATS needs one strandedness and one read type for all samples"
need "if a name then starts with a digit, prefix \`S\`"
need "1. the same sample (lanes or technical replicates; the pipeline merges their reads)"
need "2. different samples — rename them"
anchor_once "**Strandedness from an existing nf-core/rnaseq run.**"
anchor_once "**BAM input rule.**"
need "infer_strandedness() {"
need "bam_input_allowed() {"
need "RMATS_LIBTYPE=\"$(gv BAM_RMATS_LIBTYPE)\""
need "RMATS_READTYPE=\"$(gv BAM_RMATS_READTYPE)\""
need "DEXSEQ_STRAND=\"$(gv BAM_DEXSEQ_STRAND)\""
need "FC_STRAND=\"$(gv BAM_FC_STRAND)\""
need "Start from FASTQ"
need "DTU and SUPPA2 need Salmon quantification from reads"
[ "$(gv BAM_RMATS_LIBTYPE)" = none ] || need "$(gv BAM_SHEET_HEADER)"
# --- end Task 2
```
Run: `bash nfcore-rnasplice-setup/tests/check_skill.sh nfcore-rnasplice-setup/nfcore-rnasplice-setup.md`, `bash nfcore-rnasplice-setup/tests/test_strandedness.sh nfcore-rnasplice-setup/nfcore-rnasplice-setup.md`, `bash nfcore-rnasplice-setup/tests/test_bam_policy.sh nfcore-rnasplice-setup/nfcore-rnasplice-setup.md`.
Expected: the checker prints `FAIL: missing required text:` lines for the Task 2 strings and two `FAIL: anchor must appear exactly once` lines; both tests print `FAIL: ... block not found`; all exit 1.

- [ ] **Step 3: Append Step 4 to the skill**

Choose the variant by `BAM_RMATS_LIBTYPE`: **[BAM]** when it is not `none`, **[NO-BAM]** when it is `none` (the user confirmed at the checkpoint that BAM input is dropped). Substitute every `«KEY»`.

````markdown

---

## Step 4 — Input files and samplesheet rows

**Scan the working directory first:**
```bash
find {CWD} \( -name "*.fastq.gz" -o -name "*.fq.gz" \) | head -50
find {CWD} -name "*.bam" | head -20
```
[BAM] Ask (numbered), mentioning what the scan found: "Which input do you want to start from?" 1. FASTQ files (default; the pipeline aligns them with STAR and quantifies them with Salmon) · 2. Genome BAM files from a splice-aware aligner (for example the STAR BAMs of an nf-core/rnaseq run). Store `{SOURCE}` = `fastq` or `genome_bam`.
[NO-BAM] The input is FASTQ: `{SOURCE}` = `fastq`. If the user only has BAM files, explain that genome-BAM input does not work with this pipeline revision (Step 4b) and ask for the FASTQ files.

### Step 4a — FASTQ input

- FASTQ files found: list the unique directories and ask "I found FASTQ files in: `{FOUND_DIRS}`. Use this directory, or specify another?" Nothing found: ask for the full path. Store `{FASTQ_DIR}`.
- **Paired-end detection:** filenames containing `_R1_`/`_R2_`, `_1.fastq.gz`/`_2.fastq.gz`, `_1.fq.gz`/`_2.fq.gz`, or `_1_sequence`/`_2_sequence` → paired-end; otherwise single-end. Store `{LAYOUT}` = `paired` or `single`. rMATS needs one strandedness and one read type for all samples (the pipeline stops with "Cannot run rMats with mixed single and paired end samples"), so if the directory mixes paired-end and single-end files, stop and ask which files to use.
- **Sequencing date:** a leading 6-digit `YYMMDD` filename prefix → `{SEQ_DATE}` (fall back to `{TODAY_YYMMDD}`).
- **Path style:** if `{FASTQ_DIR}` is inside `{CWD}`, use paths relative to `{CWD}`; otherwise absolute paths.
- **Rows:** find each R1 and pair it with its R2 by substituting `_R1_`→`_R2_`, `_1.`→`_2.` or `_1_sequence`→`_2_sequence`; if an R2 is missing, warn and stop (one read type for all samples). Single-end rows leave `fastq_2` empty. Sample name = filename up to `_S\d+`, `_R1`, `_1.f` or `_1_sequence`. FASTQ files must end in `.fastq.gz` or `.fq.gz`; the pipeline accepts nothing else.

**Sample-name sanitisation (always, before showing the user):** replace `-`, spaces, `/`, `(`, `)` and every other character outside `A-Za-z0-9_` with `_`; if a name then starts with a digit, prefix `S` (R renames columns that start with a digit, which breaks the matching of sample names in the DEXSeq and edgeR steps). Note every substitution in the preview. Then check uniqueness. On a collision warn: "⚠️ Name collision '{NAME}': in nf-core/rnasplice, rows with the same sample name are one sample, and their reads are merged before alignment." and ask (numbered): 1. the same sample (lanes or technical replicates; the pipeline merges their reads) · 2. different samples — rename them (go to custom naming). Lane files (`X_S1_L001_R1_001`, `X_S1_L002_R1_001`) collide by design; option 1 is expected for them.

**Review and naming — order is mandatory:**
1. Show all rows as a table (sample, fastq_1, fastq_2).
2. Ask (numbered): 1. Use auto-generated names · 2. Provide custom names. For custom names: show numbered auto names next to the filenames, ask for a plain-language description, build the mapping, show it as an auto → new table and ask "Does this mapping look correct?" Custom names follow the same rule (letters, digits and `_`, starting with a letter); duplicates only where the user chose option 1 above.

The rows stay in memory: Step 5 adds the condition column (the samplesheet columns are `sample,fastq_1,fastq_2,strandedness,condition`) and writes the file.

**Strandedness is asked, never guessed.** rnasplice has no `auto` strandedness, and rMATS needs one strandedness and one read type for all samples. Ask (numbered): "Which strandedness does the library have?" 1. unstranded · 2. forward · 3. reverse · 4. I don't know. Explain: *reverse* = read 1 comes from the strand opposite to the transcript (dUTP libraries such as Illumina TruSeq Stranded and NEBNext Ultra II Directional — the most common stranded libraries); *forward* = read 1 comes from the transcript strand; *unstranded* = no strand information (for example TruSeq non-stranded). A wrong value makes rMATS discard or misassign junction reads without any error. If the user answers 4, or has an nf-core/rnaseq output directory for the same samples, use the block below. Never continue with "I don't know". Store `{STRANDEDNESS}` (one value for every row).

**Strandedness from an existing nf-core/rnaseq run.** Ask for that run's output directory (`{RNASEQ_OUTDIR}`), list its RSeQC files with `find {RNASEQ_OUTDIR} -name "*.infer_experiment.txt" | head -50`, define the function below in the Bash tool and run it on each file (small text files; light work on the login node). Show a table file → result.
```bash
infer_strandedness() {
  # usage: infer_strandedness <file.infer_experiment.txt>; prints forward, reverse, unstranded or unclear
  awk -F': ' '
    /explained by "1\+\+,1--,2\+-,2-\+"|explained by "\+\+,--"/ { f = $2 + 0; nf = 1 }
    /explained by "1\+-,1-\+,2\+\+,2--"|explained by "\+-,-\+"/ { r = $2 + 0; nr = 1 }
    END {
      if (!nf || !nr) { print "unclear"; exit }
      if (f >= 0.8) print "forward"
      else if (r >= 0.8) print "reverse"
      else if (f >= 0.3 && f <= 0.7 && r >= 0.3 && r <= 0.7) print "unstranded"
      else print "unclear"
    }' "$1"
}
```
- All files give the same `forward`, `reverse` or `unstranded`: propose it and ask the user to confirm.
- Any file gives `unclear`, or the files disagree: show the fractions (`grep "explained by" FILE`), explain that the library type cannot be read reliably, and ask the user to check the kit; do not continue until they answer 1, 2 or 3.
- The user answered 1-3 and a result disagrees: show both and ask which is right before continuing.

### Step 4b — Genome BAM input

[BAM] Genome BAMs must come from a splice-aware aligner (STAR, as in nf-core/rnaseq: `{RNASEQ_OUTDIR}/star_salmon/{sample}.markdup.sorted.bam`), be coordinate-sorted, and be aligned to the same FASTA and GTF that Step 7 gives the pipeline (same contig names). [if `«BAM_NEEDS_BAI»` is `yes`] Each BAM needs its `.bai` index next to it; check with `ls {BAM}.bai`. Find them with `find {DIR} -name "*.bam" ! -name "*toTranscriptome*"`, derive sample names by removing `.markdup.sorted.bam`, `.sorted.bam` or `.bam`, apply the sanitisation rule of Step 4a, and use one BAM per sample (BAM rows are never merged). Ask the strandedness question of Step 4a (with the nf-core/rnaseq block when the BAMs come from such a run) and the read type (numbered): 1. paired-end · 2. single-end; store `{STRANDEDNESS}` and `{LAYOUT}`. Then apply the BAM input rule.
[NO-BAM] Genome-BAM input does not work with this pipeline revision (see the README); the rule below refuses every request, and the user starts from FASTQ.

**BAM input rule.** Define this function in the Bash tool and run `bam_input_allowed {STRANDEDNESS} {LAYOUT}`:
```bash
bam_input_allowed() {
  # usage: bam_input_allowed <unstranded|forward|reverse> <paired|single>
  # The four values below are what nf-core/rnasplice «VERSION_TAG» does with genome-BAM input (recorded by this skill's verification run).
  local RMATS_LIBTYPE="«BAM_RMATS_LIBTYPE»" RMATS_READTYPE="«BAM_RMATS_READTYPE»" DEXSEQ_STRAND="«BAM_DEXSEQ_STRAND»" FC_STRAND="«BAM_FC_STRAND»"
  local lt ds fc
  case "$1" in
    unstranded) lt=fr-unstranded; ds=no; fc=0 ;;
    forward) lt=fr-secondstrand; ds=yes; fc=1 ;;
    reverse) lt=fr-firststrand; ds=reverse; fc=2 ;;
    *) echo "REFUSED: strandedness must be unstranded, forward or reverse (got '$1')"; return 1 ;;
  esac
  case "$2" in paired|single) ;; *) echo "REFUSED: read type must be paired or single (got '$2')"; return 1 ;; esac
  case "$RMATS_LIBTYPE" in
    none) echo "REFUSED: genome-BAM input does not work with this pipeline revision. Start from FASTQ."; return 1 ;;
    from_sheet) echo "ALLOWED: strandedness and read type are written to the BAM samplesheet"; return 0 ;;
  esac
  if [ "$lt" = "$RMATS_LIBTYPE" ] && [ "$2" = "$RMATS_READTYPE" ] \
     && { [ "$DEXSEQ_STRAND" = n/a ] || [ "$ds" = "$DEXSEQ_STRAND" ]; } \
     && { [ "$FC_STRAND" = n/a ] || [ "$fc" = "$FC_STRAND" ]; }; then
    echo "ALLOWED: this revision runs rMATS on BAM input as $RMATS_LIBTYPE, $RMATS_READTYPE-end, which matches the library"; return 0
  fi
  echo "REFUSED: with BAM input this revision runs rMATS as $RMATS_LIBTYPE, $RMATS_READTYPE-end (DEXSeq -s $DEXSEQ_STRAND, featureCounts -s $FC_STRAND), but the library is $1, $2-end. Start from FASTQ instead."
  return 1
}
```
- `ALLOWED`: continue. The samplesheet header is `«BAM_SHEET_HEADER»`; [if it has `strandedness`/`single_end` columns] write `{STRANDEDNESS}` and `true` (single-end) or `false` in every row.
- `REFUSED`: show the message and ask (numbered): 1. Start from the FASTQ files instead (go to Step 4a) · 2. Stop here.

With BAM input only rMATS, DEXSeq exon usage and edgeR exon usage can run: DTU and SUPPA2 need Salmon quantification from reads, so Step 8 switches them off. The read length is asked in Step 6.
````

- [ ] **Step 4: Run the checker and both tests; expected PASS**

Run the three commands of Step 2. Expected: `PASS`, `STRANDEDNESS PASS`, `BAM POLICY PASS`.

- [ ] **Step 5: Mutation rows and RED proof**

Append to `mutations.tsv` (the last two rows only when `BAM_RMATS_LIBTYPE` is a fixed value, not `none` or `from_sheet`):
```
T2-thr	strand	skill	s/if (f >= 0.8) print "forward"/if (f >= 0.6) print "forward"/	FAIL: case pe_ambiguous
T2-swap	strand	skill	s/if (f >= 0.8) print "forward"/if (f >= 0.8) print "reverse"/	FAIL: case pe_forward
T2-se	strand	skill	s/|explained by "\\+\\+,--"//	FAIL: case se_reverse
T2-auto	check	skill	s/Strandedness is asked, never guessed\./Strandedness defaults to auto./	missing required text: Strandedness is asked, never guessed.
T2-bam-lt	check	skill	s/RMATS_LIBTYPE="@BAM_RMATS_LIBTYPE@"/RMATS_LIBTYPE="fr-wrong"/	missing required text: RMATS_LIBTYPE="@BAM_RMATS_LIBTYPE@"
T2-bam-layout	bam	skill	s/ && \[ "\$2" = "\$RMATS_READTYPE" \]//	FAIL: case
```
Run: `bash nfcore-rnasplice-setup/tests/prove_mutations.sh nfcore-rnasplice-setup/nfcore-rnasplice-setup.md` → `MUTATIONS PASS (<n>)`; `bash nfcore-rnasplice-setup/tests/prove_red.sh <BASE>` → `RED PASS (...)`.

- [ ] **Step 6: Commit**

```bash
git status && git diff --stat && git diff
git add nfcore-rnasplice-setup
git commit -m "nfcore-rnasplice-setup: Step 4 input rows, strandedness, BAM rule; tests"
```

---

### Task 3: Design — conditions, contrasts, paired design, writing the sheets (Step 5)

**Files:**
- Create: `D/tests/test_validate_sheets.sh`
- Modify: `D/tests/check_skill.sh`, `D/tests/mutations.tsv`, `D/nfcore-rnasplice-setup.md` (append Step 5)

**Interfaces:**
- Consumes: rows, `{SOURCE}`, `{STRANDEDNESS}`, `{LAYOUT}`, `{SEQ_DATE}` (Task 2); `test_env.sh`.
- Produces: `{CONDITIONS}` (ordered, reference first), `{PAIRED_DESIGN}` (`true`|`false`), `{SHEET_PREFIX}`, `{SAMPLESHEET_CSV}` = `{SHEET_PREFIX}_samplesheet.csv`, `{CONTRASTS_CSV}` = `{SHEET_PREFIX}_contrasts.csv`; skill function `validate_rnasplice_sheets <samplesheet> <contrasts> <fastq|genome_bam> <0|1>` (prints `SHEETS OK` rc 0, or `ERROR: ...` lines rc 1); anchor `**Sheet validation.**`. Step 8 uses the number of distinct samples, the smallest compared condition and, per contrast, treatment+control sample counts.

- [ ] **Step 1: Write the failing test**

Create `D/tests/test_validate_sheets.sh`:
```bash
#!/bin/bash
# Usage: test_validate_sheets.sh <skill.md>
# Runs validate_rnasplice_sheets from the skill (block after "**Sheet validation.**") on good and bad sheets.
# Prints "VALIDATE PASS" or exits 1.
set -u
source "$(cd "$(dirname "$0")" && pwd)/test_env.sh" "$@"
SKILL=${1:?usage: test_validate_sheets.sh <skill.md>}
HERE=$(cd "$(dirname "$0")" && pwd); FIX="$HERE/fixtures"
gv() { awk -F'\t' -v k="$1" '$1 == k {print $2; exit}' "$FIX/gate_values.tsv"; }
CODE=$(bash "$HERE/cut_block.sh" "$SKILL" "**Sheet validation.**") || { echo "FAIL: sheet validation block not found"; exit 1; }
eval "$CODE"
declare -F validate_rnasplice_sheets >/dev/null || { echo "FAIL: validate_rnasplice_sheets is not defined by the block"; exit 1; }
S="$TEST_TMP/s.csv"; C="$TEST_TMP/c.csv"; fail=0
ok()  { local out; out=$(validate_rnasplice_sheets "$S" "$C" "$2" "$3"); [ $? -eq 0 ] && [ "$out" = "SHEETS OK" ] || { echo "FAIL: case $1: expected SHEETS OK, got: $out"; fail=1; }; }
bad() { local out; out=$(validate_rnasplice_sheets "$S" "$C" "$2" "$3"); [ $? -ne 0 ] && grep -qF -- "$4" <<< "$out" || { echo "FAIL: case $1: expected an error containing '$4', got: $out"; fail=1; }; }
good_fastq() { printf '%s\n' "sample,fastq_1,fastq_2,strandedness,condition" "WT_1,a_R1.fastq.gz,a_R2.fastq.gz,reverse,WT" "WT_2,b_R1.fastq.gz,b_R2.fastq.gz,reverse,WT" "KO_1,c_R1.fastq.gz,c_R2.fastq.gz,reverse,KO" "KO_2,d_R1.fastq.gz,d_R2.fastq.gz,reverse,KO" > "$S"; }
good_con() { printf '%s\n' "contrast,treatment,control" "KO_vs_WT,KO,WT" > "$C"; }

cp "$FIX/test_samplesheet.csv" "$S"; cp "$FIX/test_contrastsheet.csv" "$C"; ok upstream_test_sheets fastq 0
good_fastq; good_con; ok good fastq 0
good_fastq; printf '%s\n' "contrast,treatment,control" "KO_vs_WT,K0,WT" > "$C"; bad label_typo_treatment fastq 0 'treatment "K0" of contrast KO_vs_WT is not a value of the condition column'
good_fastq; printf '%s\n' "contrast,treatment,control" "KO_vs_WT,KO,wt" > "$C"; bad label_typo_control fastq 0 'control "wt" of contrast KO_vs_WT is not a value of the condition column'
good_fastq; printf '%s\n' "KO_vs_WT,KO,WT" > "$C"; bad contrasts_no_header fastq 0 'contrasts header is'
good_fastq; printf '%s\n' "contrast,treatment,control" "KO_vs_KO,KO,KO" > "$C"; bad self_contrast fastq 0 'compares KO with itself'
good_fastq; printf '%s\n' "contrast,treatment,control" "KO_vs_WT,KO,WT" "KO_vs_WT,WT,KO" > "$C"; bad duplicate_contrast fastq 0 'contrast name KO_vs_WT is used twice'
good_fastq; printf '%s\n' "contrast,treatment,control" > "$C"; bad no_contrast_rows fastq 0 'has no contrast rows'
good_fastq; sed -i '2s/,reverse,/,forward,/' "$S"; good_con; bad mixed_strandedness fastq 0 'mixed strandedness'
good_fastq; sed -i 's/,reverse,/,auto,/' "$S"; good_con; bad strandedness_auto fastq 0 'strandedness "auto"'
good_fastq; sed -i '3s/,b_R2.fastq.gz,/,,/' "$S"; good_con; bad mixed_layout fastq 0 'mixed single-end and paired-end'
good_fastq; sed -i '5d' "$S"; good_con; bad one_sample_condition fastq 0 'condition KO has 1 sample'
good_fastq; sed -i '2s/^WT_1/WT-1/' "$S"; good_con; bad dash_in_name fastq 0 'sample name "WT-1"'
good_fastq; sed -i '2s/^WT_1/1WT/' "$S"; good_con; bad digit_first fastq 0 'sample name "1WT"'
good_fastq; sed -i '2s/,WT$/,wild-type/' "$S"; good_con; bad bad_condition fastq 0 'condition "wild-type"'
good_fastq; printf '%s\n' "WT_1,a2_R1.fastq.gz,a2_R2.fastq.gz,reverse,WT" >> "$S"; good_con; ok tech_replicates fastq 0
good_fastq; printf '%s\n' "WT_1,a2_R1.fastq.gz,a2_R2.fastq.gz,reverse,KO" >> "$S"; good_con; bad tech_replicate_conflict fastq 0 'rows of sample WT_1 have different conditions'
good_fastq; sed -i '3s/^WT_2/WT_1/' "$S"; good_con; bad replicates_hide_size fastq 0 'condition WT has 1 sample'
good_fastq; sed -i '2s/a_R1.fastq.gz/a_R1.fq/' "$S"; good_con; bad not_gz fastq 0 'does not end in .fastq.gz or .fq.gz'
good_fastq; sed -i '2s/,a_R2.fastq.gz,/,"a_R2.fastq.gz",/' "$S"; good_con; bad quoted fastq 0 'contains a double quote'
good_fastq; good_con; sed -i 's/$/\r/' "$S" "$C"; ok crlf fastq 0
good_fastq; sed -i '1s/^sample,/Sample,/' "$S"; good_con; bad samplesheet_header fastq 0 'samplesheet header is'
good_fastq; good_con; ok paired_two_by_two fastq 1
good_fastq; printf '%s\n' "HET_1,e_R1.fastq.gz,e_R2.fastq.gz,reverse,HET" "HET_2,f_R1.fastq.gz,f_R2.fastq.gz,reverse,HET" >> "$S"; good_con; bad paired_three_conditions fastq 1 'exactly two conditions'
good_fastq; printf '%s\n' "KO_3,e_R1.fastq.gz,e_R2.fastq.gz,reverse,KO" >> "$S"; good_con; bad paired_unequal fastq 1 'same number of samples in both conditions'
: > "$S"; good_con; bad empty_samplesheet fastq 0 'missing or empty'
# genome-BAM sheets with the header recorded by the verification run
BH=$(gv BAM_SHEET_HEADER)
bam_rows() { awk -F',' -v h="$BH" 'BEGIN { n = split(h, c, ","); print h }
  { for (i = 1; i <= n; i++) { v = (c[i] == "sample") ? $1 : (c[i] == "condition") ? $2 : (c[i] == "genome_bam") ? $3 : (c[i] == "strandedness") ? "reverse" : (c[i] == "single_end") ? "false" : "x"; printf "%s%s", v, (i < n ? "," : "\n") } }'; }
printf '%s\n' "WT_1,WT,a.bam" "WT_2,WT,b.bam" "KO_1,KO,c.bam" "KO_2,KO,d.bam" | bam_rows > "$S"; good_con; ok bam_good genome_bam 0
printf '%s\n' "WT_1,WT,a.bam" "WT_1,WT,b.bam" "KO_1,KO,c.bam" "KO_2,KO,d.bam" | bam_rows > "$S"; good_con; bad bam_duplicate genome_bam 0 'appears twice'
printf '%s\n' "WT_1,WT,a.cram" "WT_2,WT,b.bam" "KO_1,KO,c.bam" "KO_2,KO,d.bam" | bam_rows > "$S"; good_con; bad bam_not_bam genome_bam 0 'does not end in .bam'
forbidden_ran && fail=1
[ $fail -eq 0 ] && echo "VALIDATE PASS" || exit 1
```
If the case `upstream_test_sheets` fails because the upstream sheet violates a rule of this skill (for example a sample name with `-`), the implementer reports `BLOCKED` with the error line; the controller rules (never by relaxing the validator silently).

- [ ] **Step 2: Append the Task 3 checker block; run the test and the checker to see them fail**

```bash
# --- Task 3 (Step 5)
need "## Step 5 — Conditions, contrasts and the two sheets"
need "contrast,treatment,control"
need "$(head -n1 "$FIX/test_contrastsheet.csv" | tr -d '\r')"
need "{treatment}_vs_{control}"
need "1. All pairwise comparisons · 2. Only the comparisons I list"
need "1. No, the samples are independent (default) · 2. Yes, every sample has a partner"
need "the pipeline default \`rmats_paired_stats: true\` assumes a pairing"
need "Every condition used in a contrast needs at least 2 samples"
anchor_once "**Sheet validation.**"
need "validate_rnasplice_sheets() {"
need "bamhdr=\"$(gv BAM_SHEET_HEADER)\""
need "{SHEET_PREFIX}_samplesheet.csv"
need "{SHEET_PREFIX}_contrasts.csv"
need "never edit the validator"
need "1. overwrite · 2. choose another filename"
forbid "rm -rf"
case "$(gv RMATS_BAMLIST_ORDER)" in
  sheet) need "write the rows ordered by subject within each condition" ;;
  sorted_by_name) need "propose names \`{SUBJECT}_{CONDITION}\`" ;;
  unordered) need "a paired rMATS test cannot be set up safely" ;;
esac
# --- end Task 3
```
Expected: `test_validate_sheets.sh` → `FAIL: sheet validation block not found`; checker → Task 3 `FAIL` lines.

- [ ] **Step 3: Append Step 5 to the skill**

Write only the paired-design paragraph that matches `RMATS_BAMLIST_ORDER`; substitute `«BAM_SHEET_HEADER»`.

````markdown

---

## Step 5 — Conditions, contrasts and the two sheets

**Conditions.** Ask the user to describe the groups in plain language (for example "samples 1-3 are wild type, 4-6 are Rbpms2 knockout"), interpret it, and show a table sample → condition; ask "Is this correct?" until confirmed. Condition labels use letters, digits and `_` and start with a letter (suggest short labels such as `WT`, `KO`); rows of one sample (technical replicates) share its condition. Every condition used in a contrast needs at least 2 samples (distinct sample names), because rMATS, DEXSeq and edgeR estimate variability from replicates; if one has fewer, say so and stop until the design is changed.

Then ask for the reference: with two conditions, "Which condition is the reference (control)?" (numbered); with more, ask for an order with the reference first. Store `{CONDITIONS}` = the labels in that order.

**Contrasts.** Ask (numbered): 1. All pairwise comparisons · 2. Only the comparisons I list.
- Option 1: for every pair of conditions, treatment = the later and control = the earlier one in `{CONDITIONS}`; contrast name `{treatment}_vs_{control}`. With two conditions this is one contrast.
- Option 2: ask for lines "treatment vs control", one per contrast; names as above.
Show the contrasts as a table (contrast, treatment, control). rMATS runs one prep/post pair per contrast; with more than 6 contrasts tell the user that the run takes correspondingly longer.

**Paired design.** Ask only when there are exactly two conditions with the same number of samples; otherwise set `{PAIRED_DESIGN}` = `false` without asking. Ask (numbered): "Are the samples paired (each sample of one condition has a partner from the same individual, litter or batch in the other)?" 1. No, the samples are independent (default) · 2. Yes, every sample has a partner. Explain: the pipeline default `rmats_paired_stats: true` assumes a pairing; on independent samples it gives wrong statistics, so this skill writes `false` unless the user confirms pairing. Answer 1: `{PAIRED_DESIGN}` = `false`.
[sheet] Answer 2: ask for the partner of each sample (a subject label), check that each subject has exactly one sample in each condition, and write the rows ordered by subject within each condition (rMATS pairs the i-th sample of one group with the i-th of the other, in samplesheet order). `{PAIRED_DESIGN}` = `true`.
[sorted_by_name] Answer 2: ask for the partner of each sample (a subject label) and check that sorting the sample names of each condition gives the subjects in the same order (this pipeline revision sorts the samples of each group by name before rMATS pairs them by position); if not, propose names `{SUBJECT}_{CONDITION}` and confirm them. `{PAIRED_DESIGN}` = `true`.
[unordered] Answer 2: tell the user that this pipeline revision does not keep a fixed sample order for rMATS, so a paired rMATS test cannot be set up safely; `{PAIRED_DESIGN}` = `false` (SUPPA2's paired test is switched off for the same reason).
`{PAIRED_DESIGN}` sets both `rmats_paired_stats` and SUPPA2's `diffsplice_paired` in Step 8.

**File names** (numbered): 1. `{SEQ_DATE}_{WD_NAME}` · 2. `{TODAY_YYMMDD}_{WD_NAME}` · 3. Custom prefix. Store `{SHEET_PREFIX}`; `{SAMPLESHEET_CSV}` = `{SHEET_PREFIX}_samplesheet.csv` and `{CONTRASTS_CSV}` = `{SHEET_PREFIX}_contrasts.csv`.

**Sheet validation.** Both sheets are written to a scratch directory first, validated there, and only then moved into `{CWD}`. Define this function in the Bash tool (awk only; light work on the login node):
```bash
validate_rnasplice_sheets() {
  # usage: validate_rnasplice_sheets <samplesheet.csv> <contrasts.csv> <fastq|genome_bam> <paired design: 0|1>
  # prints "SHEETS OK", or one "ERROR: ..." line per problem and returns 1
  [ -s "$1" ] && [ -s "$2" ] || { echo "ERROR: samplesheet or contrasts file missing or empty"; return 1; }
  awk -F',' -v src="$3" -v paired="$4" \
      -v fqhdr="sample,fastq_1,fastq_2,strandedness,condition" -v bamhdr="«BAM_SHEET_HEADER»" '
    function err(m) { print "ERROR: " m; bad = 1 }
    FNR == 1 { file++ }
    { sub(/\r$/, "") }
    /"/ { err("file " file " line " FNR " contains a double quote; write plain comma-separated values"); next }
    file == 1 && FNR == 1 {
      want = (src == "fastq") ? fqhdr : bamhdr
      if ($0 != want) err("samplesheet header is \"" $0 "\", expected \"" want "\"")
      ncol = split(want, h, ","); for (i = 1; i <= ncol; i++) col[h[i]] = i
      next
    }
    file == 1 && $0 == "" { next }
    file == 1 {
      if (NF != ncol) { err("samplesheet line " FNR " has " NF " fields, expected " ncol); next }
      s = $col["sample"]; c = $col["condition"]
      if (s !~ /^[A-Za-z][A-Za-z0-9_]*$/) err("sample name \"" s "\" (line " FNR "): letters, digits and _ only, starting with a letter")
      if (c !~ /^[A-Za-z][A-Za-z0-9_]*$/) err("condition \"" c "\" of " s ": letters, digits and _ only, starting with a letter")
      if (s in cond) {
        if (src != "fastq") err("sample " s " appears twice; a BAM samplesheet has one row per sample")
        else if (cond[s] != c) err("rows of sample " s " have different conditions (" cond[s] ", " c "); rows with one sample name are merged into one sample")
      } else nsamp[c]++
      cond[s] = c; conds[c] = 1
      if (src == "fastq") {
        if ($col["fastq_1"] !~ /\.f(ast)?q\.gz$/) err("fastq_1 of " s " does not end in .fastq.gz or .fq.gz")
        if ($col["fastq_2"] != "" && $col["fastq_2"] !~ /\.f(ast)?q\.gz$/) err("fastq_2 of " s " does not end in .fastq.gz or .fq.gz")
        lay = ($col["fastq_2"] == "") ? "single-end" : "paired-end"
      } else {
        if ($col["genome_bam"] !~ /\.bam$/) err("genome_bam of " s " does not end in .bam")
        lay = "n/a"
        if ("single_end" in col) {
          if ($col["single_end"] !~ /^(true|false)$/) err("single_end of " s " must be true or false")
          lay = ($col["single_end"] == "true") ? "single-end" : "paired-end"
        }
      }
      if (lay != "n/a") { if (layout == "") layout = lay; else if (lay != layout) err("mixed single-end and paired-end samples (" s "); rMATS needs one read type for all samples") }
      if ("strandedness" in col) {
        st = $col["strandedness"]
        if (st !~ /^(unstranded|forward|reverse)$/) err("strandedness \"" st "\" of " s ": use unstranded, forward or reverse")
        else if (strand == "") strand = st; else if (st != strand) err("mixed strandedness (" strand ", " st "); rMATS needs one strandedness for all samples")
      }
      next
    }
    file == 2 && FNR == 1 { if ($0 != "contrast,treatment,control") err("contrasts header is \"" $0 "\", expected \"contrast,treatment,control\""); next }
    file == 2 && $0 == "" { next }
    file == 2 {
      if (NF != 3) { err("contrasts line " FNR " has " NF " fields, expected 3"); next }
      ncon++
      if ($1 !~ /^[A-Za-z][A-Za-z0-9_-]*$/) err("contrast name \"" $1 "\": letters, digits, _ and - only, starting with a letter")
      if ($1 in seen) err("contrast name " $1 " is used twice"); seen[$1] = 1
      if (!($2 in conds)) err("treatment \"" $2 "\" of contrast " $1 " is not a value of the condition column")
      if (!($3 in conds)) err("control \"" $3 "\" of contrast " $1 " is not a value of the condition column")
      if ($2 == $3) err("contrast " $1 " compares " $2 " with itself")
      used[$2] = 1; used[$3] = 1
      next
    }
    END {
      if (file < 2) err("contrasts file was not read")
      if (ncon == 0) err("contrasts file has no contrast rows")
      for (c in used) if ((c in nsamp) && nsamp[c] < 2) err("condition " c " has " nsamp[c] " sample; each compared condition needs at least 2")
      if (paired == 1) {
        k = 0; for (c in conds) { k++; sz[k] = nsamp[c] }
        if (k != 2) err("a paired design needs exactly two conditions, found " k)
        else if (sz[1] != sz[2]) err("a paired design needs the same number of samples in both conditions (" sz[1] ", " sz[2] ")")
      }
      if (!bad) print "SHEETS OK"
      exit bad
    }' "$1" "$2"
}
```
Procedure:
1. `T=$(mktemp -d)`; write `$T/samplesheet.csv` (header `sample,fastq_1,fastq_2,strandedness,condition` for FASTQ with `{STRANDEDNESS}` in every row, or `«BAM_SHEET_HEADER»` for BAM; the sanitised names) and `$T/contrasts.csv` (header `contrast,treatment,control`).
2. Run `validate_rnasplice_sheets "$T/samplesheet.csv" "$T/contrasts.csv" {SOURCE} 0` (`1` instead of `0` when `{PAIRED_DESIGN}` is `true`). On `ERROR` lines, fix the cause with the user — never edit the validator — and repeat.
3. Check that every FASTQ or BAM path in the sheet exists: `test -s` on each (relative paths from `{CWD}`).
4. For each of `{SAMPLESHEET_CSV}` and `{CONTRASTS_CSV}`: if it already exists in `{CWD}`, ask (numbered): 1. overwrite · 2. choose another filename. Then move both files into `{CWD}` under their names and remove the scratch directory with `rmdir "$T"` (it is empty after the moves).
5. Show both files in full.
````

- [ ] **Step 4: Run test and checker; expected PASS**

Run: `bash nfcore-rnasplice-setup/tests/test_validate_sheets.sh nfcore-rnasplice-setup/nfcore-rnasplice-setup.md` → `VALIDATE PASS`; checker → `PASS`.

- [ ] **Step 5: Mutation rows, RED proof, commit**

Append to `mutations.tsv`:
```
T3-label	validate	skill	/err("treatment /d	FAIL: case label_typo_treatment
T3-header	validate	skill	s/if ($0 != "contrast,treatment,control")/if (0)/	FAIL: case contrasts_no_header
T3-min	validate	skill	s/nsamp\[c\] < 2/nsamp[c] < 1/	FAIL: case one_sample_condition
T3-strand	validate	skill	s/else if (st != strand) err/else if (0) err/	FAIL: case mixed_strandedness
T3-paired	validate	skill	s/if (k != 2) err/if (0) err/	FAIL: case paired_three_conditions
T3-rmdir	check	skill	s/remove the scratch directory with `rmdir "\$T"`/remove the scratch directory with `rm -rf "$T"`/	forbidden text present: rm -rf
```
```bash
bash nfcore-rnasplice-setup/tests/prove_mutations.sh nfcore-rnasplice-setup/nfcore-rnasplice-setup.md
bash nfcore-rnasplice-setup/tests/prove_red.sh <BASE>
git status && git diff --stat && git diff
git add nfcore-rnasplice-setup
git commit -m "nfcore-rnasplice-setup: Step 5 conditions, contrasts, paired design, validated sheets"
```
Expected: `MUTATIONS PASS (<n>)`, `RED PASS (...)`.

---

### Task 4: Read length and genome (Steps 6-7)

**Files:**
- Create: `D/tests/test_read_length.sh`
- Modify: `D/tests/check_skill.sh`, `D/tests/mutations.tsv`, `D/nfcore-rnasplice-setup.md` (append Steps 6-7)

**Interfaces:**
- Consumes: `{SOURCE}`, rows (Task 2); `test_env.sh`.
- Produces: `{READ_LENGTH}` (integer), `{ORGANISM}`, `{GENOME_DIR}`, `{FASTA_PATH}`, `{GTF_PATH}`, `{REF_TAG}` (`{ASSEMBLY}_ens{ENS_VERSION}` or `custom_{WD_NAME}`), `{STAR_INDEX}` and `{SALMON_INDEX}` (directory or empty), `{ENSEMBL_FASTA_URL}`, `{ENSEMBL_GTF_URL}` (only when a download is needed); skill function `detect_read_length <fastq.gz>...` (prints `<length> <reads>`); anchor `**Read-length detection.**`.

- [ ] **Step 1: Write the failing test**

Create `D/tests/test_read_length.sh`:
```bash
#!/bin/bash
# Usage: test_read_length.sh <skill.md>
# Runs detect_read_length from the skill (block after "**Read-length detection.**") on generated gzipped FASTQ.
# Prints "READ LENGTH PASS" or exits 1.
set -u
source "$(cd "$(dirname "$0")" && pwd)/test_env.sh" "$@"
SKILL=${1:?usage: test_read_length.sh <skill.md>}
HERE=$(cd "$(dirname "$0")" && pwd)
CODE=$(bash "$HERE/cut_block.sh" "$SKILL" "**Read-length detection.**") || { echo "FAIL: read-length block not found"; exit 1; }
eval "$CODE"
declare -F detect_read_length >/dev/null || { echo "FAIL: detect_read_length is not defined by the block"; exit 1; }
reads() { awk -v L="$1" -v n="$2" 'BEGIN { s = ""; for (i = 0; i < L; i++) s = s "A"; q = s; gsub(/A/, "I", q); for (r = 1; r <= n; r++) printf "@r%d\n%s\n+\n%s\n", r, s, q }'; }
fq() { reads "$2" "$3" | gzip -c > "$TEST_TMP/$1"; }
fail=0
check() { local name=$1 want=$2 got; shift 2; got=$(detect_read_length "$@"); [ "$got" = "$want" ] || { echo "FAIL: case $name: expected '$want', got '$got'"; fail=1; }; }
fq a.fastq.gz 101 30; fq b.fastq.gz 101 30; fq c.fastq.gz 76 50; fq d.fq.gz 151 10; fq e.fastq.gz 76 30
{ reads 60 1000; reads 90 2000; } | gzip -c > "$TEST_TMP/big.fastq.gz"
check single "101 30" "$TEST_TMP/a.fastq.gz"
check mixed "101 60" "$TEST_TMP/a.fastq.gz" "$TEST_TMP/b.fastq.gz" "$TEST_TMP/c.fastq.gz"
check tie "101 30" "$TEST_TMP/a.fastq.gz" "$TEST_TMP/e.fastq.gz"
check fq_gz "151 10" "$TEST_TMP/d.fq.gz"
check first1000 "60 1000" "$TEST_TMP/big.fastq.gz"
forbidden_ran && fail=1
[ $fail -eq 0 ] && echo "READ LENGTH PASS" || exit 1
```

- [ ] **Step 2: Append the Task 4 checker block; see both fail**

```bash
# --- Task 4 (Steps 6-7)
need "## Step 6 — Read length for rMATS"
anchor_once "**Read-length detection.**"
need "detect_read_length() {"
need "the pipeline default of 40 is wrong for almost all data"
need "\`--variable-read-length\` and \`--allow-clipping\`"
need "an integer between 20 and 1000"
need "## Step 7 — Organism and genome files"
need "{genome_base}/{organism}/{assembly}_ens{version}/"
need "Ask the base directory only for option 1"
need "custom_{WD_NAME}"
need "versionGenome"
need "$(gv STAR_VERSION_GENOME)"
need "\"indexVersion\""
need "$(gv SALMON_INDEX_VERSION)"
need "1. Reuse (default) · 2. Let the pipeline build its own"
need "\`gencode: false\` is written either way"
need "download_genome_{REF_TAG}.sh"
need "Never type a URL from memory"
# --- end Task 4
```

- [ ] **Step 3: Append Steps 6-7 to the skill**

````markdown

---

## Step 6 — Read length for rMATS

rMATS needs the read length (`rmats_read_len`), and the pipeline default of 40 is wrong for almost all data, so the length is always set here and written to the params file. rnasplice has no other read-length setting (the STAR index it builds does not depend on it), and rMATS always runs with `--variable-read-length` and `--allow-clipping`, so trimmed reads of other lengths are still counted.

**Read-length detection.** FASTQ input: define this function in the Bash tool and run it on the first FASTQ of up to 5 different samples, one file at a time (report per sample), then on all of them together (the value used). It reads only the first 1000 reads of each file, which is light work on the login node.
```bash
detect_read_length() {
  # usage: detect_read_length <fastq.gz> ...
  # prints "<length> <reads>": the most common read length among the first 1000 reads of each file (a tie goes to the longer length)
  local f
  for f in "$@"; do zcat -f "$f" | head -n 4000 | awk 'NR % 4 == 2 {print length($0)}'; done \
    | sort -n | uniq -c | sort -k1,1nr -k2,2nr | head -1 | awk '{print $2, $1}'
}
```
`{READ_LENGTH}` = the first number of the combined result. If the per-sample values differ, show them and say that `{READ_LENGTH}` is the most common length, which is what rMATS expects for variable read lengths. Tell the user: "Read length for rMATS: {READ_LENGTH} bp (`rmats_read_len: {READ_LENGTH}`)."

BAM input: ask "What is the read length of the sequencing (for example 100 or 150)?" — an integer between 20 and 1000; store it as `{READ_LENGTH}`.

---

## Step 7 — Organism and genome files

Ask, in this order:
1. "What organism is this data from? (e.g. mouse, human)"
2. (numbered): 1. Ensembl release in the standard folder (default) · 2. Custom reference — I already have a FASTA and GTF.
3. The option-specific questions below. Ask the base directory only for option 1.

- **Option 1 (Ensembl):** ask "What is the base directory where genome files and indexes are stored?" (`{genome_base}`). Folder convention, shared with `/nfcore-rnaseq-setup` so that FASTA, GTF and indexes are reused, never downloaded twice:
  ```
  {genome_base}/{organism}/{assembly}_ens{version}/
  ├── {FASTA}.fa            ← primary assembly FASTA
  ├── {GTF}.gtf             ← annotation GTF
  └── index/
      ├── star/             ← STAR index (as built by /nfcore-rnaseq-setup); reused only if compatible
      └── salmon/           ← Salmon index (as built by /nfcore-rnaseq-setup); reused only if compatible
  ```
  Mouse: GRCm39, FASTA `Mus_musculus.GRCm39.dna.primary_assembly.fa`, GTF `Mus_musculus.GRCm39.{version}.gtf`, directory `{genome_base}/mouse/mm39_ens{version}/`. Human: GRCh38, FASTA `Homo_sapiens.GRCh38.dna.primary_assembly.fa`, GTF `Homo_sapiens.GRCh38.{version}.gtf`, directory `{genome_base}/human/hg38_ens{version}/`. Other organisms: option 2. Version: the highest existing `{assembly}_ens{N}` directory unless the user asks otherwise; if none exists, the latest Ensembl release from `https://ftp.ensembl.org/pub/current_README` (`WebFetch`). Store `{GENOME_DIR}`, `{FASTA_PATH}`, `{GTF_PATH}`, `{ORGANISM}`, `{ASSEMBLY}`, `{ENS_VERSION}` and `{REF_TAG}` = `{ASSEMBLY}_ens{ENS_VERSION}`.
  If the FASTA or the GTF is missing: resolve both download URLs at run time with `WebFetch` on the Ensembl FTP listing of that release (`https://ftp.ensembl.org/pub/release-{version}/fasta/{species}/dna/` and `https://ftp.ensembl.org/pub/release-{version}/gtf/{species}/`), check each with a HEAD request (`curl -sI URL | head -1` must show 200), show them to the user, store `{ENSEMBL_FASTA_URL}` and `{ENSEMBL_GTF_URL}`, and generate `download_genome_{REF_TAG}.sh` (Step 11). Never type a URL from memory.
- **Option 2 (Custom reference):** ask for the FASTA path and the GTF path; `{GENOME_DIR}` = the directory of the FASTA, `{REF_TAG}` = `custom_{WD_NAME}`; no download. The pipeline accepts `.fa`, `.fasta`, `.fna` and the same with `.gz`, and `.gtf` or `.gtf.gz`. Check that both files exist (`test -s`). Ask (numbered) whether STAR or Salmon indexes built from exactly this FASTA and GTF exist: 1. No (default) · 2. Yes — ask for their directories and apply the compatibility checks below.

**GTF source.** `zcat -f {GTF_PATH} | grep -v "^#" | head -3`: gene IDs with a version suffix (`ENSG00000000003.15`) indicate GENCODE, without one Ensembl. Report it. `gencode: false` is written either way: the pipeline builds the transcript FASTA from the GTF itself, so it never receives a GENCODE-format transcript FASTA, the only case `gencode: true` is for.

**Existing indexes — reused only when compatible; otherwise the pipeline builds its own.** (BAM input uses no index: skip this part.)
- STAR (`{GENOME_DIR}/index/star/` for option 1): present when `SA`, `Genome`, `sjdbList.out.tab` and `genomeParameters.txt` exist. Read `grep '^versionGenome' {DIR}/genomeParameters.txt` and `grep '^sjdbOverhang' {DIR}/genomeParameters.txt`. Offer reuse — (numbered) 1. Reuse (default) · 2. Let the pipeline build its own — only when `versionGenome` is `«STAR_VERSION_GENOME»` (the index format of the STAR inside this pipeline revision); otherwise say why and let the pipeline build it. Report the `sjdbOverhang` (a value near the read length minus 1, or 100, works well). Reuse: `{STAR_INDEX}` = that directory; otherwise `{STAR_INDEX}` is empty.
- Salmon (`{GENOME_DIR}/index/salmon/` for option 1; FASTQ input only): present when `versionInfo.json` exists; offer reuse the same way only when its `"indexVersion"` is `«SALMON_INDEX_VERSION»`. Reuse: `{SALMON_INDEX}` = that directory; otherwise empty.
When the pipeline builds the STAR index of a human or mouse genome, STAR_GENOMEGENERATE needs about 32 GB of memory and an hour or more; the `nextflow.config` of Step 10 gives it 64 GB and 8 h. Indexes built by the pipeline are not kept (`save_reference: false`), so a later run builds them again.
````

- [ ] **Step 4: Run test and checker; expected PASS**

`bash nfcore-rnasplice-setup/tests/test_read_length.sh nfcore-rnasplice-setup/nfcore-rnasplice-setup.md` → `READ LENGTH PASS`; checker → `PASS`.

- [ ] **Step 5: Mutation rows, RED proof, commit**

```
T4-tie	readlen	skill	s/sort -k1,1nr -k2,2nr/sort -k1,1nr -k2,2n/	FAIL: case tie
T4-mode	readlen	skill	s/-k2,2nr | head -1 |/-k2,2nr | tail -1 |/	FAIL: case mixed
T4-head	readlen	skill	s/head -n 4000 | awk/awk/	FAIL: case first1000
T4-star	check	skill	s/@STAR_VERSION_GENOME@/9.9.9z/g	missing required text: @STAR_VERSION_GENOME@
```
```bash
bash nfcore-rnasplice-setup/tests/prove_mutations.sh nfcore-rnasplice-setup/nfcore-rnasplice-setup.md
bash nfcore-rnasplice-setup/tests/prove_red.sh <BASE>
git status && git diff --stat && git diff
git add nfcore-rnasplice-setup
git commit -m "nfcore-rnasplice-setup: Steps 6-7 read length for rMATS, genome folder, index reuse"
```

---

### Task 5: Modules, their options and the parameter policy (Steps 8-9)

**Files:**
- Modify: `D/tests/check_skill.sh`, `D/tests/mutations.tsv`, `D/nfcore-rnasplice-setup.md` (append Steps 8-9)

**Interfaces:**
- Consumes: `{READ_LENGTH}` (Task 4), `{PAIRED_DESIGN}`, sample and condition counts (Task 3), `{SOURCE}` (Task 2); `$config_kv`, `$schema_names` (Task 1).
- Produces: `{RUN_RMATS}`, `{RUN_DEXSEQ_EXON}`, `{RUN_EDGER_EXON}`, `{RUN_DEXSEQ_DTU}`, `{RUN_SUPPA}`, `{RMATS_NOVEL}` (each `true`|`false`), `{MIN_SAMPS_GENE_EXPR}`, `{MIN_SAMPS_FEATURE}` (integers); the yaml block after the anchor `**Module and option keys (written to the params file).**` (Task 6 inserts it for `{MODULE_PARAMS}`); in the checker: `DEV` (deliberate deviations) and `defaults_equal "<block text>"` (Task 6 reuses it).

- [ ] **Step 1: Append the Task 5 checker block; see it fail**

```bash
# --- Task 5 (Steps 8-9): module block
need "## Step 8 — Analyses (modules) and their settings"
need "every analysis module is switched on, so a module that is not mentioned in the params file runs anyway"
need "Default (empty answer): 1 only"
need "This skill always writes \`sashimi_plot: false\`"
need "it is the same analysis twice"
need "## Step 9 — Trimming and QC"
anchor_once "**Module and option keys (written to the params file).**"
MOD=$(bash "$CUT" "$SKILL" "**Module and option keys (written to the params file).**" 2>/dev/null)
[ -n "$MOD" ] || { echo "FAIL: module block not found"; fail=1; }
for sw in rmats dexseq_exon edger_exon dexseq_dtu suppa sashimi_plot; do
  [ "$(printf '%s\n' "$MOD" | grep -cE "^$sw: (true|false|\{RUN_[A-Z_]+\})$")" -eq 1 ] \
    && [ "$(printf '%s\n' "$MOD" | grep -c "^$sw:")" -eq 1 ] \
    || { echo "FAIL: module switch $sw must be written exactly once with an explicit value"; fail=1; }
done
printf '%s\n' "$MOD" | grep -qx 'sashimi_plot: false' || { echo "FAIL: sashimi_plot must be written as false"; fail=1; }
for kv in 'rmats: {RUN_RMATS}' 'dexseq_exon: {RUN_DEXSEQ_EXON}' 'edger_exon: {RUN_EDGER_EXON}' 'dexseq_dtu: {RUN_DEXSEQ_DTU}' 'suppa: {RUN_SUPPA}' \
          'rmats_read_len: {READ_LENGTH}' 'rmats_paired_stats: {PAIRED_DESIGN}' 'diffsplice_paired: {PAIRED_DESIGN}' 'rmats_novel_splice_site: {RMATS_NOVEL}' \
          'min_samps_gene_expr: {MIN_SAMPS_GENE_EXPR}' 'min_samps_feature_expr: {MIN_SAMPS_FEATURE}' 'min_samps_feature_prop: {MIN_SAMPS_FEATURE}' \
          'min_gene_expr: 10' 'min_feature_expr: 10' 'min_feature_prop: 0.1' 'dtu_txi: "dtuScaledTPM"' 'rmats_splice_diff_cutoff: 0.0001'; do
  printf '%s\n' "$MOD" | grep -qxF -- "$kv" || { echo "FAIL: module block must contain the line: $kv"; fail=1; }
done
case "$(gv SALMON_ROUTE)" in
  pseudo_only) r1='aligner: "star"'; r2='pseudo_aligner: "salmon"' ;;
  star_salmon_only) r1='aligner: "star_salmon"'; r2="$(gv PSEUDO_OFF_LINE)" ;;
  star_salmon_both) r1='aligner: "star_salmon"'; r2='pseudo_aligner: "salmon"' ;;
esac
for kv in "$r1" "$r2"; do printf '%s\n' "$MOD" | grep -qxF -- "$kv" || { echo "FAIL: Salmon route $(gv SALMON_ROUTE) needs the line: $kv"; fail=1; }; done
for sw in isoformswitchanalyzer leafcutter; do
  if echo "$schema_names" | grep -qx "$sw"; then printf '%s\n' "$MOD" | grep -qx "$sw: false" || { echo "FAIL: $sw is a parameter of this revision and must be written as false"; fail=1; }; fi
done
case "$(gv DTU_FILTER_SCOPE)" in
  all_samples) need "the number of samples in the samplesheet (distinct names), because the filter is applied once to all samples" ;;
  per_contrast) need "the smallest number of samples in a contrast (treatment plus control), because the filter is applied per contrast" ;;
esac
# literal values equal the recorded pipeline config, except the deliberate deviations in DEV
DEV=" aligner pseudo_aligner rmats dexseq_exon edger_exon dexseq_dtu suppa sashimi_plot isoformswitchanalyzer leafcutter max_cpus max_memory max_time "
defaults_equal() {
  local line k v cv
  while IFS= read -r line; do
    k=${line%%:*}; v=${line#*: }
    case "$v" in *"{"*) continue ;; esac
    case "$DEV" in *" $k "*) continue ;; esac
    v=${v#\"}; v=${v%\"}
    cv=$(printf '%s\n' "$config_kv" | awk -F'\t' -v k="$k" '$1 == k {print $2; f = 1} END {exit !f}') \
      || { echo "FAIL: key $k has no default in the recorded pipeline config"; fail=1; continue; }
    [ "$v" = "$cv" ] || { echo "FAIL: key $k = $v differs from the pipeline config default $cv (a deliberate deviation needs a written reason and an entry in DEV)"; fail=1; }
  done < <(printf '%s\n' "$1" | grep -E '^[a-z_0-9]+: ')
}
defaults_equal "$MOD"
# --- end Task 5
```
Run the checker. Expected: Task 5 `FAIL` lines (`module block not found`, missing needs), exit 1.

- [ ] **Step 2: Append Steps 8-9 to the skill**

Use the variants for `SALMON_ROUTE`, `DTU_FILTER_SCOPE` and `GATE_OUTCOME`; substitute `«KEY»`. Before writing the yaml block, compare each literal value with `D/tests/fixtures/config_params.txt`: the values below are the 1.0.4 defaults from the research report; where the recorded config differs, write the recorded value (the checker enforces it). [B] Add `isoformswitchanalyzer: false` and `leafcutter: false` after `sashimi_plot: false` when they are parameters in `rnasplice_schema.json`.

````markdown

---

## Step 8 — Analyses (modules) and their settings

In the pipeline's own configuration every analysis module is switched on, so a module that is not mentioned in the params file runs anyway. This skill therefore writes every switch explicitly.

Ask (numbered; several can be chosen, for example "1, 3"): "Which analyses should run?"
1. **rMATS** — differential splicing events (skipped exon, alternative 5' and 3' splice sites, mutually exclusive exons, retained intron) from junction reads, for each contrast. Recommended; usual cut-offs are FDR < 0.05 and |IncLevelDifference| > 0.1 (Akerberg et al. 2022 used these with rMATS).
2. **SUPPA2** — PSI of local events (SE, SS, MX, RI, FL) and of isoforms from Salmon transcript abundance, with differential splicing between conditions. FASTQ input only.
3. **DEXSeq exon usage** — differential usage of exon bins within a gene (exon counts relative to the gene).
4. **edgeR exon usage** — the same question with edgeR on featureCounts exon counts (edgeR function in this revision: «EDGER_DEU_FUNCTION»).
5. **DEXSeq transcript usage (DTU)** — changes in the share of each transcript within its gene (DRIMSeq filter, DEXSeq, stageR), from Salmon. FASTQ input only. ⚠️ `/bulk-rnaseq-pipeline` has the same DTU workflow (DRIMSeq → DEXSeq → stageR) on the Salmon output of an nf-core/rnaseq run; if it ran (or will run) there for these samples, do not choose it here: it is the same analysis twice.

Default (empty answer): 1 only. At least one analysis must be chosen. With BAM input, options 2 and 5 are not offered (they need Salmon quantification from reads). Set `{RUN_RMATS}`, `{RUN_SUPPA}`, `{RUN_DEXSEQ_EXON}`, `{RUN_EDGER_EXON}` and `{RUN_DEXSEQ_DTU}` to `true` for the chosen analyses and `false` for all others.

**MISO is not used.** The pipeline can draw sashimi plots with MISO (`sashimi_plot`), but MISO is archived (Python 2, no changes since 2019), and in rnasplice it only draws plots for a list of genes (by default three human Ensembl IDs); it computes no PSI or Bayes factors here. This skill always writes `sashimi_plot: false`. For sashimi plots, use rmats2sashimiplot or ggsashimi on the BAMs afterwards.

**rMATS settings** (asked only when rMATS is chosen; otherwise the values below are written unchanged):
- Novel splice sites (numbered): 1. Annotated splice sites only (default) · 2. Also detect unannotated splice sites. `{RMATS_NOVEL}` = `false` or `true`. With 2, rMATS also uses `rmats_min_intron_len` (50) and `rmats_max_exon_len` (500), written at their defaults.
- `rmats_splice_diff_cutoff` stays at 0.0001: it is the threshold of rMATS's null hypothesis (an inclusion difference larger than this counts as differential), not a reporting cut-off; filter the results by FDR and |IncLevelDifference| afterwards.
- `rmats_read_len` = `{READ_LENGTH}` (Step 6); `rmats_paired_stats` = `{PAIRED_DESIGN}` (Step 5).

**DTU filter values** (always written; used when DTU runs). The pipeline's schema (6/0/0) and its configuration (4/2/2) disagree on `min_samps_gene_expr`, `min_samps_feature_expr` and `min_samps_feature_prop`, and neither fits every design. This skill uses the rule of Love et al. 2018 (as `/bulk-rnaseq-pipeline` does) and tells the user the values:
[all_samples] `{MIN_SAMPS_GENE_EXPR}` = the number of samples in the samplesheet (distinct names), because the filter is applied once to all samples.
[per_contrast] `{MIN_SAMPS_GENE_EXPR}` = the smallest number of samples in a contrast (treatment plus control), because the filter is applied per contrast.
`{MIN_SAMPS_FEATURE}` = the size of the smallest condition used in a contrast (for both `min_samps_feature_expr` and `min_samps_feature_prop`); `min_gene_expr` 10, `min_feature_expr` 10, `min_feature_prop` 0.1, `dtu_txi` `dtuScaledTPM` (the pipeline default).

**SUPPA2, DEXSeq and edgeR settings:** pipeline defaults, written explicitly; SUPPA2's paired test `diffsplice_paired` = `{PAIRED_DESIGN}` (the pipeline default `true` assumes paired samples).

**Salmon route.**
[pseudo_only] `aligner: "star"` and `pseudo_aligner: "salmon"`: the STAR alignments feed rMATS, DEXSeq and edgeR; Salmon, run on the reads, feeds DTU and SUPPA2 once. (The pipeline's own default, `star_salmon`, quantifies with Salmon twice and runs DTU and SUPPA2 on both.) Salmon runs even when neither DTU nor SUPPA2 is chosen; it is quick.
[star_salmon_only] `aligner: "star_salmon"` and `«PSEUDO_OFF_LINE»`: Salmon quantifies the STAR alignments and feeds DTU and SUPPA2 once; the separate Salmon run on the reads is switched off.
[star_salmon_both] `aligner: "star_salmon"` and `pseudo_aligner: "salmon"`: this pipeline revision always runs both Salmon branches, so DTU and SUPPA2 results appear twice; use the ones under the `star_salmon` directory (Step 12 says so).

**Module and option keys (written to the params file).** Step 11 inserts this block, with its placeholders filled, into the params file:
```yaml
# alignment and quantification
aligner: "star"
pseudo_aligner: "salmon"
# analysis modules: every switch explicit (the pipeline default is true for all of them)
rmats: {RUN_RMATS}
dexseq_exon: {RUN_DEXSEQ_EXON}
edger_exon: {RUN_EDGER_EXON}
dexseq_dtu: {RUN_DEXSEQ_DTU}
suppa: {RUN_SUPPA}
sashimi_plot: false
# rMATS
rmats_read_len: {READ_LENGTH}
rmats_paired_stats: {PAIRED_DESIGN}
rmats_splice_diff_cutoff: 0.0001
rmats_novel_splice_site: {RMATS_NOVEL}
rmats_min_intron_len: 50
rmats_max_exon_len: 500
# DEXSeq exon usage
alignment_quality: 10
aggregation: true
save_dexseq_annotation: false
save_dexseq_plot: true
n_dexseq_plot: 10
# edgeR exon usage
save_edger_plot: true
n_edger_plot: 10
# DEXSeq transcript usage (DRIMSeq filter, DEXSeq, stageR)
dtu_txi: "dtuScaledTPM"
min_samps_gene_expr: {MIN_SAMPS_GENE_EXPR}
min_samps_feature_expr: {MIN_SAMPS_FEATURE}
min_samps_feature_prop: {MIN_SAMPS_FEATURE}
min_gene_expr: 10
min_feature_expr: 10
min_feature_prop: 0.1
# SUPPA2
suppa_per_local_event: true
suppa_per_isoform: true
generateevents_pool_genes: true
generateevents_event_type: "SE SS MX RI FL"
generateevents_boundary: "S"
generateevents_threshold: 10
generateevents_exon_length: 100
psiperevent_total_filter: 0
diffsplice_local_event: true
diffsplice_isoform: true
diffsplice_method: "empirical"
diffsplice_area: 1000
diffsplice_lower_bound: 0
diffsplice_gene_correction: true
diffsplice_paired: {PAIRED_DESIGN}
diffsplice_alpha: 0.05
diffsplice_median: false
diffsplice_tpm_threshold: 0
diffsplice_nan_threshold: 0
clusterevents_local_event: true
clusterevents_isoform: true
clusterevents_dpsithreshold: 0.05
clusterevents_eps: 0.05
clusterevents_metric: "euclidean"
clusterevents_min_pts: 20
clusterevents_method: "DBSCAN"
```
(The first two key lines are the [pseudo_only] route; for [star_salmon_only] they are `aligner: "star_salmon"` and `«PSEUDO_OFF_LINE»`; for [star_salmon_both] `aligner: "star_salmon"` and `pseudo_aligner: "salmon"`.)

[A/C] Never written, because they are not parameters of this revision: `isoformswitchanalyzer`, `leafcutter`, `ignore_tx_version`, `rmats_variable_read_len` (rMATS always uses variable read lengths here); the usage page's `local_events` is not a parameter either (the event types are `generateevents_event_type`). No iGenomes `genome` key: FASTA and GTF are always given.
[B] Never written: `rmats_variable_read_len` and `local_events` (not parameters); no iGenomes `genome` key. `isoformswitchanalyzer` and `leafcutter` are parameters of this development revision and are written `false` (this skill does not offer them).

---

## Step 9 — Trimming and QC

Trimming (Trim Galore) and QC use the pipeline defaults; this skill does not change them and writes no trimming key. If the user wants other trimming settings, they edit the params file after Step 11 and should know that the run is then not the configuration this skill verified.
````

(The parenthesis after the yaml block is plan guidance, not skill text: write only the route's two lines in the block.)

- [ ] **Step 3: Run the checker; expected PASS**

Run: `bash nfcore-rnasplice-setup/tests/check_skill.sh nfcore-rnasplice-setup/nfcore-rnasplice-setup.md` → `PASS`. A `differs from the pipeline config default` failure means the yaml value is wrong: take the recorded value; never add a key to `DEV` without a written reason in this step's report.

- [ ] **Step 4: Mutation rows, RED proof, commit**

```
T5-switch	check	skill	/^edger_exon: {RUN_EDGER_EXON}$/d	module switch edger_exon must be written exactly once
T5-dup	check	skill	s/^suppa: {RUN_SUPPA}$/suppa: {RUN_SUPPA}\nsuppa: true/	module switch suppa must be written exactly once
T5-sashimi	check	skill	s/^sashimi_plot: false$/sashimi_plot: true/	sashimi_plot must be written as false
T5-readlen	check	skill	s/^rmats_read_len: {READ_LENGTH}$/rmats_read_len: 40/	module block must contain the line: rmats_read_len: {READ_LENGTH}
T5-paired	check	skill	s/^rmats_paired_stats: {PAIRED_DESIGN}$/rmats_paired_stats: true/	module block must contain the line: rmats_paired_stats: {PAIRED_DESIGN}
T5-suppa-paired	check	skill	s/^diffsplice_paired: {PAIRED_DESIGN}$/diffsplice_paired: true/	module block must contain the line: diffsplice_paired: {PAIRED_DESIGN}
T5-dtu	check	skill	/^min_samps_feature_prop: /d	module block must contain the line: min_samps_feature_prop: {MIN_SAMPS_FEATURE}
T5-default	check	skill	s/^n_dexseq_plot: 10$/n_dexseq_plot: 11/	key n_dexseq_plot = 11 differs from the pipeline config default
T5-route	check	skill	/^pseudo_aligner: /d	Salmon route
```
[A/C only] also:
```
T5-never	check	skill	s/^sashimi_plot: false$/sashimi_plot: false\nleafcutter: false/	params key not in the recorded schema: leafcutter
```
(For route `star_salmon_only` the `T5-route` program is `/^@PSEUDO_OFF_LINE@$/d`.) Then:
```bash
bash nfcore-rnasplice-setup/tests/prove_mutations.sh nfcore-rnasplice-setup/nfcore-rnasplice-setup.md
bash nfcore-rnasplice-setup/tests/prove_red.sh <BASE>
git status && git diff --stat && git diff
git add nfcore-rnasplice-setup
git commit -m "nfcore-rnasplice-setup: Steps 8-9 explicit module switches, rMATS/DTU/SUPPA settings, parameter policy"
```

---

### Task 6: Config, params file, submission script, helpers (Steps 10-11)

**Files:**
- Create: `D/tests/render_params.sh`, `D/tests/test_render_params.sh`, `D/tests/dry_run_submit.sh`, `D/tests/dry_run_helpers.sh`
- Modify: `D/tests/check_skill.sh`, `D/tests/mutations.tsv`, `D/nfcore-rnasplice-setup.md` (append Steps 10-11)

**Interfaces:**
- Consumes: every placeholder above; the module block (Task 5); `defaults_equal`, `DEV` (Task 5); `test_env.sh` (Task 2); gate values `NEXTFLOW_MIN`, `NEXTFLOW_MAX_EXCL`, `NEXTFLOW_TESTED`, `HAS_MAX_PARAMS`, `HAS_RESOURCE_LIMITS`, `GATE_OUTCOME`, `VERSION_TAG`; fixture `trace_process_names.txt`.
- Produces: `{MULTIQC_TITLE}`, `{OUTDIR}`, `{PARAMS_YAML}` = `{SHEET_PREFIX}_params.yaml`; files written for the user: `nextflow.config`, `{PARAMS_YAML}`, `nf-core_rnasplice_{VERSION_TAG}.sh`, `download_genome_{REF_TAG}.sh` (only when Step 7 found a missing Ensembl file), [C] `create_nextflow_env_«NEXTFLOW_TESTED».sh`; anchors `**nextflow.config template.**`, `**Params file template.**`, `**Submission script.**`, `**Genome download helper.**`, [C] `**Nextflow environment helper.**`; `tests/render_params.sh <skill> <values.tsv> <out.yaml>` (values `NAME<TAB>VALUE`, prints `RENDERED <out>`), used by Task 9.

**Selector table** (Planner decision 17). For each row, take every distinct last `:`-component of the names in `D/tests/fixtures/trace_process_names.txt` that matches the row's pattern (ERE `^<row name>` on the last component) and write one `withName: '.*:<component>'` block with the row's resources; rows with no match are left out.

| Row | cpus | memory | time |
|---|---|---|---|
| STAR_GENOMEGENERATE | 8 | '64 GB' | '8h' |
| STAR_ALIGN | 8 | '48 GB' | '8h' |
| RMATS_PREP | 4 | '16 GB' | '8h' |
| RMATS_POST | 8 | '32 GB' | '16h' |
| DEXSEQ_COUNT | 2 | '8 GB' | '8h' |
| DEXSEQ_EXON | 8 | '32 GB' | '8h' |
| DEXSEQ_DTU | 8 | '32 GB' | '8h' |
| SALMON_QUANT | 8 | '16 GB' | '4h' |

- [ ] **Step 1: Write render_params.sh and the three failing tests**

Create `D/tests/render_params.sh`:
```bash
#!/bin/bash
# Usage: render_params.sh <skill.md> <values.tsv> <out.yaml>
# Builds the params file exactly as Step 11 describes: the Step 11 template (block after "**Params file template.**") with its
# {MODULE_PARAMS} line replaced by the Step 8 block (after "**Module and option keys (written to the params file).**"), then
# every {NAME} replaced by its value from values.tsv (NAME<TAB>VALUE). The star_index and salmon_index lines are deleted when
# their value is empty; any other missing or empty value, a value containing {...}, or a repeated key is an error, and then
# no output file is written. Prints "RENDERED <out.yaml>" or exits 1.
set -u
SKILL=${1:?usage: render_params.sh <skill.md> <values.tsv> <out.yaml>}
VALS=${2:?usage: render_params.sh <skill.md> <values.tsv> <out.yaml>}
OUT=${3:?usage: render_params.sh <skill.md> <values.tsv> <out.yaml>}
HERE=$(cd "$(dirname "$0")" && pwd)
T=$(mktemp -d) || exit 1
case "$T" in /tmp/tmp.?*) ;; *) echo "unexpected temp dir: $T"; exit 1 ;; esac
trap 'rm -f "$T/base" "$T/mod" "$T/joined" "$T/out"; rmdir "$T"' EXIT
bash "$HERE/cut_block.sh" "$SKILL" "**Params file template.**" > "$T/base" || { echo "ERROR: params file template not found"; exit 1; }
bash "$HERE/cut_block.sh" "$SKILL" "**Module and option keys (written to the params file).**" > "$T/mod" || { echo "ERROR: module block not found"; exit 1; }
[ "$(grep -cx '{MODULE_PARAMS}' "$T/base")" -eq 1 ] || { echo "ERROR: the template must contain exactly one {MODULE_PARAMS} line"; exit 1; }
awk -v mod="$T/mod" '$0 == "{MODULE_PARAMS}" { while ((getline l < mod) > 0) print l; next } { print }' "$T/base" > "$T/joined"
awk -F'\t' '
  NR == FNR {
    if ($1 != "") { has[$1] = 1; v[$1] = $2; if ($2 ~ /\{[A-Z_]+\}/) { print "ERROR: value of " $1 " contains a placeholder" > "/dev/stderr"; bad = 1 } }
    next
  }
  /^(star_index|salmon_index): "\{(STAR_INDEX|SALMON_INDEX)\}"$/ { k = $0; sub(/^[^{]*\{/, "", k); sub(/\}.*$/, "", k); if (!(k in has) || v[k] == "") next }
  { line = $0
    while (match(line, /\{[A-Z_]+\}/)) {
      k = substr(line, RSTART + 1, RLENGTH - 2)
      if (!(k in has) || v[k] == "") { print "ERROR: no value for {" k "}" > "/dev/stderr"; bad = 1; line = substr(line, 1, RSTART - 1) substr(line, RSTART + RLENGTH); continue }
      line = substr(line, 1, RSTART - 1) v[k] substr(line, RSTART + RLENGTH)
    }
    if (match(line, /^[a-z_0-9]+:/)) { key = substr(line, 1, RLENGTH - 1); if (key in seen) { print "ERROR: repeated key " key > "/dev/stderr"; bad = 1 }; seen[key] = 1 }
    print line
  }
  END { exit bad }' "$VALS" "$T/joined" > "$T/out" || exit 1
mv "$T/out" "$OUT" && echo "RENDERED $OUT"
```

Create `D/tests/test_render_params.sh`:
```bash
#!/bin/bash
# Usage: test_render_params.sh <skill.md>
# Renders the params file with render_params.sh from test values and checks the result and the error cases.
# Prints "RENDER PASS" or exits 1.
set -u
source "$(cd "$(dirname "$0")" && pwd)/test_env.sh" "$@"
SKILL=${1:?usage: test_render_params.sh <skill.md>}
HERE=$(cd "$(dirname "$0")" && pwd)
gv() { awk -F'\t' -v k="$1" '$1 == k {print $2; exit}' "$HERE/fixtures/gate_values.tsv"; }
V="$TEST_TMP/values.tsv"; O="$TEST_TMP/out.yaml"; RL=$(gv TEST_READ_LENGTH); fail=0
values() { printf '%s\t%s\n' VERSION "$(gv PIPELINE_REVISION)" SAMPLESHEET_CSV test_samplesheet.csv CONTRASTS_CSV test_contrasts.csv \
  SOURCE fastq OUTDIR results/2026-10-01_test MULTIQC_TITLE 261001 FASTA_PATH /data/genome.fa GTF_PATH /data/genes.gtf \
  STAR_INDEX "" SALMON_INDEX /data/index/salmon RUN_RMATS true RUN_DEXSEQ_EXON false RUN_EDGER_EXON false RUN_DEXSEQ_DTU true \
  RUN_SUPPA false READ_LENGTH "$RL" PAIRED_DESIGN false RMATS_NOVEL false MIN_SAMPS_GENE_EXPR 4 MIN_SAMPS_FEATURE 2 > "$V"; }
render() { rm -f "$O"; bash "$HERE/render_params.sh" "$SKILL" "$V" "$O" > "$TEST_TMP/log" 2>&1; }
values
if render; then
  grep -q '{' "$O" && { echo "FAIL: case full: placeholder left"; fail=1; }
  grep -q '^star_index:' "$O" && { echo "FAIL: case full: the empty star_index line was not deleted"; fail=1; }
  for l in "# nf-core/rnasplice $(gv PIPELINE_REVISION) parameters — generated by /nfcore-rnasplice-setup" \
           'salmon_index: "/data/index/salmon"' 'multiqc_title: "261001"' 'source: "fastq"' "rmats_read_len: $RL" \
           'rmats: true' 'dexseq_exon: false' 'edger_exon: false' 'dexseq_dtu: true' 'suppa: false' 'sashimi_plot: false' \
           'rmats_paired_stats: false' 'diffsplice_paired: false' 'min_samps_gene_expr: 4' 'min_samps_feature_expr: 2' 'min_samps_feature_prop: 2'; do
    grep -qxF -- "$l" "$O" || { echo "FAIL: case full: line missing: $l"; fail=1; }
  done
  [ "$(gv HAS_MAX_PARAMS)" = yes ] && { grep -qxF 'max_memory: "64.GB"' "$O" || { echo "FAIL: case full: max_memory line missing"; fail=1; }; }
  [ -z "$(grep -oE '^[a-z_0-9]+:' "$O" | sort | uniq -d)" ] || { echo "FAIL: case full: repeated key"; fail=1; }
else echo "FAIL: case full: render failed: $(cat "$TEST_TMP/log")"; fail=1; fi
values; sed -i 's#^STAR_INDEX\t$#STAR_INDEX\t/data/index/star#' "$V"
render && grep -qxF 'star_index: "/data/index/star"' "$O" || { echo "FAIL: case star_index_set"; fail=1; }
values; sed -i '/^RUN_SUPPA\t/d' "$V"
if render || [ -e "$O" ] || ! grep -qF 'no value for {RUN_SUPPA}' "$TEST_TMP/log"; then echo "FAIL: case missing_value"; fail=1; fi
values; sed -i 's#^OUTDIR\t.*#OUTDIR\t#' "$V"
if render || [ -e "$O" ] || ! grep -qF 'no value for {OUTDIR}' "$TEST_TMP/log"; then echo "FAIL: case empty_required_value"; fail=1; fi
values; sed -i 's#^OUTDIR\t.*#OUTDIR\t{X}#' "$V"
if render || [ -e "$O" ] || ! grep -qF 'contains a placeholder' "$TEST_TMP/log"; then echo "FAIL: case placeholder_in_value"; fail=1; fi
forbidden_ran && fail=1
[ $fail -eq 0 ] && echo "RENDER PASS" || exit 1
```

Create `D/tests/dry_run_submit.sh`:
```bash
#!/bin/bash
# Usage: dry_run_submit.sh <skill.md>
# Cuts the Step 11 submission script (block after "**Submission script.**"), fills its placeholders as the wizard would, and
# runs it in scenarios against stubs of module, conda, singularity and nextflow (the real condainit is never sourced).
# Prints "DRY RUN SUBMIT PASS" or exits 1.
set -u
source "$(cd "$(dirname "$0")" && pwd)/test_env.sh" "$@"
SKILL=${1:?usage: dry_run_submit.sh <skill.md>}
HERE=$(cd "$(dirname "$0")" && pwd)
gv() { awk -F'\t' -v k="$1" '$1 == k {print $2; exit}' "$HERE/fixtures/gate_values.tsv"; }
REV=$(gv PIPELINE_REVISION); TAG=$(gv VERSION_TAG); NF_T=$(gv NEXTFLOW_TESTED); NF_MIN=$(gv NEXTFLOW_MIN); NF_MAX=$(gv NEXTFLOW_MAX_EXCL)
S="$TEST_TMP/scen"; R="$TEST_TMP/run"; mkdir -p "$S" "$R" "$TEST_TMP/nonf" || exit 1
export STUB_STATE="$TEST_TMP/state"; mkdir -p "$STUB_STATE" || exit 1
cat > "$S/module" <<'EOF'
#!/bin/bash
echo "module $*" >> "$STUB_STATE/calls"
grep -qxF -- "$2" "$STUB_STATE/module_fail" 2>/dev/null && exit 1
exit 0
EOF
cat > "$S/conda" <<'EOF'
#!/bin/bash
echo "conda $*" >> "$STUB_STATE/calls"
[ -e "$STUB_STATE/conda_fail" ] && exit 1
exit 0
EOF
cat > "$S/singularity" <<'EOF'
#!/bin/bash
echo "singularity $*" >> "$STUB_STATE/calls"; exit 0
EOF
cat > "$S/nextflow" <<'EOF'
#!/bin/bash
case "$1" in
  -version) printf '\n      N E X T F L O W\n      version %s build 5940\n      created 01-01-2026 00:00 UTC\n' "$(cat "$STUB_STATE/nf_version")" ;;
  pull) echo "nextflow $* NXF_OFFLINE=${NXF_OFFLINE:-unset}" >> "$STUB_STATE/calls"; [ -e "$STUB_STATE/pull_fail" ] && exit 1 ;;
  *) echo "nextflow $* NXF_OFFLINE=${NXF_OFFLINE:-unset}" >> "$STUB_STATE/calls" ;;
esac
exit 0
EOF
: > "$S/condainit"
chmod +x "$S/module" "$S/conda" "$S/singularity" "$S/nextflow"
cp "$S/module" "$S/conda" "$S/singularity" "$TEST_TMP/nonf/"
for c in module conda singularity nextflow; do [ "$(PATH="$S:$PATH" type -P "$c")" = "$S/$c" ] || { echo "stub $c does not resolve to $S/$c"; exit 1; }; done
[ -z "$(PATH="$TEST_TMP/nonf:/usr/bin:/bin" type -P nextflow)" ] || { echo "a real nextflow is on /usr/bin:/bin; the missing-nextflow case cannot be tested"; exit 1; }
bash "$HERE/cut_block.sh" "$SKILL" "**Submission script.**" > "$TEST_TMP/tmpl.sh" || { echo "FAIL: submission script block not found"; exit 1; }
[ "$(grep -c '^source /home/software/conda/miniconda3/bin/condainit ' "$TEST_TMP/tmpl.sh")" -eq 1 ] || { echo "FAIL: expected exactly one condainit source line"; exit 1; }
[ "$(grep -c '^[[:space:]]*source ' "$TEST_TMP/tmpl.sh")" -eq 1 ] || { echo "FAIL: unexpected extra source line"; exit 1; }
sed -e "s#^source /home/software/conda/miniconda3/bin/condainit #source \"$S/condainit\" #" \
    -e "s#{USER_EMAIL}#test@example.org#g; s#{CONDA_ENV}#nf-test#g; s#{VERSION_TAG}#$TAG#g; s#{VERSION}#$REV#g; s#{PARAMS_YAML}#T_params.yaml#g" \
    "$TEST_TMP/tmpl.sh" > "$R/submit.sh"
! grep -nE '(^|[^$])\{[A-Z_]+\}' "$R/submit.sh" || { echo "FAIL: placeholder left in the submission script"; exit 1; }
! grep -q '/home/software/conda' "$R/submit.sh" || { echo "FAIL: the real condainit would be sourced"; exit 1; }
fail=0
scenario() { # scenario <name> <nextflow version> <expected exit 0|1> <expect one run: yes|no> [module_fail=<module> conda_fail pull_fail no_nextflow]
  local name=$1 ver=$2 erc=$3 erun=$4 p="$S:$PATH" w rc runs; shift 4
  rm -f "$STUB_STATE/calls" "$STUB_STATE/module_fail" "$STUB_STATE/conda_fail" "$STUB_STATE/pull_fail"
  printf '%s\n' "$ver" > "$STUB_STATE/nf_version"
  for w in "$@"; do case "$w" in
    module_fail=*) echo "${w#module_fail=}" >> "$STUB_STATE/module_fail" ;;
    conda_fail) : > "$STUB_STATE/conda_fail" ;;
    pull_fail) : > "$STUB_STATE/pull_fail" ;;
    no_nextflow) p="$TEST_TMP/nonf:/usr/bin:/bin" ;;
  esac; done
  ( cd "$R" && NXF_OFFLINE=TRUE PATH="$p" bash ./submit.sh ) > "$TEST_TMP/out" 2>&1; rc=$?
  [ "$rc" -eq "$erc" ] || { echo "FAIL: case $name: exit $rc, expected $erc"; sed 's/^/    /' "$TEST_TMP/out"; fail=1; }
  runs=$(grep -c '^nextflow run ' "$STUB_STATE/calls" 2>/dev/null); runs=${runs:-0}
  if [ "$erun" = yes ]; then
    [ "$runs" -eq 1 ] && grep -qxF "nextflow run nf-core/rnasplice -r $REV -c nextflow.config -profile slurm,singularity -params-file T_params.yaml NXF_OFFLINE=TRUE" "$STUB_STATE/calls" \
      || { echo "FAIL: case $name: expected exactly one 'nextflow run' with the exact arguments"; cat "$STUB_STATE/calls" 2>/dev/null; fail=1; }
  else
    [ "$runs" -eq 0 ] || { echo "FAIL: case $name: nextflow run was called"; fail=1; }
  fi
}
scenario tested "$NF_T" 0 yes
grep -qxF "nextflow pull nf-core/rnasplice -r $REV NXF_OFFLINE=false" "$STUB_STATE/calls" || { echo "FAIL: case pull_offline_override: the pull must run with NXF_OFFLINE=false"; fail=1; }
awk '/^module add miniconda3/ {a = NR} /^conda activate/ {b = NR} /^module add singularity/ {c = NR} END {exit !(a && b && c && a < b && b < c)}' "$STUB_STATE/calls" \
  || { echo "FAIL: case order: module miniconda3, then conda activate, then module singularity"; fail=1; }
scenario below_min 0.0.1 1 no
grep -q "older than $NF_MIN" "$TEST_TMP/out" || { echo "FAIL: case below_min: message"; fail=1; }
if [ "$NF_MAX" = none ]; then scenario far_future 99.0.0 0 yes
else scenario at_max "$NF_MAX" 1 no; scenario far_future 99.0.0 1 no; fi
scenario pull_fails "$NF_T" 0 yes pull_fail
grep -q "WARNING: nextflow pull failed" "$TEST_TMP/out" || { echo "FAIL: case pull_fails: warning"; fail=1; }
scenario singularity_module_fails "$NF_T" 1 no module_fail=singularity/3.10.4
! grep -q '^nextflow' "$STUB_STATE/calls" || { echo "FAIL: case singularity_module_fails: nextflow was called"; fail=1; }
scenario miniconda_module_fails "$NF_T" 1 no module_fail=miniconda3/v4
! grep -q '^conda' "$STUB_STATE/calls" || { echo "FAIL: case miniconda_module_fails: conda was called"; fail=1; }
scenario conda_fails "$NF_T" 1 no conda_fail
scenario nextflow_missing "$NF_T" 1 no no_nextflow
grep -q "nextflow is not in conda environment" "$TEST_TMP/out" || { echo "FAIL: case nextflow_missing: message"; fail=1; }
forbidden_ran && fail=1
[ $fail -eq 0 ] && echo "DRY RUN SUBMIT PASS" || exit 1
```

Create `D/tests/dry_run_helpers.sh`:
```bash
#!/bin/bash
# Usage: dry_run_helpers.sh <skill.md>
# Runs the genome download helper (block after "**Genome download helper.**") against a wget stub, with the real gzip, in fresh,
# re-run, corrupt, truncated, failing and relative-path scenarios; for gate outcome C also the Nextflow environment helper
# against a conda stub. Prints "DRY RUN HELPERS PASS" or exits 1.
set -u
source "$(cd "$(dirname "$0")" && pwd)/test_env.sh" "$@"
SKILL=${1:?usage: dry_run_helpers.sh <skill.md>}
HERE=$(cd "$(dirname "$0")" && pwd)
gv() { awk -F'\t' -v k="$1" '$1 == k {print $2; exit}' "$HERE/fixtures/gate_values.tsv"; }
S="$TEST_TMP/scen"; mkdir -p "$S" || exit 1
export STUB_STATE="$TEST_TMP/state"; mkdir -p "$STUB_STATE" || exit 1
cat > "$S/wget" <<'EOF'
#!/bin/bash
# wget stub: wget -c -O <out> <url>; per-URL mode in $STUB_STATE/wget_mode ("<url> fail" or "<url> truncated")
out=""; url=""
while [ $# -gt 0 ]; do case "$1" in -O) out=$2; shift 2 ;; -*) shift ;; *) url=$1; shift ;; esac; done
echo "wget $url" >> "$STUB_STATE/calls"
mode=$(awk -v u="$url" '$1 == u {print $2}' "$STUB_STATE/wget_mode" 2>/dev/null)
[ "$mode" = fail ] && exit 4
case "$url" in
  *.fa.gz) printf '>1\nACGTACGT\n' ;;
  *.gtf.gz) printf '1\ttest\texon\t1\t8\t.\t+\t.\tgene_id "G1";\n' ;;
esac | gzip -c > "$out"
if [ "$mode" = truncated ]; then head -c 10 "$out" > "$out.t" && mv "$out.t" "$out"; fi
exit 0
EOF
chmod +x "$S/wget"; export PATH="$S:$PATH"
[ "$(type -P wget)" = "$S/wget" ] || { echo "stub wget does not resolve to $S/wget"; exit 1; }
case "$(type -P gzip)" in /usr/bin/gzip|/bin/gzip) ;; *) echo "gzip is not the system gzip"; exit 1 ;; esac
bash "$HERE/cut_block.sh" "$SKILL" "**Genome download helper.**" > "$TEST_TMP/helper.tmpl" || { echo "FAIL: download helper block not found"; exit 1; }
G="$TEST_TMP/genome"; FU=https://example.invalid/genome.fa.gz; GU=https://example.invalid/genes.gtf.gz
mk_helper() { # mk_helper <GENOME_DIR as written into the script>
  sed -e "s#{GENOME_DIR}#$1#g; s#{FASTA_PATH}#$G/genome.fa#g; s#{GTF_PATH}#$G/genes.gtf#g" \
      -e "s#{ENSEMBL_FASTA_URL}#$FU#g; s#{ENSEMBL_GTF_URL}#$GU#g; s#{USER_EMAIL}#test@example.org#g; s#{REF_TAG}#test#g" \
      "$TEST_TMP/helper.tmpl" > "$TEST_TMP/helper.sh"
  ! grep -nE '(^|[^$])\{[A-Z_]+\}' "$TEST_TMP/helper.sh" || { echo "FAIL: placeholder left in the helper"; exit 1; }
}
runh() { rm -f "$STUB_STATE/calls"; ( cd "$TEST_TMP" && bash ./helper.sh ) > "$TEST_TMP/out" 2>&1; }
nwget() { local n; n=$(grep -c '^wget ' "$STUB_STATE/calls" 2>/dev/null); echo "${n:-0}"; }
fail=0
mk_helper "$G"; mkdir -p "$G"; echo keep > "$G/sentinel"
runh; rc=$?
{ [ $rc -eq 0 ] && grep -qx '>1' "$G/genome.fa" && grep -q 'gene_id "G1"' "$G/genes.gtf" && [ -s "$G/genome.fa.gz" ] && [ -s "$G/genes.gtf.gz" ] \
  && [ "$(nwget)" -eq 2 ] && [ -z "$(ls "$G" | grep '\.part$')" ] && [ -s "$G/sentinel" ]; } || { echo "FAIL: case fresh"; cat "$TEST_TMP/out"; fail=1; }
runh; rc=$?; { [ $rc -eq 0 ] && [ "$(nwget)" -eq 0 ]; } || { echo "FAIL: case rerun"; fail=1; }
rm -f "$G/genome.fa"; runh; rc=$?; { [ $rc -eq 0 ] && [ "$(nwget)" -eq 0 ] && grep -qx '>1' "$G/genome.fa"; } || { echo "FAIL: case gz_present"; fail=1; }
rm -f "$G/genes.gtf"; printf 'garbage' > "$G/genes.gtf.gz"; runh; rc=$?
{ [ $rc -eq 0 ] && [ "$(nwget)" -eq 1 ] && grep -q 'gene_id "G1"' "$G/genes.gtf"; } || { echo "FAIL: case corrupt_gz"; fail=1; }
rm -f "$G/genes.gtf" "$G/genes.gtf.gz"; echo "$GU truncated" > "$STUB_STATE/wget_mode"; runh; rc=$?
{ [ $rc -ne 0 ] && grep -q 'is not a valid gzip file' "$TEST_TMP/out" && [ ! -e "$G/genes.gtf" ] && [ ! -e "$G/genes.gtf.gz" ] && [ ! -e "$G/genes.gtf.gz.part" ]; } \
  || { echo "FAIL: case truncated_download"; cat "$TEST_TMP/out"; fail=1; }
echo "$GU fail" > "$STUB_STATE/wget_mode"; runh; rc=$?
{ [ $rc -ne 0 ] && grep -q "download failed for $GU" "$TEST_TMP/out" && [ ! -e "$G/genes.gtf" ] && grep -qx '>1' "$G/genome.fa"; } || { echo "FAIL: case wget_fails"; fail=1; }
rm -f "$STUB_STATE/wget_mode"
mk_helper genome_rel; runh; rc=$?
{ [ $rc -ne 0 ] && grep -q 'must be an absolute path' "$TEST_TMP/out" && [ "$(nwget)" -eq 0 ] && [ ! -e "$TEST_TMP/genome_rel" ]; } || { echo "FAIL: case relative_dir"; fail=1; }
[ -s "$G/sentinel" ] || { echo "FAIL: a file the helper did not create was removed"; fail=1; }
if [ "$(gv GATE_OUTCOME)" = C ]; then
  NFV=$(gv NEXTFLOW_TESTED)
  bash "$HERE/cut_block.sh" "$SKILL" "**Nextflow environment helper.**" > "$TEST_TMP/env.tmpl" || { echo "FAIL: environment helper block not found"; exit 1; }
  sed -e "s#^source /home/software/conda/miniconda3/bin/condainit #source /dev/null #" -e "s#{USER_EMAIL}#test@example.org#g" "$TEST_TMP/env.tmpl" > "$TEST_TMP/env.sh"
  ! grep -q '/home/software/conda' "$TEST_TMP/env.sh" || { echo "FAIL: the real condainit would be sourced"; exit 1; }
  printf '#!/bin/bash\necho "module $*" >> "$STUB_STATE/calls"; exit 0\n' > "$S/module"
  printf '#!/bin/bash\necho "conda $*" >> "$STUB_STATE/calls"\n[ "$1" = env ] && cat "$STUB_STATE/envs" 2>/dev/null\nexit 0\n' > "$S/conda"
  chmod +x "$S/module" "$S/conda"
  for c in module conda; do [ "$(type -P "$c")" = "$S/$c" ] || { echo "stub $c does not resolve"; exit 1; }; done
  rm -f "$STUB_STATE/calls" "$STUB_STATE/envs"; ( cd "$TEST_TMP" && bash ./env.sh ) > "$TEST_TMP/out" 2>&1
  grep -qF "conda create -y -n nf-$NFV -c conda-forge -c bioconda nextflow=$NFV" "$STUB_STATE/calls" || { echo "FAIL: case env_create"; fail=1; }
  printf 'base  /x\nnf-%s  /y\n' "$NFV" > "$STUB_STATE/envs"; rm -f "$STUB_STATE/calls"; ( cd "$TEST_TMP" && bash ./env.sh ) > "$TEST_TMP/out" 2>&1
  ! grep -q "conda create" "$STUB_STATE/calls" || { echo "FAIL: case env_exists: conda create ran again"; fail=1; }
fi
forbidden_ran && fail=1
[ $fail -eq 0 ] && echo "DRY RUN HELPERS PASS" || exit 1
```

Run the three tests on the current skill. Expected: each prints `FAIL: ... block not found` (or `ERROR: params file template not found`) and exits 1.

- [ ] **Step 2: Append the Task 6 checker block; see it fail**

```bash
# --- Task 6 (Steps 10-11)
need "## Step 10 — MultiQC title, output directory, nextflow.config"
need "## Step 11 — Params file, submission script and helper scripts"
for a in "**nextflow.config template.**" "**Params file template.**" "**Submission script.**" "**Genome download helper.**"; do anchor_once "$a"; done
need "1. \`results/{TODAY_ISO}_{WD_NAME}\` (default)"
need "If it exists, do not overwrite it"
need "\`multiqc_title\` is always double-quoted"
need "{SHEET_PREFIX}_params.yaml"
need "sbatch --dependency=afterok:"
need "#SBATCH -t 48:00:00"
SUB=$(bash "$CUT" "$SKILL" "**Submission script.**" 2>/dev/null)
launch=$(printf '%s\n' "$SUB" | awk '/^nextflow run nf-core\/rnasplice/ {p = 1} p {l = l $0; if ($0 !~ /\\[ \t]*$/) {print l; exit}}' | sed 's/\\[ \t]*/ /g')
[ "$launch" = "nextflow run nf-core/rnasplice -r {VERSION} -c nextflow.config -profile slurm,singularity -params-file {PARAMS_YAML}" ] \
  || { echo "FAIL: launch line must be exactly 'nextflow run nf-core/rnasplice -r {VERSION} -c nextflow.config -profile slurm,singularity -params-file {PARAMS_YAML}', found: '$launch'"; fail=1; }
all_nf=$(grep -cE '(^|[;&|(] *)(NXF_[A-Z_]+=[^ ]+ +)?nextflow +(run|pull|-version)' "$SKILL")
sub_nf=$(printf '%s\n' "$SUB" | grep -cE '(^|[;&|(] *)(NXF_[A-Z_]+=[^ ]+ +)?nextflow +(run|pull|-version)')
{ [ "$all_nf" -eq "$sub_nf" ] && [ "$sub_nf" -ge 3 ]; } \
  || { echo "FAIL: nextflow is invoked outside the submission script ($all_nf lines in the skill, $sub_nf in the script); the wizard must never run nextflow on the login node"; fail=1; }
while IFS= read -r l; do
  case "$l" in *"||"*) ;; *) echo "FAIL: unguarded module load (a failed module add silently breaks later ones): $l"; fail=1 ;; esac
done < <(grep -E '^[ \t]*module (add|load) ' "$SKILL")
printf '%s\n' "$SUB" | grep -qxF 'NXF_OFFLINE=false nextflow pull nf-core/rnasplice -r {VERSION} \' || { echo "FAIL: the submission script must pull the pipeline with NXF_OFFLINE=false"; fail=1; }
printf '%s\n' "$SUB" | grep -qxF "NF_MIN=\"$(gv NEXTFLOW_MIN)\"; NF_MAX_EXCL=\"$(gv NEXTFLOW_MAX_EXCL)\"" || { echo "FAIL: the submission script must carry the verified Nextflow range"; fail=1; }
CFG=$(bash "$CUT" "$SKILL" "**nextflow.config template.**" 2>/dev/null)
sels=$(printf '%s\n' "$CFG" | sed -n "s/^[ \t]*withName: '\.\*:\([A-Z0-9_]*\)' {$/\1/p")
[ -n "$sels" ] || { echo "FAIL: no withName selectors in the nextflow.config template"; fail=1; }
for s in $sels; do grep -qE "(^|:)$s$" "$FIX/trace_process_names.txt" || { echo "FAIL: selector $s matches no process of the verification runs"; fail=1; }; done
for pat in STAR_GENOMEGENERATE STAR_ALIGN RMATS_PREP RMATS_POST DEXSEQ_COUNT DEXSEQ_EXON DEXSEQ_DTU SALMON_QUANT; do
  for comp in $(sed 's/.*://' "$FIX/trace_process_names.txt" | grep -E "^$pat" | sort -u); do
    echo "$sels" | grep -qx -- "$comp" || { echo "FAIL: process $comp ran in the verification runs but has no withName selector"; fail=1; }
  done
done
[ "$(printf '%s\n' "$CFG" | grep -cE '^(timeline|report|trace|dag) +\{ enabled = true; overwrite = true;')" -eq 4 ] || { echo "FAIL: overwrite = true must be set for timeline, report, trace and dag"; fail=1; }
printf '%s\n' "$CFG" | grep -qF "fields = 'task_id,hash,native_id,name,status,exit,cpus,memory,time,realtime,peak_rss'" || { echo "FAIL: trace fields must include the requested cpus, memory and time"; fail=1; }
if [ "$(gv HAS_RESOURCE_LIMITS)" = yes ]; then
  printf '%s\n' "$CFG" | grep -qF "resourceLimits = [ cpus: 16, memory: '64 GB', time: '24h' ]" || { echo "FAIL: the resourceLimits line is missing"; fail=1; }
else
  ! printf '%s\n' "$CFG" | grep -q resourceLimits || { echo "FAIL: resourceLimits is not known to Nextflow $(gv NEXTFLOW_TESTED)"; fail=1; }
fi
BASEY=$(bash "$CUT" "$SKILL" "**Params file template.**" 2>/dev/null)
[ "$(printf '%s\n' "$BASEY" | grep -cx '{MODULE_PARAMS}')" -eq 1 ] || { echo "FAIL: the params template must contain one {MODULE_PARAMS} line"; fail=1; }
for kv in 'input: "{SAMPLESHEET_CSV}"' 'contrasts: "{CONTRASTS_CSV}"' 'source: "{SOURCE}"' 'outdir: "{OUTDIR}"' 'multiqc_title: "{MULTIQC_TITLE}"' \
          'fasta: "{FASTA_PATH}"' 'gtf: "{GTF_PATH}"' 'star_index: "{STAR_INDEX}"' 'salmon_index: "{SALMON_INDEX}"' 'gencode: false' 'save_reference: false'; do
  printf '%s\n' "$BASEY" | grep -qxF -- "$kv" || { echo "FAIL: the params template must contain the line: $kv"; fail=1; }
done
if [ "$(gv HAS_MAX_PARAMS)" = yes ]; then
  for kv in 'max_cpus: 16' 'max_memory: "64.GB"' 'max_time: "24.h"'; do printf '%s\n' "$BASEY" | grep -qxF -- "$kv" || { echo "FAIL: the params template must contain the line: $kv"; fail=1; }; done
fi
defaults_equal "$BASEY"
if [ "$(gv GATE_OUTCOME)" = C ]; then anchor_once "**Nextflow environment helper.**"; need "create_nextflow_env_$(gv NEXTFLOW_TESTED).sh"; fi
# --- end Task 6
```
Run the checker: expected Task 6 `FAIL` lines, exit 1.

- [ ] **Step 3: Append Steps 10-11 to the skill**

The `withName` blocks of the config template are the rows of the selector table that match the fixture (rule above); the example below shows all eight rows as they would appear if every name matched exactly. [HAS_RESOURCE_LIMITS = no]: leave out the `process { resourceLimits ... }` block and its sentence. [HAS_MAX_PARAMS = no]: leave out the three `max_` lines of the params template and the sentence about them. Substitute `«KEY»`.

````markdown

---

## Step 10 — MultiQC title, output directory, nextflow.config

**MultiQC title** (numbered): 1. `{SEQ_DATE}_{WD_NAME}` · 2. `{TODAY_YYMMDD}_{WD_NAME}` · 3. Custom. Store `{MULTIQC_TITLE}`.
**Output directory** (numbered): 1. `results/{TODAY_ISO}_{WD_NAME}` (default) · 2. `results/{TODAY_ISO}_{SHEET_PREFIX}` · 3. Custom. Store `{OUTDIR}`. If it exists and is not empty, ask (numbered): 1. use it anyway (the pipeline adds to it and may overwrite files) · 2. choose another.

Check for an existing config: `ls nextflow.config`. If it exists, do not overwrite it: tell the user which selectors matter for rnasplice (the `withName` lines below) so they can compare. If it does not exist, write it from the template.

**nextflow.config template.**
```nextflow
// nextflow.config — nf-core/rnasplice on SLURM + Singularity (generated by /nfcore-rnasplice-setup)
profiles {
    slurm {
        process {
            executor = 'slurm'
            queue = 'bcc'
            cpus = 2
            memory = '8 GB'
            time = '4h'

            withName: '.*:STAR_GENOMEGENERATE' {
                cpus = 8
                memory = '64 GB'
                time = '8h'
            }
            withName: '.*:STAR_ALIGN' {
                cpus = 8
                memory = '48 GB'
                time = '8h'
            }
            withName: '.*:RMATS_PREP' {
                cpus = 4
                memory = '16 GB'
                time = '8h'
            }
            withName: '.*:RMATS_POST' {
                cpus = 8
                memory = '32 GB'
                time = '16h'
            }
            withName: '.*:DEXSEQ_COUNT' {
                cpus = 2
                memory = '8 GB'
                time = '8h'
            }
            withName: '.*:DEXSEQ_EXON' {
                cpus = 8
                memory = '32 GB'
                time = '8h'
            }
            withName: '.*:DEXSEQ_DTU' {
                cpus = 8
                memory = '32 GB'
                time = '8h'
            }
            withName: '.*:SALMON_QUANT' {
                cpus = 8
                memory = '16 GB'
                time = '4h'
            }
        }
        executor {
            queueSize = 10
            submitRateLimit = '10/1min'
            pollInterval = '30s'
        }
    }
    singularity {
        singularity {
            enabled = true
            autoMounts = true
        }
    }
}

process {
    resourceLimits = [ cpus: 16, memory: '64 GB', time: '24h' ]
}

timeline { enabled = true; overwrite = true; file = "${params.outdir}/pipeline_info/execution_timeline.html" }
report   { enabled = true; overwrite = true; file = "${params.outdir}/pipeline_info/execution_report.html"   }
trace    { enabled = true; overwrite = true; file = "${params.outdir}/pipeline_info/execution_trace.txt"; fields = 'task_id,hash,native_id,name,status,exit,cpus,memory,time,realtime,peak_rss' }
dag      { enabled = true; overwrite = true; file = "${params.outdir}/pipeline_info/pipeline_dag.svg"        }
```

The selectors (`'.*:NAME'`, matching the process whatever its workflow prefix) use the process names of this skill's verification run of nf-core/rnasplice «VERSION_TAG». After the first real run, compare them with `{OUTDIR}/pipeline_info/execution_trace.txt`, whose `cpus`, `memory` and `time` columns show what each task requested. `resourceLimits` caps every task at 16 CPUs, 64 GB and 24 h; the pipeline also caps tasks through `max_cpus`, `max_memory` and `max_time`, which the params file sets to the same limits. `overwrite = true` on the four report files lets a run that failed early be started again.

---

## Step 11 — Params file, submission script and helper scripts

**Params file.** `{PARAMS_YAML}` = `{SHEET_PREFIX}_params.yaml` in `{CWD}`; if it exists, ask (numbered): 1. overwrite · 2. choose another filename. Build it from the template below: replace the `{MODULE_PARAMS}` line by the Step 8 block, fill every placeholder, delete the `star_index` line when `{STAR_INDEX}` is empty and the `salmon_index` line when `{SALMON_INDEX}` is empty. No other line is deleted, and no placeholder may remain.

**Params file template.**
```yaml
# nf-core/rnasplice {VERSION} parameters — generated by /nfcore-rnasplice-setup
input: "{SAMPLESHEET_CSV}"
contrasts: "{CONTRASTS_CSV}"
source: "{SOURCE}"
outdir: "{OUTDIR}"
multiqc_title: "{MULTIQC_TITLE}"
fasta: "{FASTA_PATH}"
gtf: "{GTF_PATH}"
gencode: false
star_index: "{STAR_INDEX}"
salmon_index: "{SALMON_INDEX}"
save_reference: false
max_cpus: 16
max_memory: "64.GB"
max_time: "24.h"
{MODULE_PARAMS}
```

Typing rules: paths and strings are double-quoted; numbers and booleans are bare (`rmats_read_len: 101`, `rmats: true`); `multiqc_title` is always double-quoted, so that a numeric-looking title (for example `261001`) stays a string. Why a params file: under Nextflow 26.04, values given on the command line (a double-dash option followed by a value) reach the parameter validator as strings, and numeric parameters were rejected ("Value is [string] but should be [number]"); a params file keeps the YAML types. So every pipeline parameter is in this file, and the launch line carries only `-params-file`.

**Submission script.** Write `nf-core_rnasplice_{VERSION_TAG}.sh` in `{CWD}` (if it exists, ask: 1. overwrite · 2. choose another filename):
```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=8G
#SBATCH -t 48:00:00
#SBATCH -p bcc
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={USER_EMAIL}
#SBATCH -o nf-core_rnasplice_{VERSION_TAG}.%j.log

module add miniconda3/v4 || { echo "ERROR: cannot load module miniconda3/v4" >&2; exit 1; }
source /home/software/conda/miniconda3/bin/condainit || { echo "ERROR: cannot source condainit" >&2; exit 1; }
conda activate {CONDA_ENV} || { echo "ERROR: cannot activate conda environment {CONDA_ENV}" >&2; exit 1; }
module add singularity/3.10.4 || { echo "ERROR: cannot load module singularity/3.10.4" >&2; exit 1; }
command -v singularity >/dev/null || { echo "ERROR: singularity is not on PATH" >&2; exit 1; }
command -v nextflow >/dev/null || { echo "ERROR: nextflow is not in conda environment {CONDA_ENV}" >&2; exit 1; }

# Nextflow version range verified for nf-core/rnasplice «VERSION_TAG» by this skill
NF_MIN="«NEXTFLOW_MIN»"; NF_MAX_EXCL="«NEXTFLOW_MAX_EXCL»"
NF_VER=$(nextflow -version 2>/dev/null | awk '/version/ {for (i = 1; i <= NF; i++) if ($i ~ /^[0-9]+\.[0-9]+\.[0-9]+/) {print $i; exit}}')
ver_ge() { [ "$(printf '%s\n%s\n' "$1" "$2" | sort -V | head -n1)" = "$2" ]; }
[ -n "$NF_VER" ] || { echo "ERROR: could not read the Nextflow version" >&2; exit 1; }
ver_ge "$NF_VER" "$NF_MIN" || { echo "ERROR: Nextflow $NF_VER is older than $NF_MIN, the minimum for this pipeline revision" >&2; exit 1; }
if [ "$NF_MAX_EXCL" != none ] && ver_ge "$NF_VER" "$NF_MAX_EXCL"; then
  echo "ERROR: Nextflow $NF_VER is $NF_MAX_EXCL or newer; this pipeline revision was verified only with older versions (see the skill README)" >&2; exit 1
fi
echo "Nextflow $NF_VER (verified with «NEXTFLOW_TESTED»)"

# ~/.bashrc may export NXF_OFFLINE=TRUE, which blocks the first download of the pipeline: pull with it switched off
NXF_OFFLINE=false nextflow pull nf-core/rnasplice -r {VERSION} \
  || echo "WARNING: nextflow pull failed; using the cached copy of nf-core/rnasplice {VERSION} if there is one" >&2

nextflow run nf-core/rnasplice -r {VERSION} -c nextflow.config -profile slurm,singularity -params-file {PARAMS_YAML}
```
The head job only coordinates the pipeline (2 CPUs, 8 GB) but must outlive every task, hence 48 h.

**Genome download helper.** Only when Step 7 found a missing Ensembl FASTA or GTF: write `download_genome_{REF_TAG}.sh` (URLs verified in this session, Step 7). It keeps the `.gz` files next to the decompressed ones, can be re-run safely, and removes only its own partial files:
```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH -t 4:00:00
#SBATCH -p bcc
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={USER_EMAIL}
#SBATCH -o download_genome_{REF_TAG}.%j.log
GENOME_DIR="{GENOME_DIR}"
case "$GENOME_DIR" in /*) ;; *) echo "ERROR: GENOME_DIR must be an absolute path: '$GENOME_DIR'" >&2; exit 1 ;; esac
mkdir -p "$GENOME_DIR" || { echo "ERROR: cannot create $GENOME_DIR" >&2; exit 1; }
fetch() {  # fetch <url of a .gz file> <decompressed target path>
  local url=$1 out=$2
  if [ -s "$out" ]; then echo "present: $out"; return 0; fi
  if ! { [ -s "$out.gz" ] && gzip -t "$out.gz" 2>/dev/null; }; then
    wget -c -O "$out.gz.part" "$url" || { echo "ERROR: download failed for $url (partial file kept for resuming: $out.gz.part)" >&2; exit 1; }
    gzip -t "$out.gz.part" 2>/dev/null || { rm -f "$out.gz.part"; echo "ERROR: $url is not a valid gzip file; partial file removed, run the helper again" >&2; exit 1; }
    mv "$out.gz.part" "$out.gz" || { echo "ERROR: cannot move $out.gz.part to $out.gz" >&2; exit 1; }
  fi
  gunzip -c "$out.gz" > "$out.part" && [ -s "$out.part" ] && mv "$out.part" "$out" \
    || { rm -f "$out.part"; echo "ERROR: decompression failed for $out.gz" >&2; exit 1; }
  echo "ready: $out"
}
fetch "{ENSEMBL_FASTA_URL}" "{FASTA_PATH}"
fetch "{ENSEMBL_GTF_URL}" "{GTF_PATH}"
echo "genome files ready in $GENOME_DIR"
```

Show every written file in full and tell the user:
```
Params file written: {PARAMS_YAML}
Script written: nf-core_rnasplice_{VERSION_TAG}.sh
To submit:  sbatch nf-core_rnasplice_{VERSION_TAG}.sh
```
If a helper was written, the pipeline must wait for it: `jid=$(sbatch --parsable download_genome_{REF_TAG}.sh)` and then `sbatch --dependency=afterok:$jid nf-core_rnasplice_{VERSION_TAG}.sh`, listing every helper job in the dependency.
````

[C only] Also append, after the genome download helper:
````markdown
**Nextflow environment helper.** Only when the user has no conda environment with Nextflow «NEXTFLOW_TESTED» (Step 2): write `create_nextflow_env_«NEXTFLOW_TESTED».sh`; it runs conda on a compute node and does nothing if the environment exists:
```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH -t 1:00:00
#SBATCH -p bcc
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={USER_EMAIL}
#SBATCH -o create_nextflow_env_«NEXTFLOW_TESTED».%j.log
module add miniconda3/v4 || { echo "ERROR: cannot load module miniconda3/v4" >&2; exit 1; }
source /home/software/conda/miniconda3/bin/condainit || { echo "ERROR: cannot source condainit" >&2; exit 1; }
if conda env list | awk '{print $1}' | grep -qx "nf-«NEXTFLOW_TESTED»"; then echo "environment nf-«NEXTFLOW_TESTED» already exists"; exit 0; fi
conda create -y -n nf-«NEXTFLOW_TESTED» -c conda-forge -c bioconda nextflow=«NEXTFLOW_TESTED» || { echo "ERROR: conda create failed" >&2; exit 1; }
echo "environment nf-«NEXTFLOW_TESTED» ready"
```
Submit it first and add its job id to the pipeline's `--dependency=afterok:` list.
````

- [ ] **Step 4: Run the checker and the three tests; expected PASS**

```bash
for t in check_skill.sh test_render_params.sh dry_run_submit.sh dry_run_helpers.sh; do bash nfcore-rnasplice-setup/tests/$t nfcore-rnasplice-setup/nfcore-rnasplice-setup.md; done
```
Expected: `PASS`, `RENDER PASS`, `DRY RUN SUBMIT PASS`, `DRY RUN HELPERS PASS`.

- [ ] **Step 5: Mutation rows, RED proof, commit**

Append (the `T6-sel` row only if `RMATS_POST` is a component in `trace_process_names.txt`; `T6-selmiss` only if `STAR_ALIGN` is):
```
T6-flag	check	skill	s/-params-file {PARAMS_YAML}$/-params-file {PARAMS_YAML} --rmats_read_len 101/	launch line must be exactly
T6-cont	check	skill	s/^nextflow run nf-core\/rnasplice -r {VERSION} -c nextflow.config -profile slurm,singularity -params-file {PARAMS_YAML}$/nextflow run nf-core\/rnasplice -r {VERSION} -c nextflow.config -profile slurm,singularity \\\n  --input x.csv -params-file {PARAMS_YAML}/	launch line must be exactly
T6-lmod	check	skill	s/^module add singularity\/3.10.4 || .*$/module add singularity\/3.10.4/	unguarded module load
T6-confine	check	skill	$s/$/\nnextflow -version/	nextflow is invoked outside the submission script
T6-sel	check	skill	s/withName: '\.\*:RMATS_POST' {/withName: '.*:RMATS_POSTX' {/	selector RMATS_POSTX matches no process
T6-selmiss	check	skill	/withName: '\.\*:STAR_ALIGN' {/,+4d	process STAR_ALIGN ran in the verification runs but has no withName selector
T6-over	check	skill	s/^trace    { enabled = true; overwrite = true;/trace    { enabled = true;/	overwrite = true must be set
T6-nfmin	submit	skill	s/^NF_MIN="[^"]*"/NF_MIN="0.0.0"/	FAIL: case below_min
T6-guard	submit	skill	s/^conda activate {CONDA_ENV} || .*$/conda activate {CONDA_ENV}/	FAIL: case conda_fails
T6-pull	submit	skill	s/^NXF_OFFLINE=false nextflow pull/nextflow pull/	FAIL: case pull_offline_override
T6-render	render	skill	s/^{MODULE_PARAMS}$/{MODULE_PARAMS}\noutdir: "x"/	FAIL: case full
T6-gzipt	helper	skill	/gzip -t "\$out.gz.part" 2>\/dev\/null ||/d	FAIL: case truncated_download
T6-abs	helper	skill	/case "\$GENOME_DIR" in/d	FAIL: case relative_dir
```
```bash
bash nfcore-rnasplice-setup/tests/prove_mutations.sh nfcore-rnasplice-setup/nfcore-rnasplice-setup.md
bash nfcore-rnasplice-setup/tests/prove_red.sh <BASE>
git status && git diff --stat && git diff
git add nfcore-rnasplice-setup
git commit -m "nfcore-rnasplice-setup: Steps 10-11 config with verified selectors, params file, submission script, helpers; render and dry-run tests"
```

---

### Task 7: Final summary, hand-off note, "rnasplice or bulk-rnaseq-pipeline", Notes (Step 12)

**Files:**
- Modify: `D/tests/check_skill.sh`, `D/tests/mutations.tsv`, `D/nfcore-rnasplice-setup.md` (append Step 12 and Notes)

**Interfaces:**
- Consumes: `{OUTDIR}`, `{VERSION}`, the module choices (Task 5); fixture `output_tree.txt`; gate values `TEST_CONTRAST`, `RMATS_B1_GROUP`, `SALMON_ROUTE`, `PIPELINE_REVISION`.
- Produces: anchor `**Hand-off note.**` whose block lists output paths relative to `{OUTDIR}` in the form `  <path>  <description>` (two leading spaces, then at least two spaces before the description; `{CONTRAST}` is the only placeholder allowed in a path). Task 8 copies these paths into the README; Task 9 checks they exist in a real run.

- [ ] **Step 1: Append the Task 7 checker block; see it fail**

```bash
# --- Task 7 (Step 12, Notes)
need "## Step 12 — Summary and hand-off"
anchor_once "**Hand-off note.**"
HO=$(bash "$CUT" "$SKILL" "**Hand-off note.**" 2>/dev/null)
n=0
while IFS= read -r p; do
  n=$((n + 1)); q=${p//\{CONTRAST\}/$(gv TEST_CONTRAST)}; q=${q%/}
  grep -qxF -- "./$q" "$FIX/output_tree.txt" || { echo "FAIL: hand-off path $p not in the recorded output tree"; fail=1; }
done < <(printf '%s\n' "$HO" | sed -n 's/^  \([^ ]\{1,\}\)  .*$/\1/p')
[ "$n" -ge 6 ] || { echo "FAIL: the hand-off note lists $n output paths (at least 6 expected)"; fail=1; }
! printf '%s\n' "$HO" | grep -q '<' || { echo "FAIL: the hand-off note still contains a <...> placeholder"; fail=1; }
need ".MATS.JC.txt"
need ".MATS.JCEC.txt"
case "$(gv RMATS_B1_GROUP)" in
  treatment) need "the treatment samples are rMATS's \`b1\`" ;;
  control) need "this pipeline passes the control samples as rMATS's \`b1\`" ;;
esac
need "**rnasplice or \`/bulk-rnaseq-pipeline\`?**"
need "do not run both on the same samples"
need "This skill does not analyse the results"
need "## Notes for the assistant"
need "Raw FASTQ/BAM files are read-only"
need "Never pass pipeline parameters on the \`nextflow run\` line"
need "Never run Nextflow, conda, Java, Python or R on the login node"
[ "$(gv SALMON_ROUTE)" = star_salmon_both ] && need "use the \`star_salmon\` ones"
# --- end Task 7
```

- [ ] **Step 2: Fill the hand-off lines from the fixture**

Each line of the hand-off block is a real path of `D/tests/fixtures/output_tree.txt` (written relative to `{OUTDIR}`, without the leading `./`, the test contrast name replaced by `{CONTRAST}`, directories ending in `/`). Find them with:
```bash
F=nfcore-rnasplice-setup/tests/fixtures/output_tree.txt; C=$(awk -F'\t' '$1 == "TEST_CONTRAST" {print $2}' nfcore-rnasplice-setup/tests/fixtures/gate_values.tsv)
grep -E '/SE\.MATS\.JC\.txt$' $F | head -3                     # rMATS: its directory is the per-contrast rMATS result directory
grep -iE 'dexseq' $F | grep -viE 'dtu|annotation' | head -10   # DEXSeq exon usage results
grep -iE 'edger' $F | head -10                                 # edgeR exon usage results
grep -iE 'dtu|stager|drimseq' $F | head -10                    # DTU results
grep -iE 'suppa' $F | head -10                                 # SUPPA2 results
grep -iE 'multiqc.*\.html$' $F                                 # MultiQC report
echo "$C"
```
Lines, in this order (one per row; leave a row out only if the fixture has no such output, and say so in the report): rMATS per-contrast directory; its `SE.MATS.JC.txt`; DEXSeq exon usage directory; edgeR exon usage directory; DTU directory; SUPPA2 directory; MultiQC report file; `pipeline_info/` is not in the tree (the tree excludes it), so it is written in the prose sentence after the block, not as a path line.

- [ ] **Step 3: Append Step 12 and the Notes to the skill**

The block shows the line format; the `<...>` parts are replaced by the fixture paths of Step 2 (the skill contains no `<` placeholders). Use the variant for `RMATS_B1_GROUP` and the [star_salmon_both] sentence only for that route.

````markdown

---

## Step 12 — Summary and hand-off

Show the user what was written (samplesheet, contrasts sheet, params file, `nextflow.config` if it was written, the submission script, any helper) and how to submit (Step 11 order). Then print the hand-off note, keeping only the lines of the analyses that run. If `{VERSION}` is not «PIPELINE_REVISION», first confirm the directory names against the pipeline's `docs/output.md` for `{VERSION}` with `WebFetch`; for «PIPELINE_REVISION» they come from this skill's verification run.

**Hand-off note.**
```text
Outputs under {OUTDIR}/ :
  <rMATS per-contrast directory>/  rMATS, one directory per contrast: for each event type (SE, A5SS, A3SS, MXE, RI) a .MATS.JC.txt and a .MATS.JCEC.txt table
  <rMATS per-contrast directory>/SE.MATS.JC.txt  example: skipped exons, junction reads only
  <DEXSeq exon usage directory>/  DEXSeq differential exon usage per contrast
  <edgeR exon usage directory>/  edgeR differential exon usage per contrast
  <DTU directory>/  DEXSeq differential transcript usage with stageR
  <SUPPA2 directory>/  SUPPA2 event and isoform PSI and differential splicing
  <MultiQC report file>  MultiQC report
```
`pipeline_info/` holds the execution report, the timeline and the trace (requested and used resources per task).

**Reading the rMATS tables.** Each event type (SE, A5SS, A3SS, MXE, RI) has a `.MATS.JC.txt` table (junction-spanning reads only) and a `.MATS.JCEC.txt` table (junction and exon-body reads).
[treatment] `IncLevelDifference` = inclusion level of the treatment group minus that of the control group (in this pipeline the treatment samples are rMATS's `b1`); positive values mean more inclusion in the treatment.
[control] `IncLevelDifference` = inclusion level of the control group minus that of the treatment group (this pipeline passes the control samples as rMATS's `b1`); positive values mean more inclusion in the control, so flip the sign to read it as treatment minus control.
Usual filters: FDR < 0.05 and |IncLevelDifference| > 0.1, with read support in every replicate (no missing `IJC_SAMPLE_*`/`SJC_SAMPLE_*` values).

**Which output answers which question.**

| Question | Output |
|---|---|
| Which splice events (cassette exons, alternative 5'/3' splice sites, mutually exclusive exons, retained introns) change between conditions? | rMATS (from junction reads); SUPPA2 local events (from transcript abundance) |
| Which parts (exon bins) of a gene are used more or less? | DEXSeq and edgeR exon usage |
| Which transcripts of a gene change their share of the gene's expression? | DEXSeq DTU (stageR-confirmed); SUPPA2 isoform PSI |

[star_salmon_both] DTU and SUPPA2 results appear twice (under the STAR+Salmon and the Salmon directories); use the `star_salmon` ones.

**rnasplice or `/bulk-rnaseq-pipeline`?**
- Gene-level differential expression, GSEA and transcript usage on the Salmon output of an nf-core/rnaseq run: `/bulk-rnaseq-pipeline`. Its optional DTU module is the same workflow (DRIMSeq filter → DEXSeq → stageR) as the DEXSeq DTU of rnasplice; do not run both on the same samples.
- Which exons or splice events change (cassette exons, alternative 5'/3' splice sites, retained introns, mutually exclusive exons), for example after knocking out a splicing factor: this pipeline with rMATS (optionally SUPPA2); usage of exon bins: DEXSeq or edgeR exon usage.
- This skill does not analyse the results; a report for the rMATS tables is a separate, later skill.

---

## Notes for the assistant

- **`gh` CLI is not available on this HPC cluster.** Use `WebFetch` for all GitHub API calls.
- **Always present finite-choice questions as numbered lists.** Use open questions only for the email, the conda environment, paths and the plain-language description of the groups.
- Never run Nextflow, conda, Java, Python or R on the login node. The wizard runs only light commands there (`find`, `ls`, `test`, `zcat | head`, `grep`, `awk`, `mktemp`, `mv`); all pipeline work goes through `sbatch`.
- Raw FASTQ/BAM files are read-only; never modify them.
- Never pass pipeline parameters on the `nextflow run` line; every parameter is in the params file.
- Every module switch is written explicitly; `sashimi_plot` is always `false`.
- Strandedness is asked, never guessed; one value for all samples.
- `rmats_read_len` is always written (detected for FASTQ, asked for BAM).
- `rmats_paired_stats` and `diffsplice_paired` are `true` only for a confirmed paired design.
- Both sheets are validated with `validate_rnasplice_sheets` before they are written; never edit the validator to make a sheet pass.
- Never overwrite an existing `nextflow.config`; ask before overwriting any other file.
- Never embed a download URL that was not verified in this session.
- Every `module add` in a generated script is guarded (`|| { ...; exit 1; }`): under Lmod a failed `module add` silently breaks every later one.
````

- [ ] **Step 4: Run the checker; expected PASS**

- [ ] **Step 5: Mutation rows, RED proof, commit**

```
T7-path	check	skill	s/^Outputs under {OUTDIR}\/ :$/Outputs under {OUTDIR}\/ :\n  no_such_dir\/  bogus/	hand-off path no_such_dir/ not in the recorded output tree
T7-few	check	skill	/^  [^ ]*SE\.MATS\.JC\.txt  /d	output paths (at least 6 expected)
```
(`T7-few` applies only when exactly 6 path lines are written; with 7 lines its program is `/^  [^ ]*SE\.MATS\.JC\.txt  /d;/MultiQC report$/d`.)
```bash
bash nfcore-rnasplice-setup/tests/prove_mutations.sh nfcore-rnasplice-setup/nfcore-rnasplice-setup.md
bash nfcore-rnasplice-setup/tests/prove_red.sh <BASE>
git status && git diff --stat && git diff
git add nfcore-rnasplice-setup
git commit -m "nfcore-rnasplice-setup: Step 12 hand-off from the verified output tree, rMATS sign, bulk-rnaseq-pipeline guidance, notes"
```

---

### Task 8: Checker consolidation, README, root README

**Files:**
- Create: `D/README.md`, `D/tests/run_all_tests.sh`
- Modify: `D/tests/check_skill.sh` (README block; tidy), `D/tests/mutations.tsv`, `README.md` (root: row and intro bullet)

**Interfaces:**
- Consumes: every earlier task; gate values (`GATE_DATE`, `RUN_DIRS`, `NEXTFLOW_TESTED`, `PIPELINE_REVISION`, `BAM_*`, `SALMON_ROUTE`, `EDGER_DEU_FUNCTION`).
- Produces: `tests/run_all_tests.sh` (prints `ALL TESTS PASS`, used by Tasks 9-10); the README line `Cluster acceptance (nf-core test data): PENDING` that Task 9 replaces by `Cluster acceptance (nf-core test data): DONE (YYYY-MM-DD)`.

- [ ] **Step 1: Create run_all_tests.sh**

```bash
#!/bin/bash
# Usage: run_all_tests.sh   — runs every static test of the skill; prints "ALL TESTS PASS" or exits 1
set -u
HERE=$(cd "$(dirname "$0")" && pwd); SK="$HERE/../nfcore-rnasplice-setup.md"; rc=0
for t in check_skill.sh test_strandedness.sh test_bam_policy.sh test_validate_sheets.sh test_read_length.sh test_render_params.sh dry_run_submit.sh dry_run_helpers.sh prove_mutations.sh; do
  echo "== $t"; bash "$HERE/$t" "$SK" || rc=1
done
grep -qF '| nf-core/rnasplice setup | `/nfcore-rnasplice-setup` |' "$HERE/../../README.md" || { echo "FAIL: root README row missing"; rc=1; }
! grep -qF '(in development, not yet validated)' "$HERE/../../README.md" || { echo "FAIL: root README row still marked in development"; rc=1; }
[ $rc -eq 0 ] && echo "ALL TESTS PASS" || exit 1
```

- [ ] **Step 2: Append the README checker block; see it fail**

```bash
# --- Task 8 (README)
needr "# \`/nfcore-rnasplice-setup\` — nf-core/rnasplice Differential Splicing Setup Skill"
needr "cp nfcore-rnasplice-setup.md ~/.claude/commands/"
needr "## Validation status"
needr "nf-core/rnasplice $(gv PIPELINE_REVISION)"
needr "Nextflow $(gv NEXTFLOW_TESTED)"
needr "4 paired-end human chrX samples"
needr "No real biological data"
needr "## Known limitations"
needr "## rnasplice or \`/bulk-rnaseq-pipeline\`"
needr "sashimi_plot"
needr "tests/run_all_tests.sh"
needr "## Where the results are"
# no plan variant marker may survive in the skill or the README
VM='\[(A/C|A/B|B|C|BAM|NO-BAM|fixed|from_sheet|none|pseudo_only|star_salmon_only|star_salmon_both|treatment|control|sheet|sorted_by_name|unordered|all_samples|per_contrast)\]|\[C: '
! grep -nE "$VM" "$SKILL" || { echo "FAIL: a plan variant marker is left in the skill (lines above)"; fail=1; }
! grep -nE "$VM" "$README" 2>/dev/null || { echo "FAIL: a plan variant marker is left in the README (lines above)"; fail=1; }
grep -qE '^Cluster acceptance \(nf-core test data\): (PENDING|DONE \([0-9]{4}-[0-9]{2}-[0-9]{2}\))' "$README" 2>/dev/null \
  || { echo "FAIL: README must have the line 'Cluster acceptance (nf-core test data): PENDING' or '... DONE (YYYY-MM-DD)'"; fail=1; }
while IFS= read -r p; do
  needr "$p"
done < <(printf '%s\n' "$HO" | sed -n 's/^  \([^ ]\{1,\}\)  .*$/\1/p' | head -3)
# --- end Task 8
```
(This block goes after the Task 7 block, so `$HO` is set.)

- [ ] **Step 3: Write `D/README.md`**

Substitute `«KEY»`; use the BAM sentence that matches `BAM_RMATS_LIBTYPE` (fixed value / `from_sheet` / `none`):

````markdown
# `/nfcore-rnasplice-setup` — nf-core/rnasplice Differential Splicing Setup Skill

An interactive Claude Code skill that sets up and submits [nf-core/rnasplice](https://nf-co.re/rnasplice) for differential alternative splicing of bulk RNA-seq on an HPC cluster with SLURM and Singularity: rMATS (splicing events), SUPPA2 (event and isoform PSI), DEXSeq and edgeR differential exon usage, and DEXSeq differential transcript usage (DRIMSeq, DEXSeq, stageR). It follows the style of `/nfcore-rnavar-setup`: every pipeline parameter goes into one params file, and every analysis switch is written explicitly. It sets up the run; it does not analyse the results.

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
| Nextflow «NEXTFLOW_MIN» or newer [C: and older than «NEXTFLOW_MAX_EXCL»] | In a conda environment; verified with Nextflow «NEXTFLOW_TESTED» (env `«CONDA_ENV_TESTED»`). The submission script checks the version on the compute node |
| Internet on compute nodes | The job pulls the pipeline and the containers (depot.galaxyproject.org). If `~/.bashrc` exports `NXF_OFFLINE=TRUE`, the script still pulls the pipeline (`NXF_OFFLINE=false nextflow pull`) |
| Internet on the login node | The wizard uses `WebFetch` (GitHub API, Ensembl) |

---

## What the skill does

| Step | Topic | Asked or detected |
|------|-------|-------------------|
| 0 | Working directory | Detected (`pwd`) |
| 1 | Email for SLURM notifications | Asked |
| 2 | Conda environment with Nextflow | Asked (version checked later, on the compute node) |
| 3 | Pipeline version | Detected (GitHub API); default nf-core/rnasplice «PIPELINE_REVISION», the verified revision |
| 4 | Input and samplesheet rows | FASTQ or genome BAM; files found, paired-end detected, names sanitised; strandedness asked |
| 5 | Conditions, contrasts, paired design | Asked; both sheets validated before they are written |
| 6 | Read length for rMATS | Detected from the FASTQ (asked for BAM) |
| 7 | Organism and genome | Shared genome folder or custom FASTA/GTF; compatible STAR/Salmon indexes reused |
| 8 | Analyses and their settings | Asked: rMATS (default), SUPPA2, DEXSeq/edgeR exon usage, DEXSeq DTU |
| 9 | Trimming and QC | Pipeline defaults |
| 10 | MultiQC title, output directory, `nextflow.config` | Asked; config written only if absent |
| 11 | Params file, submission script, helper | Generated |
| 12 | Summary and hand-off | Where each result is and which question it answers |

### Key design points

- **Every analysis switch is explicit.** In the pipeline's own configuration every module (rMATS, DEXSeq, edgeR, DTU, SUPPA2, MISO sashimi plots) is on, so a switch left out of the params file runs that module anyway. The params file always contains all of them.
- **MISO is not used** (`sashimi_plot: false`): MISO is archived (Python 2) and rnasplice uses it only for sashimi plots of a fixed gene list.
- **Strandedness is asked**, never guessed (rnasplice has no `auto`), one value for all samples (rMATS needs it); the wizard can read the RSeQC results of an nf-core/rnaseq run of the same samples.
- **`rmats_read_len` is always set** from the reads (the pipeline default of 40 is wrong for almost all data).
- **Paired rMATS statistics only when confirmed:** `rmats_paired_stats` and SUPPA2's `diffsplice_paired` default to `true` in the pipeline; the skill writes `false` unless the user confirms a paired design.
- **DTU filters are explicit** (the pipeline's schema and configuration disagree); the skill uses the Love et al. 2018 rule.
- **Genome BAM input:** [fixed] accepted only for «BAM_RMATS_LIBTYPE» / «BAM_RMATS_READTYPE»-end libraries, because this revision cannot be told the strandedness or read type of BAM input; other libraries are pointed to FASTQ input. [from_sheet] strandedness and read type are written to the BAM samplesheet. [none] not available: genome-BAM input did not work with this revision in the verification run.
- **Salmon route:** «SALMON_ROUTE» (see the skill's Step 8).

---

## Output files

| File | Description |
|------|-------------|
| `{prefix}_samplesheet.csv` | Samplesheet (`sample,fastq_1,fastq_2,strandedness,condition`, or the genome-BAM form) |
| `{prefix}_contrasts.csv` | Contrasts (`contrast,treatment,control`) |
| `{prefix}_params.yaml` | Every pipeline parameter of the run, passed with `-params-file` |
| `nextflow.config` | SLURM + Singularity profiles with process selectors (written only if absent) |
| `nf-core_rnasplice_{version_tag}.sh` | Pipeline submission script (checks the Nextflow version, pulls the pipeline, runs it) |
| `download_genome_{ref_tag}.sh` | Ensembl FASTA/GTF download (only if missing) |
[C: | `create_nextflow_env_«NEXTFLOW_TESTED».sh` | Creates a conda environment with Nextflow «NEXTFLOW_TESTED» (only if needed) |]

## Where the results are

Relative to the output directory, as recorded by the verification run (`{CONTRAST}` = the contrast name); copy the path lines of the skill's Step 12 hand-off note into this table, one row each, with the same descriptions:

| Path | Content |
|------|---------|
| (each hand-off path, verbatim) | (its description) |

`pipeline_info/` holds the execution report, timeline and trace.

---

## Typical workflow

```bash
# 1. Only if a helper was generated
jid=$(sbatch --parsable download_genome_GRCm39_ens115.sh)
# 2. The pipeline, after the helper
sbatch --dependency=afterok:$jid nf-core_rnasplice_«VERSION_TAG».sh
```
(The helper name is an example.)

---

## Validation status

- Static tests: `tests/run_all_tests.sh` runs the checker (every parameter against the recorded schema of nf-core/rnasplice «PIPELINE_REVISION», module switches, defaults against the recorded pipeline configuration, process selectors against the processes of the verification runs, the launch line, guarded module loads), the tests that run the skill's own code blocks (sheet validation, strandedness from RSeQC, BAM rule, read length, params rendering) and the stub dry runs of the submission script and the download helper, and the mutation proofs of every structural check.
- Verification gate («GATE_DATE»): nf-core/rnasplice «PIPELINE_REVISION» with Nextflow «NEXTFLOW_TESTED» on the cluster (runs: «RUN_DIRS»): the pipeline's own test profile, a run with a fully typed params file and all modules, an rMATS-only run, and a genome-BAM run. The test data are 4 paired-end human chrX samples (about 30 MB of FASTQ). Details: `tests/fixtures/gate_report.md`.
Cluster acceptance (nf-core test data): PENDING
- No real biological data was run. Resources in `nextflow.config` are judgement, not measured on a real genome.

---

## Known limitations

- nf-core/rnasplice «PIPELINE_REVISION» is the verified revision; a newer release is offered only after its schema has been checked, and is not verified by this skill.
- Genome-BAM input: see Key design points.
- rMATS needs one strandedness and one read type for all samples.
- Paired designs: available only when this revision keeps a fixed sample order for rMATS (verification result: «RMATS_BAMLIST_ORDER»).
- Trimming and SUPPA2 clustering use the pipeline defaults.
- edgeR exon usage uses «EDGER_DEU_FUNCTION» in this revision.
- MISO, IsoformSwitchAnalyzeR and LeafCutter are not offered.
- The skill does not analyse the results.

---

## rnasplice or `/bulk-rnaseq-pipeline`

`/bulk-rnaseq-pipeline` does gene-level differential expression, GSEA and, optionally, differential transcript usage (DRIMSeq → DEXSeq → stageR) on the Salmon output of an nf-core/rnaseq run. rnasplice answers which exons and splice events change (rMATS, SUPPA2, DEXSeq/edgeR exon usage). Its DEXSeq DTU is the same workflow as the bulk skill's DTU module: run it in one place only.
````

- [ ] **Step 4: Root README**

In `README.md` (root): replace the stub row by
```
| nf-core/rnasplice setup | `/nfcore-rnasplice-setup` | Interactive setup wizard for nf-core/rnasplice differential alternative splicing (rMATS, SUPPA2, DEXSeq/edgeR exon usage, DEXSeq DTU; every module switch explicit) on SLURM + Singularity |
```
and change the "Pipeline setup" bullet to end with `nf-core/rnavar (RNA-seq variant calling) and nf-core/rnasplice (differential alternative splicing).`

- [ ] **Step 5: Mutation rows, consolidated run, RED proof, commit**

```
T8-marker	check	readme	s/^Cluster acceptance (nf-core test data): /Cluster acceptance: /	README must have the line
T8-data	check	readme	s/No real biological data/Real data/	README missing required text: No real biological data
```
Tidy `check_skill.sh` only by removing exact duplicate `need` lines (keep the first occurrence); no rule is weakened or removed. Then:
```bash
bash nfcore-rnasplice-setup/tests/run_all_tests.sh
bash nfcore-rnasplice-setup/tests/prove_red.sh <BASE>
git status && git diff --stat && git diff
git add nfcore-rnasplice-setup README.md
git commit -m "nfcore-rnasplice-setup: README, run_all_tests, root README registration"
```
Expected: `ALL TESTS PASS`, `RED PASS (...)`.

---

### Task 9: Cluster acceptance on the nf-core test data (controller-run)

The controller follows the skill text of the branch by hand (it is not installed into `~/.claude/commands/`), with the answers below, in a fresh directory, and submits the generated files. Every Nextflow job runs through `sbatch`; the login node only writes files and reads results.

**Files:**
- Create (scratch): `/net/bmc-lab3/data/bcc/yannvrb/rnasplice_accept/` (`data/` copied from `/net/bmc-lab3/data/bcc/yannvrb/rnasplice_gate/data/`)
- Modify: `D/README.md` (acceptance record), possibly skill and tests (defects found here go through a fix dispatch with tests, not controller edits)

**Interfaces:**
- Consumes: the whole skill; `tests/render_params.sh`; gate fixtures.
- Produces: README line `Cluster acceptance (nf-core test data): DONE (YYYY-MM-DD)` with the run record.

- [ ] **Step 1: Static gate**

`bash nfcore-rnasplice-setup/tests/run_all_tests.sh` → `ALL TESTS PASS`.

- [ ] **Step 2: Run (a), all modules, by following the wizard**

In `/net/bmc-lab3/data/bcc/yannvrb/rnasplice_accept/a/` follow Steps 0-12 with: email `yannvrb@mit.edu`; conda env `«CONDA_ENV_TESTED»`; FASTQ input from `../data/` (copies of the gate data); strandedness 1 (unstranded; the test data are unstranded); auto names; conditions from the upstream sheet (GBR, YRI), reference YRI; all pairwise (one contrast `GBR_vs_YRI`); not paired; custom reference (`../data/X.fa.gz`, `../data/genes_chrX.gtf`), no existing indexes; analyses 1, 2, 3, 4, 5; defaults elsewhere. Check before submitting:
```bash
grep -v '^#SBATCH' nf-core_rnasplice_*.sh | grep -c ' --'     # expected 0: no double-dash option outside the #SBATCH header
grep -E '^(rmats|dexseq_exon|edger_exon|dexseq_dtu|suppa|sashimi_plot|rmats_read_len|rmats_paired_stats):' *_params.yaml
```
Expected: all five switches `true`, `sashimi_plot: false`, `rmats_read_len: «TEST_READ_LENGTH»`, `rmats_paired_stats: false`. Render the same answers with `tests/render_params.sh` (values file written from the answers) and `diff` it with the wizard's params file: identical. Submit `sbatch nf-core_rnasplice_«VERSION_TAG».sh`; wait with a Monitor until-loop on `squeue -j <id>`.

Expected after the run: the log has `Nextflow «NEXTFLOW_TESTED» (verified with «NEXTFLOW_TESTED»)` and `Pipeline completed successfully`; in `pipeline_info/execution_trace.txt` every task `COMPLETED`; tasks of RMATS_PREP/RMATS_POST, DEXSEQ_EXON, EDGER_EXON, DEXSEQ_DTU, STAGER and SUPPA present; no MISO task. Selector check:
```bash
TR=$(ls -d results/*/pipeline_info)/execution_trace.txt
for s in $(sed -n "s/^[ \t]*withName: '\.\*:\([A-Z0-9_]*\)' {$/\1/p" nextflow.config); do
  awk -F'\t' -v s="$s" 'NR > 1 {n = $4; sub(/ \(.*$/, "", n); m = n; sub(/.*:/, "", m); if (m == s) {print s, $7, $8, $9; f = 1}} END {if (!f) print s, "NO TASK"}' "$TR" | sort -u
done
```
Expected: every selector has tasks, with the cpus/memory/time of its row (memory as shown by Nextflow, e.g. `48 GB`). And the rMATS read length: the RMATS_PREP `.command.sh` (work dir from the trace hash) contains `--readLength «TEST_READ_LENGTH»`. Every hand-off path of Step 12 exists under the output directory (with `{CONTRAST}` = `GBR_vs_YRI`).

- [ ] **Step 3: Run (b), rMATS only**

In `../b/`, the same answers but analyses `1` only (default). Expected: `rmats: true` and the four other switches `false` in the params file; the run completes; `awk -F'\t' 'NR > 1 {print $4}' <trace> | grep -cE 'DEXSEQ|EDGER|DRIMSEQ|STAGER|SUPPA|MISO'` = 0; the rMATS directory has the JC/JCEC tables for all five event types.

- [ ] **Step 4: Run (c), genome-BAM input (only if `BAM_RMATS_LIBTYPE` is not `none`)**

In `../c/`: input 2 (genome BAM) with the BAMs published by run (a) (`find ../a/results -name "*.bam" ! -name "*toTranscriptome*"`), strandedness unstranded, paired-end. Expected: `bam_input_allowed unstranded paired` prints `ALLOWED` when the gate's values are `fr-unstranded`/`paired` (or `from_sheet`); analyses 1, 3, 4 offered (2 and 5 not); run completes; trace has RMATS, DEXSEQ_EXON and EDGER_EXON tasks. Negative check (no run): `bam_input_allowed reverse paired` prints `REFUSED: ... Start from FASTQ instead.` when the values are fixed.

- [ ] **Step 5: Record and commit**

Replace the README line `Cluster acceptance (nf-core test data): PENDING` by `Cluster acceptance (nf-core test data): DONE (YYYY-MM-DD)` followed by one bullet per run: job id, wall time, task count, verdict, what was checked (selectors, read length, switches, hand-off paths), and what was not exercised (real data, Ensembl download helper end to end, [C] the environment helper, paired designs). Run `bash nfcore-rnasplice-setup/tests/run_all_tests.sh` → `ALL TESTS PASS`.
```bash
git status && git diff --stat && git diff
git add nfcore-rnasplice-setup/README.md
git commit -m "nfcore-rnasplice-setup: record cluster acceptance on the nf-core test data"
```
A defect found here (a run fails, a selector does not apply, a hand-off path is missing) is recorded in the ledger and fixed by a dispatched fix task that adds a test or checker line proven RED first; the run is repeated.

---

### Task 10: Whole-branch final review and finishing

- [ ] **Step 1: Final review**

Dispatch a fresh reviewer on the most capable available model (earlier final reviews used `fable`) with: the spec, this plan, the ledger, `git diff <plan commit>..HEAD`, the gate report, and these focus points: silent wrong results from pipeline defaults (module switches, paired stats, read length, DTU filters, Salmon route), unvalidated user input reaching a sheet or script, unsafe `rm` paths, stale text (gate-dependent sentences vs `gate_values.tsv`, README vs skill), honest validation status, login-node work in the wizard text, Lmod guards. The reviewer runs `tests/run_all_tests.sh` and mutates at least five checks of its own choosing.

- [ ] **Step 2: One fix wave**

The controller rules each finding (`Ruling:` + `Cost if wrong:` in the ledger); one fresh implementer fixes all accepted Critical/Important findings and the cheap Minors in one wave, each fix with a test or checker line proven RED on the pre-fix commit (`prove_red.sh`) and mutation rows where structural. `ALL TESTS PASS` afterwards.

- [ ] **Step 3: Scoped re-review**

A fresh reviewer checks only the fix-wave diff against the findings. Remaining Minors are parked and listed for the user.

- [ ] **Step 4: Finishing**

Use superpowers:finishing-a-development-branch. Show `git status` and `git log --oneline <base>..HEAD`; ask the user (numbered): 1. merge locally into master · 2. keep the branch · 3. other. Install (`cp nfcore-rnasplice-setup/nfcore-rnasplice-setup.md ~/.claude/commands/`, keeping a `.bak` of any existing file) and merge only after the user approves. Never `git push` without the user's explicit approval.

---

## Self-review (planner)

1. **Spec coverage.** Purpose and decisions table: Tasks 1-8 (standalone skill, params file, explicit switches, modules incl. DTU overlap warning in Step 8 and Step 12, MISO dropped, FASTQ and BAM input, relationship to bulk-rnaseq-pipeline in Step 12 and README). Architecture (folder, README, checker, fixtures, root row, install by copy): Tasks 0, 1, 8, 10. Task 0 items 1-5: G1/G2 (launch under 26.04.6, fallbacks), G3/G8 (process names, resources, output tree, STAR/edgeR facts), G3 (typing), G5 (BAM strandedness), G6/G7 (containers, miso_genes); stop condition: G11. Wizard Steps 0-11: Tasks 1-7 (Steps 0-3, 4, 5, 6-7, 8-9, 10-11, 12). Parameter policy (typing rules, never-emitted keys, schema check): Tasks 1, 5, 6. Testing and acceptance (checker in rnavar style, need/forbid, no `--flag`, mutation proofs; cluster runs (a)/(b)/(c); honest README): Tasks 1-9. Risks (old release, BAM, read type, resources, output names): gate values and README limitations. Out of scope respected (no MISO analysis, no report skill, no `salmon_results`/`transcriptome_bam`, no iGenomes, no IsoformSwitchAnalyzeR/LeafCutter offered, no trimming/clustering tuning). One deliberate reading: the spec's "`gencode: true` for GENCODE transcript FASTA" becomes "`gencode: false` always", because the skill never passes a transcript FASTA (Planner decision 16).
2. **Placeholder scan.** No TBD/TODO. Gate-dependent text uses `«KEY»` with every key defined in the Task 0 table and verified by the checker; variant paragraphs are given in full for each outcome. The hand-off paths (Task 7) are read from the fixture with the given commands and verified by the checker.
3. **Type consistency.** Placeholders: `{SOURCE}`, `{LAYOUT}`, `{STRANDEDNESS}`, `{PAIRED_DESIGN}`, `{SHEET_PREFIX}`, `{SAMPLESHEET_CSV}`, `{CONTRASTS_CSV}`, `{READ_LENGTH}`, `{STAR_INDEX}`, `{SALMON_INDEX}`, `{RUN_*}`, `{RMATS_NOVEL}`, `{MIN_SAMPS_GENE_EXPR}`, `{MIN_SAMPS_FEATURE}`, `{MULTIQC_TITLE}`, `{OUTDIR}`, `{PARAMS_YAML}`, `{VERSION}`, `{VERSION_TAG}` are produced once and consumed with the same names; `render_params.sh` and `test_render_params.sh` use exactly the template's placeholder set. Anchors are named identically in skill text, checker and tests. Mutation runner names match `prove_mutations.sh`.
4. **Review Focus.** The five items each name the owning task and its test or checker line; all five are pinned (Tasks 2, 3, 4, 5, 9).
