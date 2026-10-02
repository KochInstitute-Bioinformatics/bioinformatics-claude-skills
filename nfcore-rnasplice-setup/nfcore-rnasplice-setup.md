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

Tell the user that this pipeline revision needs Nextflow 26.04.0 or newer (verified with Nextflow 26.04.6).

Do not run `nextflow` here to check the version: Nextflow runs only on compute nodes. The submission script (Step 11) checks the version when the job starts and stops with a clear message if it is outside the verified range.

---

## Step 3 — nf-core/rnasplice version

This skill was verified against nf-core/rnasplice 1b447239488097651d8eac44bca2c1556865eb0f with Nextflow 26.04.6; that is a commit of the development branch (`dev-1b44723`), pinned because release 1.0.4 does not run under that Nextflow (see the README).

Fetch the latest release with the `WebFetch` tool against:
```
https://api.github.com/repos/nf-core/rnasplice/releases/latest
```
Extract `tag_name` (strip a leading `v`) as `{TAG}`.
- Set `{VERSION}` = 1b447239488097651d8eac44bca2c1556865eb0f (`nextflow run -r` accepts a commit).
- If the fetch fails (no network, rate limit, or no `tag_name` in the answer): tell the user that the latest release could not be checked and keep the pinned commit.
- If `{TAG}` is 1.0.4: tell the user that the pinned commit is used (release 1.0.4 does not run under Nextflow 26.04.6).
- If `{TAG}` is newer than 1.0.4: tell the user, fetch `https://raw.githubusercontent.com/nf-core/rnasplice/{TAG}/nextflow_schema.json` with `WebFetch`, list the keys this skill writes (the Step 8 module block and the Step 11 template) that are missing from it, and ask (numbered): 1. Use the pinned commit (verified; default) · 2. Use release {TAG} (not verified by this skill; offer it only when no key is missing). If the user picks the release, `{VERSION}` = {TAG}, and warn the user that the Nextflow requirement (Step 2) and the Nextflow version range the submission script checks (Step 11) were verified for the pinned commit only, not for release {TAG}.

Set `{VERSION_TAG}` = `{VERSION}` when it is a release tag, or `dev-` followed by the first 7 characters when it is a commit; it is used in file names.

**Important:** Do NOT use the `gh` CLI — it is not installed on this HPC cluster. Always use `WebFetch`.
