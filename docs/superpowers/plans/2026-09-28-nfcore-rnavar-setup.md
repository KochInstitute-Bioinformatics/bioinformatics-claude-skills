# nfcore-rnavar-setup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an interactive `nfcore-rnavar-setup` Claude Code skill that sets up and submits nf-core/rnavar (RNA-seq short-variant calling) on SLURM + Singularity, mirroring `nfcore-rnaseq-setup`.

**Architecture:** The deliverable is one self-contained Markdown skill file (`nfcore-rnavar-setup/nfcore-rnavar-setup.md`) plus a README. Because a skill is a prompt, not code, correctness is enforced by a shell checker (`nfcore-rnavar-setup/tests/check_skill.sh`) that (a) validates every `--parameter` in the skill against the real nf-core/rnavar 1.3.0 `nextflow_schema.json`, (b) asserts required strings are present, and (c) asserts forbidden strings are absent. Each task extends the checker first (failing), then the skill text (passing). A final task runs the skill's output for real on a compute node.

**Tech Stack:** Markdown skill file, bash + grep + curl (no python on this cluster), nf-core/rnavar 1.3.0, Nextflow, Singularity, SLURM.

**Spec:** `docs/superpowers/specs/2026-09-28-nfcore-rnavar-setup-design.md`

## Global Constraints

- Pipeline: `nf-core/rnavar`, verified against tag **1.3.0** (2026-06-03). The skill fetches the latest tag at run time via `WebFetch` of `https://api.github.com/repos/nf-core/rnavar/releases/latest`; it never uses the `gh` CLI (not installed on this cluster).
- Every `--parameter` the skill emits must exist in the rnavar `nextflow_schema.json`. **Not** in the schema, never emit: `--annotation_cache`, `--gencode`, `--strandedness`.
- `--read_length` is always passed (schema default 150 is wrong for other read lengths). `--seq_platform illumina` is passed. `aligner` defaults to `star` and is omitted.
- Flags equal to their pipeline default are omitted from the generated script.
- Known sites are never assumed: either `--dbsnp/--dbsnp_tbi/--known_indels/--known_indels_tbi` are supplied or `--skip_baserecalibration` is set (rnavar does not skip BQSR automatically and errors late).
- Resource URLs (GATK bundle, Mouse Genomes Project, annotation caches) are resolved at run time and confirmed by the user; the skill embeds no URL it has not verified.
- `nextflow.config` is generated only if absent, never overwritten. `resourceLimits` cap: 64 GB memory, 24 h, 16 CPUs.
- Email (Step 1) and conda environment (Step 2) are asked in separate messages; neither is pre-filled.
- All finite-choice questions are numbered lists. No heavy computation on the login node; raw FASTQ/BAM are read-only.
- The skill file must be self-contained: no references to steps in other skills.
- Git: show `git status` and `git diff` before each commit; never `git push` without explicit user approval.

## Review Focus

1. **Known sites missing** → user gets a wizard choice, and `--skip_baserecalibration` is emitted when chosen (never a late BQSR failure). Tested in Task 3.
2. **Mixed read lengths across samples** → skill reports the distribution and uses the most common length with a warning, rather than silently using the first file. Tested in Task 2.
3. **Sample names with `-`, spaces, or duplicates** → sanitised; duplicates warned as lane-merging. Tested in Task 2.
4. **Existing STAR index with a different `sjdbOverhang`** → not reused; a separate `star_rnavar_sjdb{N-1}` index is built. Tested in Task 3.
5. **Annotation chosen but cache not on disk / no internet on compute nodes** → skill refuses to emit `--download_cache` silently and offers a pre-download job. Tested in Task 4.
6. **Existing `nextflow.config`** → left untouched; user told which process selectors to verify. Tested in Task 5.

---

## File Structure

| File | Responsibility |
|---|---|
| `nfcore-rnavar-setup/nfcore-rnavar-setup.md` | The skill prompt (Steps 0–13, Notes). Created in Task 1, extended through Task 5. |
| `nfcore-rnavar-setup/tests/check_skill.sh` | Static checker: schema-validates flags, required/forbidden strings. Created Task 1, extended each task. |
| `nfcore-rnavar-setup/README.md` | User-facing README in the same shape as `nfcore-rnaseq-setup/README.md`. Task 6. |
| `README.md` (repo root) | Add one row to the skills table. Task 6. |
| `~/.claude/commands/nfcore-rnavar-setup.md` | Installed copy for the user (same place `nfcore-rnaseq-setup.md` lives). Task 6. |

All repo paths are relative to `/net/bmc-lab3/data/bcc/yannvrb/Github_repos/bioinformatics-claude-skills`.

---

### Task 1: Checker scaffold and skill Steps 0–3

**Files:**
- Create: `nfcore-rnavar-setup/tests/check_skill.sh`
- Create: `nfcore-rnavar-setup/nfcore-rnavar-setup.md`

**Interfaces:**
- Produces: `check_skill.sh <skill.md>` exits 0 iff all checks pass; prints `FAIL: ...` lines. Helper functions `need "<string>"` and `forbid "<string>"` that later tasks append calls to. Skill file header + Steps 0–3 that later tasks append after.

- [ ] **Step 1: Write the checker**

Create `nfcore-rnavar-setup/tests/check_skill.sh`:

```bash
#!/bin/bash
# Static checks for nfcore-rnavar-setup.md. Usage: check_skill.sh <skill.md> [schema.json]
set -u
SKILL=${1:?usage: check_skill.sh <skill.md> [schema.json]}
SCHEMA=${2:-/tmp/rnavar_schema_1.3.0.json}
fail=0

if [ ! -s "$SKILL" ]; then echo "FAIL: skill file missing or empty: $SKILL"; exit 1; fi
if [ ! -s "$SCHEMA" ]; then
  curl -sfL https://raw.githubusercontent.com/nf-core/rnavar/1.3.0/nextflow_schema.json -o "$SCHEMA" \
    || { echo "cannot fetch schema"; exit 2; }
fi

need()   { grep -qF -- "$1" "$SKILL" || { echo "FAIL: missing required text: $1"; fail=1; }; }
forbid() { ! grep -qF -- "$1" "$SKILL" || { echo "FAIL: forbidden text present: $1"; fail=1; }; }

# (a) every --flag in the skill is a schema parameter or an allowlisted non-rnavar flag
schema_names=$(grep -oE '"[a-z_0-9]+": *\{' "$SCHEMA" | sed -E 's/"([a-z_0-9]+)".*/\1/' | sort -u)
allow="runMode genomeDir genomeFastaFiles sjdbGTFfile sjdbOverhang runThreadN mail-type mail-user mem"
for flag in $(grep -oE '(^|[ `(=])--[A-Za-z_0-9-]+' "$SKILL" | sed -E 's/^[^-]*--//' | sort -u); do
  if ! echo "$schema_names $allow" | tr ' ' '\n' | grep -qx -- "$flag"; then
    echo "FAIL: flag not in rnavar schema or allowlist: --$flag"; fail=1
  fi
done

# (b)/(c) content checks, appended per task -------------------------------
# --- Task 1
need "# nf-core/rnavar Pipeline Setup Skill"
need "## Step 0"
need "## Step 3"
need "api.github.com/repos/nf-core/rnavar/releases/latest"
need "Do NOT use the \`gh\` CLI"
forbid "--annotation_cache"
forbid "--gencode"
forbid "--strandedness"
# --- end Task 1

[ $fail -eq 0 ] && echo "PASS" || exit 1
```

- [ ] **Step 2: Run the checker to verify it fails**

Run: `bash nfcore-rnavar-setup/tests/check_skill.sh nfcore-rnavar-setup/nfcore-rnavar-setup.md`
Expected: `FAIL: skill file missing or empty: nfcore-rnavar-setup/nfcore-rnavar-setup.md`, exit 1.

- [ ] **Step 3: Write skill Steps 0–3**

Create `nfcore-rnavar-setup/nfcore-rnavar-setup.md`:

````markdown
# nf-core/rnavar Pipeline Setup Skill

You are helping the user set up and submit an nf-core/rnavar RNA-seq variant-calling pipeline on an HPC cluster using SLURM and Singularity. rnavar implements the GATK best-practices RNA-seq short-variant workflow: 2-pass STAR → MarkDuplicates → SplitNCigarReads → base recalibration (BQSR) → HaplotypeCaller → soft variant filtering → optional annotation (SnpEff / VEP). Walk through each step below in order, asking the user only what you need and performing automated steps silently.

---

## Step 0 — Establish working directory

Before anything else, run `pwd` to record the current working directory (`{CWD}`). All scripts, the samplesheet CSV, `nextflow.config`, and the output directory will be written here. Inform the user: "Working directory: `{CWD}`"

Also derive:
- `{WD_NAME}` = basename of `{CWD}`
- `{TODAY_YYMMDD}` = today's date formatted as `YYMMDD`

---

## Step 1 — Email address

Ask the user: "What email address should SLURM use for pipeline notifications (completion/failure)?"

Do **not** pre-suggest or pre-fill any email address. Wait for the answer before Step 2.

---

## Step 2 — Conda environment

Ask the user: "Which conda environment should be activated for Nextflow?"

Do **not** pre-suggest any environment name. Ask this as a separate question after Step 1 — never combine Steps 1 and 2 in one message.

---

## Step 3 — Latest nf-core/rnavar version

Fetch the latest release with the `WebFetch` tool against:
```
https://api.github.com/repos/nf-core/rnavar/releases/latest
```
Extract `tag_name` (strip any leading `v`) and use it as `{VERSION}`. Tell the user which version you found. The parameter names in this skill were verified against release 1.3.0; if `{VERSION}` is a newer major or minor version, fetch `https://raw.githubusercontent.com/nf-core/rnavar/{VERSION}/nextflow_schema.json` and confirm every parameter this skill emits still exists before writing any script.

**Important:** Do NOT use the `gh` CLI — it is not installed on this HPC cluster. Always use `WebFetch`.
````

- [ ] **Step 4: Run the checker to verify it passes**

Run: `bash nfcore-rnavar-setup/tests/check_skill.sh nfcore-rnavar-setup/nfcore-rnavar-setup.md`
Expected: `PASS`. (If curl cannot reach GitHub from the login node, run the check via `srun -p bcc --mem=2G -t 00:05:00 bash ...`.)

- [ ] **Step 5: Commit**

```bash
git status && git diff --stat
git add nfcore-rnavar-setup
git commit -m "Add nfcore-rnavar-setup skill skeleton (Steps 0-3) and static checker"
```
(commit message ends with the attribution line given in the session's system-reminder)

---

### Task 2: Samplesheet and read length (Steps 4–5)

**Files:**
- Modify: `nfcore-rnavar-setup/tests/check_skill.sh` (append Task 2 checks)
- Modify: `nfcore-rnavar-setup/nfcore-rnavar-setup.md` (append Steps 4–5)

**Interfaces:**
- Consumes: `need`/`forbid` helpers; placeholders `{CWD}`, `{WD_NAME}`, `{TODAY_YYMMDD}` from Task 1.
- Produces: `{SEQ_DATE}`, `{SAMPLESHEET_CSV}`, `{READ_LENGTH}` used by Tasks 3 and 5.

- [ ] **Step 1: Append failing checks**

Insert before the final `[ $fail ...` line of `check_skill.sh`:

```bash
# --- Task 2
need "## Step 4"
need "sample,fastq_1,fastq_2"
need "sample,bam,bai"
need "sample,cram,crai"
need "Supplying FASTQ files and a BAM/CRAM file for the same sample"
need "Replace every \`-\` with \`_\`"
need "merged before alignment"
need "## Step 5"
need "--read_length"
need "most common read length"
need "sjdbOverhang"
forbid "strandedness,"
# --- end Task 2
```

- [ ] **Step 2: Run to verify it fails**

Run: `bash nfcore-rnavar-setup/tests/check_skill.sh nfcore-rnavar-setup/nfcore-rnavar-setup.md`
Expected: multiple `FAIL: missing required text:` lines, exit 1.

- [ ] **Step 3: Append Steps 4–5 to the skill**

````markdown

---

## Step 4 — Raw data location and samplesheet

**Before asking the user, scan the working directory for input files:**

```bash
find {CWD} -name "*.fastq.gz" -o -name "*.fq.gz" | head -30
find {CWD} -name "*.bam" -o -name "*.cram" | head -10
```

- FASTQ files found: list the unique directories and ask "I found FASTQ files in: `{FOUND_DIRS}`. Use this directory, or specify another?"
- Only BAM/CRAM found: ask whether these are STAR-aligned, duplicate-marked BAMs (rnavar accepts them as input and skips alignment) and use the BAM samplesheet form below.
- Nothing found: ask for the full path to the input files.

**Paired-end detection:** filenames containing `_R1_`/`_R2_`, `_1.fastq.gz`/`_2.fastq.gz`, or `_1_sequence`/`_2_sequence` → paired-end. Otherwise single-end: warn "Single-end data is accepted by rnavar but gives weaker variant calls than paired-end; continue?" Detect the sequencing date from a leading 6-digit `YYMMDD` filename prefix as `{SEQ_DATE}` (fall back to today's date).

**Path style:** if the input directory is inside `{CWD}`, use paths relative to `{CWD}`; otherwise absolute paths.

**Samplesheet forms (exactly one file type per sample):**
- FASTQ: `sample,fastq_1,fastq_2` (rnavar has no `strandedness` column). Find each R1, pair it with its R2 by substituting `_R1_`→`_R2_`, `_1.`→`_2.`, or `_1_sequence`→`_2_sequence`; warn and leave `fastq_2` empty if R2 is missing. Sample name = filename up to `_S\d+`, `_R1`, or `_1_sequence`.
- BAM: `sample,bam,bai`. CRAM: `sample,cram,crai`.
- **Never mix types for one sample.** Supplying FASTQ files and a BAM/CRAM file for the same sample makes the pipeline error; check for this and stop with a clear message.

**Sample-name sanitisation (always, before showing the user):**
- Replace every `-` with `_`; replace spaces, `/`, `(`, `)` and other special characters with `_`. Note substitutions in the preview.
- Check uniqueness after sanitisation. On collision warn: "⚠️ Name collision '{NAME}': in nf-core/rnavar, rows with the same sample name are treated as lanes of one sample and their reads are merged before alignment. Provide distinct names if these are different samples." Then go to custom naming.

**Review and naming — order is mandatory:**
1. Show the full samplesheet (all rows) as a table.
2. Ask about names (numbered): 1. Use auto-generated names · 2. Provide custom names. For custom names, show numbered auto names next to filenames, ask for a plain-language description, build the mapping, show it as an auto→new table and ask "Does this mapping look correct?" Validate custom names: no `-`; no duplicates (same warning as above).
3. Ask for the samplesheet filename (numbered): 1. `{SEQ_DATE}_{WD_NAME}_samplesheet.csv` · 2. `{TODAY_YYMMDD}_{WD_NAME}_samplesheet.csv` · 3. Custom.
4. Write the file only after names and filename are confirmed. Store as `{SAMPLESHEET_CSV}`.

---

## Step 5 — Read length

rnavar uses `--read_length` to set STAR's `sjdbOverhang` (`read_length − 1`), and its default is 150, which is wrong for most other libraries. **Always detect and always pass it.**

```bash
zcat {FASTQ_FILE} | awk 'NR%4==2 {print length($0)}' | head -n 1000 | sort -n | uniq -c | sort -rn | head -3
```

Run this on the first FASTQ of several different samples (up to 5). If lengths differ between samples, report the distribution, use the most common read length as `{READ_LENGTH}`, and warn that `sjdbOverhang` is tuned to it. For BAM/CRAM input, ask the user for the read length. Tell the user: "Detected read length {READ_LENGTH} bp → sjdbOverhang {READ_LENGTH − 1}." Store `{READ_LENGTH}`.
````

- [ ] **Step 4: Run to verify it passes**

Run: `bash nfcore-rnavar-setup/tests/check_skill.sh nfcore-rnavar-setup/nfcore-rnavar-setup.md`
Expected: `PASS`.

- [ ] **Step 5: Commit**

```bash
git status && git diff --stat
git add nfcore-rnavar-setup
git commit -m "nfcore-rnavar-setup: samplesheet and read-length steps"
```

---

### Task 3: Genome, STAR index and known sites (Steps 6–7)

**Files:**
- Modify: `nfcore-rnavar-setup/tests/check_skill.sh`
- Modify: `nfcore-rnavar-setup/nfcore-rnavar-setup.md`

**Interfaces:**
- Consumes: `{READ_LENGTH}` (Task 2).
- Produces: `{GENOME_DIR}`, `{FASTA_PATH}`, `{GTF_PATH}`, `{STAR_INDEX}`, `{ORGANISM}`, `{ASSEMBLY}`, `{ENS_VERSION}`, `{KNOWN_SITES_LINES}` (a block of flag lines or a single `--skip_baserecalibration \`), `{SJDB_OVERHANG}` used by Task 5.

- [ ] **Step 1: Append failing checks**

```bash
# --- Task 3
need "## Step 6"
need "star_rnavar_sjdb"
need "GTF source"
need "no flag is emitted"
need "## Step 7"
need "--dbsnp"
need "--dbsnp_tbi"
need "--known_indels"
need "--known_indels_tbi"
need "--skip_baserecalibration"
need "does not skip base recalibration automatically"
need "resolve the resource URLs at run time"
need "--star_index"
# --- end Task 3
```

- [ ] **Step 2: Run to verify it fails**

Run: `bash nfcore-rnavar-setup/tests/check_skill.sh nfcore-rnavar-setup/nfcore-rnavar-setup.md`
Expected: `FAIL: missing required text:` lines for the Task 3 strings.

- [ ] **Step 3: Append Steps 6–7 to the skill**

````markdown

---

## Step 6 — Organism and genome files

Ask:
1. "What organism is this data from? (e.g. mouse, human)"
2. "What is the base directory where genome files and indexes are stored?"

Use this folder convention (shared with other nf-core skills so FASTA/GTF are reused, never re-downloaded):

```
{genome_base}/{organism}/{assembly}_ens{version}/
├── {FASTA}.fa                    ← primary assembly FASTA
├── {GTF}.gtf                     ← annotation GTF
└── index/
    └── star_rnavar_sjdb{N-1}/    ← STAR index built for THIS read length
```

- Mouse: assembly GRCm39, FASTA `Mus_musculus.GRCm39.dna.primary_assembly.fa`, GTF `Mus_musculus.GRCm39.{version}.gtf`, directory `{genome_base}/mouse/mm39_ens{version}/`.
- Human: assembly GRCh38, FASTA `Homo_sapiens.GRCh38.dna.primary_assembly.fa`, GTF `Homo_sapiens.GRCh38.{version}.gtf`, directory `{genome_base}/human/hg38_ens{version}/`.
- Other organisms: ask for the FASTA and GTF paths; skip the checks below.

**Existing FASTA/GTF:** if present, report the paths and reuse them. If missing, fetch the latest Ensembl release from `https://ftp.ensembl.org/pub/current_README`, and generate the download commands in the helper script (Step 12).

**GTF source:** inspect `grep -v "^#" {GTF_PATH} | head -3`. Gene IDs with a version suffix (`ENSG00000000003.15`) indicate GENCODE; without one, Ensembl. Report which it is, but **no flag is emitted** — rnavar's schema has no `gencode` parameter.

**STAR index — always its own index per read length.** `sjdbOverhang` is fixed when the index is built, and an index made for another read length (for example one built by another pipeline) would be wrong here. Set `{SJDB_OVERHANG}` = `{READ_LENGTH} − 1` and look for `{GENOME_DIR}/index/star_rnavar_sjdb{SJDB_OVERHANG}/`:
- Present and non-empty (`SA`, `Genome`, `sjdbList.out.tab` exist): use it, `{STAR_INDEX}` = that path.
- Missing: add the index build to the helper script (Step 12) and tell the user it must run first.

`--star_index '{STAR_INDEX}'` is always passed so rnavar never rebuilds the index inside the workflow.

---

## Step 7 — Known sites for base recalibration

rnavar passes the known-sites files directly to GATK BaseRecalibrator and **does not skip base recalibration automatically** — if they are missing the run fails late, after alignment. So always resolve this now. Ask (numbered):

1. **Use known-sites VCFs** — add `--dbsnp {DBSNP}`, `--dbsnp_tbi {DBSNP}.tbi`, `--known_indels {INDELS}`, `--known_indels_tbi {INDELS}.tbi`.
   - Human: the GATK resource-bundle dbSNP and Mills/1000G known-indels VCFs for the matching assembly.
   - Mouse: Mouse Genomes Project variants (SNPs and indels) for GRCm39.
   - Check whether the files already exist under `{GENOME_DIR}/known_sites/`; if so reuse them and verify the `.tbi` indexes exist.
   - If they do not exist, **resolve the resource URLs at run time**: use `WebFetch` on the current GATK resource-bundle page (human) or the Mouse Genomes Project / Ensembl variation FTP listing (mouse), show the exact URLs to the user, and download only after they confirm. Never type a URL from memory. Add the download and `tabix -p vcf` steps to the helper script (Step 12).
2. **Skip base recalibration** — add `--skip_baserecalibration`. Tell the user the trade-off: base qualities are not recalibrated, which is slightly less accurate but is the right choice for organisms without a curated variant set, or to get a first result quickly.

Store the result as `{KNOWN_SITES_LINES}`: the four `--dbsnp`/`--known_indels` lines for option 1, or the single line `--skip_baserecalibration \` for option 2.
````

- [ ] **Step 4: Run to verify it passes**

Run: `bash nfcore-rnavar-setup/tests/check_skill.sh nfcore-rnavar-setup/nfcore-rnavar-setup.md`
Expected: `PASS`.

- [ ] **Step 5: Commit**

```bash
git status && git diff --stat
git add nfcore-rnavar-setup
git commit -m "nfcore-rnavar-setup: genome, per-read-length STAR index, known sites"
```

---

### Task 4: Variant options and annotation (Steps 8–9)

**Files:**
- Modify: `nfcore-rnavar-setup/tests/check_skill.sh`
- Modify: `nfcore-rnavar-setup/nfcore-rnavar-setup.md`

**Interfaces:**
- Produces: `{VARIANT_LINES}` (possibly empty) and `{ANNOTATION_LINES}` (possibly empty), `{ANNOTATION_TOOL}` ∈ {none, snpeff, vep, merge}, used by Task 5.

- [ ] **Step 1: Append failing checks**

```bash
# --- Task 4
need "## Step 8"
need "--remove_duplicates"
need "--star_twopass"
need "--gatk_hc_call_conf"
need "--gatk_vf_qd_filter"
need "--gatk_vf_fs_filter"
need "--gatk_vf_window_size"
need "--gatk_vf_cluster_size"
need "--skip_variantfiltration"
need "--generate_gvcf"
need "--bam_csi_index"
need "## Step 9"
need "--tools"
need "--snpeff_cache"
need "--vep_cache"
need "--snpeff_db"
need "--vep_genome"
need "--vep_species"
need "--vep_cache_version"
need "--download_cache"
need "needs internet from compute nodes"
need "not a parameter in the rnavar schema"
# --- end Task 4
```

- [ ] **Step 2: Run to verify it fails**

Run: `bash nfcore-rnavar-setup/tests/check_skill.sh nfcore-rnavar-setup/nfcore-rnavar-setup.md`
Expected: `FAIL: missing required text:` for the Task 4 strings.

- [ ] **Step 3: Append Steps 8–9 to the skill**

````markdown

---

## Step 8 — Variant calling options

Ask each as a numbered choice. **Omit any flag whose value equals the pipeline default** — only emit a flag when the user changes it.

**a) Duplicates.** 1. Keep duplicates marked (default — omit) · 2. Remove duplicates — add `--remove_duplicates`.

**b) STAR two-pass.** Two-pass mapping (`--star_twopass`) is on by default and recommended for calling; keep it. Emit `--star_twopass false` only if the user explicitly asks to disable it.

**c) Calling and filtering thresholds.** Ask: 1. Pipeline defaults (recommended; emits nothing) · 2. Customise. If customising, ask for each and emit only values that differ from the default: `--gatk_hc_call_conf` (default 20), `--gatk_vf_qd_filter` (2), `--gatk_vf_fs_filter` (30), `--gatk_vf_window_size` (35), `--gatk_vf_cluster_size` (3). Offer `--skip_variantfiltration` if the user wants unfiltered calls.

**d) gVCFs.** Ask: "Will you jointly call variants across samples later?" 1. No (omit) · 2. Yes — add `--generate_gvcf`.

**e) Large chromosomes.** Only if the genome has chromosomes longer than 512 Mb (not human or mouse): add `--bam_csi_index` and tell the user it disables variant filtration.

**f) Save intermediates.** 1. No (omit) · 2. Yes — add `--save_align_intermeds` (recommended if the BAMs will feed allele-specific expression analysis).

Collect the emitted lines as `{VARIANT_LINES}`.

---

## Step 9 — Optional variant annotation

Ask (numbered): 1. No annotation (default — omit `--tools`) · 2. SnpEff · 3. VEP · 4. Both merged. Emit `--tools snpeff`, `--tools vep`, or `--tools merge` respectively. Set `{ANNOTATION_TOOL}`.

If annotation is chosen, the caches must be on disk **before** submission:
- Ask for existing cache directories → `--snpeff_cache '{DIR}'` and/or `--vep_cache '{DIR}'`, plus the matching identifiers: `--snpeff_db`, `--vep_genome`, `--vep_species`, `--vep_cache_version` (ask; do not guess versions).
- If a cache is missing, do **not** silently add `--download_cache`: that option needs internet from compute nodes, which may not be available, and the job then stalls without a clear error. Instead tell the user this, and offer to generate a pre-download helper script (Step 12) that they run where internet is available.

Note for the assistant: `annotation_cache` appears on the rnavar usage page but is not a parameter in the rnavar schema — never emit it. Use `--snpeff_cache`, `--vep_cache` and `--download_cache` only.

Collect the emitted lines as `{ANNOTATION_LINES}`.
````

- [ ] **Step 4: Run to verify it passes**

Run: `bash nfcore-rnavar-setup/tests/check_skill.sh nfcore-rnavar-setup/nfcore-rnavar-setup.md`
Expected: `PASS`. (The note about `annotation_cache` deliberately omits the leading dashes, so it passes both the flag check and `forbid "--annotation_cache"`.)

- [ ] **Step 5: Commit**

```bash
git status && git diff --stat
git add nfcore-rnavar-setup
git commit -m "nfcore-rnavar-setup: variant options and optional annotation"
```

---

### Task 5: Config, submission script, helpers, hand-off (Steps 10–13 and Notes)

**Files:**
- Modify: `nfcore-rnavar-setup/tests/check_skill.sh`
- Modify: `nfcore-rnavar-setup/nfcore-rnavar-setup.md`

**Interfaces:**
- Consumes: `{VERSION}`, `{SAMPLESHEET_CSV}`, `{FASTA_PATH}`, `{GTF_PATH}`, `{STAR_INDEX}`, `{READ_LENGTH}`, `{KNOWN_SITES_LINES}`, `{VARIANT_LINES}`, `{ANNOTATION_LINES}`, `{SJDB_OVERHANG}`, `{GENOME_DIR}`.
- Produces: `nf-core_rnavar_{VERSION}.sh`, `nextflow.config` (if absent), helper scripts, hand-off text.

- [ ] **Step 1: Append failing checks**

```bash
# --- Task 5
need "## Step 10"
need "do not overwrite it"
need "resourceLimits"
need "'.*:STAR_ALIGN'"
need "'.*:GATK4_HAPLOTYPECALLER'"
need "'.*:GATK4_BASERECALIBRATOR'"
need "'.*:GATK4_SPLITNCIGARREADS'"
need "## Step 11"
need "nextflow run nf-core/rnavar -r {VERSION}"
need "--seq_platform illumina"
need "--read_length {READ_LENGTH}"
need "--star_index '{STAR_INDEX}'"
need "## Step 12"
need "build_star_index_rnavar"
need "gunzip -c"
need "## Step 13"
need "## Notes for the assistant"
need "read-only"
# --- end Task 5
```

- [ ] **Step 2: Run to verify it fails**

Run: `bash nfcore-rnavar-setup/tests/check_skill.sh nfcore-rnavar-setup/nfcore-rnavar-setup.md`
Expected: `FAIL: missing required text:` for the Task 5 strings.

- [ ] **Step 3: Append Steps 10–13 and Notes to the skill**

````markdown

---

## Step 10 — MultiQC title, output directory, nextflow.config

**MultiQC title** (numbered): 1. `{SEQ_DATE}_{WD_NAME}` · 2. `{TODAY_YYMMDD}_{WD_NAME}` · 3. Custom. Store as `{MULTIQC_TITLE}`.
**Output directory** (numbered): same three options. Store as `{OUTDIR}`.

Check for an existing config: `ls nextflow.config`. **If it exists, do not overwrite it** — instead tell the user the selectors that matter for rnavar (below) and that they can compare them. If it does not exist, write:

```nextflow
// nextflow.config — nf-core/rnavar on SLURM + Singularity
profiles {
    slurm {
        process {
            executor = 'slurm'
            queue = 'bcc'
            cpus = 2
            memory = '8 GB'
            time = '4h'

            withName: '.*:STAR_ALIGN' {
                cpus = 8
                memory = '64 GB'
                time = '8h'
            }
            withName: '.*:GATK4_SPLITNCIGARREADS' {
                cpus = 2
                memory = '16 GB'
                time = '8h'
            }
            withName: '.*:GATK4_BASERECALIBRATOR' {
                cpus = 2
                memory = '16 GB'
                time = '8h'
            }
            withName: '.*:GATK4_HAPLOTYPECALLER' {
                cpus = 2
                memory = '16 GB'
                time = '8h'
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

params {
    max_cpus   = 16
    max_memory = '64 GB'
    max_time   = '24h'
}

process {
    resourceLimits = [
        cpus:   params.max_cpus,
        memory: params.max_memory,
        time:   params.max_time
    ]
}

timeline { enabled = true; file = "${params.outdir}/pipeline_info/execution_timeline.html" }
report   { enabled = true; file = "${params.outdir}/pipeline_info/execution_report.html"   }
trace    { enabled = true; file = "${params.outdir}/pipeline_info/execution_trace.txt"     }
dag      { enabled = true; file = "${params.outdir}/pipeline_info/pipeline_dag.svg"        }
```

rnavar's own `base.config` defines label-based resources only (`process_medium` = 6 CPU/36 GB/8 h, `process_high` = 12 CPU/72 GB/16 h), so the `withName` overrides above use regex selectors (`'.*:NAME'`) that do not depend on the workflow-name prefix. After the first run, compare the selectors with the process names in `{OUTDIR}/pipeline_info/execution_trace.txt` and adjust if any did not match.

---

## Step 11 — Generate the submission script

Write `nf-core_rnavar_{VERSION}.sh` in `{CWD}`:

```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -p bcc
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={USER_EMAIL}

module add miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit
conda activate {CONDA_ENV}
module add singularity/3.10.4

nextflow run nf-core/rnavar -r {VERSION} -c nextflow.config -profile slurm,singularity \
--input {SAMPLESHEET_CSV} \
--fasta {FASTA_PATH} \
--gtf {GTF_PATH} \
--star_index '{STAR_INDEX}' \
--read_length {READ_LENGTH} \
--seq_platform illumina \
{KNOWN_SITES_LINES}\
{VARIANT_LINES}\
{ANNOTATION_LINES}\
--multiqc_title {MULTIQC_TITLE} \
--outdir {OUTDIR}
```

Each `{..._LINES}` placeholder is either empty or one-or-more complete lines ending in ` \`. The final line has no trailing backslash. Show the full file and instruct:
```
Script written: nf-core_rnavar_{VERSION}.sh
To submit:  sbatch nf-core_rnavar_{VERSION}.sh
```
If any helper script (Step 12) was generated, list the order: helpers first, then the pipeline.

---

## Step 12 — Helper scripts (only for missing resources)

Generate only what is missing. Each is an `sbatch` script (`#SBATCH -N 1 -n 8 --mem=64G -t 8:00:00 -p bcc --mail-type=END,FAIL`), run on a compute node — never on the login node. Always use `gunzip -c file.gz > file` (never `gunzip -k`; not available on CentOS 7).

**`build_star_index_rnavar_{ASSEMBLY}_ens{VERSION_ENS}.sh`** — downloads (`wget -c`) and decompresses the FASTA and GTF if absent, then:

```bash
module add star/2.7.9a
mkdir -p "{GENOME_DIR}/index/star_rnavar_sjdb{SJDB_OVERHANG}"
STAR \
    --runMode genomeGenerate \
    --genomeDir "{GENOME_DIR}/index/star_rnavar_sjdb{SJDB_OVERHANG}" \
    --genomeFastaFiles "{FASTA_PATH}" \
    --sjdbGTFfile "{GTF_PATH}" \
    --sjdbOverhang {SJDB_OVERHANG} \
    --runThreadN 8
```

**`prepare_known_sites_{ASSEMBLY}.sh`** — for the URLs the user confirmed in Step 7: `wget -c` each file into `{GENOME_DIR}/known_sites/`, then `module add htslib` (or the site's tabix module) and `tabix -p vcf` each VCF that lacks a `.tbi`.

**`prepare_annotation_cache_{TOOL}.sh`** — only if the user has no cache: a script the user runs where internet is available, using the tool's own cache installer (`vep_install` or `snpEff download`) into the directory passed to `--vep_cache` / `--snpeff_cache`.

---

## Step 13 — Hand-off note

Print where the results will be, so later analyses can find them:
```
Outputs under {OUTDIR}/ :
  variant_calling/   filtered VCFs (per sample) and, with --generate_gvcf, gVCFs
  preprocessing/     recalibrated / duplicate-marked BAMs
  multiqc/           MultiQC report
These VCFs and BAMs are the inputs expected by the ase-pipeline skill (allele-specific expression).
```
Before printing, confirm the actual directory names against the pipeline's `docs/output` page for `{VERSION}` with `WebFetch`, and use the real names.

---

## Notes for the assistant

- **`gh` CLI is not available on this HPC cluster.** Use `WebFetch` for all GitHub API calls.
- **Always present finite-choice questions as numbered lists.** Use open questions only when no reasonable discrete set exists (email, conda env, custom paths).
- Never run heavy computation on the login node; all work through `sbatch`.
- Raw FASTQ/BAM files are read-only; never modify them.
- Omit flags that equal the pipeline default.
- Never overwrite an existing `nextflow.config`.
- `--read_length` is always emitted; a STAR index made for a different read length is never reused.
- Known sites are never assumed — either supply all four files or `--skip_baserecalibration`.
- Never embed a download URL that was not verified in this session.
- Not rnavar parameters, never emit: `annotation_cache`, `gencode`, `strandedness`.
````

- [ ] **Step 4: Run to verify it passes**

Run: `bash nfcore-rnavar-setup/tests/check_skill.sh nfcore-rnavar-setup/nfcore-rnavar-setup.md`
Expected: `PASS`. If a flag-validity failure appears for a helper-script flag, add that flag to the `allow` list in the checker (only for non-rnavar tools such as STAR) — never to silence an rnavar parameter.

- [ ] **Step 5: Commit**

```bash
git status && git diff --stat
git add nfcore-rnavar-setup
git commit -m "nfcore-rnavar-setup: config, submission script, helpers, hand-off"
```

---

### Task 6: Real-run acceptance, README, install

**Files:**
- Create: `nfcore-rnavar-setup/README.md`
- Modify: `README.md` (repo root; add one table row)
- Create (install): `/home/yannvrb/.claude/commands/nfcore-rnavar-setup.md`
- Possibly modify: `nfcore-rnavar-setup/nfcore-rnavar-setup.md` (fixes found while testing)

**Interfaces:**
- Consumes: the finished skill from Tasks 1–5.

- [ ] **Step 1: Install and dry-run the wizard**

Run: `cp nfcore-rnavar-setup/nfcore-rnavar-setup.md /home/yannvrb/.claude/commands/nfcore-rnavar-setup.md`
Then, in a fresh scratch working directory containing 2–3 small paired-end FASTQ files, invoke `/nfcore-rnavar-setup` and answer the wizard. Expected: it asks email then conda env in separate messages; finds the FASTQs; builds `sample,fastq_1,fastq_2`; reports the read length; asks the known-sites question; and writes `nf-core_rnavar_{VERSION}.sh` and `nextflow.config`. Save the generated script and config.

- [ ] **Step 2: Validate every generated flag against the schema**

Run: `grep -oE -- '--[a-z_0-9]+' nf-core_rnavar_*.sh | sort -u` and confirm each name is in `/tmp/rnavar_schema_1.3.0.json` (`grep -c '"NAME":' /tmp/rnavar_schema_1.3.0.json` ≥ 1). Expected: no unknown flag.

- [ ] **Step 3: Run the pipeline's own test profile on a compute node**

Write a scratch sbatch script (not on the login node) that runs the same conda/Singularity block as the generated script, then `nextflow run nf-core/rnavar -r {VERSION} -profile test,slurm,singularity -c nextflow.config --outdir test_out`. Submit with `sbatch -p bcc`. Expected: workflow completes; `test_out/pipeline_info/execution_trace.txt` exists.

- [ ] **Step 4: Verify the config selectors against the trace**

Run: `cut -f4 test_out/pipeline_info/execution_trace.txt | sort -u | grep -E 'STAR_ALIGN|SPLITNCIGAR|BASERECALIBRATOR|HAPLOTYPECALLER'`
Expected: each of the four selectors in `nextflow.config` matches at least one real process name. If any does not match (for example the SplitNCigarReads process has a different name), fix the selector in the skill's Step 10 template and in `check_skill.sh`'s `need` lines, re-run the checker, and repeat this step.

- [ ] **Step 5: Exercise the two failure paths**

(a) Re-run the wizard choosing "use known-sites" with a nonexistent path and confirm it stops or offers the download/skip choice instead of writing a script that would fail late. (b) Give FASTQs of two different read lengths and confirm the skill reports the distribution and passes the most common `--read_length`.

- [ ] **Step 6: Write the READMEs**

Create `nfcore-rnavar-setup/README.md` in the shape of `nfcore-rnaseq-setup/README.md` (read that file first): what it does, prerequisites, usage (`/nfcore-rnavar-setup`), the wizard steps in one table, outputs (`nf-core_rnavar_{version}.sh`, `nextflow.config`, samplesheet, helper scripts), and Known limitations (annotation caches must be pre-downloaded; known sites required or BQSR skipped; STAR index built per read length). Add to the root `README.md` skills table, after the nfcore-scrnaseq row:

```
| nf-core/rnavar setup | `/nfcore-rnavar-setup` | Interactive setup wizard for nf-core/rnavar RNA-seq variant calling (GATK best practices: 2-pass STAR, SplitNCigarReads, BQSR, HaplotypeCaller, optional SnpEff/VEP) on SLURM + Singularity |
```

- [ ] **Step 7: Final check and commit**

```bash
bash nfcore-rnavar-setup/tests/check_skill.sh nfcore-rnavar-setup/nfcore-rnavar-setup.md
git status && git diff --stat
git add nfcore-rnavar-setup README.md
git commit -m "Add nfcore-rnavar-setup README and register in skills table"
```
Expected: checker prints `PASS`. **Do not push** — ask the user first.

---

## Self-review notes

- Spec coverage: Steps 0–3 (Task 1), Step 4/5 (Task 2), Steps 6/7 (Task 3), Steps 8/9 (Task 4), Steps 10–13 + notes (Task 5), acceptance tests 1–4 of the spec (Task 6 Steps 1–5). The spec's roadmap is context, not built here.
- Consistency: placeholders defined once and reused with the same names: `{READ_LENGTH}`, `{SJDB_OVERHANG}`, `{STAR_INDEX}`, `{KNOWN_SITES_LINES}`, `{VARIANT_LINES}`, `{ANNOTATION_LINES}`.
- Known deviation from the spec, deliberate: the spec's `--tools` values list included `bcfann`; the skill omits it because it needs extra BCFtools annotation inputs and was not requested.
