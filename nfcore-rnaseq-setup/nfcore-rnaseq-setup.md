# nf-core/rnaseq Pipeline Setup Skill

You are helping the user set up and submit an nf-core/rnaseq bulk RNA-seq pipeline on an HPC cluster using SLURM and Singularity. Walk through each step below in order, asking the user only what you need and performing automated steps silently.

---

## Step 0 — Establish working directory

Before anything else, run `pwd` to record the current working directory (`{CWD}`). All scripts, the samplesheet CSV, `nextflow.config`, and the output directory will be written here. Inform the user: "Working directory: `{CWD}`"

Also derive:
- `{WD_NAME}` = basename of `{CWD}` (e.g. `Test3` from `/data/bcc/user/Test3`)
- `{TODAY_YYMMDD}` = today's date formatted as `YYMMDD` (e.g. `260513` for 2026-05-13)

This also serves as the base for relative path resolution in Step 4.

---

## Step 1 — Email address

Ask the user: "What email address should SLURM use for pipeline notifications (completion/failure)?"

Do **not** pre-suggest or pre-fill any specific email address. Wait for the user's answer before proceeding to Step 2.

---

## Step 2 — Conda environment

Ask the user: "Which conda environment should be activated for Nextflow?"

Do **not** pre-suggest or pre-fill any specific environment name. Ask this as a separate question after Step 1 — never combine Steps 1 and 2 in the same message.

---

## Step 3 — Latest nf-core/rnaseq version

Fetch the latest release using the `WebFetch` tool against this URL:
```
https://api.github.com/repos/nf-core/rnaseq/releases/latest
```
Extract the `tag_name` field (strip any leading `v`). Use this version in the script. Inform the user of the version you found.

**Important:** Do NOT use the `gh` CLI — it is not installed on this HPC cluster. Always use `WebFetch` directly.

---

## Step 4 — Raw data location and samplesheet

**Before asking the user, scan the current working directory for FASTQ files:**

```bash
find {CWD} -name "*.fastq.gz" | head -30
```

- If `.fastq.gz` files are found, identify the unique directories containing them and present them as suggestions. Ask: "I found FASTQ files in the following location(s): `{FOUND_DIRS}`. Is this the directory you'd like to use, or would you like to specify a different path?"
- If no files are found, ask: "Where do your raw FASTQ files live? Please provide the full path to the directory containing the `.fastq.gz` files."

Once you have the path, **auto-detect** layout, sequencing date, and read length by scanning the directory:

```bash
find {FASTQ_DIR} -name "*.fastq.gz" | head -30
```

**Paired-end detection:** If filenames contain `_R1_`, `_R2_`, `_1.fastq.gz`, `_2.fastq.gz`, or `_1_sequence.fastq.gz` / `_2_sequence.fastq.gz` patterns → paired-end.

**Single-end detection:** If no R1/R2 patterns are found → single-end.

**Sequencing date detection:** Inspect filenames for a leading 6-digit prefix matching `YYMMDD` (e.g. `250430` in `250430Yil_D25-176064_NA_sequence.fastq.gz`). Store as `{SEQ_DATE}`. Fall back to today's date if not found.

Report detected layout and ask the user to confirm. Proceed with their answer if they correct you.

**Path style in the samplesheet:**
- If the FASTQ directory is inside `{CWD}` (i.e. `{FASTQ_DIR}` starts with `{CWD}/`), use **paths relative to `{CWD}`** (e.g. `data/sample_1_sequence.fastq.gz`).
- If the FASTQ directory is outside `{CWD}`, use full absolute paths.

**Samplesheet construction:**

*Paired-end:*
- Find all R1 files recursively. For each R1, find its R2 counterpart by substituting `_R1_` → `_R2_`, `_1.` → `_2.`, or `_1_sequence` → `_2_sequence`. If R2 is missing for a sample, warn the user and leave `fastq_2` empty.
- Derive the sample name by keeping everything up to `_S\d+`, `_R1`, or `_1_sequence` in the filename.
- CSV columns: `sample,fastq_1,fastq_2,strandedness`

*Single-end:*
- Find all `.fastq.gz` files recursively.
- Derive the sample name by stripping lane/index suffixes from the filename.
- CSV columns: `sample,fastq_1,strandedness`

In both cases use full absolute paths and set strandedness to `auto`.

**Sample name sanitisation (always apply before showing the user):**
- Replace every `-` with `_` in auto-generated sample names. Dashes cause parsing errors in R. Do this silently and note any substitutions in the preview.
- Replace any other special characters (spaces, `/`, `(`, `)`, etc.) with `_`.
- Check that all auto-generated names are unique after sanitisation. If any collisions occur, warn: "⚠️ Auto-generated name collision: '{NAME}' was derived from multiple files. In nf-core/rnaseq, samples with the same name are merged as sequencing lanes. Please provide custom names to distinguish them." Then proceed to the custom naming flow.

**Samplesheet review and naming — order is mandatory:**

1. **Show the full samplesheet** (all rows, not just a preview) in a table so the user can review every sample.

2. **Ask about sample names** — present as numbered options:
   1. Use auto-generated names — proceed
   2. Provide custom names

   - If the user chooses **custom names**: display the full numbered list of auto-generated names alongside their FASTQ filenames. Ask the user to describe the naming in plain language (e.g. "samples 1–5 are CTRL replicates, 6–9 are TMZ"). Interpret the description and construct the full mapping. If ambiguous or incomplete, ask a follow-up for only the unresolved samples.

     Display the mapping as a two-column table (auto → new name) and ask: "Does this mapping look correct?" Accept corrections until confirmed.

     **Validation (apply to all custom names before writing):**
     - Any name containing `-`: warn "Sample names cannot contain dashes — this causes issues in R. Please replace `-` with `_`."
     - Any duplicate name: warn "⚠️ Duplicate sample name '{NAME}' detected. In nf-core/rnaseq, samples sharing the same name are treated as technical replicates from different lanes and their reads will be **merged before alignment**. If these are distinct biological samples, this will silently combine their data. Please provide unique names."

3. **Ask for the samplesheet filename** — present as numbered options:
   1. `{SEQ_DATE}_{WD_NAME}_samplesheet.csv` — sequencing date (YYMMDD)
   2. `{TODAY_YYMMDD}_{WD_NAME}_samplesheet.csv` — today's date (YYMMDD)
   3. Custom filename — ask the user to type it

4. **Write the file** to the current working directory only after names and filename are confirmed.

---

## Step 5 — Trimming

Present as numbered options:
1. Data is already trimmed — add `--skip_trimming`; skip trimmer and min-reads questions
2. Trim with `fastp` — faster, modern, recommended
3. Trim with `trimgalore` — classic, widely used default

If trimming (options 2 or 3): add `--trimmer '{TRIMMER}'`.

Then ask: "What is the minimum number of reads a sample must retain after trimming?" Present as numbered options:
1. Default (10,000) — omit the flag
2. Custom value — add `--min_trimmed_reads {N}`

---

## Step 5b — rRNA removal

Ask: "Would you like to remove ribosomal RNA reads before alignment?" (Useful when the library is **not** ribo-depleted at the wet lab stage.) Present as numbered options:
1. No — omit
2. Yes — add `--remove_ribo_rna`

---

## Step 6 — 3' DGE / UMI detection

**If paired-end layout was detected in Step 4: skip this entire step.** Paired-end data is treated as standard RNA-seq — omit all 3' DGE flags and proceed to Step 7.

**If single-end layout was detected**, run the following checks silently.

### Detection procedure

**Check 1 — Filename keywords**
Scan all file paths for: `quantseq`, `umi`, `3dge`, `3prime`, `brbseq`, `brb-seq`, `marsseq`, `mars-seq`, `dprime`, `lexogen` (case-insensitive). Match = strong evidence.

**Check 2 — Read length** (measure now if not yet detected):
```bash
zcat {FASTQ_FILE} | awk 'NR%4==2 {print length($0)}' | head -n 100 | sort | uniq -c | sort -rn | head -3
```
- ≤ 20 bp: strong evidence.
- 20–50 bp: mild evidence.
- > 50 bp: no strong signal (QuantSeq can produce 76–100 bp reads).

**Check 3 — PolyT content**
```bash
zcat {FASTQ_FILE} | awk 'NR%4==2' | head -n 50
```
≥ 30% of reads starting with `^T{4,}` → strong evidence.

### Scoring and user confirmation

Always ask the user to confirm:
- **High confidence** (2+ strong signals): "This looks like a **3' DGE experiment** ([reason]). Is that correct?"
- **Possible** (1 weak signal): "I noticed [reason] — is this standard RNA-seq or a 3' DGE experiment (e.g. QuantSeq, BRB-seq, MARS-seq)?"
- **No signals**: "I found no 3' DGE indicators. Is this standard RNA-seq, or a 3' DGE experiment?"

### Outcome

**If confirmed standard RNA-seq:**
Omit all 3' DGE flags. Check `nextflow.config` for misconfigured length correction:
```bash
grep "extra_salmon_quant_args" nextflow.config
```
- Active and uncommented → warn: "⚠️ `--noLengthCorrection` is enabled but this is standard RNA-seq. This biases TPM toward longer transcripts. Commenting it out now." Comment it out automatically.
- Absent or commented → no action.

**If confirmed 3' DGE:**

1. Always add `--skip_stringtie` (skip Step 11).

2. Check and fix `nextflow.config`:
   - Active and uncommented → tell user it is already correct.
   - Commented out → uncomment automatically. Tell user: "Re-enabled `--noLengthCorrection` in `nextflow.config`."
   - Absent entirely → add `extra_salmon_quant_args = '--noLengthCorrection'` to `params {}` automatically. Tell user: "Added `--noLengthCorrection` to `nextflow.config`."

3. Ask: "Does this library include UMIs?" Present as numbered options:
   1. No UMIs — omit UMI flags
   2. Yes — proceed to UMI questions

   If yes:
   a. "How are the UMIs structured?" Present as numbered options:
      1. Separate UMI read (R1 = UMI-only, R2 = RNA): add `--with_umi --umi_discard_read 1`
      2. Embedded inline in the read (single-end): ask for pattern (e.g. `NNNNNN`), add `--with_umi --umitools_bc_pattern '{PATTERN}' --skip_markduplicates`
   b. "Which UMI deduplication tool?" Present as numbered options (omit flag if choosing default):
      1. `umitools` — default
      2. `umicollapse` — faster for large datasets
   c. "Which UMI grouping method?" Present as numbered options (omit flag if choosing default):
      1. `directional` — default; handles UMI sequencing errors
      2. `unique` — strict exact match
      3. `cluster` — groups by edit distance
      4. `adjacency` — frequency-aware clustering

---

## Step 7 — Organism and genome files

Ask:
1. "What organism is this data from? (e.g. mouse, human, rat, zebrafish)"
2. "What is the base directory where genome files and indexes should be stored?"

**GTF source — auto-detect, do not ask the user:**
If a GTF already exists in the genome directory, inspect it:
```bash
grep -v "^#" {GTF_PATH} | head -3
```
- Gene IDs like `ENSG...` / `ENSMUSG...` (no version suffix) → Ensembl, no flag needed.
- Gene IDs like `ENSG00000000003.15` (with version suffix) → GENCODE, add `--gencode`.
- GTF not yet downloaded → Ensembl by default (no flag).

### Unified folder structure

```
{genome_base}/{organism}/{assembly}_ens{version}/
├── {FASTA}.fa                  ← primary assembly FASTA
├── {FASTA}.fa.gz               ← compressed original (kept)
├── {GTF}.gtf                   ← annotation GTF
├── {GTF}.gtf.gz
├── {CDNA}.cdna.all.fa          ← raw Ensembl cDNA FASTA
├── {CDNA}.cdna.all.fa.gz
├── {CDNA}.cdna.all.noversionnum.fa   ← version-stripped cDNA
├── {CDNA}.cdna.gtf_filtered.fa       ← GTF-filtered + version-stripped (used for Salmon index)
├── gtf_transcript_ids.txt      ← transcript IDs from GTF (used for filtering)
├── decoys.txt
├── gentrome.fa
└── index/
    ├── star/
    ├── salmon/
    ├── bowtie2/     ← only if aligner is bowtie2_salmon
    └── kallisto/    ← only if pseudo_aligner is kallisto
```

### Mouse (Mus musculus):
- Assembly: GRCm39 | Directory: `{genome_base}/mouse/mm39_ens{version}/`
- FASTA: `Mus_musculus.GRCm39.dna.primary_assembly.fa`
- GTF: `Mus_musculus.GRCm39.{version}.gtf`
- cDNA: `Mus_musculus.GRCm39.cdna.all.fa`

### Human (Homo sapiens):
- Assembly: GRCh38 | Directory: `{genome_base}/human/hg38_ens{version}/`
- FASTA: `Homo_sapiens.GRCh38.dna.primary_assembly.fa`
- GTF: `Homo_sapiens.GRCh38.{version}.gtf`
- cDNA: `Homo_sapiens.GRCh38.cdna.all.fa`

### Genome version check procedure:
1. Fetch latest Ensembl release: `https://ftp.ensembl.org/pub/current_README`
2. Check: `ls {genome_base}/{organism}/` — find highest `ens{N}` present.
3. If latest version **present** and all files + index subdirs exist: report paths and proceed.
4. If **missing or outdated**: detect read length now (needed for `sjdbOverhang` in the STAR index build script — `sjdbOverhang = read_length − 1`):
   ```bash
   zcat {FASTQ_FILE} | awk 'NR%4==2 {print length($0)}' | head -n 100 | sort | uniq -c | sort -rn | head -3
   ```
   Store as `{READ_LENGTH}`. Default to 99 if detection fails. Inform the user: "Detected read length: {READ_LENGTH} bp → sjdbOverhang will be set to {READ_LENGTH − 1}."
   Then generate the index-build script (Step 14b) and create dirs with `mkdir -p`.

   **Note:** Read length detection is deferred to this point because if the index already exists it is irrelevant — the `sjdbOverhang` is baked into the index at build time and cannot be changed without a full rebuild. Re-using an existing index across runs with different read lengths works fine.

Store resolved paths for the submission script.

---

## Step 7b — nextflow.config check and generation

Check whether `nextflow.config` exists in the current working directory:

```bash
ls nextflow.config
```

- **If it exists**: do not overwrite it. The `--noLengthCorrection` check is handled in Step 6.
- **If it does not exist**: generate one automatically using the template below and inform the user.

Set `{SALMON_MEMORY}` based on organism:
- Mouse → `16 GB`
- Human → `40 GB`
- Other / unknown → `40 GB` (conservative default for large genomes)

```nextflow
// nextflow.config

profiles {
    slurm {
        process {
            executor = 'slurm'
            queue = 'bcc'
            cpus = 2
            memory = '8 GB'
            time = '4h'

            withName: 'NFCORE_RNASEQ:RNASEQ:.*:STAR_ALIGN' {
                cpus = 8
                memory = '64 GB'
                time = '8h'
            }

            withName: 'NFCORE_RNASEQ:RNASEQ:.*:SALMON_QUANT' {
                cpus = 4
                memory = '{SALMON_MEMORY}'
                time = '8h'
            }

            withName: 'NFCORE_RNASEQ:RNASEQ:.*:FASTQC' {
                cpus = 2
                memory = '4 GB'
                time = '2h'
            }

            withName: 'NFCORE_RNASEQ:RNASEQ:.*:MULTIQC' {
                cpus = 1
                memory = '8 GB'
                time = '2h'
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
    max_memory = '128 GB'
    max_time   = '24h'

    // extra_salmon_quant_args = '--noLengthCorrection'  // 3' DGE only — do not enable for standard RNA-seq
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

---

## Step 8 — MultiQC report title

Ask: "What title would you like for the MultiQC report?" Present as numbered options:
1. `{SEQ_DATE}_{WD_NAME}` — sequencing date (YYMMDD)
2. `{TODAY_YYMMDD}_{WD_NAME}` — today's date (YYMMDD)
3. Custom title — ask the user to type it

---

## Step 9 — Aligner

Options:
1. `star_salmon` — STAR + Salmon (recommended)
2. `star_rsem` — STAR + RSEM (ENCODE-compatible)
3. `hisat2` — HISAT2 only (lower memory; skip Step 10)
4. `bowtie2_salmon` — Bowtie2 + Salmon (lower memory than STAR)

---

## Step 10 — Pseudo-aligner

*(Skip if aligner is `hisat2`.)*

Options:
1. `salmon` — recommended alongside STAR
2. `kallisto`
3. None — omit `--pseudo_aligner`

---

## Step 11 — Skip StringTie

*(Skip entirely if experiment is 3' DGE — always omitted in that case.)*

Ask: "Would you like to skip StringTie transcript assembly?" Present as numbered options:
1. Yes, skip StringTie (recommended for standard quantification) — add `--skip_stringtie`
2. No, run StringTie — omit flag

---

## Step 12 — Additional output and QC options

**a) BAM files** — "Would you like to save intermediate BAM files?" Present as numbered options:
1. No — omit
2. Yes — add `--save_align_intermeds`

**b) QC steps** — "Would you like to skip all QC steps except MultiQC?" Present as numbered options:
1. No, run all QC steps — omit
2. Yes, skip QC — add `--skip_qc` (skips FastQC, RSeQC, dupRadar, Qualimap, preseq, DESeq2 QC, biotype QC)

---

## Step 13 — Output directory

Ask: "What should the output directory be named?" Present as numbered options:
1. `{SEQ_DATE}_{WD_NAME}` — sequencing date (YYMMDD)
2. `{TODAY_YYMMDD}_{WD_NAME}` — today's date (YYMMDD)
3. Custom name — ask the user to type it

---

## Step 14a — Generate the submission script

Write `nf-core_rnaseq_{version}.sh` in the current working directory:

```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -p bcc
#SBATCH --mail-type=END
#SBATCH --mail-user={USER_EMAIL}

module add miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit
conda activate {CONDA_ENV}
module add singularity/3.10.4

nextflow run nf-core/rnaseq -r {VERSION} -c nextflow.config -profile slurm,singularity \
--input {SAMPLESHEET_CSV} \
--fasta {FASTA_PATH} \
--gtf {GTF_PATH} \
{INDEX_LINES}\
--multiqc_title {MULTIQC_TITLE} \
--aligner '{ALIGNER}' \
{PSEUDO_ALIGNER_LINE}\
{TRIMMER_LINE}\
{MIN_TRIMMED_READS_LINE}\
{SKIP_TRIMMING_LINE}\
{REMOVE_RIBO_RNA_LINE}\
{UMI_LINES}\
{SKIP_STRINGTIE_LINE}\
{GENCODE_LINE}\
{SAVE_ALIGN_INTERMEDS_LINE}\
{SKIP_QC_LINE}\
--outdir {OUTDIR}
```

**Index lines by aligner:**
- `star_salmon`: `--star_index '{GENOME_DIR}/index/star' \ --salmon_index '{GENOME_DIR}/index/salmon' \`
- `star_rsem`: `--star_index '{GENOME_DIR}/index/star' \ --rsem_index '{GENOME_DIR}/index/rsem' \`
- `hisat2`: `--hisat2_index '{GENOME_DIR}/index/hisat2' \`
- `bowtie2_salmon`: `--bowtie2_index '{GENOME_DIR}/index/bowtie2' \ --salmon_index '{GENOME_DIR}/index/salmon' \`

**Other substitutions:**
- Omit any flag that matches its pipeline default (e.g. `--umi_dedup_tool umitools`, `--umitools_grouping_method directional`, `--min_trimmed_reads 10000`).
- `{UMI_LINES}`: `--with_umi` + one of `--umi_discard_read 1` or `--umitools_bc_pattern '{PATTERN}'` + `--skip_markduplicates` + optional tool/method flags. Omit entirely if no UMIs.

After writing, show the full file contents and instruct:
```
Script written: nf-core_rnaseq_{version}.sh
To submit:  sbatch nf-core_rnaseq_{version}.sh
```

If indexes needed building:
```
NOTE: Genome indexes are not yet built. Build them first:
  sbatch build_genome_index_{assembly}_ens{version}.sh
```

---

## Step 14b — Genome index build script (only when indexes are missing)

Generate `build_genome_index_{assembly}_ens{version}.sh`. Rules:

- **`gunzip`**: always use `gunzip -c file.gz > file`, never `gunzip -k` (not available on CentOS 7).
- **Modules**: load `star` and `salmon` directly — do not use nf-env conda or Singularity.
- **Salmon index must use a GTF-filtered, version-stripped cDNA FASTA.** Two issues must both be fixed:
  1. *Version mismatch*: Ensembl cDNA IDs have version suffixes (e.g. `ENSMUST00000082392.1`) but the GTF does not. Strip them with `awk`, or `CUSTOM_TX2GENE` fails with "No attribute in GTF matching transcripts".
  2. *Extra transcripts*: Ensembl cDNA includes transcripts absent from the GTF (patches, alt sequences). These get quantified but have no tx2gene entry, causing the R tximport step to fail with "No column contains all vector entries". Filter the cDNA to GTF-annotated transcripts only.
- **`sjdbOverhang`** = `{READ_LENGTH} − 1` (detected in Step 7 when index was found missing). Default to 99 if unknown.

```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=64G
#SBATCH -t 8:00:00
#SBATCH -p bcc
#SBATCH --mail-type=END
#SBATCH --mail-user={USER_EMAIL}
#SBATCH -o build_genome_index_{assembly}_ens{version}.log

module add star/2.7.9a
module add salmon/1.10.0
# module add bowtie2/...  ← uncomment if using bowtie2_salmon
# module add kallisto/... ← uncomment if using kallisto pseudo-aligner

GENOME_DIR="{GENOME_DIR}"
FASTA="${GENOME_DIR}/{FASTA_FILENAME}"
GTF="${GENOME_DIR}/{GTF_FILENAME}"
CDNA="${GENOME_DIR}/{SPECIES}.cdna.all.fa"
CDNA_NOVERSION="${GENOME_DIR}/{SPECIES}.cdna.all.noversionnum.fa"
CDNA_FILTERED="${GENOME_DIR}/{SPECIES}.cdna.gtf_filtered.fa"

# --- Step 1: Download (wget -c skips already-downloaded files) ---
cd "${GENOME_DIR}"
wget -c "{ENSEMBL_FASTA_URL}" -O {FASTA_FILENAME}.gz
wget -c "{ENSEMBL_GTF_URL}"   -O {GTF_FILENAME}.gz
wget -c "{ENSEMBL_CDNA_URL}"  -O {SPECIES}.cdna.all.fa.gz

# --- Step 2: Decompress (keep .gz originals) ---
gunzip -c {FASTA_FILENAME}.gz        > "${FASTA}"
gunzip -c {GTF_FILENAME}.gz          > "${GTF}"
gunzip -c {SPECIES}.cdna.all.fa.gz  > "${CDNA}"

# --- Step 3: Build STAR index ---
STAR \
    --runMode genomeGenerate \
    --genomeDir "${GENOME_DIR}/index/star" \
    --genomeFastaFiles "${FASTA}" \
    --sjdbGTFfile "${GTF}" \
    --sjdbOverhang {SJDB_OVERHANG} \
    --runThreadN 32

# --- Step 4: Build Salmon decoy-aware index ---
# 4a: Strip Ensembl version numbers from cDNA transcript IDs to match the GTF
#     e.g. ENSMUST00000082392.1 → ENSMUST00000082392
awk '/^>/ {sub(/\.[0-9]+/, "")} {print}' "${CDNA}" > "${CDNA_NOVERSION}"

# 4b: Filter to only transcripts annotated in the GTF
#     Prevents tximport failure from cDNA transcripts absent in the tx2gene table
grep -v "^#" "${GTF}" | grep -oP 'transcript_id "\K[^"]+' | sort -u \
    > "${GENOME_DIR}/gtf_transcript_ids.txt"

awk 'BEGIN { while ((getline line < "'"${GENOME_DIR}/gtf_transcript_ids.txt"'") > 0) ids[line] = 1 }
     /^>/ { split($1, a, " "); id = substr(a[1], 2); keep = (id in ids) }
     keep { print }' "${CDNA_NOVERSION}" > "${CDNA_FILTERED}"

grep "^>" "${FASTA}" | cut -d " " -f 1 | sed 's/>//' > "${GENOME_DIR}/decoys.txt"
cat "${CDNA_FILTERED}" "${FASTA}" > "${GENOME_DIR}/gentrome.fa"

salmon index \
    -t "${GENOME_DIR}/gentrome.fa" \
    -d "${GENOME_DIR}/decoys.txt" \
    -i "${GENOME_DIR}/index/salmon" \
    -p 32

echo "Index build complete."
```

---

## Notes for the assistant

- **`gh` CLI is not available on this HPC cluster.** Never attempt `gh api` or any `gh` command. Use `WebFetch` for all GitHub API calls.
- **Always present questions with a finite set of answers as a numbered list of options.** Only use open-ended text questions when no reasonable discrete set exists (e.g. free-text email, conda env name, or custom paths the user must supply themselves).
- Never run heavy computation on the login node. All SLURM work goes through `sbatch`.
- Raw FASTQ files are read-only; never modify them.
- For organisms other than mouse or human, ask the user to supply FASTA and GTF paths manually and skip the genome check step.
- For 3' DGE: always add `--skip_stringtie` without asking (Step 11 is skipped).
- `--umi_discard_read 1` is only valid for **paired-end** data where R1 is a UMI-only read. Never apply it to single-end data (it discards the only read).
- The pseudo-aligner (Step 10) is skipped when aligner is `hisat2`.
- `--skip_markduplicates` must always accompany `--with_umi`.
- Omit flags that match pipeline defaults to keep scripts clean.
- SALMON_QUANT memory: mouse = 16 GB, human = 40 GB, other = 40 GB. The generated `nextflow.config` reflects this.
