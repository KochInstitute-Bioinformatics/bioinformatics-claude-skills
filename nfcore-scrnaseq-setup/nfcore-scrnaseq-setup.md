# nf-core/scrnaseq Pipeline Setup Skill

You are helping the user set up and submit an nf-core/scrnaseq single-cell RNA-seq pipeline on an HPC cluster using SLURM and Singularity. Walk through each step below in order, asking the user only what you need and performing automated steps silently.

This skill sets up **one aligner per run** (like the bulk nf-core/rnaseq skill). If the user wants to compare aligners, run the skill once per aligner — but always from a **separate launch directory per run** (see Notes), because concurrent Nextflow runs launched from the same directory collide on the `.nextflow` session lock.

---

## Step 0 — Establish working directory

Before anything else, run `pwd` to record the current working directory (`{CWD}`). All scripts, the samplesheet CSV, `nextflow.config`, and the output directory will be written here. Inform the user: "Working directory: `{CWD}`"

Also derive:
- `{WD_NAME}` = basename of `{CWD}` (e.g. `Test3` from `/data/bcc/user/Test3`)
- `{TODAY_YYMMDD}` = today's date formatted as `YYMMDD` (e.g. `260707` for 2026-07-07)

---

## Step 1 — Email address

Ask the user: "What email address should SLURM use for pipeline notifications (completion/failure)?"

Do **not** pre-suggest or pre-fill any specific email address. Wait for the user's answer before proceeding to Step 2.

---

## Step 2 — Conda environment

Ask the user: "Which conda environment should be activated for Nextflow?"

Do **not** pre-suggest or pre-fill any specific environment name. Ask this as a separate question after Step 1 — never combine Steps 1 and 2 in the same message.

---

## Step 3 — Latest nf-core/scrnaseq version

Fetch the latest release using the `WebFetch` tool against this URL:
```
https://api.github.com/repos/nf-core/scrnaseq/releases/latest
```
Extract the `tag_name` field (strip any leading `v`). Use this version in the script. Inform the user of the version you found.

**Version requirements — enforce silently:**
- **10x 3' v4 (GEM-X) chemistry requires scrnaseq ≥ 3.0.0.** Earlier releases have no `10XV4` protocol and mis-assign barcodes (<10% valid). If the detected chemistry is v4 (Step 5) and the resolved version is < 3.0.0, warn the user and stop.
- The `simpleaf` aligner was named **`alevin`** before 3.0.0. On ≥ 3.0.0 use `simpleaf`; if the user is pinned to an older version, use `alevin` instead.

**Important:** Do NOT use the `gh` CLI — it is not installed on this HPC cluster. Always use `WebFetch` directly.

---

## Step 4 — Raw data location and samplesheet

**Before asking the user, scan the current working directory for FASTQ files:**

```bash
find {CWD} -name "*.fastq.gz" | head -30
```

- If `.fastq.gz` files are found, identify the unique directories containing them and present them as suggestions. Ask: "I found FASTQ files in the following location(s): `{FOUND_DIRS}`. Is this the directory you'd like to use, or would you like to specify a different path?"
- If no files are found, ask: "Where do your raw FASTQ files live? Please provide the full path to the directory containing the `.fastq.gz` files."

### 10x read-role detection (barcode vs cDNA)

10x Chromium data has, per sample and lane, a **barcode+UMI read (R1, short)** and a **cDNA read (R2, long)**, plus optional index reads (`I1`/`I2`) which are **not used**. In the nf-core/scrnaseq samplesheet, `fastq_1` **must be R1 (barcode+UMI)** and `fastq_2` **must be R2 (cDNA)** — never swap them.

Assign roles by **measured read length**, not by filename alone (filenames are usually reliable but confirm at least once):

```bash
# for a candidate R1 and R2 file
zcat {FILE} | awk 'NR%4==2 {print length($0)}' | head -100 | sort | uniq -c | sort -rn | head -3
```
- ~26–28 bp read → **barcode read → `fastq_1`** (16 bp cell barcode + 10/12 bp UMI)
- ~90 bp read (or any read markedly longer than the barcode read) → **cDNA read → `fastq_2`**
- Exclude any `_I1_` / `_I2_` index files.

Store both read lengths — the barcode-read length feeds chemistry detection in Step 5.

### Sequencing date detection

Inspect filenames for a leading 6-digit `YYMMDD` prefix (e.g. `250430` in `250430Yil_...`). Store as `{SEQ_DATE}`. Fall back to today's date if not found.

### Samplesheet construction

CSV columns: `sample,fastq_1,fastq_2` (optionally add `expected_cells` — see Step 9).

- Find all R1 (barcode) files recursively. For each, find its R2 (cDNA) counterpart by substituting `_R1_` → `_R2_` (or `_1.` → `_2.`). If R2 is missing, warn the user and do **not** proceed with that sample.
- Derive the sample name by keeping everything up to `_S\d+`, `_R1`, or the lane token `_L\d+` in the filename.

**Path style:**
- If the FASTQ directory is inside `{CWD}` (starts with `{CWD}/`), use paths **relative to `{CWD}`**.
- Otherwise use full absolute paths.

**Multi-lane handling — note the difference from bulk RNA-seq:**
10x samples are frequently split across lanes (`_L001_`, `_L002_`, …). In nf-core/scrnaseq, **multiple rows sharing the same `sample` name are the correct, intended way to feed lanes** — they are concatenated per sample before counting. This is the *desired* behaviour here, the opposite of the "accidental merge" warning in the bulk skill. Emit one row per lane, all with the same sample name. Do **not** warn about duplicate sample names when they are distinct lanes of the same library; only warn if two genuinely different biological samples resolve to the same name.

**Sample name sanitisation (always apply before showing the user):**
- Replace every `-` with `_` (dashes break downstream R / Seurat). Do this silently and note substitutions in the preview.
- Replace spaces, `/`, `(`, `)`, and other special characters with `_`.

**Samplesheet review and naming — order is mandatory:**

1. **Show the full samplesheet** (all rows) in a table so the user can review every sample and lane.
2. **Ask about sample names** — numbered options:
   1. Use auto-generated names — proceed
   2. Provide custom names — display the numbered list of auto-generated names alongside their FASTQ filenames, ask the user to describe the naming in plain language, construct and display the mapping (auto → new) as a two-column table, and accept corrections until confirmed. Validate: reject any name containing `-`; for genuinely distinct libraries reject duplicates (but remember same-name lanes are legitimate).
3. **Ask for the samplesheet filename** — numbered options:
   1. `{SEQ_DATE}_{WD_NAME}_samplesheet.csv` — sequencing date (YYMMDD)
   2. `{TODAY_YYMMDD}_{WD_NAME}_samplesheet.csv` — today's date (YYMMDD)
   3. Custom filename
4. **Write the file** to the current working directory only after names and filename are confirmed.

---

## Step 5 — Chemistry / protocol detection

Auto-detect the 10x chemistry, then **always ask the user to confirm**. Store the result as `{PROTOCOL}`.

### Detection procedure

**Barcode-read length** (from Step 4):
- **26 bp** (16 CB + 10 UMI) → 10x 3' **v2** → `10XV2`
- **28 bp** (16 CB + 12 UMI) → 10x 3' **v3 or v4** — geometry is identical, so length alone cannot separate them. Disambiguate by whitelist membership (below).
- Other lengths / 5' libraries → present the explicit list and ask.

**v3 vs v4 disambiguation (whitelist membership):**
v3 and v4 share barcode geometry but use different barcode whitelists (v3 = `3M-february-2018`, ~6.8M barcodes; v4/GEM-X = a larger ~7.4M-barcode list). Sample the first 16 bp of ~10,000 barcodes from the R1 file and measure the fraction present in each whitelist; the higher match wins.

```bash
# first 16 bp of barcodes
zcat {R1_FILE} | awk 'NR%4==2 {print substr($0,1,16)}' | head -10000 | sort -u > /tmp/bc_sample.txt
# compare against each whitelist you have on disk (adjust paths)
comm -12 <(sort /tmp/bc_sample.txt) <(sort {V3_WHITELIST}) | wc -l
comm -12 <(sort /tmp/bc_sample.txt) <(sort {V4_WHITELIST}) | wc -l
```
Both whitelists ship gzip'd in the pipeline assets and can be extracted the same way (see Step 11b), e.g.:
```bash
git -C ~/.nextflow/assets/nf-core/scrnaseq show {VERSION}:assets/whitelist/10x_V3_barcode_whitelist.txt.gz | gunzip -c > {CWD}/refs/10x_V3_whitelist.txt
git -C ~/.nextflow/assets/nf-core/scrnaseq show {VERSION}:assets/whitelist/10x_V4_barcode_whitelist.txt.gz | gunzip -c > {CWD}/refs/10x_V4_whitelist.txt
```
If neither whitelist is available, fall back to asking, defaulting the suggestion to **v4 (GEM-X)** for recently sequenced libraries.

### Confirm with the user

Report the finding and ask to confirm, e.g.:
> "Detected **10x 3' v4 (GEM-X)** — 16 bp cell barcode + 12 bp UMI, barcodes match the v4 whitelist (`{N}`/10000) far better than v3 (`{M}`/10000). I'll use `--protocol 10XV4`. Correct? (If not, pick: 1. 10XV2  2. 10XV3  3. 10XV4  4. 5-prime / other)"

### Protocol value per aligner (resolved in Step 6)

- **cellranger**: always use `--protocol auto`. Cell Ranger auto-detects chemistry internally and is the most robust; do not force a protocol for it.
- **star / simpleaf**: use the detected value directly (`10XV2` / `10XV3` / `10XV4`).
- **kallisto** with **v4**: kb-python (0.28.2, shipped with the pipeline) does **not** know `10XV4`. Run kallisto as **`--protocol 10XV3`** (identical geometry) and inject the v4 whitelist via `kb count -w` through a small extra config (Step 11b). For v2/v3, use the detected value normally.

---

## Step 6 — Aligner

Present as numbered options. Include the one-line guidance so the user can choose deliberately:

1. **`cellranger`** — 10x's reference counter (STAR + EmptyDrops). Includes intronic reads by default (v7+), builds its own biotype-filtered reference (drops rRNA / Mt-rRNA / IG / TR / many pseudogenes). Best when you want results directly comparable to 10x Cloud / Loupe or to prior Cell Ranger runs. **Default / recommended.**
2. **`star`** — STARsolo. Fast, exonic `Gene` counts by default (no intron inclusion), full Ensembl gene universe. Good open-source alternative; expect ~fewer genes/UMIs per cell than Cell Ranger and a higher apparent mito-%.
3. **`simpleaf`** — alevin-fry (simpleaf). USA mode counts spliced + unspliced + ambiguous → velocity-ready and tracks Cell Ranger's per-cell counts closely. Detects the most genes overall (retains rRNA / IG). *(On scrnaseq < 3.0.0 this aligner is named `alevin`.)*
4. **`kallisto`** — kallisto | bustools. Fast pseudoalignment, spliced counts (tracks STARsolo). Note the v4 whitelist caveat (handled automatically in Step 11b).

Store as `{ALIGNER}`. Resolve `{PROTOCOL}` per the table in Step 5.

---

## Step 7 — CellBender ambient-RNA removal

CellBender (`remove-background`) removes ambient / empty-droplet RNA. It is GPU-oriented; **this cluster has no GPUs**, so it runs on CPU (~150 epochs, ≈ 9–14 h per sample). The injected config gives it a 24 h budget and retries on OOM.

Ask — numbered options (default = Yes):
1. **Yes — run CellBender** (default). Omit `--skip_cellbender`; ensure the `nextflow.config` includes the `CELLBENDER_REMOVEBACKGROUND` resource block (Step 8b).
2. **No — skip it.** Add `--skip_cellbender` to the submission script.

Inform the user of the ~9–14 h/sample CPU cost when they choose Yes, so the head-job walltime (2 days) makes sense.

---

## Step 8 — Organism and reference files

Ask:
1. "What organism is this data from? (e.g. mouse, human)"
2. "What is the base directory where genome files should be stored?"

Unlike the bulk skill, **you do not pre-build aligner indexes** — nf-core/scrnaseq builds the aligner-specific reference internally (Cell Ranger `mkref`, STAR index, salmon index, or kb ref) from the FASTA + GTF, and `--save_reference` caches it under the output directory. You only need the **primary-assembly FASTA** and the **GTF** present.

### Folder structure

```
{genome_base}/{organism}/{assembly}_ens{version}/
├── {SPECIES}.dna.primary_assembly.fa      ← primary assembly FASTA
├── {SPECIES}.dna.primary_assembly.fa.gz   ← compressed original (kept)
├── {SPECIES}.{version}.gtf                 ← annotation GTF
└── {SPECIES}.{version}.gtf.gz
```

### Mouse (Mus musculus)
- Assembly GRCm39 | Directory `{genome_base}/mouse/mm39_ens{version}/`
- FASTA `Mus_musculus.GRCm39.dna.primary_assembly.fa`
- GTF `Mus_musculus.GRCm39.{version}.gtf`

### Human (Homo sapiens)
- Assembly GRCh38 | Directory `{genome_base}/human/hg38_ens{version}/`
- FASTA `Homo_sapiens.GRCh38.dna.primary_assembly.fa`
- GTF `Homo_sapiens.GRCh38.{version}.gtf`

### Procedure
1. Fetch latest Ensembl release: `https://ftp.ensembl.org/pub/current_README` (via `WebFetch`).
2. `ls {genome_base}/{organism}/` — find the highest `ens{N}` present.
3. If the FASTA and GTF already exist, report their paths and proceed.
4. If missing, generate a small download script (`download_reference_{assembly}_ens{version}.sh`) using `wget -c` and `gunzip -c` (never `gunzip -k` — unavailable on CentOS 7), create the directory with `mkdir -p`, and tell the user to run it before submitting:

```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH -t 2:00:00
#SBATCH -p bcc
#SBATCH -o download_reference_{assembly}_ens{version}.log

GENOME_DIR="{genome_base}/{organism}/{assembly}_ens{version}"
mkdir -p "${GENOME_DIR}"; cd "${GENOME_DIR}"

wget -c "{ENSEMBL_FASTA_URL}" -O {FASTA_FILENAME}.gz
wget -c "{ENSEMBL_GTF_URL}"   -O {GTF_FILENAME}.gz
gunzip -c {FASTA_FILENAME}.gz > {FASTA_FILENAME}
gunzip -c {GTF_FILENAME}.gz   > {GTF_FILENAME}
echo "Reference download complete."
```

**Reference reuse (optional):** if a previous run saved a reference (`--save_reference` writes it under `<outdir>/reference/`), it can be reused by pointing the aligner-specific flag (`--cellranger_index` / `--star_index` / `--salmon_index` / `--kallisto_index`) at the saved directory to skip the rebuild. Default behaviour is to build + save.

For organisms other than mouse/human, ask the user to supply FASTA and GTF paths manually and skip the Ensembl lookup.

---

## Step 8b — nextflow.config check and generation

Check whether `nextflow.config` exists in the current working directory (`ls nextflow.config`).

- **If it exists**: do not overwrite it. If CellBender is enabled (Step 7), verify it contains a `CELLBENDER_REMOVEBACKGROUND` block; if not, tell the user and offer to add it.
- **If it does not exist**: generate one using the template below and inform the user.

```nextflow
// nextflow.config

profiles {
    slurm {
        process {
            executor = 'slurm'
            queue    = 'bcc'
            cpus     = 2
            memory   = '8 GB'
            time     = '4h'

            withLabel: process_low         { cpus = 2;  memory = '12 GB'; time = '4h'  }
            withLabel: process_medium      { cpus = 6;  memory = '36 GB'; time = '8h'  }
            withLabel: process_high        { cpus = 12; memory = '72 GB'; time = '16h' }
            withLabel: process_high_memory {           memory = '120 GB'                }

            // Heavy aligner / counting steps
            withName: '.*CELLRANGER_COUNT'        { cpus = 12; memory = '90 GB';  time = '24h' }
            withName: '.*STAR_ALIGN'              { cpus = 12; memory = '72 GB';  time = '16h' }
            withName: '.*SIMPLEAF_QUANT'          { cpus = 12; memory = '72 GB';  time = '16h' }
            withName: '.*KALLISTOBUSTOOLS_COUNT'  { cpus = 12; memory = '72 GB';  time = '16h' }

            // CellBender on CPU (no GPU on this cluster): long walltime + OOM retries
            withName: '.*CELLBENDER_REMOVEBACKGROUND' {
                cpus          = 4
                memory        = { 32.GB * task.attempt }
                time          = '24.h'
                errorStrategy = { task.exitStatus in [104,134,137,139,140,143,247] ? 'retry' : 'finish' }
                maxRetries    = 3
            }
        }

        executor {
            queueSize       = 20
            submitRateLimit = '10/1min'
            pollInterval    = '30s'
        }
    }

    singularity {
        singularity {
            enabled    = true
            autoMounts = true
            cacheDir   = '{HOME}/.singularity/cache'
        }
    }
}

params {
    max_cpus   = 16
    max_memory = '128 GB'
    max_time   = '48h'
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

Substitute `{HOME}` with the user's home (`echo $HOME`). If CellBender is skipped (Step 7 → No), you may leave the `CELLBENDER_REMOVEBACKGROUND` block in place (harmless) or remove it.

---

## Step 9 — Expected cells and other options

**a) Expected cells** — ask: "Do you know the approximate number of cells loaded/recovered per sample? This improves cell calling." Numbered options:
1. No — omit (aligners fall back to their knee/EmptyDrops defaults)
2. Yes, same for all samples — add an `expected_cells` column to the samplesheet with that value for every row
3. Yes, per-sample — collect values and fill the column per sample

**b) MultiQC title** — ask: "What title for the MultiQC report?" Numbered options:
1. `{SEQ_DATE}_{WD_NAME}`
2. `{TODAY_YYMMDD}_{WD_NAME}`
3. Custom

---

## Step 10 — Output directory

Ask: "What should the output directory be named?" Numbered options:
1. `{SEQ_DATE}_{WD_NAME}` — sequencing date (YYMMDD)
2. `{TODAY_YYMMDD}_{WD_NAME}` — today's date (YYMMDD)
3. Custom name

Store as `{OUTDIR}`.

---

## Step 11 — Generate the submission script

Write `nf-core_scrnaseq_{version}_{aligner}.sh` in the current working directory. The head job is intentionally **tiny** — with `executor = slurm` the head process only orchestrates and submits child jobs, so it needs almost no resources but must stay alive for the whole run (hence the 2-day walltime).

```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH -t 2-00:00:00
#SBATCH -p bcc
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={USER_EMAIL}
#SBATCH -o nf-core_scrnaseq_{version}_{aligner}.log

module add miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit
conda activate {CONDA_ENV}
module add singularity/3.10.4

nextflow run nf-core/scrnaseq -r {VERSION} \
-profile slurm,singularity \
-c nextflow.config {KALLISTO_CONFIG_FLAG}\
--input {SAMPLESHEET_CSV} \
--outdir {OUTDIR} \
--aligner {ALIGNER} \
--protocol {PROTOCOL} \
--fasta {FASTA_PATH} \
--gtf {GTF_PATH} \
--save_reference \
{CELLBENDER_LINE}--multiqc_title {MULTIQC_TITLE} \
-resume
```

**Substitutions:**
- `{KALLISTO_CONFIG_FLAG}`: `-c kallisto_v4.config ` **only** when `{ALIGNER}` is kallisto AND chemistry is v4 (see Step 11b). Otherwise empty.
- `{PROTOCOL}`: `auto` for cellranger; detected value for star/simpleaf; `10XV3` for kallisto-on-v4 (with the whitelist config); detected value otherwise.
- `{CELLBENDER_LINE}`: empty when CellBender is enabled; `--skip_cellbender \` (with trailing continuation) when skipped.
- Drop any flag that matches a pipeline default to keep the script clean.

After writing, show the full file and instruct:
```
Script written: nf-core_scrnaseq_{version}_{aligner}.sh
To submit:  sbatch nf-core_scrnaseq_{version}_{aligner}.sh
```
If the reference still needs downloading (Step 8):
```
NOTE: Reference FASTA/GTF not present yet. Download it first:
  sbatch download_reference_{assembly}_ens{version}.sh
```

---

## Step 11b — kallisto v4 whitelist (only when aligner = kallisto AND chemistry = v4)

kb-python 0.28.2 (bundled with the pipeline) has no built-in v4 onlist, so kallisto must be told the v4 whitelist explicitly. Two artifacts are needed.

**1. Extract the v4 barcode whitelist** (it ships gzip'd inside the pipeline repo, which Nextflow caches under `~/.nextflow/assets/nf-core/scrnaseq`):

```bash
mkdir -p {CWD}/refs
git -C ~/.nextflow/assets/nf-core/scrnaseq show \
  {VERSION}:assets/whitelist/10x_V4_barcode_whitelist.txt.gz \
  | gunzip -c > {CWD}/refs/10x_V4_whitelist.txt
wc -l {CWD}/refs/10x_V4_whitelist.txt   # expect 7,372,800 barcodes
```
(If the pipeline hasn't been pulled yet, run `nextflow pull nf-core/scrnaseq -r {VERSION}` first.) This same file feeds the v3/v4 disambiguation in Step 5.

**2. Write `kallisto_v4.config`** in the working directory, injecting the whitelist via `kb count -w`:

```nextflow
process {
    withName: '.*KALLISTOBUSTOOLS_COUNT' {
        ext.args = '-w {CWD}/refs/10x_V4_whitelist.txt'
    }
}
```

The submission script then adds `-c kallisto_v4.config` and runs kallisto with `--protocol 10XV3` (identical geometry to v4).

---

## Notes for the assistant

- **`gh` CLI is not available on this HPC cluster.** Never attempt `gh api`. Use `WebFetch` for all GitHub API calls.
- **Always present questions with a finite set of answers as a numbered list of options.** Only use open-ended text questions when no reasonable discrete set exists (email, conda env name, custom paths/values).
- Never run heavy computation on the login node. All pipeline work goes through `sbatch`; the tiny head job orchestrates SLURM child jobs.
- Raw FASTQ files are read-only; never modify them.
- **`fastq_1` = R1 = barcode+UMI; `fastq_2` = R2 = cDNA. Never swap them.** Verify by read length at least once. Exclude `I1`/`I2` index reads.
- **Same sample name across multiple rows = lane concatenation, which is correct for 10x.** Do not warn about it (opposite of the bulk skill).
- **10x 3' v4 (GEM-X) needs scrnaseq ≥ 3.0.0 and `--protocol 10XV4`.** v3 and v4 share geometry (16 bp CB + 12 bp UMI); separate them by whitelist membership, not read length.
- **cellranger → `--protocol auto`.** It auto-detects chemistry and is the most robust; do not force a protocol for it.
- **kallisto + v4**: run as `--protocol 10XV3` + inject the v4 whitelist through `kallisto_v4.config` (kb-python can't parse `10XV4`).
- **The `simpleaf` aligner is named `alevin` on scrnaseq < 3.0.0.**
- **CellBender has no GPU here** → CPU run ~9–14 h/sample; the `CELLBENDER_REMOVEBACKGROUND` block gives it 24 h and OOM retries. Skip it with `--skip_cellbender` if not needed.
- **`gunzip -k` is unavailable on CentOS 7** — always use `gunzip -c file.gz > file`.
- **Running several aligners to compare?** Launch each from its **own directory** (e.g. `runs/{aligner}/`) so the concurrent Nextflow runs don't collide on the shared `.nextflow` session lock. Point `--outdir` at a shared `results/` and re-use `--save_reference` output where possible.
- Downstream QC thresholds are **not** portable across aligners: the biotype/intron differences drive large gaps in per-cell gene/UMI counts and mito-% (e.g. ~3.5% mito in Cell Ranger vs 11–17% in STARsolo/kallisto on the same sample). Re-tune `nFeature` / `nCount` / `percent.mt` per aligner before Seurat.
```
