# `/nfcore-rnaseq-setup` — nf-core/rnaseq Pipeline Setup Skill

An interactive Claude Code skill that walks you through setting up and submitting an [nf-core/rnaseq](https://nf-co.re/rnaseq) bulk RNA-seq pipeline on an HPC cluster running SLURM and Singularity.

---

## Installation

```bash
cp nfcore-rnaseq-setup.md ~/.claude/commands/
```

Then invoke it in Claude Code:

```
/nfcore-rnaseq-setup
```

---

## Prerequisites

| Requirement | Notes |
|-------------|-------|
| SLURM scheduler | Tested on bcc queue |
| Singularity ≥ 3.10 | Loaded via `module add singularity/3.10.4` |
| Nextflow | Available in a conda environment (e.g. `nf-env`) |
| Internet access | Required from login node for GitHub API + Ensembl downloads |

---

## What the skill does

The skill is a guided wizard that collects your settings one step at a time and generates two ready-to-submit SLURM scripts. Here is what each step covers:

### Step 0 — Working directory
Detects the current directory. All output files (samplesheet, config, scripts) are written here.

### Step 1 — Email
SLURM `--mail-user` address for job completion/failure notifications.

### Step 2 — Conda environment
Name of the conda environment that has Nextflow installed.

### Step 3 — nf-core/rnaseq version
Fetches the latest release tag automatically from the GitHub API. No manual lookup needed.

### Step 4 — FASTQ files and samplesheet
- Scans the working directory for `.fastq.gz` files automatically.
- Auto-detects **paired-end** vs **single-end** layout from filename patterns (`_R1_`/`_R2_`, `_1_sequence`/`_2_sequence`, etc.).
- Detects the **sequencing date** from leading `YYMMDD` filename prefixes.
- Builds the samplesheet CSV with sanitised sample names (dashes replaced with underscores to prevent R parsing errors).
- Presents the full samplesheet for review and offers custom naming if needed.
- Uses relative paths when FASTQ files are inside the working directory, absolute paths otherwise.

### Step 5 — Trimming
Choose between:
- Already trimmed (skip trimming)
- `fastp` (recommended)
- `trimgalore`

Optionally set a minimum reads-after-trimming threshold.

### Step 5b — rRNA removal
Optionally add `--remove_ribo_rna` for libraries that are not ribo-depleted.

### Step 6 — 3' DGE / UMI detection *(single-end only)*
For single-end data, the skill checks three signals:
1. Filename keywords (quantseq, lexogen, brbseq, marsseq, …)
2. Read length (≤20 bp = strong evidence, 20–50 bp = mild)
3. PolyT content at read starts

Prompts for confirmation, then configures `--noLengthCorrection` in `nextflow.config` and optional UMI flags (`--with_umi`, `--umitools_bc_pattern`, etc.) as appropriate.

### Step 7 — Organism and genome
Supports **mouse** (GRCm39) and **human** (GRCh38) with automatic Ensembl version detection.

Folder structure created under your genome base directory:
```
{genome_base}/{organism}/{assembly}_ens{version}/
├── *.fa / *.gtf / *.cdna.all.fa       ← downloaded from Ensembl
├── *.cdna.all.noversionnum.fa          ← version-stripped for Salmon
├── *.cdna.gtf_filtered.fa              ← GTF-filtered cDNA (used for Salmon index)
├── decoys.txt / gentrome.fa
└── index/
    ├── star/
    └── salmon/
```

If the genome index is missing or outdated, the skill generates an index-build script (Step 14b) and creates the required directories.

### Step 7b — nextflow.config
Checks for an existing `nextflow.config`. If none is found, generates one with SLURM resource profiles for STAR, Salmon, FastQC, and MultiQC.

### Step 8 — MultiQC title
Names the MultiQC report (sequencing date + project, today's date + project, or custom).

### Step 9 — Aligner
- `star_salmon` (default, recommended)
- `star_rsem`
- `hisat2`
- `bowtie2_salmon`

### Step 10 — Pseudo-aligner
- `salmon` (recommended alongside STAR)
- `kallisto`
- None

### Step 11 — StringTie
Optionally skip StringTie transcript assembly (recommended for standard differential expression).

### Step 12 — Additional options
- Save intermediate BAM files
- Skip all QC steps except MultiQC

### Step 13 — Output directory
Names the output directory (sequencing date + project, today's date + project, or custom).

### Step 14a — Pipeline submission script
Writes `nf-core_rnaseq_{version}.sh` — a ready-to-submit SLURM script with all chosen flags.

### Step 14b — Genome index build script *(only when indexes are missing)*
Writes `build_genome_index_{assembly}_ens{version}.sh`, which:
1. Downloads FASTA, GTF, and cDNA from Ensembl via `wget -c`
2. Decompresses with `gunzip -c` (compatible with CentOS 7 — no `-k` flag)
3. Builds a STAR index with the correct `sjdbOverhang` (= read length − 1, auto-detected)
4. Strips Ensembl version suffixes from cDNA transcript IDs
5. Filters cDNA to GTF-annotated transcripts only (prevents tximport failures)
6. Builds a Salmon decoy-aware index from `gentrome.fa`

---

## Output files

After running the skill you will have:

| File | Description |
|------|-------------|
| `{date}_{project}_samplesheet.csv` | Input samplesheet for nf-core/rnaseq |
| `nextflow.config` | SLURM + Singularity resource profiles |
| `nf-core_rnaseq_{version}.sh` | Pipeline SLURM submission script |
| `build_genome_index_{assembly}_ens{version}.sh` | Index build script (if genome was missing) |

---

## Typical workflow

```bash
# 1. (If genome indexes are missing) build them first — takes ~2–4 hours
sbatch build_genome_index_GRCh38_ens115.sh

# 2. Once indexes are ready, submit the pipeline
sbatch nf-core_rnaseq_3.26.0.sh
```

---

## Cluster-specific notes

- **Queue**: `bcc`
- **Singularity module**: `singularity/3.10.4`
- **Conda init**: `source /home/software/conda/miniconda3/bin/condainit`
- **`gh` CLI is not available** — the skill uses the GitHub REST API via HTTP instead
- **`gunzip -k` is not available** on CentOS 7 — use `gunzip -c file.gz > file`

---

## Known limitations

- Organisms other than mouse and human require manual FASTA/GTF paths
- `--umi_discard_read 1` is only valid for paired-end data (R1 = UMI-only read)
- Re-using a STAR index built with a different read length works fine — `sjdbOverhang` is baked in at index build time
