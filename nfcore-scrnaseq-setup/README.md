# `/nfcore-scrnaseq-setup` — nf-core/scrnaseq Pipeline Setup Skill

An interactive Claude Code skill that walks you through setting up and submitting an [nf-core/scrnaseq](https://nf-co.re/scrnaseq) single-cell RNA-seq pipeline on an HPC cluster running SLURM and Singularity. It is the single-cell companion to [`/nfcore-rnaseq-setup`](../nfcore-rnaseq-setup) and follows the same wizard structure.

It sets up **one aligner per run**. To compare aligners, run the skill once per aligner from a separate launch directory (see *Cluster-specific notes*).

---

## Installation

```bash
cp nfcore-scrnaseq-setup.md ~/.claude/commands/
```

Then invoke it in Claude Code:

```
/nfcore-scrnaseq-setup
```

---

## Prerequisites

| Requirement | Notes |
|-------------|-------|
| SLURM scheduler | Tested on `bcc` queue |
| Singularity ≥ 3.10 | Loaded via `module add singularity/3.10.4` |
| Nextflow | Available in a conda environment (e.g. `nf-env`) |
| nf-core/scrnaseq ≥ 3.0.0 | **Required for 10x 3' v4 (GEM-X)** and for the `simpleaf` aligner name |
| Internet access | Required from login node for GitHub API + Ensembl downloads |

---

## What the skill does

A guided wizard that collects your settings one step at a time and generates ready-to-submit SLURM scripts.

### Step 0 — Working directory
Detects the current directory. All output files (samplesheet, config, scripts) are written here.

### Step 1 — Email
SLURM `--mail-user` address for job completion/failure notifications.

### Step 2 — Conda environment
Name of the conda environment that has Nextflow installed.

### Step 3 — nf-core/scrnaseq version
Fetches the latest release tag automatically from the GitHub API. Enforces ≥ 3.0.0 when v4 chemistry is detected, and maps the `simpleaf`/`alevin` aligner name to the version.

### Step 4 — FASTQ files and samplesheet
- Scans the working directory for `.fastq.gz` files.
- Assigns **10x read roles by measured read length**: R1 (16 bp CB + 10/12 bp UMI, ~26–28 bp) → `fastq_1`; R2 (cDNA, ~90 bp) → `fastq_2`. Index reads (`I1`/`I2`) are excluded.
- Detects the **sequencing date** from leading `YYMMDD` filename prefixes.
- Builds a `sample,fastq_1,fastq_2` samplesheet with sanitised names (dashes → underscores).
- **Treats same-name rows as lane concatenation** — the correct behaviour for multi-lane 10x libraries (opposite of the bulk skill's duplicate warning).

### Step 5 — Chemistry / protocol detection
Auto-detects 10x chemistry, then asks you to confirm:
- Barcode-read length → v2 (26 bp) vs v3/v4 (28 bp).
- **v3 vs v4 disambiguation** by barcode-whitelist membership (v3 `3M-february-2018` ~6.8M vs v4/GEM-X ~7.4M).
- Resolves the right `--protocol` per aligner (`auto` for cellranger; detected value for star/simpleaf; `10XV3` + whitelist for kallisto-on-v4).

### Step 6 — Aligner
- `cellranger` (default) — reference counter, intron-inclusive, biotype-filtered reference
- `star` — STARsolo, exonic counts, full gene universe
- `simpleaf` — alevin-fry, USA mode (velocity-ready), tracks Cell Ranger closely *(named `alevin` on < 3.0.0)*
- `kallisto` — kallisto|bustools, spliced counts

Each option carries a one-line "when to use it" note.

### Step 7 — CellBender
Ambient-RNA removal, **enabled by default**. Because this cluster has no GPUs it runs on CPU (~9–14 h/sample); the generated config gives it a 24 h budget with OOM retries. Choose *No* to add `--skip_cellbender`.

### Step 8 — Organism and reference
Supports **mouse** (GRCm39) and **human** (GRCh38) with automatic Ensembl version detection. No pre-built index is required — the pipeline builds the aligner reference internally from FASTA + GTF and caches it with `--save_reference`. If the FASTA/GTF are missing, the skill writes a small `download_reference_*.sh` script.

### Step 8b — nextflow.config
Generates a SLURM + Singularity config with label-based resources, per-aligner `withName` blocks (Cell Ranger / STARsolo / simpleaf / kallisto), and the CPU-tuned `CELLBENDER_REMOVEBACKGROUND` block (24 h walltime, memory scaling on retry).

### Step 9 — Expected cells & MultiQC title
Optional `expected_cells` column (improves cell calling) and the MultiQC report title.

### Step 10 — Output directory
Names the output directory (sequencing date + project, today's date + project, or custom).

### Step 11 — Submission script
Writes `nf-core_scrnaseq_{version}_{aligner}.sh` — a tiny orchestration head job (executor = slurm submits the real child jobs) with all chosen flags.

### Step 11b — kallisto v4 whitelist *(only for kallisto + v4)*
Extracts the v4 barcode whitelist from the cached pipeline repo and writes `kallisto_v4.config` to inject it via `kb count -w`, since kb-python 0.28.2 can't parse `10XV4`.

---

## Output files

| File | Description |
|------|-------------|
| `{date}_{project}_samplesheet.csv` | Input samplesheet for nf-core/scrnaseq |
| `nextflow.config` | SLURM + Singularity resource profiles (incl. CellBender block) |
| `nf-core_scrnaseq_{version}_{aligner}.sh` | Pipeline SLURM submission script |
| `download_reference_{assembly}_ens{version}.sh` | Reference download script (if FASTA/GTF were missing) |
| `kallisto_v4.config` + `refs/10x_V4_whitelist.txt` | Only for the kallisto + v4 case |

---

## Typical workflow

```bash
# 1. (If reference is missing) download FASTA + GTF first
sbatch download_reference_GRCm39_ens115.sh

# 2. Submit the pipeline
sbatch nf-core_scrnaseq_4.1.0_cellranger.sh
```

---

## Cluster-specific notes

- **Queue**: `bcc`
- **Singularity module**: `singularity/3.10.4`
- **Conda init**: `source /home/software/conda/miniconda3/bin/condainit`
- **`gh` CLI is not available** — the skill uses the GitHub REST API via HTTP instead
- **`gunzip -k` is not available** on CentOS 7 — use `gunzip -c file.gz > file`
- **No GPUs** — CellBender runs on CPU (~9–14 h/sample)
- **Comparing aligners?** Run the skill once per aligner and launch each from its **own directory** (e.g. `runs/{aligner}/`) so concurrent Nextflow runs don't collide on the shared `.nextflow` session lock.

---

## Known limitations

- Organisms other than mouse and human require manual FASTA/GTF paths.
- **10x 3' v4 (GEM-X) requires nf-core/scrnaseq ≥ 3.0.0.** Earlier versions mis-assign v4 barcodes.
- kallisto on v4 runs as `10XV3` geometry with an injected v4 whitelist (kb-python 0.28.2 has no native `10XV4`).
- Downstream QC thresholds (`nFeature` / `nCount` / `percent.mt`) are **not** portable across aligners — the intron-inclusion and reference-biotype differences shift per-cell counts and mito-% substantially. Re-tune per aligner before Seurat.
```
