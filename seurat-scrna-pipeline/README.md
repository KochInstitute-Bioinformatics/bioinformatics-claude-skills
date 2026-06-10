# `/seurat-scrna-pipeline` — scRNA-seq Downstream Analysis Skill

An interactive Claude Code skill that generates a complete single-cell RNA-seq downstream analysis pipeline (Seurat v5 / Harmony) as a series of R Markdown reports and SLURM submission scripts, designed to run inside a Singularity container on an HPC cluster.

---

## Installation

```bash
cp seurat-scrna-pipeline.md ~/.claude/commands/
```

Then invoke it in Claude Code:

```
/seurat-scrna-pipeline
```

---

## Prerequisites

| Requirement | Notes |
|-------------|-------|
| SLURM scheduler | Tested on `bcc` queue |
| Singularity ≥ 3.5 | Loaded via `module load singularity/3.5.0` |
| Container image | `yannvrb56/r4_6_0_python3_singlecell_visium:latest` (SIF on shared storage) |
| Starting data | Per-sample Seurat RDS, CellRanger output, or a merged/integrated RDS |
| Internet access | Required from compute node for scType, Azimuth references, and msigdbr |

The container bundles Seurat v5, Harmony, scDblFinder, SingleR, scType, DESeq2, clusterProfiler, enrichplot, msigdbr (26.x), LIANA, limma, and pandoc.

---

## What the skill does

The skill is a guided wizard that detects your input data, collects analysis settings, and writes one R Markdown report per processing stage plus a matching SLURM script for each. Stages:

### Rmd 01 — Data Import & QC
Per-sample processing: Ensembl→symbol gene name conversion (if needed), QC metrics, adaptive MAD-based per-sample filtering (with a mito floor for CellBender data), and scDblFinder doublet scoring. Saves one cleaned RDS per sample.

### Rmd 02 — Merge & QC Review
Merges all samples, manual cell-cycle scoring, rough clustering **for QC inspection only** (no integration), and per-cell filtering (doublet removal — never whole-cluster exclusion). Saves the filtered object.

### Rmd 03 — Integration
SCTransform on split layers → PCA → **Harmony** integration → clustering → UMAP → FindAllMarkers. Generates the **ACT** cluster input file for the CellVision online annotation tool. Includes checkpoints for the expensive SCT+PCA and FindAllMarkers steps.

### Rmd 04 — Cell Type Assignment
Consensus annotation across **SingleR** (one or two references), **scType** (organ-level tissue context), **Azimuth**, and **ACT**. Produces a per-cluster comparison table and consensus vote, plus an auto-generated `RenameIdents` scaffold for review.

### Rmd 05 — Final Annotation & Differential Abundance
Applies the user's reviewed cell type labels, annotated UMAPs, marker dot plots, proportion tables, and **differential abundance testing** (propeller, implemented via limma).

### Rmd 06 — Differential Gene Expression *(optional)*
Pseudobulk **DESeq2** per cell type (recommended for ≥2 replicates/condition) or Wilcoxon `FindMarkers`. Includes targeted gene-list dot plots for curated pathways.

### Rmd 07 — Gene Set Enrichment Analysis *(optional)*
**GSEA** (clusterProfiler / fgsea) on DEG results, ranked by DESeq2 Wald statistic. **MSigDB collections are chosen per project based on biological context** (Hallmark, canonical pathways, GO:BP, immunologic signatures, cell-type signatures, oncogenic, TF targets).

### Rmd 08 — Ligand–Receptor Analysis *(optional)*
Intercellular communication inference with **LIANA** (MouseConsensus / Consensus resource) on the full integrated object, filtered to sender/receiver cell types of interest.

### Summary HTML report
A self-contained Bootstrap page summarising study design, QC, integration, and annotation, with links to each individual report and the final Seurat object.

---

## Output files

| File | Description |
|------|-------------|
| `{date}_{project}_01_dataImport.Rmd` … `_08_LR_LIANA.Rmd` | One R Markdown per analysis stage |
| `run_01_dataImport.sh` … `run_08_LR_LIANA.sh` | Matching SLURM submission scripts |
| `qc_rds/` | Per-sample QC-filtered Seurat objects |
| `results/{date}_{project}/` | Rendered HTML reports, RDS objects, xlsx tables, figures |
| `{project}_summary_report.html` | Standalone summary with links to all outputs |

---

## Typical workflow

```bash
sbatch -p bcc run_01_dataImport.sh
sbatch -p bcc run_02_qc.sh                 # review QC plots, verify cell counts
sbatch -p bcc run_03_integration.sh        # then submit ACT input to thecellvision.org
sbatch -p bcc run_04_celltypeAssignment.sh # review annotation_consensus_per_cluster.xlsx
# fill in consensus cell type labels, then:
sbatch -p bcc run_05_finalAnnotation.sh
# optional downstream analyses:
sbatch -p bcc run_06_DEG.sh
sbatch -p bcc run_07_<celltype>_GSEA.sh
sbatch -p bcc run_08_LR_LIANA.sh
```

---

## Cluster-specific notes

- **Queue**: `bcc` — always submit with `sbatch -p bcc script.sh`; the `#SBATCH -p bcc` header alone silently fails on this cluster.
- **Singularity module**: `module load singularity/3.5.0`.
- **Never run R or heavy computation on the login node** — use `sbatch` (or `srun`/`salloc` for debugging).
- **`gh` CLI is not available**.
- Packages are baked into the container image — never `install.packages()` at runtime.

---

## Key gotchas the skill handles

These are hard-won lessons baked into the skill so you don't rediscover them:

- **Ensembl IDs as gene names** silently break `percent.mt`, SingleR, scType, ACT, and Azimuth — converted up front.
- **Seurat v5 API**: `slot=`→`layer=`; `CellCycleScoring`, `as.SingleCellExperiment()`, and `DietSeurat()` call defunct internals — manual equivalents used.
- **Harmony ≥ 1.0**: use `RunHarmony()` directly, not `IntegrateLayers(HarmonyIntegration)`.
- **Azimuth corrupts gene names** (mouse→human, in place) — downstream gene-expression code loads the pre-Azimuth object.
- **propeller**: `speckle` is not installed — a manual limma reimplementation (arcsine-sqrt + moderated t-test) is used; run per sorted cell population.
- **msigdbr 26.x**: uses `collection=`/`subcollection=` (not the deprecated `category=`); for mouse, ortholog-maps from the human DB (`db_species = "HS"`).
- **LIANA v0.1.14** needs a Seurat v3-style assay — the v5 object is rebuilt before running.

---

## Known limitations

- Azimuth `pbmcref` (human PBMC) annotates immune clusters well but cannot label epithelial/stromal cells — rely on scType + markers for those.
- Consensus annotation labels are provisional; manual validation against canonical markers is required before publication.
- The mouse-native MSigDB (`db_species = "MM"`) is not bundled in the current image — GSEA uses human→mouse ortholog mapping.
- Word/PowerPoint deliverables are produced via pandoc (`officer` and `python-pptx` are not in the image).
