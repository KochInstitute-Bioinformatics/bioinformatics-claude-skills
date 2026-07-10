# `/bulk-rnaseq-pipeline` — Bulk RNA-seq Downstream Analysis Skill

An interactive Claude Code skill that generates a complete bulk RNA-seq downstream analysis pipeline (tximport → DESeq2 / edgeR → GSEA) as a series of R Markdown reports and SLURM submission scripts, starting from **nf-core/rnaseq `star_salmon` output** and running inside a Singularity container on an HPC cluster. It is the bulk companion to [`/seurat-scrna-pipeline`](../seurat-scrna-pipeline) and picks up where [`/nfcore-rnaseq-setup`](../nfcore-rnaseq-setup) leaves off.

---

## Installation

```bash
cp bulk-rnaseq-pipeline.md ~/.claude/commands/
```

Then invoke it in Claude Code:

```
/bulk-rnaseq-pipeline
```

---

## Prerequisites

| Requirement | Notes |
|-------------|-------|
| SLURM scheduler | Tested on `bcc` queue |
| Singularity ≥ 3.5 | Loaded via `module load singularity/3.5.0` |
| Container image | `yannvrb56/bulkrnaseq:latest` (SIF at `/net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif`) |
| Starting data | nf-core/rnaseq `star_salmon` output — per-sample `quant.sf` + multiqc general stats |
| Internet access | Required from compute node for `msigdbr` 26.x (fetches gene-set data at runtime) if GSEA is used |

The container bundles tximport, DESeq2, edgeR, apeglm, fgsea, msigdbr (26.x), GseaVis, openxlsx, tidyverse, matrixStats, ggrepel, showtext, and pandoc.

---

## What the skill does

The skill is a guided wizard that detects your nf-core/rnaseq output, **confirms the library mode**, collects analysis settings, and writes one R Markdown report per stage plus a matching SLURM script. Stages:

### Rmd 01 — Import, Annotation & QC
`tximport` of the per-sample `quant.sf` files to gene level, annotation from nf-core's `salmon.merged.tx2gene.tsv` + a streaming-awk `gene_biotype` map, a curated multiqc QC join, expressed-protein-coding filtering (MPG + `LowExpLowVar`), QC boxplots, and exploratory PCA coloured by both biology and each QC metric. Saves a checkpoint RDS.

### Rmd 02 — Differential Expression
Builds `dds` (mode-branched — see below), runs DESeq2, per-contrast `lfcShrink(apeglm)` with the Wald `stat` recovered for ranking, volcano plots (top-N labelled by `abs(stat)`), a DEG summary table, a multi-sheet supplemental workbook, and GEO-ready processed-data matrices.

### Rmd 03 — GSEA *(optional)*
In-R `fgsea` ranked by the DESeq2 Wald statistic, over user-chosen **MSigDB collections** (per-project biology choice — never hardcoded Hallmark), NES plots, leading-edge tables, and optional `GseaVis::gseaNb` enrichment curves.

### Summary / client report
A self-contained summary HTML linking every stage report and data deliverable, plus an optional **consolidated client report** Rmd that re-loads saved outputs and embeds all figures/tables as a single portable file (base64, zero external `src`).

---

## The central distinction: paired-end vs 3′ DGE

The single most error-prone part of bulk RNA-seq downstream analysis — the skill detects it from multiqc, **asks you to confirm**, and branches everything on it:

| | Paired-end (standard mRNA) | 3′ DGE (3' tag) |
|---|---|---|
| Length normalisation | **Yes** — transcript length matters | **No** — 3' tag count is length-independent |
| Expression units | **TPM** (`txi$abundance`) → `l2tpm` | **CPM** (edgeR `cpm()`) → `l2cpm` |
| DESeq2 construction | `DESeqDataSetFromTximport(txi, …)` (uses length offset) | `DESeqDataSetFromMatrix(round(txi$counts), …)` (no offset) |
| edgeR used? | No | Yes (`DGEList` → `cpm`) |

Getting the mode wrong silently corrupts every downstream fold change, so it is confirmed before any Rmd is generated.

---

## Output files

| File | Description |
|------|-------------|
| `{date}_{project}_01_import_qc.Rmd` … `_03_gsea.Rmd` | One R Markdown per analysis stage |
| `run_01_import_qc.sh` … `run_03_gsea.sh` | Matching SLURM submission scripts |
| `ref/tx2gene.tsv`, `ref/gene2typesym.tsv` | Staged annotation (transcript→gene, gene→biotype/symbol) |
| `results/{date}_{project}/` | Rendered HTML reports, checkpoint RDS, xlsx tables, GEO matrices, figures |
| `*_supplemental.xlsx`, `*_DEG_summary.xlsx`, `*_GSEA_fgsea.xlsx` | Per-gene table, DEG counts, GSEA leading-edge sheets |
| `{project}_summary_report.html` | Standalone summary with links to all outputs |

---

## Typical workflow

```bash
sbatch -p bcc run_01_import_qc.sh
# Review 01_import_qc.html — boxplots, PCA (biology vs QC overlays), sample grouping
sbatch -p bcc run_02_deg.sh
# Review volcano plots and *_supplemental.xlsx
sbatch -p bcc run_03_gsea.sh          # if requested
```

---

## Cluster-specific notes

- **Queue**: `bcc` — always submit with `sbatch -p bcc script.sh`; the `#SBATCH -p bcc` header alone silently fails on this cluster.
- **Singularity module**: `module load singularity/3.5.0`.
- **Never run R or Singularity on the login node** — use `sbatch` (or `srun`/`salloc` for debugging).
- **`gh` CLI is not available**.
- Packages are baked into the `yannvrb56/bulkrnaseq` image — never `install.packages()` at runtime; rebuild the image to add packages.
- Raw `quant.sf` and nf-core outputs are read-only; all outputs go to `results/{date}_{project}/`.

---

## Key gotchas the skill handles

Hard-won lessons baked into the skill so you don't rediscover them:

- **Library mode** (paired-end vs 3′ DGE) drives both the expression units and the DESeq2 constructor — confirmed up front.
- **`quant.sf` counts are not integers** — always `round(txi$counts, 0)` for DESeq2 / GEO matrices.
- **Prefer nf-core's `salmon.merged.tx2gene.tsv`** for `tx2gene` (exact-matched to the quantified transcripts, no version mismatch).
- **Filtered nf-core GTFs have no top-level `gene` lines** — `gene_biotype` is extracted from `transcript` lines with streaming awk, never `rtracklayer::import` (multi-GB → OOM risk).
- **QC join** drops per-read ` Read 1`/` Read 2` sub-rows and asserts non-NA to catch sample-name mismatches loudly.
- **apeglm drops the `stat` column** — the Wald statistic is recovered from `dds@rowRanges` and used as the GSEA ranking metric.
- **Volcano labelling** ranks by `abs(stat)` (top ~20), not a LFC/padj threshold, to avoid an unreadable wall of labels on large DEG sets.
- **msigdbr 26.x** uses `collection=`/`subcollection=` and fetches data at runtime (needs internet); for mouse it ortholog-maps from the human DB (`db_species = "HS"`).
- **Namespace masking** — load DESeq2/Bioc before tidyverse so dplyr verbs win; `cache = FALSE` in every Rmd (chunks mutate a shared `assembleDat`).

---

## Known limitations

- Starts specifically from nf-core/rnaseq `star_salmon` output (`quant.sf`); other quantifiers require adapting the import step.
- If no GTF is reachable, the `gene_biotype` protein-coding filter is skipped (with a warning).
- GSEA on offline-compute clusters requires pre-seeding the `msigdbdf` cache (or baking it into the image) before submitting `run_03`.
- The mouse-native MSigDB is not bundled — GSEA uses human→mouse ortholog mapping.
