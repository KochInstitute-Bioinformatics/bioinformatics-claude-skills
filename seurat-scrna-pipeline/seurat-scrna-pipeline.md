# scRNAseq Downstream Analysis Skill (Seurat v5 / Harmony / HPC)

You are helping the user run single-cell RNA-seq downstream analysis on an HPC cluster using SLURM and Singularity. The container is `yannvrb56/r4_6_0_python3_singlecell_visium:latest` (SIF at `/net/bmc-lab3/data/bcc/shared/singularity_images/r4_6_0_python3_singlecell_visium_latest.sif`). Generate R Markdown files (one per processing stage) and a SLURM submission script for each. Walk through the steps below in order. Perform all automated steps (file scans, directory detection) silently before presenting results.

---

## Step 0 — Working directory and constants

Run `pwd` to record `{CWD}`. Derive:
- `{WD_NAME}` = basename of `{CWD}`
- `{TODAY_YYMMDD}` = today's date as `YYMMDD`
- `{RESULTS_DIR}` = `{CWD}/results/{TODAY_YYMMDD}_{WD_NAME}`

Inform the user:
```
Working directory : {CWD}
Results will go to: {RESULTS_DIR}
```

---

## Step 1 — Email address

Ask: "What email address should SLURM use for job notifications?"

---

## Step 2 — Verify Singularity image

```bash
ls /net/bmc-lab3/data/bcc/shared/singularity_images/r4_6_0_python3_singlecell_visium_latest.sif 2>/dev/null
```

- **Found**: proceed silently.
- **Not found**: stop and tell the user to pull it on a compute node first.

---

## Step 3 — Detect input data structure

**Auto-detect silently:**

```bash
# Pre-converted Seurat RDS per sample (most common starting point)
find {CWD} -maxdepth 8 -name "*filtered_matrix.seurat.rds" 2>/dev/null | head -20

# CellRanger multi per_sample_outs
find {CWD} -maxdepth 8 -name "sample_filtered_feature_bc_matrix" -type d 2>/dev/null | head -10

# Standard CellRanger count
find {CWD} -maxdepth 8 -name "filtered_feature_bc_matrix" -type d 2>/dev/null | head -10

# Features file for gene symbol mapping
find {CWD} -maxdepth 10 -name "features.tsv.gz" 2>/dev/null | head -3
```

Tell the user what was found. Ask: "Where would you like to start?"

1. **From per-sample Seurat RDS files** (CellRanger Multi output, already converted)
2. **From raw CellRanger output** (filtered_feature_bc_matrix directories)
3. **From a merged/QC-filtered RDS** — skip to integration
4. **From an integrated RDS** — skip to cell type annotation

Store as `{ENTRY_POINT}`.

**If starting from Seurat RDS (option 1):** Note the path to the `features.tsv.gz` file — it is needed for Ensembl ID → gene symbol conversion (Step 7).

---

## Step 4 — Sample discovery and metadata

**Auto-detect samples** based on entry point. Show the user a table of detected samples. Ask for:
1. The condition/group for each sample
2. Which cell populations were sorted (e.g. CD45+, EpCAM+, GFP+) — this affects annotation reference choice

Validate: conditions cannot have spaces. Warn if all samples are in one condition (no DEG possible).

Store as `{SAMPLE_META}`.

---

## Step 5 — Organism

Ask: "What organism?"
1. Mouse (*Mus musculus*) — MT prefix: `^mt-`
2. Human (*Homo sapiens*) — MT prefix: `^MT-`
3. Other — ask for MT prefix

Store as `{MT_PATTERN}` and `{ORGANISM}`.

---

## Step 6 — Clustering resolution

Ask: "What clustering resolution(s)?"
1. Default (0.4)
2. Single custom value
3. Multiple values (e.g. 0.3, 0.5, 0.8) — generate `clustree` plot

Store as `{RESOLUTION}`.

---

## Step 7 — Annotation references

All annotation methods are always run (SingleR, scType, Azimuth, ACT). The only choices are which references to use.

**SingleR reference — always ask two things:**

1. Primary immune reference: `ImmGenData()` (mouse) or `HumanPrimaryCellAtlasData()` (human)
2. If the dataset contains **non-immune populations** (epithelial, fibroblasts, endothelial etc.): also run `MouseRNAseqData()` (mouse) or `BlueprintEncodeData()` (human) as a second SingleR reference — these cover broader cell types that ImmGen/HPCA miss.

**scType tissue context:** Ask: "What tissue context should scType use? (e.g. Intestine, Lung, Liver, Brain, Immune system)" — do NOT use "Mouse" or "Human" as the tissue; use organ-level context. Full list in the ScTypeDB.

**Azimuth reference:** Ask which reference to use. Common: `pbmcref` (PBMC/immune), `bonemarrowref`, `lungref`, `kidneyref`, `heartref`, `mousecortexref`. Say `none` to skip. Note: `pbmcref` is a human PBMC reference and is poor for epithelial/stromal cell types.

Store as `{SINGLER_REF1}`, `{SINGLER_REF2}` (or NULL), `{SCTYPE_TISSUE}`, `{AZIMUTH_REFERENCE}`.

---

## Step 8 — Downstream analysis menu

**Differential abundance (propeller) is always run** in Rmd 05 — no need to ask. It needs ≥2 biological replicates per condition; warn if not met.

Ask: "Which DEG comparisons?" and "Which method: pseudobulk DESeq2 or Wilcoxon FindMarkers?"
Recommend pseudobulk DESeq2 when ≥2 biological replicates per condition.

Ask which optional downstream modules to generate (multi-select):
1. **GSEA** (clusterProfiler / fgsea / msigdbr) on DEG results — ranked by Wald statistic. Ask which cell type(s) to rank, **and which MSigDB collections to test** — this depends on the biology and should be chosen per project. Offer (multi-select), with guidance:
   - **H** — Hallmark, 50 broad programs. Good default baseline for any study.
   - **C2:CP** — canonical pathways (Reactome, KEGG, WikiPathways, PID, BioCarta). Best for signalling / metabolic / pathway-level questions.
   - **C5 GO:BP** — Gene Ontology Biological Process. Broad functional coverage.
   - **C7:IMMUNESIGDB** — immunologic perturbation signatures. For immune-focused studies.
   - **C8** — cell-type signatures. For lineage / differentiation questions.
   - **C6** — oncogenic signatures. For cancer / tumour studies.
   - **C3:TFT** — transcription-factor target sets. For regulatory / TF-activity questions.
   Store the chosen collections as `{GSEA_COLLECTIONS}` (a list of `collection`/`subcollection` pairs).
2. **GSVA** (per-cell or pseudobulk gene-set scores) — Hallmark / custom gene sets.
3. **Ligand-receptor (LIANA)** — intercellular communication. Ask which sender and receiver cell types are of interest.
4. **Population subclustering** — extract one lineage (e.g. myeloid) and re-cluster.

Note: **`msigdbr` must be present in the container image** for GSEA/GSVA. If the user selects either, confirm it is installed (`Rscript -e 'library(msigdbr)'` inside the SIF) before generating — it is not in the base image by default and requires a Docker rebuild, not a runtime `install.packages()`.

---

## Step 9 — SLURM resources

Ask: "How much memory/time?"
1. Small (< 10k cells) — 32G, 4h, 8 CPUs
2. Medium (10k–50k cells) — 64G, 8h, 16 CPUs
3. Large (> 50k cells) — 128G, 12h, 32 CPUs
4. Custom

For Rmd 02 QC (rough clustering only, no integration): can use half the resources of the integration step.

---

## Step 10 — Project name / author

Ask for author name and project title for Rmd headers.

---

## Step 11 — Generate Rmd 01: Data Import & QC

Write `{TODAY_YYMMDD}_{WD_NAME}_01_dataImport.Rmd`:

Key sections:
1. **Gene symbol conversion** (if starting from Seurat RDS with Ensembl IDs): load `features.tsv.gz` to build an Ensembl → symbol map, then after loading each RDS extract the counts matrix, rename rows from Ensembl IDs to gene symbols using `make.unique(feat_map$gene_symbol[match(ens_ids, feat_map$ensembl_id)])`, and rebuild with `CreateSeuratObject(counts = counts_mat)`. Do this BEFORE calculating percent.mt — otherwise `^mt-` will match nothing and all cells will show 0% MT.
2. Add metadata: `orig.ident`, `Condition`, `CellPop`, `percent.mt`, `l2.n_feature`, `l2.n_count`
3. **MAD-based per-sample QC filter** with 5% floor on max_mt: `max_mt <- max(min(100, med_mt + 3 * mad_mt), 5)` — the floor prevents CellBender-corrected data (near-zero MT) from filtering all cells
4. **scDblFinder** via direct SCE construction (never use `as.SingleCellExperiment()` or `DietSeurat()` — both call defunct `PackageCheck()`):
   ```r
   counts_mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
   sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts_mat))
   sce <- scDblFinder(sce)
   ```
5. Save each sample to `qc_rds/{sample_id}.seurat.rds` using exact sample IDs as filenames.
6. Save `QC_adaptive_thresholds.csv` with n_before, n_after, min_feat, max_feat, max_mt per sample.

---

## Step 12 — Generate Rmd 02: Merge & QC Review

**This Rmd is for QC review only. No Harmony integration here.**

Write `{TODAY_YYMMDD}_{WD_NAME}_02_qc.Rmd`:

Key sections:
1. **Load samples**: Use `list.files(file.path(CWD, "qc_rds"), pattern = "\\.seurat\\.rds$")`. After loading, immediately verify the count matches expected: `cat("Loaded", length(seurat_objects), "samples — expected {N_SAMPLES}\n")`. If the count is wrong, stop and investigate before continuing (stale files from previous runs can contaminate the merge).
2. **Cell cycle scoring**: Manual approach — `CellCycleScoring` calls defunct `slot=` API internally. Use:
   ```r
   norm_data <- GetAssayData(merged, assay = "RNA", layer = "data")
   s_present <- intersect(s.genes, rownames(norm_data))
   merged$S.Score <- colMeans(norm_data[s_present, , drop = FALSE])
   ```
3. **Rough clustering** for QC visualisation: `split(RNA, orig.ident)` → `SCTransform(glmGamPoi)` → `RunPCA` → `RunUMAP` → `FindClusters(resolution = 0.8)`. This is for inspection only — the integration Rmd will redo SCT from scratch.
4. **QC plots**: DimPlot by rough.clusters, FeaturePlot for percent.mt / nFeature / doublet score, VlnPlot, cluster QC table with doublet % per cluster.
5. **Filtering**: Remove doublets (`scDblFinder.class == "singlet"`) only. **Never remove whole clusters** — only per-cell QC filtering is appropriate. The cluster plots are for informational review, not for deciding which clusters to exclude.
6. **Gene-level filter**: keep if `rowSums(cm > 0) >= 3` OR `rowMeans(cm) >= 0.005`.
7. Save `seurat.filtered.rds`.

---

## Step 13 — Generate Rmd 03: Integration

Write `{TODAY_YYMMDD}_{WD_NAME}_03_integration.Rmd`:

Key sections:
1. Load `seurat.filtered.rds`.
2. **SCT+PCA checkpoint** — saves the expensive SCTransform+PCA so re-runs (e.g. to change resolution) skip straight to Harmony. Auto-invalidates if cell count changed (i.e. user re-ran Rmd 02 with different filters):
   ```r
   CHECKPOINT <- file.path(OUT, "checkpoint_sct_pca.rds")
   filtered_n <- ncol(seurat.filtered)
   if (file.exists(CHECKPOINT)) {
     chk <- readRDS(CHECKPOINT)
     if (ncol(chk) == filtered_n) {
       seurat.filtered <- chk
     } else {
       cat("Stale checkpoint — recomputing\n"); file.remove(CHECKPOINT)
     }
     rm(chk)
   }
   if (!file.exists(CHECKPOINT)) {
     DefaultAssay(seurat.filtered) <- "RNA"
     seurat.filtered[["SCT"]] <- NULL
     seurat.filtered[["RNA"]] <- split(seurat.filtered[["RNA"]], f = seurat.filtered$orig.ident)
     seurat.filtered <- SCTransform(seurat.filtered, method = "glmGamPoi", verbose = FALSE)
     seurat.filtered <- RunPCA(seurat.filtered, npcs = 50, verbose = FALSE)
     saveRDS(seurat.filtered, CHECKPOINT)
   }
   ```
3. **Harmony integration** using `RunHarmony` directly (do NOT use `IntegrateLayers(HarmonyIntegration)` — it calls `harmony::HarmonyMatrix` which is not exported in harmony >= 1.0):
   ```r
   cat("PCA rows:", nrow(Embeddings(seurat.filtered, "pca")),
       "| Object cells:", ncol(seurat.filtered), "\n")
   set.seed(42)
   seurat.integrated.harmony <- RunHarmony(
     seurat.filtered,
     group.by.vars  = "orig.ident",
     reduction      = "pca",
     reduction.save = "harmony",
     verbose        = FALSE
   )
   ```
   The diagnostic print confirms the PCA and metadata dimensions match before calling RunHarmony. If they don't match, it indicates a data integrity problem (e.g. a stale file was loaded in Rmd 02) — fix the data, not the code.
4. `JoinLayers(assay = "RNA")` → `FindNeighbors(reduction = "harmony")` → `FindClusters(resolution = {RESOLUTION})` → `RunUMAP(reduction = "harmony", reduction.name = "umap.harmony")`.
5. UMAP plots: by cluster, by Condition, by CellPop, by orig.ident (sample).
6. **FindAllMarkers checkpoint** — saves the slow FindAllMarkers result for re-use when tweaking ACT format or gene count:
   ```r
   MARKERS_CHECKPOINT <- file.path(OUT, "all.markers.rds")
   if (file.exists(MARKERS_CHECKPOINT)) {
     all.markers <- readRDS(MARKERS_CHECKPOINT)
   } else {
     DefaultAssay(seurat.integrated.harmony) <- "RNA"
     all.markers <- FindAllMarkers(seurat.integrated.harmony, only.pos = TRUE,
                                   min.pct = 0.25, logfc.threshold = 0.25,
                                   test.use = "wilcox", verbose = FALSE)
     saveRDS(all.markers, MARKERS_CHECKPOINT)
   }
   ```
7. **ACT input file** — format: one line per cluster, `cluster_X: gene1, gene2, ...` (up to 200 genes, sorted by significance). This is what the ACT online tool at https://www.thecellvision.org/celltype/ expects:
   ```r
   act_lines <- all.markers %>%
     filter(p_val_adj < 0.05) %>%
     group_by(cluster) %>%
     arrange(p_val_adj, desc(avg_log2FC)) %>%
     slice_head(n = 200) %>%
     ungroup() %>%
     mutate(Cluster = paste0("cluster_", cluster)) %>%
     group_by(Cluster) %>%
     summarise(genes = paste(gene, collapse = ", "), .groups = "drop") %>%
     arrange(as.integer(sub("cluster_", "", Cluster))) %>%
     mutate(line = paste0(Cluster, ": ", genes))
   writeLines(act_lines$line, file.path(OUT, "ACT_cluster_input_harmony.txt"))
   ```
8. Save `seurat.filtered.rds` and `seurat.integrated.harmony.rds`.

---

## Step 14 — Generate Rmd 04: Cell Type Assignment

Write `{TODAY_YYMMDD}_{WD_NAME}_04_celltypeAssignment.Rmd`:

**Annotation order: SingleR → scType → Azimuth (last)**. Azimuth's `RunAzimuth` converts mouse→human gene names in-place, corrupting SCT scale.data that scType needs.

**Critical**: `RunAzimuth` permanently overwrites the RNA assay rownames (e.g. `Ptprc` → `PTPRC`). This means `seurat.annotated.rds` (saved after Azimuth) has human gene names, while `all.markers.rds` (from FindAllMarkers in Rmd 03) retains the original gene names. Any downstream Rmd that needs gene expression (DotPlot, FeaturePlot, FindMarkers) must load `seurat.integrated.harmony.rds`, not `seurat.annotated.rds`.

Key sections:
1. **SingleR** — build SCE directly (not via DietSeurat):
   ```r
   counts_mat    <- GetAssayData(seurat.integrated.harmony, assay = "RNA", layer = "counts")
   logcounts_mat <- log1p(t(t(counts_mat) / colSums(counts_mat)) * 1e4)
   sce <- SingleCellExperiment::SingleCellExperiment(
     assays = list(counts = counts_mat, logcounts = logcounts_mat))
   ```
   Run SingleR once per reference (ImmGen/HPCA for immune cells; MouseRNAseqData/BlueprintEncode if non-immune populations present). Store results as separate metadata columns: `singler_immgen_main`, `singler_immgen_fine`, `singler_broad_label`.

2. **scType** — use `GetAssayData(assay = "SCT", layer = "scale.data")` for the score matrix. Set tissue to `{SCTYPE_TISSUE}` (organ-level, not "Mouse"/"Human"). Wrap in `tryCatch` since it requires internet access.

3. **Azimuth** — `RunAzimuth(reference = "{AZIMUTH_REFERENCE}")`. Wrap in `tryCatch`. Adds `predicted.celltype.l1`, `predicted.celltype.l2`, `predicted.celltype.l2.score`.

4. **ACT output loading** — load `ACT_Annotationresults_top10.txt` from `{RESULTS_DIR}` if present; graceful fallback if missing. Parse: `filter(Order == 1)`, strip `cluster_` prefix from Cluster column.

5. **Comparison table** per cluster — one row per cluster, majority label from each tool.

6. **Consensus vote** across all tools that ran. Flag clusters as `uncertain` if `n_agree < 3`. Write `annotation_consensus_per_cluster.xlsx`.

7. **RenameIdents scaffold** — auto-generate and `cat()` a ready-to-paste `RenameIdents()` call from consensus labels. Leave the actual call commented out (`# TODO`) for user review.

8. Cell type proportions per condition.

9. Save `seurat.annotated.rds`.

---

## Step 15 — Generate Rmd 05: Final Annotation (always)

After the user has reviewed `annotation_consensus_per_cluster.xlsx` from Rmd 04, they provide a tab-separated file (e.g. `harmony.clusters_consensus.txt`) with two columns — `harmony.clusters` and `consensus` — mapping each cluster number to their chosen cell type label.

Write `{TODAY_YYMMDD}_{WD_NAME}_05_finalAnnotation.Rmd`:

Key sections:
1. **Load `seurat.integrated.harmony.rds`** — NOT `seurat.annotated.rds`. The integrated object retains original gene names (pre-Azimuth), which are required for DotPlot and marker gene lookup. The annotated object has Azimuth-converted human gene names and cannot be used for gene expression visualisation.
2. Read the consensus txt file: `read.table(..., header = TRUE, sep = "\t")`.
3. Apply `RenameIdents` from the consensus mapping and store as `seurat[["celltype"]]`.
4. **UMAP plots**: overall by celltype, split by Condition, split by CellPop.
5. **Top 4 markers per cell type** from `all.markers.rds` (saved in Rmd 03):
   - Map cluster numbers → cell types using the consensus file
   - `group_by(celltype) %>% arrange(p_val_adj, desc(avg_log2FC)) %>% slice_head(n = 4)`
   - **Always filter marker genes to those present in the object** before calling DotPlot:
     `marker_genes <- marker_genes[marker_genes %in% rownames(seurat)]`
     Genes from `all.markers.rds` that were removed by the gene-level QC filter in Rmd 02 will cause `FetchData()` to abort with "None of the requested variables were found" if not filtered.
6. **DotPlot** — order cell types meaningfully (immune → epithelial → stromal). Use `scale_colour_distiller(palette = "RdBu", direction = -1)`.
7. Cell type proportions per condition (stacked bar). Also write per-sample, per-condition, and group-mean proportion tables to xlsx.
8. **Differential abundance — propeller via limma** (always run if ≥2 replicates/condition). **`speckle` is NOT in the image; `limma` is** — use a manual reimplementation of Phipson et al. (2022): arcsine-sqrt-transform per-sample proportions, then a limma moderated t-test.
   ```r
   library(limma)
   propeller_limma <- function(metadata, celltype_col, sample_col, group_col, ref_group) {
     # as.character() to drop unused factor levels — subsetting to a CellPop keeps
     # phantom levels that create empty rows and break the design matrix.
     samp  <- as.character(metadata[[sample_col]])
     ctype <- as.character(metadata[[celltype_col]])
     grp_v <- as.character(metadata[[group_col]])
     prop  <- table(samp, ctype); prop <- prop / rowSums(prop)
     grp   <- factor(grp_v[match(rownames(prop), samp)],
                     levels = c(ref_group, setdiff(unique(grp_v), ref_group)))
     y      <- asin(sqrt(t(prop)))            # cell types x samples (limma = genes x samples)
     fit    <- eBayes(lmFit(y, model.matrix(~grp)), trend = TRUE, robust = TRUE)
     topTable(fit, coef = 2, number = Inf, sort.by = "P")  # add prop_ratio etc. for output
   }
   ```
   Run **separately per sorted CellPop** (e.g. CD45 and EpCAM each) — proportions are compositional within a fraction, so mixing fractions is meaningless. Output columns to match: `celltype, prop_ref, prop_test, prop_ratio, t_stat, P.Value, FDR`; write to `propeller_results.xlsx`.
   Plot as a volcano with **`log2(prop_ratio)`** on the x-axis (there is no logFC column — `prop_ratio` is test/ref proportion) and `-log10(FDR)` on the y-axis; colour by `FDR < 0.05`.
9. Save as `seurat.annotated.rds` (overwrites the Azimuth-corrupted version with the clean annotated object).

SLURM resources: 64G, 8 CPUs, 2h — no heavy computation, just loading and plotting.

---

## Step 16 — Generate Rmd 06: DEG (if requested)

*(Only if user did not skip DEG in Step 8.)*

Write `{TODAY_YYMMDD}_{WD_NAME}_06_DEG.Rmd`. Load `seurat.annotated.rds`, set `DefaultAssay = "RNA"`, `JoinLayers(assay = "RNA")` first.

**Pseudobulk DESeq2 — the recipe that works** (loop over each cell type):
1. Aggregate with `group.by = "orig.ident"` **only**, then join `Condition` from metadata separately. Do NOT use `group.by = c("orig.ident","Condition")` — it mangles column names and couples the two factors:
   ```r
   bulk <- AggregateExpression(sub, group.by = "orig.ident", assays = "RNA",
                               return.seurat = FALSE)$RNA
   colnames(bulk) <- gsub("-", "_", colnames(bulk))  # AggregateExpression turns _ into -
   ```
2. **`round(bulk)`** before DESeq2 — `countData` must be integers. These are CellBender-corrected counts, not raw (state this in any methods text).
3. Require **≥10 cells per sample** for that cell type; drop samples below. Skip the cell type entirely if fewer than 4 samples (or <2 per group) survive.
4. Gene pre-filter inside DESeq2: `dds <- dds[rowSums(counts(dds) >= 5) >= 2, ]`.
5. `design = ~Condition`, `contrast = c("Condition", "<test>", "<ref>")` with `factor(levels = c(ref, test))` so the reference is explicit.
6. Write one CSV per cell type: `DEG/{CellPop}__{celltype}_DESeq2.csv` (keep `gene`, `log2FoldChange`, `stat`, `pvalue`, `padj`). The `stat` (Wald statistic) column is needed for GSEA ranking.
7. Significance: `padj < 0.05 & abs(log2FoldChange) > 0.25`. Write a `DEG_summary.xlsx` (n_up / n_down per cell type).

**Wilcoxon alternative** (if <2 replicates): `FindMarkers(test.use = "wilcox")` per cell type, `DefaultAssay = "RNA"`, `JoinLayers` first.

**Targeted gene-list visualisation** (optional, very useful): for curated pathway gene lists, `DotPlot(sub, features = genes, group.by = "Condition") + coord_flip()`, filtering `genes <- genes[genes %in% rownames(sub)]` first.

---

## Step 17 — Generate Rmd 07: GSEA (if requested)

Write `{TODAY_YYMMDD}_{WD_NAME}_07_<celltype>_GSEA.Rmd`. Requires `msigdbr`, `clusterProfiler`, `enrichplot` in the image (confirm in Step 8).

1. Load the cell type's DEG CSV. **Rank by Wald statistic (`stat`)**, not log2FC — it captures effect size and significance together:
   ```r
   ranked <- deg %>% filter(!is.na(stat), !is.na(gene)) %>%
     distinct(gene, .keep_all = TRUE) %>% arrange(desc(stat))
   gene_ranks <- setNames(ranked$stat, ranked$gene)
   ```
2. Gene sets via `msigdbr`. **Use the current API** (the image has msigdbr 26.x): `collection=` / `subcollection=`, NOT the deprecated `category=` / `subcategory=`. **For mouse, use `db_species = "HS", species = "Mus musculus"`** — the mouse-native database (`db_species = "MM"`) is not installed, so msigdbr ortholog-maps the human DB to mouse symbols (returns the `gene_symbol` column in mouse symbols). Loop over the user's `{GSEA_COLLECTIONS}`:
   ```r
   get_sets <- function(coll, subcoll = NULL) {
     msigdbr(db_species = "HS", species = "{ORGANISM}",
             collection = coll, subcollection = subcoll) %>%
       dplyr::select(gs_name, gene_symbol)
   }
   # e.g. H -> get_sets("H"); GO:BP -> get_sets("C5","GO:BP");
   #      immune -> get_sets("C7","IMMUNESIGDB"); pathways -> get_sets("C2","CP:REACTOME")
   ```
   (For human data, drop `species`/use `species = "Homo sapiens"`.)
3. `GSEA(gene_ranks, TERM2GENE = sets, minGSSize = 10, maxGSSize = 500, pAdjustMethod = "BH")` once per collection. **Do not set `nPerm`** — clusterProfiler defaults to `method = "multilevel"` (fgsea), which ignores `nPerm` and is both faster and more accurate than the obsolete fixed-permutation approach.
4. Filter significant at `p.adjust < 0.05`. NES dot plots (colour by direction), and `gseaplot2` enrichment curves for top hits.
5. **`gseaplot2` chunks need `fig.width >= 14`** or the p-value table is cropped on the right.
6. Export leading-edge genes (`core_enrichment`) per significant set to xlsx, one sheet per collection.

---

## Step 18 — Generate Rmd 08: Ligand-Receptor / LIANA (if requested)

Write `{TODAY_YYMMDD}_{WD_NAME}_08_LR_LIANA.Rmd`. Load `seurat.annotated.rds`.

1. **LIANA v0.1.14 requires a Seurat v3 `Assay`, not v5 `Assay5`.** `as(obj[["RNA"]], "Assay")` fails on split layers ("subscript out of bounds"). Rebuild the object instead:
   ```r
   seurat <- JoinLayers(seurat)
   options(Seurat.object.assay.version = "v3")
   seurat <- CreateSeuratObject(
     counts    = GetAssayData(seurat, assay = "RNA", layer = "counts"),
     meta.data = seurat@meta.data)
   seurat <- NormalizeData(seurat)
   Idents(seurat) <- seurat$celltype
   ```
2. **Run on the full integrated object** (all cell types), then filter results to sender/receiver cell types of interest. Do NOT subset before running.
3. Resource: `"MouseConsensus"` for mouse (native symbols, no homolog conversion) or `"Consensus"` for human. The ligand column in the resource is **`source_genesymbol`** (not `genesymbol_intercell_source`).
4. `liana_wrap(seurat, resource = "MouseConsensus", expr_prop = 0.1)` then `liana_aggregate()`. Lower `aggregate_rank` = more methods agree the interaction is active.
5. Filter `source %in% SENDERS, target %in% RECEIVERS`; `liana_dotplot(...)`. Run per condition separately if comparing.
6. Annotate prioritised interactions by downstream programme where biologically known (e.g. Il1b→Il1r1 = NF-κB; Il6/Osm→Il6st = STAT3; Areg/Egf→Egfr = proliferation).

---

## Step 19 — Generate SLURM scripts

One script per Rmd. SLURM template:

```bash
#!/bin/bash
#SBATCH -N 1
#SBATCH -n {CPUS}
#SBATCH --mem={MEM}
#SBATCH -t {TIME}
#SBATCH -p bcc
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user={USER_EMAIL}
#SBATCH -o {LOG_FILE}

CWD={CWD}
OUT={RESULTS_DIR}
SIF=/net/bmc-lab3/data/bcc/shared/singularity_images/r4_6_0_python3_singlecell_visium_latest.sif
R_LIBS=/net/bmc-lab3/data/bcc/shared/azimuth_refs

mkdir -p ${OUT}
module load singularity/3.5.0
export SINGULARITYENV_R_LIBS_USER=${R_LIBS}

singularity exec \
    --bind ${CWD}:${CWD} \
    --bind ${R_LIBS}:${R_LIBS} \
    ${SIF} \
    Rscript -e "rmarkdown::render(
        '${CWD}/{RMD_FILE}',
        output_dir    = '${OUT}',
        knit_root_dir = '${CWD}'
    )"
```

Scripts: `run_01_dataImport.sh`, `run_02_qc.sh`, `run_03_integration.sh`, `run_04_celltypeAssignment.sh`, `run_05_finalAnnotation.sh`, and (if requested) `run_06_DEG.sh`, `run_07_<celltype>_GSEA.sh`, `run_08_LR_LIANA.sh`.

Rmd 02 QC can use half the memory/CPUs of the integration step (no Harmony needed).
Rmd 05 Final Annotation needs only 64G / 8 CPUs / 2h — no heavy computation.
Rmd 07 GSEA and Rmd 08 LIANA: 64G / 8 CPUs / 2–4h is ample.

---

## Step 20 — Generate summary HTML report

After all Rmds are written, also generate a standalone `{WD_NAME}_summary_report.html` in `{RESULTS_DIR}`. This is a self-contained Bootstrap HTML page (no external dependencies) that:
- Summarises the study design, QC statistics, integration results, and cluster annotation consensus table
- Describes the final Seurat object in detail: all `@meta.data` columns with types/values, assays (RNA + SCT), reductions (pca, harmony, umap.harmony), and a quick-start R code block
- Links to each individual HTML report and to `seurat.annotated.rds` using the base URL pattern provided by the user (or ask if not known)

Ask: "Is there a web-accessible URL prefix for the results directory (e.g. https://bmc-data.mit.edu/BCC/Oli/{project}/)? This will be used to create clickable links in the summary report."

---

## Final submission instructions

Display after all files are written:

```
Files written:
  {TODAY_YYMMDD}_{WD_NAME}_01_dataImport.Rmd
  {TODAY_YYMMDD}_{WD_NAME}_02_qc.Rmd
  {TODAY_YYMMDD}_{WD_NAME}_03_integration.Rmd
  {TODAY_YYMMDD}_{WD_NAME}_04_celltypeAssignment.Rmd
  run_01_dataImport.sh
  run_02_qc.sh
  run_03_integration.sh
  run_04_celltypeAssignment.sh

Submission order:
  sbatch -p bcc run_01_dataImport.sh
  # Wait for completion, then:
  sbatch -p bcc run_02_qc.sh
  # Review 02_qc.html — check QC plots, verify cell counts look right
  # Then:
  sbatch -p bcc run_03_integration.sh
  # While 03 runs: submit ACT_cluster_input_harmony.txt to
  #   https://www.thecellvision.org/celltype/
  # Save ACT result as ACT_Annotationresults_top10.txt in results dir
  sbatch -p bcc run_04_celltypeAssignment.sh
  # Review annotation_consensus_per_cluster.xlsx
  # Fill in harmony.clusters_consensus.txt with your final cell type labels
  sbatch -p bcc run_05_finalAnnotation.sh
  # Review propeller_results.xlsx and the annotated UMAP. Then any requested:
  sbatch -p bcc run_06_DEG.sh
  sbatch -p bcc run_07_<celltype>_GSEA.sh
  sbatch -p bcc run_08_LR_LIANA.sh
```

---

## Notes for the assistant

### Cluster-based analysis
- **Never exclude whole clusters from QC.** Only per-cell filtering is appropriate (percent.mt, nFeature, doublet score). The cluster QC plots in Rmd 02 are for informational review only.
- **BAD_CLUSTERS parameter is banned** from the pipeline. Remove it entirely.

### Gene naming
- Pre-existing Seurat RDS files from CellRanger often use Ensembl IDs (e.g. `ENSMUSG00000019987`) as gene names, not symbols. This causes: (a) `PercentageFeatureSet(pattern = "^mt-")` to return 0 for all cells, (b) SingleR, scType, ACT, and Azimuth to fail silently or return no matches.
- Always check gene names after loading the first RDS. If they look like Ensembl IDs, add the conversion step (see Step 11). The `features.tsv.gz` from any sample's CellRanger output contains the Ensembl→symbol mapping.
- Use `make.unique()` after mapping to handle duplicate gene symbols (common with mitochondrial or pseudogenes).

### Data integrity before integration
- After loading `qc_rds/`, always print and verify the number of samples matches the expected count. Stale files from previous pipeline attempts (e.g. `SeuratProject.seurat.rds` from a pool-level combined object) can be inadvertently picked up by `list.files()`, inflating the cell count by 2–3× and causing all downstream failures. This is easier to catch early than to diagnose from a failed Harmony run.

### Harmony integration
- Use `RunHarmony(obj, group.by.vars = "orig.ident", reduction = "pca", reduction.save = "harmony", verbose = FALSE)` directly.
- Do NOT use `IntegrateLayers(method = HarmonyIntegration, ...)` — it calls `harmony::HarmonyMatrix` internally which is not exported in harmony >= 1.0.
- Do NOT pass `assay.use =` to `RunHarmony` — removed in harmony v1.2+. Harmony operates on the PCA embedding directly.
- Do NOT call `JoinLayers(assay = "SCT")` — `SCTAssay` is a v4-style class that does not support `JoinLayers`.
- Before calling `RunHarmony`, print `nrow(Embeddings(obj, "pca"))` vs `ncol(obj)` to verify they match. If they don't, the root cause is always a data integrity problem (wrong files loaded, wrong object), not an API issue.

### SeuratObject v5 API
- `slot=` argument is removed everywhere → use `layer=`.
- `CellCycleScoring` calls defunct `slot=` internally. Use manual `GetAssayData(layer = "data")` + `colMeans` approach.
- `as.SingleCellExperiment()` and `DietSeurat()` call defunct `PackageCheck()`. Build SCE directly from counts matrix.
- `JoinLayers(assay = "RNA")` is needed before `NormalizeData`, `FindAllMarkers`, and `FindMarkers`.
- Seurat plotting functions (`DimPlot`, `VlnPlot`, `FeaturePlot`) work correctly in the r4_6_0 container — use them normally.

### SCTransform on split layers
- `split(RNA, f = orig.ident)` → `SCTransform` → `RunPCA` is the correct Seurat v5 pattern for Harmony. Do this in the integration Rmd (not the QC Rmd).
- The resulting SCT assay stores scale.data correctly for all cells as an `SCTAssay` object.
- The rough clustering in Rmd 02 QC can also use this pattern, or use standard `NormalizeData` → `ScaleData` — either works for the QC step.

### ACT file format
- The ACT online tool (https://www.thecellvision.org/celltype/) expects: one line per cluster, format `cluster_X: gene1, gene2, gene3, ...` Up to 200 genes per cluster, sorted by significance. This is NOT a tabular format.
- Always use `ungroup()` before `select()` in dplyr chains that had a `group_by()` — otherwise the grouping variable persists as an extra column in the output.
- When submitting to ACT, select the relevant tissue (e.g. "Colon" for colorectal data) and the correct species.

### Annotation references
- For datasets with **mixed immune + epithelial/stromal populations**, always run two SingleR references: one immune-focused (ImmGen / HPCA) and one broad (MouseRNAseqData / BlueprintEncode). ImmGen alone will force-assign epithelial clusters to the nearest immune label.
- scType tissue should be organ-level (e.g. "Intestine", "Lung") not species-level ("Mouse") — the organ context dramatically improves epithelial and stromal cell type resolution.
- Azimuth `pbmcref` is appropriate for immune cells but cannot annotate epithelial, fibroblast, or endothelial populations — set user expectation accordingly.
- Consensus labels from automated tools are always provisional — manual validation against canonical marker genes is required before publication.

### Azimuth gene name corruption
- `RunAzimuth` permanently overwrites the RNA assay `rownames` by converting species-specific gene symbols to human symbols (e.g. mouse `Ptprc` → human `PTPRC`). This happens in-place and cannot be undone without reloading.
- Consequence: `seurat.annotated.rds` (saved after Azimuth) has human gene names. Any code that calls `DotPlot`, `FeaturePlot`, `FindMarkers`, or `GetAssayData` by gene name will fail or silently return nothing if run on `seurat.annotated.rds` with original (mouse) gene lists.
- Always load `seurat.integrated.harmony.rds` (pre-Azimuth) for gene expression visualisation and marker analysis. Use `seurat.annotated.rds` only for metadata (cell type labels, Azimuth predictions, SingleR labels) and UMAP coordinates.
- In Rmd 05 Final Annotation: always load `seurat.integrated.harmony.rds`, not `seurat.annotated.rds`.
- When building DotPlot from `all.markers.rds`: always filter features to `marker_genes[marker_genes %in% rownames(seurat)]` before passing to `DotPlot()` — genes removed by the gene-level QC filter in Rmd 02 will cause an abort otherwise.

### Checkpoints
- Add two checkpoints in the integration Rmd:
  1. After SCTransform + RunPCA (most expensive step: ~30–60 min on 50k cells). Auto-invalidate by comparing cell count.
  2. After FindAllMarkers (can take 30–90 min). This allows tweaking ACT format or gene count thresholds without re-running the marker test.
- Document in each checkpoint how to force a re-run (delete the `.rds` file).

### Differential abundance (propeller)
- **`speckle` is NOT in the image** — do not call `speckle::propeller()`. Use the manual limma reimplementation (arcsine-sqrt + `eBayes(trend=TRUE, robust=TRUE)`); `limma` is installed.
- Use `as.character()` on the metadata columns before `table()` — subsetting to a CellPop leaves phantom factor levels that break the design matrix.
- Run **per sorted CellPop separately** (proportions are compositional within a fraction).
- Output has **no logFC column** — use `log2(prop_ratio)` (test/ref) for volcano x-axes. Columns: `celltype, prop_ref, prop_test, prop_ratio, t_stat, P.Value, FDR`.

### Pseudobulk DESeq2
- Aggregate with `group.by = "orig.ident"` only, join `Condition` separately. `AggregateExpression` converts `_` → `-` in column names — fix with `gsub`.
- `round()` the aggregated matrix before DESeq2 (integer requirement). These are CellBender-corrected counts, not raw — say so in methods.
- Require ≥10 cells/sample for a cell type; skip the cell type if <4 samples (or <2/group) survive.
- Keep the `stat` (Wald) column in output CSVs — it is the correct GSEA ranking metric.

### GSEA (clusterProfiler / fgsea)
- `msigdbr` (26.x), `clusterProfiler`, `enrichplot` are in the image (verified). `msigdbr` is added via Docker rebuild, never a runtime `install.packages()`.
- **The collection(s) to test are a per-project choice driven by biology** — always ask (Step 8), never hardcode Hallmark/C7/GO:BP.
- **msigdbr 26.x API:** use `collection=` / `subcollection=` — `category=` / `subcategory=` are deprecated and will error. Columns are `gs_name` / `gene_symbol` / `gs_collection`.
- **Mouse:** use `db_species = "HS", species = "Mus musculus"` — the mouse-native DB (`db_species = "MM"`) is NOT installed; msigdbr ortholog-maps human → mouse symbols automatically.
- Rank genes by DESeq2 Wald statistic (`stat`), not log2FC.
- Do not set `nPerm` — `GSEA()` uses `method = "multilevel"` (fgsea) by default, which ignores it and is more accurate than the obsolete fixed-permutation `nPerm=1000` approach.
- `gseaplot2` needs `fig.width >= 14` or the p-value table is cropped.

### LIANA ligand-receptor
- v0.1.14 requires a Seurat v3 `Assay`. `as(obj[["RNA"]], "Assay")` fails on v5 split-layer objects ("subscript out of bounds"). Rebuild: `JoinLayers` → `options(Seurat.object.assay.version="v3")` → `CreateSeuratObject(counts = GetAssayData(...,layer="counts"), meta.data=...)` → `NormalizeData`.
- Run on the full integrated object, filter results to senders/receivers afterward — never subset first.
- Use `"MouseConsensus"` resource for mouse (native symbols). The ligand column is `source_genesymbol`.
- Lower `aggregate_rank` = higher-confidence interaction (more methods agree).

### Figure conventions (recurring user preferences)
- x-axis labels: `theme(axis.text.x = element_text(angle = 45, hjust = 1))`.
- Variable-size dot plots: always add `guides(size = guide_legend(override.aes = list(...)))` or the legend collapses to a single dot size.
- Prefer `facet_wrap(~ var, nrow = 1)` over vertically stacked `facet_grid` for compactness; scale figure height to the tallest facet, not the sum.

### Deliverables (Word / PowerPoint / methods)
- `pandoc` 3.9 is in the container: `pandoc methods.md --to docx -o methods.docx` and `pandoc slides.md --to pptx -o slides.pptx` (`##` headings become slides; `::: notes :::` blocks become speaker notes). `--resource-path={CWD}` so embedded image paths resolve.
- `officer` (R) and `python-pptx` are NOT in the image — use pandoc, not those.
- Export ggplot figures to PNG (`ggsave(..., dpi = 300, bg = "white")`) for slide assembly rather than embedding HTML.
- Tool versions for methods text: Cell Ranger in `pipeline_info/nf_core_scrnaseq_software_mqc_versions.yml`; CellBender only in `cellrangermulti/*/cellbender_removebackground/*.log`.

### SLURM
- Always submit with `sbatch -p bcc script.sh` on the command line. The `#SBATCH -p bcc` header alone silently fails on this cluster.
- Never run R, Singularity, or heavy computation on the login node.
- `module load singularity/3.5.0` (not `module add`).

### General
- R line continuation: never use `\` — R expressions continue when a line ends with an open operator (`+`, `,`, `|>`, `(`). Keep `if/else` on one line.
- `set.seed(42)` immediately before every `RunUMAP` and `RunHarmony`.
- FindAllMarkers: always `DefaultAssay = "RNA"` with `JoinLayers` first.
- Raw CellRanger / h5ad files are read-only. Write all outputs to `{RESULTS_DIR}`.
- Always install packages in the Docker image — never suggest `install.packages()` at runtime.
- `gh` CLI is not available on this cluster.
