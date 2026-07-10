# Bulk RNAseq Downstream Analysis Skill (tximport / DESeq2 / edgeR / HPC)

You are helping the user run **bulk RNA-seq downstream analysis** on an HPC cluster using SLURM and Singularity, starting from **nf-core/rnaseq `star_salmon` output** (`quant.sf` files). The container is `docker://yannvrb56/bulkrnaseq` (pull to SIF at `/net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif`). Generate R Markdown files (one per processing stage) and a SLURM submission script for each. Walk through the steps below in order. Perform all automated steps (file scans, mode detection, QC parsing) silently before presenting results.

**This pipeline supports two library types and the difference between them is the single most error-prone part of the workflow — get the mode right before generating anything:**

| | Paired-end (standard mRNA) | 3′ DGE (3' tag) |
|---|---|---|
| Length normalisation | **Yes** — transcript length matters | **No** — 3' tag count is length-independent |
| Expression units | **TPM** (`txi$abundance`) → `l2tpm` | **CPM** (edgeR `cpm()`) → `l2cpm` |
| DESeq2 construction | `DESeqDataSetFromTximport(txi, ...)` (uses length offset) | `DESeqDataSetFromMatrix(countData = round(txi$counts), ...)` (no offset) |
| edgeR used? | No | Yes (`DGEList` → `cpm`) |

**Violating this distinction silently corrupts every downstream result.** A 3′ DGE dataset analysed with length-aware tximport offsets will produce wrong fold changes; a paired-end dataset analysed as CPM-only loses the length correction. See the **Notes for the assistant → Library mode** section.

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

## Step 2 — Verify / pull Singularity image

```bash
ls /net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif 2>/dev/null
```

- **Found**: verify it has the required packages (run silently inside the SIF):
  ```bash
  module load singularity/3.5.0
  singularity exec /net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif \
    Rscript -e 'for (p in c("tximport","DESeq2","edgeR","apeglm","fgsea","msigdbr","GseaVis","openxlsx","tidyverse","matrixStats","ggrepel","showtext")) cat(sprintf("%-12s %s\n", p, requireNamespace(p, quietly=TRUE)))'
  ```
  `tximport`, `DESeq2`, `edgeR`, `apeglm`, `fgsea` are mandatory. `msigdbr` is needed only if GSEA is requested; `GseaVis` only for `gseaNb` enrichment curves; `showtext` only for unicode-label PCA figures. (`rtracklayer` is no longer required — annotation is staged as TSVs, not parsed from the GTF in R.) If a mandatory one is `FALSE`, stop and tell the user the image must be rebuilt (packages go in the Docker image, never `install.packages()` at runtime).
- **Not found**: tell the user to pull it on a compute node first (not the login node):
  ```bash
  cd /net/bmc-lab3/data/bcc/shared/singularity_images
  singularity pull bulkrnaseq_latest.sif docker://yannvrb56/bulkrnaseq:latest
  ```

---

## Step 3 — Detect input data structure

**Auto-detect silently** — find the nf-core/rnaseq salmon quantifications and the multiqc stats:

```bash
# Per-sample salmon quantifications (the import unit)
find {CWD} -maxdepth 8 -name "quant.sf" 2>/dev/null | head -40

# nf-core multiqc general stats (for QC join + mode detection)
find {CWD} -maxdepth 8 -path "*multiqc*star_salmon*" -name "multiqc_general_stats.txt" 2>/dev/null | head -3

# nf-core's own tx2gene (PREFERRED annotation source — exact match to quant.sf)
find {CWD} -maxdepth 8 -name "salmon.merged.tx2gene.tsv" 2>/dev/null | grep "star_salmon" | head -3

# Reference GTF (only needed for gene_biotype, which the tx2gene.tsv lacks)
find {CWD} /Genomes -maxdepth 4 -name "*.gtf" 2>/dev/null | head -10
```

Each `quant.sf` lives in a sample folder (`.../star_salmon/<SAMPLE>/quant.sf`). The parent folder name is the sample id. Tell the user how many quantified samples were found and where. If no `quant.sf` is found, stop — this pipeline starts from nf-core/rnaseq `star_salmon` output.

Store the multiqc stats path as `{MULTIQC_STATS}` and the salmon base dir as `{SALMON_DIR}`.

---

## Step 4 — Detect library mode (paired-end vs 3′ DGE) and confirm

**Auto-detect silently** from `{MULTIQC_STATS}`:

```bash
# Paired-end signals: "Read 1"/"Read 2" sub-rows, non-zero properly_paired, ISR/ISF library type
grep -E " Read [12]\t| Read [12]$" {MULTIQC_STATS} | head -2   # present => paired-end
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i ~ /properly_paired_percent/)pp=i; if($i ~ /salmon-library_types/)lt=i}} NR>1 && $1 !~ /Read [12]$/{print $1"\t"$pp"\t"$lt; if(NR>4)exit}' {MULTIQC_STATS}
```

**Decision rule:** `Read 2` rows present OR `properly_paired_percent > 0` → **paired-end**. Otherwise (single-end, library type `SR`/`U`, properly_paired empty/0) → **3′ DGE**.

Present the inferred mode and the evidence, then **ask the user to confirm**:
"Detected **{MODE}** (evidence: …). Is that correct? [paired-end / 3pDGE]"

Store as `{MODE}` ∈ {`paired-end`, `3pDGE`}. This drives Step 9 (import) and Step 12 (DESeq2 construction).

---

## Step 5 — Organism & annotation

Ask: "What organism / genome build?" (mouse `mm39`, human `GRCh38`, etc.) → store `{ORGANISM}`.

**Salmon `quant.sf` quantifies at the transcript level — `tx2gene` is required for gene-level summarisation; do not skip it.** Two pieces of annotation are needed:

1. **`tx2gene` (transcript_id → gene_id):** **prefer nf-core's own `star_salmon/salmon.merged.tx2gene.tsv`** (detected in Step 3) over re-parsing the GTF. It is exact-matched to the quantified transcripts (no transcript-version mismatch) and reads in milliseconds. Stage it: `cp .../star_salmon/salmon.merged.tx2gene.tsv {CWD}/ref/tx2gene.tsv`. Columns: `transcript_id, gene_id, gene_name` (header present); tximport uses the first two.
2. **`gene2typesym` (gene_id → gene_biotype, gene_name):** the nf-core tx2gene **lacks `gene_biotype`**, which the `protein_coding` filter needs — so derive it from the GTF. **Do NOT `rtracklayer::import` the GTF** (a full mouse/human GTF is multi-GB and risks OOM in a 32 G job). Extract it with a streaming awk pass over **transcript lines** (filtered nf-core GTFs have **no top-level `gene` feature lines** — keying on `$3=="gene"` silently yields an empty map):
   ```bash
   awk -F'\t' '$3=="transcript"{
     gid=gbt=gnm="";
     if (match($9,/gene_id "[^"]+"/))      gid=substr($9,RSTART+9, RLENGTH-10);
     if (match($9,/gene_biotype "[^"]+"/)) gbt=substr($9,RSTART+14,RLENGTH-15);
     if (match($9,/gene_name "[^"]+"/))    gnm=substr($9,RSTART+11,RLENGTH-12);
     if(!(gid in s)){s[gid]=1; print gid"\t"gbt"\t"gnm}
   }' {GTF_PATH} > {CWD}/ref/gene2typesym.tsv   # no header: gene_id, gene_biotype, gene_name
   ```
   Verify it is non-empty and has a sane biotype spread (`cut -f2 ref/gene2typesym.tsv | sort | uniq -c`) — expect ~20k `protein_coding` for mouse/human.

Store `{TX2GENE_TSV}` and `{GENE2TYPESYM_TSV}`. The GTF copy under nf-core `work/` is transient — stage the two small TSVs into `{CWD}/ref/` so the Rmd does not depend on `work/`. If no GTF is reachable at all, fall back to building both from the GTF, or proceed without the `gene_biotype` filter (warn the user that protein-coding filtering will be skipped).

---

## Step 6 — Sample sheet (auto-scaffold + metadata)

**Auto-generate a samples scaffold** from the detected `quant.sf` folders and the QC join, then have the user fill in experimental factors.

1. Build a draft `{WD_NAME}_samples.xlsx` with one row per detected sample:
   - `Sample` — the sample folder name (must match the multiqc `Sample` column exactly)
   - `Folder` — relative path to the folder containing `quant.sf`
   - `Condition` — **blank, user fills in** (the primary DE factor)
   - plus the curated QC columns from Step 7, pre-joined
2. Show the table. Ask the user to fill `Condition` (and any extra factors: `Cell`, `Experiment`, `Genotype`, `Sex`, etc.).
3. Validate: `Condition` values cannot contain spaces; warn if a project has only one condition (no DE possible) or fewer than 2 replicates per condition (DE underpowered).

Store as `{SAMPLE_META}`. If the user already has a hand-made samples xlsx (as in the example Rmds), offer to use that instead and just join QC onto it.

---

## Step 7 — QC metric join (curated)

Parse `{MULTIQC_STATS}` and join a **curated subset** of QC metrics onto `samples` (keyed on `Sample`), so they are available for colouring exploratory PCA. **Filter out the per-read sub-rows first** (`Sample` ending in ` Read 1` / ` Read 2`) — they have mostly empty columns and will corrupt the join.

Curated columns (rename to short, plottable names):

| multiqc column | samples column |
|---|---|
| `salmon-percent_mapped` | `PctMapped` |
| `star-uniquely_mapped_percent` | `UniqMapPct` |
| `picard_mark_duplicates-PERCENT_DUPLICATION` | `PctDup` |
| `sortmerna-rRNA_pct` | `rRNApct` |
| `qualimap_rnaseq-5_3_bias` | `Bias5to3` |
| `qualimap_rnaseq-reads_aligned` | `ReadsAligned` |

R snippet (goes in Rmd 01):
```r
qc <- read.delim("{MULTIQC_STATS}", check.names = FALSE)
qc <- qc[!grepl(" Read [12]$", qc$Sample), ]          # drop per-read sub-rows
qc <- qc %>% dplyr::transmute(
  Sample       = Sample,
  PctMapped    = `salmon-percent_mapped`,
  UniqMapPct   = `star-uniquely_mapped_percent`,
  PctDup       = `picard_mark_duplicates-PERCENT_DUPLICATION`,
  rRNApct      = `sortmerna-rRNA_pct`,
  Bias5to3     = `qualimap_rnaseq-5_3_bias`,
  ReadsAligned = `qualimap_rnaseq-reads_aligned`)
samples <- dplyr::left_join(samples, qc, by = "Sample")
stopifnot(!any(is.na(samples$PctMapped)))             # catches Sample-name mismatches early
```
The `stopifnot` is important: if the join key doesn't match (e.g. the samples `Sample` differs from the multiqc `Sample`), every QC column becomes `NA` silently. Fail loudly instead. Some nf-core columns are occasionally absent depending on the run — guard each `select` with `any(colnames(qc) == ...)` and drop missing ones with a warning rather than erroring.

These columns become PCA colour overlays in Step 11: `plotPCA(vsd, intgroup = "PctMapped")`, `"PctDup"`, `"rRNApct"`, etc. — a fast way to spot a PC driven by a technical artefact rather than biology.

---

## Step 8 — DE comparisons & downstream menu

Ask: "Which condition is the reference (baseline) level?" and "Which contrasts do you want?" (e.g. each treatment vs reference). For multi-factor designs, ask for the design formula (default `~ Condition`).

Ask which optional downstream modules to generate (multi-select):
1. **GSEA** — in-R `fgsea`, ranked by the DESeq2 Wald statistic (`stat`). Ask which **MSigDB collection(s)** to test (per-project biology choice; same menu as the seurat skill: H / C2:CP / C5 GO:BP / C7:IMMUNESIGDB / C8 / C6 / C3:TFT). Store as `{GSEA_COLLECTIONS}`.
2. **GseaVis enrichment plots** — `gseaNb` curves for named gene sets of interest.
3. **Counts-of-interest plots** — `plotCounts` jitter plots for a user-supplied gene list (e.g. sex genes `Xist`/`Ddx3y` for sex-check, or pathway genes).
4. **Heatmap** — `ComplexHeatmap` of top variable / top DE genes.

`msigdbr` + `fgsea` must be in the image for module 1; `GseaVis` for module 2 (verified in Step 2).

---

## Step 9 — SLURM resources

Bulk RNAseq is light compared with scRNA. Defaults:
1. Standard (≤ 24 samples) — 32G, 2h, 4 CPUs
2. Large (> 24 samples, or many contrasts/GSEA collections) — 64G, 4h, 8 CPUs
3. Custom

Never run R or Singularity on the login node.

---

## Step 10 — Project name / author

Ask for author name and project title for the Rmd headers.

---

## Step 11 — Generate Rmd 01: Import, Annotation & QC

Write `{TODAY_YYMMDD}_{WD_NAME}_01_import_qc.Rmd`. **Self-contained — do NOT `source('Rcode/*.R')`; inline every function used.** Set `options(scipen = 9)`.

Key sections:
1. **Load libraries** (mode-dependent): always `tximport, DESeq2, tidyverse, matrixStats, openxlsx, ggrepel`; add `edgeR` only for 3′ DGE.
2. **Load annotation** by reading the two staged TSVs (Step 5) — `tx2gene` from nf-core's `salmon.merged.tx2gene.tsv`, `gene2typesym` from the awk-extracted biotype map. Do NOT `rtracklayer::import` the GTF in the Rmd:
   ```r
   tx2gene <- read.delim(TX2GENE_TSV, header = TRUE) %>%
     dplyr::select(transcript_id, gene_id) %>% distinct() %>% filter(!is.na(transcript_id))
   gene2typesym <- read.delim(GENE2TYPESYM_TSV, header = FALSE,
     col.names = c("gene_id","gene_biotype","gene_name")) %>% distinct()
   ```
3. **Import samples + QC join** (Steps 6–7).
4. **tximport — MODE-BRANCHED (the critical fork):**
   ```r
   files <- file.path(dir, samples$Folder, "quant.sf")
   names(files) <- samples$Sample
   stopifnot(all(file.exists(files)))
   txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

   intCt <- round(txi$counts, 0)                 # ALWAYS round: quant.sf counts are not integers
   colnames(intCt) <- paste0(colnames(intCt), ".intCt")
   ```
   Then, **paired-end**:
   ```r
   l2tpm <- log2(txi$abundance + 1)              # TPM = length-normalised
   colnames(l2tpm) <- paste0(colnames(l2tpm), ".l2tpm")
   dumpDat <- merge(intCt, l2tpm, by = 0, all = TRUE)
   ```
   or **3′ DGE**:
   ```r
   y     <- edgeR::DGEList(txi$counts)           # NO length normalisation for 3' tag data
   cpms  <- edgeR::cpm(y)
   l2cpm <- log2(cpms + 1)
   colnames(l2cpm) <- paste0(colnames(l2cpm), ".l2cpm")
   dumpDat <- merge(intCt, l2cpm, by = 0, all = TRUE)
   ```
   Use the expression suffix (`.l2tpm` vs `.l2cpm`) consistently everywhere downstream.
5. **Annotate** `dumpDat` → `assembleDat`: merge `gene2typesym`; add `Avg = rowMeans(.l2tpm/.l2cpm)`, `Var = rowVars(...)`, `LowExpLowVar = (Avg <= 0.1 & Var <= 0.1)`.
6. **MPG (max-per-gene representative)** — loop over unique `gene_name`, mark the highest-`Avg` geneid `"Yes"` (first on ties). Used to deduplicate symbols for plotting/ranking.
7. **QC boxplots** — all genes, then expressed protein-coding only (`MPG == "Yes" & gene_biotype == "protein_coding" & LowExpLowVar == "No"`).
8. **Exploratory PCA** — `vsd <- vst(dds, blind = TRUE)`, `plotPCA(vsd, intgroup = "Condition", ntop = 500)`, plus one PCA per curated QC column (`PctMapped`, `PctDup`, `rRNApct`, `Bias5to3`, …) to catch technical drivers. (Build `dds` here per Step 12 so the same object is reused.) **Lay the panel out compactly** — put the `dds`/`vst` build in its own chunk, then a plotting chunk with `fig.width=7, fig.height=6, out.width="50%", fig.show="hold"` so the ~7 PCAs tile 2-per-row (larger panels, still compact) instead of stacking full-height with whitespace between them. `print()` each plot (including the first) inside the loop so `fig.show="hold"` collects them all. (Tune density via `out.width`: `"50%"` = 2/row, `"33%"` = 3/row smaller.)
9. Save a checkpoint: `saveRDS(list(txi = txi, samples = samples, assembleDat = assembleDat), file.path(OUT, "import_checkpoint.rds"))` and write `assembleDat` to xlsx.

---

## Step 12 — Generate Rmd 02: Differential Expression

Write `{TODAY_YYMMDD}_{WD_NAME}_02_deg.Rmd`. Load the checkpoint from Rmd 01.

Key sections:
1. **Set factors / reference level:** `samples$Condition <- factor(samples$Condition, levels = c("{REF}", ...))`.
2. **Build `dds` — MODE-BRANCHED (the second critical fork):**
   - **paired-end:**
     ```r
     dds <- DESeqDataSetFromTximport(txi = txi, colData = samples, design = ~ Condition)
     ```
     (passes `txi` so DESeq2 uses the average-transcript-length offset — correct for full-length libraries.)
   - **3′ DGE:**
     ```r
     dds <- DESeqDataSetFromMatrix(countData = intCt, colData = samples, design = ~ Condition)
     ```
     (uses the rounded integer counts directly, **no** length offset — correct for 3' tag libraries. `intCt` here must be `round(txi$counts, 0)` without the `.intCt` colname suffix.)
   - `dds$Condition <- relevel(dds$Condition, "{REF}")`.
3. **`dds <- DESeq(dds)`**; `resultsNames(dds)`.
4. **Per-contrast results** — `lfcShrink(dds, coef = i, type = "apeglm")`; recover the Wald `stat` from `dds@rowRanges@elementMetadata[, paste0("WaldStatistic_", i)]` (apeglm drops `stat`; the `stat`/Wald column is the GSEA ranking metric). Suffix columns per comparison and merge into `assembleDat`.
5. **Write outputs:** supplemental workbook `*_supplemental.xlsx` — a multi-sheet file bundling the per-gene `assembled` table (counts, l2tpm/l2cpm, DE stats, **plus a per-contrast `.sig` direction column** = Up/Down/ns by `LFC_SIG`/`PADJ_SIG`, **restricted to the same expressed-protein-coding-MPG universe** as the summary/volcano so all three report identical counts — genes outside that universe are `ns` by design) and the `DEG_summary` sheet (Step 5b). Plus GEO processed-data matrices (`.intCt`, `.l2tpm`/`.l2cpm`) as tab-delimited txt keyed by ENSID. Keep the `.sig`/summary cutoffs identical to the volcano.
5b. **DEG summary table** — per-contrast counts of Up / Down / Total over the expressed protein-coding MPG universe. **Define `LFC_SIG` / `PADJ_SIG` once as shared constants (default `1` / `0.05`) and use them for BOTH this summary AND the volcano colouring/guides** so the table and the figure can never report different significant-gene sets. Write `*_DEG_summary.xlsx`.
   ```r
   data.frame(Comparison = cn,
              Up   = sum(lfc >  LFC_SIG & padj < PADJ_SIG),
              Down = sum(lfc < -LFC_SIG & padj < PADJ_SIG),
              Total_DE = sum(abs(lfc) > LFC_SIG & padj < PADJ_SIG))
   ```
6. **Volcano plots** — filter to expressed protein-coding MPG genes; `NA padj → 1`, `NA log2FC → 0`, ceiling small padj (`< 1e-20 → 1e-20`); inline `ggplot`, one plot per contrast. **Conventional orientation: `log2FC` on x, `-log10(padj)` on y** (the `1e-20` floor caps the y-axis at `-log10(1e-20) = 20` — this is the auto axis cap, not a hardcoded `ylim`). **Colour by direction × significance**: Up = `LFC > 1 & padj < 0.05` (red), Down = `LFC < -1 & padj < 0.05` (blue), else grey; dashed guides at `x = ±1` and `y = -log10(0.05)`:
   ```r
   d$sig <- dplyr::case_when(d[[col_lfc]] >  1 & d[[col_padj]] < 0.05 ~ "Up",
                             d[[col_lfc]] < -1 & d[[col_padj]] < 0.05 ~ "Down",
                             TRUE ~ "ns")
   # ... aes(x = LFC, y = -log10(padj), color = sig) +
   #     scale_color_manual(c(Up="red", Down="blue", ns="grey70"), breaks=c("Up","Down"))
   ```
   **Label only the top N (default 20) most significant genes per contrast — rank by `abs(.stat)` (Wald), not a LFC/padj threshold.** A fixed threshold over-labels large DEG sets (thousands pass `abs(LFC)≥1 & padj<0.05` → unreadable wall of labels + `ggrepel: N unlabeled data points (too many overlaps)`). `abs(stat)` ranks cleanly with no ties at the padj floor:
   ```r
   label_idx <- head(order(abs(d[[paste0(i, ".stat")]]), decreasing = TRUE), N_LABEL)  # N_LABEL <- 20
   d$toLabel <- NA_character_; d$toLabel[label_idx] <- d$gene_name[label_idx]
   ```
   `geom_label_repel(color = "black")` on `toLabel` only; x-limits from `max(abs(LFC))`.
7. **Write `.rnk` files** for the GSEA module: `gene_name` + Wald `stat`, MPG protein-coding only, `complete.cases`, into `{RESULTS_DIR}/gsea/`.

---

## Step 13 — Generate Rmd 03: GSEA (if requested)

Write `{TODAY_YYMMDD}_{WD_NAME}_03_gsea.Rmd`. **Use in-R `fgsea` — self-contained, no external Java GSEA.**

1. Build ranked vector per contrast from the DESeq2 Wald `stat` (dedup by gene, drop NA, `arrange(desc(stat))`).
2. Gene sets via `msigdbr` over `{GSEA_COLLECTIONS}`. **Mouse: `msigdbr(db_species = "HS", species = "Mus musculus", collection = ..., subcollection = ...)`** then split `gene_symbol` by `gs_name` into a named list for fgsea. (Use current API: `collection`/`subcollection`, not deprecated `category`/`subcategory`.)
   - **Caveat — `msigdbr` 26.x fetches its gene-set data at runtime** (via the `msigdbdf` data package) instead of bundling it; the first `msigdbr()` call downloads tens of MB. The render therefore **needs internet on the compute node**. On an offline-compute cluster the GSEA step will hang — pre-seed the `msigdbdf` cache (run `msigdbr::msigdbr()` once on a node with internet, or bake `msigdbdf` into the image) before submitting `run_03`.
3. `fgsea(pathways = sets, stats = gene_ranks, minSize = 10, maxSize = 500)` per collection per contrast; `padj < 0.05`.
4. NES bar/dot plots (colour by sign); leading-edge genes per significant set to xlsx (one sheet per collection).
5. If GseaVis requested: `gseaNb()` enrichment curves for named sets of interest; **`fig.width >= 14`** or the p-value table is cropped.

---

## Step 14 — Generate SLURM scripts

One script per Rmd. Template:

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
SIF=/net/bmc-lab3/data/bcc/shared/singularity_images/bulkrnaseq_latest.sif

mkdir -p ${OUT}
module load singularity/3.5.0

singularity exec \
    --bind ${CWD}:${CWD} \
    ${SIF} \
    Rscript -e "rmarkdown::render(
        '${CWD}/{RMD_FILE}',
        output_dir    = '${OUT}',
        knit_root_dir = '${CWD}'
    )"
```

Scripts: `run_01_import_qc.sh`, `run_02_deg.sh`, and (if requested) `run_03_gsea.sh`.

---

## Step 15 — Generate summary report

After the reports are rendered, write a standalone `{WD_NAME}_summary_report.html` into `{RESULTS_DIR}` — a self-contained page (inline CSS, **no external dependencies**, no R needed) that ties the deliverables together:
- Study-design block: mode (paired-end/3pDGE), samples/conditions, organism/annotation, contrast(s), significance cutoff.
- QC summary table (the curated multiqc metrics per sample).
- **Card links to each stage HTML** (`_01_import_qc.html`, `_02_deg.html`, `_03_gsea.html`) and to the data deliverables (supplemental/DEG-summary/GSEA xlsx, GEO matrices).
- Headline results: DEG up/down counts and top GSEA programs.

**Use relative links** (bare filenames) so the page works both over a web mount and from the local filesystem as long as it sits in `{RESULTS_DIR}` beside the reports — no need to know a web URL prefix. If the user has a web-accessible prefix (e.g. `https://bmc-data.mit.edu/BCC/...`) and wants absolute links, ask for it; otherwise default to relative. Verify every `href` target exists before finishing.

**Two report flavours — pick by audience:**
- *Summary index* (above) — a lightweight page that **links out** to the three stage HTMLs. For the analyst / internal use.
- *Consolidated client report* — **one self-contained HTML that embeds the actual results** (figures + tables from all three stages) for handing to the scientist. Build it as a separate Rmd (`{...}_report.Rmd`) that **loads the saved outputs** (`import_checkpoint.rds`, `*_supplemental.xlsx`, `*_GSEA_fgsea.xlsx`) and re-plots — do NOT re-run DESeq. Use `echo=FALSE` globally (hide all code), narrative section intros in plain language (define "up = higher in <test>"), `output: html_document` with `toc_float`, `theme: flatly`, `df_print: kable`, `self_contained: true`. Render with a clean `output_file` like `{PROJECT}_Results_Report.html`. Verify it embeds figures as base64 with zero external `src` (truly portable single file).

## Final submission instructions

```
Submission order:
  sbatch -p bcc run_01_import_qc.sh
  # Review 01_import_qc.html — check boxplots, PCA (biology vs QC overlays), sample grouping
  sbatch -p bcc run_02_deg.sh
  # Review volcano plots and *_assembled.xlsx
  sbatch -p bcc run_03_gsea.sh     # if requested
```

---

## Notes for the assistant

### Library mode (paired-end vs 3′ DGE) — the central distinction
- **Decide and confirm the mode before generating any Rmd** (Step 4). It changes the import normalisation AND the DESeq2 constructor.
- **Paired-end:** `tximport` → `txi$abundance` is TPM (length-normalised). Expression = `l2tpm = log2(TPM + 1)`. DE via `DESeqDataSetFromTximport(txi, ...)`, which applies the average-transcript-length offset. This is correct because full-length library coverage scales with transcript length.
- **3′ DGE:** counts come from a single 3' tag per molecule, so count is **independent of transcript length** — length normalisation is wrong. Do **not** use TPM and do **not** pass `txi` to DESeq2. Expression = `l2cpm` via edgeR `DGEList` → `cpm`. DE via `DESeqDataSetFromMatrix(countData = round(txi$counts), ...)` on integer counts with no offset.
- The `tximport()` call itself is identical for both modes; the fork is entirely downstream (units + DESeq2 constructor + whether edgeR is used).

### Integer counts
- `quant.sf` / `txi$counts` are estimated counts and are **not integers**. DESeq2 requires integer `countData`. Always `round(txi$counts, 0)`.
- For paired-end, `DESeqDataSetFromTximport` handles this internally from `txi` (do not pre-round the object you pass it) — but the `intCt` table for the dump/GEO files still needs rounding.
- For 3′ DGE, you pass `round(txi$counts, 0)` directly as `countData`.

### tximport / annotation
- Salmon quantifies transcripts; `tx2gene` (transcript_id → gene_id) is mandatory for gene-level summarisation.
- **Prefer nf-core's `star_salmon/salmon.merged.tx2gene.tsv` for `tx2gene`** — it is exact-matched to the quantified transcripts (no transcript-version mismatch with `quant.sf`) and avoids parsing the GTF. Stage it to `{CWD}/ref/`.
- **`gene_biotype` is NOT in the nf-core tx2gene** — derive `gene2typesym` from the GTF. **Do not `rtracklayer::import` a full GTF** (multi-GB → OOM risk in a 32 G job); use the streaming awk extract in Step 5.
- **Filtered nf-core GTFs have no top-level `gene` feature lines** — extract biotype from `$3=="transcript"` lines and dedup by `gene_id`. Keying on `$3=="gene"` silently produces an empty map (the failure mode that motivated this approach).
- Stage both as small TSVs in `{CWD}/ref/`; the GTF under nf-core `work/` is transient and gets cleaned.
- The annotation may map several gene_ids to one gene_name. The **MPG** ("max-per-gene") flag selects one representative geneid per symbol (highest mean expression) so plots/rankings aren't duplicated. First-on-ties.
- `LowExpLowVar` (Avg ≤ 0.1 & Var ≤ 0.1) + `gene_biotype == "protein_coding"` + `MPG == "Yes"` is the standard "expressed, informative" filter for boxplots, volcanoes, and rnk files.

### QC join from multiqc
- Source: `multiqc/star_salmon/multiqc_report_data/multiqc_general_stats.txt` (60 columns).
- **Filter out `Sample` rows ending in ` Read 1` / ` Read 2`** before joining — they are per-read FastQC sub-rows with mostly empty columns and will pollute the merge.
- Join key is the multiqc `Sample` column = the salmon sample folder name. If `samples$Sample` differs, the join silently fills `NA`; assert non-NA (`stopifnot`) to fail loudly.
- Curated columns: `PctMapped, UniqMapPct, PctDup, rRNApct, Bias5to3, ReadsAligned`. Guard each against absence (`any(colnames(qc) == ...)`) since not every nf-core run emits all modules.
- Purpose: colour exploratory PCA by each QC metric (`plotPCA(vsd, intgroup = "PctDup")` etc.) to distinguish biological from technical PCs.

### apeglm / Wald statistic
- `lfcShrink(type = "apeglm")` needs `coef = <name from resultsNames(dds)>` (not a `contrast`), and it **drops the `stat` column**. Recover the Wald statistic from `dds@rowRanges@elementMetadata[, paste0("WaldStatistic_", coef)]` and keep it — it is the correct GSEA ranking metric (captures effect size + significance), not log2FC.

### GSEA
- Prefer in-R `fgsea` ranked by Wald `stat` over the external Java preranked GSEA — self-contained and reproducible inside the container.
- Collections are a **per-project biology choice** — always ask, never hardcode Hallmark.
- **msigdbr current API:** `collection=` / `subcollection=` (not `category`/`subcategory`). **Mouse:** `db_species = "HS", species = "Mus musculus"` — msigdbr ortholog-maps the human DB to mouse symbols (the mouse-native DB is typically not installed). Columns: `gs_name`, `gene_symbol`.
- **msigdbr 26.x downloads gene-set data at runtime** (the `msigdbdf` data package), so the GSEA render needs internet on the compute node. On offline-compute clusters, pre-seed the `msigdbdf` cache or bake it into the image before submitting `run_03`, or it will hang on the first `msigdbr()` call.
- `GseaVis::gseaNb` chunks need `fig.width >= 14` or the p-value table is cropped.

### Self-contained Rmds
- Do **not** `source('Rcode/data_formatting_tools.R')` etc. — those IGB helpers are not assumed present. Inline only the functions actually used (MPG loop, volcano builder). The example Rmds source them but barely call them; everything needed is standard `tximport`/`DESeq2`/`edgeR`/`ggplot`.

### Figures
- **Exploratory PCA panel**: build `dds`/`vst` in one chunk, plot in a second chunk with `fig.width=5, fig.height=4, out.width="33%", fig.show="hold"` so the Condition + per-QC-metric PCAs tile 3-per-row. Stacking them at full `fig.height` leaves large vertical gaps between plots.
- Volcano: conventional orientation (`log2FC` on x, `-log10(padj)` on y); `NA padj → 1`, `NA log2FC → 0`, ceiling tiny padj (`< 1e-20 → 1e-20`, which auto-caps the y-axis at 20 — not a manual `ylim`). Colour Up=`LFC>1 & padj<0.05` red / Down=`LFC<-1 & padj<0.05` blue / ns grey, dashed guides at `±1` and `-log10(0.05)`. **Label only the top N≈20 genes by `abs(.stat)`** (Wald), never a LFC/padj threshold — thresholds over-label large DEG sets into an unreadable mass.
- Unicode labels (e.g. Δ for deletion genotypes) need `library(showtext); showtext_auto()` before `ggsave`.
- `ggsave(..., dpi = 300)` for figures destined for slides/supplements.

### SLURM / environment
- Always submit with `sbatch -p bcc script.sh` — the `#SBATCH -p bcc` header alone silently fails on this cluster.
- `module load singularity/3.5.0` (not `module add`). Never run R or Singularity on the login node.
- Packages live in the `yannvrb56/bulkrnaseq` Docker image — never `install.packages()` at runtime; rebuild the image to add packages.
- Raw `quant.sf` and nf-core outputs are read-only. Write all outputs to `{RESULTS_DIR}`.

### General
- **Namespace masking:** in any Rmd mixing Bioconductor + tidyverse, **load `DESeq2`/Bioc BEFORE `tidyverse`** so dplyr verbs win. `S4Vectors` exports `rename`/`first`/`setdiff`/`union`/`intersect`; loading it after dplyr shadows `dplyr::rename` and a bare `rename(new = old)` then errors with `object 'old' not found`. (The stage Rmds load tidyverse first but only call `dplyr::`-prefixed verbs; the consolidated report calls bare verbs, so order matters there.)
- **Set `knitr::opts_chunk$set(cache = FALSE)`** in every generated Rmd. These Rmds mutate a shared `assembleDat` across chunks (DE merge → `.sig` annotation → supplemental write); with `cache = TRUE`, knitr serves a stale cached copy and downstream sheets silently disagree with the live data (observed: `.sig` counts not matching the DEG summary until the cache was cleared). The pipeline is fast (tximport/DESeq are seconds) so caching buys nothing here. If a cached Rmd ever shows stale numbers, delete its `*_cache/` dir and re-render.
- Figure layout density is tunable via `out.width` (PCA panel: `"50%"` = 2/row default).
- `options(scipen = 9)` so numbers display normally.
- R line continuation: never `\`; lines continue when ending on an open operator (`+`, `,`, `|>`, `(`).
- Keep the expression-unit suffix (`.l2tpm` vs `.l2cpm`) consistent with the chosen mode across all Rmds.
