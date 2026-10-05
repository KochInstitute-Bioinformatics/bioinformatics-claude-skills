# Usage (inside the R image, on a compute node only; tests/r/run_r_tests.sh does this through sbatch):
#   Rscript test_usage_cross.R <skill.md> <tests/fixtures dir>
# Extracts the Rmd blocks under **File `04_usage.Rmd`.** and **File `05_cross_tool.Rmd`.** of the skill, evaluates their `functions`
# chunks, then renders 01_import.Rmd, 04_usage.Rmd and 05_cross_tool.Rmd on the fixture tree and compares the counts, the gene sets,
# the intersections and the directions with an independent recomputation from the raw tables (base R, no code of the skill).
# The data frames SYNTHETIC_* hold three tools with hand-computed intersections and directions: they test the counting logic only,
# never the signs of a real tool. The sign of the real tables is recomputed here from the rules (DEXSeq and SUPPA2 flipped).
# Prints one "ok" or "FAIL:" line per check, then "R TESTS PASS (<n>)" or "R TESTS FAIL (<failed>/<n>)" (exit 1).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: test_usage_cross.R <skill.md> <fixtures dir>")
skill <- args[1]; fx <- args[2]; rt <- file.path(fx, "real_testdata")
n_ok <- 0; n_fail <- 0
check <- function(name, cond) {
  r <- tryCatch({ do.call(stopifnot, stats::setNames(list(isTRUE(all(cond))), name)); TRUE },
                error = function(e) { cat("FAIL:", conditionMessage(e), "\n"); FALSE })
  if (r) { cat("ok  ", name, "\n"); n_ok <<- n_ok + 1 } else n_fail <<- n_fail + 1
  invisible(r)
}
finish <- function() {
  if (n_fail == 0) { cat(sprintf("R TESTS PASS (%d)\n", n_ok)); quit(status = 0) }
  cat(sprintf("R TESTS FAIL (%d/%d)\n", n_fail, n_ok + n_fail)); quit(status = 1)
}
work <- tempfile("rsr_usage_"); dir.create(work)
txt <- readLines(skill, warn = FALSE)

# --- extract an Rmd block by its file name
get_rmd <- function(name) {
  h <- which(txt == sprintf("**File `%s`.**", name))
  if (!check(sprintf("the skill has the heading **File `%s`.** once", name), length(h) == 1)) finish()
  o <- h + which(txt[(h + 1):length(txt)] == "````rmd")[1]
  e <- o + which(txt[(o + 1):length(txt)] == "````")[1]
  if (!check(sprintf("a ````rmd block follows the heading of %s", name), !is.na(o) && !is.na(e) && e > o + 1)) finish()
  txt[(o + 1):(e - 1)]
}
functions_env <- function(rmd, name) {
  fb <- grep("^```\\{r functions", rmd)
  if (!check(sprintf("%s has one functions chunk", name), length(fb) == 1)) finish()
  fe <- fb + which(rmd[(fb + 1):length(rmd)] == "```")[1]
  e <- new.env(parent = globalenv())
  err <- tryCatch({ eval(parse(text = rmd[(fb + 1):(fe - 1)]), envir = e); "" }, error = function(x) conditionMessage(x))
  if (!check(sprintf("the functions chunk of %s evaluates (%s)", name, err), err == "")) finish()
  e
}
rmd4 <- get_rmd("04_usage.Rmd"); rmd5 <- get_rmd("05_cross_tool.Rmd")
E4 <- functions_env(rmd4, "04_usage.Rmd"); E5 <- functions_env(rmd5, "05_cross_tool.Rmd")
need_fun <- function(e, nm) all(vapply(nm, exists, logical(1), envir = e, inherits = FALSE))
if (!check("04 defines its functions", need_fun(E4, c("sig_usage", "none_message", "usage_wording", "top_by_gene", "nearest_rows", "usage_counts",
    "bin_number", "stager_alpha", "read_gene_table", "read_stager", "read_dtu_groups", "volcano_usage", "exon_usage_plot", "dtu_group_plot")))) finish()
if (!check("05 defines its functions", need_fun(E5, c("strip_version", "expand_genes", "top_by_gene", "calls_from_events", "calls_from_gene_table",
    "calls_from_stager", "membership_table", "intersection_counts", "cross_table", "concordance_table", "upset_engine", "make_upset", "TOOL_LABELS")))) finish()
txt_of <- function(rmd) paste(rmd, collapse = "\n")
SIGN <- "Every effect is stated as **treatment minus control**"
check("04 states the sign rule before its first function", regexpr(SIGN, txt_of(rmd4), fixed = TRUE) < regexpr("```{r functions", txt_of(rmd4), fixed = TRUE) &&
      regexpr(SIGN, txt_of(rmd4), fixed = TRUE) > 0)
check("05 states the sign rule before its first function", regexpr(SIGN, txt_of(rmd5), fixed = TRUE) < regexpr("```{r functions", txt_of(rmd5), fixed = TRUE) &&
      regexpr(SIGN, txt_of(rmd5), fixed = TRUE) > 0)
check("05 states that directions only are compared, never effect sizes across scales",
      grepl("directions compared only, never effect sizes across scales", txt_of(rmd5), fixed = TRUE))
check("05 states the background of the comparison", grepl("Background:", txt_of(rmd5), fixed = TRUE))
check("04 names the stageR columns gene and transcript", grepl("`gene`", txt_of(rmd4), fixed = TRUE) && grepl("`transcript`", txt_of(rmd4), fixed = TRUE) && grepl("stageR", txt_of(rmd4), fixed = TRUE))

# --- 1. SYNTHETIC, Rmd 04: counting logic
SYNTHETIC_feat <- data.frame(
  gene_id = c("G1", "G1", "G2", "G2", "G3", "G4"), feature_id = c("G1:E001", "G1:E002", "G2:E001", "G2:E002", "G3:E001", "G4:E001"),
  padj = c(0.04, 0.001, 0.02, 0.2, NA, 0.5), effect = c(2, -3, 1, 4, NA, 1),
  direction = c("treatment_up", "control_up", "treatment_up", "treatment_up", "none", "treatment_up"), stringsAsFactors = FALSE)
check("sig_usage: padj strictly below the FDR, NA never significant", identical(E4$sig_usage(SYNTHETIC_feat, 0.05), c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)) &&
      identical(E4$sig_usage(data.frame(padj = c(0.05, 0.049)), 0.05), c(FALSE, TRUE)))
tg <- E4$top_by_gene(SYNTHETIC_feat, E4$sig_usage(SYNTHETIC_feat, 0.05))
check("top_by_gene: one row per gene with a significant feature (G1, G2)", identical(sort(tg$gene_id), c("G1", "G2")))
check("top_by_gene: G1 takes the feature with the SMALLEST padj (E002, control_up), not the first or the largest effect",
      identical(tg$top_feature[tg$gene_id == "G1"], "G1:E002") && identical(tg$direction[tg$gene_id == "G1"], "control_up") && tg$n_sig_features[tg$gene_id == "G1"] == 2L)
check("top_by_gene: G2 has one significant feature (E001, treatment_up)", identical(tg$top_feature[tg$gene_id == "G2"], "G2:E001") && identical(tg$direction[tg$gene_id == "G2"], "treatment_up"))
check("top_by_gene: a zero effect gives direction NA", is.na(E4$top_by_gene(data.frame(gene_id = "Z", feature_id = "Z:E1", padj = 0.01, effect = 0, direction = "none", stringsAsFactors = FALSE), TRUE)$direction))
check("usage_wording: positive is more used in treatment, negative more used in control",
      identical(E4$usage_wording("treatment_up"), "more used in treatment") && identical(E4$usage_wording("control_up"), "more used in control"))
check("none_message names the FDR", identical(E4$none_message(0.05), "none at FDR 0.05") && identical(E4$none_message(0.5), "none at FDR 0.5"))
uc <- E4$usage_counts("exon_dexseq", "exon_bin", SYNTHETIC_feat$padj, E4$sig_usage(SYNTHETIC_feat, 0.05), SYNTHETIC_feat$direction)
check("usage_counts: 5 tested, 3 significant, 2 more used in treatment, 1 more used in control",
      identical(as.integer(uc[1, c("tested", "significant", "more used in treatment", "more used in control")]), c(5L, 3L, 2L, 1L)))
check("usage_counts without directions leaves the direction columns NA", all(is.na(E4$usage_counts("exon_dexseq", "gene", c(0.01, 0.9), c(TRUE, FALSE))[1, c("more used in treatment", "more used in control")])))
check("nearest_rows gives the three smallest padj", identical(E4$nearest_rows(data.frame(gene_id = letters[1:5], padj = c(0.9, 0.1, NA, 0.2, 0.05), stringsAsFactors = FALSE), 3)$gene_id, c("e", "b", "d")))
check("bin_number reads the exon number after the last colon (merged gene groups included)",
      identical(E4$bin_number(c("ENSG1:E001", "ENSG1+ENSG2:E024", "ENSG3:E100")), c(1L, 24L, 100L)))
# stageR table: columns geneID txID gene transcript; NA where not screened
stf <- file.path(work, "stager.tsv")
writeLines(c("geneID\ttxID\tgene\ttranscript", "GA\tTA1\t0.01\t0.02", "GA\tTA2\t0.01\t1", "GB\tTB1\t0.5\tNA", "GC\tTC1\t0.03\t0.04"), stf)
st <- E4$read_stager(stf)
check("read_stager reads gene and transcript as the screened adjusted p-values, NA where not screened",
      identical(st$gene, c(0.01, 0.01, 0.5, 0.03)) && identical(st$transcript, c(0.02, 1, NA, 0.04)) && identical(st$tx_id, c("TA1", "TA2", "TB1", "TC1")))
check("stager_alpha reads the overall alpha of the real stageRTx object (0.05)", isTRUE(all.equal(E4$stager_alpha(file.path(rt, "salmon/dexseq_dtu/results/stager/stageRTx.GBR-YRI.rds")), 0.05)))
check("stager_alpha is NA for a missing file", is.na(E4$stager_alpha(file.path(work, "absent.rds"))))
bad_stf <- file.path(work, "stager_bad.tsv"); writeLines(c("geneID\ttxID\tgene", "GA\tTA1\t0.01"), bad_stf)
check("read_stager stops naming the file and the missing column transcript",
      grepl("stager_bad.tsv", tryCatch({ E4$read_stager(bad_stf); "" }, error = function(e) conditionMessage(e)), fixed = TRUE) &&
      grepl("transcript", tryCatch({ E4$read_stager(bad_stf); "" }, error = function(e) conditionMessage(e)), fixed = TRUE))
# DTU group columns are named by the condition names
dtf <- file.path(work, "dtu.tsv")
writeLines(c("groupID\tfeatureID\tpadj\tGBR-1\tYRI_2", "GA\tTA1\t0.01\t1.5\t2.5", "GA\tTA2\t0.02\t3\t4"), dtf)
dg <- E4$read_dtu_groups(dtf, "GBR-1", "YRI_2")
check("read_dtu_groups takes the treatment and control values from the columns named by the conditions (GBR-1, YRI_2)",
      identical(dg$treatment, c(1.5, 3)) && identical(dg$control, c(2.5, 4)) && identical(dg$feature_id, c("TA1", "TA2")))
check("read_dtu_groups stops naming the file and the missing condition column",
      grepl("column not found: GBR", tryCatch({ E4$read_dtu_groups(dtf, "GBR", "YRI_2"); "" }, error = function(e) conditionMessage(e)), fixed = TRUE))

# --- 2. SYNTHETIC, Rmd 05: three tools, intersections and directions computed by hand
SYNTHETIC_rmats <- data.frame(
  contrast = "SYN", tool = "rmats", gene_id = c("G1", "G1", "G2", "G3", "G4", "G5", "G6", "G11"), gene_symbol = NA_character_, feature_type = "SE",
  feature_id = paste0("SE:", 1:8), effect = c(0.3, -0.5, -0.2, 0.4, 0.05, 0.3, 0.2, 0.35), effect_scale = "dPSI",
  padj = c(0.01, 0.04, 0.02, 0.2, 0.001, 0.001, 0.03, 0.01), direction = c("treatment_up", "control_up", "control_up", "treatment_up", "treatment_up", "treatment_up", "treatment_up", "treatment_up"),
  n_support = 100, low_support = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE), stringsAsFactors = FALSE)
# significant at FDR 0.05, |dPSI| 0.1: SE:1 (G1), SE:2 (G1), SE:3 (G2), SE:7 (G6); G3 fails FDR, G4 fails dPSI, G5 is low support
SYNTHETIC_suppa <- data.frame(
  contrast = "SYN", tool = "suppa", gene_id = c("G1", "G2", "G3", "G7", "G11"), gene_symbol = NA_character_, feature_type = "SE",
  feature_id = paste0("SE:", 11:15), effect = c(0.25, 0.3, -0.2, -0.4, 0.4), effect_scale = "dPSI", padj = c(0.02, 0.01, 0.01, 0.01, 0.02),
  direction = c("treatment_up", "treatment_up", "control_up", "control_up", "treatment_up"), n_support = NA_real_, low_support = FALSE, stringsAsFactors = FALSE)
SYNTHETIC_dex_feat <- data.frame(
  contrast = "SYN", tool = "exon_dexseq", gene_id = c("G1", "G3", "G3", "G6", "G9+G10"), gene_symbol = NA_character_, feature_type = "exon_bin",
  feature_id = c("G1:E001", "G3:E002", "G3:E003", "G6:E001", "G9+G10:E001"), effect = c(2, 3, -4, -1, 5), effect_scale = "log2FC",
  padj = c(0.001, 0.01, 0.04, 0.5, 0.002), direction = c("treatment_up", "treatment_up", "control_up", "control_up", "treatment_up"),
  n_support = NA_real_, low_support = FALSE, stringsAsFactors = FALSE)
SYNTHETIC_dex_gene <- data.frame(gene_id = c("G1", "G2", "G3", "G6", "G8", "G9+G10"), padj = c(0.01, 0.2, 0.03, 0.01, 0.04, 0.02), stringsAsFactors = FALSE)
cr <- E5$calls_from_events(SYNTHETIC_rmats, 0.05, 0.1)
cs <- E5$calls_from_events(SYNTHETIC_suppa, 0.05, 0.1)
cd <- E5$calls_from_gene_table(SYNTHETIC_dex_gene, SYNTHETIC_dex_feat, 0.05)
check("SYNTHETIC rMATS set is G1, G2, G6, G11 (not G3: FDR, not G4: dPSI, not G5: low support; no non-significant gene)", identical(sort(cr$genes), sort(c("G1", "G2", "G6", "G11"))))
check("SYNTHETIC SUPPA2 set is G1, G2, G3, G7, G11", identical(sort(cs$genes), sort(c("G1", "G2", "G3", "G7", "G11"))))
check("SYNTHETIC DEXSeq set is the genes of the per-gene table below the FDR (G1, G3, G6, G8, G9, G10): the exon padj of G6 (0.5) and G8 (no bin) do not matter, G9+G10 is split",
      identical(sort(cd$genes), c("G1", "G10", "G3", "G6", "G8", "G9")))
dirof <- function(cl, g) { i <- match(g, cl$dir$gene_id); if (is.na(i)) NA_character_ else cl$dir$direction[i] }
check("SYNTHETIC rMATS directions: G1 takes SE:1 (smallest padj, treatment_up), G2 control_up, G6 treatment_up",
      identical(dirof(cr, "G1"), "treatment_up") && identical(dirof(cr, "G2"), "control_up") && identical(dirof(cr, "G6"), "treatment_up"))
check("SYNTHETIC DEXSeq directions: G1 treatment_up, G3 takes E002 (smallest padj among the significant bins: treatment_up), G6 and G8 have no significant bin (NA), G9 and G10 inherit the bin",
      identical(dirof(cd, "G1"), "treatment_up") && identical(dirof(cd, "G3"), "treatment_up") && is.na(dirof(cd, "G6")) && is.na(dirof(cd, "G8")) &&
      identical(dirof(cd, "G9"), "treatment_up") && identical(dirof(cd, "G10"), "treatment_up"))
calls_syn <- list(rmats = cr, suppa = cs, exon_dexseq = cd)
ic <- E5$intersection_counts(calls_syn)
icv <- stats::setNames(ic$n_genes, ic$combination)
check("SYNTHETIC UpSet intersections equal the hand-computed ones (exact membership: rmats&suppa&exon_dexseq 1 (G1), rmats&suppa 2 (G2, G11), suppa&exon_dexseq 1 (G3), rmats&exon_dexseq 1 (G6), suppa only 1 (G7), exon_dexseq only 3 (G8, G9, G10), rmats only none)",
      identical(icv[["rmats&suppa&exon_dexseq"]], 1L) && identical(icv[["rmats&suppa"]], 2L) && identical(icv[["suppa&exon_dexseq"]], 1L) &&
      identical(icv[["rmats&exon_dexseq"]], 1L) && identical(icv[["suppa"]], 1L) && identical(icv[["exon_dexseq"]], 3L) && !"rmats" %in% names(icv) && sum(icv) == 9L)
mt <- E5$membership_table(calls_syn)
check("SYNTHETIC membership table: one row per gene of the union (9), one logical column per tool", nrow(mt) == 9L && all(c("gene_id", "rmats", "suppa", "exon_dexseq") %in% names(mt)) &&
      identical(mt$exon_dexseq[mt$gene_id == "G6"], TRUE) && identical(mt$suppa[mt$gene_id == "G6"], FALSE))
xt <- E5$cross_table(calls_syn)
check("SYNTHETIC cross table: exactly the genes called by two or more tools (G1, G2, G3, G6, G11)", setequal(xt$gene_id, c("G1", "G2", "G3", "G6", "G11")) && nrow(xt) == 5L && all(xt$n_tools >= 2))
xd <- function(g, tool) xt[[paste0("direction_", tool)]][xt$gene_id == g]
check("SYNTHETIC cross table: G2 is control_up in rMATS and treatment_up in SUPPA2, G6 has direction NA in DEXSeq and none for SUPPA2 (not called)",
      identical(xd("G2", "rmats"), "control_up") && identical(xd("G2", "suppa"), "treatment_up") && is.na(xd("G6", "exon_dexseq")) && is.na(xd("G6", "suppa")))
cc <- E5$concordance_table(calls_syn)
pair <- function(a, b) cc[cc$tool_a == a & cc$tool_b == b, ]
check("SYNTHETIC concordance rmats-suppa: 3 genes in both, G1 and G11 same direction, G2 opposite directions (counted discordant): asymmetric, so a swap of same and opposite shows",
      identical(as.integer(pair("rmats", "suppa")[, c("genes_in_both", "same_direction", "opposite_direction", "no_direction")]), c(3L, 2L, 1L, 0L)))
check("SYNTHETIC concordance rmats-exon_dexseq: G1 same direction, G6 has no direction in DEXSeq (not counted as either)",
      identical(as.integer(pair("rmats", "exon_dexseq")[, c("genes_in_both", "same_direction", "opposite_direction", "no_direction")]), c(2L, 1L, 0L, 1L)))
check("SYNTHETIC concordance suppa-exon_dexseq: G1 same direction, G3 opposite (SUPPA2 control_up, DEXSeq treatment_up)",
      identical(as.integer(pair("suppa", "exon_dexseq")[, c("genes_in_both", "same_direction", "opposite_direction", "no_direction")]), c(2L, 1L, 1L, 0L)))
# directions are compared across scales: dPSI of 0.2 and log2FC of 3 in the same direction are concordant (a comparison of the numbers would not be)
sx_a <- list(genes = "GX", tested = "GX", dir = data.frame(gene_id = "GX", n_sig_features = 1L, top_feature = "a", top_effect = 0.2, top_padj = 0.01, direction = "treatment_up", stringsAsFactors = FALSE))
sx_b <- list(genes = "GX", tested = "GX", dir = data.frame(gene_id = "GX", n_sig_features = 1L, top_feature = "b", top_effect = 3, top_padj = 0.01, direction = "treatment_up", stringsAsFactors = FALSE))
cc2 <- E5$concordance_table(list(rmats = sx_a, exon_dexseq = sx_b))
check("a dPSI of 0.2 and a log2 fold change of 3 with the same sign are concordant (directions only, never effect sizes across scales)",
      identical(as.integer(cc2[, c("same_direction", "opposite_direction")]), c(1L, 0L)))
check("a pair of tools with no gene in common gives a row with zero counts", {
  z <- E5$concordance_table(list(rmats = sx_a, suppa = list(genes = "GY", tested = "GY", dir = sx_a$dir[0, ]))); identical(as.integer(z[, c("genes_in_both", "same_direction", "opposite_direction")]), c(0L, 0L, 0L)) })
check("expand_genes splits merged DEXSeq gene groups and strip_version removes a version suffix only when present",
      identical(E5$expand_genes(data.frame(gene_id = c("A+B", "C"), v = 1:2, stringsAsFactors = FALSE))$gene_id, c("A", "B", "C")) &&
      identical(E5$strip_version(c("ENSG00000147403.5", "ENSG00000147403", "ENST00000406022.2")), c("ENSG00000147403", "ENSG00000147403", "ENST00000406022.2")))
check("only one tool: the cross table has no gene and the concordance no row",
      nrow(E5$cross_table(list(rmats = cr))) == 0 && nrow(E5$concordance_table(list(rmats = cr))) == 0)
# the UpSet plot: both engines drawn, the one the image has by default
check("upset_engine is ComplexUpset or UpSetR (both are in the image)", E5$upset_engine() %in% c("ComplexUpset", "UpSetR"))
for (eng in c("ComplexUpset", "UpSetR")) {
  f <- file.path(work, paste0("upset_", eng, ".png"))
  drawn <- tryCatch(E5$make_upset(calls_syn, f, engine = eng), error = function(e) { cat("upset error:", conditionMessage(e), "\n"); FALSE })
  check(sprintf("make_upset draws the SYNTHETIC sets with %s", eng), isTRUE(drawn) && file.exists(f) && file.size(f) > 1000)
}
check("make_upset draws nothing when fewer than two tools have a gene", identical(E5$make_upset(list(rmats = cr, suppa = list(genes = character(0), tested = character(0), dir = cs$dir[0, ])), file.path(work, "none.png")), FALSE) &&
      !file.exists(file.path(work, "none.png")))

# --- 3. render 01, 04 and 05 on the fixture tree
base_params <- list(outdir = rt, report_dir = "", contrasts = file.path(rt, "261003_a_contrasts.csv"),
                    samplesheet = file.path(rt, "261003_a_samplesheet.csv"), modules = "rmats suppa dexseq_exon edger_exon dtu", rmats_type = "JC",
                    fdr = 0.05, dpsi = 0.1, min_reads = 10, events = TRUE, cross = TRUE, topgene = TRUE, enrich = TRUE,
                    organism = "human", gtf = "", enrich_dbs = "GO", paired_rmats = FALSE, paired_suppa = TRUE)
render <- function(rd, rmd_name, ...) {
  pr <- utils::modifyList(base_params, c(list(report_dir = rd), list(...)))
  tryCatch({ rmarkdown::render(file.path(rd, rmd_name), output_dir = rd, knit_root_dir = rd, params = pr, quiet = TRUE, envir = new.env()); TRUE },
           error = function(e) { cat("render error:", conditionMessage(e), "\n"); FALSE })
}
rmd1 <- get_rmd("01_import.Rmd")
new_dir <- function(nm) { d <- file.path(work, nm); dir.create(d); writeLines(rmd1, file.path(d, "01_import.Rmd"))
  writeLines(rmd4, file.path(d, "04_usage.Rmd")); writeLines(rmd5, file.path(d, "05_cross_tool.Rmd")); d }
html_of <- function(rd, nm) gsub("[[:space:]]+", " ", paste(readLines(file.path(rd, nm), warn = FALSE), collapse = " "))
has <- function(html, s) grepl(s, html, fixed = TRUE)
nonempty <- function(f) file.exists(f) && file.size(f) > 0
num <- function(x) suppressWarnings(as.numeric(ifelse(x %in% c("nan", "NA", ""), NA, x)))
C <- "GBR_vs_YRI"
rd_csv <- function(f) utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
rd_tsv <- function(f) utils::read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
dir_txt <- function(x) ifelse(is.na(x) | x == 0, NA_character_, ifelse(x > 0, "treatment_up", "control_up"))

# --- raw recomputation (base R from the raw fixture tables; the sign rules are applied here: DEXSeq and SUPPA2 flipped)
raw_usage <- function(fdr) {
  dex <- rd_csv(file.path(rt, "star/dexseq_exon/results/DEXSeqResults.GBR_vs_YRI.csv")); pg <- rd_csv(file.path(rt, "star/dexseq_exon/results/perGeneQValue.GBR_vs_YRI.csv"))
  eex <- rd_csv(file.path(rt, "star/edger/contrast_GBR_vs_YRI.usage.exon.csv")); egn <- rd_csv(file.path(rt, "star/edger/contrast_GBR_vs_YRI.usage.gene.csv"))
  esm <- rd_csv(file.path(rt, "star/edger/contrast_GBR_vs_YRI.usage.simes.csv"))
  dtx <- rd_tsv(file.path(rt, "salmon/dexseq_dtu/results/dexseq/DEXSeqResults.GBR-YRI.tsv")); stg <- rd_tsv(file.path(rt, "salmon/dexseq_dtu/results/stager/getAdjustedPValues.GBR-YRI.tsv"))
  sg <- function(p) !is.na(p) & p < fdr
  row <- function(tool, level, p, sig, eff = NULL) list(tool = tool, level = level, tested = sum(!is.na(p)), significant = sum(sig),
    up = if (is.null(eff)) NA_integer_ else sum(sig & eff > 0, na.rm = TRUE), down = if (is.null(eff)) NA_integer_ else sum(sig & eff < 0, na.rm = TRUE))
  e_dex <- -dex$log2fold_YRI_GBR; e_dtu <- -dtx$log2fold_YRI_GBR
  txs <- sg(stg$transcript); tx_eff <- e_dtu[match(stg$txID[txs], dtx$featureID)]
  g_tested <- unique(stg$geneID[!is.na(stg$gene)]); g_sig <- unique(stg$geneID[sg(stg$gene)])
  rows <- list(
    row("exon_dexseq", "gene", pg$padj, sg(pg$padj)), row("exon_dexseq", "exon_bin", dex$padj, sg(dex$padj), e_dex),
    row("exon_edger", "gene_simes", esm$FDR, sg(esm$FDR)), row("exon_edger", "gene_ftest", egn$FDR, sg(egn$FDR)), row("exon_edger", "exon", eex$FDR, sg(eex$FDR), eex$logFC),
    list(tool = "dtu", level = "gene_stager", tested = length(g_tested), significant = length(g_sig), up = NA_integer_, down = NA_integer_),
    list(tool = "dtu", level = "transcript_stager", tested = sum(!is.na(stg$transcript)), significant = sum(txs), up = sum(tx_eff > 0), down = sum(tx_eff < 0)),
    row("dtu", "transcript_dexseq", dtx$padj, sg(dtx$padj), e_dtu))
  list(counts = do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE)), pg_genes = pg$groupID[sg(pg$padj)], esm_genes = esm$Geneid[sg(esm$FDR)],
       dtu_genes = g_sig, dtu_tx = stg$txID[txs], dtu_tx_dir = stats::setNames(dir_txt(tx_eff), stg$txID[txs]), dex_bins = paste0(dex$groupID, ":", dex$featureID)[sg(dex$padj)],
       dex_bin_dir = stats::setNames(dir_txt(e_dex[sg(dex$padj)]), paste0(dex$groupID, ":", dex$featureID)[sg(dex$padj)]))
}
check_counts <- function(rd, fdr, tag) {
  ex <- raw_usage(fdr)$counts
  cn <- rd_csv(file.path(rd, "tables", paste0(C, ".usage.counts.csv")))
  k <- match(paste(ex$tool, ex$level), paste(cn$tool, cn$level))
  check(sprintf("the counts table has the %d rows (tool, level) of the recomputation at %s", nrow(ex), tag), !anyNA(k) && nrow(cn) == nrow(ex))
  if (anyNA(k)) return(invisible(NULL))
  cn <- cn[k, ]
  check(sprintf("tested counts equal the raw recomputation at %s (%s)", tag, paste(ex$tested, collapse = ",")), identical(as.integer(cn$tested), as.integer(ex$tested)))
  check(sprintf("significant counts equal the raw recomputation at %s (%s: %s)", tag, paste(ex$level, collapse = "/"), paste(ex$significant, collapse = ",")),
        identical(as.integer(cn$significant), as.integer(ex$significant)))
  d <- !is.na(ex$up)
  check(sprintf("direction counts (more used in treatment / in control) equal the raw recomputation, DEXSeq flipped, at %s", tag),
        identical(as.integer(cn[["more used in treatment"]][d]), as.integer(ex$up[d])) && identical(as.integer(cn[["more used in control"]][d]), as.integer(ex$down[d])) &&
        all(is.na(cn[["more used in treatment"]][!d])))
}
tab_ids <- function(rd, nm, col) { f <- file.path(rd, "tables", paste0(C, ".", nm, ".csv")); if (!file.exists(f)) return(NULL); utils::read.csv(f, stringsAsFactors = FALSE, colClasses = "character")[[col]] }

# 3a. Rmd 04 at the defaults (FDR 0.05)
rd1 <- new_dir("rd_default")
ok <- render(rd1, "01_import.Rmd")
check("01_import.Rmd renders for the five modules", ok && all(file.exists(file.path(rd1, "rds", paste0(C, c(".exon_dexseq", ".exon_edger", ".dtu"), ".rds")))))
if (!ok) finish()
ok4 <- render(rd1, "04_usage.Rmd")
check("04_usage.Rmd renders on the fixtures with the default thresholds", ok4 && file.exists(file.path(rd1, "04_usage.html")))
if (ok4) {
  h <- html_of(rd1, "04_usage.html"); ex <- raw_usage(0.05)
  check("fixture sanity: the real tables hold exon bins, genes and DTU rows below the FDR (so the counts below are not all zero)",
        ex$counts$significant[ex$counts$level == "exon_bin"] > 0 && ex$counts$significant[ex$counts$level == "gene_stager"] > 0 && ex$counts$significant[ex$counts$level == "transcript_stager"] > 0)
  check("04 prints the sign rule at the top: treatment minus control, more used in treatment and in control",
        has(h, "treatment minus control") && has(h, "more used in treatment") && has(h, "more used in control"))
  check_counts(rd1, 0.05, "FDR 0.05")
  check("tables/exon_dexseq.genes.csv holds exactly the genes of perGeneQValue below the FDR", setequal(tab_ids(rd1, "exon_dexseq.genes", "gene_id"), ex$pg_genes) && length(tab_ids(rd1, "exon_dexseq.genes", "gene_id")) == length(ex$pg_genes))
  check("tables/exon_edger.genes.csv exists and holds exactly the Simes genes below the FDR", !is.null(tab_ids(rd1, "exon_edger.genes", "gene_id")) && setequal(tab_ids(rd1, "exon_edger.genes", "gene_id"), ex$esm_genes) && length(tab_ids(rd1, "exon_edger.genes", "gene_id")) == length(ex$esm_genes))
  check("tables/dtu.genes.csv holds exactly the genes with stageR gene below the FDR", setequal(tab_ids(rd1, "dtu.genes", "gene_id"), ex$dtu_genes) && length(tab_ids(rd1, "dtu.genes", "gene_id")) == length(ex$dtu_genes))
  check("tables/dtu.transcripts.csv holds exactly the transcripts with stageR transcript below the FDR, with the direction of the flipped DEXSeq column",
        setequal(tab_ids(rd1, "dtu.transcripts", "transcript_id"), ex$dtu_tx) && identical(unname(stats::setNames(tab_ids(rd1, "dtu.transcripts", "direction"), tab_ids(rd1, "dtu.transcripts", "transcript_id"))[names(ex$dtu_tx_dir)]), unname(ex$dtu_tx_dir)))
  check("tables/exon_dexseq.bins.csv holds exactly the exon bins with padj below the FDR, with the direction of the flipped column",
        setequal(tab_ids(rd1, "exon_dexseq.bins", "feature_id"), ex$dex_bins) && identical(unname(stats::setNames(tab_ids(rd1, "exon_dexseq.bins", "direction"), tab_ids(rd1, "exon_dexseq.bins", "feature_id"))[names(ex$dex_bin_dir)]), unname(ex$dex_bin_dir)))
  check("04 prints 'none at FDR 0.05' for the empty edgeR results", length(gregexpr("none at FDR 0.05", h, fixed = TRUE)[[1]]) >= 2 && has(h, "three nearest"))
  check("04 prints the stageR alpha read from the stageRTx object (0.05) and names the gene and transcript columns",
        has(h, "overall alpha of 0.05") && has(h, "<code>gene</code>") && has(h, "<code>transcript</code>"))
  check("the volcano plots of the three tools exist and are not empty",
        all(sapply(c("exon_dexseq", "exon_edger", "dtu"), function(t) nonempty(file.path(rd1, "figures", sprintf("%s.usage.volcano.%s.png", C, t))))))
  check("the exon usage plot of the top DEXSeq genes and the isoform plot of the top DTU genes exist",
        nonempty(file.path(rd1, "figures", paste0(C, ".usage.exon_top.png"))) && nonempty(file.path(rd1, "figures", paste0(C, ".usage.dtu_top.png"))))
  check("04 states that the exon and transcript usage tools use the sample sheet design as it is", has(h, "design of the sample sheet"))
}
# 3b. a second threshold: a literal threshold in the Rmd would fail one of the two
rd2 <- new_dir("rd_fdr5")
render(rd2, "01_import.Rmd", min_reads = 0)
ok4b <- render(rd2, "04_usage.Rmd", fdr = 0.9)
check("04_usage.Rmd renders at FDR 0.9", ok4b)
if (ok4b) {
  check_counts(rd2, 0.9, "FDR 0.9")
  ex <- raw_usage(0.9)
  check("at FDR 0.9 the gene tables equal the recomputation (perGeneQValue, Simes, stageR)", setequal(tab_ids(rd2, "exon_dexseq.genes", "gene_id"), ex$pg_genes) &&
        setequal(tab_ids(rd2, "exon_edger.genes", "gene_id"), ex$esm_genes) && setequal(tab_ids(rd2, "dtu.genes", "gene_id"), ex$dtu_genes) && length(ex$esm_genes) > 0)
  check("at FDR 0.9 the page names the FDR 0.9", has(html_of(rd2, "04_usage.html"), "FDR 0.9"))
}
# 3c. events = FALSE and topgene = FALSE leave out tables, volcano plots and top-gene plots, the counts stay
rd3 <- new_dir("rd_off"); dir.create(file.path(rd3, "rds")); file.copy(list.files(file.path(rd1, "rds"), full.names = TRUE), file.path(rd3, "rds"))
ok4c <- render(rd3, "04_usage.Rmd", events = FALSE, topgene = FALSE)
check("04 renders with events and topgene switched off", ok4c)
if (ok4c) {
  check("with events and topgene off no volcano, no top-gene plot and no gene table is written, the counts table stays",
        length(list.files(file.path(rd3, "figures"), pattern = "volcano|_top")) == 0 && length(list.files(file.path(rd3, "tables"), pattern = "genes|bins|transcripts")) == 0 &&
        file.exists(file.path(rd3, "tables", paste0(C, ".usage.counts.csv"))))
  check("04 says the events section is switched off", has(html_of(rd3, "04_usage.html"), "switched off"))
}
# 3d. a missing RDS is reported, not an error; the same for Rmd 05
rdm <- new_dir("rd_missing"); dir.create(file.path(rdm, "rds"))
okm4 <- render(rdm, "04_usage.Rmd"); okm5 <- render(rdm, "05_cross_tool.Rmd")
check("04 and 05 render when no RDS of the contrast exists", okm4 && okm5)
if (okm4 && okm5) {
  h4 <- html_of(rdm, "04_usage.html"); h5 <- html_of(rdm, "05_cross_tool.html")
  check("04 says that the DEXSeq, edgeR and DTU RDS were not found", has(h4, paste0(C, ".exon_dexseq.rds")) && has(h4, paste0(C, ".exon_edger.rds")) && has(h4, paste0(C, ".dtu.rds")) && has(h4, "not found"))
  check("05 says that the RDS of every tool was not found and that there is no comparison", has(h5, paste0(C, ".rmats.rds")) && has(h5, paste0(C, ".dtu.rds")) && has(h5, "not found") && has(h5, "fewer than two"))
}
# 3e. only the gene-level raw files missing (outdir without them): the exon-level result is still shown
fb <- file.path(work, "fb_out"); dir.create(file.path(fb, "star/dexseq_exon/results"), recursive = TRUE)
rdf <- new_dir("rd_nopergene"); dir.create(file.path(rdf, "rds")); file.copy(list.files(file.path(rd1, "rds"), full.names = TRUE), file.path(rdf, "rds"))
okf <- render(rdf, "04_usage.Rmd", outdir = fb)
check("04 renders when the per-gene, Simes and stageR files of the run are absent, and names the missing files",
      okf && has(html_of(rdf, "04_usage.html"), "perGeneQValue.GBR_vs_YRI.csv") && has(html_of(rdf, "04_usage.html"), "was not found"))

# 3f. Rmd 05 on the fixtures: sets, intersections, directions
raw5 <- function(fdr, dpsi, min_reads) {
  pick <- function(df) { df <- df[order(df$padj, 1 - abs(df$eff), df$feat), ]; df <- df[!duplicated(df$gene), ]; stats::setNames(dir_txt(df$eff), df$gene) }
  # rMATS
  rm <- do.call(rbind, lapply(c("SE", "A5SS", "A3SS", "MXE", "RI"), function(ty) {
    r <- utils::read.delim(file.path(rt, sprintf("star/rmats/GBR_vs_YRI/%s.MATS.JC.txt", ty)), quote = "\"", stringsAsFactors = FALSE, colClasses = "character")
    sup <- function(a, b) sapply(seq_along(a), function(i) mean(as.numeric(strsplit(a[i], ",")[[1]]) + as.numeric(strsplit(b[i], ",")[[1]])))
    s <- pmin(sup(r$IJC_SAMPLE_1, r$SJC_SAMPLE_1), sup(r$IJC_SAMPLE_2, r$SJC_SAMPLE_2))
    data.frame(gene = r$GeneID, feat = paste0(ty, ":", r$ID), padj = num(r$FDR), eff = num(r$IncLevelDifference), low = s < min_reads, stringsAsFactors = FALSE) }))
  rm_sig <- !is.na(rm$padj) & !is.na(rm$eff) & rm$padj < fdr & abs(rm$eff) >= dpsi & !rm$low
  # SUPPA2: BH within each file over the tested rows; the effect is minus the dPSI of the file
  sp <- do.call(rbind, lapply(c("per_local_event/local_GBR-YRI.dpsi", "per_isoform/transcript_GBR-YRI.dpsi"), function(f) {
    t <- do.call(rbind, strsplit(readLines(file.path(rt, "salmon/suppa/diffsplice", f))[-1], "\t", fixed = TRUE))
    d <- num(t[, 2]); p <- num(t[, 3]); q <- rep(NA_real_, length(d)); ok <- !is.na(d) & !is.na(p); q[ok] <- stats::p.adjust(p[ok], method = "BH")
    data.frame(gene = sub(";.*$", "", t[, 1]), feat = sub("^[^;]*;", "", t[, 1]), padj = q, eff = -d, stringsAsFactors = FALSE) }))
  sp_sig <- !is.na(sp$padj) & !is.na(sp$eff) & sp$padj < fdr & abs(sp$eff) >= dpsi
  # DEXSeq exon: genes from perGeneQValue, bins from the exon table (merged groups A+B count for A and for B)
  ex_ <- function(df) { parts <- strsplit(df$gene, "+", fixed = TRUE); df[rep(seq_len(nrow(df)), lengths(parts)), ] -> o; o$gene <- unlist(parts); o }
  pg <- rd_csv(file.path(rt, "star/dexseq_exon/results/perGeneQValue.GBR_vs_YRI.csv")); dx <- rd_csv(file.path(rt, "star/dexseq_exon/results/DEXSeqResults.GBR_vs_YRI.csv"))
  dxf <- ex_(data.frame(gene = dx$groupID, feat = paste0(dx$groupID, ":", dx$featureID), padj = dx$padj, eff = -dx$log2fold_YRI_GBR, stringsAsFactors = FALSE))
  dxs <- !is.na(dxf$padj) & dxf$padj < fdr
  pgx <- ex_(data.frame(gene = pg$groupID, padj = pg$padj, stringsAsFactors = FALSE))
  # edgeR: Simes genes, exon bins by FDR
  es <- rd_csv(file.path(rt, "star/edger/contrast_GBR_vs_YRI.usage.simes.csv")); ee <- rd_csv(file.path(rt, "star/edger/contrast_GBR_vs_YRI.usage.exon.csv"))
  eef <- data.frame(gene = ee$Geneid, feat = paste0(ee$Geneid, ":", ee$Chr, ":", ee$Start, "-", ee$End), padj = ee$FDR, eff = ee$logFC, stringsAsFactors = FALSE)
  # DTU: genes and transcripts of stageR; the transcript with the smallest DEXSeq padj among the confirmed ones gives the direction
  dt <- rd_tsv(file.path(rt, "salmon/dexseq_dtu/results/dexseq/DEXSeqResults.GBR-YRI.tsv")); sg <- rd_tsv(file.path(rt, "salmon/dexseq_dtu/results/stager/getAdjustedPValues.GBR-YRI.tsv"))
  dtf <- data.frame(gene = dt$groupID, feat = dt$featureID, padj = dt$padj, eff = -dt$log2fold_YRI_GBR, stringsAsFactors = FALSE)
  conf <- paste(sg$geneID, sg$txID)[!is.na(sg$transcript) & sg$transcript < fdr]
  list(rmats = list(genes = unique(rm$gene[rm_sig]), dir = pick(rm[rm_sig, ])), suppa = list(genes = unique(sp$gene[sp_sig]), dir = pick(sp[sp_sig, ])),
       exon_dexseq = list(genes = unique(pgx$gene[!is.na(pgx$padj) & pgx$padj < fdr]), dir = pick(dxf[dxs, ])),
       exon_edger = list(genes = unique(es$Geneid[!is.na(es$FDR) & es$FDR < fdr]), dir = pick(eef[!is.na(eef$padj) & eef$padj < fdr, ])),
       dtu = list(genes = unique(sg$geneID[!is.na(sg$gene) & sg$gene < fdr]), dir = pick(dtf[paste(dtf$gene, dtf$feat) %in% conf, ])))
}
TOOLS <- c("rmats", "suppa", "exon_dexseq", "exon_edger", "dtu")
check_cross <- function(rd, fdr, dpsi, min_reads, tag) {
  ex <- raw5(fdr, dpsi, min_reads)
  mt <- rd_csv(file.path(rd, "tables", paste0(C, ".cross_tool.membership.csv")))
  check(sprintf("the membership table has one column per tool and the gene sets equal the recomputation at %s (%s genes)", tag, paste(sapply(ex, function(x) length(x$genes)), collapse = ",")),
        all(TOOLS %in% names(mt)) && all(sapply(TOOLS, function(t) setequal(mt$gene_id[mt[[t]] %in% c(TRUE, "TRUE")], ex[[t]]$genes))) && !anyDuplicated(mt$gene_id))
  univ <- unique(unlist(lapply(ex, function(x) x$genes)))
  check(sprintf("the union of the sets is the membership table (%d genes) at %s", length(univ), tag), setequal(mt$gene_id, univ))
  member <- sapply(TOOLS, function(t) univ %in% ex[[t]]$genes)
  combo <- apply(member, 1, function(r) paste(TOOLS[r], collapse = "&")); tb <- table(combo)
  ic <- rd_csv(file.path(rd, "tables", paste0(C, ".cross_tool.intersections.csv")))
  check(sprintf("the intersections (exact memberships) equal the recomputation at %s", tag), identical(stats::setNames(as.integer(ic$n_genes), ic$combination)[names(tb)], stats::setNames(as.integer(tb), names(tb))) && sum(ic$n_genes) == length(univ))
  xt <- utils::read.csv(file.path(rd, "tables", paste0(C, ".cross_tool.genes.csv")), stringsAsFactors = FALSE, colClasses = "character")
  two <- univ[rowSums(member) >= 2]
  check(sprintf("the cross-tool table holds exactly the %d genes called by two or more tools at %s", length(two), tag), setequal(xt$gene_id, two) && length(xt$gene_id) == length(two))
  if (length(two) > 0) {
    dirs_ok <- sapply(TOOLS, function(t) { col <- xt[[paste0("direction_", t)]]; want <- unname(ex[[t]]$dir[xt$gene_id]); want[!xt$gene_id %in% ex[[t]]$genes] <- NA
      identical(is.na(col), is.na(want)) && all(col[!is.na(col)] == want[!is.na(want)]) })
    check(sprintf("each tool's direction in the cross-tool table equals the recomputation (sign of the smallest-padj significant feature; DEXSeq and SUPPA2 flipped) at %s", tag), all(dirs_ok))
    cc <- rd_csv(file.path(rd, "tables", paste0(C, ".cross_tool.concordance.csv")))
    pr <- utils::combn(TOOLS, 2)
    expd <- do.call(rbind, lapply(seq_len(ncol(pr)), function(j) { a <- pr[1, j]; b <- pr[2, j]; both <- intersect(ex[[a]]$genes, ex[[b]]$genes)
      da <- unname(ex[[a]]$dir[both]); db <- unname(ex[[b]]$dir[both]); k <- !is.na(da) & !is.na(db)
      data.frame(a = a, b = b, both = length(both), same = sum(k & da == db), opp = sum(k & da != db)) }))
    kk <- match(paste(expd$a, expd$b), paste(cc$tool_a, cc$tool_b))
    check(sprintf("the concordance counts (same / opposite direction) equal the recomputation for every pair of tools at %s", tag),
          !anyNA(kk) && identical(as.integer(cc$genes_in_both[kk]), as.integer(expd$both)) && identical(as.integer(cc$same_direction[kk]), as.integer(expd$same)) && identical(as.integer(cc$opposite_direction[kk]), as.integer(expd$opp)))
  }
  ex
}
rd5a <- new_dir("rd5_default"); ok <- render(rd5a, "01_import.Rmd"); ok5 <- ok && render(rd5a, "05_cross_tool.Rmd")
check("05_cross_tool.Rmd renders on the fixtures with the default thresholds", ok5 && file.exists(file.path(rd5a, "05_cross_tool.html")))
if (ok5) {
  h <- html_of(rd5a, "05_cross_tool.html")
  check("05 prints the sign rule at the top, the direction rule and the statement that only directions are compared", has(h, "treatment minus control") && has(h, "directions compared only, never effect sizes across scales") && has(h, "Background:"))
  ex <- check_cross(rd5a, 0.05, 0.1, 10, "FDR 0.05, |dPSI| 0.1, min_reads 10")
  check("at the defaults rMATS has no significant gene (all low support) and DTU and SUPPA2 have some (the sets are not empty by accident)", length(ex$rmats$genes) == 0 && length(ex$dtu$genes) > 0)
  check("the UpSet figure exists when two or more tools have a gene", nonempty(file.path(rd5a, "figures", paste0(C, ".upset.png"))))
  check("05 reports which UpSet engine was used", has(h, "ComplexUpset") || has(h, "UpSetR"))
}
rd5b <- new_dir("rd5_perm"); ok <- render(rd5b, "01_import.Rmd", min_reads = 0); ok5b <- ok && render(rd5b, "05_cross_tool.Rmd", fdr = 0.9, dpsi = 0.05, min_reads = 0)
check("05_cross_tool.Rmd renders at FDR 0.9, |dPSI| 0.05, min_reads 0", ok5b)
if (ok5b) {
  ex <- check_cross(rd5b, 0.9, 0.05, 0, "FDR 0.9, |dPSI| 0.05, min_reads 0")
  two <- sum(table(unlist(lapply(ex, function(x) x$genes))) >= 2)
  check(sprintf("at the permissive thresholds some genes are called by two or more tools (%d), so directions and concordance are exercised", two), two > 0)
  check("the UpSet figure exists at the permissive thresholds", nonempty(file.path(rd5b, "figures", paste0(C, ".upset.png"))))
}
# a single tool in the modules: no comparison, the page says so
rd5c <- new_dir("rd5_one"); ok <- render(rd5c, "01_import.Rmd", modules = "dtu"); ok5c <- ok && render(rd5c, "05_cross_tool.Rmd", modules = "dtu")
check("05 renders with one tool and says that fewer than two tools are available", ok5c && has(html_of(rd5c, "05_cross_tool.html"), "fewer than two") && !file.exists(file.path(rd5c, "figures", paste0(C, ".upset.png"))))
finish()
