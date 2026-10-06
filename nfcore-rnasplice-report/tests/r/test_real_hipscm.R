# Usage (inside the R image, on a compute node only; tests/r/run_r_tests.sh does this through sbatch):
#   Rscript test_real_hipscm.R <skill.md> <tests/fixtures dir>
# Tests on the trimmed real tables of the human hiPSC-CM run (tests/fixtures/real_hipscm/, contrast KO_vs_WT, treatment KO, control WT):
#   1. the code that 04, 05 and 06 share is the same function in each (direction of a gene, calls, the DTU flag);
#   2. SYNTHETIC ties of the gene direction: two features of opposite sign with the same smallest padj (the larger absolute effect wins,
#      in every row order), a three-way tie, a tie with equal absolute effects (the first feature identifier wins), a missing effect;
#   3. the real rows of sign_rows.tsv (fixture real_hipscm): the raw value is the value of the fixture table, and harmonise_effect gives
#      the expected sign;
#   4. the whole chain 01 to 07 rendered on the fixtures at the defaults (FDR 0.05, |dPSI| 0.1, min_reads 10; 06 on the DEXSeq exon
#      usage genes only, see MODS_06 below), then: the RDS effect of
#      every real sign row, the counts per type and level against a base-R recomputation from the trimmed tables, TEAD3 (rMATS SE 199981
#      and MXE 42852 both FDR 0: the larger change wins), RBPMS2 (control_up in rMATS, SUPPA2, DEXSeq and edgeR), the DTU flag
#      `effect_flag` (absolute log2 fold change above 10), the order of the DEXSeq bins table, and the two gene counts of 06.
# Prints one "ok" or "FAIL:" line per check, then "R TESTS PASS (<n>)" or "R TESTS FAIL (<failed>/<n>)" (exit 1).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: test_real_hipscm.R <skill.md> <fixtures dir>")
skill <- args[1]; fx <- args[2]; rt <- file.path(fx, "real_hipscm")
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
work <- tempfile("rsr_real_"); dir.create(work)
txt <- readLines(skill, warn = FALSE)
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
NAMES <- c("01_import.Rmd", "02_rmats.Rmd", "03_suppa.Rmd", "04_usage.Rmd", "05_cross_tool.Rmd", "06_enrichment.Rmd", "07_summary.Rmd")
RMD <- stats::setNames(lapply(NAMES, get_rmd), NAMES)
E1 <- functions_env(RMD[["01_import.Rmd"]], "01"); E4 <- functions_env(RMD[["04_usage.Rmd"]], "04")
E5 <- functions_env(RMD[["05_cross_tool.Rmd"]], "05"); E6 <- functions_env(RMD[["06_enrichment.Rmd"]], "06")
num <- function(x) suppressWarnings(as.numeric(ifelse(x %in% c("nan", "NA", ""), NA, x)))
rd_csv <- function(f) utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
rd_chr <- function(f) utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character", na.strings = character(0))
rd_tsv <- function(f) utils::read.delim(f, stringsAsFactors = FALSE, check.names = FALSE, quote = "")

# --- 1. the shared code is one function in every report that uses it
same <- function(nm, envs) { d <- lapply(envs, function(e) if (exists(nm, envir = e, inherits = FALSE)) deparse(get(nm, envir = e)) else NULL)
  !any(vapply(d, is.null, logical(1))) && all(vapply(d[-1], identical, logical(1), d[[1]])) }
check("top_by_gene (the direction of a gene) is the same function in 04, 05 and 06", same("top_by_gene", list(E4, E5, E6)))
for (nm in c("strip_version", "strip_tx_version", "expand_genes", "tidy_ids", "sig_usage", "sig_events", "make_calls", "calls_from_events",
             "calls_from_gene_table", "calls_from_stager"))
  check(sprintf("%s is the same function in 05 and 06", nm), same(nm, list(E5, E6)))
check("effect_flag and ARTEFACT_LOG2FC (10) are the same in 04 and 05", same("effect_flag", list(E4, E5)) && identical(E4$ARTEFACT_LOG2FC, 10) && identical(E5$ARTEFACT_LOG2FC, 10))
check("effect_flag marks an absolute effect above 10 only (10 itself, NA and 9.9 are not marked)",
      identical(E5$effect_flag(c(-21.5, 22.6, 10, 9.9, NA, -10.01)), c("zero-count artefact", "zero-count artefact", "", "", "", "zero-count artefact")))

# --- 2. SYNTHETIC ties of the gene direction, in every Rmd that has top_by_gene, in every row order
SYNTHETIC_tie <- data.frame(
  gene_id    = c("T1", "T1", "T1", "T2", "T2", "T2", "T3", "T3", "T4", "T4"),
  feature_id = c("T1:a", "T1:b", "T1:c", "T2:a", "T2:b", "T2:c", "T3:b", "T3:a", "T4:a", "T4:b"),
  padj       = c(0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0),
  effect     = c(-0.3, 0.5, -0.9, 0.2, -0.7, 0.6, 0.4, -0.4, NA, 0.2), stringsAsFactors = FALSE)
SYNTHETIC_tie$direction <- ifelse(is.na(SYNTHETIC_tie$effect), "none", ifelse(SYNTHETIC_tie$effect > 0, "treatment_up", "control_up"))
perms <- c(list(seq_len(nrow(SYNTHETIC_tie)), rev(seq_len(nrow(SYNTHETIC_tie)))), local({ set.seed(42); lapply(1:20, function(i) sample(nrow(SYNTHETIC_tie))) }))
tops <- function(E) lapply(perms, function(p) { d <- SYNTHETIC_tie[p, ]; tg <- E$top_by_gene(d, rep(TRUE, nrow(d))); tg[order(tg$gene_id), c("gene_id", "top_feature", "direction")] })
for (nm in c("04", "05", "06")) {
  E <- list("04" = E4, "05" = E5, "06" = E6)[[nm]]; tt <- tops(E)
  pick <- function(g) unique(vapply(tt, function(x) paste(x$top_feature[x$gene_id == g], x$direction[x$gene_id == g]), character(1)))
  check(sprintf("SYNTHETIC tie (%s): two features with padj 0 and opposite signs: the larger absolute effect wins (T1:b treatment_up, not T1:a, not T1:c with a larger padj) in all %d row orders", nm, length(perms)),
        identical(pick("T1"), "T1:b treatment_up"))
  check(sprintf("SYNTHETIC tie (%s): a three-way tie at padj 0: the largest absolute effect wins (T2:b control_up) in all row orders", nm), identical(pick("T2"), "T2:b control_up"))
  check(sprintf("SYNTHETIC tie (%s): equal padj and equal absolute effects: the first feature identifier wins (T3:a control_up) in all row orders", nm), identical(pick("T3"), "T3:a control_up"))
  check(sprintf("SYNTHETIC tie (%s): a feature without an effect loses a tie (T4:b treatment_up)", nm), identical(pick("T4"), "T4:b treatment_up"))
}

# --- 3. the real rows of sign_rows.tsv: raw value from the fixture table, harmonised sign
sr <- read.delim(file.path(fx, "sign_rows.tsv"), comment.char = "#", quote = "", stringsAsFactors = FALSE, colClasses = "character")
sr <- sr[sr$fixture == "real_hipscm", , drop = FALSE]
tool_of <- c(rMATS = "rmats", SUPPA2 = "suppa", DEXSeq_exon = "exon_dexseq", DEXSeq_DTU = "dtu", edgeR_exon = "exon_edger")
want <- c("+" = 1, "-" = -1)
check(sprintf("sign_rows.tsv has at least 6 real_hipscm rows for each of the five tools (%s)", paste(table(sr$tool), collapse = ",")),
      all(names(tool_of) %in% sr$tool) && all(table(sr$tool)[names(tool_of)] >= 6) && all(sr$treatment == "KO") && all(sr$control == "WT") && all(nzchar(sr$verification)))
P <- function(...) file.path(rt, ...)
rm_tab <- lapply(stats::setNames(nm = c("SE", "A5SS", "A3SS", "MXE", "RI")), function(ty)
  utils::read.delim(P("star/rmats/KO_vs_WT", sprintf("%s.MATS.JC.txt", ty)), quote = "\"", stringsAsFactors = FALSE, colClasses = "character", check.names = FALSE))
sup_lines <- function(f) do.call(rbind, strsplit(readLines(f)[-1], "\t", fixed = TRUE))
su_loc <- sup_lines(P("salmon/suppa/diffsplice/per_local_event/local_KO-WT.dpsi")); su_iso <- sup_lines(P("salmon/suppa/diffsplice/per_isoform/transcript_KO-WT.dpsi"))
dex <- rd_csv(P("star/dexseq_exon/results/DEXSeqResults.KO_vs_WT.csv")); dtu <- rd_tsv(P("salmon/dexseq_dtu/results/dexseq/DEXSeqResults.KO-WT.tsv"))
edx <- rd_csv(P("star/edger/contrast_KO_vs_WT.usage.exon.csv"))
raw_of <- function(tool, feature) {
  switch(tool,
    rMATS = { ty <- sub("\\..*$", "", feature); id <- sub("^[^=]*=([0-9]+):.*$", "\\1", feature); t <- rm_tab[[ty]]; t$IncLevelDifference[t[[1]] == id] },
    SUPPA2 = su_loc[su_loc[, 1] == feature, 2],
    DEXSeq_exon = dex$log2fold_WT_KO[paste0(dex$groupID, ":", dex$featureID) == feature],
    DEXSeq_DTU = dtu$log2fold_WT_KO[paste0(dtu$groupID, ":", dtu$featureID) == feature],
    edgeR_exon = edx$logFC[paste0(edx$Geneid, ":", edx$Chr, ":", edx$Start, "-", edx$End) == feature])
}
for (tl in names(tool_of)) {
  s <- sr[sr$tool == tl, , drop = FALSE]
  found <- vapply(seq_len(nrow(s)), function(i) { v <- raw_of(tl, s$feature[i]); length(v) == 1 && isTRUE(all.equal(as.numeric(v), as.numeric(s$raw_value[i]))) }, logical(1))
  check(sprintf("real hipscm: every %s sign row is one row of the fixture table with the raw value of the row (%d rows)", tl, nrow(s)), all(found))
  check(sprintf("real hipscm: harmonise_effect gives the expected sign for every %s sign row", tl),
        all(sign(E1$harmonise_effect(tool_of[[tl]], as.numeric(s$raw_value))) == unname(want[s$expected_sign])))
}

# --- 4. the chain 01 to 07 on the fixtures, at the defaults
rd <- file.path(work, "report"); dir.create(rd)
for (nm in NAMES) writeLines(RMD[[nm]], file.path(rd, nm))
pr <- list(outdir = rt, report_dir = rd, contrasts = P("261006_project_contrasts.csv"), samplesheet = P("261006_project_samplesheet.csv"),
           modules = "rmats suppa dexseq_exon edger_exon dtu", rmats_type = "JC", fdr = 0.05, dpsi = 0.1, min_reads = 10, events = TRUE, cross = TRUE,
           topgene = TRUE, enrich = TRUE, organism = "human", gtf = "", enrich_dbs = "GO", paired_rmats = FALSE, paired_suppa = FALSE)
# 06 runs on the DEXSeq exon usage genes only: every GO enrichment takes 6 to 25 s in the image, and 06 with all five tools (45
# enrichments) takes about 9 minutes on these tables, which every R mutation row would pay again. The five-tool run of 06 on the
# whole real output was part of the real-data run; the gene sets of all five tools are tested here through 05 (same functions).
MODS_06 <- "dexseq_exon"
ok_all <- TRUE
for (nm in NAMES) {
  of <- if (nm == "07_summary.Rmd") "report_summary.html" else NULL
  p_nm <- if (nm == "06_enrichment.Rmd") utils::modifyList(pr, list(modules = MODS_06)) else pr
  t0 <- Sys.time()
  ok <- tryCatch({ rmarkdown::render(file.path(rd, nm), output_file = of, output_dir = rd, knit_root_dir = rd, params = p_nm, quiet = TRUE, envir = new.env()); TRUE },
                 error = function(e) { cat("render error in", nm, ":", conditionMessage(e), "\n"); FALSE })
  cat(sprintf("     %s rendered in %.0f s\n", nm, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  check(sprintf("real hipscm: %s renders on the trimmed real tables", nm), ok)
  ok_all <- ok_all && ok
}
if (!ok_all) finish()
C <- "KO_vs_WT"
TB <- function(nm) file.path(rd, "tables", paste0(C, ".", nm, ".csv"))
html_of <- function(nm) gsub("[[:space:]]+", " ", paste(readLines(file.path(rd, nm), warn = FALSE), collapse = " "))
has <- function(h, s) grepl(s, h, fixed = TRUE)
rds <- function(tool) readRDS(file.path(rd, "rds", sprintf("%s.%s.rds", C, tool)))

# 4a. the RDS effect of every real sign row: the expected sign and the size of the raw value
for (tl in names(tool_of)) {
  s <- sr[sr$tool == tl, , drop = FALSE]; g <- rds(tool_of[[tl]])
  fid <- switch(tl, rMATS = paste0(sub("\\..*$", "", s$feature), ":", sub("^[^=]*=([0-9]+):.*$", "\\1", s$feature)), SUPPA2 = sub("^[^;]*;", "", s$feature),
                DEXSeq_DTU = sub("^[^:]*:", "", s$feature), s$feature)
  gid <- if (tl == "SUPPA2") sub(";.*$", "", s$feature) else NULL
  k <- if (is.null(gid)) match(fid, g$feature_id) else match(paste(gid, fid), paste(g$gene_id, g$feature_id))
  check(sprintf("real hipscm: the RDS effect of every %s sign row has the expected sign and the size of the raw value", tl),
        !anyNA(k) && all(sign(g$effect[k]) == unname(want[s$expected_sign])) && isTRUE(all.equal(abs(g$effect[k]), abs(as.numeric(s$raw_value)))) &&
        all(g$direction[k] == ifelse(s$expected_sign == "+", "treatment_up", "control_up")))
}

# 4b. counts per type and level against a base-R recomputation from the trimmed tables
mean_list <- function(x) vapply(strsplit(x, ",", fixed = TRUE), function(v) mean(as.numeric(v)), numeric(1))
rm_all <- do.call(rbind, lapply(names(rm_tab), function(ty) { t <- rm_tab[[ty]]
  s <- pmin(mean_list(t$IJC_SAMPLE_1) + mean_list(t$SJC_SAMPLE_1), mean_list(t$IJC_SAMPLE_2) + mean_list(t$SJC_SAMPLE_2))
  data.frame(type = ty, gene = t$GeneID, id = paste0(ty, ":", t[[1]]), fdr = num(t$FDR), d = num(t$IncLevelDifference), low = s < 10, stringsAsFactors = FALSE) }))
rm_all$sig <- !is.na(rm_all$fdr) & rm_all$fdr < 0.05 & !is.na(rm_all$d) & abs(rm_all$d) >= 0.1 & !rm_all$low
exp_rm <- do.call(rbind, lapply(names(rm_tab), function(ty) { i <- rm_all$type == ty
  c(events = sum(i), tested = sum(i & !is.na(rm_all$fdr)), low_support = sum(i & rm_all$low), significant = sum(i & rm_all$sig), up = sum(i & rm_all$sig & rm_all$d > 0), down = sum(i & rm_all$sig & rm_all$d < 0)) }))
cn <- rd_csv(TB("rmats.counts"))
check(sprintf("real hipscm: rMATS counts per type (events, tested, low_support, significant, by direction) equal the recomputation (significant %s)", paste(exp_rm[, "significant"], collapse = "/")),
      identical(cn$type, names(rm_tab)) && identical(dim(as.matrix(cn[, 2:7])), dim(exp_rm)) && all(as.matrix(cn[, 2:7]) == exp_rm) &&
      sum(exp_rm[, "significant"]) > 0 && sum(exp_rm[, "low_support"]) > 0)
bh <- function(d, p) { ok <- !is.na(d) & !is.na(p); q <- rep(NA_real_, length(p)); q[ok] <- stats::p.adjust(p[ok], method = "BH"); q }
su_frame <- function(m, kind) { d <- num(m[, 2]); p <- num(m[, 3]); ev <- sub("^[^;]*;", "", m[, 1])
  data.frame(gene = sub(";.*$", "", m[, 1]), id = ev, group = if (kind == "local") sub(":.*$", "", ev) else "isoform", eff = -d, padj = bh(d, p), stringsAsFactors = FALSE) }
su_all <- rbind(su_frame(su_loc, "local"), su_frame(su_iso, "isoform"))
su_all$sig <- !is.na(su_all$padj) & su_all$padj < 0.05 & !is.na(su_all$eff) & abs(su_all$eff) >= 0.1
cs <- rd_csv(TB("suppa.counts"))
exp_su <- t(vapply(cs$type, function(ty) { i <- su_all$group == ty
  c(sum(i), sum(i & is.na(su_all$eff)), sum(i & !is.na(su_all$padj)), sum(i & su_all$sig), sum(i & su_all$sig & su_all$eff > 0), sum(i & su_all$sig & su_all$eff < 0)) }, numeric(6)))
check(sprintf("real hipscm: SUPPA2 counts per type (events, untested, tested, significant, by direction) equal the recomputation, BH per file (types %s)", paste(cs$type, collapse = ",")),
      all(unique(su_all$group) %in% cs$type) && identical(dim(as.matrix(cs[, 2:7])), dim(exp_su)) && all(as.matrix(cs[, 2:7]) == exp_su) && sum(exp_su[, 4]) > 0)
pg <- rd_csv(P("star/dexseq_exon/results/perGeneQValue.KO_vs_WT.csv")); esm <- rd_csv(P("star/edger/contrast_KO_vs_WT.usage.simes.csv"))
egn <- rd_csv(P("star/edger/contrast_KO_vs_WT.usage.gene.csv")); stg <- rd_tsv(P("salmon/dexseq_dtu/results/stager/getAdjustedPValues.KO-WT.tsv"))
sg <- function(p) !is.na(p) & p < 0.05
txs <- sg(stg$transcript); e_dtu <- -dtu$log2fold_WT_KO; tx_eff <- e_dtu[match(paste(stg$geneID, stg$txID)[txs], paste(dtu$groupID, dtu$featureID))]
e_dex <- -dex$log2fold_WT_KO
exp_u <- rbind(
  c("exon_dexseq", "exon_bin", sum(!is.na(dex$padj)), sum(sg(dex$padj)), sum(sg(dex$padj) & e_dex > 0, na.rm = TRUE), sum(sg(dex$padj) & e_dex < 0, na.rm = TRUE)),
  c("exon_dexseq", "gene", sum(!is.na(pg$padj)), sum(sg(pg$padj)), NA, NA),
  c("exon_edger", "exon", sum(!is.na(edx$FDR)), sum(sg(edx$FDR)), sum(sg(edx$FDR) & edx$logFC > 0), sum(sg(edx$FDR) & edx$logFC < 0)),
  c("exon_edger", "gene_ftest", sum(!is.na(egn$FDR)), sum(sg(egn$FDR)), NA, NA), c("exon_edger", "gene_simes", sum(!is.na(esm$FDR)), sum(sg(esm$FDR)), NA, NA),
  c("dtu", "gene_stager", length(unique(stg$geneID[!is.na(stg$gene)])), length(unique(stg$geneID[sg(stg$gene)])), NA, NA),
  c("dtu", "transcript_stager", sum(!is.na(stg$transcript)), sum(txs), sum(tx_eff > 0), sum(tx_eff < 0)),
  c("dtu", "transcript_dexseq", sum(!is.na(dtu$padj)), sum(sg(dtu$padj)), sum(sg(dtu$padj) & e_dtu > 0, na.rm = TRUE), sum(sg(dtu$padj) & e_dtu < 0, na.rm = TRUE)))
cu <- rd_chr(TB("usage.counts")); k <- match(paste(exp_u[, 1], exp_u[, 2]), paste(cu$tool, cu$level))
got_u <- if (anyNA(k)) NULL else cbind(cu$tool[k], cu$level[k], cu$tested[k], cu$significant[k], cu[["more used in treatment"]][k], cu[["more used in control"]][k])
if (!is.null(got_u)) got_u[got_u %in% c("", "NA")] <- NA
check(sprintf("real hipscm: usage counts per tool and level (tested, significant, by direction) equal the recomputation (significant %s)", paste(exp_u[, 4], collapse = "/")),
      !is.null(got_u) && nrow(cu) == nrow(exp_u) && identical(unname(got_u), unname(exp_u)))
h4 <- html_of("04_usage.html")
check("real hipscm: stageR's alpha of the run (0.05) is read from the stageRTx object and the run's table is used as it is", has(h4, "overall alpha of 0.05") && has(h4, "equal to the FDR of this report"))

# 4c. TEAD3 (ENSG00000007866): SE 199981 (-0.467) and MXE 42852 (+0.4) both have FDR 0; the larger change gives the direction
xt <- rd_chr(TB("cross_tool.genes")); xr <- function(g) xt[xt$gene_id == g, , drop = FALSE]
te <- rm_all[gsub("\"", "", rm_all$gene) == "ENSG00000007866" & rm_all$sig, , drop = FALSE]
tied <- te[te$fdr == min(te$fdr), , drop = FALSE]
check("real hipscm fixture: TEAD3 has significant rMATS events tied at FDR 0 with opposite signs, among them SE:199981 (-0.467) and MXE:42852 (+0.4)",
      min(te$fdr) == 0 && all(c("SE:199981", "MXE:42852") %in% tied$id) && any(tied$d > 0) && any(tied$d < 0))
best <- tied[order(-abs(tied$d), tied$id), ][1, ]
check(sprintf("real hipscm: TEAD3 rMATS direction is that of the tied event with the largest change (%s, %s): control_up, effect -0.467", best$id, best$d),
      nrow(xr("ENSG00000007866")) == 1 && identical(xr("ENSG00000007866")$direction_rmats, "control_up") && isTRUE(all.equal(as.numeric(xr("ENSG00000007866")$effect_rmats), best$d)) &&
      isTRUE(all.equal(best$d, -0.467)))
g_rm <- rds("rmats"); set.seed(7)
dirs <- vapply(1:5, function(i) { cl <- E5$calls_from_events(g_rm[sample(nrow(g_rm)), ], 0.05, 0.1); paste(cl$dir[cl$dir$gene_id == "ENSG00000007866", c("top_feature", "direction")], collapse = " ") }, character(1))
check("real hipscm: the TEAD3 direction does not change when the rows of the rMATS RDS are shuffled (5 orders)", length(unique(dirs)) == 1 && identical(dirs[1], "SE:199981 control_up"))
# ENSG00000164754: two stageR-confirmed transcripts with DEXSeq padj 0 (DTU); the larger absolute change gives the direction
dt <- data.frame(gene = dtu$groupID, tx = dtu$featureID, padj = dtu$padj, eff = e_dtu, stringsAsFactors = FALSE)
conf <- paste(stg$geneID, stg$txID)[txs]
d64 <- dt[dt$gene == "ENSG00000164754" & paste(dt$gene, dt$tx) %in% conf, , drop = FALSE]; d64 <- d64[d64$padj == min(d64$padj), , drop = FALSE]
b64 <- d64[order(-abs(d64$eff), d64$tx), ][1, ]
check(sprintf("real hipscm: ENSG00000164754 (DTU, %d confirmed transcripts tied at padj 0) takes the larger change (%s, %.2f)", nrow(d64), b64$tx, b64$eff),
      nrow(d64) >= 2 && min(d64$padj) == 0 && nrow(xr("ENSG00000164754")) == 1 && isTRUE(all.equal(as.numeric(xr("ENSG00000164754")$effect_dtu), b64$eff)))

# 4d. RBPMS2 (ENSG00000166831): the exon skipped in KO; every exon-level tool says control_up
rb <- xr("ENSG00000166831")
check("real hipscm: RBPMS2 is control_up in rMATS (-0.89), SUPPA2, DEXSeq exon and edgeR exon (-1.588), as on the whole run",
      nrow(rb) == 1 && identical(rb$direction_rmats, "control_up") && isTRUE(all.equal(as.numeric(rb$effect_rmats), -0.89)) && identical(rb$direction_suppa, "control_up") &&
      identical(rb$direction_exon_dexseq, "control_up") && identical(rb$direction_exon_edger, "control_up") && abs(as.numeric(rb$effect_exon_edger) + 1.5879454015777) < 1e-9)

# 4e. the DTU flag: effect_flag is 'zero-count artefact' exactly when the absolute log2 fold change is above 10
tx <- rd_chr(TB("dtu.transcripts")); gt <- rd_chr(TB("dtu.genes"))
flag_ok <- function(eff, fl) identical(fl, ifelse(!is.na(eff) & abs(eff) > 10, "zero-count artefact", ""))
check("real hipscm: dtu.transcripts.csv has effect_flag, 'zero-count artefact' exactly where abs(effect) > 10 (ENST00000958519 -21.54 flagged, ENST00000883954 5.68 not)",
      "effect_flag" %in% names(tx) && flag_ok(num(tx$effect), tx$effect_flag) && identical(tx$effect_flag[tx$transcript_id == "ENST00000958519"], "zero-count artefact") &&
      identical(tx$effect_flag[tx$transcript_id == "ENST00000883954"], "") && isTRUE(all.equal(num(tx$effect[tx$transcript_id == "ENST00000958519"]), 21.5415458514456)))
check("real hipscm: dtu.genes.csv has effect_flag on the effect of the top transcript", "effect_flag" %in% names(gt) && flag_ok(num(gt$top_effect), gt$effect_flag) && any(gt$effect_flag != ""))
check("real hipscm: cross_tool.genes.csv has effect_flag for effect_dtu (RBPMS2 DTU +22.63 flagged), and effect_dtu is not capped",
      "effect_flag" %in% names(xt) && flag_ok(num(xt$effect_dtu), xt$effect_flag) && identical(rb$effect_flag, "zero-count artefact") && max(abs(num(xt$effect_dtu)), na.rm = TRUE) > 20)
check("real hipscm: 04 says that the DTU threshold of 10 is a heuristic", has(h4, "zero-count artefact") && has(h4, "heuristic"))

# 4f. the DEXSeq bins table: bins with an effect first, then padj, then the larger absolute effect
bn <- rd_csv(TB("exon_dexseq.bins")); hasE <- !is.na(bn$effect)
nae <- dex[sg(dex$padj) & is.na(dex$log2fold_WT_KO), , drop = FALSE]
check("real hipscm fixture: some significant bins of a group of overlapping genes have no effect and a smaller padj than significant bins with an effect (ordered by padj alone they would come first)",
      nrow(nae) > 0 && min(nae$padj) < max(dex$padj[sg(dex$padj) & !is.na(dex$log2fold_WT_KO)]))
e1 <- bn[hasE, , drop = FALSE]
check("real hipscm: the DEXSeq bins table lists every bin with an effect before the bins without one, by padj and then the larger absolute effect",
      sum(hasE) > 0 && all(hasE[seq_len(sum(hasE))]) && !is.unsorted(e1$padj) && all(vapply(split(abs(e1$effect), e1$padj), function(v) !is.unsorted(-v), logical(1))) && any(!hasE))

# 4g. 06 prints the pipeline's own count of DEXSeq groups and the genes after the split, and says what genes means
grp <- unique(pg$groupID[sg(pg$padj)]); split_genes <- unique(unlist(strsplit(grp, "+", fixed = TRUE)))
h6 <- html_of("06_enrichment.html"); sets <- rd_csv(TB("cross_tool.sets"))
check(sprintf("real hipscm fixture: DEXSeq calls groups of overlapping genes (%d groups, %d genes after the split)", length(grp), length(split_genes)), length(split_genes) > length(grp))
check(sprintf("real hipscm: 06 prints both DEXSeq numbers (%d identifiers as the tool writes them, %d genes after the split)", length(grp), length(split_genes)),
      has(h6, sprintf("Called by DEXSeq exon: %d identifiers as the tool writes them (the pipeline", length(grp))) && has(h6, "own count of <code>perGeneQValue</code> groups below the FDR;") &&
      has(h6, sprintf("and %d genes after splitting the groups", length(split_genes))) && has(h6, "genes after splitting the groups of overlapping genes"))
check("real hipscm: the sets table has the identifiers before the split next to the genes (DEXSeq)",
      sets$significant_ids_before_split[sets$tool == "exon_dexseq"] == length(grp) && sets$significant_genes[sets$tool == "exon_dexseq"] == length(split_genes))
sm <- rd_chr(TB("enrich.summary"))
check("real hipscm: the enrichment summary counts genes after the split (DEXSeq group all), and the GO enrichment of that group ran",
      as.integer(sm$genes[sm$tool == "exon_dexseq" & sm$group == "all"]) == length(split_genes) && identical(sm$status[sm$tool == "exon_dexseq" & sm$group == "all"], "ran"))
check("real hipscm: every report and the summary page exist", all(file.exists(file.path(rd, c(sub("Rmd$", "html", NAMES[1:6]), "report_summary.html")))))
finish()
