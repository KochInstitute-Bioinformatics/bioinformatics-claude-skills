# Usage (inside the R image, on a compute node only; tests/r/run_r_tests.sh does this through sbatch):
#   Rscript test_harmonise.R <skill.md> <tests/fixtures dir>
# Extracts the Rmd block under the heading **File `01_import.Rmd`.** of the skill, evaluates its `functions` chunk and checks
# the sign rules, the readers and the QC counts on the real fixtures; then renders the whole Rmd on the fixture tree and
# checks the RDS files against an independent recomputation from the raw tables.
# Prints one "ok" or "FAIL:" line per check, then "R TESTS PASS (<n>)" or "R TESTS FAIL (<failed>/<n>)" (exit 1).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: test_harmonise.R <skill.md> <fixtures dir>")
skill <- args[1]; fx <- args[2]; rt <- file.path(fx, "real_testdata")
n_ok <- 0; n_fail <- 0
# check(name, cond): stopifnot with a named message; a failure is printed and counted, the later checks still run
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
err_of <- function(expr) tryCatch({ force(expr); "" }, error = function(e) conditionMessage(e))
COLS <- c("contrast", "tool", "gene_id", "gene_symbol", "feature_type", "feature_id", "effect", "effect_scale", "padj",
          "direction", "n_support")
work <- tempfile("rsr_test_"); dir.create(work)

# --- extract the Rmd block and its functions chunk
txt <- readLines(skill, warn = FALSE)
h <- which(txt == "**File `01_import.Rmd`.**")
if (!check("the skill has the heading **File `01_import.Rmd`.** once", length(h) == 1)) finish()
o <- h + which(txt[(h + 1):length(txt)] == "````rmd")[1]
e <- o + which(txt[(o + 1):length(txt)] == "````")[1]
if (!check("a ````rmd block follows the heading", !is.na(o) && !is.na(e) && e > o + 1)) finish()
rmd <- txt[(o + 1):(e - 1)]
fb <- grep("^```\\{r functions", rmd)
if (!check("the Rmd has one functions chunk", length(fb) == 1)) finish()
fe <- fb + which(rmd[(fb + 1):length(rmd)] == "```")[1]
eval(parse(text = rmd[(fb + 1):(fe - 1)]), envir = globalenv())
if (!check("the functions chunk defines the readers", all(sapply(c("harmonise_effect", "direction_of", "read_rmats_raw",
    "read_rmats", "read_suppa", "read_dexseq", "read_edger_exon", "read_contrasts", "dexseq_fold_column", "bh_adjust"),
    exists, envir = globalenv())))) finish()

# --- 1. the five sign rules and the direction labels
check("harmonise_effect suppa 0.3 gives -0.3", harmonise_effect("suppa", 0.3) == -0.3)
check("harmonise_effect rmats 0.3 gives 0.3", harmonise_effect("rmats", 0.3) == 0.3)
check("harmonise_effect exon_dexseq 0.3 gives -0.3", harmonise_effect("exon_dexseq", 0.3) == -0.3)
check("harmonise_effect dtu 0.3 gives -0.3", harmonise_effect("dtu", 0.3) == -0.3)
check("harmonise_effect exon_edger 0.3 gives 0.3", harmonise_effect("exon_edger", 0.3) == 0.3)
check("harmonise_effect stops on an unknown tool", grepl("unknown tool: x", err_of(harmonise_effect("x", 1)), fixed = TRUE))
check("direction_of labels +, -, 0, NA", identical(direction_of(c(1, -1, 0, NA)), c("treatment_up", "control_up", "none", "none")))

# --- 2. the real rows of sign_rows.tsv
sr <- read.delim(file.path(fx, "sign_rows.tsv"), comment.char = "#", quote = "", stringsAsFactors = FALSE, colClasses = "character")
tool_of <- c(rMATS = "rmats", SUPPA2 = "suppa", DEXSeq_exon = "exon_dexseq", DEXSeq_DTU = "dtu", edgeR_exon = "exon_edger")
want <- c("+" = 1, "-" = -1, "0" = 0)
check("sign_rows.tsv has rows of all five tools", setequal(unique(sr$tool), names(tool_of)))
for (i in seq_len(nrow(sr)))
  check(sprintf("sign row %s %s: harmonised sign %s", sr$tool[i], sr$feature[i], sr$expected_sign[i]),
        sign(harmonise_effect(tool_of[[sr$tool[i]]], as.numeric(sr$raw_value[i]))) == want[[sr$expected_sign[i]]])

# --- 3. the rMATS reader: duplicate ID, quotes, one row per event, padj = FDR
se <- file.path(rt, "star/rmats/GBR_vs_YRI/SE.MATS.JC.txt")
raw_q <- read.delim(se, quote = "\"", stringsAsFactors = FALSE, colClasses = "character")   # quotes removed by read.delim
rr <- read_rmats_raw(se)
check("rMATS raw table has one ID column and the second one renamed ID_2", sum(names(rr) == "ID") == 1 && "ID_2" %in% names(rr))
check("rMATS raw table: ID and ID_2 hold the same values", identical(rr$ID, rr$ID_2))
r1 <- read_rmats(se, type = "SE", contrast = "GBR_vs_YRI")
check("rMATS harmonised table has the contract columns", identical(names(r1), COLS))
check("rMATS gene_id has no quote characters", !any(grepl("\"", r1$gene_id)))
check("rMATS gene_symbol has no quote characters", !any(grepl("\"", r1$gene_symbol)))
check("rMATS gene_id equals the unquoted GeneID", identical(r1$gene_id, raw_q$GeneID))
check("rMATS gene_symbol equals the unquoted geneSymbol", identical(r1$gene_symbol, raw_q$geneSymbol))
check("rMATS one row per event (202 events, unique feature_id)", nrow(r1) == 202 && nrow(r1) == nrow(raw_q) && !anyDuplicated(r1$feature_id))
check("rMATS padj equals the file's FDR", isTRUE(all.equal(r1$padj, as.numeric(raw_q$FDR))))
check("rMATS effect equals IncLevelDifference (no flip)", isTRUE(all.equal(r1$effect, as.numeric(raw_q$IncLevelDifference))))
check("rMATS effect_scale is dPSI", all(r1$effect_scale == "dPSI"))
sup <- function(a, b) sapply(seq_along(a), function(i) mean(as.numeric(strsplit(a[i], ",")[[1]]) + as.numeric(strsplit(b[i], ",")[[1]])))
ns <- pmin(sup(raw_q$IJC_SAMPLE_1, raw_q$SJC_SAMPLE_1), sup(raw_q$IJC_SAMPLE_2, raw_q$SJC_SAMPLE_2))
check("rMATS n_support is the smaller group mean of IJC + SJC", isTRUE(all.equal(r1$n_support, ns)))

# --- 4. SUPPA2: nan rows kept with effect NA, BH on the tested p-values
lf <- file.path(rt, "salmon/suppa/diffsplice/per_local_event/local_GBR-YRI.dpsi")
lt <- do.call(rbind, strsplit(readLines(lf)[-1], "\t", fixed = TRUE))
dpsi <- suppressWarnings(as.numeric(ifelse(lt[, 2] == "nan", NA, lt[, 2]))); pv <- suppressWarnings(as.numeric(ifelse(lt[, 3] == "nan", NA, lt[, 3])))
s1 <- read_suppa(lf, kind = "local", treatment = "GBR", control = "YRI", contrast = "GBR_vs_YRI")
check("SUPPA2 local: every row kept (3919)", nrow(s1) == nrow(lt) && nrow(s1) == 3919)
check("SUPPA2 local: a nan row has effect NA and direction none", any(is.na(s1$effect)) && all(s1$direction[is.na(s1$effect)] == "none"))
check("SUPPA2 local: effect NA exactly where dPSI is nan (2823 rows)", identical(is.na(s1$effect), is.na(dpsi)) && sum(is.na(dpsi)) == 2823)
check("SUPPA2 local: gene_id is the text before the first ;", identical(s1$gene_id, sub(";.*$", "", lt[, 1])))
tested <- !is.na(dpsi) & !is.na(pv)
bh <- rep(NA_real_, length(pv)); bh[tested] <- p.adjust(pv[tested], method = "BH")
check("SUPPA2 padj is BH of the p-val column over the tested rows", isTRUE(all.equal(s1$padj, bh)))
check("SUPPA2 padj differs from the raw p-value (BH applied)", any(abs(s1$padj - pv) > 1e-12, na.rm = TRUE))
check("bh_adjust equals p.adjust BH on the non-NA values", isTRUE(all.equal(bh_adjust(c(0.01, NA, 0.04, 0.03)), c(p.adjust(c(0.01, 0.04, 0.03), "BH")[1], NA, p.adjust(c(0.01, 0.04, 0.03), "BH")[2:3]))))

# --- 5. header-only files: zero rows, the full column set
ho <- function(src, name) { f <- file.path(work, name); writeLines(readLines(src, n = 1), f); f }
z1 <- read_rmats(ho(se, "SE.MATS.JC.txt"), type = "SE", contrast = "GBR_vs_YRI")
check("header-only rMATS file gives zero rows and the contract columns", nrow(z1) == 0 && identical(names(z1), COLS))
z2 <- read_suppa(ho(lf, "local_GBR-YRI.dpsi"), kind = "local", treatment = "GBR", control = "YRI", contrast = "GBR_vs_YRI")
check("header-only SUPPA2 file gives zero rows and the contract columns", nrow(z2) == 0 && identical(names(z2), COLS))
dx <- file.path(rt, "star/dexseq_exon/results/DEXSeqResults.GBR_vs_YRI.csv")
z3 <- read_dexseq(ho(dx, "DEXSeqResults.GBR_vs_YRI.csv"), tool = "exon_dexseq", treatment = "GBR", control = "YRI", contrast = "GBR_vs_YRI")
check("header-only DEXSeq file gives zero rows and the contract columns", nrow(z3) == 0 && identical(names(z3), COLS))
ed <- file.path(rt, "star/edger/contrast_GBR_vs_YRI.usage.exon.csv")
z4 <- read_edger_exon(ho(ed, "contrast_GBR_vs_YRI.usage.exon.csv"), contrast = "GBR_vs_YRI")
check("header-only edgeR file gives zero rows and the contract columns", nrow(z4) == 0 && identical(names(z4), COLS))

# --- 6. a missing column stops with the file and the column named
mc <- file.path(work, "renamed_SE.MATS.JC.txt"); l <- readLines(se); l[1] <- sub("IncLevelDifference", "IncLevelDiff_renamed", l[1]); writeLines(l, mc)
m6 <- err_of(read_rmats(mc, type = "SE", contrast = "GBR_vs_YRI"))
check("missing column stop names the file and IncLevelDifference", grepl(mc, m6, fixed = TRUE) && grepl("IncLevelDifference", m6, fixed = TRUE))
ms <- file.path(work, "renamed_local.dpsi"); l <- readLines(lf); l[1] <- sub("_p-val", "_pvalue", l[1]); writeLines(l, ms)
m6b <- err_of(read_suppa(ms, kind = "local", treatment = "GBR", control = "YRI", contrast = "GBR_vs_YRI"))
check("missing SUPPA2 column stop names the file and local_GBR-local_YRI_p-val", grepl(ms, m6b, fixed = TRUE) && grepl("local_GBR-local_YRI_p-val", m6b, fixed = TRUE))
me <- file.path(work, "renamed_edger.csv"); l <- readLines(ed); l[1] <- sub("logFC", "log_FC", l[1]); writeLines(l, me)
m6c <- err_of(read_edger_exon(me, contrast = "GBR_vs_YRI"))
check("missing edgeR column stop names the file and logFC", grepl(me, m6c, fixed = TRUE) && grepl("logFC", m6c, fixed = TRUE))

# --- 7. SYNTHETIC contrast with - and _ in the group names: treatment GBR-1, control YRI_2
SYNTHETIC <- file.path(work, "SYNTHETIC_contrasts.csv")
writeLines(c("contrast,treatment,control", "GBR-1_vs_YRI_2,GBR-1,YRI_2"), SYNTHETIC)
sc <- read_contrasts(SYNTHETIC)
check("SYNTHETIC contrasts: treatment GBR-1 and control YRI_2 read whole", identical(sc$treatment, "GBR-1") && identical(sc$control, "YRI_2") && identical(sc$contrast, "GBR-1_vs_YRI_2"))
check("SYNTHETIC: the DEXSeq column is log2fold_YRI_2_GBR-1", identical(dexseq_fold_column(sc$treatment, sc$control), "log2fold_YRI_2_GBR-1"))
sd <- file.path(work, "SYNTHETIC_DEXSeqResults.csv"); l <- readLines(dx); l[1] <- sub("log2fold_YRI_GBR", "log2fold_YRI_2_GBR-1", l[1], fixed = TRUE); writeLines(l, sd)
d7 <- read_dexseq(sd, tool = "exon_dexseq", treatment = "GBR-1", control = "YRI_2", contrast = "GBR-1_vs_YRI_2")
check("SYNTHETIC: the DEXSeq reader finds log2fold_YRI_2_GBR-1", nrow(d7) == length(l) - 1)
m7 <- err_of(read_dexseq(dx, tool = "exon_dexseq", treatment = "GBR-1", control = "YRI_2", contrast = "GBR-1_vs_YRI_2"))
check("SYNTHETIC: a missing log2fold_YRI_2_GBR-1 stops naming the file and that column", grepl(dx, m7, fixed = TRUE) && grepl("log2fold_YRI_2_GBR-1", m7, fixed = TRUE))

# --- 8. render the whole Rmd on the fixture tree; the RDS files against an independent recomputation from the raw tables
rd <- file.path(work, "report"); dir.create(rd)
writeLines(rmd, file.path(rd, "01_import.Rmd"))
pr <- list(outdir = rt, report_dir = rd, contrasts = file.path(rt, "261003_a_contrasts.csv"),
           samplesheet = file.path(rt, "261003_a_samplesheet.csv"), modules = "rmats suppa dexseq_exon edger_exon dtu",
           rmats_type = "JC", fdr = 0.05, dpsi = 0.1, min_reads = 10, events = TRUE, cross = TRUE, topgene = TRUE,
           enrich = TRUE, organism = "human", gtf = "", enrich_dbs = "GO", paired_rmats = FALSE, paired_suppa = TRUE)
rm_ok <- tryCatch({ rmarkdown::render(file.path(rd, "01_import.Rmd"), output_dir = rd, knit_root_dir = rd, params = pr,
                                      quiet = TRUE, envir = new.env()); TRUE },
                  error = function(e) { cat("render error:", conditionMessage(e), "\n"); FALSE })
check("01_import.Rmd renders on the fixtures", rm_ok && file.exists(file.path(rd, "01_import.html")))
if (!rm_ok) finish()
rds <- function(tool) file.path(rd, "rds", sprintf("GBR_vs_YRI.%s.rds", tool))
for (tool in c("rmats", "suppa", "exon_dexseq", "exon_edger", "dtu"))
  check(sprintf("rds %s exists with the contract columns", tool), file.exists(rds(tool)) && identical(names(readRDS(rds(tool))), COLS))
qc <- readRDS(file.path(rd, "rds", "GBR_vs_YRI.qc.rds"))
check("qc$suppa_untested counts the nan rows of the local file (2823)", identical(as.integer(qc$suppa_untested[["local"]]), 2823L))
check("qc$suppa_untested counts the nan rows of the isoform file (3173)", identical(as.integer(qc$suppa_untested[["isoform"]]), 3173L))
# expected direction from the raw value by the sign rule, written out here (not through harmonise_effect)
dir_raw <- function(v, flip) { v <- if (flip) -v else v; ifelse(is.na(v) | v == 0, "none", ifelse(v > 0, "treatment_up", "control_up")) }
cmp_dir <- function(tool, got, ids, exp) {
  m <- match(got$feature_id, ids)
  check(sprintf("rds %s: every row maps to a raw row", tool), !anyNA(m))
  check(sprintf("rds %s direction equals the raw sign rule, row by row", tool), identical(got$direction, exp[m]))
  check(sprintf("rds %s direction counts equal the raw recomputation (both directions present)", tool),
        identical(as.vector(table(factor(got$direction, c("treatment_up", "control_up", "none")))),
                  as.vector(table(factor(exp[m], c("treatment_up", "control_up", "none"))))) &&
        all(c("treatment_up", "control_up") %in% got$direction))
}
# rMATS: every type of the JC fixtures, events with the smaller group mean of IJC + SJC below 10 removed
g <- readRDS(rds("rmats")); ids <- character(0); ex <- character(0); n_low <- 0; n_all <- 0
for (ty in c("SE", "A5SS", "A3SS", "MXE", "RI")) {
  r <- read.delim(file.path(rt, sprintf("star/rmats/GBR_vs_YRI/%s.MATS.JC.txt", ty)), quote = "\"", stringsAsFactors = FALSE, colClasses = "character")
  s <- pmin(sup(r$IJC_SAMPLE_1, r$SJC_SAMPLE_1), sup(r$IJC_SAMPLE_2, r$SJC_SAMPLE_2)); keep <- s >= 10
  n_low <- n_low + sum(!keep); n_all <- n_all + nrow(r)
  ids <- c(ids, paste0(ty, ":", r$ID[keep])); ex <- c(ex, dir_raw(as.numeric(r$IncLevelDifference[keep]), FALSE))
}
check("rds rmats: rows = raw events minus the low-support events", nrow(g) == n_all - n_low && n_low > 0)
check("qc rmats_low_support equals the raw recomputation", identical(as.integer(qc$rmats_low_support), as.integer(n_low)))
check("rds rmats: every row has n_support >= min_reads", all(g$n_support >= 10))
cmp_dir("rmats", g, ids, ex)
# SUPPA2: local events and isoforms
g <- readRDS(rds("suppa")); itf <- file.path(rt, "salmon/suppa/diffsplice/per_isoform/transcript_GBR-YRI.dpsi")
it <- do.call(rbind, strsplit(readLines(itf)[-1], "\t", fixed = TRUE))
ids <- c(sub("^[^;]*;", "", lt[, 1]), sub("^[^;]*;", "", it[, 1]))
ex <- c(dir_raw(dpsi, TRUE), dir_raw(suppressWarnings(as.numeric(ifelse(it[, 2] == "nan", NA, it[, 2]))), TRUE))
check("rds suppa: local and isoform rows all kept", nrow(g) == nrow(lt) + nrow(it))
cmp_dir("suppa", g, ids, ex)
# DEXSeq exon usage and DTU: flip of log2fold_YRI_GBR
r <- read.csv(dx, stringsAsFactors = FALSE, check.names = FALSE)
cmp_dir("exon_dexseq", readRDS(rds("exon_dexseq")), paste0(r$groupID, ":", r$featureID), dir_raw(r$log2fold_YRI_GBR, TRUE))
r <- read.delim(file.path(rt, "salmon/dexseq_dtu/results/dexseq/DEXSeqResults.GBR-YRI.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
cmp_dir("dtu", readRDS(rds("dtu")), r$featureID, dir_raw(r$log2fold_YRI_GBR, TRUE))
# edgeR exon usage: logFC kept
r <- read.csv(ed, stringsAsFactors = FALSE, check.names = FALSE)
cmp_dir("exon_edger", readRDS(rds("exon_edger")), paste0(r$Geneid, ":", r$Chr, ":", r$Start, "-", r$End), dir_raw(r$logFC, FALSE))
html <- paste(readLines(file.path(rd, "01_import.html"), warn = FALSE), collapse = "\n")
check("the report prints the sign rule of SUPPA2 and the BH note", grepl("control minus treatment", html, fixed = TRUE) && grepl("Benjamini-Hochberg", html, fixed = TRUE))
finish()
