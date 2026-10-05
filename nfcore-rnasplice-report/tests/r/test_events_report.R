# Usage (inside the R image, on a compute node only; tests/r/run_r_tests.sh does this through sbatch):
#   Rscript test_events_report.R <skill.md> <tests/fixtures dir>
# Extracts the Rmd blocks under **File `02_rmats.Rmd`.** and **File `03_suppa.Rmd`.** of the skill, evaluates their `functions`
# chunks (counting, thresholds, wording of the directions, the PSI table), then renders 01_import.Rmd, 02_rmats.Rmd and
# 03_suppa.Rmd on the fixture tree and checks the tables and figures against an independent recomputation from the raw tables.
# The data frame SYNTHETIC holds known effects, signs and flags: it tests the counting logic only, never the signs of a real tool.
# Prints one "ok" or "FAIL:" line per check, then "R TESTS PASS (<n>)" or "R TESTS FAIL (<failed>/<n>)" (exit 1).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: test_events_report.R <skill.md> <fixtures dir>")
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
work <- tempfile("rsr_events_"); dir.create(work)
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
rmd2 <- get_rmd("02_rmats.Rmd"); rmd3 <- get_rmd("03_suppa.Rmd")
E2 <- functions_env(rmd2, "02_rmats.Rmd"); E3 <- functions_env(rmd3, "03_suppa.Rmd")
need_fun <- function(e, nm) all(vapply(nm, exists, logical(1), envir = e, inherits = FALSE))
if (!check("02 defines its functions", need_fun(E2, c("sig_flag", "split_counts", "direction_wording", "nearest_events", "top_events",
    "low_support_report", "none_message", "find_rmats_dir", "read_inclevels", "psi_long", "psi_legend", "psi_plot", "volcano_plot")))) finish()
if (!check("03 defines its functions", need_fun(E3, c("sig_flag", "split_counts", "direction_wording", "nearest_events", "top_events",
    "none_message", "suppa_event_type", "suppa_group", "volcano_plot")))) finish()
txt_of <- function(rmd) paste(rmd, collapse = "\n")
check("02 states the sign rule before its first function",
      regexpr("Every effect is stated as **treatment minus control**", txt_of(rmd2), fixed = TRUE) < regexpr("```{r functions", txt_of(rmd2), fixed = TRUE))
check("03 states the sign rule before its first function",
      regexpr("Every effect is stated as **treatment minus control**", txt_of(rmd3), fixed = TRUE) < regexpr("```{r functions", txt_of(rmd3), fixed = TRUE))

# --- 1. SYNTHETIC: known effects with signs, thresholds FDR 0.05 and |dPSI| 0.1, one low-support event
SYNTHETIC <- data.frame(
  contrast = "SYN", tool = "rmats",
  gene_id = paste0("G", 1:11), gene_symbol = paste0("g", 1:11),
  feature_type = c(rep("SE", 9), "A5SS", "A5SS"),
  feature_id = c(paste0("SE:", 1:9), "A5SS:10", "A5SS:11"),
  effect = c(0.30, -0.25, 0.05, 0.40, 0.50, 0.20, 0.10, 0.30, -0.10, -0.60, NA),
  effect_scale = "dPSI",
  padj = c(0.01, 0.02, 0.001, 0.20, 0.001, 0.04, 0.0499, 0.05, 0.02, 0.03, NA),
  direction = c("treatment_up", "control_up", "treatment_up", "treatment_up", "treatment_up", "treatment_up", "treatment_up",
                "treatment_up", "control_up", "control_up", "none"),
  n_support = 100, low_support = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE)
# SE: rows 1, 2, 6, 7 (exactly 0.10 counts), 9 (exactly -0.10 counts) are significant; 3 fails dPSI, 4 fails FDR, 5 is low support
# (would be significant), 8 has padj exactly 0.05 (not below). A5SS: row 10 is significant, row 11 has no effect and no padj.
sg <- E2$sig_flag(SYNTHETIC, 0.05, 0.1)
check("SYNTHETIC: sig_flag marks rows 1, 2, 6, 7, 9 and 10 and no other", identical(which(sg), c(1L, 2L, 6L, 7L, 9L, 10L)))
sc <- E2$split_counts(SYNTHETIC, 0.05, 0.1, c("SE", "A5SS", "MXE"))
check("SYNTHETIC: split_counts has the type, the counts and the wording of both directions as columns",
      all(c("type", "events", "low_support", "significant", "more included in treatment", "more skipped in treatment") %in% names(sc)))
check("SYNTHETIC: SE counts 9 events, 1 low support, 5 significant (3 more included in treatment, 2 more skipped in treatment)",
      identical(as.integer(sc[sc$type == "SE", c("events", "low_support", "significant", "more included in treatment", "more skipped in treatment")]),
                c(9L, 1L, 5L, 3L, 2L)))
check("SYNTHETIC: A5SS counts 2 events, 1 significant, 0 more included in treatment, 1 more skipped in treatment",
      identical(as.integer(sc[sc$type == "A5SS", c("events", "significant", "more included in treatment", "more skipped in treatment")]), c(2L, 1L, 0L, 1L)))
check("SYNTHETIC: a type without events has all counts 0", all(sc[sc$type == "MXE", -1] == 0))
check("direction_wording of a positive effect is 'more included in treatment'",
      identical(E2$direction_wording("treatment_up"), "more included in treatment"))
check("direction_wording of a negative effect is 'more skipped in treatment'",
      identical(E2$direction_wording("control_up"), "more skipped in treatment"))
check("SYNTHETIC: the direction of a positive effect is counted under 'more included in treatment' (03 as well)",
      identical(E3$direction_wording("treatment_up"), "more included in treatment") && identical(E3$direction_wording("control_up"), "more skipped in treatment"))
lr <- E2$low_support_report(SYNTHETIC, 0.05, 0.1)
check("SYNTHETIC: low_support_report gives 1 flagged and 1 of them otherwise significant", identical(as.integer(lr[c("flagged", "would_be_significant")]), c(1L, 1L)))
check("SYNTHETIC: none_message names the FDR and the dPSI", identical(E2$none_message(0.05, 0.1), "none at FDR 0.05, |dPSI| 0.1"))
check("SYNTHETIC: nearest_events gives the three smallest padj without the low-support event (SE:3, SE:1, SE:2)",
      identical(E2$nearest_events(SYNTHETIC, 3)$feature_id, c("SE:3", "SE:1", "SE:2")))
check("SYNTHETIC: top_events puts the significant events first, by padj (SE:1, SE:2, SE:9)",
      identical(E2$top_events(SYNTHETIC, 0.05, 0.1, 3)$feature_id, c("SE:1", "SE:2", "SE:9")))
syn3 <- SYNTHETIC; syn3$tool <- "suppa"; syn3$low_support <- FALSE
check("SYNTHETIC: 03 sig_flag and split_counts count the same rows", identical(which(E3$sig_flag(syn3, 0.05, 0.1)), c(1L, 2L, 5L, 6L, 7L, 9L, 10L)))
check("suppa_event_type parses the type after the FIRST ; only",
      identical(unname(E3$suppa_event_type(c("ENSG1;SE:X:1-2:3-4:+", "ENSG1;A5:X:1-2:3-4:-", "ENSG1;MX:X:1-2;3-4:+", "SE:X:1-2:3-4:+", "ENST0001"))),
                c("SE", "A5", "MX", "SE", "ENST0001")))
check("suppa_group is isoform for an isoform row and the parsed type for a local event",
      identical(unname(E3$suppa_group(data.frame(feature_type = c("SE", "isoform", "AF"), feature_id = c("SE:X:1-2:3-4:+", "ENST0001", "AF:X:1-2:3-4:+"),
                                                 stringsAsFactors = FALSE))), c("SE", "isoform", "AF")))

# --- 2. the PSI table: sample order from the sample sheet, group 1 = treatment
inc <- E2$read_inclevels(file.path(rt, "star/rmats/GBR_vs_YRI"), c("SE:3", "SE:5"), "JC")
check("read_inclevels reads the two events by ID (the second ID column does not matter)", identical(inc$feature_id, c("SE:3", "SE:5")))
check("read_inclevels returns the comma lists of IncLevel1 and IncLevel2 of SE:3", identical(inc$IncLevel1[1], "0.771,1.0") && identical(inc$IncLevel2[1], "1.0,1.0"))
SYNTHETIC_sheet <- data.frame(sample = c("S_B", "S_A", "C_B", "C_A"), condition = c("GBR", "GBR", "YRI", "YRI"), stringsAsFactors = FALSE)
pl <- E2$psi_long(inc, SYNTHETIC_sheet, "GBR", "YRI")
p3 <- pl[pl$feature_id == "SE:3", ]
check("psi_long: the first IncLevel1 value belongs to the FIRST treatment sample of the sheet (S_B), not to the first by name",
      p3$sample[1] == "S_B" && isTRUE(all.equal(p3$psi[1], 0.771)) && p3$sample[2] == "S_A")
check("psi_long: samples come in the order of the sheet (S_B S_A C_B C_A)", identical(as.character(p3$sample), c("S_B", "S_A", "C_B", "C_A")))
check("psi_long: IncLevel1 samples have the role treatment and the condition GBR, IncLevel2 samples control and YRI",
      identical(as.character(p3$role), c("treatment", "treatment", "control", "control")) && identical(as.character(p3$condition), c("GBR", "GBR", "YRI", "YRI")))
check("psi_long stops when the sheet has another number of samples than the IncLevel list",
      grepl("samples", tryCatch({ E2$psi_long(inc, SYNTHETIC_sheet[-1, ], "GBR", "YRI"); "" }, error = function(e) conditionMessage(e)), fixed = TRUE))
check("psi_legend names the treatment first and the control second", identical(E2$psi_legend("GBR", "YRI"), c("treatment (GBR)", "control (YRI)")))
lab <- data.frame(feature_id = c("SE:3", "SE:5"), gene_symbol = c("VAMP7", "TMLHE"), stringsAsFactors = FALSE)
pp <- E2$psi_plot(pl, lab, "GBR", "YRI", "GBR_vs_YRI")
check("psi_plot: the legend levels are treatment (GBR) then control (YRI)", identical(levels(pp$data$role), c("treatment (GBR)", "control (YRI)")))
check("psi_plot: the points of S_B carry the treatment label, those of C_A the control label",
      all(pp$data$role[pp$data$sample == "S_B"] == "treatment (GBR)") && all(pp$data$role[pp$data$sample == "C_A"] == "control (YRI)"))
check("psi_plot: the x axis keeps the sample sheet order", identical(levels(pp$data$sample), c("S_B", "S_A", "C_B", "C_A")))
# rMATS directory: the contrast's own, or the one with the suffix _paired; both present: the pairing switch decides
fo <- file.path(work, "fx_out"); for (d in c("star/rmats/X_paired", "star/rmats/Y", "star/rmats/Y_paired", "star/rmats/Z"))
  dir.create(file.path(fo, d), recursive = TRUE)
for (d in c("X_paired", "Y", "Y_paired")) file.copy(file.path(rt, "star/rmats/GBR_vs_YRI/SE.MATS.JC.txt"), file.path(fo, "star/rmats", d, "SE.MATS.JC.txt"))
check("find_rmats_dir finds the _paired directory when it is the only one", identical(E2$find_rmats_dir(fo, "X", FALSE, "JC"), file.path(fo, "star/rmats/X_paired")))
check("find_rmats_dir with both: unpaired run takes the plain directory, paired run the _paired one",
      identical(E2$find_rmats_dir(fo, "Y", FALSE, "JC"), file.path(fo, "star/rmats/Y")) && identical(E2$find_rmats_dir(fo, "Y", TRUE, "JC"), file.path(fo, "star/rmats/Y_paired")))
check("find_rmats_dir gives NA when no table of the type exists", is.na(E2$find_rmats_dir(fo, "Z", FALSE, "JC")) && is.na(E2$find_rmats_dir(fo, "Y", FALSE, "JCEC")))

# --- 3. render 01, 02 and 03 on the fixture tree
base_params <- list(outdir = rt, report_dir = "", contrasts = file.path(rt, "261003_a_contrasts.csv"),
                    samplesheet = file.path(rt, "261003_a_samplesheet.csv"), modules = "rmats suppa", rmats_type = "JC",
                    fdr = 0.05, dpsi = 0.1, min_reads = 10, events = TRUE, cross = TRUE, topgene = TRUE, enrich = TRUE,
                    organism = "human", gtf = "", enrich_dbs = "GO", paired_rmats = FALSE, paired_suppa = TRUE)
render <- function(rd, rmd_name, ...) {
  pr <- utils::modifyList(base_params, c(list(report_dir = rd), list(...)))
  tryCatch({ rmarkdown::render(file.path(rd, rmd_name), output_dir = rd, knit_root_dir = rd, params = pr, quiet = TRUE, envir = new.env()); TRUE },
           error = function(e) { cat("render error:", conditionMessage(e), "\n"); FALSE })
}
rmd1 <- get_rmd("01_import.Rmd")
new_dir <- function(nm) { d <- file.path(work, nm); dir.create(d); writeLines(rmd1, file.path(d, "01_import.Rmd"))
  writeLines(rmd2, file.path(d, "02_rmats.Rmd")); writeLines(rmd3, file.path(d, "03_suppa.Rmd")); d }
# the HTML is wrapped by pandoc: runs of white space (also line breaks) count as one space
html_of <- function(rd, nm) gsub("[[:space:]]+", " ", paste(readLines(file.path(rd, nm), warn = FALSE), collapse = " "))
has <- function(html, s) grepl(s, html, fixed = TRUE)
num <- function(x) suppressWarnings(as.numeric(ifelse(x %in% c("nan", "NA", ""), NA, x)))
C <- "GBR_vs_YRI"

# --- raw recomputation of rMATS, per type: events, low support, significant, directions (no code of the skill involved)
sup <- function(a, b) sapply(seq_along(a), function(i) mean(as.numeric(strsplit(a[i], ",")[[1]]) + as.numeric(strsplit(b[i], ",")[[1]])))
TYPES <- c("SE", "A5SS", "A3SS", "MXE", "RI")
raw_rmats <- function(fdr, dpsi, min_reads) {
  out <- list()
  for (ty in TYPES) {
    r <- read.delim(file.path(rt, sprintf("star/rmats/GBR_vs_YRI/%s.MATS.JC.txt", ty)), quote = "\"", stringsAsFactors = FALSE, colClasses = "character")
    s <- pmin(sup(r$IJC_SAMPLE_1, r$SJC_SAMPLE_1), sup(r$IJC_SAMPLE_2, r$SJC_SAMPLE_2)); lo <- s < min_reads
    fd <- num(r$FDR); d <- num(r$IncLevelDifference)
    sig <- !is.na(fd) & !is.na(d) & fd < fdr & abs(d) >= dpsi & !lo
    would <- !is.na(fd) & !is.na(d) & fd < fdr & abs(d) >= dpsi & lo
    out[[ty]] <- list(events = nrow(r), low = sum(lo), sig = sum(sig), up = sum(sig & d > 0), down = sum(sig & d < 0), would = sum(would),
                      ids = paste0(ty, ":", r$ID[sig]))
  }
  out
}
sum_of <- function(x, f) sum(vapply(x, function(y) y[[f]], numeric(1)))
nonempty <- function(f) file.exists(f) && file.size(f) > 0

# 3a. defaults: min_reads 10 leaves no significant rMATS event (all the ones passing the thresholds are low support)
rd10 <- new_dir("rd10")
ok <- render(rd10, "01_import.Rmd")
check("01_import.Rmd renders (rmats and suppa) for the events reports", ok && file.exists(file.path(rd10, "rds", paste0(C, ".rmats.rds"))))
if (!ok) finish()
ok2 <- render(rd10, "02_rmats.Rmd")
check("02_rmats.Rmd renders on the fixtures with the default thresholds", ok2 && file.exists(file.path(rd10, "02_rmats.html")))
if (ok2) {
  h <- html_of(rd10, "02_rmats.html"); ex <- raw_rmats(0.05, 0.1, 10)
  check("raw recomputation: no significant rMATS event at min_reads 10, and some at min_reads 0 (the test data)",
        sum_of(ex, "sig") == 0 && sum_of(raw_rmats(0.05, 0.1, 0), "sig") > 0)
  check("02 prints the sign rule at the top: treatment minus control, more included in treatment",
        has(h, "treatment minus control") && has(h, "more included in treatment") && has(h, "more skipped in treatment"))
  check("02 prints 'none at FDR 0.05, |dPSI| 0.1' once for each of the five types", length(gregexpr("none at FDR 0.05, |dPSI| 0.1", h, fixed = TRUE)[[1]]) == 5)
  check("02 prints the three nearest events under each empty result", length(gregexpr("three nearest events by adjusted p-value", h, fixed = TRUE)[[1]]) == 5)
  check("02 reports the flagged low_support events and how many would otherwise be significant",
        has(h, sprintf("flagged low_support (read support below 10): %d; of these %d would be significant", sum_of(ex, "low"), sum_of(ex, "would"))) && sum_of(ex, "would") > 0)
  cn <- read.csv(file.path(rd10, "tables", paste0(C, ".rmats.counts.csv")), check.names = FALSE, stringsAsFactors = FALSE)
  check("the count table has the five types and zero significant events", identical(cn$type, TYPES) && all(cn$significant == 0))
  check("the count table's low_support column equals the raw recomputation", identical(as.integer(cn$low_support), as.integer(sapply(ex, function(y) y$low))))
  check("the significant tables of the five types exist and hold no row",
        all(sapply(TYPES, function(ty) { f <- file.path(rd10, "tables", sprintf("%s.rmats.%s.significant.csv", C, ty)); file.exists(f) && nrow(read.csv(f)) == 0 })))
  check("the volcano of every type and the top-event PSI plot exist and are not empty",
        all(sapply(TYPES, function(ty) nonempty(file.path(rd10, "figures", sprintf("%s.rmats.volcano.%s.png", C, ty))))) &&
        nonempty(file.path(rd10, "figures", paste0(C, ".rmats.psi_top.png"))))
  check("02 says the rMATS statistics are unpaired for this run", has(h, "unpaired statistics"))
}
ok3 <- render(rd10, "03_suppa.Rmd")
check("03_suppa.Rmd renders on the fixtures with the default thresholds", ok3 && file.exists(file.path(rd10, "03_suppa.html")))

# 3b. min_reads 0: the counts equal the raw recomputation for two threshold pairs (a literal threshold in the Rmd would fail the second)
rd0 <- new_dir("rd0")
ok <- render(rd0, "01_import.Rmd", min_reads = 0)
check("01_import.Rmd renders with min_reads 0", ok)
if (!ok) finish()
for (th in list(c(0.05, 0.1), c(0.9, 0.05))) {
  fdr <- th[1]; dp <- th[2]; tag <- sprintf("FDR %s, dPSI %s", fdr, dp)
  ok2 <- render(rd0, "02_rmats.Rmd", fdr = fdr, dpsi = dp, min_reads = 0)
  check(sprintf("02_rmats.Rmd renders with min_reads 0 at %s", tag), ok2)
  if (!ok2) next
  h <- html_of(rd0, "02_rmats.html"); ex <- raw_rmats(fdr, dp, 0)
  cn <- read.csv(file.path(rd0, "tables", paste0(C, ".rmats.counts.csv")), check.names = FALSE, stringsAsFactors = FALSE)
  check(sprintf("rMATS significant counts per type equal the raw recomputation at %s (%s)", tag, paste(sapply(ex, function(y) y$sig), collapse = ",")),
        identical(cn$type, TYPES) && identical(as.integer(cn$significant), as.integer(sapply(ex, function(y) y$sig))) && sum_of(ex, "sig") > 0)
  check(sprintf("rMATS direction split equals the raw recomputation at %s", tag),
        identical(as.integer(cn[["more included in treatment"]]), as.integer(sapply(ex, function(y) y$up))) &&
        identical(as.integer(cn[["more skipped in treatment"]]), as.integer(sapply(ex, function(y) y$down))))
  check(sprintf("the significant tables hold exactly the raw significant events at %s", tag),
        all(sapply(TYPES, function(ty) { f <- file.path(rd0, "tables", sprintf("%s.rmats.%s.significant.csv", C, ty))
          file.exists(f) && setequal(read.csv(f, stringsAsFactors = FALSE)$feature_id, ex[[ty]]$ids) && nrow(read.csv(f)) == ex[[ty]]$sig })))
  check(sprintf("no flagged event at min_reads 0 and no empty-result message at %s", tag),
        all(cn$low_support == 0) && !has(h, "none at FDR") && has(h, sprintf("FDR %s, abs(dPSI) %s", fdr, dp)))
  check(sprintf("volcano and PSI figures exist and are not empty at %s", tag),
        all(sapply(TYPES, function(ty) nonempty(file.path(rd0, "figures", sprintf("%s.rmats.volcano.%s.png", C, ty))))) &&
        nonempty(file.path(rd0, "figures", paste0(C, ".rmats.psi_top.png"))) && nonempty(file.path(rd0, "figures", paste0(C, ".rmats.effect_hist.png"))))
}

# 3c. SUPPA2: counts per event type and for isoforms equal the raw recomputation (BH per file over the tested rows; effect = minus dPSI)
sdir <- file.path(rt, "salmon/suppa/diffsplice")
raw_suppa <- function(fdr, dpsi) {
  rd_ <- function(f) { t <- do.call(rbind, strsplit(readLines(f)[-1], "\t", fixed = TRUE)); data.frame(lab = t[, 1], d = num(t[, 2]), p = num(t[, 3]), stringsAsFactors = FALSE) }
  one <- function(t, grp) {
    ok <- !is.na(t$d) & !is.na(t$p); q <- rep(NA_real_, nrow(t)); q[ok] <- p.adjust(t$p[ok], method = "BH")
    sig <- !is.na(q) & q < fdr & abs(t$d) >= dpsi
    lapply(split(seq_len(nrow(t)), grp), function(i) list(events = length(i), untested = sum(is.na(t$d[i])), tested = sum(!is.na(q[i])), sig = sum(sig[i]),
           up = sum(sig[i] & t$d[i] < 0), down = sum(sig[i] & t$d[i] > 0), ids = sub("^[^;]*;", "", t$lab[i][sig[i]])))
  }
  l <- rd_(file.path(sdir, "per_local_event/local_GBR-YRI.dpsi")); i <- rd_(file.path(sdir, "per_isoform/transcript_GBR-YRI.dpsi"))
  c(one(l, sub(":.*$", "", sub("^[^;]*;", "", l$lab))), list(isoform = one(i, rep("isoform", nrow(i)))$isoform))
}
STY <- c("SE", "A5", "A3", "MX", "RI", "AF", "AL", "isoform")
for (th in list(c(0.05, 0.1), c(0.9, 0.05))) {
  fdr <- th[1]; dp <- th[2]; tag <- sprintf("FDR %s, dPSI %s", fdr, dp)
  rdx <- if (fdr == 0.05) rd10 else rd0
  ok3 <- render(rdx, "03_suppa.Rmd", fdr = fdr, dpsi = dp, min_reads = if (fdr == 0.05) 10 else 0)
  check(sprintf("03_suppa.Rmd renders at %s", tag), ok3)
  if (!ok3) next
  h <- html_of(rdx, "03_suppa.html"); ex <- raw_suppa(fdr, dp)
  check(sprintf("raw recomputation: the SUPPA2 groups are %s", paste(STY, collapse = " ")), setequal(names(ex), STY))
  cn <- read.csv(file.path(rdx, "tables", paste0(C, ".suppa.counts.csv")), check.names = FALSE, stringsAsFactors = FALSE)
  check(sprintf("SUPPA2 count table has one row per event type and the isoforms, in order, at %s", tag), identical(cn$type, STY))
  check(sprintf("SUPPA2 significant counts per type equal the raw recomputation at %s (%s)", tag, paste(sapply(ex[STY], function(y) y$sig), collapse = ",")),
        identical(as.integer(cn$significant), as.integer(sapply(ex[STY], function(y) y$sig))))
  check(sprintf("SUPPA2 events, untested and tested counts equal the raw recomputation at %s", tag),
        identical(as.integer(cn$events), as.integer(sapply(ex[STY], function(y) y$events))) && identical(as.integer(cn$untested), as.integer(sapply(ex[STY], function(y) y$untested))) &&
        identical(as.integer(cn$tested), as.integer(sapply(ex[STY], function(y) y$tested))))
  check(sprintf("SUPPA2 direction split equals the raw recomputation (effect is minus dPSI) at %s", tag),
        identical(as.integer(cn[["more included in treatment"]]), as.integer(sapply(ex[STY], function(y) y$up))) &&
        identical(as.integer(cn[["more skipped in treatment"]]), as.integer(sapply(ex[STY], function(y) y$down))))
  check(sprintf("SUPPA2 significant tables hold exactly the raw significant events at %s", tag),
        all(sapply(STY, function(ty) { f <- file.path(rdx, "tables", sprintf("%s.suppa.%s.significant.csv", C, ty))
          file.exists(f) && setequal(read.csv(f, stringsAsFactors = FALSE)$feature_id, ex[[ty]]$ids) && nrow(read.csv(f)) == ex[[ty]]$sig })))
  nz <- sum(sapply(ex[STY], function(y) y$events > 0 && y$sig == 0))
  check(sprintf("03 prints the empty-result message once per group without a significant event (%d) at %s", nz, tag),
        length(gregexpr("none at FDR", h, fixed = TRUE)[[1]][gregexpr("none at FDR", h, fixed = TRUE)[[1]] > 0]) == nz)
  check(sprintf("volcano of every group with tested events, the histogram and the top plot exist at %s", tag),
        all(sapply(STY[sapply(ex[STY], function(y) y$tested > 0)], function(ty) nonempty(file.path(rdx, "figures", sprintf("%s.suppa.volcano.%s.png", C, ty))))) &&
        nonempty(file.path(rdx, "figures", paste0(C, ".suppa.effect_hist.png"))) && nonempty(file.path(rdx, "figures", paste0(C, ".suppa.top.png"))))
  check(sprintf("03 states the Benjamini-Hochberg correction, the sign rule and the pairing at %s", tag),
        has(h, "Benjamini-Hochberg") && has(h, "treatment minus control") && has(h, "more included in treatment") && has(h, "paired statistics") && has(h, sprintf("FDR %s, abs(dPSI) %s", fdr, dp)))
}

# --- 4. a missing RDS is reported, not an error; events = FALSE and topgene = FALSE leave out tables and plots
rdm <- new_dir("rd_missing"); dir.create(file.path(rdm, "rds"))
okm2 <- render(rdm, "02_rmats.Rmd"); okm3 <- render(rdm, "03_suppa.Rmd")
check("02 and 03 render when the RDS of the contrast is missing", okm2 && okm3)
if (okm2 && okm3) {
  check("02 says that the rMATS RDS of the contrast was not found", has(html_of(rdm, "02_rmats.html"), paste0(C, ".rmats.rds")) && has(html_of(rdm, "02_rmats.html"), "not found"))
  check("03 says that the SUPPA2 RDS of the contrast was not found", has(html_of(rdm, "03_suppa.html"), paste0(C, ".suppa.rds")) && has(html_of(rdm, "03_suppa.html"), "not found"))
}
rdo <- new_dir("rd_off"); dir.create(file.path(rdo, "rds")); file.copy(list.files(file.path(rd0, "rds"), full.names = TRUE), file.path(rdo, "rds"))
oko2 <- render(rdo, "02_rmats.Rmd", min_reads = 0, events = FALSE, topgene = FALSE); oko3 <- render(rdo, "03_suppa.Rmd", min_reads = 0, events = FALSE, topgene = FALSE)
check("02 and 03 render with events and topgene switched off", oko2 && oko3)
check("with events and topgene off no volcano, no PSI plot, no top plot and no significant table is written",
      length(list.files(file.path(rdo, "figures"), pattern = "volcano|psi_top|\\.top\\.")) == 0 && length(list.files(file.path(rdo, "tables"), pattern = "significant")) == 0)
check("02 says the events section is switched off", oko2 && has(html_of(rdo, "02_rmats.html"), "switched off"))
finish()
