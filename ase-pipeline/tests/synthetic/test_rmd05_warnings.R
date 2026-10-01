# Usage: Rscript test_rmd05_warnings.R <generated Rmd 05> <work dir>
# Per-sample dispersion-warning path of Rmd 05: the Rmd's code is run (knitr::purl) on a copy of its inputs in <work dir>, with
# bb_estimate_rho_trim wrapped so that it warns on the data of ONE sample (the second in sorted order). The Rmd must print the
# warning with that sample and keep it in rho_warning of that sample only. Needs the Rmd's RESULTS_DIR to hold ase_checkpoint.rds,
# ase_imbalance_checkpoint.rds and phaser/. Prints WARNING ATTRIBUTION PASS/FAIL; exits 1 on FAIL. Run it in an sbatch job (R).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: Rscript test_rmd05_warnings.R <generated Rmd 05> <work dir>", call. = FALSE)
RMD <- normalizePath(args[1]); WORK <- args[2]
L0 <- readLines(RMD)
rl <- grep('^RESULTS_DIR <- "', L0, value = TRUE); if (length(rl) != 1) stop("no RESULTS_DIR line in ", RMD, call. = FALSE)
OLD <- sub('^RESULTS_DIR <- "([^"]*)".*', "\\1", rl); NEW <- file.path(normalizePath(WORK, mustWork = FALSE), "results")
unlink(NEW, recursive = TRUE); dir.create(NEW, recursive = TRUE)
file.copy(file.path(OLD, c("ase_checkpoint.rds", "ase_imbalance_checkpoint.rds")), NEW); file.copy(file.path(OLD, "phaser"), NEW, recursive = TRUE)
src <- knitr::purl(RMD, output = file.path(WORK, "rmd05_purl.R"), quiet = TRUE)
L <- readLines(src)
L <- gsub(paste0('RESULTS_DIR <- "', OLD, '"'), paste0('RESULTS_DIR <- "', NEW, '"'), L, fixed = TRUE)
e <- grep("^tst <- genes\\$totalCount >= MIN_DEPTH", L); if (length(e) != 1) stop("the test chunk of Rmd 05 was not found", call. = FALSE)
inj <- c("WARN_SAMPLE <- sort(unique(genes$sample))[2]",
         "bb_estimate_rho_trim_orig <- bb_estimate_rho_trim",
         "bb_estimate_rho_trim <- function(x, n, ...) {   # warns on the data of WARN_SAMPLE only (same rows as hap_gene_test uses)",
         "  if (isTRUE(all.equal(sum(n), sum(genes$totalCount[genes$sample == WARN_SAMPLE & tst])))) warning('injected test warning')",
         "  bb_estimate_rho_trim_orig(x, n, ...) }")
L <- append(L, inj, after = e)
writeLines(L, file.path(WORK, "rmd05_warn.R"))
out <- suppressWarnings(system2("Rscript", file.path(WORK, "rmd05_warn.R"), stdout = TRUE, stderr = TRUE))
st <- attr(out, "status"); st <- if (is.null(st)) 0 else st
cat(grep("WARNING \\(rho estimation\\)", out, value = TRUE), sep = "\n")
s <- if (st == 0) readRDS(file.path(NEW, "ase_phaser_checkpoint.rds"))$samples else NULL
ws <- sort(unique(s$sample))[2]
if (!is.null(s)) print(s[, c("sample", "n_genes", "rho_used", "rho_warning")])
ok <- st == 0 && !is.null(s) && identical(s$rho_warning == "", s$sample != ws) && grepl("injected test warning", s$rho_warning[s$sample == ws]) &&
  any(grepl(paste0("WARNING \\(rho estimation\\), sample ", ws, " : injected test warning"), out))
if (st != 0) cat(utils::tail(out, 20), sep = "\n")
cat(if (ok) "WARNING ATTRIBUTION PASS\n" else "WARNING ATTRIBUTION FAIL\n")
quit(status = if (ok) 0 else 1)
