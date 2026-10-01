# Usage: Rscript test_rmd05_stops.R <generated Rmd 05> <work dir>
# Rmd 05 must stop with a clear message on bad phASER gene tables. For each case a copy of the Rmd's inputs is made under
# <work dir>/<case>, one sample's table (the third in sorted order) is damaged, and the Rmd (RESULTS_DIR redirected to the copy) is
# rendered in its own R process; the render must fail and its output must contain the expected message. Needs the Rmd's
# RESULTS_DIR to hold ase_checkpoint.rds, ase_imbalance_checkpoint.rds and phaser/. Prints one STOP TEST line per case and
# ALL STOP TESTS PASS; exits 1 otherwise. Run it in an sbatch job (R, rmarkdown).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: Rscript test_rmd05_stops.R <generated Rmd 05> <work dir>", call. = FALSE)
RMD <- normalizePath(args[1]); WORK <- normalizePath(args[2], mustWork = FALSE); dir.create(WORK, showWarnings = FALSE, recursive = TRUE)
L0 <- readLines(RMD)
rl <- grep('^RESULTS_DIR <- "', L0, value = TRUE); if (length(rl) != 1) stop("no RESULTS_DIR line in ", RMD, call. = FALSE)
OLD <- sub('^RESULTS_DIR <- "([^"]*)".*', "\\1", rl)
smp <- sort(sub("\\.gene_ae\\.txt$", "", list.files(file.path(OLD, "phaser"), pattern = "\\.gene_ae\\.txt$")))
if (length(smp) < 3) stop("need at least 3 phASER tables in ", file.path(OLD, "phaser"), call. = FALSE)
S <- smp[3]
bad <- 0
case <- function(name, damage, expect) {
  d <- file.path(WORK, name); R <- file.path(d, "results"); unlink(d, recursive = TRUE); dir.create(R, recursive = TRUE)
  file.copy(file.path(OLD, c("ase_checkpoint.rds", "ase_imbalance_checkpoint.rds")), R); file.copy(file.path(OLD, "phaser"), R, recursive = TRUE)
  damage(file.path(R, "phaser", paste0(S, ".gene_ae.txt")))
  rmd <- file.path(d, basename(RMD)); writeLines(gsub(paste0('RESULTS_DIR <- "', OLD, '"'), paste0('RESULTS_DIR <- "', R, '"'), L0, fixed = TRUE), rmd)
  out <- suppressWarnings(system2("Rscript", c("-e", shQuote(sprintf("rmarkdown::render('%s', output_dir = '%s', knit_root_dir = '%s', quiet = TRUE)", rmd, R, d))),
                                  stdout = TRUE, stderr = TRUE))
  st <- attr(out, "status"); st <- if (is.null(st)) 0 else st
  txt <- paste(out, collapse = " ")
  ok <- st != 0 && all(vapply(expect, function(e) grepl(e, txt, fixed = TRUE), logical(1)))
  cat(sprintf("STOP TEST %-14s %s: render exit %d; expected %s\n", name, if (ok) "PASS" else "FAIL", st, paste(shQuote(expect), collapse = " and ")))
  if (!ok) { bad <<- 1; cat(utils::tail(out, 15), sep = "\n") }
}
rw <- function(f, fun) { L <- readLines(f); writeLines(fun(L), f) }
case("missing_table", function(f) file.remove(f), paste("phASER gene counts missing for sample(s)", S))
case("zero_counts", function(f) rw(f, function(L) { p <- strsplit(L[-1], "\t"); L[-1] <- vapply(p, function(x) { x[5:7] <- "0"; paste(x, collapse = "\t") }, ""); L }),
     paste("no gene with haplotype reads in the phASER table of sample(s)", S))
case("header_only", function(f) rw(f, function(L) L[1]), paste("no gene with haplotype reads in the phASER table of sample(s)", S))
case("empty_count", function(f) rw(f, function(L) { x <- strsplit(L[2], "\t")[[1]]; x[5] <- ""; L[2] <- paste(x, collapse = "\t"); L }),
     c("missing counts (aCount, bCount or totalCount) in the phASER gene tables", paste0(S, " / ")))
case("wrong_columns", function(f) rw(f, function(L) sub("\t[^\t]*\t([^\t]*)$", "\t\\1", L)),
     c("unexpected columns in", "expected (columns of phaser_gene_ae at the commit pinned in Step 18)"))
cat(if (bad == 0) "ALL STOP TESTS PASS\n" else "SOME STOP TESTS FAIL\n")
quit(status = bad)
