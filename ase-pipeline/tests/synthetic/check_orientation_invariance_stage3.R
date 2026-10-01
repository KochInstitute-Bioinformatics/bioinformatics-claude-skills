# Usage: Rscript check_orientation_invariance_stage3.R <RESULTS_DIR a> <RESULTS_DIR b>
# Two Rmd 05 renders of the same unphased counts in which some genes have aCount and bCount swapped (for example the emulated
# phASER tables and the same tables with the emulated flips undone): the haplotype labels are arbitrary, so the test must give
# identical p, padj, sig, major_frac and rho_used. Prints ORIENTATION INVARIANCE PASS/FAIL; exits 1 on FAIL or if nothing was swapped.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: Rscript check_orientation_invariance_stage3.R <RESULTS_DIR a> <RESULTS_DIR b>", call. = FALSE)
a <- readRDS(file.path(args[1], "ase_phaser_checkpoint.rds"))$gene; b <- readRDS(file.path(args[2], "ase_phaser_checkpoint.rds"))$gene
same_rows <- identical(a$gene_id, b$gene_id) && identical(a$sample, b$sample)
swapped <- if (same_rows) sum(a$aCount != b$aCount) else NA
ok <- same_rows && isTRUE(swapped > 0) && identical(a$totalCount, b$totalCount) && isTRUE(all.equal(a$p, b$p)) &&
  isTRUE(all.equal(a$padj, b$padj)) && identical(a$sig, b$sig) && isTRUE(all.equal(a$major_frac, b$major_frac)) && isTRUE(all.equal(a$rho_used, b$rho_used))
cat("rows", nrow(a), "; rows with aCount and bCount swapped", swapped, "\n")
cat(if (ok) "ORIENTATION INVARIANCE PASS\n" else "ORIENTATION INVARIANCE FAIL\n")
quit(status = if (ok) 0 else 1)
