# Usage: Rscript prove_evaluate_stage3.R <phased RESULTS_DIR> <unphased RESULTS_DIR> <truth dir (outbred_phase)> <work dir>
# Non-vacuity proofs for evaluate_stage3.R (rmd mode). The two RESULTS_DIRs are Rmd 05 renders on the direct phASER tables
# (phased: direct_gene_ae_phased with PHASED_GT 1; unphased: direct_gene_ae_unphased with PHASED_GT 0); each must hold
# ase_phaser_checkpoint.rds, summary_numbers_phaser.tsv, the Rmd 05 HTML and phaser/<sample>.gene_ae.txt.
# 1. Both untampered results must make the evaluator exit 0.
# 2. Each tampered copy (written under <work dir>/tamper/<name>) must make it exit 1 with the FAIL line of the named gate tag.
# Every gate tag of evaluate_stage3.R ("[tag]") appears below at least once (check_skill.sh checks this).
# Ends with "ALL TAMPER PROOFS OK" and exit 0, or "SOME TAMPER PROOFS BAD" and exit 1. Run it in an sbatch job (R).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: Rscript prove_evaluate_stage3.R <phased RESULTS_DIR> <unphased RESULTS_DIR> <truth dir> <work dir>", call. = FALSE)
base <- c(phased = args[1], unphased = args[2]); TRUTH <- args[3]; WORK <- args[4]
self <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
EV <- file.path(dirname(normalizePath(self)), "evaluate_stage3.R")
if (!file.exists(EV)) stop("evaluate_stage3.R not found next to this script: ", EV, call. = FALSE)
tg <- read.delim(file.path(TRUTH, "truth_genes.tsv"), stringsAsFactors = FALSE)
cls <- setNames(tg$class, tg$gene_id)
bad <- 0
evaluate <- function(dir, gt) {
  out <- suppressWarnings(system2("Rscript", c(EV, "rmd", dir, TRUTH, gt), stdout = TRUE, stderr = TRUE))
  st <- attr(out, "status"); list(out = out, status = if (is.null(st)) 0 else st)
}
# 1. untampered results: exit 0, and the two_block numbers (observed, oracle) for the bound tampers below
tb <- list()
for (gt in names(base)) {
  r <- evaluate(base[[gt]], gt)
  ln <- grep("^PASS \\[two_block\\]", r$out, value = TRUE)
  cat(sprintf("UNTAMPERED %-8s exit %d; %s\n", gt, r$status, if (length(ln)) ln else "no PASS [two_block] line"))
  if (r$status != 0 || length(ln) != 1) { bad <- 1; next }
  tb[[gt]] <- as.integer(strsplit(sub(".*observed ([0-9]+), oracle \\(true rho [0-9.]+, input counts\\) ([0-9]+),.*", "\\1 \\2", ln), " ")[[1]])
  rho_true <- as.numeric(sub(".*oracle \\(true rho ([0-9.]+),.*", "\\1", ln))
}
if (bad) { cat("SOME TAMPER PROOFS BAD (an untampered result does not pass)\n"); quit(status = 1) }
# 2. tampered copies
run <- function(name, gt, expect, edit = function(ck) ck, edit_tsv = NULL, edit_dir = NULL, gt_arg = gt, after = NULL) {
  d <- file.path(WORK, "tamper", name); unlink(d, recursive = TRUE); dir.create(d, recursive = TRUE)
  src <- base[[gt]]
  file.copy(c(file.path(src, c("ase_phaser_checkpoint.rds", "summary_numbers_phaser.tsv")),
              list.files(src, pattern = "_05_phaser\\.html$", full.names = TRUE)), d)
  file.copy(file.path(src, "phaser"), d, recursive = TRUE)
  ck <- readRDS(file.path(d, "ase_phaser_checkpoint.rds")); ck <- edit(ck); saveRDS(ck, file.path(d, "ase_phaser_checkpoint.rds"))
  if (!is.null(edit_tsv)) {
    s <- read.delim(file.path(d, "summary_numbers_phaser.tsv"), stringsAsFactors = FALSE); s <- edit_tsv(s)
    write.table(s, file.path(d, "summary_numbers_phaser.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  if (!is.null(edit_dir)) edit_dir(d)
  r <- evaluate(d, gt_arg)
  hit <- r$out[startsWith(r$out, paste("FAIL", expect))]
  ok <- r$status == 1 && length(hit) > 0
  cat(sprintf("%s tamper %-24s exit %d; expected 'FAIL %s': %s\n", if (ok) "PROOF OK " else "PROOF BAD", name, r$status, expect,
              if (length(hit)) hit[1] else paste("not found; FAIL lines:", paste(grep("^FAIL", r$out, value = TRUE), collapse = " | "))))
  if (!ok) bad <<- 1
  if (!is.null(after)) after(ck)
}
plt <- function(g) g$sig & !cls[g$gene_id] %in% c("null", "null_linked")
drop_calls <- function(ck, cl, n) {   # make the first n significant cells of class cl not significant
  i <- which(ck$gene$sig & cls[ck$gene$gene_id] == cl)[seq_len(n)]; ck$gene$sig[i] <- FALSE; ck$gene$direction[i] <- "none"; ck }
below <- function(x) x[1] - (ceiling(0.9 * x[2]) - 1)   # calls to remove so that ceiling(0.9 x oracle) - 1 remain
run("inverted_sign", "phased", "[phased_orient]", function(ck) {
  i <- plt(ck$gene); a <- ck$gene$aCount[i]; ck$gene$aCount[i] <- ck$gene$bCount[i]; ck$gene$bCount[i] <- a
  ck$gene$hap_A_frac[i] <- 1 - ck$gene$hap_A_frac[i]; ck })
run("swapped_labels", "phased", "[phased_label]", function(ck) {
  i <- plt(ck$gene); ck$gene$direction[i] <- ifelse(ck$gene$direction[i] == "haplotype A higher", "haplotype B higher", "haplotype A higher"); ck })
run("not_gw_phased", "phased", "[phased_orient]", function(ck) { i <- which(plt(ck$gene))[1]; ck$gene$gw_phased[i] <- FALSE; ck })
run("strong_missed", "phased", "[hap_strong]", function(ck) drop_calls(ck, "hap_strong", 1))
run("two_block_below_phased", "phased", "[two_block]", function(ck) drop_calls(ck, "two_block", below(tb$phased)))
run("two_block_below_unphased", "unphased", "[two_block]", function(ck) drop_calls(ck, "two_block", below(tb$unphased)))
run("lowdepth_missed", "phased", "[hap_lowdepth]", function(ck) drop_calls(ck, "hap_lowdepth", sum(ck$gene$sig & cls[ck$gene$gene_id] == "hap_lowdepth")))
nullcalls <- function(ck) {
  i <- which(!ck$gene$sig & cls[ck$gene$gene_id] == "null")[1:7]; ck$gene$sig[i] <- TRUE
  ck$gene$direction[i] <- if (ck$constants$PHASED_GT == 1) ifelse(ck$gene$hap_A_frac[i] > 0.5, "haplotype A higher", "haplotype B higher") else
    "no direction (haplotype labels arbitrary)"; ck }
run("null_calls_phased", "phased", "[null]", nullcalls)
run("null_calls_unphased", "unphased", "[null]", nullcalls)
run("lowdepth_category", "phased", "[lowdepth_category]", function(ck) {
  i <- cls[ck$comparison$gene_id] == "hap_lowdepth" & ck$comparison$sig; ck$comparison$category[i] <- "phASER only"; ck })
run("p_below_min_depth", "phased", "[rows]", function(ck) {
  i <- which(ck$gene$totalCount < ck$constants$MIN_DEPTH)[1]; if (is.na(i)) stop("no row below MIN_DEPTH to tamper"); ck$gene$p[i] <- 0.5; ck })
run("zero_count_row", "phased", "[rows]", function(ck) {
  r <- ck$gene[1, ]; r$aCount <- 0; r$bCount <- 0; r$totalCount <- 0; r$gene_id <- "SYNZERO"; ck$gene <- rbind(ck$gene, r); ck })
run("phased_gt_flag", "phased", "[phased_gt]", function(ck) { ck$constants$PHASED_GT <- 0; ck })
run("unphased_as_phased", "unphased", "[phased_gt]", gt_arg = "phased")
run("summary_count", "phased", "[summary]", edit_tsv = function(s) { s$sig_genes[1] <- s$sig_genes[1] + 1; s })
run("summary_row_missing", "phased", "[summary]", edit_tsv = function(s) s[-1, ])
run("unphased_direction", "unphased", "[unphased_direction]", function(ck) { i <- which(ck$gene$sig)[1]; ck$gene$direction[i] <- "haplotype A higher"; ck })
run("nonsig_direction", "unphased", "[nonsig_none]", function(ck) {
  i <- which(!ck$gene$sig)[1]; ck$gene$direction[i] <- "no direction (haplotype labels arbitrary)"; ck })
run("counts_altered", "phased", "[input_complete]", function(ck) { i <- which(plt(ck$gene))[1]; ck$gene$aCount[i] <- ck$gene$aCount[i] + 1; ck })
run("row_count_log", "phased", "[row_counts]", edit_dir = function(d) {
  h <- list.files(d, pattern = "_05_phaser\\.html$", full.names = TRUE); L <- readLines(h, warn = FALSE)
  k <- grep("with haplotype counts;", L)
  L[k] <- sub("([0-9]+) with haplotype counts;", "999999 with haplotype counts;", L[k]); writeLines(L, h) })
# dropped rows: two significant two_block cells and five null cells removed from the gene table. The completeness gate must fail,
# and the previous oracle (computed on the Rmd's own gene table) would have passed this copy: shown after the run.
old_oracle <- function(ck) {
  g <- ck$gene; MIN_D <- ck$constants$MIN_DEPTH; u <- g[g$totalCount >= MIN_D, ]
  orp <- function(x, n, rho = rho_true) { k <- 0:n; a <- 0.5 * (1 - rho) / rho; lp <- lchoose(n, k) + lbeta(k + a, n - k + a) - lbeta(a, a)
    min(1, sum(exp(lp[lp <= lp[x + 1] + 1e-9]))) }
  u$padj_or <- ave(mapply(orp, u$aCount, u$totalCount), u$sample, FUN = function(p) p.adjust(p, method = "BH"))
  u$sig_or <- u$padj_or < ck$constants$FDR_SIG & pmax(u$aCount, u$bCount) / u$totalCount - 0.5 >= ck$constants$ABS_DEV_SIG
  tbk <- cls[g$gene_id] == "two_block"; n_or <- sum(u$sig_or[cls[u$gene_id] == "two_block"]); obs <- sum(g$sig[tbk])
  pass <- obs >= 0.9 * n_or
  cat(sprintf("%s old oracle (on the Rmd gene table) for this copy: observed %d, oracle %d, bound >= %.1f: the old gate would %s\n",
              if (pass) "PROOF OK " else "PROOF BAD", obs, n_or, 0.9 * n_or, if (pass) "PASS (blind to dropped rows)" else "FAIL"))
  if (!pass) bad <<- 1
}
run("dropped_rows", "phased", "[input_complete]", function(ck) {
  g <- ck$gene; i <- c(which(g$sig & cls[g$gene_id] == "two_block")[1:2], which(cls[g$gene_id] == "null")[1:5])
  ck$gene <- g[-i, ]; ck$comparison <- ck$comparison[!paste(ck$comparison$sample, ck$comparison$gene_id) %in% paste(g$sample[i], g$gene_id[i]), ]; ck },
  edit_tsv = NULL, after = old_oracle)
cat(if (bad == 0) "ALL TAMPER PROOFS OK\n" else "SOME TAMPER PROOFS BAD\n")
quit(status = bad)
