# Usage: Rscript prove_evaluate_stage3_phaser.R <phased RESULTS_DIR> <unphased RESULTS_DIR> <truth dir (outbred_phase)> <work dir>
#                [<phased ASEReadCounter table dir> <unphased ASEReadCounter table dir>]
# Non-vacuity proofs for evaluate_stage3.R, phaser mode (gates [ph_*]). The two RESULTS_DIRs hold the raw outputs of phaser_count.sh
# (phaser/<sample>.haplotypic_counts.txt and .gene_ae.txt) on aligned synthetic data: phased = PHASED_GT 1 with the phased VCFs,
# unphased = PHASED_GT 0 with the unphased VCFs. The ASEReadCounter tables (<sample>.table; default <RESULTS_DIR>/ase_counts) give the
# oracle of the phased two_block gate; they are copied into every tampered copy and passed to the evaluator as its 5th argument.
# 1. Both untampered results must make the evaluator exit 0.
# 2. Each tampered copy (written under <work dir>/tamper_phaser/<name>; only the first sample's tables are edited) must make it exit 1
#    with the FAIL line of the named gate tag.
# Every [ph_*] gate tag of evaluate_stage3.R appears below at least once (check_skill.sh checks this).
# Ends with "ALL PHASER TAMPER PROOFS OK" and exit 0, or "SOME PHASER TAMPER PROOFS BAD" and exit 1. Run it in an sbatch job (R).
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(4, 6)) stop("usage: Rscript prove_evaluate_stage3_phaser.R <phased RESULTS_DIR> <unphased RESULTS_DIR> <truth dir> <work dir> [<phased ASE dir> <unphased ASE dir>]", call. = FALSE)
base <- c(phased = args[1], unphased = args[2]); TRUTH <- args[3]; WORK <- args[4]
ase <- if (length(args) == 6) c(phased = args[5], unphased = args[6]) else c(phased = file.path(args[1], "ase_counts"), unphased = file.path(args[2], "ase_counts"))
self <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
EV <- file.path(dirname(normalizePath(self)), "evaluate_stage3.R")
if (!file.exists(EV)) stop("evaluate_stage3.R not found next to this script: ", EV, call. = FALSE)
rd <- function(f, ...) read.delim(file.path(TRUTH, f), stringsAsFactors = FALSE, ...)
tg <- rd("truth_genes.tsv"); tig <- rd("truth_individual_genes.tsv")
tp <- rd("truth_phase.tsv", colClasses = c(gt = "character", alt_on = "character"))
sm <- read.csv(file.path(TRUTH, "samples.csv"), stringsAsFactors = FALSE)
S1 <- sm$sample[1]; I1 <- sm$individual[1]
ph1 <- tp[tp$individual == I1 & tp$gt == "0/1", ]
cls <- setNames(tg$class, tg$gene_id)
bad <- 0
ts <- rd("truth_snps.tsv")
evaluate <- function(dir, gt, ase_dir) {
  out <- suppressWarnings(system2("Rscript", c(EV, "phaser", dir, TRUTH, gt, ase_dir), stdout = TRUE, stderr = TRUE))
  st <- attr(out, "status"); list(out = out, status = if (is.null(st)) 0 else st)
}
for (gt in names(base)) {
  r <- evaluate(base[[gt]], gt, ase[[gt]])
  cat(sprintf("UNTAMPERED %-8s exit %d; %d PASS lines, %d FAIL lines\n", gt, r$status, sum(startsWith(r$out, "PASS")), sum(startsWith(r$out, "FAIL"))))
  if (r$status != 0 || !any(startsWith(r$out, "PASS [ph_"))) bad <- 1
}
if (bad) { cat("SOME PHASER TAMPER PROOFS BAD (an untampered result does not pass)\n"); quit(status = 1) }
rt <- function(f) read.delim(f, colClasses = "character", quote = "", na.strings = character(0), check.names = FALSE)
wt <- function(d, f) write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE)
pos_of <- function(v) as.integer(vapply(strsplit(strsplit(v, ",")[[1]], "_"), `[`, "", 2))
run <- function(name, gt, expect, hc = NULL, ga = NULL, ac = NULL) {   # hc / ga / ac: functions that edit the first sample's tables
  d <- file.path(WORK, "tamper_phaser", name); unlink(d, recursive = TRUE); dir.create(file.path(d, "phaser"), recursive = TRUE)
  dir.create(file.path(d, "ase_counts"))
  file.copy(list.files(file.path(base[[gt]], "phaser"), pattern = "\\.(haplotypic_counts|gene_ae)\\.txt$", full.names = TRUE), file.path(d, "phaser"))
  file.copy(file.path(ase[[gt]], paste0(sm$sample, ".table")), file.path(d, "ase_counts"))
  fh <- file.path(d, "phaser", paste0(S1, ".haplotypic_counts.txt")); fg <- file.path(d, "phaser", paste0(S1, ".gene_ae.txt"))
  fa <- file.path(d, "ase_counts", paste0(S1, ".table"))
  if (!is.null(hc)) wt(hc(rt(fh)), fh)
  if (!is.null(ga)) wt(ga(rt(fg)), fg)
  if (!is.null(ac)) wt(ac(rt(fa)), fa)
  r <- evaluate(d, gt, file.path(d, "ase_counts"))
  hit <- r$out[startsWith(r$out, paste("FAIL", expect))]
  ok <- r$status == 1 && length(hit) > 0
  cat(sprintf("%s tamper %-22s exit %d; expected 'FAIL %s': %s\n", if (ok) "PROOF OK " else "PROOF BAD", name, r$status, expect,
              if (length(hit)) hit[1] else paste("not found; FAIL lines:", paste(grep("^FAIL", r$out, value = TRUE), collapse = " | "))))
  if (!ok) bad <<- 1
}
multi <- function(h) which(as.integer(h$variantCount) >= 2)
is_trans <- function(v) length(unique(ph1$alt_on[match(pos_of(v), ph1$pos)])) > 1
# [ph_blocks]: one multi-SNP block phased against the truth (its first allele swapped between haplotype A and B)
run("block_phase_wrong", "unphased", "[ph_blocks]", hc = function(h) {
  j <- multi(h)[1]; A <- strsplit(h$haplotypeA[j], ",")[[1]]; B <- strsplit(h$haplotypeB[j], ",")[[1]]
  t <- A[1]; A[1] <- B[1]; B[1] <- t; h$haplotypeA[j] <- paste(A, collapse = ","); h$haplotypeB[j] <- paste(B, collapse = ","); h })
# [ph_blocks]: only 19 multi-SNP blocks left
run("too_few_blocks", "phased", "[ph_blocks]", hc = function(h) { m <- multi(h); h[-m[-(1:19)], ] })
# [ph_trans]: every block that joins SNPs in trans removed
run("no_trans_block", "unphased", "[ph_trans]", hc = function(h) h[!vapply(seq_len(nrow(h)), function(j) as.integer(h$variantCount[j]) >= 2 && is_trans(h$variants[j]), logical(1)), ])
# [ph_counts]: one gene with more reads than the simulated fragments covering its het SNPs
run("count_inflated", "phased", "[ph_counts]", ga = function(g) {
  j <- which(as.numeric(g$totalCount) > 0)[1]
  g$aCount[j] <- as.character(as.numeric(g$aCount[j]) + 100000); g$totalCount[j] <- as.character(as.numeric(g$totalCount[j]) + 100000); g })
# [ph_two_block]: unphased run claiming both blocks; phased run with one block only
run("two_block_unphased", "unphased", "[ph_two_block]", ga = function(g) { g$n_variants[cls[g$name] %in% "two_block"] <- "6"; g })
run("two_block_phased", "phased", "[ph_two_block]", ga = function(g) { g$n_variants[cls[g$name] %in% "two_block"] <- "3"; g })
# phased two_block gate (controller ruling): oracle = het SNPs with >= 1 read in the sample's ASEReadCounter table
vid <- function(p) { k <- match(p, ts$pos); paste(ts$chrom[k], ts$pos[k], ts$ref[k], ts$alt[k], sep = "_") }
tb_edit <- function(g, f) {   # f(positions used) -> new positions; rewrites variants and n_variants of the two_block genes
  for (j in which(cls[g$name] %in% "two_block")) { p <- sort(f(pos_of(g$variants[j]), g$name[j]))
    g$variants[j] <- paste(vid(p), collapse = ","); g$n_variants[j] <- as.character(length(p)) }
  g }
# [ph_two_block_over]: 6 variants reported where ASEReadCounter covers only 5 (the uncovered het SNP added)
run("two_block_overreport", "phased", "[ph_two_block_over]", ga = function(g) tb_edit(g, function(p, x) union(p, ph1$pos[ph1$gene_id == x])))
# [ph_two_block_under]: a covered SNP dropped from each two_block gene (still both clusters)
run("two_block_dropped_snp", "phased", "[ph_two_block_under]", ga = function(g) tb_edit(g, function(p, x) p[-length(p)]))
# [ph_two_block_over] again: the oracle is read from the ASEReadCounter table passed in: one covered two_block SNP removed from it
run("ase_table_dropped_snp", "phased", "[ph_two_block_over]", ac = function(a) {
  tbp <- ph1$pos[ph1$gene_id %in% names(cls)[cls == "two_block"]]; a[-which(as.integer(a$position) %in% tbp)[1], ] })
# [ph_two_block]: variants of one cluster only (count unchanged)
run("two_block_one_cluster", "phased", "[ph_two_block]", ga = function(g) {
  for (j in which(cls[g$name] %in% "two_block")) { p <- pos_of(g$variants[j]); c2 <- ts$cluster[match(p, ts$pos)] == 2
    g$variants[j] <- paste(vid(c(p[c2], p[c2][seq_len(sum(!c2))])), collapse = ",") }
  g })
# [ph_gw_phased]: one covered gene not genome-wide phased in the phased run
run("not_gw_phased", "phased", "[ph_gw_phased]", ga = function(g) {
  j <- which(as.numeric(g$totalCount) > 0 & !cls[g$name] %in% "two_block")[1]; g$gw_phased[j] <- "0"; g })
# [ph_orient]: aCount and bCount swapped in one planted gene whose fraction deviates by more than 0.1
h1 <- setNames(tig$h1[tig$individual == I1], tig$gene_id[tig$individual == I1])
run("orient_swapped", "phased", "[ph_orient]", ga = function(g) {
  a <- as.numeric(g$aCount); n <- as.numeric(g$totalCount)
  j <- which(n > 0 & !is.na(h1[g$name]) & h1[g$name] != 0.5 & abs(a / n - 0.5) > 0.1)[1]
  t <- g$aCount[j]; g$aCount[j] <- g$bCount[j]; g$bCount[j] <- t; g })
# [ph_orient_both]: the planted genes whose haplotype 1 carries ALT at the first het SNP lose their counts, so the orientation check
# would only see one kind of gene (a REF- or ALT-based rule could then pass it)
first_alt <- vapply(split(ph1, ph1$gene_id), function(p) p$alt_on[which.min(p$pos)] == "1", logical(1))
run("orient_one_kind", "phased", "[ph_orient_both]", ga = function(g) {
  j <- !is.na(h1[g$name]) & h1[g$name] != 0.5 & g$name %in% names(first_alt)[first_alt]
  g$aCount[j] <- "0"; g$bCount[j] <- "0"; g$totalCount[j] <- "0"; g })
cat(if (bad == 0) "ALL PHASER TAMPER PROOFS OK\n" else "SOME PHASER TAMPER PROOFS BAD\n")
quit(status = bad)
