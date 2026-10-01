# Usage: Rscript evaluate_stage3.R rmd    <RESULTS_DIR> <truth dir (outbred_phase)> <phased|unphased>
#        Rscript evaluate_stage3.R phaser <RESULTS_DIR> <truth dir (outbred_phase)> <phased|unphased>
# rmd: Rmd 05 results (ase_phaser_checkpoint.rds, summary_numbers_phaser.tsv) against the planted truth.
# phaser: the raw phASER outputs (RESULTS_DIR/phaser/<sample>.*) against the truth (added by the phaser_count task).
# Prints PASS / FAIL / REPORT lines; exits 1 if any criterion fails. REPORT lines are informational.
# Orientation (truth): haplotype 1 is the left allele of the phased genotype (alt_on 1 is written 1|0), and phASER's
# aCount is haplotype 1 in genome-wide phased genes, so h1 > 0.5 means "haplotype A higher".
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4 || !args[1] %in% c("rmd", "phaser") || !args[4] %in% c("phased", "unphased"))
  stop("usage: Rscript evaluate_stage3.R <rmd|phaser> <RESULTS_DIR> <truth dir> <phased|unphased>", call. = FALSE)
what <- args[1]; RES <- args[2]; TRUTH <- args[3]; GT <- args[4]
fail <- 0
crit <- function(cond, msg) { cat(if (isTRUE(cond)) "PASS" else "FAIL", msg, "\n"); if (!isTRUE(cond)) fail <<- 1 }
rd <- function(f, ...) read.delim(file.path(TRUTH, f), stringsAsFactors = FALSE, ...)
tg <- rd("truth_genes.tsv"); tig <- rd("truth_individual_genes.tsv"); tf <- rd("truth_fragments.tsv")
tp <- rd("truth_phase.tsv", colClasses = c(gt = "character", alt_on = "character")); ts <- rd("truth_snps.tsv")
sm <- read.csv(file.path(TRUTH, "samples.csv"), stringsAsFactors = FALSE)
cell <- merge(merge(tf, sm[, c("sample", "individual")], by = "sample"), tig, by = c("individual", "gene_id"))
cell <- merge(cell, tg[, c("gene_id", "class")], by = "gene_id")
if (what == "phaser") stop("the phaser mode is not written yet (phaser_count task); nothing was checked", call. = FALSE)
if (what == "rmd") {
  ck <- readRDS(file.path(RES, "ase_phaser_checkpoint.rds"))
  crit(identical(as.numeric(ck$constants$PHASED_GT), if (GT == "phased") 1 else 0), sprintf("checkpoint PHASED_GT matches '%s'", GT))
  g <- ck$gene; MIN_D <- ck$constants$MIN_DEPTH
  crit(nrow(g) > 0 && all(g$totalCount > 0) && all(is.na(g$p[g$totalCount < MIN_D])),
       "no gene row with totalCount 0, and every row below MIN_DEPTH has p NA")
  m <- merge(cell, g[, c("sample", "gene_id", "aCount", "bCount", "totalCount", "gw_phased", "p", "padj", "sig", "hap_A_frac",
                         "major_frac", "direction")], by = c("sample", "gene_id"), all.x = TRUE)
  m$sig[is.na(m$sig)] <- FALSE
  cnt <- function(cls) c(sum(m$sig[m$class == cls]), sum(m$class == cls))
  x <- cnt("hap_strong"); crit(x[2] > 0 && x[1] == x[2], sprintf("hap_strong cells significant: %d of %d (all)", x[1], x[2]))
  # two_block (controller ruling, Task 4): per-cell power at about 100 reads and h1 = 0.8 is well below 1, so "all 12" is not a
  # justified gate. Bound: observed calls >= 0.9 x the calls of an ORACLE test on the SAME counts: exact two-sided beta-binomial
  # against 0.5 with the TRUE dispersion PHI_BIO = 0.005 (simulate_ase_stage3.R, line 18: each cell's haplotype-1 probability is
  # drawn from a beta with mean h1 and rho PHI_BIO, then the read pairs binomially), BH within sample over the same tested rows
  # (totalCount >= MIN_DEPTH), and the same sig rule (padj < FDR_SIG and major fraction - 0.5 >= ABS_DEV_SIG).
  PHI_TRUE <- 0.005
  or_p <- function(x, n, rho) {   # same definition as bb_pvalue of the skill, written out here so the oracle does not depend on it
    k <- 0:n; a <- 0.5 * (1 - rho) / rho
    lp <- lchoose(n, k) + lbeta(k + a, n - k + a) - lbeta(a, a)
    min(1, sum(exp(lp[lp <= lp[x + 1] + 1e-9])))
  }
  gt_rows <- g[g$totalCount >= MIN_D, c("sample", "gene_id", "aCount", "totalCount")]
  gt_rows$p_or <- mapply(or_p, gt_rows$aCount, gt_rows$totalCount, MoreArgs = list(rho = PHI_TRUE))
  gt_rows$padj_or <- ave(gt_rows$p_or, gt_rows$sample, FUN = function(p) p.adjust(p, method = "BH"))
  gt_rows$sig_or <- gt_rows$padj_or < ck$constants$FDR_SIG &
    pmax(gt_rows$aCount, gt_rows$totalCount - gt_rows$aCount) / gt_rows$totalCount - 0.5 >= ck$constants$ABS_DEV_SIG
  m <- merge(m, gt_rows[, c("sample", "gene_id", "padj_or", "sig_or")], by = c("sample", "gene_id"), all.x = TRUE)
  m$sig_or[is.na(m$sig_or)] <- FALSE
  x <- cnt("two_block"); n_or <- sum(m$sig_or[m$class == "two_block"])
  crit(x[2] > 0 && n_or > 0 && x[1] >= 0.9 * n_or,
       sprintf("two_block cells significant: observed %d, oracle (true rho %.3f, same counts) %d, bound >= %.1f, of %d cells", x[1], PHI_TRUE, n_or, 0.9 * n_or, x[2]))
  for (cls in c("hap_strong", "hap_moderate", "hap_lowdepth", "null", "null_linked"))
    cat(sprintf("REPORT oracle (true rho %.3f) %-12s: observed %d, oracle %d of %d cells\n", PHI_TRUE, cls, cnt(cls)[1], sum(m$sig_or[m$class == cls]), cnt(cls)[2]))
  x <- cnt("hap_lowdepth")
  if (GT == "phased") {
    crit(x[2] == 18 && x[1] >= 15, sprintf("hap_lowdepth cells significant (phased genotypes): %d of %d (>= 15)", x[1], x[2]))
  } else {
    cat(sprintf("REPORT hap_lowdepth cells significant (unphased genotypes, single best block): %d of %d\n", x[1], x[2]))
  }
  cmp <- ck$comparison
  lo <- merge(m[m$class == "hap_lowdepth", c("sample", "gene_id", "sig")], cmp[, c("sample", "gene_id", "category")],
              by = c("sample", "gene_id"), all.x = TRUE)
  n_only <- sum(lo$sig & lo$category %in% "phASER only (no SNP tested unphased)")
  crit(n_only == sum(lo$sig), sprintf("significant hap_lowdepth cells labelled 'phASER only (no SNP tested unphased)': %d of %d", n_only, sum(lo$sig)))
  mod <- tg$gene_id[tg$class == "hap_moderate"]
  x <- cnt("hap_moderate")
  cat(sprintf("REPORT hap_moderate cells significant: phASER %d of %d; unphased ACAT (Rmd 02) %d\n", x[1], x[2],
              sum(cmp$sig_unphased[cmp$gene_id %in% mod])))
  nul <- m$class %in% c("null", "null_linked")
  crit(sum(nul) > 0 && sum(m$sig[nul]) <= 6, sprintf("null cells significant: %d of %d tested (limit 6)", sum(m$sig[nul]), sum(nul & !is.na(m$p))))
  pl <- m$sig & !nul
  if (GT == "phased") {
    agree <- sign(m$hap_A_frac[pl] - 0.5) == sign(m$h1[pl] - 0.5)
    crit(sum(pl) > 0 && all(m$gw_phased[pl]) && all(agree),
         sprintf("phased: significant planted cells genome-wide phased, haplotype A = haplotype 1 (left GT allele): %d of %d", sum(agree), sum(pl)))
    lab <- ifelse(m$h1[pl] > 0.5, "haplotype A higher", "haplotype B higher")
    crit(sum(pl) > 0 && all(m$direction[pl] == lab),
         sprintf("phased: direction label of significant planted cells follows haplotype 1: %d of %d", sum(m$direction[pl] == lab), sum(pl)))
  } else {
    crit(sum(g$sig) > 0 && all(g$direction[g$sig] == "no direction (haplotype labels arbitrary)"), "unphased: no significant row carries a direction")
  }
  crit(all(g$direction[!g$sig] == "none"), "rows that are not significant have the direction 'none'")
  s <- read.delim(file.path(RES, "summary_numbers_phaser.tsv"), stringsAsFactors = FALSE)
  crit(identical(names(s), c("sample", "genotypes", "genes_with_counts", "genes_tested", "sig_genes", "gw_phased_genes", "sig_unphased", "sig_both",
                             "sig_phaser_only", "sig_phaser_only_untested_unphased", "sig_unphased_only", "rho_used", "rho_own", "rho_cohort")) &&
       nrow(s) == nrow(sm) && all(s$sig_genes == vapply(s$sample, function(x) sum(g$sig[g$sample == x]), integer(1))),
       "summary_numbers_phaser.tsv: exact columns, one row per sample, sig_genes equal to the gene table")
  # informational: the phASER test next to the unphased ACAT test of Rmd 02, per planted class (sample-gene cells)
  cc <- merge(cell[, c("sample", "gene_id", "class")], cmp[, c("sample", "gene_id", "sig", "sig_unphased", "category")],
              by = c("sample", "gene_id"), all.x = TRUE)
  cc$sig[is.na(cc$sig)] <- FALSE; cc$sig_unphased[is.na(cc$sig_unphased)] <- FALSE
  for (cls in c("hap_strong", "hap_moderate", "hap_lowdepth", "two_block", "null", "null_linked")) {
    i <- cc$class == cls
    cat(sprintf("REPORT %-12s cells %3d: phASER sig %3d, unphased ACAT sig %3d, both %3d, phASER only %3d (no SNP tested unphased %3d), unphased only %3d\n",
                cls, sum(i), sum(cc$sig[i]), sum(cc$sig_unphased[i]), sum(cc$sig[i] & cc$sig_unphased[i]), sum(cc$sig[i] & !cc$sig_unphased[i]),
                sum(i & cc$category %in% "phASER only (no SNP tested unphased)"), sum(!cc$sig[i] & cc$sig_unphased[i])))
  }
  cat(sprintf("REPORT null false-positive rate: phASER %d of %d null cells tested (%.3f); unphased ACAT %d of %d null cells\n",
              sum(m$sig[nul]), sum(nul & !is.na(m$p)), sum(m$sig[nul]) / max(1, sum(nul & !is.na(m$p))),
              sum(cc$sig_unphased[cc$class %in% c("null", "null_linked")]), sum(cc$class %in% c("null", "null_linked"))))
  print(m[!nul, c("sample", "gene_id", "class", "h1", "totalCount", "gw_phased", "hap_A_frac", "padj", "sig", "direction")], row.names = FALSE)
}
quit(status = fail)
