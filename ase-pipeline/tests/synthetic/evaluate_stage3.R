# Usage: Rscript evaluate_stage3.R rmd    <RESULTS_DIR> <truth dir (outbred_phase)> <phased|unphased>
#        Rscript evaluate_stage3.R phaser <RESULTS_DIR> <truth dir (outbred_phase)> <phased|unphased> [ASEReadCounter table dir]
# rmd: Rmd 05 results (ase_phaser_checkpoint.rds, summary_numbers_phaser.tsv, the Rmd 05 HTML) against the planted truth and against
#      the input phASER gene tables the Rmd read (RESULTS_DIR/phaser/<sample>.gene_ae.txt).
# phaser: the raw phASER outputs of phaser_count.sh (RESULTS_DIR/phaser/<sample>.haplotypic_counts.txt and .gene_ae.txt) against the
#      truth: block phasing, trans blocks, count bounds, two_block genes, genome-wide phase and orientation (gates [ph_*]).
# Prints PASS / FAIL / REPORT lines; exits 1 if any criterion fails. REPORT lines are informational. Every gate message starts with
# a [tag]; prove_evaluate_stage3.R (rmd mode) and prove_evaluate_stage3_phaser.R (phaser mode, tags [ph_*]) hold a tampered-copy
# proof for every tag (check_skill.sh checks that).
# Orientation (truth): haplotype 1 is the left allele of the phased genotype (alt_on 1 is written 1|0), and phASER's aCount is
# haplotype 1 in genome-wide phased genes, so h1 > 0.5 means "haplotype A higher".
# Not exercised by the synthetic data: in the phased tables every gene with counts is genome-wide phased, so the rule "rows that are
# not genome-wide phased get no direction under PHASED_GT = 1" is not tested here.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 4:5 || !args[1] %in% c("rmd", "phaser") || !args[4] %in% c("phased", "unphased"))
  stop("usage: Rscript evaluate_stage3.R <rmd|phaser> <RESULTS_DIR> <truth dir> <phased|unphased> [ASEReadCounter table dir]", call. = FALSE)
what <- args[1]; RES <- args[2]; TRUTH <- args[3]; GT <- args[4]
# phaser mode, phased two_block gate: the ASEReadCounter tables <sample>.table (default RESULTS_DIR/ase_counts) give the oracle
ASE_DIR <- if (length(args) == 5) args[5] else file.path(RES, "ase_counts")
fail <- 0
crit <- function(cond, msg) { cat(if (isTRUE(cond)) "PASS" else "FAIL", msg, "\n"); if (!isTRUE(cond)) fail <<- 1 }
rd <- function(f, ...) read.delim(file.path(TRUTH, f), stringsAsFactors = FALSE, ...)
tg <- rd("truth_genes.tsv"); tig <- rd("truth_individual_genes.tsv"); tf <- rd("truth_fragments.tsv")
# truth_phase.tsv and truth_snps.tsv are read for the phaser mode; the rmd mode does not use them
tp <- rd("truth_phase.tsv", colClasses = c(gt = "character", alt_on = "character")); ts <- rd("truth_snps.tsv")
sm <- read.csv(file.path(TRUTH, "samples.csv"), stringsAsFactors = FALSE)
cell <- merge(merge(tf, sm[, c("sample", "individual")], by = "sample"), tig, by = c("individual", "gene_id"))
cell <- merge(cell, tg[, c("gene_id", "class")], by = "gene_id")
# true dispersion of the generator: PHI_BIO in simulate_ase_stage3.R next to this script (each cell's haplotype-1 probability is
# drawn from a beta with mean h1 and rho PHI_BIO, then the read pairs binomially)
self <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
gen <- file.path(dirname(normalizePath(self)), "simulate_ase_stage3.R")
if (!file.exists(gen)) stop("the generator ", gen, " is missing: the true dispersion (PHI_BIO) cannot be read", call. = FALSE)
pl_line <- grep("PHI_BIO <- [0-9.eE-]+", readLines(gen), value = TRUE)
if (length(pl_line) != 1) stop("expected exactly one 'PHI_BIO <- <number>' assignment in ", gen, call. = FALSE)
PHI_TRUE <- suppressWarnings(as.numeric(sub(".*PHI_BIO <- ([0-9.eE-]+).*", "\\1", pl_line)))
if (!is.finite(PHI_TRUE) || PHI_TRUE <= 0) stop("PHI_BIO in ", gen, " is not a positive number", call. = FALSE)
if (what == "phaser") {
  # raw phASER outputs of phaser_count.sh (RESULTS_DIR/phaser/<sample>.haplotypic_counts.txt and .gene_ae.txt) against the truth
  fin <- file.path(RES, "phaser", paste0(rep(sm$sample, each = 2), c(".haplotypic_counts.txt", ".gene_ae.txt")))
  if (!all(file.exists(fin))) stop("phASER outputs missing: ", paste(fin[!file.exists(fin)], collapse = ", "), call. = FALSE)
  pos_of <- function(v) as.integer(vapply(strsplit(strsplit(v, ",")[[1]], "_"), `[`, "", 2))
  for (i in seq_len(nrow(sm))) {
    s <- sm$sample[i]; ind <- sm$individual[i]
    hc <- read.delim(file.path(RES, "phaser", paste0(s, ".haplotypic_counts.txt")), stringsAsFactors = FALSE, quote = "",
                     colClasses = c(variants = "character", haplotypeA = "character", haplotypeB = "character", blockGWPhase = "character"))
    ga <- read.delim(file.path(RES, "phaser", paste0(s, ".gene_ae.txt")), stringsAsFactors = FALSE, quote = "",
                     colClasses = c(name = "character", gw_phased = "character", variants = "character", log2_aFC = "character"))
    ph <- tp[tp$individual == ind & tp$gt == "0/1", ]
    h1allele <- setNames(ifelse(ph$alt_on == "1", ts$alt[match(ph$pos, ts$pos)], ts$ref[match(ph$pos, ts$pos)]), ph$pos)
    multi <- hc[hc$variantCount >= 2, ]
    okb <- vapply(seq_len(nrow(multi)), function(j) {
      A <- strsplit(multi$haplotypeA[j], ",")[[1]]; t1 <- unname(h1allele[as.character(pos_of(multi$variants[j]))])
      !anyNA(t1) && (all(A == t1) || all(A != t1)) }, logical(1))
    ntr <- sum(vapply(seq_len(nrow(multi)), function(j) length(unique(ph$alt_on[match(pos_of(multi$variants[j]), ph$pos)])) > 1, logical(1)))
    crit(nrow(multi) >= 20 && mean(okb) >= 0.99,
         sprintf("[ph_blocks] %s: %d multi-SNP blocks (>= 20), %.4f phased as in the truth (>= 0.99)", s, nrow(multi), mean(okb)))
    crit(ntr >= 1, sprintf("[ph_trans] %s: blocks joining SNPs in trans: %d (>= 1)", s, ntr))
    g <- merge(ga[ga$totalCount > 0, ], tf[tf$sample == s, ], by.x = "name", by.y = "gene_id")
    crit(nrow(g) > 0 && all(g$totalCount <= g$cover_h1 + g$cover_h2),
         sprintf("[ph_counts] %s: no gene count exceeds the simulated fragments covering its het SNPs (%d genes; median ratio %.3f)", s, nrow(g),
                 median(g$totalCount / (g$cover_h1 + g$cover_h2))))
    tb <- g[g$name %in% tg$gene_id[tg$class == "two_block"], ]
    if (GT == "unphased") {
      crit(nrow(tb) == 2 && all(tb$n_variants <= 3), sprintf("[ph_two_block] %s: two_block genes use one block (n_variants %s <= 3)", s, paste(tb$n_variants, collapse = ",")))
    } else {
      # controller ruling (Task 5): the oracle is the number of the gene's het SNPs with >= 1 read in THIS sample's ASEReadCounter
      # table (an independent tool, not phASER's allelic_counts); the first SNP of each two_block gene sits at transcript offset 2
      # and the WASP filter removes the few reads there, so 6 is not observable on the synthetic data
      act <- read.delim(file.path(ASE_DIR, paste0(s, ".table")), stringsAsFactors = FALSE)
      cov_key <- paste(act$contig, act$position)[act$totalCount >= 1]   # contig AND position
      orc <- vapply(tb$name, function(x) sum(paste(ph$chrom, ph$pos)[ph$gene_id == x] %in% cov_key), integer(1))
      both <- vapply(tb$variants, function(v) length(unique(ts$cluster[match(pos_of(v), ts$pos)])) == 2, logical(1))
      crit(nrow(tb) == 2 && all(tb$n_variants == orc & both & tb$gw_phased == "1"),
           sprintf("[ph_two_block] %s: two_block genes use both blocks: n_variants %s = het SNPs with reads in ASEReadCounter %s, both clusters %s, genome-wide phased %s",
                   s, paste(tb$n_variants, collapse = ","), paste(orc, collapse = ","), paste(both, collapse = ","), paste(tb$gw_phased, collapse = ",")))
      crit(nrow(tb) == 2 && all(tb$n_variants <= orc), sprintf("[ph_two_block_over] %s: two_block genes report no more variants than ASEReadCounter covers (%s <= %s)",
                                                                s, paste(tb$n_variants, collapse = ","), paste(orc, collapse = ",")))
      crit(nrow(tb) == 2 && all(tb$n_variants >= orc), sprintf("[ph_two_block_under] %s: two_block genes use every het SNP that ASEReadCounter covers (%s >= %s)",
                                                                s, paste(tb$n_variants, collapse = ","), paste(orc, collapse = ",")))
      crit(all(g$gw_phased == "1"), sprintf("[ph_gw_phased] %s: every covered gene genome-wide phased (%d of %d)", s, sum(g$gw_phased == "1"), nrow(g)))
      h <- merge(g, tig[tig$individual == ind, c("gene_id", "h1")], by.x = "name", by.y = "gene_id")
      pl <- h[h$h1 != 0.5 & abs(h$aCount / h$totalCount - 0.5) > 0.1, ]
      ag <- sign(pl$aCount / pl$totalCount - 0.5) == sign(pl$h1 - 0.5)
      crit(nrow(pl) > 0 && all(ag), sprintf("[ph_orient] %s: aCount is haplotype 1 (left GT allele) in %d of %d planted genes", s, sum(ag), nrow(pl)))
      # carried item (Task 2): cis/trans are unbalanced per class and individual, so the genes checked above must hold both phases:
      # genes whose haplotype 1 carries the ALT allele at the first het SNP and genes where it carries REF (a rule "aCount = the
      # ALT (or REF) haplotype of the first SNP" agrees with "aCount = haplotype 1" on one kind only), and at least one trans gene
      first_alt <- vapply(pl$name, function(x) { p <- ph[ph$gene_id == x, ]; p$alt_on[which.min(p$pos)] == "1" }, logical(1))
      trans_g <- vapply(pl$name, function(x) length(unique(ph$alt_on[ph$gene_id == x])) > 1, logical(1))
      crit(sum(first_alt) >= 1 && sum(!first_alt) >= 1 && sum(trans_g) >= 1,
           sprintf("[ph_orient_both] %s: planted genes checked for orientation: %d with ALT on haplotype 1 at the first het SNP, %d with REF; %d trans (each >= 1)",
                   s, sum(first_alt), sum(!first_alt), sum(trans_g)))
    }
    lo <- g[g$name %in% tg$gene_id[tg$class == "hap_lowdepth"], ]
    cat(sprintf("REPORT %s: hap_lowdepth genes, SNPs used: %s of 8\n", s, paste(lo$n_variants, collapse = ",")))
    cat(sprintf("REPORT %s: %d blocks (%d multi-SNP, %d with a trans pair), %d of %d multi-SNP blocks phased as in the truth; %d genes with counts (%d genome-wide phased)\n",
                s, nrow(hc), nrow(multi), ntr, sum(okb), nrow(multi), nrow(g), sum(g$gw_phased == "1")))
    if (GT == "phased") {
      pc <- merge(cell[cell$sample == s, c("gene_id", "class", "h1", "frag_h1", "frag_h2")], g[, c("name", "aCount", "bCount", "totalCount")],
                  by.x = "gene_id", by.y = "name")
      pc <- pc[pc$h1 != 0.5, ]
      cat(sprintf("REPORT %s: planted genes recovered: %s\n", s, paste(sprintf("%s(%s h1 %.2f: truth %.3f, phASER aCount/total %.3f)", pc$gene_id, pc$class,
                  pc$h1, pc$frag_h1 / (pc$frag_h1 + pc$frag_h2), pc$aCount / pc$totalCount), collapse = "; ")))
    }
  }
}
if (what == "rmd") {
  ck <- readRDS(file.path(RES, "ase_phaser_checkpoint.rds"))
  crit(identical(as.numeric(ck$constants$PHASED_GT), if (GT == "phased") 1 else 0), sprintf("[phased_gt] checkpoint PHASED_GT matches '%s'", GT))
  g <- ck$gene; MIN_D <- ck$constants$MIN_DEPTH
  crit(nrow(g) > 0 && all(g$totalCount > 0) && all(is.na(g$p[g$totalCount < MIN_D])),
       "[rows] no gene row with totalCount 0, and every row below MIN_DEPTH has p NA")
  # the input phASER gene tables (what the Rmd read): independent of the Rmd's own output
  cls <- c(contig = "character", start = "numeric", stop = "numeric", name = "character", aCount = "numeric", bCount = "numeric",
           totalCount = "numeric", log2_aFC = "character", n_variants = "numeric", variants = "character", gw_phased = "character",
           bam = "character")
  fin <- file.path(RES, "phaser", paste0(sm$sample, ".gene_ae.txt"))
  if (!all(file.exists(fin)))
    stop("input phASER tables missing in ", file.path(RES, "phaser"), ": ", paste(sm$sample[!file.exists(fin)], collapse = ", "), call. = FALSE)
  inp <- do.call(rbind, lapply(seq_along(fin), function(i) {
    d <- read.delim(fin[i], colClasses = cls, quote = "", na.strings = character(0)); d$sample <- rep(sm$sample[i], nrow(d)); d }))
  kept <- inp[inp$totalCount > 0, c("sample", "name", "aCount", "bCount", "totalCount")]
  key_in <- paste(kept$sample, kept$name); key_g <- paste(g$sample, g$gene_id)
  mi <- match(key_in, key_g)
  same <- !is.na(mi)
  same[same] <- g$aCount[mi[same]] == kept$aCount[same] & g$bCount[mi[same]] == kept$bCount[same] & g$totalCount[mi[same]] == kept$totalCount[same]
  crit(nrow(kept) > 0 && !anyDuplicated(key_g) && all(same) && all(key_g %in% key_in),
       sprintf("[input_complete] every input row with counts is in the gene table once with identical counts, and nothing else: %d of %d input rows matched; %d gene rows, %d duplicated, %d not in the input",
               sum(same), nrow(kept), nrow(g), sum(duplicated(key_g)), sum(!key_g %in% key_in)))
  htm <- list.files(RES, pattern = "_05_phaser\\.html$", full.names = TRUE)
  ln <- if (length(htm) == 1) grep("[0-9]+ sample-gene rows read; [0-9]+ with haplotype counts;", readLines(htm, warn = FALSE), value = TRUE) else character(0)   # the printed line, not the echoed code
  hit <- if (length(ln) == 1) regmatches(ln, regexpr("[0-9]+ sample-gene rows read; [0-9]+ with haplotype counts; [0-9]+ at or above MIN_DEPTH", ln)) else character(0)
  num <- if (length(hit) == 1) as.numeric(strsplit(gsub("[^0-9]+", " ", hit), " ")[[1]]) else numeric(0)
  num <- num[!is.na(num)]
  want <- c(nrow(inp), nrow(kept), sum(kept$totalCount >= MIN_D))
  crit(length(htm) == 1 && length(num) == 3 && all(num == want),
       sprintf("[row_counts] the Rmd 05 log reports rows read / with counts / at MIN_DEPTH equal to the input tables: reported %s, input %s",
               paste(num, collapse = "/"), paste(want, collapse = "/")))
  m <- merge(cell, g[, c("sample", "gene_id", "aCount", "bCount", "totalCount", "gw_phased", "p", "padj", "sig", "hap_A_frac",
                         "major_frac", "direction")], by = c("sample", "gene_id"), all.x = TRUE)
  m$sig[is.na(m$sig)] <- FALSE
  cnt <- function(cls) c(sum(m$sig[m$class == cls]), sum(m$class == cls))
  x <- cnt("hap_strong"); crit(x[2] > 0 && x[1] == x[2], sprintf("[hap_strong] hap_strong cells significant: %d of %d (all)", x[1], x[2]))
  # two_block (controller ruling): observed calls >= 0.9 x the calls of an ORACLE test on the INPUT counts: exact two-sided
  # beta-binomial against 0.5 with the TRUE dispersion PHI_BIO, BH within sample over the input rows with totalCount >= MIN_DEPTH,
  # and the same sig rule (padj < FDR_SIG and major fraction - 0.5 >= ABS_DEV_SIG). The input tables serve both runs: for the
  # unphased run they hold the single-block counts the Rmd sees, which the truth cover counts (all blocks) are not.
  or_p <- function(x, n, rho) {   # same definition as bb_pvalue of the skill, written out here so the oracle does not depend on it
    k <- 0:n; a <- 0.5 * (1 - rho) / rho
    lp <- lchoose(n, k) + lbeta(k + a, n - k + a) - lbeta(a, a)
    min(1, sum(exp(lp[lp <= lp[x + 1] + 1e-9])))
  }
  orc <- kept[kept$totalCount >= MIN_D, ]
  orc$p_or <- mapply(or_p, orc$aCount, orc$totalCount, MoreArgs = list(rho = PHI_TRUE))
  orc$padj_or <- ave(orc$p_or, orc$sample, FUN = function(p) p.adjust(p, method = "BH"))
  orc$sig_or <- orc$padj_or < ck$constants$FDR_SIG &
    pmax(orc$aCount, orc$totalCount - orc$aCount) / orc$totalCount - 0.5 >= ck$constants$ABS_DEV_SIG
  m <- merge(m, data.frame(sample = orc$sample, gene_id = orc$name, sig_or = orc$sig_or), by = c("sample", "gene_id"), all.x = TRUE)
  m$sig_or[is.na(m$sig_or)] <- FALSE
  x <- cnt("two_block"); n_or <- sum(m$sig_or[m$class == "two_block"])
  crit(x[2] > 0 && n_or > 0 && x[1] >= 0.9 * n_or,
       sprintf("[two_block] two_block cells significant: observed %d, oracle (true rho %.3f, input counts) %d, bound >= %.1f, of %d cells", x[1], PHI_TRUE, n_or, 0.9 * n_or, x[2]))
  for (cls in c("hap_strong", "hap_moderate", "hap_lowdepth", "null", "null_linked"))
    cat(sprintf("REPORT oracle (true rho %.3f) %-12s: observed %d, oracle %d of %d cells\n", PHI_TRUE, cls, cnt(cls)[1], sum(m$sig_or[m$class == cls]), cnt(cls)[2]))
  x <- cnt("hap_lowdepth")
  if (GT == "phased") {
    crit(x[2] == 18 && x[1] >= 15, sprintf("[hap_lowdepth] hap_lowdepth cells significant (phased genotypes): %d of %d (>= 15)", x[1], x[2]))
  } else {
    cat(sprintf("REPORT hap_lowdepth cells significant (unphased genotypes, single best block): %d of %d\n", x[1], x[2]))
  }
  cmp <- ck$comparison
  lo <- merge(m[m$class == "hap_lowdepth", c("sample", "gene_id", "sig")], cmp[, c("sample", "gene_id", "category")],
              by = c("sample", "gene_id"), all.x = TRUE)
  n_only <- sum(lo$sig & lo$category %in% "phASER only (no SNP tested unphased)")
  crit(n_only == sum(lo$sig), sprintf("[lowdepth_category] significant hap_lowdepth cells labelled 'phASER only (no SNP tested unphased)': %d of %d", n_only, sum(lo$sig)))
  mod <- tg$gene_id[tg$class == "hap_moderate"]
  x <- cnt("hap_moderate")
  cat(sprintf("REPORT hap_moderate cells significant: phASER %d of %d; unphased ACAT (Rmd 02) %d\n", x[1], x[2],
              sum(cmp$sig_unphased[cmp$gene_id %in% mod])))
  nul <- m$class %in% c("null", "null_linked")
  crit(sum(nul) > 0 && sum(m$sig[nul]) <= 6, sprintf("[null] null cells significant: %d of %d tested (limit 6)", sum(m$sig[nul]), sum(nul & !is.na(m$p))))
  pl <- m$sig & !nul
  if (GT == "phased") {
    agree <- sign(m$hap_A_frac[pl] - 0.5) == sign(m$h1[pl] - 0.5)
    crit(sum(pl) > 0 && all(m$gw_phased[pl]) && all(agree),
         sprintf("[phased_orient] phased: significant planted cells genome-wide phased, haplotype A = haplotype 1 (left GT allele): %d of %d", sum(agree), sum(pl)))
    lab <- ifelse(m$h1[pl] > 0.5, "haplotype A higher", "haplotype B higher")
    crit(sum(pl) > 0 && all(m$direction[pl] == lab),
         sprintf("[phased_label] phased: direction label of significant planted cells follows haplotype 1: %d of %d", sum(m$direction[pl] == lab), sum(pl)))
  } else {
    crit(sum(g$sig) > 0 && all(g$direction[g$sig] == "no direction (haplotype labels arbitrary)"), "[unphased_direction] unphased: no significant row carries a direction")
  }
  crit(all(g$direction[!g$sig] == "none"), "[nonsig_none] rows that are not significant have the direction 'none'")
  s <- read.delim(file.path(RES, "summary_numbers_phaser.tsv"), stringsAsFactors = FALSE)
  crit(identical(names(s), c("sample", "genotypes", "genes_with_counts", "genes_tested", "sig_genes", "gw_phased_genes", "sig_unphased", "sig_both",
                             "sig_phaser_only", "sig_phaser_only_untested_unphased", "sig_unphased_only", "rho_used", "rho_own", "rho_cohort")) &&
       nrow(s) == nrow(sm) && all(s$sig_genes == vapply(s$sample, function(x) sum(g$sig[g$sample == x]), integer(1))),
       "[summary] summary_numbers_phaser.tsv: exact columns, one row per sample, sig_genes equal to the gene table")
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
