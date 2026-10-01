#!/usr/bin/env Rscript
# Usage: Rscript check_phase_reads_stage3.R <outdir>
# From the FASTQs alone: every read pair that carries alleles at two or more heterozygous SNPs of one gene (exact 41 bp allele
# windows, one window per SNP and allele) must carry alleles of ONE haplotype of truth_phase.tsv (cis or trans as planted).
# Partition: in every sample at least 99% of such pairs (correct data >= 99.6%; one wrong alt_on in one gene: 98.1%) agree with the truth, and at least 50 pairs link SNPs in trans.
# Orientation: in every sample, for every planted cell (h1 != 0.5: hap_strong, hap_moderate, hap_lowdepth, two_block) with >= 20 read pairs carrying one haplotype of the
# gene, the majority haplotype equals the one favoured by h1 (h1 > 0.5: haplotype 1 = left allele of the phased GT); at least 8 of the 13 planted cells must have >= 20 pairs.
# Orientation of null and null_linked genes (h1 = 0.5) cannot be observed from reads or tables: tests of haplotype direction must use planted genes only.
suppressPackageStartupMessages(library(Biostrings))
O <- file.path(commandArgs(trailingOnly = TRUE)[1], "outbred_phase"); W <- 20L; MINREADS <- 20L; MINCELLS <- 8L
gseq <- as.character(readDNAStringSet(file.path(O, "../genome/genome.fa"))[[1]])
ex <- read.delim(file.path(O, "../genome/genome.gtf"), header = FALSE, quote = "", stringsAsFactors = FALSE)
ex <- ex[ex$V3 == "exon", ]; ex$gene_id <- sub('.*gene_id "([^"]+)".*', "\\1", ex$V9)
ts <- read.delim(file.path(O, "truth_snps.tsv"), stringsAsFactors = FALSE)
tp <- read.delim(file.path(O, "truth_phase.tsv"), stringsAsFactors = FALSE, colClasses = c(gt = "character", alt_on = "character"))
sm <- read.csv(file.path(O, "samples.csv"), stringsAsFactors = FALSE)
tg <- read.delim(file.path(O, "truth_genes.tsv"), stringsAsFactors = FALSE)
th <- read.delim(file.path(O, "truth_individual_genes.tsv"), stringsAsFactors = FALSE)
spl <- lapply(split(ex, ex$gene_id), function(e) { e <- e[order(e$V4), ]; strsplit(paste(substring(gseq, e$V4, e$V5), collapse = ""), "")[[1]] })
allok <- TRUE
for (i in seq_len(nrow(sm))) {
  ph <- tp[tp$individual == sm$individual[i] & tp$gt == "0/1", c("pos", "alt_on")]
  s <- merge(ts, ph, by = "pos")
  w <- do.call(rbind, lapply(seq_len(nrow(s)), function(k) {
    sp <- spl[[s$gene_id[k]]]; o <- s$tx_off[k]
    if (o <= W || o + W > length(sp)) return(NULL)
    hr <- sp; hr[o] <- s$ref[k]; ha <- sp; ha[o] <- s$alt[k]
    data.frame(snp = k, allele = c("ref", "alt"), seq = c(paste(hr[(o - W):(o + W)], collapse = ""), paste(ha[(o - W):(o + W)], collapse = "")),
               stringsAsFactors = FALSE)
  }))
  pd <- PDict(DNAStringSet(w$seq))
  h1 <- vwhichPDict(pd, readDNAStringSet(sm$fastq_1[i], format = "fastq"))
  h2 <- vwhichPDict(pd, reverseComplement(readDNAStringSet(sm$fastq_2[i], format = "fastq")))
  hits <- mapply(function(x, y) unique(c(x, y)), h1, h2, SIMPLIFY = FALSE)
  total <- 0; agree <- 0; trans <- 0
  for (q in which(lengths(hits) >= 2)) {
    hw <- w[hits[[q]], ]
    if (anyDuplicated(hw$snp) || length(unique(s$gene_id[hw$snp])) != 1) next   # both alleles of one SNP (error) or two genes
    hap <- ifelse(hw$allele == "alt", s$alt_on[hw$snp], ifelse(s$alt_on[hw$snp] == "1", "2", "1"))
    total <- total + 1; agree <- agree + (length(unique(hap)) == 1)
    if (length(unique(s$alt_on[hw$snp])) > 1) trans <- trans + 1
  }
  # orientation: for each planted cell (h1 != 0.5) the haplotype that most read pairs carry (via alt_on: ALT allele = haplotype alt_on,
  # REF allele = the other) must be the one the truth h1 favours. A pair counts when all its allele windows are in one gene and on one haplotype.
  hapw <- ifelse(w$allele == "alt", s$alt_on[w$snp], ifelse(s$alt_on[w$snp] == "1", "2", "1"))
  pi <- rep(seq_along(hits), lengths(hits)); wi <- unlist(hits)
  d <- unique(data.frame(pi = pi, gene = s$gene_id[w$snp[wi]], hap = hapw[wi], stringsAsFactors = FALSE))
  d <- d[!(d$pi %in% d$pi[duplicated(d$pi)]), ]
  cnt <- table(factor(d$gene, levels = tg$gene_id), factor(d$hap, levels = c("1", "2")))
  hh <- th[th$individual == sm$individual[i], ]; hh <- hh[match(tg$gene_id, hh$gene_id), ]
  cells <- which(tg$class %in% c("hap_strong", "hap_moderate", "hap_lowdepth", "two_block"))
  nn <- cnt[cells, 1] + cnt[cells, 2]; el <- nn >= MINREADS
  want <- ifelse(hh$h1[cells] > 0.5, 1L, 2L); got <- ifelse(cnt[cells, 1] > cnt[cells, 2], 1L, 2L)
  oagree <- sum(el & want == got); oel <- sum(el)
  ook <- oel >= MINCELLS && oagree == oel
  cat(sprintf("%s (%s): orientation %d of %d planted cells with >= %d read pairs agree with sign(h1 - 0.5) (%d cells below %d pairs, smallest cell %d): %s\n",
              sm$sample[i], sm$individual[i], oagree, oel, MINREADS, sum(!el), MINREADS, min(nn), if (ook) "orientation ok" else "orientation FAIL"))
  okk <- total > 0 && agree / total >= 0.99 && trans >= 50
  cat(sprintf("%s (%s): %d read pairs link two or more SNPs, %.4f agree with truth_phase, %d link SNPs in trans: %s\n",
              sm$sample[i], sm$individual[i], total, if (total > 0) agree / total else NA, trans, if (okk) "partition ok" else "partition FAIL"))
  allok <- allok && okk && ook
}
cat(if (allok) "PHASE READS PASS\n" else "PHASE READS FAIL\n"); quit(status = if (allok) 0 else 1)
