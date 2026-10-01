#!/usr/bin/env Rscript
# Usage: Rscript check_phase_reads_stage3.R <outdir>
# From the FASTQs alone: every read pair that carries alleles at two or more heterozygous SNPs of one gene (exact 41 bp allele
# windows, one window per SNP and allele) must carry alleles of ONE haplotype of truth_phase.tsv (cis or trans as planted).
# Pass: in every sample at least 99% of such pairs (correct data: >= 99.6%; one wrong alt_on in one gene: 98.1%) agree with the truth, and at least 50 pairs link SNPs in trans.
suppressPackageStartupMessages(library(Biostrings))
O <- file.path(commandArgs(trailingOnly = TRUE)[1], "outbred_phase"); W <- 20L
gseq <- as.character(readDNAStringSet(file.path(O, "../genome/genome.fa"))[[1]])
ex <- read.delim(file.path(O, "../genome/genome.gtf"), header = FALSE, quote = "", stringsAsFactors = FALSE)
ex <- ex[ex$V3 == "exon", ]; ex$gene_id <- sub('.*gene_id "([^"]+)".*', "\\1", ex$V9)
ts <- read.delim(file.path(O, "truth_snps.tsv"), stringsAsFactors = FALSE)
tp <- read.delim(file.path(O, "truth_phase.tsv"), stringsAsFactors = FALSE, colClasses = c(gt = "character", alt_on = "character"))
sm <- read.csv(file.path(O, "samples.csv"), stringsAsFactors = FALSE)
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
  okk <- total > 0 && agree / total >= 0.99 && trans >= 50
  cat(sprintf("%s (%s): %d read pairs link two or more SNPs, %.4f agree with truth_phase, %d link SNPs in trans: %s\n",
              sm$sample[i], sm$individual[i], total, if (total > 0) agree / total else NA, trans, if (okk) "ok" else "FAIL"))
  allok <- allok && okk
}
cat(if (allok) "PHASE READS PASS\n" else "PHASE READS FAIL\n"); quit(status = if (allok) 0 else 1)
