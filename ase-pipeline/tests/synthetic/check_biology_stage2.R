#!/usr/bin/env Rscript
# Usage: Rscript check_biology_stage2.R <outdir>   (synthetic genome, plus-strand genes; F1 reciprocal x condition data)
# Re-derives, from the 12 F1 FASTQs alone, the ALT fraction at each gene's SNPs by exact matching of
# 41 bp ref/alt haplotype windows centred on each SNP, and compares it with 1 - the realised strain-A
# fragment fraction (truth_fragments.tsv pooled over all samples). Pass: |observed - expected| <= 0.05 per gene.
suppressPackageStartupMessages({ library(Biostrings) })
D <- file.path(commandArgs(trailingOnly = TRUE)[1], "f1_recip"); W <- 20L
gseq <- as.character(readDNAStringSet(file.path(D, "../genome/genome.fa"))[[1]])
ex <- read.delim(file.path(D, "../genome/genome.gtf"), header = FALSE, quote = "", stringsAsFactors = FALSE)
ex <- ex[ex$V3 == "exon", ]; ex$gene_id <- sub('.*gene_id "([^"]+)".*', "\\1", ex$V9)
tg <- read.delim(file.path(D, "truth_genes.tsv")); ts <- read.delim(file.path(D, "truth_snps.tsv"), stringsAsFactors = FALSE)
tf <- read.delim(file.path(D, "truth_fragments.tsv"), stringsAsFactors = FALSE)
sm <- read.csv(file.path(D, "samples.csv"))
wins <- list()
for (g in tg$gene_id) {
  e <- ex[ex$gene_id == g, ]; e <- e[order(e$V4), ]
  stopifnot(all(e$V7 == "+"))
  sp <- strsplit(paste(substring(gseq, e$V4, e$V5), collapse = ""), "")[[1]]
  s <- ts[ts$gene_id == g, ]; off <- s$tx_off
  stopifnot(all(sp[off] == s$ref))
  hr <- sp; hr[off] <- s$ref; ha <- sp; ha[off] <- s$alt
  for (k in seq_along(off)) if (off[k] > W && off[k] + W <= length(sp)) {
    wins[[length(wins) + 1]] <- data.frame(gene = g, allele = "ref", seq = paste(hr[(off[k] - W):(off[k] + W)], collapse = ""))
    wins[[length(wins) + 1]] <- data.frame(gene = g, allele = "alt", seq = paste(ha[(off[k] - W):(off[k] + W)], collapse = ""))
  }
}
wins <- do.call(rbind, wins)
pd <- PDict(DNAStringSet(wins$seq))
tot <- numeric(nrow(wins))
for (i in seq_len(nrow(sm))) {
  r <- c(readDNAStringSet(sm$fastq_1[i], format = "fastq"), reverseComplement(readDNAStringSet(sm$fastq_2[i], format = "fastq")))
  tot <- tot + rowSums(vcountPDict(pd, r))
}
wins$n <- tot
res <- do.call(rbind, lapply(tg$gene_id, function(g) {
  w <- wins[wins$gene == g, ]; a <- sum(w$n[w$allele == "alt"]); r <- sum(w$n[w$allele == "ref"])
  f <- tf[tf$gene_id == g, ]; realised_A <- sum(f$frag_A) / sum(f$frag_A + f$frag_B)
  data.frame(gene_id = g, class = tg$class[tg$gene_id == g], expected_alt = round(1 - realised_A, 3), ref_reads = r, alt_reads = a,
             observed = round(a / (a + r), 3))
}))
res$abs_diff <- round(abs(res$observed - res$expected_alt), 3)
print(res, row.names = FALSE)
ok <- nrow(res) == 60 && !anyNA(res$abs_diff) && all(res$abs_diff <= 0.05)
cat(if (ok) "BIOLOGY PASS (all |observed - realised| <= 0.05)\n" else "BIOLOGY FAIL\n")
quit(status = if (ok) 0 else 1)
