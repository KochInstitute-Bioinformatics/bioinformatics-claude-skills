#!/usr/bin/env Rscript
# Usage: Rscript check_biology.R <outdir> <mode: f1|outbred-unused>  (synthetic genome, plus-strand genes)
# Re-derives, from the F1 FASTQs alone, the ALT fraction at each gene's SNPs by exact matching of
# 41 bp ref/alt haplotype windows centred on each SNP, and compares it with truth_genes.tsv p_alt.
suppressPackageStartupMessages({ library(Biostrings); library(rtracklayer) })
D <- commandArgs(trailingOnly = TRUE)[1]; W <- 20L
gseq <- as.character(readDNAStringSet(file.path(D, "genome/genome.fa"))[[1]])
ex <- as.data.frame(import(file.path(D, "genome/genome.gtf"))); ex <- ex[ex$type == "exon", ]
tg <- read.delim(file.path(D, "f1/truth_genes.tsv")); ts <- read.delim(file.path(D, "f1/truth_snps.tsv"), stringsAsFactors = FALSE)
sm <- read.csv(file.path(D, "f1/samples.csv"))
wins <- list()   # each: gene, allele, seq
for (g in tg$gene_id) {
  e <- ex[ex$gene_id == g, ]; e <- e[order(e$start), ]
  stopifnot(all(e$strand == "+"))
  sp <- strsplit(paste(substring(gseq, e$start, e$end), collapse = ""), "")[[1]]
  cum <- c(0, cumsum(e$end - e$start + 1)); s <- ts[ts$gene_id == g, ]
  off <- vapply(s$pos, function(p) { i <- which(p >= e$start & p <= e$end); cum[i] + p - e$start[i] + 1 }, numeric(1))
  hr <- sp; hr[off] <- s$ref; ha <- sp; ha[off] <- s$alt
  for (k in seq_along(off)) if (off[k] > W && off[k] + W <= length(sp)) {
    wins[[length(wins) + 1]] <- data.frame(gene = g, allele = "ref", seq = paste(hr[(off[k] - W):(off[k] + W)], collapse = ""))
    wins[[length(wins) + 1]] <- data.frame(gene = g, allele = "alt", seq = paste(ha[(off[k] - W):(off[k] + W)], collapse = ""))
  }
}
wins <- do.call(rbind, wins)
pd <- PDict(DNAStringSet(wins$seq))
tot <- matrix(0, nrow(wins), 1)
for (i in seq_len(nrow(sm))) {
  r <- c(readDNAStringSet(sm$fastq_1[i], format = "fastq"), reverseComplement(readDNAStringSet(sm$fastq_2[i], format = "fastq")))
  # mate1 forward reads; mate2 reverse-complemented gives forward orientation too
  tot[, 1] <- tot[, 1] + rowSums(vcountPDict(pd, r))
}
wins$n <- tot[, 1]
res <- do.call(rbind, lapply(tg$gene_id, function(g) {
  w <- wins[wins$gene == g, ]; a <- sum(w$n[w$allele == "alt"]); r <- sum(w$n[w$allele == "ref"])
  data.frame(gene_id = g, class = tg$class[tg$gene_id == g], planted = tg$p_alt[tg$gene_id == g], ref_reads = r, alt_reads = a,
             observed = round(a / (a + r), 3), se = round(sqrt(0.25 / (a + r)), 3))
}))
res$abs_diff <- round(abs(res$observed - res$planted), 3)
print(res, row.names = FALSE)
ok <- all(res$abs_diff <= 0.1)
cat(if (ok) "BIOLOGY PASS (all |observed - planted| <= 0.1)\n" else "BIOLOGY FAIL\n")
quit(status = if (ok) 0 else 1)
