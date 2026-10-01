#!/usr/bin/env Rscript
# Stage 3 synthetic ground-truth ASE data for phASER: six outbred individuals with cis/trans phase configurations, low-depth genes,
# two-block genes, phased and unphased VCFs, direct count tables and emulated phASER gene tables (outbred_phase/).
# Usage: Rscript simulate_ase_stage3.R <outdir> <seed>
# Conventions: haplotype 1 / 2 are the left / right allele of the phased GT. At a het SNP alt_on = 1 puts ALT on haplotype 1 (GT 1|0),
# alt_on = 2 puts ALT on haplotype 2 (GT 0|1). h1 in truth_individual_genes.tsv is the EXPRESSION fraction of haplotype 1 (not of
# ALT or REF). truth_fragments: frag_* are fragments drawn per haplotype; cover_* are those covering >= 1 het SNP of the individual.
# Direct counts (ASEReadCounter columns: refCount = REF-allele fragments, altCount = ALT-allele fragments) and the emulated phASER
# gene tables (aCount = haplotype-1 side, bCount = haplotype-2 side for the phased table) ignore sequencing errors: they test the
# Rmds, not the aligner or phASER. Deterministic for a given seed (no dependence on time or locale).
suppressPackageStartupMessages(library(Biostrings))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: simulate_ase_stage3.R <outdir> <seed>")
outdir <- args[1]; seed <- as.integer(args[2]); set.seed(seed)
READLEN <- 100L; ERR <- 0.002; PHI_BIO <- 0.005; NGENES <- 80L; GLEN <- 1500000L; MIN_DEPTH <- 10L
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
writeLines(as.character(seed), file.path(outdir, "seed.txt"))
abs_out <- normalizePath(outdir)

# ---- genome and genes (Stage 2 block, 80 genes of 4-5 exons)
gdir <- file.path(outdir, "genome"); dir.create(gdir, showWarnings = FALSE)
sq <- sample(c("A", "C", "G", "T"), GLEN, replace = TRUE, prob = c(0.275, 0.225, 0.225, 0.275))
rows <- character(0); cur <- 2000L; genes <- list()
for (i in 1:NGENES) {
  ne <- sample(4:5, 1); el <- sample(150:300, ne, replace = TRUE); il <- sample(300:2000, ne - 1, replace = TRUE)
  st <- cur + c(0, cumsum(el[-ne] + il)); en <- st + el - 1L
  stopifnot(max(en) < GLEN - 2000L)
  gid <- sprintf("SYN%06d", i); tid <- paste0(gid, ".1")
  a <- function(extra = "") paste0('gene_id "', gid, '"; ', extra)
  rows <- c(rows,
    paste("chr1", "synthetic", "gene", min(st), max(en), ".", "+", ".", a(), sep = "\t"),
    paste("chr1", "synthetic", "transcript", min(st), max(en), ".", "+", ".", a(paste0('transcript_id "', tid, '";')), sep = "\t"),
    paste("chr1", "synthetic", "exon", st, en, ".", "+", ".", a(paste0('transcript_id "', tid, '"; exon_number "', seq_len(ne), '";')), sep = "\t"))
  genes[[gid]] <- list(id = gid, start = st, end = en, strand = "+", len = sum(en - st + 1))
  cur <- max(en) + sample(2000:4000, 1)
}
writeLines(rows, file.path(gdir, "genome.gtf"))
s <- paste(sq, collapse = "")
writeLines(c(">chr1", substring(s, seq(1, GLEN, 60), pmin(seq(60, GLEN + 59, 60), GLEN))), file.path(gdir, "genome.fa"))
gseq <- s; chrom <- "chr1"; gene_ids <- names(genes)
stopifnot(all(vapply(genes, function(x) x$len, 1) >= 600))

alt_of <- function(ref) sample(setdiff(c("A", "C", "G", "T"), ref), 1)
spliced <- function(g) paste(substring(gseq, genes[[g]]$start, genes[[g]]$end), collapse = "")
revcomp <- function(v) as.character(reverseComplement(DNAStringSet(v)))
add_err <- function(r) {
  ne <- rbinom(length(r), READLEN, ERR)
  for (i in which(ne > 0)) {
    p <- sample.int(READLEN, ne[i]); ch <- strsplit(r[i], "")[[1]]
    ch[p] <- vapply(ch[p], alt_of, character(1)); r[i] <- paste(ch, collapse = "")
  }
  r
}
write_fq <- function(path, name, seq, mate) {
  con <- gzfile(path, "wb"); on.exit(close(con))
  q <- strrep("I", READLEN)
  writeLines(paste0("@", name, "/", mate, "\n", seq, "\n+\n", q), con)
}
vcf_hdr <- function(sample = NULL) c("##fileformat=VCFv4.2", paste0("##contig=<ID=", chrom, ",length=", nchar(gseq), ">"),
  if (!is.null(sample)) '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')

wtab <- function(x, f) write.table(x, f, sep = "\t", quote = FALSE, row.names = FALSE)

frags2 <- function(hseq, n, snp_offs) {
  L <- nchar(hseq)
  if (n == 0) return(list(r1 = character(0), r2 = character(0), cover = matrix(FALSE, 0, length(snp_offs))))
  fl <- pmin(pmax(round(rnorm(n, 250, 30)), 150), 350, L); st <- floor(runif(n) * (L - fl + 1)) + 1
  r1 <- substring(hseq, st, st + READLEN - 1); r2 <- revcomp(substring(hseq, st + fl - READLEN, st + fl - 1))
  cover <- (outer(st, snp_offs, "<=") & outer(st + READLEN - 1, snp_offs, ">=")) |
           (outer(st + fl - READLEN, snp_offs, "<=") & outer(st + fl - 1, snp_offs, ">="))
  list(r1 = r1, r2 = r2, cover = cover)
}

# Block emulation for unphased phaser_gene_ae (copied from the plan)
emulate_blocks <- function(cover1, cover2, het) {   # fragment x SNP logical matrices of haplotype 1 / 2; het: logical over the gene's SNPs
  cv <- rbind(cover1, cover2)[, het, drop = FALSE]; hs <- c(rep(1L, nrow(cover1)), rep(2L, nrow(cover2)))
  cov_snp <- which(colSums(cv) > 0); if (length(cov_snp) == 0) return(NULL)
  comp <- seq_along(cov_snp)                        # union-find: SNPs joined when one fragment covers both
  find <- function(i) { while (comp[i] != i) i <- comp[i]; i }
  sub <- cv[, cov_snp, drop = FALSE]
  for (f in which(rowSums(sub) >= 2)) { s <- which(sub[f, ]); r <- find(s[1]); for (x in s[-1]) comp[find(x)] <- r }
  root <- vapply(seq_along(cov_snp), find, integer(1))
  blocks <- lapply(split(seq_along(cov_snp), root), function(ix) {
    anyc <- rowSums(sub[, ix, drop = FALSE]) > 0
    list(snps = which(het)[cov_snp[ix]], a = sum(anyc & hs == 1L), b = sum(anyc & hs == 2L)) })
  tot <- vapply(blocks, function(x) x$a + x$b, numeric(1)); first <- vapply(blocks, function(x) min(x$snps), numeric(1))
  list(n_blocks = length(blocks), best = blocks[[order(-tot, first)[1]]])
}

# ---- classes
glen <- vapply(genes, function(x) x$len, 1)
ord <- sample(gene_ids)
tb <- ord[glen[ord] >= 1100][1:2]
lw <- setdiff(ord, tb); lw <- lw[glen[lw] >= 1000][1:3]
if (anyNA(tb) || anyNA(lw)) stop("not enough long genes for two_block / hap_lowdepth")
rest <- setdiff(ord, c(tb, lw))
cls <- setNames(rep("null", NGENES), gene_ids)
cls[tb] <- "two_block"; cls[lw] <- "hap_lowdepth"; cls[rest[1:4]] <- "hap_strong"; cls[rest[5:8]] <- "hap_moderate"; cls[rest[9:14]] <- "null_linked"
spec <- data.frame(class = c("hap_strong", "hap_moderate", "hap_lowdepth", "two_block", "null_linked", "null"),
                   h_abs = c(0.3, 0.17, 0.4, 0.3, 0, 0), n_snps = c(4, 4, 8, 6, 4, 4), lambda = c(450, 450, 30, 450, 450, 450),
                   layout = c("cluster", "cluster", "spread", "two_cluster", "cluster", "random"), stringsAsFactors = FALSE)
rownames(spec) <- spec$class
tg <- data.frame(gene_id = gene_ids, class = unname(cls[gene_ids]), stringsAsFactors = FALSE)
tg$h_abs <- vapply(tg$class, function(k) format(spec[k, "h_abs"]), ""); tg$n_snps <- spec[tg$class, "n_snps"]
tg$lambda <- spec[tg$class, "lambda"]; tg$layout <- spec[tg$class, "layout"]

# ---- SNPs (transcript offsets == spliced offsets: plus-strand genes)
pairwise40 <- function(o) all(diff(sort(o)) >= 40)
snps <- list()
for (g in gene_ids) {
  x <- genes[[g]]; len <- x$len; lay <- spec[cls[[g]], "layout"]; cl <- rep(0L, spec[cls[[g]], "n_snps"])
  if (lay == "cluster") {
    repeat { w0 <- sample.int(len - 299, 1); off <- sort(w0 + sample.int(300, 4) - 1); if (pairwise40(off)) break }
  } else if (lay == "spread") {
    o0 <- sample(61:(len - 900), 1); off <- o0 + 120L * (0:7)
  } else if (lay == "two_cluster") {
    repeat { o1 <- sort(sample.int(150, 3)); if (pairwise40(o1)) break }
    repeat { o2 <- sort(len - 150L + sample.int(150, 3)); if (pairwise40(o2)) break }
    off <- c(o1, o2); cl <- rep(1:2, each = 3)
  } else {
    repeat { off <- sort(sample.int(len, 4)); if (pairwise40(off)) break }
  }
  cum <- c(0, cumsum(x$end - x$start + 1))
  pos <- vapply(off, function(o) { i <- max(which(cum < o)); x$start[i] + o - cum[i] - 1 }, numeric(1))
  ref <- substring(gseq, pos, pos)
  alt <- vapply(ref, alt_of, character(1))
  snps[[g]] <- data.frame(chrom = chrom, pos = as.integer(pos), ref = unname(ref), alt = unname(alt), gene_id = g, tx_off = off, cluster = cl, stringsAsFactors = FALSE)
}
snp <- do.call(rbind, snps); rownames(snp) <- NULL
stopifnot(!anyDuplicated(snp$pos))
snp <- snp[order(snp$pos), ]; rownames(snp) <- NULL
gidx <- lapply(setNames(gene_ids, gene_ids), function(g) which(snp$gene_id == g))   # ascending pos == ascending offset
stopifnot(all(vapply(genes, function(x) x$strand, "") == "+"))
haplo <- function(g, alleles) { s <- strsplit(spliced(g), "")[[1]]; s[snp$tx_off[gidx[[g]]]] <- alleles; paste(s, collapse = "") }

# ---- individuals: genotypes, phase (alt_on), h1
inds <- paste0("ind", 1:6); NI <- length(inds)
gtm <- matrix(NA_character_, nrow(snp), NI); aon <- matrix(NA_integer_, nrow(snp), NI)
h1m <- matrix(0.5, NGENES, NI, dimnames = list(gene_ids, inds)); cluster_gene <- gene_ids[cls %in% c("hap_strong", "hap_moderate", "null_linked")]
is_cis <- matrix(NA, NGENES, NI, dimnames = list(gene_ids, inds))
for (k in seq_len(NI)) for (g in gene_ids) {
  i <- gidx[[g]]; n <- length(i)
  gt <- if (cls[[g]] == "null") sample(c("0/1", "0/0", "1/1"), n, replace = TRUE, prob = c(0.6, 0.2, 0.2)) else rep("0/1", n)
  if (g %in% cluster_gene) {
    cis <- runif(1) < 0.5; is_cis[g, k] <- cis
    if (cis) a <- rep(sample(1:2, 1), n) else repeat { a <- sample(1:2, n, replace = TRUE); if (length(unique(a)) == 2) break }
  } else a <- sample(1:2, n, replace = TRUE)
  gtm[i, k] <- gt; aon[i, k] <- ifelse(gt == "0/1", a, NA_integer_)
  hab <- spec[cls[[g]], "h_abs"]
  h1m[g, k] <- if (hab > 0) 0.5 + hab * sample(c(-1, 1), 1) else 0.5
}
stopifnot(all(colSums(is_cis[cluster_gene, , drop = FALSE]) >= 2), all(colSums(!is_cis[cluster_gene, , drop = FALSE]) >= 2))

od <- file.path(outdir, "outbred_phase"); dir.create(file.path(od, "direct_counts"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(od, "direct_gene_ae_phased"), showWarnings = FALSE); dir.create(file.path(od, "direct_gene_ae_unphased"), showWarnings = FALSE)
wtab(tg, file.path(od, "truth_genes.tsv")); wtab(snp, file.path(od, "truth_snps.tsv"))
tig <- do.call(rbind, lapply(seq_len(NI), function(k) do.call(rbind, lapply(gene_ids, function(g)
  data.frame(individual = inds[k], gene_id = g, h1 = format(h1m[g, k], digits = 12, trim = TRUE), n_het = sum(gtm[gidx[[g]], k] == "0/1"), stringsAsFactors = FALSE)))))
wtab(tig, file.path(od, "truth_individual_genes.tsv"))
tph <- do.call(rbind, lapply(seq_len(NI), function(k) data.frame(individual = inds[k], chrom = snp$chrom, pos = snp$pos, gene_id = snp$gene_id,
  gt = gtm[, k], alt_on = ifelse(is.na(aon[, k]), "NA", as.character(aon[, k])), stringsAsFactors = FALSE)))
wtab(tph, file.path(od, "truth_phase.tsv"))

for (k in seq_len(NI)) {
  gt <- gtm[, k]; a <- aon[, k]
  ph <- ifelse(gt == "0/0", "0|0", ifelse(gt == "1/1", "1|1", ifelse(a == 2, "0|1", "1|0")))   # left allele = haplotype 1
  for (pf in c("", ".phased")) {
    vc <- data.frame(CHROM = snp$chrom, POS = snp$pos, ID = ".", REF = snp$ref, ALT = snp$alt, QUAL = ".", FILTER = ".", INFO = ".", FORMAT = "GT",
                     S = if (pf == "") gt else ph, stringsAsFactors = FALSE)
    names(vc)[10] <- inds[k]
    for (nm in c("all", "het")) {
      f <- file.path(od, paste0(inds[k], ".", nm, pf, ".vcf"))
      writeLines(c(vcf_hdr(inds[k]), paste0("#", paste(names(vc), collapse = "\t"))), f)
      write.table(if (nm == "het") vc[gt == "0/1", ] else vc, f, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    }
  }
}

# ---- samples: fragments, reads, direct counts, emulated phASER gene tables
lg2 <- function(a, b) if (b == 0 && a == 0) "inf" else if (b == 0) "inf" else if (a == 0) "-inf" else format(log2(a / b), digits = 6)
gene_row <- function(g, a, b, nv, vars, gwp, s) {
  paste(c("chr1", sprintf("%d", as.integer(min(genes[[g]]$start) - 1)), sprintf("%d", as.integer(max(genes[[g]]$end))), g, a, b, a + b, lg2(a, b), nv, vars, gwp, s), collapse = "\t")
}
GH <- "contig\tstart\tstop\tname\taCount\tbCount\ttotalCount\tlog2_aFC\tn_variants\tvariants\tgw_phased\tbam"
ss <- paste0("obp_s", 1:6); tfa <- list(); flips <- list()
for (k in seq_len(NI)) {
  s <- ss[k]; gt <- gtm[, k]; a <- aon[, k]
  r1 <- character(0); r2 <- character(0); refc <- numeric(nrow(snp)); altc <- numeric(nrow(snp)); tr <- list()
  rp <- GH; ru <- GH; fl <- list()
  for (g in gene_ids) {
    i <- gidx[[g]]; het <- gt[i] == "0/1"; h1 <- h1m[g, k]; lam <- spec[cls[[g]], "lambda"]; offs <- snp$tx_off[i]
    alleles1 <- ifelse(gt[i] == "1/1", snp$alt[i], ifelse(gt[i] == "0/0", snp$ref[i], ifelse(a[i] %in% 1, snp$alt[i], snp$ref[i])))
    alleles2 <- ifelse(gt[i] == "1/1", snp$alt[i], ifelse(gt[i] == "0/0", snp$ref[i], ifelse(a[i] %in% 1, snp$ref[i], snp$alt[i])))
    hs1 <- haplo(g, alleles1); hs2 <- haplo(g, alleles2)
    ps <- rbeta(1, h1 * (1 - PHI_BIO) / PHI_BIO, (1 - h1) * (1 - PHI_BIO) / PHI_BIO)
    n <- rpois(1, lam); n1 <- rbinom(1, n, ps)
    tries <- 0
    repeat {
      f1 <- frags2(hs1, n1, offs); f2 <- frags2(hs2, n - n1, offs)
      if (cls[[g]] != "hap_lowdepth") break
      tries <- tries + 1; if (tries > 1000) stop("hap_lowdepth redraw failed for ", g)
      if (all(colSums(f1$cover) + colSums(f2$cover) < MIN_DEPTH) && sum(rowSums(f1$cover) > 0) + sum(rowSums(f2$cover) > 0) >= 18) break
    }
    r1 <- c(r1, f1$r1, f2$r1); r2 <- c(r2, f1$r2, f2$r2)
    c1 <- colSums(f1$cover); c2 <- colSums(f2$cover)
    refc[i] <- ifelse(het, ifelse(a[i] %in% 1, c2, c1), 0); altc[i] <- ifelse(het, ifelse(a[i] %in% 1, c1, c2), 0)
    ch1 <- sum(rowSums(f1$cover[, het, drop = FALSE]) > 0); ch2 <- sum(rowSums(f2$cover[, het, drop = FALSE]) > 0)
    msd <- if (any(het)) max((c1 + c2)[het]) else 0
    tr[[g]] <- data.frame(sample = s, gene_id = g, frag_h1 = n1, frag_h2 = n - n1, cover_h1 = ch1, cover_h2 = ch2, max_snp_depth = msd, stringsAsFactors = FALSE)
    vname <- function(ix) paste0("chr1_", snp$pos[i][ix], "_", snp$ref[i][ix], "_", snp$alt[i][ix], collapse = ",")
    cov_ix <- which(het & (c1 + c2) > 0)
    if (ch1 + ch2 == 0) { u <- gene_row(g, 0, 0, 0, "", 1, s); rp <- c(rp, u); ru <- c(ru, u); fl[[g]] <- data.frame(sample = s, gene_id = g, n_blocks = 0, best_block_n_snps = 0, flipped = "NA") ; next }
    rp <- c(rp, gene_row(g, ch1, ch2, length(cov_ix), vname(cov_ix), 1, s))
    eb <- emulate_blocks(f1$cover, f2$cover, het); flipped <- runif(1) < 0.5
    ba <- eb$best$a; bb <- eb$best$b; if (flipped) { t <- ba; ba <- bb; bb <- t }
    ru <- c(ru, gene_row(g, ba, bb, length(eb$best$snps), vname(eb$best$snps), 0, s))
    fl[[g]] <- data.frame(sample = s, gene_id = g, n_blocks = eb$n_blocks, best_block_n_snps = length(eb$best$snps), flipped = ifelse(flipped, "TRUE", "FALSE"))
  }
  r1 <- add_err(r1); r2 <- add_err(r2); idx <- seq_along(r1)
  write_fq(file.path(od, paste0(s, "_1.fastq.gz")), paste0(s, "_", idx), r1, 1)
  write_fq(file.path(od, paste0(s, "_2.fastq.gz")), paste0(s, "_", idx), r2, 2)
  tot <- refc + altc; kp <- tot >= 1
  counts <- data.frame(contig = snp$chrom, position = snp$pos, variantID = ".", refAllele = snp$ref, altAllele = snp$alt,
                       refCount = refc, altCount = altc, totalCount = tot, lowMAPQDepth = 0, lowBaseQDepth = 0, rawDepth = tot,
                       otherBases = 0, improperPairs = 0, stringsAsFactors = FALSE)[kp, ]
  tf <- file.path(od, "direct_counts", paste0(s, ".table")); wtab(counts, tf)
  file.copy(tf, file.path(od, "direct_counts", paste0(s, ".unfiltered.table")))
  writeLines(c("vW\tvA\tn", paste("1", "1", sum(counts$refCount), sep = "\t"), paste("1", "2", sum(counts$altCount), sep = "\t"), "none\t-\t0"),
             file.path(od, "direct_counts", paste0(s, ".wasp_stats.tsv")))
  writeLines(rp, file.path(od, "direct_gene_ae_phased", paste0(s, ".gene_ae.txt")))
  writeLines(ru, file.path(od, "direct_gene_ae_unphased", paste0(s, ".gene_ae.txt")))
  tfa[[k]] <- do.call(rbind, tr); flips[[k]] <- do.call(rbind, fl)
}
wtab(do.call(rbind, tfa), file.path(od, "truth_fragments.tsv")); wtab(do.call(rbind, flips), file.path(od, "direct_flips.tsv"))
write.csv(data.frame(sample = ss, fastq_1 = file.path(abs_out, "outbred_phase", paste0(ss, "_1.fastq.gz")),
  fastq_2 = file.path(abs_out, "outbred_phase", paste0(ss, "_2.fastq.gz")), condition = "ctrl", cross_direction = "NA", individual = inds),
  file.path(od, "samples.csv"), row.names = FALSE, quote = FALSE)
cat("done\n")
