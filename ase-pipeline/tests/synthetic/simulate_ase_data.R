#!/usr/bin/env Rscript
# Synthetic ground-truth ASE data generator (F1 cross and outbred designs).
# Usage: Rscript simulate_ase_data.R <genome.fasta|SYNTHETIC> <genome.gtf|SYNTHETIC> <outdir> <seed>
# If both genome arguments are the literal SYNTHETIC, a deterministic synthetic genome
# (chr1, 300 kb, ~45% GC, 20 plus-strand genes) is generated into <outdir>/genome/ and used.
suppressPackageStartupMessages({ library(Biostrings); library(rtracklayer) })
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: simulate_ase_data.R <genome.fasta|SYNTHETIC> <genome.gtf|SYNTHETIC> <outdir> <seed>")
fasta <- args[1]; gtf <- args[2]; outdir <- args[3]; seed <- as.integer(args[4])
set.seed(seed)
SYNTH <- identical(fasta, "SYNTHETIC") && identical(gtf, "SYNTHETIC")
MIN_LEN <- if (SYNTH) 400L else 300L   # nf-core test genome has too few loci for 400
MAX_GENES <- 12L; READLEN <- 100L; ERR <- 0.002; LAMBDA <- 450
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
writeLines(as.character(seed), file.path(outdir, "seed.txt"))

if (SYNTH) {
  gdir <- file.path(outdir, "genome"); dir.create(gdir, showWarnings = FALSE)
  GLEN <- 300000L
  sq <- sample(c("A", "C", "G", "T"), GLEN, replace = TRUE, prob = c(0.275, 0.225, 0.225, 0.275))
  rows <- character(0); cur <- 2000L
  for (i in 1:20) {
    ne <- sample(3:5, 1); el <- sample(120:300, ne, replace = TRUE); il <- sample(300:2000, ne - 1, replace = TRUE)
    st <- cur + c(0, cumsum(el[-ne] + il)); en <- st + el - 1L
    stopifnot(max(en) < GLEN - 2000L)
    gid <- sprintf("SYN%06d", i); tid <- paste0(gid, ".1")
    a <- function(extra = "") paste0('gene_id "', gid, '"; ', extra)
    rows <- c(rows,
      paste("chr1", "synthetic", "gene", min(st), max(en), ".", "+", ".", a(), sep = "\t"),
      paste("chr1", "synthetic", "transcript", min(st), max(en), ".", "+", ".", a(paste0('transcript_id "', tid, '";')), sep = "\t"),
      paste("chr1", "synthetic", "exon", st, en, ".", "+", ".", a(paste0('transcript_id "', tid, '"; exon_number "', seq_len(ne), '";')), sep = "\t"))
    cur <- max(en) + sample(2000:4000, 1)
  }
  writeLines(rows, file.path(gdir, "genome.gtf"))
  s <- paste(sq, collapse = "")
  writeLines(c(">chr1", substring(s, seq(1, GLEN, 60), pmin(seq(60, GLEN + 59, 60), GLEN))), file.path(gdir, "genome.fa"))
  fasta <- file.path(gdir, "genome.fa"); gtf <- file.path(gdir, "genome.gtf")
}


genome <- readDNAStringSet(fasta)
chrom <- names(genome)[1]; chrom <- sub("\\s.*$", "", chrom)
gseq <- as.character(genome[[1]])
ex <- as.data.frame(import(gtf))
ex <- ex[ex$type == "exon", ]

# ---- 1. genes: first transcript (file order), spliced length filter
genes <- list()
for (g in unique(ex$gene_id)) {
  e <- ex[ex$gene_id == g, ]
  e <- e[e$transcript_id == e$transcript_id[1], ]
  e <- e[order(e$start), ]
  len <- sum(e$end - e$start + 1)
  if (len >= MIN_LEN) genes[[g]] <- list(id = g, start = e$start, end = e$end, strand = as.character(e$strand[1]), len = len)
}
gene_ids <- names(genes)
# drop genes whose exons overlap an earlier kept gene (their reads/SNPs cannot be separated)
keep <- character(0)
for (g in gene_ids) {
  ov <- FALSE
  for (h in keep) ov <- ov || any(outer(genes[[g]]$start, genes[[h]]$end, "<=") & outer(genes[[g]]$end, genes[[h]]$start, ">="))
  if (!ov) keep <- c(keep, g) else cat("dropped (exon overlap):", g, "\n")
}
gene_ids <- head(keep, MAX_GENES)
genes <- genes[gene_ids]
cat("kept genes:", length(genes), "\n")
stopifnot(length(genes) >= 6)

# positions covered by exons of more than one kept gene are not eligible for SNPs
cov <- unlist(lapply(genes, function(x) unlist(Map(seq.int, x$start, x$end))))
ambiguous <- as.integer(names(which(table(cov) > 1)))
comp <- c(A = "T", C = "G", G = "C", T = "A")
alt_of <- function(ref) sample(setdiff(c("A", "C", "G", "T"), ref), 1)

# ---- classes: seeded assignment
ord <- sample(gene_ids)
bias_ok <- ord[vapply(ord, function(g) any(genes[[g]]$end - genes[[g]]$start + 1 >= 60), logical(1))]
bias_gene <- bias_ok[1]
rest <- setdiff(ord, bias_gene)
imb_p <- if (length(gene_ids) >= 12) c(0.70, 0.70, 0.85, 0.30) else c(0.70, 0.30)  # full quota needs 12 genes
n_imb <- length(imb_p)
imb_genes <- rest[seq_len(n_imb)]
cls <- setNames(rep("null", length(gene_ids)), gene_ids); pal <- setNames(rep(0.5, length(gene_ids)), gene_ids)
cls[bias_gene] <- "bias"
cls[imb_genes] <- "imbalanced"; pal[imb_genes] <- imb_p

# ---- 2. SNPs
snps <- list()
for (g in gene_ids) {
  x <- genes[[g]]
  pos_all <- setdiff(unlist(Map(seq.int, x$start, x$end)), ambiguous)
  if (g == bias_gene) {
    big <- which(x$end - x$start + 1 >= 60)
    i <- big[sample.int(length(big), 1)]
    w0 <- x$start[i] + sample.int(x$end[i] - x$start[i] + 1 - 59, 1) - 1
    pos <- sort(w0 + sample.int(60, 5) - 1)
  } else {
    repeat {
      pos <- sort(sample(pos_all, 4))
      if (all(diff(pos) >= 40)) break
    }
  }
  ref <- substring(gseq, pos, pos)
  alt <- vapply(ref, alt_of, character(1))
  snps[[g]] <- data.frame(chrom = chrom, pos = pos, ref = unname(ref), alt = unname(alt), gene_id = g, stringsAsFactors = FALSE)
}
snp <- do.call(rbind, snps); rownames(snp) <- NULL
stopifnot(!anyDuplicated(snp$pos))
snp <- snp[order(snp$pos), ]; rownames(snp) <- NULL

# spliced sequence (genomic + orientation) for a gene and the offset of each SNP in it
spliced <- function(g) paste(substring(gseq, genes[[g]]$start, genes[[g]]$end), collapse = "")
snp_off <- function(g) {
  x <- genes[[g]]; s <- snp[snp$gene_id == g, ]
  cum <- c(0, cumsum(x$end - x$start + 1))
  vapply(s$pos, function(p) { i <- which(p >= x$start & p <= x$end); cum[i] + p - x$start[i] + 1 }, numeric(1))
}
revcomp <- function(v) as.character(reverseComplement(DNAStringSet(v)))
# haplotype transcript in transcript orientation, given alleles (character vector per SNP of gene)
haplo <- function(g, alleles) {
  s <- strsplit(spliced(g), "")[[1]]
  s[snp_off(g)] <- alleles
  s <- paste(s, collapse = "")
  if (genes[[g]]$strand == "-") revcomp(s) else s
}

# ---- 6. fragments -> read pairs
frags <- function(hseq, n) {
  L <- nchar(hseq)
  if (n == 0) return(list(r1 = character(0), r2 = character(0)))
  fl <- pmin(pmax(round(rnorm(n, 250, 30)), 150), 350, L)
  st <- floor(runif(n) * (L - fl + 1)) + 1
  r1 <- substring(hseq, st, st + READLEN - 1)
  r2 <- revcomp(substring(hseq, st + fl - READLEN, st + fl - 1))
  list(r1 = r1, r2 = r2)
}
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

# simulate one sample; hapA/hapB: lists (per gene) of haplotype transcripts (REF-side, ALT-side); has_alt: gene can be imbalanced
simulate_sample <- function(sample, hapR, hapA, p_use, dirp) {
  r1 <- character(0); r2 <- character(0)
  for (g in gene_ids) {
    n <- rpois(1, LAMBDA); na <- rbinom(1, n, p_use[[g]])
    fa <- frags(hapA[[g]], na); fr <- frags(hapR[[g]], n - na)
    r1 <- c(r1, fr$r1, fa$r1); r2 <- c(r2, fr$r2, fa$r2)
  }
  r1 <- add_err(r1); r2 <- add_err(r2)
  idx <- seq_along(r1)
  write_fq(file.path(dirp, paste0(sample, "_1.fastq.gz")), paste0(sample, "_", idx), r1, 1)
  write_fq(file.path(dirp, paste0(sample, "_2.fastq.gz")), paste0(sample, "_", idx), r2, 2)
}

truth_genes <- data.frame(gene_id = gene_ids, n_snps = as.integer(table(snp$gene_id)[gene_ids]),
                          p_alt = unname(pal[gene_ids]), class = unname(cls[gene_ids]), stringsAsFactors = FALSE)
vcf_hdr <- function(sample = NULL) c("##fileformat=VCFv4.2", paste0("##contig=<ID=", chrom, ",length=", nchar(gseq), ">"),
  if (!is.null(sample)) '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
write_common <- function(d) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  write.table(truth_genes, file.path(d, "truth_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(snp, file.path(d, "truth_snps.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
abs_out <- normalizePath(outdir)

# ---- 4. F1
d <- file.path(outdir, "f1"); write_common(d)
sites <- data.frame(CHROM = snp$chrom, POS = snp$pos, ID = ".", REF = snp$ref, ALT = snp$alt, QUAL = ".", FILTER = ".", INFO = ".")
writeLines(c(vcf_hdr(), paste0("#", paste(names(sites), collapse = "\t"))), file.path(d, "parental_snps.vcf"))
write.table(sites, file.path(d, "parental_snps.vcf"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
hapR <- lapply(setNames(gene_ids, gene_ids), function(g) haplo(g, snp$ref[snp$gene_id == g]))
hapA <- lapply(setNames(gene_ids, gene_ids), function(g) haplo(g, snp$alt[snp$gene_id == g]))
ss <- paste0("f1_s", 1:6)
cond <- rep(c("ctrl", "treat"), each = 3)
for (s in ss) simulate_sample(s, hapR, hapA, as.list(pal), d)
write.csv(data.frame(sample = ss, fastq_1 = file.path(abs_out, "f1", paste0(ss, "_1.fastq.gz")),
  fastq_2 = file.path(abs_out, "f1", paste0(ss, "_2.fastq.gz")), condition = cond, cross_direction = "AxB", individual = ""),
  file.path(d, "samples.csv"), row.names = FALSE, quote = FALSE)

# ---- 5. outbred
d <- file.path(outdir, "outbred"); write_common(d)
inds <- paste0("ind", 1:4); os <- paste0("outbred_s", 1:4); ocond <- c("ctrl", "ctrl", "treat", "treat")
for (k in seq_along(inds)) {
  gt <- sample(c("0/1", "0/0", "1/1"), nrow(snp), replace = TRUE, prob = c(0.6, 0.2, 0.2))
  vc <- data.frame(CHROM = snp$chrom, POS = snp$pos, ID = ".", REF = snp$ref, ALT = snp$alt, QUAL = ".", FILTER = ".", INFO = ".",
                   FORMAT = "GT", S = gt)
  names(vc)[10] <- inds[k]
  for (nm in c("all", "het")) {
    f <- file.path(d, paste0(inds[k], ".", nm, ".vcf"))
    writeLines(c(vcf_hdr(inds[k]), paste0("#", paste(names(vc), collapse = "\t"))), f)
    write.table(if (nm == "het") vc[vc[[10]] == "0/1", ] else vc,
                f, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  }
  # haplotypes: REF-side carries ref at het, hom alleles as genotyped; ALT-side carries alt at het
  hR <- list(); hA <- list(); pu <- list()
  for (g in gene_ids) {
    i <- which(snp$gene_id == g); o <- order(snp$pos[i]); i <- i[o]
    a_ref <- ifelse(gt[i] == "1/1", snp$alt[i], snp$ref[i])
    a_alt <- ifelse(gt[i] == "0/0", snp$ref[i], snp$alt[i])
    hR[[g]] <- haplo(g, a_ref); hA[[g]] <- haplo(g, a_alt)
    pu[[g]] <- if (any(gt[i] == "0/1")) pal[[g]] else 0.5
  }
  simulate_sample(os[k], hR, hA, pu, d)
}
write.csv(data.frame(sample = os, fastq_1 = file.path(abs_out, "outbred", paste0(os, "_1.fastq.gz")),
  fastq_2 = file.path(abs_out, "outbred", paste0(os, "_2.fastq.gz")), condition = ocond, cross_direction = "", individual = inds),
  file.path(d, "samples.csv"), row.names = FALSE, quote = FALSE)
cat("done\n")
