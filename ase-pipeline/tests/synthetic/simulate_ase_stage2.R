#!/usr/bin/env Rscript
# Stage 2 synthetic ground-truth ASE data: reciprocal x two-condition F1 (f1_recip/) and paired two-condition outbred (outbred_diff/).
# Usage: Rscript simulate_ase_stage2.R <outdir> <seed>
# Direct counts (ASEReadCounter format) are counted from the simulated fragments, ignoring sequencing errors.
suppressPackageStartupMessages(library(Biostrings))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: simulate_ase_stage2.R <outdir> <seed>")
outdir <- args[1]; seed <- as.integer(args[2]); set.seed(seed)
READLEN <- 100L; ERR <- 0.002; LAMBDA <- 450; PHI_BIO <- 0.005; NGENES <- 60L; GLEN <- 1000000L
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
writeLines(as.character(seed), file.path(outdir, "seed.txt"))
abs_out <- normalizePath(outdir)

# ---- genome and genes (Stage 1 synthetic block, 60 genes of 4-5 exons)
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

# ---- classes
ord <- sample(gene_ids)
bias_gene <- ord[1]; dense_genes <- ord[2:3]
plan <- data.frame(
  class = c("strain_A_high", "strain_A_high", "strain_B_high", "strain_B_high", "maternal", "maternal", "paternal", "paternal",
            "strain_and_maternal", "diff_up", "diff_up", "diff_down", "diff_on_strain"),
  fA = c(0.75, 0.70, 0.25, 0.30, 0.5, 0.5, 0.5, 0.5, 0.70, 0.5, 0.5, 0.5, 0.70),
  fM = c(0.5, 0.5, 0.5, 0.5, 0.85, 0.95, 0.15, 0.10, 0.80, 0.5, 0.5, 0.5, 0.5),
  fC = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.75, 0.75, 0.25, 0.30),
  expect_strain = c("A_higher", "A_higher", "B_higher", "B_higher", "none", "none", "none", "none", "A_higher", "any", "any", "any", "any"),
  expect_poo = c("none", "none", "none", "none", "maternal", "maternal", "paternal", "paternal", "maternal", "none", "none", "none", "none"),
  expect_diff = c(rep("none", 9), "up", "up", "down", "down"), stringsAsFactors = FALSE)
f1 <- data.frame(gene_id = gene_ids, class = "null", b0 = 0, b1 = 0, bc = 0, expect_strain = "none", expect_poo = "none",
                 expect_diff = "none", stringsAsFactors = FALSE)
rownames(f1) <- gene_ids
f1[bias_gene, "class"] <- "bias"; f1[dense_genes, "class"] <- "dense_null"
pg <- ord[3 + seq_len(nrow(plan))]
f1[pg, "class"] <- plan$class; f1[pg, "b0"] <- qlogis(plan$fA); f1[pg, "b1"] <- qlogis(plan$fM); f1[pg, "bc"] <- qlogis(plan$fC)
f1[pg, "expect_strain"] <- plan$expect_strain; f1[pg, "expect_poo"] <- plan$expect_poo; f1[pg, "expect_diff"] <- plan$expect_diff

# ---- SNPs (transcript offsets == spliced offsets: plus-strand genes)
snps <- list()
for (g in gene_ids) {
  x <- genes[[g]]; len <- x$len
  if (g == bias_gene) {
    w0 <- sample.int(len - 59, 1); off <- sort(w0 + sample.int(60, 5) - 1)
  } else if (g %in% dense_genes) {
    repeat { w0 <- sample.int(len - 299, 1); off <- sort(w0 + sample.int(300, 8) - 1); if (all(diff(off) >= 20)) break }
  } else {
    repeat { off <- sort(sample.int(len, 4)); if (all(diff(off) >= 40)) break }
  }
  cum <- c(0, cumsum(x$end - x$start + 1))
  pos <- vapply(off, function(o) { i <- max(which(cum < o)); x$start[i] + o - cum[i] - 1 }, numeric(1))
  ref <- substring(gseq, pos, pos)
  alt <- vapply(ref, alt_of, character(1))
  snps[[g]] <- data.frame(chrom = chrom, pos = pos, ref = unname(ref), alt = unname(alt), gene_id = g, tx_off = off, stringsAsFactors = FALSE)
}
snp <- do.call(rbind, snps); rownames(snp) <- NULL
stopifnot(!anyDuplicated(snp$pos))
snp <- snp[order(snp$pos), ]; rownames(snp) <- NULL
gidx <- lapply(setNames(gene_ids, gene_ids), function(g) which(snp$gene_id == g))   # ascending pos == ascending offset

haplo <- function(g, alleles) { s <- strsplit(spliced(g), "")[[1]]; s[snp$tx_off[gidx[[g]]]] <- alleles; paste(s, collapse = "") }

frags2 <- function(hseq, n, snp_offs) {
  L <- nchar(hseq)
  if (n == 0) return(list(r1 = character(0), r2 = character(0), cover = matrix(FALSE, 0, length(snp_offs))))
  fl <- pmin(pmax(round(rnorm(n, 250, 30)), 150), 350, L); st <- floor(runif(n) * (L - fl + 1)) + 1
  r1 <- substring(hseq, st, st + READLEN - 1); r2 <- revcomp(substring(hseq, st + fl - READLEN, st + fl - 1))
  cover <- (outer(st, snp_offs, "<=") & outer(st + READLEN - 1, snp_offs, ">=")) |
           (outer(st + fl - READLEN, snp_offs, "<=") & outer(st + fl - 1, snp_offs, ">="))
  list(r1 = r1, r2 = r2, cover = cover)
}

# simulate one sample. hR/hA: per-gene REF-side / ALT-side haplotypes; p_ref: per-gene probability of a REF-side fragment
# (before biological overdispersion); counted: per-gene logical over the gene's SNPs (which SNPs get direct counts).
# Returns list(truth = per-gene fragment counts, counts = ASEReadCounter table).
sim_sample <- function(sample, dirp, hR, hA, p_ref, counted) {
  r1 <- character(0); r2 <- character(0); refc <- numeric(nrow(snp)); altc <- numeric(nrow(snp)); tr <- list()
  for (g in gene_ids) {
    p <- p_ref[[g]]
    ps <- if (p <= 0 || p >= 1) p else rbeta(1, p * (1 - PHI_BIO) / PHI_BIO, (1 - p) * (1 - PHI_BIO) / PHI_BIO)
    n <- rpois(1, LAMBDA); nA <- rbinom(1, n, ps)
    offs <- snp$tx_off[gidx[[g]]]
    fr <- frags2(hR[[g]], nA, offs); fa <- frags2(hA[[g]], n - nA, offs)
    r1 <- c(r1, fr$r1, fa$r1); r2 <- c(r2, fr$r2, fa$r2)
    refc[gidx[[g]]] <- colSums(fr$cover) * counted[[g]]; altc[gidx[[g]]] <- colSums(fa$cover) * counted[[g]]
    tr[[g]] <- data.frame(sample = sample, gene_id = g, frag_A = nA, frag_B = n - nA, stringsAsFactors = FALSE)
  }
  r1 <- add_err(r1); r2 <- add_err(r2); idx <- seq_along(r1)
  write_fq(file.path(dirp, paste0(sample, "_1.fastq.gz")), paste0(sample, "_", idx), r1, 1)
  write_fq(file.path(dirp, paste0(sample, "_2.fastq.gz")), paste0(sample, "_", idx), r2, 2)
  tot <- refc + altc; k <- tot >= 1
  counts <- data.frame(contig = snp$chrom, position = snp$pos, variantID = ".", refAllele = snp$ref, altAllele = snp$alt,
                       refCount = refc, altCount = altc, totalCount = tot, lowMAPQDepth = 0, lowBaseQDepth = 0, rawDepth = tot,
                       otherBases = 0, improperPairs = 0, stringsAsFactors = FALSE)[k, ]
  list(truth = do.call(rbind, tr), counts = counts)
}
wtab <- function(x, f) write.table(x, f, sep = "\t", quote = FALSE, row.names = FALSE)

# ================= F1 =================
d <- file.path(outdir, "f1_recip"); dir.create(file.path(d, "direct_counts"), recursive = TRUE, showWarnings = FALSE)
cellp <- function(j, dd, t) plogis(f1$b0[j] + f1$b1[j] * dd + f1$bc[j] * t)
tg <- data.frame(gene_id = gene_ids, n_snps = as.integer(table(snp$gene_id)[gene_ids]), class = f1$class,
                 b0 = format(f1$b0, digits = 12, trim = TRUE), b1 = format(f1$b1, digits = 12, trim = TRUE), bc = format(f1$bc, digits = 12, trim = TRUE),
                 p_A_AxB_ctrl = format(cellp(seq_along(gene_ids), 1, 0), digits = 12, trim = TRUE),
                 p_A_AxB_treat = format(cellp(seq_along(gene_ids), 1, 1), digits = 12, trim = TRUE),
                 p_A_BxA_ctrl = format(cellp(seq_along(gene_ids), -1, 0), digits = 12, trim = TRUE),
                 p_A_BxA_treat = format(cellp(seq_along(gene_ids), -1, 1), digits = 12, trim = TRUE),
                 expect_strain = f1$expect_strain, expect_poo = f1$expect_poo, expect_diff = f1$expect_diff, stringsAsFactors = FALSE)
wtab(tg, file.path(d, "truth_genes.tsv")); wtab(snp, file.path(d, "truth_snps.tsv"))
sites <- data.frame(CHROM = snp$chrom, POS = snp$pos, ID = ".", REF = snp$ref, ALT = snp$alt, QUAL = ".", FILTER = ".", INFO = ".")
writeLines(c(vcf_hdr(), paste0("#", paste(names(sites), collapse = "\t"))), file.path(d, "parental_snps.vcf"))
write.table(sites, file.path(d, "parental_snps.vcf"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
hR <- lapply(setNames(gene_ids, gene_ids), function(g) haplo(g, snp$ref[gidx[[g]]]))
hA <- lapply(setNames(gene_ids, gene_ids), function(g) haplo(g, snp$alt[gidx[[g]]]))
ss <- sprintf("f1r_s%02d", 1:12)
dirn <- rep(c("AxB", "BxA"), each = 6); cond <- rep(rep(c("ctrl", "treat"), each = 3), 2)
allc <- setNames(lapply(gene_ids, function(g) rep(TRUE, length(gidx[[g]]))), gene_ids)
tfa <- list()
for (k in seq_along(ss)) {
  dd <- if (dirn[k] == "AxB") 1 else -1; t <- as.numeric(cond[k] == "treat")
  p <- setNames(as.list(plogis(f1$b0 + f1$b1 * dd + f1$bc * t)), gene_ids)
  r <- sim_sample(ss[k], d, hR, hA, p, allc)
  tfa[[k]] <- r$truth; wtab(r$counts, file.path(d, "direct_counts", paste0(ss[k], ".table")))
}
wtab(do.call(rbind, tfa), file.path(d, "truth_fragments.tsv"))
write.csv(data.frame(sample = ss, fastq_1 = file.path(abs_out, "f1_recip", paste0(ss, "_1.fastq.gz")),
  fastq_2 = file.path(abs_out, "f1_recip", paste0(ss, "_2.fastq.gz")), condition = cond, cross_direction = dirn, individual = ss),
  file.path(d, "samples.csv"), row.names = FALSE, quote = FALSE)

# ================= outbred, paired two-condition =================
d <- file.path(outdir, "outbred_diff"); dir.create(file.path(d, "direct_counts"), recursive = TRUE, showWarnings = FALSE)
inds <- paste0("ind", 1:4)
ord2 <- sample(setdiff(gene_ids, bias_gene))
oc <- data.frame(gene_id = gene_ids, class = "null", f_ctrl = 0.5, f_treat = 0.5, consistent = FALSE, expect_diff = "none", stringsAsFactors = FALSE)
rownames(oc) <- gene_ids
oc[bias_gene, "class"] <- "bias"
bi <- ord2[1:3]; dp <- ord2[4:6]; dc <- ord2[7:8]
oc[bi, c("class", "f_ctrl", "f_treat")] <- list("base_imbalanced", 0.75, 0.75)
oc[dp, c("class", "f_ctrl", "f_treat")] <- list("diff_phase", 0.5, 0.8)
oc[dc, c("class", "f_ctrl", "f_treat")] <- list("diff_consistent", 0.5, 0.75); oc[dc, "consistent"] <- TRUE; oc[dc, "expect_diff"] <- "up"
GT <- function(n) matrix(sample(c("0/1", "0/0", "1/1"), n * 4, replace = TRUE, prob = c(0.6, 0.2, 0.2)), n, 4)
gtm <- matrix(NA_character_, nrow(snp), 4); phase <- matrix(1, length(gene_ids), 4, dimnames = list(gene_ids, inds))
for (g in gene_ids) {
  ph <- if (g %in% dp) c(1, 1, -1, -1) else sample(c(1, -1), 4, replace = TRUE)
  planted <- !(oc[g, "class"] %in% c("null", "bias"))
  repeat {
    m <- GT(length(gidx[[g]])); nh <- colSums(m == "0/1")
    if (!planted) break
    if (sum(nh > 0) >= 3 && (!(g %in% dp) || (any(ph[nh > 0] == 1) && any(ph[nh > 0] == -1)))) break
  }
  gtm[gidx[[g]], ] <- m; phase[g, ] <- ph
}
tg2 <- data.frame(gene_id = gene_ids, n_snps = as.integer(table(snp$gene_id)[gene_ids]), class = oc$class, f_ctrl = oc$f_ctrl, f_treat = oc$f_treat,
                  consistent = ifelse(oc$consistent, "TRUE", "FALSE"), expect_diff = oc$expect_diff, stringsAsFactors = FALSE)
wtab(tg2, file.path(d, "truth_genes.tsv")); wtab(snp, file.path(d, "truth_snps.tsv"))
tig <- list(); hap <- list(); pa <- list()
for (k in 1:4) {
  gt <- gtm[, k]
  vc <- data.frame(CHROM = snp$chrom, POS = snp$pos, ID = ".", REF = snp$ref, ALT = snp$alt, QUAL = ".", FILTER = ".", INFO = ".", FORMAT = "GT", S = gt)
  names(vc)[10] <- inds[k]
  for (nm in c("all", "het")) {
    f <- file.path(d, paste0(inds[k], ".", nm, ".vcf"))
    writeLines(c(vcf_hdr(inds[k]), paste0("#", paste(names(vc), collapse = "\t"))), f)
    write.table(if (nm == "het") vc[vc[[10]] == "0/1", ] else vc, f, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  }
  hR <- list(); hA <- list(); pc <- list(); pt <- list(); cnt <- list()
  for (g in gene_ids) {
    i <- gidx[[g]]
    hR[[g]] <- haplo(g, ifelse(gt[i] == "1/1", snp$alt[i], snp$ref[i]))
    hA[[g]] <- haplo(g, ifelse(gt[i] == "0/0", snp$ref[i], snp$alt[i]))
    nhet <- sum(gt[i] == "0/1"); cnt[[g]] <- gt[i] == "0/1"
    pal <- function(f) if (nhet == 0) 0.5 else if (oc[g, "consistent"] || phase[g, k] == 1) f else 1 - f
    pc[[g]] <- pal(oc[g, "f_ctrl"]); pt[[g]] <- pal(oc[g, "f_treat"])
    tig[[length(tig) + 1]] <- data.frame(individual = inds[k], gene_id = g, phase = phase[g, k], n_het = nhet, p_alt_ctrl = pc[[g]], p_alt_treat = pt[[g]])
  }
  hap[[k]] <- list(hR = hR, hA = hA, pc = pc, pt = pt, cnt = cnt)
}
wtab(do.call(rbind, tig), file.path(d, "truth_individual_genes.tsv"))
os <- paste0("ob_s", 1:8); oind <- rep(inds, 2); ocond <- rep(c("ctrl", "treat"), each = 4)
tfa <- list()
for (j in 1:8) {
  h <- hap[[match(oind[j], inds)]]
  palt <- if (ocond[j] == "ctrl") h$pc else h$pt
  r <- sim_sample(os[j], d, h$hR, h$hA, lapply(palt, function(x) 1 - x), h$cnt)   # p_alt is the ALT-side fraction
  tfa[[j]] <- r$truth
  wtab(r$counts, file.path(d, "direct_counts", paste0(os[j], ".table")))
  file.copy(file.path(d, "direct_counts", paste0(os[j], ".table")), file.path(d, "direct_counts", paste0(os[j], ".unfiltered.table")))
  writeLines(c("vW\tvA\tn", paste("1", "1", sum(r$truth$frag_A), sep = "\t"), paste("1", "2", sum(r$truth$frag_B), sep = "\t"), "none\t-\t0"),
             file.path(d, "direct_counts", paste0(os[j], ".wasp_stats.tsv")))
}
wtab(do.call(rbind, tfa), file.path(d, "truth_fragments.tsv"))
write.csv(data.frame(sample = os, fastq_1 = file.path(abs_out, "outbred_diff", paste0(os, "_1.fastq.gz")),
  fastq_2 = file.path(abs_out, "outbred_diff", paste0(os, "_2.fastq.gz")), condition = ocond, cross_direction = "NA", individual = oind),
  file.path(d, "samples.csv"), row.names = FALSE, quote = FALSE)
cat("done\n")
