# Stage 3 statistics tests (phASER gene-level haplotype test); hap_gene_test is extracted from the shipped skill text.
# Usage: Rscript test_ase_stats_stage3.R <ase-pipeline.md>     (needs >= 8 CPUs: parallel::mclapply, mc.cores = 8)
args <- commandArgs(trailingOnly = TRUE); skill <- args[1]
txt <- readLines(skill)
b <- grep("^# --- ase-stats-begin", txt); e <- grep("^# --- ase-stats-end", txt)
stopifnot(length(b) == 1, length(e) == 1, e > b)
eval(parse(text = txt[(b + 1):(e - 1)]))
fail <- 0
ok <- function(cond, msg) { if (!isTRUE(cond)) { cat("FAIL:", msg, "\n"); fail <<- 1 } else cat("ok  ", msg, "\n") }
if (!exists("hap_gene_test")) { cat("FAIL: missing from the stats block: hap_gene_test\n"); quit(status = 1) }
const <- function(name) {   # a Step 8 default, read from the skill text (never hard-coded here)
  l <- grep(paste0("^\\| `", name, "` \\|"), txt, value = TRUE); stopifnot(length(l) == 1)
  as.numeric(trimws(strsplit(l, "\\|")[[1]][3]))
}
RHO_MIN <- const("RHO_MIN"); MIN_DEPTH <- const("MIN_DEPTH")
ok(RHO_MIN == 0.01 && MIN_DEPTH == 10, sprintf("Step 8 defaults read from the skill: RHO_MIN %g, MIN_DEPTH %g", RHO_MIN, MIN_DEPTH))
MC <- 8
seeds_rbind <- function(seeds, f) {
  r <- parallel::mclapply(seeds, f, mc.cores = MC)
  stopifnot(!any(vapply(r, function(x) inherits(x, "try-error"), logical(1))))
  do.call(rbind, r)
}
gate <- function(m, col, label, pooled = 0.065, worst = 0.09)
  ok(mean(m[, col]) <= pooled && max(m[, col]) <= worst,
     sprintf("%s: pooled %.4f (<= %.3f), worst seed %.3f (<= %.2f)", label, mean(m[, col]), pooled, max(m[, col]), worst))
report <- function(m, label) cat(sprintf("     %s: %s\n", label, paste(sprintf("%s %.4f", colnames(m), colMeans(m, na.rm = TRUE)), collapse = ", ")))
depth_ln <- function(m, med = 200, sdlog = 1) round(exp(rnorm(m, log(med), sdlog))) + 5
draw <- function(n, h, rho) { pr <- if (rho > 0) rbeta(length(n), h * (1 - rho) / rho, (1 - h) * (1 - rho) / rho) else h; rbinom(length(n), n, pr) }

# ---- 1. interface, symmetry (haplotype labels are arbitrary), unusable rows, few genes, the RHO_MIN floor
set.seed(1); n <- depth_ln(400); a <- draw(n, 0.5, 0.02); smp <- rep(c("s1", "s2"), each = 200)
r <- hap_gene_test(a, n - a, smp, RHO_MIN)
ok(identical(names(r), c("p", "rho_used", "samples")) &&
   identical(names(r$samples), c("sample", "n_genes", "rho_own", "central_frac", "rho_cohort", "rho_used")),
   "hap_gene_test returns p, rho_used and samples (sample, n_genes, rho_own, central_frac, rho_cohort, rho_used)")
ok(length(r$p) == 400 && all(is.finite(r$p)) && all(r$p >= 0 & r$p <= 1), "one p per row, every p in [0, 1]")
r_sw <- hap_gene_test(n - a, a, smp, RHO_MIN)
ok(isTRUE(all.equal(r$p, r_sw$p, tolerance = 1e-12)) && isTRUE(all.equal(r$samples, r_sw$samples, tolerance = 1e-12)),
   "swapping the haplotype labels (a <-> b) in every row leaves p and the dispersions unchanged")
fl <- runif(400) < 0.5
r_fl <- hap_gene_test(ifelse(fl, n - a, a), ifelse(fl, a, n - a), smp, RHO_MIN)
ok(isTRUE(all.equal(r$p, r_fl$p, tolerance = 1e-12)) && isTRUE(all.equal(r$samples, r_fl$samples, tolerance = 1e-12)),
   "flipping the labels of a random half of the rows leaves p and the dispersions unchanged")
ae <- c(a, 0, 5, NA, 7, Inf, 30); be <- c(n - a, 0, 0, 3, NA, 2, 0); se <- c(smp, rep("s1", 6))
re <- hap_gene_test(ae, be, se, RHO_MIN)
ok(all(is.na(re$p[400 + c(1, 3, 4, 5)])), "zero-coverage, NA and Inf rows: p NA")
ok(is.finite(re$p[402]) && re$p[406] < 1e-3, "monoallelic rows are tested: 5 of 5 gives a finite p, 30 of 30 gives p < 1e-3")
keep <- c(1:400, 402, 406)
ok(isTRUE(all.equal(re$p[keep], hap_gene_test(ae[keep], be[keep], se[keep], RHO_MIN)$p, tolerance = 1e-12)),
   "NA, Inf and zero rows do not enter the dispersion fit (same p as without them)")
ok(identical(as.integer(re$samples$n_genes), c(202L, 200L)), "n_genes counts the usable rows per sample")
r19 <- hap_gene_test(a[1:19], (n - a)[1:19], rep("s1", 19), RHO_MIN)
ok(all(is.na(r19$p)) && is.na(r19$samples$rho_cohort), "no sample with 20 usable genes: nothing tested (p NA, rho_cohort NA)")
ix <- c(1:15, 201:400)
mix <- hap_gene_test(a[ix], (n - a)[ix], rep(c("s1", "s2"), c(15, 200)), RHO_MIN)
ok(is.na(mix$samples$rho_own[1]) && all(is.finite(mix$p[1:15])) &&
   isTRUE(all.equal(mix$samples$rho_used[1], max(mix$samples$rho_cohort[1], RHO_MIN))),
   "a sample with fewer than 20 genes is tested with max(rho_cohort, RHO_MIN)")
floor_run <- function(s) { set.seed(100 + s); nn <- depth_ln(12000); aa <- draw(nn, 0.5, 0)
  hap_gene_test(aa, nn - aa, rep(paste0("s", 1:6), each = 2000), RHO_MIN)$samples$rho_used }
fl0 <- unlist(parallel::mclapply(1:10, floor_run, mc.cores = MC))
ok(length(fl0) == 60 && all(fl0 == RHO_MIN), sprintf("binomial data: rho_used sits at the RHO_MIN floor in all 60 samples (max %.4f)", max(fl0)))

# ---- 2. cohorts: null size with real imbalance present, power against an oracle, few genes, heterogeneous dispersion
cohort_run <- function(seed, rho, frac_imb, S = 6, G = 2000, h_imb = 0.75, rho_s = NULL, G_s = NULL) {
  set.seed(seed)
  Gs <- if (is.null(G_s)) rep(G, S) else G_s; rs <- if (is.null(rho_s)) rep(rho, S) else rho_s; S <- length(Gs)
  d <- do.call(rbind, lapply(seq_len(S), function(s) {
    nn <- depth_ln(Gs[s]); imb <- runif(Gs[s]) < frac_imb
    h <- ifelse(imb, ifelse(runif(Gs[s]) < 0.5, h_imb, 1 - h_imb), 0.5)
    data.frame(sample = paste0("s", s), a = draw(nn, h, rs[s]), n = nn, imb = imb, rho_true = rs[s], stringsAsFactors = FALSE)
  }))
  r <- hap_gene_test(d$a, d$n - d$a, d$sample, RHO_MIN)
  po <- vapply(seq_len(nrow(d)), function(i) bb_pvalue(d$a[i], d$n[i], d$rho_true[i]), numeric(1))   # oracle: the true dispersion
  sz <- function(p, sel) mean(p[sel] < 0.05, na.rm = TRUE)
  c(size = sz(r$p, !d$imb), size_s1 = sz(r$p, !d$imb & d$sample == "s1"),
    power = if (frac_imb > 0) sz(r$p, d$imb) else NA_real_, oracle_power = if (frac_imb > 0) sz(po, d$imb) else NA_real_,
    oracle_size = sz(po, !d$imb), rho_cohort = r$samples$rho_cohort[1], rho_used_mean = mean(r$samples$rho_used))
}
k <- 0
for (rho in c(0.005, 0.02, 0.05)) for (fi in c(0, 0.2, 0.4)) {
  k <- k + 1
  m <- seeds_rbind(1:10, function(s) cohort_run(1000 * k + s, rho, fi))
  lab <- sprintf("cohort 6 x 2000 genes, true rho %.3f, %d%% imbalanced", rho, round(100 * fi))
  report(m, lab); gate(m, "size", paste("null size,", lab))
  if (rho == 0.02 && fi == 0.2)
    ok(mean(m[, "power"]) >= 0.7 && mean(m[, "power"]) >= 0.9 * mean(m[, "oracle_power"]),
       sprintf("power, h 0.75, rho 0.02, 20%% imbalanced: %.4f (>= 0.7 and >= 0.9 x oracle %.4f)", mean(m[, "power"]), mean(m[, "oracle_power"])))
}
avg_runs <- function(seed0, reps, ...) { x <- t(vapply(seq_len(reps), function(j) cohort_run(seed0 + j, ...), numeric(7))); colMeans(x, na.rm = TRUE) }
for (G in c(20, 30, 60)) {
  m <- seeds_rbind(1:10, function(s) avg_runs(20000 + 1000 * G + 100 * s, 30, rho = 0.02, frac_imb = 0.2, G = G))
  report(m, sprintf("few genes: 30 cohorts of 6 samples x %d genes, rho 0.02, 20%% imbalanced", G))
  gate(m, "size", sprintf("null size with %d genes per sample", G))
}
m <- seeds_rbind(1:10, function(s) avg_runs(300000 + 1000 * s, 100, rho = 0.02, frac_imb = 0.2, G_s = c(15, rep(300, 5))))
report(m, "one sample with 15 genes in a cohort of 6 (tested with max(rho_cohort, RHO_MIN))")
gate(m, "size_s1", "null size of the 15-gene sample")
m <- seeds_rbind(1:10, function(s) cohort_run(40000 + s, 0.02, 0.2, rho_s = c(0.08, 0.05, 0.02, 0.02, 0.01, 0.01)))
report(m, "heterogeneous cohort (true rho 0.08, 0.05, 0.02, 0.02, 0.01, 0.01), 20% imbalanced")
gate(m, "size", "heterogeneous cohort, null size over all samples"); gate(m, "size_s1", "heterogeneous cohort, null size of the rho 0.08 sample")

# ---- 3. fragment level: phASER-style haplotype counts (each read pair once) against the unphased per-SNP path + ACAT (Rmd 02)
frag_gene <- function(pos, h, nf, span, rl = 100) {   # one gene and sample: fragments, haplotype-1 indicator, SNPs each fragment covers
  L <- round(runif(nf, 150, 350)); st <- floor(runif(nf) * (span - L + 1)) + 1; A <- runif(nf) < h
  cv <- (outer(st, pos, "<=") & outer(st + rl - 1, pos, ">=")) | (outer(st + L - rl, pos, "<=") & outer(st + L - 1, pos, ">="))
  anyc <- rowSums(cv) > 0
  list(snp_a = colSums(cv & A), snp_n = colSums(cv), hap_a = sum(anyc & A), hap_n = sum(anyc))
}
frag_run <- function(seed, lambda, S = 6, G = 1500, rho_bio = 0.02, frac_imb = 0.3, span = 1500) {
  set.seed(seed)
  k <- sample(c(1, 2, 4, 8), G, TRUE); imb <- runif(G) < frac_imb
  pos <- lapply(k, function(kk) sort(sample(400:1000, kk)))             # the SNPs of a gene lie within 600 transcript bases
  out <- lapply(seq_len(S), function(s) {
    h <- ifelse(imb, ifelse(runif(G) < 0.5, 0.75, 0.25), 0.5)
    hs <- rbeta(G, h * (1 - rho_bio) / rho_bio, (1 - h) * (1 - rho_bio) / rho_bio)
    lapply(seq_len(G), function(g) frag_gene(pos[[g]], hs[g], rpois(1, lambda), span))
  })
  acat_p <- matrix(NA_real_, G, S); hap_a <- matrix(NA_real_, G, S); hap_n <- matrix(NA_real_, G, S)
  for (s in seq_len(S)) {   # unphased path of Rmd 01/02: sites with depth >= MIN_DEPTH, trimmed rho floored at RHO_MIN, bb_pvalue, ACAT
    sa <- unlist(lapply(out[[s]], `[[`, "snp_a")); sn <- unlist(lapply(out[[s]], `[[`, "snp_n")); sg <- rep(seq_len(G), k)
    use <- sn >= MIN_DEPTH
    rt <- bb_estimate_rho_trim(sa[use], sn[use])[["rho"]]; rh <- if (is.na(rt)) NA_real_ else max(rt, RHO_MIN)
    ps <- rep(NA_real_, length(sa))
    if (!is.na(rh)) ps[use] <- vapply(which(use), function(i) bb_pvalue(sa[i], sn[i], rh), numeric(1))
    acat_p[, s] <- vapply(seq_len(G), function(g) acat(ps[sg == g]), numeric(1))
    hap_a[, s] <- vapply(out[[s]], `[[`, numeric(1), "hap_a"); hap_n[, s] <- vapply(out[[s]], `[[`, numeric(1), "hap_n")
  }
  okn <- hap_n >= MIN_DEPTH
  r <- hap_gene_test(ifelse(okn, hap_a, NA), ifelse(okn, hap_n - hap_a, NA), rep(paste0("s", seq_len(S)), each = G), RHO_MIN)
  hp <- matrix(r$p, G, S); I <- matrix(imb, G, S); big <- matrix(k >= 4, G, S); K8 <- matrix(k == 8, G, S)
  det <- function(p, sel) mean(!is.na(p[sel]) & p[sel] < 0.05)          # not tested counts as not detected
  c(hap_size = mean(hp[!I] < 0.05, na.rm = TRUE), acat_size = mean(acat_p[!I] < 0.05, na.rm = TRUE),
    hap_power_k4 = det(hp, I & big), acat_power_k4 = det(acat_p, I & big),
    hap_tested_k8 = mean(!is.na(hp[K8])), acat_tested_k8 = mean(!is.na(acat_p[K8])))
}
for (lam in c(300, 40)) {
  m <- seeds_rbind(1:10, function(s) frag_run(50000 + 100 * lam + s, lam)); report(m, sprintf("fragment level, %d fragments per gene", lam))
  gate(m, "hap_size", sprintf("fragment level, lambda %d, haplotype test null size", lam))
  ok(mean(m[, "hap_power_k4"]) >= mean(m[, "acat_power_k4"]),
     sprintf("fragment level, lambda %d, genes with >= 4 SNPs: haplotype-test power %.4f >= unphased ACAT power %.4f", lam,
             mean(m[, "hap_power_k4"]), mean(m[, "acat_power_k4"])))
  if (lam == 40) ok(mean(m[, "hap_tested_k8"]) >= 0.9,
     sprintf("lambda 40, 8-SNP genes: tested by the haplotype test %.4f (>= 0.9); by the unphased path %.4f (reported)",
             mean(m[, "hap_tested_k8"]), mean(m[, "acat_tested_k8"])))
}

# ---- 4. runtime and memory: one sample of genome scale on one core (Rmd 05 runs the samples one after the other)
set.seed(900); nn <- pmin(round(exp(rnorm(20000, log(300), 1.2))) + 5, 50000); aa <- draw(nn, 0.5, 0.02)
invisible(gc(reset = TRUE))
tt <- system.time(rr <- hap_gene_test(aa, nn - aa, rep("s1", 20000), RHO_MIN))[["elapsed"]]
g <- gc(); mx <- sum(g[, ncol(g)])
ok(tt <= 120 && mx <= 4096, sprintf("runtime: one sample, 20,000 genes (depth median 300, max 50,000): %.0f s (<= 120), max memory %.0f MB (<= 4096); projected 100 samples %.2f h on one core",
                                    tt, mx, tt * 100 / 3600))
quit(status = fail)
