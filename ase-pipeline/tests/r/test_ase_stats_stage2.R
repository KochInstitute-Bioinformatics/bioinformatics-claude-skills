# Stage 2 statistics tests (reciprocal F1, differential ASE); the functions are extracted from the shipped skill text.
# Usage: Rscript test_ase_stats_stage2.R <ase-pipeline.md>     (needs >= 8 CPUs: parallel::mclapply, mc.cores = 8)
args <- commandArgs(trailingOnly = TRUE); skill <- args[1]
txt <- readLines(skill)
b <- grep("^# --- ase-stats-begin", txt); e <- grep("^# --- ase-stats-end", txt)
stopifnot(length(b) == 1, length(e) == 1, e > b)
eval(parse(text = txt[(b + 1):(e - 1)]))
fail <- 0
ok <- function(cond, msg) { if (!isTRUE(cond)) { cat("FAIL:", msg, "\n"); fail <<- 1 } else cat("ok  ", msg, "\n") }
fns <- c("chrom_class", "thin_snps", "bb_ll_rows", "bb_glm_fit", "bb_glm_lrt", "bb_moment_phi", "ase_glm_test", "bb_pair_phi_trim", "ase_paired_test")
miss <- fns[!vapply(fns, exists, logical(1))]
if (length(miss) > 0) { cat("FAIL: missing from the stats block:", paste(miss, collapse = ", "), "\n"); quit(status = 1) }
const <- function(name) {   # a Step 8 default, read from the skill text (never hard-coded here)
  l <- grep(paste0("^\\| `", name, "` \\|"), txt, value = TRUE); stopifnot(length(l) == 1)
  as.numeric(trimws(strsplit(l, "\\|")[[1]][3]))
}
RHO_MIN <- const("RHO_MIN"); THIN_BP <- const("THIN_BP")
ok(RHO_MIN == 0.01 && THIN_BP == 500, sprintf("Step 8 defaults read from the skill: RHO_MIN %g, THIN_BP %g", RHO_MIN, THIN_BP))
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

# ---- 1. helpers
ok(identical(chrom_class(c("1", "chr2", "X", "chrX", "Y", "MT", "chrM", "M", "x")),
             c("autosome", "autosome", "X", "X", "Y", "MT", "MT", "MT", "X")), "chrom_class maps Ensembl and chr-prefixed names")
ok(identical(thin_snps(c(100, 150, 700, 720, 1500), c(10, 50, 30, 30, 5), 500), c(FALSE, TRUE, TRUE, FALSE, TRUE)),
   "thin_snps keeps the deepest SNP per window (ties broken by position)")
ok(identical(thin_snps(42, 7, 500), TRUE), "thin_snps keeps a single SNP")
ok(identical(thin_snps(numeric(0), numeric(0), 500), logical(0)), "thin_snps with no SNP returns logical(0)")
em <- tryCatch({ thin_snps(c(100, NA, 900), c(10, 20, 30), 500); "" }, error = function(e) conditionMessage(e))
ok(grepl("thin_snps: 1 SNP(s) with a missing position or depth", em, fixed = TRUE), "thin_snps stops with a clear message on a missing position")
ok(all(thin_snps(c(1, 2000, 4000), c(1, 1, 1), 500)), "SNPs further apart than the window are all kept")
set.seed(11); ps <- sort(sample(1:5000, 60)); kp <- ps[thin_snps(ps, rpois(60, 50), THIN_BP)]
ok(all(diff(kp) >= THIN_BP), "kept SNPs are at least THIN_BP apart")

# ---- 2. GLM fit: equals glm(binomial) at rho = 0 and aod::betabin with the dispersion fixed
set.seed(12); n <- rpois(12, 80) + 5; d <- rep(c(1, -1), 6); X <- cbind("(Intercept)" = 1, d = d)
y <- rbinom(12, n, plogis(0.4 + 0.8 * d))
g <- glm(cbind(y, n - y) ~ d, family = binomial); f0 <- bb_glm_fit(y, n, X, 0)
ok(max(abs(f0$beta - coef(g))) < 1e-4 && abs(f0$loglik - as.numeric(logLik(g))) < 1e-6,
   "bb_glm_fit with rho = 0 equals glm(binomial): coefficients and log-likelihood")
suppressPackageStartupMessages(library(aod))
dd <- data.frame(y = y, n = n, d = d)
a1 <- aod::betabin(cbind(y, n - y) ~ d, ~ 1, data = dd, fixpar = list(3, 0.02))
a0 <- aod::betabin(cbind(y, n - y) ~ 1, ~ 1, data = dd, fixpar = list(2, 0.02))
f2 <- bb_glm_fit(y, n, X, 0.02); f20 <- bb_glm_fit(y, n, X[, 1, drop = FALSE], 0.02)
ok(max(abs(f2$beta - a1@fixed.param)) < 1e-3, "bb_glm_fit at rho 0.02 matches aod::betabin (phi fixed at 0.02): coefficients")
ok(abs((f2$loglik - f20$loglik) - (a1@logL - a0@logL)) < 1e-3, "bb_glm_fit log-likelihood difference matches aod (same LRT)")
l2 <- bb_glm_lrt(y, n, X, "d", 0.02)
ok(abs(l2$stat - 2 * (f2$loglik - f20$loglik)) < 1e-6 && l2$df == 1, "bb_glm_lrt statistic = 2 x log-likelihood difference, df 1")
el <- tryCatch({ bb_glm_lrt(y, n, X, "cond", 0.02); "" }, error = function(e) conditionMessage(e))
ok(grepl("bb_glm_lrt: tested column(s) not in the design: cond", el, fixed = TRUE), "bb_glm_lrt stops (no silent p = 1) when the tested column is not in the design")
ra <- ase_glm_test(list(u = list(y = y, n = n, X = X)), list(condition = "cond"), 0.01)
ok(ra$status == "fit_error" && is.na(ra$p_condition), "ase_glm_test: a test of a column absent from the unit's design is fit_error with p NA")

# ---- 3. separation and driver edge cases
set.seed(13)
mk <- function(y, n, d) list(y = y, n = n, X = cbind("(Intercept)" = 1, d = d))
d6 <- c(1, 1, 1, -1, -1, -1)
base <- lapply(1:200, function(i) { n <- rpois(6, 100) + 10; mk(rbinom(6, n, plogis(rnorm(1, 0, 0.5))), n, d6) })
names(base) <- paste0("g", 1:200)
edge <- list(sep = mk(c(60, 70, 50, 0, 0, 0), c(60, 70, 50, 55, 65, 45), d6),
             onedir = mk(c(30, 32, 28), rep(50, 3), c(1, 1, 1)),
             onerow = mk(30, 50, 1),
             empty = list(y = numeric(0), n = numeric(0), X = cbind("(Intercept)" = numeric(0), d = numeric(0))),
             broken = mk(c(30, 32, 28, 20, 22, 25), c(50, 50, Inf, 50, 50, 50), d6),
             onesnp = mk(c(40, 45, 38, 12, 10, 15), rep(50, 6), d6))
r <- ase_glm_test(c(base, edge), list(strain = "(Intercept)", parent_of_origin = "d"), RHO_MIN)
ok(identical(names(r), c("unit", "n_rows", "df_resid", "depth", "phi_common", "phi_unit", "phi_used", "status", "beta_(Intercept)",
                         "beta_d", "stat_strain", "stat_parent_of_origin", "p_strain", "p_parent_of_origin")), "ase_glm_test output columns")
st <- setNames(r$status, r$unit); rr <- function(u, col) r[[col]][r$unit == u]
ok(st[["sep"]] == "ok_at_bound" && rr("sep", "p_parent_of_origin") < 1e-6 && rr("sep", "beta_d") > 5,
   "monoallelic maternal gene (complete separation): ok_at_bound, p < 1e-6, b1 large and positive, no error")
ok(all(st[c("onedir", "onerow", "empty")] == "not_estimable") && all(is.na(r$p_strain[r$unit %in% c("onedir", "onerow", "empty")])),
   "one direction only / one row / no rows: not_estimable with p NA")
ok(st[["broken"]] == "fit_error" && is.na(rr("broken", "p_strain")), "a unit whose fit fails is fit_error with p NA")
ok(all(st[names(base)] %in% c("ok", "ok_at_bound")) && all(is.finite(r$p_strain[r$unit %in% names(base)])), "a failing unit does not affect the others")
ok(st[["onesnp"]] == "ok" && rr("onesnp", "beta_d") > 0.8 && rr("onesnp", "p_parent_of_origin") < 0.05,
   "a gene with one SNP (one row per sample) is tested; maternal-rich counts give b1 > 0")
ok(all(r$phi_used[r$status %in% c("ok", "ok_at_bound")] >= RHO_MIN), "phi_used is never below RHO_MIN")
# 3b. convergence rule of bb_glm_fit (regression): optim is masked to record or limit L-BFGS-B, the fit code is unchanged
raw <- NA_integer_; optim <- function(...) { o <- stats::optim(...); raw <<- o$convergence; o }
y52 <- c(53, 55, 57, 58, 60, 65); n52 <- c(105, 108, 100, 106, 115, 126); X52 <- cbind("(Intercept)" = 1, d = d6)
f52 <- bb_glm_fit(y52, n52, X52, 0.01); rm(optim)
nm <- stats::optim(c(0, 0), function(b) -sum(bb_ll_rows(y52, n52, plogis(drop(X52 %*% b)), 0.01)), control = list(reltol = 1e-14, maxit = 5000))
ok(identical(raw, 52L) && f52$conv == 0 && abs(f52$loglik + nm$value) < 1e-6 && max(abs(f52$beta - nm$par)) < 1e-3,
   sprintf("L-BFGS-B code %s at the optimum is accepted (conv 0; log-likelihood equals a Nelder-Mead fit to %.1e)", raw, abs(f52$loglik + nm$value)))
optim <- function(...) { a <- list(...); a$control <- list(maxit = 1); do.call(stats::optim, a) }
unc <- list(y = c(10, 80, 15, 70, 5, 90), n = c(100, 300, 60, 100, 20, 150), X = cbind("(Intercept)" = 1, d = d6))
f1 <- bb_glm_fit(unc$y, unc$n, unc$X, 0.02); r1 <- ase_glm_test(list(u = unc), list(poo = "d"), RHO_MIN); rm(optim)
ok(f1$conv != 0 && r1$status == "not_converged" && is.na(r1$p_poo), "a fit stopped after 1 iteration (maxit = 1) stays not_converged with p NA")

# ---- 4. moment estimator: recovers phi, no Neyman-Scott shrinkage with 6 rows and 2 coefficients
sim_units <- function(U, d, phi, depth = function(m) round(runif(m, 30, 500))) {
  u <- lapply(seq_len(U), function(i) {
    p <- plogis(rnorm(1, 0, 0.5) + rnorm(1, 0, 0.5) * d); n <- depth(length(d))
    pr <- if (phi > 0) rbeta(length(d), p * (1 - phi) / phi, (1 - p) * (1 - phi) / phi) else p
    list(y = rbinom(length(d), n, pr), n = n, X = cbind("(Intercept)" = 1, d = d))
  })
  setNames(u, paste0("u", seq_len(U)))
}
for (ph in c(0.005, 0.02, 0.05)) {
  est <- unlist(parallel::mclapply(1:8, function(s) { set.seed(200 + s); bb_moment_phi(sim_units(1000, d6, ph)) }, mc.cores = MC))
  ok(abs(mean(est) / ph - 1) < 0.15, sprintf("bb_moment_phi recovers phi %.3f: mean %.4f over 8 seeds (within 15%%)", ph, mean(est)))
}
est0 <- unlist(parallel::mclapply(1:8, function(s) { set.seed(300 + s); bb_moment_phi(sim_units(1000, d6, 0)) }, mc.cores = MC))
ok(all(est0 < 2e-3), sprintf("binomial data: bb_moment_phi max %.2e (< 2e-3)", max(est0)))

# ---- 5. reciprocal F1 model (rows = samples, logit p_A = b0 + b1 d): size, leakage, power, sign
sim_design <- function(G, d, cond, b0, b1, bc, phi, phi_sd = 0, depth = function(m) round(runif(m, 60, 400))) {
  lapply(seq_len(G), function(g) {
    ph <- if (phi_sd > 0) min(phi * exp(rnorm(1, -phi_sd^2 / 2, phi_sd)), 0.3) else phi
    p <- plogis(b0[g] + b1[g] * d + bc[g] * cond); n <- depth(length(d))
    pr <- if (ph > 0) rbeta(length(d), p * (1 - ph) / ph, (1 - p) * (1 - ph) / ph) else p
    list(y = rbinom(length(d), n, pr), n = n)
  })
}
recip_run <- function(seed, nA, nB, phi, phi_sd = 0) {
  set.seed(seed); d <- c(rep(1, nA), rep(-1, nB))
  cls <- rep(c("null", "strain_only", "imprinted_only", "power_strain", "power_poo"), c(2000, 400, 400, 200, 200)); G <- length(cls)
  b0 <- ifelse(cls == "strain_only", sample(c(-1, 1), G, TRUE) * qlogis(0.8), ifelse(cls == "power_strain", qlogis(0.7), 0))
  b1 <- ifelse(cls == "imprinted_only", sample(c(-1, 1), G, TRUE) * qlogis(0.9), ifelse(cls == "power_poo", qlogis(0.8), 0))
  s <- sim_design(G, d, 0, b0, b1, rep(0, G), phi, phi_sd)
  u <- setNames(lapply(s, function(z) list(y = z$y, n = z$n, X = cbind("(Intercept)" = 1, d = d))), paste0("g", seq_len(G)))
  r <- ase_glm_test(u, list(strain = "(Intercept)", parent_of_origin = "d"), RHO_MIN)
  sz <- function(p, k) mean(p[cls == k] < 0.05, na.rm = TRUE)
  ss <- which(r$p_strain < 0.05 & cls %in% c("strain_only", "power_strain"))           # which(): p = NA never indexes
  sp <- which(r$p_parent_of_origin < 0.05 & cls %in% c("imprinted_only", "power_poo"))
  c(null_strain = sz(r$p_strain, "null"), null_poo = sz(r$p_parent_of_origin, "null"),
    leak_poo = sz(r$p_parent_of_origin, "strain_only"), leak_strain = sz(r$p_strain, "imprinted_only"),
    pow_strain = sz(r$p_strain, "power_strain"), pow_poo = sz(r$p_parent_of_origin, "power_poo"),
    sign_strain = mean(sign(r[["beta_(Intercept)"]][ss]) == sign(b0[ss])), sign_poo = mean(sign(r$beta_d[sp]) == sign(b1[sp])),
    bad_status = mean(!r$status %in% c("ok", "ok_at_bound")), phi_common = r$phi_common[1])
}
rc <- list(list("3+3, phi 0.005", 3, 3, 0.005, 0), list("3+3, phi 0.02", 3, 3, 0.02, 0), list("2+2, phi 0.02", 2, 2, 0.02, 0),
           list("3+2, phi 0.02", 3, 2, 0.02, 0), list("4+2, phi 0.02", 4, 2, 0.02, 0),
           list("3+3, phi 0.02 heterogeneous (lognormal sd 0.7)", 3, 3, 0.02, 0.7))
for (cs in rc) {
  m <- seeds_rbind(1:10, function(s) recip_run(4000 + s, cs[[2]], cs[[3]], cs[[4]], cs[[5]])); report(m, paste("reciprocal", cs[[1]]))
  het <- cs[[5]] > 0; pl <- if (het) 0.075 else 0.065; wo <- if (het) 0.10 else 0.09
  for (col in c("null_strain", "null_poo", "leak_poo", "leak_strain")) gate(m, col, paste("reciprocal", cs[[1]], col), pl, wo)
  ok(min(m[, "sign_strain"]) >= 0.99 && min(m[, "sign_poo"]) >= 0.99, sprintf("reciprocal %s: sign of b0 and b1 correct in >= 99%% of calls", cs[[1]]))
  ok(max(m[, "bad_status"]) <= 0.01, sprintf("reciprocal %s: at most 1%% of genes not ok (max %.4f)", cs[[1]], max(m[, "bad_status"])))
  if (cs[[1]] == "3+3, phi 0.005") ok(mean(m[, "pow_strain"]) >= 0.9 && mean(m[, "pow_poo"]) >= 0.9,
    sprintf("reciprocal 3+3 power: strain 0.7 %.3f, maternal 0.8 %.3f (>= 0.9 each)", mean(m[, "pow_strain"]), mean(m[, "pow_poo"])))
}

# ---- 6. F1 differential (logit p_A = b0 + bc cond [+ b1 d]): size with strain / imprinted genes present, power, confounding
diff_run <- function(seed, dir_ctrl, dir_treat, phi, with_d = TRUE) {
  set.seed(seed); d <- c(dir_ctrl, dir_treat); cond <- rep(c(0, 1), c(length(dir_ctrl), length(dir_treat)))
  cls <- rep(c("null", "imprinted", "power"), c(2000, 600, 300)); G <- length(cls)
  b0 <- rnorm(G, 0, 0.6)
  b1 <- ifelse(cls == "imprinted", sample(c(-1, 1), G, TRUE) * qlogis(0.9), 0)
  bc <- ifelse(cls == "power", sample(c(-1, 1), G, TRUE) * qlogis(0.7), 0)
  s <- sim_design(G, d, cond, b0, b1, bc, phi)
  X <- if (with_d && length(unique(d)) == 2) cbind("(Intercept)" = 1, cond = cond, d = d) else cbind("(Intercept)" = 1, cond = cond)
  u <- setNames(lapply(s, function(z) list(y = z$y, n = z$n, X = X)), paste0("g", seq_len(G)))
  r <- ase_glm_test(u, list(condition = "cond"), RHO_MIN); sg <- r$p_condition < 0.05; pw <- which(sg & cls == "power")
  c(null = mean(sg[cls == "null"], na.rm = TRUE), imprinted = mean(sg[cls == "imprinted"], na.rm = TRUE),
    power = mean(sg[cls == "power"], na.rm = TRUE), sign = mean(sign(r$beta_cond[pw]) == sign(bc[pw])),
    bad_status = mean(!r$status %in% c("ok", "ok_at_bound")))
}
dc <- list(list("3 vs 3, one direction", c(1, 1, 1), c(1, 1, 1), TRUE), list("2 vs 2, one direction", c(1, 1), c(1, 1), TRUE),
           list("3 vs 2, one direction", c(1, 1, 1), c(1, 1), TRUE),
           list("unbalanced directions, with d term", c(1, 1, -1), c(1, -1, -1), TRUE))
for (cs in dc) {
  m <- seeds_rbind(1:10, function(s) diff_run(5000 + s, cs[[2]], cs[[3]], 0.02, cs[[4]])); report(m, paste("differential F1", cs[[1]]))
  gate(m, "null", paste("differential F1", cs[[1]], "null genes")); gate(m, "imprinted", paste("differential F1", cs[[1]], "imprinted genes (no condition effect)"))
  ok(min(m[, "sign"]) >= 0.99 && max(m[, "bad_status"]) <= 0.01, sprintf("differential F1 %s: sign correct, <= 1%% not ok", cs[[1]]))
  if (cs[[1]] == "3 vs 3, one direction") ok(mean(m[, "power"]) >= 0.7, sprintf("differential F1 3 vs 3 power (0.5 -> 0.7): %.3f (>= 0.7)", mean(m[, "power"])))
}
m <- seeds_rbind(1:10, function(s) diff_run(5000 + s, c(1, 1, -1), c(1, -1, -1), 0.02, FALSE))
cat(sprintf("     EVIDENCE unbalanced directions WITHOUT the d term: imprinted-gene size %.4f (with the d term: gated above)\n", mean(m[, "imprinted"])))
set.seed(5100); s <- sim_design(50, c(1, 1, 1, -1, -1, -1), c(0, 0, 0, 1, 1, 1), rep(0, 50), rep(0, 50), rep(0, 50), 0.02)
uc <- setNames(lapply(s, function(z) list(y = z$y, n = z$n, X = cbind("(Intercept)" = 1, cond = c(0, 0, 0, 1, 1, 1), d = c(1, 1, 1, -1, -1, -1)))), paste0("c", 1:50))
rcf <- ase_glm_test(uc, list(condition = "cond"), RHO_MIN)
ok(all(rcf$status == "not_estimable") && all(is.na(rcf$p_condition)), "condition confounded with cross direction: every gene not_estimable (Rmd 04 stops before this)")

# ---- 7. outbred paired test: per individual LRT, summed over individuals (direction-free), ACAT per gene
sim_pairs <- function(S, I, phi_pair, frac_diff = 0, consistent = FALSE, shift = qlogis(0.8), depth = function(m) round(runif(m, 30, 200))) {
  g <- expand.grid(i = seq_len(I), s = seq_len(S)); g <- g[runif(nrow(g)) < 0.6, ]            # heterozygous individual x SNP cells
  diff_snp <- runif(S) < frac_diff; m <- nrow(g); ph <- sample(c(-1, 1), m, replace = TRUE)   # phase differs between individuals
  l0 <- ifelse(runif(m) < 0.5, ph * qlogis(0.7), 0)                                           # baseline imbalance on the high haplotype
  l1 <- l0 + ifelse(diff_snp[g$s], if (consistent) shift else ph * shift, 0)
  n <- matrix(depth(2 * m), m); p <- plogis(cbind(l0, l1))
  pr <- if (phi_pair > 0) matrix(rbeta(2 * m, p * (1 - phi_pair) / phi_pair, (1 - p) * (1 - phi_pair) / phi_pair), m) else p
  data.frame(snp = paste0("s", rep(g$s, 2)), individual = paste0("i", rep(g$i, 2)), cond = rep(c(0, 1), each = m),
             y = rbinom(2 * m, as.vector(n), as.vector(pr)), n = as.vector(n), diff = rep(diff_snp[g$s], 2), stringsAsFactors = FALSE)
}
pair_run <- function(seed, I, phi, frac_diff = 0, consistent = FALSE) {
  set.seed(seed); s <- sim_pairs(3000, I, phi, frac_diff, consistent)
  gmap <- rep(seq_len(3000), sample(1:4, 3000, TRUE))[1:3000]                                  # genes of 1-4 SNPs
  r <- ase_paired_test(s, RHO_MIN); x <- r$snp; idx <- as.integer(sub("^s", "", x$snp))
  truth <- tapply(s$diff, s$snp, any)[x$snp]; gene <- gmap[idx]
  gp <- tapply(x$p, gene, acat); gnull <- tapply(!truth, gene, all)[names(gp)]
  # sizes and power are over TESTED SNPs only (status ok); untested SNPs (too_few_individuals, p NA) are the separate `tested` fraction
  tst <- x$status == "ok"; sg <- tst & x$p < 0.05
  alt_power <- function(fixed) {   # the same test with bb_pair_phi_trim masked: "true" = the true phi (oracle), "all" = its untrimmed phi_all
    real <- bb_pair_phi_trim; on.exit(assign("bb_pair_phi_trim", real, envir = globalenv()))
    mask <- if (fixed == "true") function(y, n, cell, ...) c(phi = phi, phi_all = phi, central_frac = 1) else
      function(y, n, cell, ...) { e <- real(y, n, cell, ...); c(phi = e[["phi_all"]], phi_all = e[["phi_all"]], central_frac = 1) }
    assign("bb_pair_phi_trim", mask, envir = globalenv())
    ro <- ase_paired_test(s, RHO_MIN); xo <- ro$snp; to <- xo$status == "ok"; tr <- tapply(s$diff, s$snp, any)[xo$snp]
    c(mean((xo$p < 0.05)[to & tr]), mean(ro$phi$phi_pair))
  }
  po <- pm <- c(NA_real_, NA_real_)
  if (frac_diff > 0) { po <- alt_power("true"); pm <- alt_power("all") }
  c(snp_size = mean(sg[tst & !truth]), gene_size = mean(gp[gnull] < 0.05, na.rm = TRUE),
    power = if (frac_diff > 0) mean(sg[tst & truth]) else NA_real_, power_oracle = po[1], power_untrimmed = pm[1], phi_untrimmed = pm[2],
    sign_up = if (consistent) mean(x$mean_delta[sg & truth] > 0) else NA_real_,
    tested = mean(tst), phi = mean(r$phi$phi_pair))
}
pow_gate <- function(m, label)   # ruled gate: power >= 0.70 absolute AND >= 0.9 x the oracle power (true dispersion)
  ok(mean(m[, "power"]) >= 0.70 && mean(m[, "power"]) >= 0.9 * mean(m[, "power_oracle"]),
     sprintf("%s: power %.3f (>= 0.70 and >= 0.9 x oracle power %.3f = %.3f)", label, mean(m[, "power"]), mean(m[, "power_oracle"]), 0.9 * mean(m[, "power_oracle"])))
for (ph in c(0.01, 0.02, 0.05)) {   # raw pair dispersion (no floor), no SNP changes: recovers the truth
  e <- seeds_rbind(1:8, function(k) { set.seed(6500 + k); s <- sim_pairs(3000, 1, ph); bb_pair_phi_trim(s$y, s$n, s$snp) })
  ok(abs(mean(e[, "phi"]) / ph - 1) < 0.15, sprintf("bb_pair_phi_trim recovers pair phi %.3f: mean %.4f (all-SNP moment %.4f, kept %.3f) over 8 seeds (within 15%%)",
                                                    ph, mean(e[, "phi"]), mean(e[, "phi_all"]), mean(e[, "central_frac"])))
}
for (I in c(2, 3, 4, 6)) for (ph in c(0.005, 0.02)) {
  m <- seeds_rbind(1:10, function(s) pair_run(6000 + 10 * I + s, I, ph)); report(m, sprintf("paired outbred I = %d, phi %.3f", I, ph))
  cat(sprintf("     TESTED FRACTION I = %d, phi %.3f: %.4f of the reported SNPs have >= 2 informative individuals (sizes below are over these)\n", I, ph, mean(m[, "tested"])))
  gate(m, "snp_size", sprintf("paired outbred I = %d, phi %.3f, SNP", I, ph)); gate(m, "gene_size", sprintf("paired outbred I = %d, phi %.3f, gene (ACAT)", I, ph))
}
m <- seeds_rbind(1:10, function(s) pair_run(7000 + s, 4, 0.02, 0.1)); report(m, "paired outbred I = 4, 10% SNPs changed, phase-heterogeneous")
gate(m, "snp_size", "paired outbred with 10% changed SNPs: size on the unchanged SNPs")
pow_gate(m, "paired outbred power, phase-heterogeneous change 0.5 -> 0.8, I = 4")
untrimmed <- function(m, label) cat(sprintf("     EVIDENCE %s, untrimmed moment pair phi (phi_all): phi %.4f for a true 0.02, power %.4f (trimmed: phi %.4f, power %.4f; oracle %.4f)\n",
                                            label, mean(m[, "phi_untrimmed"]), mean(m[, "power_untrimmed"]), mean(m[, "phi"]), mean(m[, "power"]), mean(m[, "power_oracle"])))
untrimmed(m, "10% changed, phase-heterogeneous")
m <- seeds_rbind(1:10, function(s) pair_run(7100 + s, 4, 0.02, 0.1, TRUE)); report(m, "paired outbred I = 4, consistent change")
gate(m, "snp_size", "paired outbred with 10% consistently changed SNPs: size on the unchanged SNPs")
pow_gate(m, "paired outbred power, consistent change 0.5 -> 0.8, I = 4")
untrimmed(m, "10% changed, consistent")
ok(min(m[, "sign_up"]) >= 0.95, sprintf("paired outbred consistent change: mean_delta > 0 in %.3f of calls (>= 0.95)", min(m[, "sign_up"])))
m <- seeds_rbind(1:10, function(s) pair_run(7300 + s, 4, 0.02, 0.3)); report(m, "paired outbred I = 4, 30% SNPs changed, phase-heterogeneous")
cat(sprintf("     EVIDENCE 30%% of SNPs changed: power %.4f against oracle %.4f, pair phi %.4f for a true 0.02 (limit of the trimmed estimate)\n",
            mean(m[, "power"]), mean(m[, "power_oracle"]), mean(m[, "phi"])))
set.seed(7200); bg <- sim_pairs(300, 3, 0.01)
ed <- data.frame(snp = c("a", "a", "a", "a", "b", "b", "c", "c", "c"), individual = c("i1", "i1", "i2", "i2", "i1", "i1", "i3", "i3", "i1"),
                 cond = c(0, 1, 0, 1, 0, 1, 0, 0, 0), y = c(20, 40, 22, 38, 10, 12, 5, 6, 7), n = c(50, 50, 50, 50, 20, 20, 10, 10, 10), diff = FALSE)
pe <- ase_paired_test(rbind(bg, ed), RHO_MIN)
ok(identical(names(pe$snp), c("snp", "n_individuals", "n_failed", "stat", "df", "p", "mean_delta", "max_abs_delta", "n_up", "n_down", "status")) &&
   identical(names(pe$phi), c("individual", "phi_pair")), "ase_paired_test output columns")
pa <- pe$snp[pe$snp$snp == "a", ]; pb <- pe$snp[pe$snp$snp == "b", ]
ok(pa$status == "ok" && pa$n_individuals == 2 && pa$n_up == 2 && is.finite(pa$p), "SNP in 2 paired individuals: tested, both REF up")
ok(pb$status == "too_few_individuals" && is.na(pb$p), "SNP paired in 1 individual only: too_few_individuals, p NA")
ok(!"c" %in% pe$snp$snp, "a SNP without any paired individual is not reported as tested")
ok(all(pe$phi$phi_pair >= RHO_MIN), "per-individual pair dispersion never below RHO_MIN")

# ---- 8. fragment level: SNPs that share read pairs (ASEReadCounter counts a fragment at every SNP it covers)
frag_counts <- function(pos, p, nf, span, rl = 100) {
  L <- round(runif(nf, 150, 350)); st <- floor(runif(nf) * (span - L + 1)) + 1; A <- runif(nf) < p
  cv <- (outer(st, pos, "<=") & outer(st + rl - 1, pos, ">=")) | (outer(st + L - rl, pos, "<=") & outer(st + L - 1, pos, ">="))
  list(a = colSums(cv & A), n = colSums(cv))
}
frag_run <- function(seed, G = 1500, nA = 3, nB = 3, phi_bio = 0.005, span = 1500) {
  set.seed(seed); d <- c(rep(1, nA), rep(-1, nB)); S <- length(d)
  k <- sample(c(1, 2, 4, 8), G, TRUE); dense <- runif(G) < 0.3
  genes <- lapply(seq_len(G), function(g) {
    pos <- if (dense[g]) sort(sample(600:900, k[g])) else sort(sample(1:span, k[g]))       # dense: every SNP within 300 bp
    pr <- rbeta(S, 0.5 * (1 - phi_bio) / phi_bio, 0.5 * (1 - phi_bio) / phi_bio)        # null: animal-level variation only
    cnt <- lapply(seq_len(S), function(s) frag_counts(pos, pr[s], rpois(1, 150), span))
    list(pos = pos, a = sapply(cnt, `[[`, "a"), n = sapply(cnt, `[[`, "n"))               # SNP x sample matrices
  })
  mat <- function(v) if (is.matrix(v)) v else matrix(v, nrow = 1)
  unit <- function(gi, thin) {
    a <- mat(genes[[gi]]$a); n <- mat(genes[[gi]]$n)
    keep <- if (thin) thin_snps(genes[[gi]]$pos, rowSums(n), THIN_BP) else rep(TRUE, nrow(n))
    y <- colSums(a[keep, , drop = FALSE]); nn <- colSums(n[keep, , drop = FALSE]); ok_s <- nn > 0
    list(y = y[ok_s], n = nn[ok_s], X = cbind("(Intercept)" = 1, d = d)[ok_s, , drop = FALSE])
  }
  tests <- list(strain = "(Intercept)", parent_of_origin = "d")
  rt <- ase_glm_test(setNames(lapply(seq_len(G), unit, thin = TRUE), paste0("g", 1:G)), tests, RHO_MIN)
  ru <- ase_glm_test(setNames(lapply(seq_len(G), unit, thin = FALSE), paste0("g", 1:G)), tests, RHO_MIN)
  sz <- function(p, sel) mean(p[sel] < 0.05, na.rm = TRUE)
  multi <- k >= 2
  rmd02 <- function(thin, rho_thin = thin) {   # the Rmd 02 F1 gene path, per sample: free-mean rho (corrected, floored) + gene LRT
    unlist(lapply(seq_len(S), function(s) {   # thin: SNPs of the gene LRT; rho_thin: SNPs of the rho estimate
      ab <- rmd02_counts(genes, s, thin); r <- if (rho_thin == thin) ab else rmd02_counts(genes, s, rho_thin)
      e <- bb_estimate_rho_gene(r$x, r$n, r$gid); rho <- max(if (is.na(e[["corrected"]])) 0 else e[["corrected"]], RHO_MIN)
      vapply(which(multi), function(gi) { i <- ab$gid == gi; if (sum(i) < 2) NA_real_ else bb_gene_lrt(ab$x[i], ab$n[i], rho)$p }, numeric(1))
    }))
  }
  p02u <- rmd02(FALSE); p02t <- rmd02(TRUE); p02tu <- rmd02(TRUE, FALSE); dm <- rep(dense[multi], S)
  c(thin_strain = sz(rt$p_strain, TRUE), thin_poo = sz(rt$p_parent_of_origin, TRUE), thin_dense_poo = sz(rt$p_parent_of_origin, dense & multi),
    unthin_strain = sz(ru$p_strain, TRUE), unthin_poo = sz(ru$p_parent_of_origin, TRUE), unthin_dense_poo = sz(ru$p_parent_of_origin, dense & multi),
    rmd02_unthinned = mean(p02u < 0.05, na.rm = TRUE), rmd02_thinned = mean(p02t < 0.05, na.rm = TRUE),
    rmd02_dense_unthinned = mean(p02u[dm] < 0.05, na.rm = TRUE), rmd02_dense_thinned = mean(p02t[dm] < 0.05, na.rm = TRUE),
    rmd02_thinned_unthinned_rho = mean(p02tu < 0.05, na.rm = TRUE))
}
rmd02_counts <- function(genes, s, thin) {   # one sample's SNP counts (x, n, gene index), all SNPs or thinned (deepest per THIN_BP window)
  mat <- function(v) if (is.matrix(v)) v else matrix(v, nrow = 1)
  x <- integer(0); n <- integer(0); gid <- integer(0)
  for (gi in seq_along(genes)) {
    a <- mat(genes[[gi]]$a); nn <- mat(genes[[gi]]$n)
    keep <- if (thin) thin_snps(genes[[gi]]$pos, rowSums(nn), THIN_BP) else rep(TRUE, nrow(nn))
    x <- c(x, a[keep, s]); n <- c(n, nn[keep, s]); gid <- c(gid, rep(gi, sum(keep)))
  }
  list(x = x, n = n, gid = gid)
}
m <- seeds_rbind(1:10, function(s) frag_run(8000 + s)); report(m, "fragment-level null (shared read pairs)")
gate(m, "thin_strain", "fragment level, thinned sums, strain test"); gate(m, "thin_poo", "fragment level, thinned sums, parent-of-origin test")
gate(m, "thin_dense_poo", "fragment level, thinned sums, SNP-dense genes, parent-of-origin test", 0.07, 0.10)
cat(sprintf("UNTHINNED_F1_GENE_SIZE rmd02_unthinned=%.4f rmd02_thinned=%.4f rmd02_dense_unthinned=%.4f rmd02_dense_thinned=%.4f unthinned_sum_poo=%.4f\n",
            mean(m[, "rmd02_unthinned"]), mean(m[, "rmd02_thinned"]), mean(m[, "rmd02_dense_unthinned"]), mean(m[, "rmd02_dense_thinned"]), mean(m[, "unthin_poo"])))
cat(sprintf("     EVIDENCE same thinned Rmd 02 gene LRT with the UNTHINNED free-mean rho: size %.4f, worst seed %.4f (thinned rho: %.4f, worst seed %.4f)\n",
            mean(m[, "rmd02_thinned_unthinned_rho"]), max(m[, "rmd02_thinned_unthinned_rho"]), mean(m[, "rmd02_thinned"]), max(m[, "rmd02_thinned"])))

# ---- 8b. Rmd 02 rho_gene fallback: fewer than 5 genes keep 2 or more thinned SNPs, so the gene LRT on the thinned SNPs uses
#          the UNTHINNED free-mean rho (as on the Stage 1 synthetic acceptance data, where 2 of 12 genes keep 2 SNPs)
frag_fallback_run <- function(seed, G = 300, S = 6, phi_bio = 0.005, span = 1500) {
  set.seed(seed); k <- sample(c(2, 4, 8), G, TRUE); sparse <- seq_len(G) <= 3   # only 3 genes can keep 2 or more SNPs
  genes <- lapply(seq_len(G), function(g) {
    pos <- if (sparse[g]) sort(sample(1:span, 4)) else sort(sample(600:900, k[g]))   # dense: every SNP within 300 bp
    pr <- rbeta(S, 0.5 * (1 - phi_bio) / phi_bio, 0.5 * (1 - phi_bio) / phi_bio)
    cnt <- lapply(seq_len(S), function(s) frag_counts(pos, pr[s], rpois(1, 150), span))
    list(pos = pos, a = sapply(cnt, `[[`, "a"), n = sapply(cnt, `[[`, "n"))
  })
  per <- lapply(seq_len(S), function(s) {
    th <- rmd02_counts(genes, s, TRUE); un <- rmd02_counts(genes, s, FALSE)
    et <- bb_estimate_rho_gene(th$x, th$n, th$gid)[["corrected"]]; eu <- bb_estimate_rho_gene(un$x, un$n, un$gid)[["corrected"]]
    fb <- is.na(et) && !is.na(eu)                                   # step (2) of the Rmd 02 rho_gene chain is taken
    rho <- if (!is.na(et)) max(et, RHO_MIN) else if (!is.na(eu)) max(eu, RHO_MIN) else max(bb_estimate_rho(un$x, un$n), RHO_MIN)
    p <- vapply(seq_len(G), function(gi) { i <- th$gid == gi; if (sum(i) < 1) NA_real_ else bb_gene_lrt(th$x[i], th$n[i], rho)$p }, numeric(1))
    list(p = p, fb = fb, rho = rho, eu = eu)
  })
  p <- unlist(lapply(per, `[[`, "p"))
  c(fallback_size = mean(p < 0.05, na.rm = TRUE), fallback_dense_size = mean(p[rep(!sparse, S)] < 0.05, na.rm = TRUE),
    fallback_taken = mean(vapply(per, `[[`, logical(1), "fb")), rho_used = mean(vapply(per, `[[`, numeric(1), "rho")),
    rho_unthinned_corrected = mean(vapply(per, `[[`, numeric(1), "eu")))
}
m <- seeds_rbind(1:10, function(s) frag_fallback_run(8100 + s)); report(m, "Rmd 02 rho_gene fallback (unthinned rho, thinned gene LRT)")
ok(all(m[, "fallback_taken"] == 1), sprintf("fallback scenario: the unthinned rho_gene is used in every sample (fraction %.2f)", mean(m[, "fallback_taken"])))
gate(m, "fallback_size", "Rmd 02 F1 gene LRT with the rho_gene fallback to the unthinned free-mean rho: null size")
gate(m, "fallback_dense_size", "Rmd 02 F1 gene LRT with the rho_gene fallback, SNP-dense genes (one SNP left): null size")
cat(sprintf("RHO_GENE_FALLBACK fallback_size=%.4f worst=%.4f dense=%.4f dense_worst=%.4f rho_used=%.4f rho_unthinned_corrected=%.4f\n",
            mean(m[, "fallback_size"]), max(m[, "fallback_size"]), mean(m[, "fallback_dense_size"]), max(m[, "fallback_dense_size"]),
            mean(m[, "rho_used"]), mean(m[, "rho_unthinned_corrected"])))

# ---- 9. runtime (one core): projected time for 60,000 units x 12 rows must fit the 4 h run-script limit with margin
set.seed(900); ub <- sim_units(5000, rep(c(1, -1), 6), 0.02)
tt <- system.time(rb <- ase_glm_test(ub, list(strain = "(Intercept)", parent_of_origin = "d"), RHO_MIN))[["elapsed"]]
proj <- tt * 60000 / 5000 / 3600
ok(proj < 3, sprintf("runtime: 5000 units x 12 rows, 2 tests, one core: %.0f s; projected 60,000 units %.2f h (< 3 h)", tt, proj))

quit(status = fail)
