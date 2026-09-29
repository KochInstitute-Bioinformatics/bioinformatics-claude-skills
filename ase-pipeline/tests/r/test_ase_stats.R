args <- commandArgs(trailingOnly = TRUE); skill <- args[1]
txt <- readLines(skill)
b <- grep("^# --- ase-stats-begin", txt); e <- grep("^# --- ase-stats-end", txt)
stopifnot(length(b) == 1, length(e) == 1, e > b)
eval(parse(text = txt[(b + 1):(e - 1)]))
set.seed(1); fail <- 0
ok <- function(cond, msg) { if (!isTRUE(cond)) { cat("FAIL:", msg, "\n"); fail <<- 1 } else cat("ok  ", msg, "\n") }

# 1. null calibration: p-values approximately uniform under H0 (rho = 0.02, n = 50)
rho0 <- 0.02; n <- rep(50, 3000); a <- 0.5 * (1 - rho0) / rho0
x <- rbinom(3000, n, rbeta(3000, a, a))
p <- mapply(bb_pvalue, x, n, MoreArgs = list(rho = rho0))
ok(abs(mean(p < 0.05) - 0.05) < 0.02, sprintf("null size at 0.05 = %.3f", mean(p < 0.05)))
ok(ks.test(p, "punif")$p.value > 0.001 || abs(mean(p) - 0.5) < 0.05, "null p-values roughly uniform")

# 2. power: alt fraction 0.7, n = 60, low overdispersion
x1 <- rbinom(500, 60, 0.7); p1 <- sapply(x1, bb_pvalue, n = 60, rho = 0.02)
ok(mean(p1 < 0.05) > 0.5, sprintf("power at p=0.7 = %.2f", mean(p1 < 0.05)))

# 3. rho estimation recovers roughly the truth and handles rho -> 0
xx <- rbinom(2000, 50, rbeta(2000, a, a)); est <- bb_estimate_rho(xx, rep(50, 2000))
ok(est > 0.01 && est < 0.04, sprintf("rho estimate %.4f near 0.02", est))
xb <- rbinom(2000, 50, 0.5); estb <- bb_estimate_rho(xb, rep(50, 2000))
ok(estb < 0.005, sprintf("binomial data gives rho %.5f near 0", estb))
ok(is.finite(bb_pvalue(25, 50, rho = 0)), "rho = 0 falls back to the binomial")

# 4. edge cases: zero total returns NA; x = 0 and x = n give small p at high depth
ok(is.na(bb_pvalue(0, 0, rho = 0.02)), "n = 0 gives NA")
ok(bb_pvalue(0, 80, rho = 0.02) < 1e-6, "all-ref at depth 80 is significant")

# 5. ACAT combination
ok(acat(c(1, 1, 1)) >= 0.98, "acat of all ones stays near 1 (p capped at 0.99)")
ok(acat(c(1e-10, 0.5, 0.5)) < 1e-6, "acat dominated by a tiny p")
ok(acat(c(0.2, 0.3, 0.4)) > 0.1 && acat(c(0.2, 0.3, 0.4)) < 0.5, "acat of moderate p-values is moderate")

# 5b. ACAT must not collapse to ~1 when a SNP p-value is exactly 1 (bb_pvalue returns 1 at the modal count)
ok(acat(c(1e-6, 1)) < 1e-3, sprintf("acat(1e-6, 1) = %.2e stays significant", acat(c(1e-6, 1))))
ok(acat(c(1e-6, 0.999)) < 1e-3, sprintf("acat(1e-6, 0.999) = %.2e stays significant", acat(c(1e-6, 0.999))))
ok(is.na(acat(c(NA_real_, NA_real_))), "acat of only NA gives NA")

# 3b. boundary handling and exact zero
xb2 <- rbinom(2000, 50, 0.5); nb2 <- rep(50, 2000)
ok(identical(bb_estimate_rho(xb2, nb2), 0), "binomial data gives exactly 0 from bb_estimate_rho")
gb <- rep(1:500, each = 4); xg0 <- rbinom(2000, 60, 0.5)
eg0 <- bb_estimate_rho_gene(xg0, rep(60, 2000), gb)
ok(identical(unname(eg0[["free"]]), 0) && identical(unname(eg0[["corrected"]]), 0), "binomial data gives exactly 0 from bb_estimate_rho_gene")
for (xn in list(c(25, 50), c(30, 50), c(3, 10), c(0, 20), c(7, 9), c(60, 100), c(1, 2)))
  ok(abs(bb_pvalue(xn[1], xn[2], 0) - binom.test(xn[1], xn[2], 0.5)$p.value) < 1e-8, sprintf("bb_pvalue(%d, %d, 0) equals binom.test", xn[1], xn[2]))
ok(is.na(bb_pvalue(11, 10, 0.02)) && is.na(bb_pvalue(-1, 10, 0.02)) && is.na(bb_pvalue(2.5, 10, 0.02)) && is.na(bb_pvalue(NA, 10, 0.02)), "bb_pvalue guards x > n, x < 0, non-integer and NA x")

# 6. F1 free-mean overdispersion (df-corrected, floored) and gene-level LRT
set.seed(2)
RHO_MIN <- 0.01; rho_t <- 0.02; nn <- 60
a2 <- function(p) p * (1 - rho_t) / rho_t
sim <- function(k, G_null = 500, frac_imb = 0.4) {
  G_imb <- round(G_null * frac_imb / (1 - frac_imb)); pg <- c(rep(0.5, G_null), rep(0.7, G_imb)); G <- length(pg)
  list(x = unlist(lapply(pg, function(p) rbinom(k, nn, rbeta(k, a2(p), a2(1 - p))))), n = rep(nn, G * k),
       gene = rep(seq_len(G), each = k), G_null = G_null)
}
s4 <- sim(4); e4 <- bb_estimate_rho_gene(s4$x, s4$n, s4$gene)
r_h0 <- bb_estimate_rho(s4$x, s4$n)
cat(sprintf("     4-SNP mixed data: uncorrected %.4f, corrected %.4f, H0-based %.4f\n", e4[["free"]], e4[["corrected"]], r_h0))
ok(e4[["corrected"]] > 0.014 && e4[["corrected"]] < 0.03, sprintf("corrected free-mean rho %.4f within (0.014, 0.03) of 0.02 with 40%% imbalanced genes", e4[["corrected"]]))
ok(e4[["corrected"]] > e4[["free"]], "correction raises the free-mean estimate")
ok(r_h0 > 2 * e4[["corrected"]], sprintf("H0-based rho %.4f is inflated (> 2 x corrected %.4f)", r_h0, e4[["corrected"]]))
ok(all(is.na(bb_estimate_rho_gene(s4$x[1:12], rep(nn, 12), s4$gene[1:12]))), "fewer than 5 usable genes gives NA")
for (k in c(2, 4)) {
  s <- if (k == 4) s4 else sim(2); e <- if (k == 4) e4 else bb_estimate_rho_gene(s$x, s$n, s$gene)
  rg <- max(e[["corrected"]], RHO_MIN)
  pn <- sapply(seq_len(s$G_null), function(g) { i <- which(s$gene == g); bb_gene_lrt(s$x[i], s$n[i], rg)$p })
  cat(sprintf("     %d-SNP genes: uncorrected %.4f, corrected %.4f, rho_gene %.4f, null size %.3f\n", k, e[["free"]], e[["corrected"]], rg, mean(pn < 0.05)))
  ok(mean(pn < 0.05) <= 0.07, sprintf("full-path null size at 0.05 for %d-SNP genes = %.3f (<= 0.07)", k, mean(pn < 0.05)))
}
pw <- sapply(rep(c(0.7, 0.3), each = 250), function(p) bb_gene_lrt(rbinom(4, nn, rbeta(4, a2(p), a2(1 - p))), rep(nn, 4), rho_t)$p)
ok(mean(pw < 0.05) > 0.9, sprintf("gene LRT power at p = 0.7 / 0.3 = %.2f", mean(pw < 0.05)))
ok(bb_gene_lrt(c(50, 48, 52), c(60, 60, 60), rho_t)$phat > 0.5, "gene LRT phat above 0.5 for ALT-rich counts")
ok(is.na(bb_gene_lrt(c(0, 0), c(0, 0), rho_t)$p), "all-zero counts give NA")
s1 <- bb_gene_lrt(40, 60, rho_t); ok(is.finite(s1$p) && s1$p >= 0 && s1$p <= 1, "single-SNP gene works")
ok(is.finite(bb_gene_lrt(c(30, 35), c(60, 60), 0)$p), "gene LRT with rho = 0 works")
b1 <- bb_gene_lrt(c(0, 0, 0), c(60, 60, 60), rho_t); ok(is.finite(b1$p) && b1$p < 1e-6 && b1$phat < 0.01, "boundary p-hat (all-REF gene) does not error")
b2 <- bb_gene_lrt(c(60, 60), c(60, 60), 0.05); ok(is.finite(b2$p) && b2$phat > 0.99, "boundary p-hat (all-ALT gene) does not error")

quit(status = fail)
