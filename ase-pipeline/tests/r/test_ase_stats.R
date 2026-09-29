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
ok(abs(acat(c(1, 1, 1)) - 1) < 1e-6, "acat of all ones is 1")
ok(acat(c(1e-10, 0.5, 0.5)) < 1e-6, "acat dominated by a tiny p")
ok(acat(c(0.2, 0.3, 0.4)) > 0.1 && acat(c(0.2, 0.3, 0.4)) < 0.5, "acat of moderate p-values is moderate")

# 6. F1 free-mean overdispersion and gene-level LRT
set.seed(2)
G <- 300; k <- 4; nn <- 60; rho_t <- 0.02; a2 <- function(p) p * (1 - rho_t) / rho_t
pg <- ifelse(seq_len(G) <= 0.4 * G, 0.7, 0.5)
gid <- rep(seq_len(G), each = k)
xg <- unlist(lapply(pg, function(p) rbinom(k, nn, rbeta(k, a2(p), a2(1 - p)))))
r_free <- bb_estimate_rho_gene(xg, rep(nn, G * k), gid); r_h0 <- bb_estimate_rho(xg, rep(nn, G * k))
ok(r_free > 0.01 && r_free < 0.04, sprintf("free-mean rho %.4f near 0.02 with 40%% imbalanced genes", r_free))
ok(r_h0 > 2 * r_free, sprintf("H0-based rho %.4f is inflated (> 2 x free-mean %.4f)", r_h0, r_free))
ok(is.na(bb_estimate_rho_gene(xg[1:12], rep(nn, 12), gid[1:12])), "fewer than 5 usable genes gives NA")
pn <- sapply(1:500, function(i) bb_gene_lrt(rbinom(k, nn, rbeta(k, a2(0.5), a2(0.5))), rep(nn, k), rho_t)$p)
ok(abs(mean(pn < 0.05) - 0.05) < 0.02, sprintf("gene LRT null size at 0.05 = %.3f", mean(pn < 0.05)))
pw <- sapply(rep(c(0.7, 0.3), each = 250), function(p) bb_gene_lrt(rbinom(k, nn, rbeta(k, a2(p), a2(1 - p))), rep(nn, k), rho_t)$p)
ok(mean(pw < 0.05) > 0.9, sprintf("gene LRT power at p = 0.7 / 0.3 = %.2f", mean(pw < 0.05)))
ok(bb_gene_lrt(c(50, 48, 52), c(60, 60, 60), rho_t)$phat > 0.5, "gene LRT phat above 0.5 for ALT-rich counts")
ok(is.na(bb_gene_lrt(c(0, 0), c(0, 0), rho_t)$p), "all-zero counts give NA")
s1 <- bb_gene_lrt(40, 60, rho_t); ok(is.finite(s1$p) && s1$p >= 0 && s1$p <= 1, "single-SNP gene works")
ok(is.finite(bb_gene_lrt(c(30, 35), c(60, 60), 0)$p), "gene LRT with rho = 0 works")
b1 <- bb_gene_lrt(c(0, 0, 0), c(60, 60, 60), rho_t); ok(is.finite(b1$p) && b1$p < 1e-6 && b1$phat < 0.01, "boundary p-hat (all-REF gene) does not error")
b2 <- bb_gene_lrt(c(60, 60), c(60, 60), 0.05); ok(is.finite(b2$p) && b2$phat > 0.99, "boundary p-hat (all-ALT gene) does not error")

quit(status = fail)
