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
quit(status = fail)
