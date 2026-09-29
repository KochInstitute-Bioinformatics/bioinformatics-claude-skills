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
# the shipped default of RHO_MIN is read from the Step 8 table of the skill text (not hard-coded here)
rm_line <- grep("^\\| `RHO_MIN` \\|", txt, value = TRUE)
RHO_MIN <- as.numeric(trimws(strsplit(rm_line[1], "\\|")[[1]][3])); stopifnot(length(RHO_MIN) == 1, is.finite(RHO_MIN))
ok(abs(RHO_MIN - 0.01) < 1e-12, sprintf("RHO_MIN default read from the skill = %g", RHO_MIN))
# pure-binomial ML sits at the boundary only about half of the time (sd of rho-hat about 6.5e-4): over 10 seeds require
# most estimates (>= 40%) to be exactly 0 and every estimate to be small
zr <- sapply(101:110, function(s) { set.seed(s); bb_estimate_rho(rbinom(2000, 50, 0.5), rep(50, 2000)) })
ok(mean(zr == 0) >= 0.4 && all(zr < 2e-3), sprintf("binomial data, 10 seeds, bb_estimate_rho: %d exactly 0, max %.2e", sum(zr == 0), max(zr)))
gb <- rep(1:500, each = 4)
zg <- sapply(101:110, function(s) { set.seed(s); bb_estimate_rho_gene(rbinom(2000, 60, 0.5), rep(60, 2000), gb) })
ok(mean(zg["free", ] == 0) >= 0.4 && all(zg["free", ] < 2e-3), sprintf("binomial data, 10 seeds, bb_estimate_rho_gene free: %d exactly 0, max %.2e", sum(zg["free", ] == 0), max(zg["free", ])))
ok(mean(zg["corrected", ] == 0) >= 0.4 && all(zg["corrected", ] < RHO_MIN), sprintf("binomial data, 10 seeds, bb_estimate_rho_gene corrected: %d exactly 0, max %.4f (below RHO_MIN)", sum(zg["corrected", ] == 0), max(zg["corrected", ])))
for (xn in list(c(25, 50), c(30, 50), c(3, 10), c(0, 20), c(7, 9), c(60, 100), c(1, 2)))
  ok(abs(bb_pvalue(xn[1], xn[2], 0) - binom.test(xn[1], xn[2], 0.5)$p.value) < 1e-8, sprintf("bb_pvalue(%d, %d, 0) equals binom.test", xn[1], xn[2]))
ok(is.na(bb_pvalue(11, 10, 0.02)) && is.na(bb_pvalue(-1, 10, 0.02)) && is.na(bb_pvalue(2.5, 10, 0.02)) && is.na(bb_pvalue(NA, 10, 0.02)), "bb_pvalue guards x > n, x < 0, non-integer and NA x")

# 6. F1 free-mean overdispersion (df-corrected, floored) and gene-level LRT
set.seed(2)
rho_t <- 0.02; nn <- 60
a2 <- function(p) p * (1 - rho_t) / rho_t
sim <- function(k, G_null = 500, frac_imb = 0.4, nfun = function(m) rep(nn, m)) {
  G_imb <- round(G_null * frac_imb / (1 - frac_imb)); pg <- c(rep(0.5, G_null), rep(0.7, G_imb)); G <- length(pg)
  pk <- rep(pg, each = k); n <- nfun(G * k)
  list(x = rbinom(G * k, n, rbeta(G * k, a2(pk), a2(1 - pk))), n = n, gene = rep(seq_len(G), each = k), G_null = G_null)
}
s4 <- sim(4); e4 <- bb_estimate_rho_gene(s4$x, s4$n, s4$gene)
r_h0 <- bb_estimate_rho(s4$x, s4$n)
cat(sprintf("     4-SNP mixed data: uncorrected %.4f, corrected %.4f, H0-based %.4f\n", e4[["free"]], e4[["corrected"]], r_h0))
ok(e4[["corrected"]] > 0.014 && e4[["corrected"]] < 0.03, sprintf("corrected free-mean rho %.4f within (0.014, 0.03) of 0.02 with 40%% imbalanced genes", e4[["corrected"]]))
ok(e4[["corrected"]] > e4[["free"]], "correction raises the free-mean estimate")
ok(r_h0 > 2 * e4[["corrected"]], sprintf("H0-based rho %.4f is inflated (> 2 x corrected %.4f)", r_h0, e4[["corrected"]]))
ok(all(is.na(bb_estimate_rho_gene(s4$x[1:12], rep(nn, 12), s4$gene[1:12]))), "fewer than 5 usable genes gives NA")
# full path (estimate rho on the mixed data, correct, floor, gene LRT on the NULL genes only), 10 seeds x 2000 null genes
size_run <- function(seed, k, nfun) {
  set.seed(seed); s <- sim(k, 2000, 0.4, nfun); e <- bb_estimate_rho_gene(s$x, s$n, s$gene)
  rg <- max(e[["corrected"]], RHO_MIN)
  pn <- vapply(seq_len(s$G_null), function(g) { i <- ((g - 1) * k + 1):(g * k); bb_gene_lrt(s$x[i], s$n[i], rg)$p }, numeric(1))
  c(size = mean(pn < 0.05), free = e[["free"]], corrected = e[["corrected"]], rho_gene = rg)
}
cases <- list(list("2-SNP, n = 60", 2, NULL), list("4-SNP, n = 60", 4, NULL),
              list("2-SNP, heterogeneous depth 30-200", 2, function(m) round(runif(m, 30, 200))),
              list("4-SNP, heterogeneous depth 30-200", 4, function(m) round(runif(m, 30, 200))))
for (cs in cases) {
  nf <- if (is.null(cs[[3]])) function(m) rep(nn, m) else cs[[3]]
  res <- parallel::mclapply(1:10, function(s) size_run(1000 + s, cs[[2]], nf), mc.cores = 8)
  stopifnot(!any(vapply(res, function(r) inherits(r, "try-error"), logical(1))))
  m <- do.call(rbind, res)
  cat(sprintf("     %s: mean uncorrected %.4f, mean corrected %.4f, mean rho_gene %.4f, sizes %s\n", cs[[1]], mean(m[, "free"]), mean(m[, "corrected"]), mean(m[, "rho_gene"]), paste(sprintf("%.3f", m[, "size"]), collapse = " ")))
  ok(mean(m[, "size"]) <= 0.065 && max(m[, "size"]) <= 0.08, sprintf("full-path null size, %s: pooled %.4f (<= 0.065), worst seed %.3f (<= 0.08)", cs[[1]], mean(m[, "size"]), max(m[, "size"])))
}
pw <- sapply(rep(c(0.7, 0.3), each = 250), function(p) bb_gene_lrt(rbinom(4, nn, rbeta(4, a2(p), a2(1 - p))), rep(nn, 4), rho_t)$p)
ok(mean(pw < 0.05) > 0.9, sprintf("gene LRT power at p = 0.7 / 0.3 = %.2f", mean(pw < 0.05)))
ok(bb_gene_lrt(c(50, 48, 52), c(60, 60, 60), rho_t)$phat > 0.5, "gene LRT phat above 0.5 for ALT-rich counts")
ok(is.na(bb_gene_lrt(c(0, 0), c(0, 0), rho_t)$p), "all-zero counts give NA")
s1 <- bb_gene_lrt(40, 60, rho_t); ok(is.finite(s1$p) && s1$p >= 0 && s1$p <= 1, "single-SNP gene works")
ok(is.finite(bb_gene_lrt(c(30, 35), c(60, 60), 0)$p), "gene LRT with rho = 0 works")
b1 <- bb_gene_lrt(c(0, 0, 0), c(60, 60, 60), rho_t); ok(is.finite(b1$p) && b1$p < 1e-6 && b1$phat < 0.01, "boundary p-hat (all-REF gene) does not error")
b2 <- bb_gene_lrt(c(60, 60), c(60, 60), 0.05); ok(is.finite(b2$p) && b2$phat > 0.99, "boundary p-hat (all-ALT gene) does not error")

# 7. Outbred: trimmed (central-sites, truncation-corrected) rho, per-SNP bb_pvalue and gene-level ACAT, unphased data
#    (task 8 fix D1: the all-sites H0 rho was inflated by real imbalance and detected 0 of 16 planted genes)
stopifnot(exists("bb_estimate_rho_trim"))
dep_syn <- function(m) round(runif(m, 30, 200))          # depth range of the synthetic acceptance data
sim_ob <- function(G, rho_t, frac_imb, depth = dep_syn) {
  k <- sample(1:4, G, replace = TRUE); imb <- runif(G) < frac_imb
  pg <- ifelse(imb, sample(c(0.85, 0.7, 0.3), G, replace = TRUE), 0.5)
  gene <- rep(seq_len(G), k); p <- pg[gene]; p <- ifelse(runif(length(p)) < 0.5, 1 - p, p)   # unphased: orientation random per SNP
  n <- depth(length(p))
  pr <- if (rho_t > 0) rbeta(length(p), p * (1 - rho_t) / rho_t, (1 - p) * (1 - rho_t) / rho_t) else p
  list(x = rbinom(length(p), n, pr), n = n, gene = gene, imb = imb, pg = pg)
}
ob_path <- function(s, rho) {   # the Rmd 02 outbred path: floor, per-SNP test, ACAT per gene, BH within the sample
  r <- max(rho, RHO_MIN); p <- mapply(bb_pvalue, s$x, s$n, MoreArgs = list(rho = r))
  gp <- tapply(p, s$gene, acat); list(p = p, gp = gp, gpadj = p.adjust(gp, "BH"))
}
ob_run <- function(seed, rho_t, frac_imb, G = 300) {
  set.seed(seed); s <- sim_ob(G, rho_t, frac_imb)
  tr <- bb_estimate_rho_trim(s$x, s$n); h0 <- bb_estimate_rho(s$x, s$n)
  a <- ob_path(s, tr[["rho"]]); b <- ob_path(s, h0); nullsnp <- !s$imb[s$gene]
  c(rho_trim = tr[["rho"]], rho_h0 = h0, snp_size = mean(a$p[nullsnp] < 0.05), gene_size = mean(a$gp[!s$imb] < 0.05),
    pow85 = mean(a$gpadj[s$imb & s$pg == 0.85] < 0.05), pow70_30 = mean(a$gpadj[s$imb & s$pg != 0.85] < 0.05),
    h0_pow85 = mean(b$gpadj[s$imb & s$pg == 0.85] < 0.05))
}
ob_cases <- list(list("binomial, 1/3 genes imbalanced", 0, 1/3), list("rho 0.01, 1/3 imbalanced", 0.01, 1/3),
                 list("rho 0.02, 1/3 imbalanced", 0.02, 1/3), list("rho 0.02, no imbalance", 0.02, 0),
                 list("rho 0.05, no imbalance", 0.05, 0), list("rho 0.02, 10% imbalanced", 0.02, 0.1))
ob <- list()
for (cs in ob_cases) {
  res <- parallel::mclapply(1:10, function(s) ob_run(2000 + s, cs[[2]], cs[[3]]), mc.cores = 8)
  stopifnot(!any(vapply(res, function(r) inherits(r, "try-error"), logical(1))))
  m <- do.call(rbind, res); ob[[cs[[1]]]] <- m
  cat(sprintf("     outbred %s: rho_trim %.4f (%.4f-%.4f), rho_h0 %.4f, SNP size %.4f, gene size %.4f (worst %.3f), power 0.85 %.3f, power 0.70/0.30 %.3f, H0-rho power 0.85 %.3f\n",
              cs[[1]], mean(m[, "rho_trim"]), min(m[, "rho_trim"]), max(m[, "rho_trim"]), mean(m[, "rho_h0"]), mean(m[, "snp_size"]),
              mean(m[, "gene_size"]), max(m[, "gene_size"]), mean(m[, "pow85"], na.rm = TRUE), mean(m[, "pow70_30"], na.rm = TRUE), mean(m[, "h0_pow85"], na.rm = TRUE)))
  ok(mean(m[, "snp_size"]) <= 0.06 && mean(m[, "gene_size"]) <= 0.06 && max(m[, "gene_size"]) <= 0.10,
     sprintf("outbred null size, %s: SNP %.4f, gene %.4f pooled (<= 0.06), worst seed %.3f (<= 0.10)", cs[[1]], mean(m[, "snp_size"]), mean(m[, "gene_size"]), max(m[, "gene_size"])))
  if (cs[[2]] > 0) ok(min(m[, "rho_trim"]) > 0.5 * cs[[2]] && max(m[, "rho_trim"]) < 0.5,
                      sprintf("rho_trim neither collapses nor blows up, %s: %.4f-%.4f", cs[[1]], min(m[, "rho_trim"]), max(m[, "rho_trim"])))
}
m <- ob[["binomial, 1/3 genes imbalanced"]]
ok(mean(m[, "pow85"]) >= 0.95 && mean(m[, "pow70_30"]) >= 0.75,
   sprintf("outbred power at the synthetic depth (binomial): 0.85 genes %.3f (>= 0.95), 0.70/0.30 genes %.3f (>= 0.75)", mean(m[, "pow85"]), mean(m[, "pow70_30"])))
ok(mean(m[, "h0_pow85"]) < 0.1, sprintf("the all-sites H0 rho (%.3f) kills power in the same data: 0.85 genes %.3f", mean(m[, "rho_h0"]), mean(m[, "h0_pow85"])))
m <- ob[["rho 0.01, 1/3 imbalanced"]]
ok(mean(m[, "pow85"]) >= 0.9 && mean(m[, "rho_trim"]) < mean(m[, "rho_h0"]) / 3,
   sprintf("rho 0.01 with imbalance: power 0.85 genes %.3f (>= 0.9); rho_trim %.4f < rho_h0 %.4f / 3", mean(m[, "pow85"]), mean(m[, "rho_trim"]), mean(m[, "rho_h0"])))
for (nm in c("rho 0.02, no imbalance", "rho 0.05, no imbalance")) {
  m <- ob[[nm]]; tru <- if (grepl("0.05", nm)) 0.05 else 0.02
  ok(abs(mean(m[, "rho_trim"]) / tru - 1) < 0.15, sprintf("%s: mean rho_trim %.4f within 15%% of %.2f", nm, mean(m[, "rho_trim"]), tru))
}
# synthetic-acceptance-sized sample: 12 genes (4 planted 0.85/0.70/0.30/0.70), about 30 sites, binomial, 100 seeds
small <- do.call(rbind, parallel::mclapply(1:100, function(sd) { set.seed(3000 + sd)
  k <- rep(c(4, 2, 3), 4); pg <- c(0.85, 0.7, 0.3, 0.7, rep(0.5, 8)); gene <- rep(1:12, k); p <- pg[gene]
  p <- ifelse(runif(length(p)) < 0.5, 1 - p, p); n <- dep_syn(length(p)); x <- rbinom(length(p), n, p)
  tr <- bb_estimate_rho_trim(x, n); a <- ob_path(list(x = x, n = n, gene = gene), tr[["rho"]])
  c(rho = tr[["rho"]], det = sum(a$gpadj[1:4] < 0.05), fp = sum(a$gpadj[5:12] < 0.05)) }, mc.cores = 8))
cat(sprintf("     small (30 sites): rho_trim median %.4f, 95th pct %.4f, max %.4f; planted detected %.2f of 4; null FP %.3f of 8\n",
            median(small[, "rho"]), quantile(small[, "rho"], 0.95), max(small[, "rho"]), mean(small[, "det"]), mean(small[, "fp"])))
# with only about 30 sites (11 imbalanced) the trimmed fit sometimes cannot separate them and stays near the naive value
# (conservative); it must never run away to a U-shaped fit (which would make central counts significant)
ok(all(is.finite(small[, "rho"])) && median(small[, "rho"]) < 0.01 && max(small[, "rho"]) < 0.3,
   sprintf("small samples: rho_trim finite, median %.4f < 0.01, max %.4f < 0.3 (no run-away)", median(small[, "rho"]), max(small[, "rho"])))
ok(mean(small[, "det"]) >= 3 && mean(small[, "fp"]) / 8 <= 0.05,
   sprintf("small samples: %.2f of 4 planted genes detected (>= 3), null FP rate %.4f (<= 0.05)", mean(small[, "det"]), mean(small[, "fp"]) / 8))
ok(all(is.na(bb_estimate_rho_trim(rbinom(19, 50, 0.5), rep(50, 19)))), "rho_trim with fewer than 20 sites gives NA")
zt <- sapply(101:110, function(s) { set.seed(s); bb_estimate_rho_trim(rbinom(2000, 50, 0.5), rep(50, 2000))[["rho"]] })
ok(mean(zt == 0) >= 0.4 && all(zt < 3e-3), sprintf("binomial data, 10 seeds, rho_trim: %d exactly 0, max %.2e", sum(zt == 0), max(zt)))
w <- NULL; withCallingHandlers({ set.seed(7); s <- sim_ob(3000, 0.02, 0.2, function(m) round(exp(runif(m, log(10), log(1000))))); tt <- system.time(bb_estimate_rho_trim(s$x, s$n)) },
                               warning = function(cond) { w <<- conditionMessage(cond); invokeRestart("muffleWarning") })
ok(is.null(w), sprintf("rho_trim converges without a warning on %d sites at depth 10-1000 (%.1f s)", length(s$x), tt[["elapsed"]]))

quit(status = fail)
