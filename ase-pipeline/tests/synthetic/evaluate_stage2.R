# Usage: Rscript evaluate_stage2.R <reciprocal|differential_f1|differential_outbred> <RESULTS_DIR> <truth dir (f1_recip or outbred_diff)> <values.tsv>
# values.tsv: the NAME<TAB>value file the Rmds were generated from (render_from_skill.sh); THIN_BP and the strain names come from it.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("usage: Rscript evaluate_stage2.R <what> <RESULTS_DIR> <truth dir> <values.tsv>", call. = FALSE)
what <- args[1]; RES <- args[2]; TRUTH <- args[3]; VALS <- args[4]
fail <- 0
crit <- function(cond, msg) { cat(if (isTRUE(cond)) "PASS" else "FAIL", msg, "\n"); if (!isTRUE(cond)) fail <<- 1 }
tg <- read.delim(file.path(TRUTH, "truth_genes.tsv"), stringsAsFactors = FALSE)
vv <- read.delim(VALS, header = FALSE, colClasses = "character", quote = "", col.names = c("name", "value"))
val <- function(k) { v <- vv$value[vv$name == k]; if (length(v) != 1) stop("values file: need exactly one ", k, call. = FALSE); v }
if (what == "reciprocal") {
  THIN_BP <- as.numeric(val("THIN_BP")); SA <- val("STRAIN_A"); SB <- val("STRAIN_B")
  ck <- readRDS(file.path(RES, "ase_reciprocal_checkpoint.rds"))
  g <- ck$gene
  crit(isTRUE(all.equal(ck$constants$THIN_BP, THIN_BP)) && identical(ck$constants$STRAIN_A, SA) && identical(ck$constants$STRAIN_B, SB),
       sprintf("checkpoint constants match the values file (THIN_BP %s, STRAIN_A %s, STRAIN_B %s)", THIN_BP, SA, SB))
  m <- merge(tg, g, by = "gene_id", all.x = TRUE)
  for (k in c("sig_strain", "sig_parent_of_origin")) m[[k]][is.na(m[[k]])] <- FALSE
  print(m[m$class != "null", c("gene_id", "class", "expect_strain", "expect_poo", "frac_A", "maternal_frac", "padj_strain",
                               "padj_parent_of_origin", "direction_strain", "direction_parent_of_origin", "status")], row.names = FALSE)
  # every criterion below needs at least one truth gene of its class: a relabelled truth table must FAIL, not PASS as "0 of 0"
  exS <- m$expect_strain %in% c("A_higher", "B_higher")
  okS <- m$sig_strain & ((m$expect_strain == "A_higher" & m$b0_strain > 0) | (m$expect_strain == "B_higher" & m$b0_strain < 0))
  crit(sum(exS) > 0 && all(okS[exS]), sprintf("strain effect recovered with the correct sign: %d of %d", sum(okS[exS]), sum(exS)))
  exP <- m$expect_poo %in% c("maternal", "paternal")
  okP <- m$sig_parent_of_origin & ((m$expect_poo == "maternal" & m$b1_parent_of_origin > 0) | (m$expect_poo == "paternal" & m$b1_parent_of_origin < 0))
  crit(sum(exP) > 0 && all(okP[exP]), sprintf("parent-of-origin effect recovered with the correct sign: %d of %d", sum(okP[exP]), sum(exP)))
  so <- m$class %in% c("strain_A_high", "strain_B_high"); po <- m$class %in% c("maternal", "paternal")
  crit(sum(so) > 0 && sum(m$sig_parent_of_origin[so]) == 0,
       sprintf("no parent-of-origin call among strain-only genes: %d of %d", sum(m$sig_parent_of_origin[so]), sum(so)))
  crit(sum(po) > 0 && sum(m$sig_strain[po]) == 0, sprintf("no strain call among parent-of-origin-only genes: %d of %d", sum(m$sig_strain[po]), sum(po)))
  nul <- m$class %in% c("null", "dense_null", "bias")
  crit(sum(nul) > 0 && sum(m$sig_strain[nul]) <= 2 && sum(m$sig_parent_of_origin[nul]) <= 2,
       sprintf("null genes: strain false positives %d, parent-of-origin false positives %d, of %d (limit 2 each)",
               sum(m$sig_strain[nul]), sum(m$sig_parent_of_origin[nul]), sum(nul)))
  crit(sum(nul | exS | exP) > 0 && all(m$status[nul | exS | exP] %in% c("ok", "ok_at_bound")), "every planted and null gene was tested")
  su <- ck$snps_used
  dnc <- m$gene_id[m$class == "dense_null"]; dn <- su[su$gene_id %in% dnc, ]
  crit(length(dnc) > 0 && all(tapply(dn$exon_pos, dn$gene_id, function(p) length(p) < 2 || min(diff(sort(p))) >= THIN_BP)),
       sprintf("dense_null genes (%d): kept SNPs at least THIN_BP (%s) apart", length(dnc), THIN_BP))
  # direction strings of the gene table follow the signs of b0 and b1 (a swapped label would otherwise pass the sign criteria)
  expS <- ifelse(!g$sig_strain, "none", ifelse(g$b0_strain > 0, paste(SA, "higher"), paste(SB, "higher")))
  expP <- ifelse(!g$sig_parent_of_origin, "none", ifelse(g$b1_parent_of_origin > 0, "maternal higher", "paternal higher"))
  crit(any(g$sig_strain) && any(g$sig_parent_of_origin) && identical(g$direction_strain, expS) && identical(g$direction_parent_of_origin, expP),
       sprintf("direction strings match the signs of b0 and b1 in all %d genes (strain: %d mismatches, parent of origin: %d)",
               nrow(g), sum(g$direction_strain != expS), sum(g$direction_parent_of_origin != expP)))
  # summary_numbers_reciprocal.tsv against the gene table
  sf <- file.path(RES, "summary_numbers_reciprocal.tsv")
  sm <- tryCatch(read.delim(sf, stringsAsFactors = FALSE, quote = "",
                            colClasses = c(test = "character", direction_pos = "character", direction_neg = "character")),
                 error = function(e) NULL)
  des <- ck$design
  want <- function(test, sig, b, p, lp, ln) data.frame(test = test, genes_tested = sum(!is.na(p)), genes_not_tested = nrow(ck$not_tested) + sum(is.na(p)),
                                                     sig_genes = sum(sig), direction_pos = lp, n_pos = sum(sig & b > 0, na.rm = TRUE),
                                                     direction_neg = ln, n_neg = sum(sig & b < 0, na.rm = TRUE), phi_common = g$phi_common[1],
                                                     samples_AxB = sum(des$cross_direction == "AxB"), samples_BxA = sum(des$cross_direction == "BxA"),
                                                     stringsAsFactors = FALSE)
  w <- rbind(want("strain", g$sig_strain, g$b0_strain, g$p_strain, paste(SA, "higher"), paste(SB, "higher")),
             want("parent_of_origin", g$sig_parent_of_origin, g$b1_parent_of_origin, g$p_parent_of_origin, "maternal higher", "paternal higher"))
  same_sm <- !is.null(sm) && identical(names(sm), names(w)) && nrow(sm) == 2 &&
    all(vapply(names(w), function(k) if (is.numeric(w[[k]])) isTRUE(all.equal(as.numeric(sm[[k]]), w[[k]], tolerance = 1e-6))
                                     else identical(as.character(sm[[k]]), w[[k]]), logical(1)))
  crit(same_sm, "summary_numbers_reciprocal.tsv equals the counts, labels and design of the gene table")
  if (!same_sm && !is.null(sm)) { cat("summary file:\n"); print(sm); cat("expected from the gene table:\n"); print(w) }
}
if (what %in% c("differential_f1", "differential_outbred")) {
  THIN_BP <- as.numeric(val("THIN_BP")); SA <- val("STRAIN_A"); SB <- val("STRAIN_B"); RC <- val("REF_CONDITION")
  ck <- readRDS(file.path(RES, "ase_differential_checkpoint.rds"))
  crit(isTRUE(all.equal(ck$constants$THIN_BP, THIN_BP)) && identical(ck$constants$STRAIN_A, SA) && identical(ck$constants$STRAIN_B, SB) &&
         identical(ck$constants$REF_CONDITION, RC),
       sprintf("checkpoint constants match the values file (THIN_BP %s, STRAIN_A %s, STRAIN_B %s, REF_CONDITION %s)", THIN_BP, SA, SB, RC))
  CT <- paste("treat vs", RC)
  # summary_numbers_differential.tsv recomputed from the checkpoint tables (a unit listed in not_tested and with p NA counts once)
  sf <- file.path(RES, "summary_numbers_differential.tsv")
  sm <- tryCatch(read.delim(sf, stringsAsFactors = FALSE, quote = "", colClasses = c(contrast = "character", mode = "character",
                                                                                      level = "character", dispersion = "character")),
                 error = function(e) NULL)
  want_row <- function(t, lev, idcol, pcol, up, dn, disp) {
    nt <- ck$not_tested$id[ck$not_tested$contrast == CT & ck$not_tested$level == lev]
    data.frame(contrast = CT, mode = ck$constants$MODE, level = lev, tested = sum(!is.na(t[[pcol]])),
               not_tested = length(unique(c(nt, t[[idcol]][is.na(t[[pcol]])]))), sig = sum(t$sig), n_up = up, n_down = dn,
               dispersion = disp, stringsAsFactors = FALSE)
  }
  cmp_summary <- function(w) {
    ok <- !is.null(sm) && identical(names(sm), names(w)) && nrow(sm) == nrow(w) &&
      all(vapply(names(w), function(k) if (is.numeric(w[[k]]) || all(is.na(w[[k]]))) isTRUE(all.equal(as.numeric(sm[[k]]), as.numeric(w[[k]])))
                                       else identical(as.character(sm[[k]]), w[[k]]), logical(1)))
    crit(ok, "summary_numbers_differential.tsv equals the counts, directions and dispersion of the checkpoint tables (gene and SNP rows)")
    if (!ok && !is.null(sm)) { cat("summary file:\n"); print(sm); cat("expected:\n"); print(w) }
  }
}
if (what == "differential_f1") {
  g <- ck$gene[ck$gene$contrast == CT, ]; sn <- ck$snp[ck$snp$contrast == CT, ]
  m <- merge(tg, g, by = "gene_id", all.x = TRUE); m$sig[is.na(m$sig)] <- FALSE
  print(m[m$class != "null", c("gene_id", "class", "expect_diff", "frac_A_ref", "frac_A_test", "delta_frac", "padj", "direction", "status")], row.names = FALSE)
  # every criterion needs at least one truth gene of its class: a relabelled truth table must FAIL, not PASS as "0 of 0"
  ex <- m$expect_diff %in% c("up", "down")
  okD <- m$sig & ((m$expect_diff == "up" & m$delta_frac > 0) | (m$expect_diff == "down" & m$delta_frac < 0))
  crit(sum(m$expect_diff == "up") > 0 && sum(m$expect_diff == "down") > 0 && all(okD[ex]),
       sprintf("F1 differential genes recovered with the correct sign: %d of %d", sum(okD[ex]), sum(ex)))
  sp <- m$class %in% c("strain_A_high", "strain_B_high", "maternal", "paternal", "strain_and_maternal")
  crit(sum(sp) > 0 && sum(m$sig[sp]) == 0, sprintf("no differential call among strain / parent-of-origin genes (no condition effect): %d of %d", sum(m$sig[sp]), sum(sp)))
  nul <- m$class %in% c("null", "dense_null", "bias")
  crit(sum(nul) > 0 && sum(m$sig[nul]) <= 2, sprintf("null genes: %d false positives of %d (limit 2)", sum(m$sig[nul]), sum(nul)))
  crit(all(m$status[nul | ex | sp] %in% c("ok", "ok_at_bound")), "every planted and null gene was tested")
  crit(isTRUE(ck$design$cross_direction_term[1]), "the cross-direction term was used (both directions present)")
  lab <- function(t) ifelse(!t$sig, "none", ifelse(t$b_condition > 0, paste0(SA, " fraction higher in treat"), paste0(SA, " fraction lower in treat")))
  crit(any(g$sig) && any(sn$sig) && identical(g$direction, lab(g)) && identical(sn$direction, lab(sn)) &&
         all(sign(g$delta_frac[g$sig]) == sign(g$b_condition[g$sig])),
       sprintf("direction strings match the sign of b_condition (genes: %d mismatches, SNPs: %d)", sum(g$direction != lab(g)), sum(sn$direction != lab(sn))))
  cmp_summary(rbind(want_row(g, "gene", "gene_id", "p_condition", sum(g$sig & g$delta_frac > 0), sum(g$sig & g$delta_frac < 0), sprintf("phi_common %.4f", g$phi_common[1])),
                    want_row(sn, "SNP", "SNP", "p_condition", sum(sn$sig & sn$delta_frac > 0), sum(sn$sig & sn$delta_frac < 0), sprintf("phi_common %.4f", sn$phi_common[1]))))
}
if (what == "differential_outbred") {
  g <- ck$gene[ck$gene$contrast == CT, ]; sn <- ck$snp[ck$snp$contrast == CT, ]
  # truth columns only: the truth's n_snps (SNPs simulated) would collide with the gene table's n_snps (SNPs tested)
  m <- merge(tg[, c("gene_id", "class", "expect_diff")], g, by = "gene_id", all.x = TRUE); m$sig[is.na(m$sig)] <- FALSE
  print(m[m$class != "null", c("gene_id", "class", "expect_diff", "n_snps", "acat_p", "padj", "max_abs_delta", "sig")], row.names = FALSE)
  ex <- m$class %in% c("diff_phase", "diff_consistent")
  crit(sum(ex) > 0 && sum(m$sig[ex]) >= 4, sprintf("outbred differential genes detected (gene level, unphased): %d of %d (need >= 4)", sum(m$sig[ex]), sum(ex)))
  bi <- m$class == "base_imbalanced"
  crit(sum(bi) > 0 && sum(m$sig[bi]) == 0, sprintf("no call among imbalanced but unchanged genes: %d of %d", sum(m$sig[bi]), sum(bi)))
  nul <- m$class %in% c("null", "bias")
  crit(sum(nul) > 0 && sum(m$sig[nul]) <= 2, sprintf("null genes: %d false positives of %d (limit 2)", sum(m$sig[nul]), sum(nul)))
  s <- merge(sn[sn$sig, ], ck$snp_gene, by = "SNP")
  s <- merge(s, tg[, c("gene_id", "class", "expect_diff")], by = "gene_id")
  cs <- s$class == "diff_consistent"
  crit(any(cs) && all(s$direction[cs] == "REF lower in treat (all tested individuals)") && all(s$expect_diff[cs] == "down"),
       sprintf("consistent genes: significant SNPs say 'REF lower in treat (all tested individuals)' (%d of %d)",
               sum(s$direction[cs] == "REF lower in treat (all tested individuals)"), sum(cs)))
  crit(any(s$direction[s$class == "diff_phase"] == "mixed (phase differs)"),
       sprintf("phase-heterogeneous genes: at least one SNP labelled 'mixed (phase differs)' (%d of %d)",
               sum(s$direction[s$class == "diff_phase"] == "mixed (phase differs)"), sum(s$class == "diff_phase")))
  lab <- ifelse(!sn$sig, "none", ifelse(sn$n_down == 0, "REF higher in treat (all tested individuals)",
                                        ifelse(sn$n_up == 0, "REF lower in treat (all tested individuals)", "mixed (phase differs)")))
  crit(any(sn$sig) && identical(sn$direction, lab), sprintf("SNP direction strings match n_up / n_down (%d mismatches)", sum(sn$direction != lab)))
  ph <- ck$dispersion$phi_pair[ck$dispersion$contrast == CT]
  dsp <- sprintf("phi_pair %.4f-%.4f", min(ph, na.rm = TRUE), max(ph, na.rm = TRUE))
  cmp_summary(rbind(want_row(g, "gene", "gene_id", "acat_p", NA_integer_, NA_integer_, dsp),
                    want_row(sn, "SNP", "SNP", "p", sum(sn$sig & sn$n_down == 0), sum(sn$sig & sn$n_up == 0), dsp)))
}
if (!what %in% c("reciprocal", "differential_f1", "differential_outbred")) stop("unknown evaluation: ", what, call. = FALSE)
quit(status = fail)
