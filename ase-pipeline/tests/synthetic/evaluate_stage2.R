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
quit(status = fail)
