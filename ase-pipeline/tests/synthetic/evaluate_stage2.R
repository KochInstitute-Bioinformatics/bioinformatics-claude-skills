# Usage: Rscript evaluate_stage2.R <reciprocal|differential_f1|differential_outbred> <RESULTS_DIR> <truth dir (f1_recip or outbred_diff)>
args <- commandArgs(trailingOnly = TRUE); what <- args[1]; RES <- args[2]; TRUTH <- args[3]
fail <- 0
crit <- function(cond, msg) { cat(if (isTRUE(cond)) "PASS" else "FAIL", msg, "\n"); if (!isTRUE(cond)) fail <<- 1 }
tg <- read.delim(file.path(TRUTH, "truth_genes.tsv"), stringsAsFactors = FALSE)
if (what == "reciprocal") {
  g <- readRDS(file.path(RES, "ase_reciprocal_checkpoint.rds"))$gene
  m <- merge(tg, g, by = "gene_id", all.x = TRUE)
  for (k in c("sig_strain", "sig_parent_of_origin")) m[[k]][is.na(m[[k]])] <- FALSE
  print(m[m$class != "null", c("gene_id", "class", "expect_strain", "expect_poo", "frac_A", "maternal_frac", "padj_strain",
                               "padj_parent_of_origin", "direction_strain", "direction_parent_of_origin", "status")], row.names = FALSE)
  exS <- m$expect_strain %in% c("A_higher", "B_higher")
  okS <- m$sig_strain & ((m$expect_strain == "A_higher" & m$b0_strain > 0) | (m$expect_strain == "B_higher" & m$b0_strain < 0))
  crit(all(okS[exS]), sprintf("strain effect recovered with the correct sign: %d of %d", sum(okS[exS]), sum(exS)))
  exP <- m$expect_poo %in% c("maternal", "paternal")
  okP <- m$sig_parent_of_origin & ((m$expect_poo == "maternal" & m$b1_parent_of_origin > 0) | (m$expect_poo == "paternal" & m$b1_parent_of_origin < 0))
  crit(all(okP[exP]), sprintf("parent-of-origin effect recovered with the correct sign: %d of %d", sum(okP[exP]), sum(exP)))
  so <- m$class %in% c("strain_A_high", "strain_B_high"); po <- m$class %in% c("maternal", "paternal")
  crit(sum(m$sig_parent_of_origin[so]) == 0, sprintf("no parent-of-origin call among strain-only genes: %d of %d", sum(m$sig_parent_of_origin[so]), sum(so)))
  crit(sum(m$sig_strain[po]) == 0, sprintf("no strain call among parent-of-origin-only genes: %d of %d", sum(m$sig_strain[po]), sum(po)))
  nul <- m$class %in% c("null", "dense_null", "bias")
  crit(sum(m$sig_strain[nul]) <= 2 && sum(m$sig_parent_of_origin[nul]) <= 2,
       sprintf("null genes: strain false positives %d, parent-of-origin false positives %d, of %d (limit 2 each)",
               sum(m$sig_strain[nul]), sum(m$sig_parent_of_origin[nul]), sum(nul)))
  crit(all(m$status[nul | exS | exP] %in% c("ok", "ok_at_bound")), "every planted and null gene was tested")
  su <- readRDS(file.path(RES, "ase_reciprocal_checkpoint.rds"))$snps_used
  dn <- su[su$gene_id %in% m$gene_id[m$class == "dense_null"], ]
  crit(all(tapply(dn$exon_pos, dn$gene_id, function(p) length(p) < 2 || min(diff(sort(p))) >= 500)),
       "dense_null genes: kept SNPs at least THIN_BP (500) apart")
}
quit(status = fail)
