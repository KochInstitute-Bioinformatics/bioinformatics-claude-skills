# Usage (inside the R image, on a compute node only; tests/r/run_r_tests.sh does this through sbatch):
#   Rscript test_enrich_summary.R <skill.md> <tests/fixtures dir>
# Extracts the Rmd blocks under **File `06_enrichment.Rmd`.** and **File `07_summary.Rmd`.** of the skill and tests them.
# Rmd 06: renders on SYNTHETIC gene lists (human and mouse Ensembl identifiers taken from org.Hs.eg.db and org.Mm.eg.db inside this
# test: the genes of one GO term, a universe of 2000 testable genes) and compares the tables with a hypergeometric test written here
# in base R (phyper), the group sizes, the universe size, the skip rule, the organism and database messages; renders on the real
# fixtures and compares the group sizes and the universe with the cross-tool report. Rmd 07: renders on a scratch report directory
# (two reports, one figure, one table) and on a complete report directory of the fixtures, and compares what it lists and counts with
# the files and the counts tables. The SYNTHETIC data test the counting and listing logic only, never a sign of a real tool.
# Prints one "ok" or "FAIL:" line per check, then "R TESTS PASS (<n>)" or "R TESTS FAIL (<failed>/<n>)" (exit 1).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: test_enrich_summary.R <skill.md> <fixtures dir>")
skill <- args[1]; fx <- args[2]; rt <- file.path(fx, "real_testdata")
n_ok <- 0; n_fail <- 0
check <- function(name, cond) {
  r <- tryCatch({ do.call(stopifnot, stats::setNames(list(isTRUE(all(cond))), name)); TRUE },
                error = function(e) { cat("FAIL:", conditionMessage(e), "\n"); FALSE })
  if (r) { cat("ok  ", name, "\n"); n_ok <<- n_ok + 1 } else n_fail <<- n_fail + 1
  invisible(r)
}
finish <- function() {
  if (n_fail == 0) { cat(sprintf("R TESTS PASS (%d)\n", n_ok)); quit(status = 0) }
  cat(sprintf("R TESTS FAIL (%d/%d)\n", n_fail, n_ok + n_fail)); quit(status = 1)
}
work <- tempfile("rsr_enrich_"); dir.create(work)
txt <- readLines(skill, warn = FALSE)

get_rmd <- function(name) {
  h <- which(txt == sprintf("**File `%s`.**", name))
  if (!check(sprintf("the skill has the heading **File `%s`.** once", name), length(h) == 1)) finish()
  o <- h + which(txt[(h + 1):length(txt)] == "````rmd")[1]
  e <- o + which(txt[(o + 1):length(txt)] == "````")[1]
  if (!check(sprintf("a ````rmd block follows the heading of %s", name), !is.na(o) && !is.na(e) && e > o + 1)) finish()
  txt[(o + 1):(e - 1)]
}
functions_env <- function(rmd, name) {
  fb <- grep("^```\\{r functions", rmd)
  if (!check(sprintf("%s has one functions chunk", name), length(fb) == 1)) finish()
  fe <- fb + which(rmd[(fb + 1):length(rmd)] == "```")[1]
  e <- new.env(parent = globalenv())
  err <- tryCatch({ eval(parse(text = rmd[(fb + 1):(fe - 1)]), envir = e); "" }, error = function(x) conditionMessage(x))
  if (!check(sprintf("the functions chunk of %s evaluates (%s)", name, err), err == "")) finish()
  e
}
rmd6 <- get_rmd("06_enrichment.Rmd"); rmd7 <- get_rmd("07_summary.Rmd")
# the code chunks are scanned before any of them is evaluated: a network call, download or installation is never run by this test
for (nm in c("06", "07")) {
  r <- if (nm == "06") rmd6 else rmd7
  chunks <- local({ inb <- FALSE; keep <- character(0); for (l in r) { if (grepl("^```\\{r", l)) { inb <- TRUE; next }; if (inb && l == "```") { inb <- FALSE; next }; if (inb) keep <- c(keep, l) }; keep })
  if (!check(sprintf("%s has no network call, download or installation in its code (KEGGREST, download, url(, install.packages)", nm),
        !any(grepl("KEGGREST|download|url\\(|install\\.packages|curl|httr|RCurl|https?:|ftp:|BiocManager|install_|wget|socketConn|system2?\\(|source\\(", chunks)) && !any(grepl("^```\\{[^r]", r)))) finish()
}
E6 <- functions_env(rmd6, "06_enrichment.Rmd"); E7 <- functions_env(rmd7, "07_summary.Rmd")
need_fun <- function(e, nm) all(vapply(nm, exists, logical(1), envir = e, inherits = FALSE))
if (!check("06 defines its functions", need_fun(E6, c("collect_calls", "group_genes", "to_entrez", "run_go", "requested_dbs", "MIN_GENES", "ENRICH_CUTOFF", "ORGDB")))) finish()
if (!check("07 defines its functions", need_fun(E7, c("REPORTS", "read_counts", "html_warnings", "list_files")))) finish()
txt_of <- function(rmd) paste(rmd, collapse = "\n")
SIGN <- "Every effect is stated as **treatment minus control**"
for (nm in c("06", "07")) {
  r <- if (nm == "06") rmd6 else rmd7
  check(sprintf("%s states the sign rule before its first function", nm), regexpr(SIGN, txt_of(r), fixed = TRUE) > 0 &&
        regexpr(SIGN, txt_of(r), fixed = TRUE) < regexpr("```{r functions", txt_of(r), fixed = TRUE))
  used <- sub("params\\$", "", unique(unlist(regmatches(txt_of(r), gregexpr("params\\$[A-Za-z_]+", txt_of(r))))))
  check(sprintf("%s uses only the 18 parameters of the contract (%s)", nm, paste(setdiff(used, c("outdir", "report_dir", "contrasts", "samplesheet", "modules", "fdr", "dpsi", "rmats_type", "min_reads", "events", "cross", "topgene", "enrich", "organism", "enrich_dbs", "gtf", "paired_rmats", "paired_suppa")), collapse = ",")),
        all(used %in% c("outdir", "report_dir", "contrasts", "samplesheet", "modules", "fdr", "dpsi", "rmats_type", "min_reads", "events", "cross", "topgene", "enrich", "organism", "enrich_dbs", "gtf", "paired_rmats", "paired_suppa")))
}
check("06 says that the universe is the testable genes and that a group needs at least 10 genes", grepl("testable", txt_of(rmd6), fixed = TRUE) && grepl("at least 10 genes", txt_of(rmd6), fixed = TRUE))
check("the skill constants: minimum 10 genes, cutoff 0.05", identical(E6$MIN_GENES, 10) && isTRUE(all.equal(E6$ENRICH_CUTOFF, 0.05)))

suppressPackageStartupMessages({ library(org.Hs.eg.db); library(org.Mm.eg.db); library(GO.db) })
# --- SYNTHETIC data: gene identifiers drawn from the annotation packages (no identifier is typed by hand)
one_to_one <- function(egens) { # Bimap entrez -> ensembl
  l <- as.list(egens); l <- l[!vapply(l, function(x) all(is.na(x)), logical(1))]
  d <- data.frame(entrez = rep(names(l), lengths(l)), ensembl = unlist(l, use.names = FALSE), stringsAsFactors = FALSE)
  d[!(duplicated(d$entrez) | duplicated(d$entrez, fromLast = TRUE)) & !(duplicated(d$ensembl) | duplicated(d$ensembl, fromLast = TRUE)), , drop = FALSE]
}
SPECIES <- list(
  human = list(pkg = "org.Hs.eg.db", map = one_to_one(org.Hs.egENSEMBL), go2eg = as.list(org.Hs.egGO2ALLEGS), fake = "ENSG000999999%02d"),
  mouse = list(pkg = "org.Mm.eg.db", map = one_to_one(org.Mm.egENSEMBL), go2eg = as.list(org.Mm.egGO2ALLEGS), fake = "ENSMUSG000999999%02d"))
TERM_UP <- "GO:0006412"; TERM_DOWN <- "GO:0006955"   # translation; immune response
term_ens <- function(sp, term) sp$map$ensembl[sp$map$entrez %in% unique(sp$go2eg[[term]])]
C_SYN <- "SYN"
# a SYNTHETIC rMATS table (one event per gene, a second non-significant event for three genes): `up` and `down` are significant genes,
# `bg` testable and not significant, `na` not testable (padj NA)
synth_table <- function(up, down, bg, na) {
  mk <- function(g, effect, padj, dir) data.frame(contrast = C_SYN, tool = "rmats", gene_id = g, gene_symbol = NA_character_, feature_type = "SE", feature_id = NA_character_,
    effect = effect, effect_scale = "dPSI", padj = padj, direction = dir, n_support = 100, low_support = FALSE, stringsAsFactors = FALSE)
  d <- rbind(mk(up, 0.5, 0.001, "treatment_up"), mk(down, -0.5, 0.001, "control_up"), mk(bg, 0.01, 0.9, "treatment_up"), mk(na, NA_real_, NA_real_, "none"),
             mk(up[1:3], 0.01, 0.5, "treatment_up"))
  d$feature_id <- paste0("SE:", seq_len(nrow(d)))
  d
}
synth <- function(species, n_up, n_down, seed = 1) {
  sp <- SPECIES[[species]]; set.seed(seed)
  fake <- sprintf(sp$fake, 1:3)
  up <- sample(term_ens(sp, TERM_UP), n_up)
  down <- sample(setdiff(term_ens(sp, TERM_DOWN), up), n_down)
  rest <- setdiff(sp$map$ensembl, c(up, down))
  nbg <- 2000 - length(up) - length(down) - 3
  bg <- sample(rest, nbg); na <- sample(setdiff(rest, bg), 100)
  list(species = species, up = c(up, fake[1]), down = down, bg = c(bg, fake[2:3]), na = na)
}
base_params <- list(outdir = rt, report_dir = "", contrasts = "", samplesheet = file.path(rt, "261003_a_samplesheet.csv"), modules = "rmats", rmats_type = "JC",
                    fdr = 0.05, dpsi = 0.1, min_reads = 10, events = TRUE, cross = TRUE, topgene = TRUE, enrich = TRUE,
                    organism = "human", gtf = "", enrich_dbs = "GO", paired_rmats = FALSE, paired_suppa = TRUE)
render <- function(rd, rmd_name, ...) {
  pr <- utils::modifyList(base_params, c(list(report_dir = rd), list(...)))
  of <- if (rmd_name == "07_summary.Rmd") "report_summary.html" else NULL
  tryCatch({ rmarkdown::render(file.path(rd, rmd_name), output_file = of, output_dir = rd, knit_root_dir = rd, params = pr, quiet = TRUE, envir = new.env()); TRUE },
           error = function(e) { cat("render error:", conditionMessage(e), "\n"); FALSE })
}
new_dir <- function(nm) { d <- file.path(work, nm); dir.create(d)
  writeLines(rmd6, file.path(d, "06_enrichment.Rmd")); writeLines(rmd7, file.path(d, "07_summary.Rmd")); d }
get_other <- function(n) { h <- which(txt == sprintf("**File `%s`.**", n)); o <- h + which(txt[(h + 1):length(txt)] == "````rmd")[1]; e <- o + which(txt[(o + 1):length(txt)] == "````")[1]; txt[(o + 1):(e - 1)] }
RMD_OTHER <- stats::setNames(lapply(c("01_import.Rmd", "02_rmats.Rmd", "03_suppa.Rmd", "04_usage.Rmd", "05_cross_tool.Rmd"), get_other), c("01_import.Rmd", "02_rmats.Rmd", "03_suppa.Rmd", "04_usage.Rmd", "05_cross_tool.Rmd"))
new_dir_all <- function(nm) { d <- file.path(work, nm); dir.create(d); for (r in names(RMD_OTHER)) writeLines(RMD_OTHER[[r]], file.path(d, r))
  writeLines(rmd6, file.path(d, "06_enrichment.Rmd")); writeLines(rmd7, file.path(d, "07_summary.Rmd")); d }
synth_dir <- function(nm, s) {
  d <- new_dir(nm); dir.create(file.path(d, "rds"))
  saveRDS(synth_table(s$up, s$down, s$bg, s$na), file.path(d, "rds", paste0(C_SYN, ".rmats.rds")))
  writeLines(c("contrast,treatment,control", paste0(C_SYN, ",T1,C1")), file.path(d, "contrasts.csv"))
  d
}
html_of <- function(rd, nm) gsub("[[:space:]]+", " ", paste(readLines(file.path(rd, nm), warn = FALSE), collapse = " "))
has <- function(html, s) grepl(s, html, fixed = TRUE)
tab_file <- function(rd, ...) file.path(rd, "tables", paste0(C_SYN, ".enrich.", paste(..., sep = "."), ".csv"))
rd_csv <- function(f) utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE, colClasses = "character", na.strings = character(0))
summ_of <- function(rd, C = C_SYN) { s <- rd_csv(file.path(rd, "tables", paste0(C, ".enrich.summary.csv"))); s }
srow <- function(s, tool, group) s[s$tool == tool & s$group == group, , drop = FALSE]
dpath <- function(rd) file.path(rd, "contrasts.csv")

# independent hypergeometric p-value of a GO term (base R): the genes of the term come from the annotation package's GO2ALLEGS
# (the term and all its descendants), restricted to the universe and to the genes annotated to the ontology
indep_p <- function(species, ont, term, input_ens, universe_ens) {
  sp <- SPECIES[[species]]
  ent <- function(e) sort(unique(sp$map$entrez[match(e, sp$map$ensembl)][!is.na(match(e, sp$map$ensembl))]))
  inp <- ent(input_ens); uni <- ent(universe_ens)
  ont_terms <- names(sp$go2eg)[AnnotationDbi::Ontology(names(sp$go2eg)) %in% ont]
  ann <- unique(unlist(sp$go2eg[ont_terms], use.names = FALSE)); ann_uni <- intersect(uni, ann)
  tg <- intersect(unique(sp$go2eg[[term]]), ann_uni)
  n <- length(intersect(inp, ann_uni)); k <- length(intersect(tg, inp)); M <- length(tg); N <- length(ann_uni)
  stats::phyper(k - 1, M, N - M, n, lower.tail = FALSE)
}
# relative comparison of p-values (they range from 1e-35 to 1): all.equal would fall back to the absolute difference for values below its tolerance
rel_eq <- function(a, b, tol = 1e-6) length(a) == length(b) && length(a) > 0 && all(is.finite(a)) && all(is.finite(b)) && all(a > 0) && all(b > 0) && max(abs(a / b - 1)) < tol
check("rel_eq separates p-values that all.equal would call equal (1e-20 against 5e-20) and accepts equal ones", !rel_eq(1e-20, 5e-20) && isTRUE(all.equal(1e-20, 5e-20, tolerance = 1e-6)) && rel_eq(c(1e-20, 0.3), c(1e-20 * (1 + 1e-9), 0.3)))
# SYNTHETIC human: groups treatment_up (30 genes, one unmapped), control_up (12), all (42)
S_h <- synth("human", 29, 12)
rdA <- synth_dir("synth_human", S_h)
okA <- render(rdA, "06_enrichment.Rmd", contrasts = dpath(rdA))
check("06 renders on the SYNTHETIC human gene list", okA)
if (!okA) finish()
hA <- html_of(rdA, "06_enrichment.html")
sA <- summ_of(rdA)
check("SYNTHETIC human: the group sizes are 42 (all), 30 (treatment_up), 12 (control_up)",
      identical(as.integer(srow(sA, "rmats", "all")$genes), 42L) && identical(as.integer(srow(sA, "rmats", "treatment_up")$genes), 30L) && identical(as.integer(srow(sA, "rmats", "control_up")$genes), 12L))
check("SYNTHETIC human: the unmapped identifiers are counted (1 of the 30 and of the 42, none of the 12) and reported",
      identical(as.integer(srow(sA, "rmats", "treatment_up")$genes_unmapped), 1L) && identical(as.integer(srow(sA, "rmats", "treatment_up")$genes_mapped), 29L) &&
      identical(as.integer(srow(sA, "rmats", "all")$genes_unmapped), 1L) && identical(as.integer(srow(sA, "rmats", "control_up")$genes_unmapped), 0L) && has(hA, "unmapped"))
tbl <- readRDS(file.path(rdA, "rds", paste0(C_SYN, ".rmats.rds")))
n_testable <- length(unique(tbl$gene_id[!is.na(tbl$padj)]))
check("the SYNTHETIC table has 2000 testable genes (padj present), more rows than genes, and 100 more genes that are not testable", n_testable == 2000 && nrow(tbl) > n_testable && length(unique(tbl$gene_id)) == 2100)
check("universe size equals the number of genes with a non-NA padj (2000), for every group, and 3 of them are unmapped (the identifiers that are not in the annotation)",
      all(as.integer(sA$universe_genes) == n_testable) && all(as.integer(sA$universe_unmapped) == 3L) && all(as.integer(sA$universe_mapped) == 1997L))
check("the status of the three groups is ran", all(sA$status == "ran"))
for (g in c("all", "treatment_up", "control_up")) for (o in c("go_bp", "go_mf", "go_cc"))
  check(sprintf("SYNTHETIC human: table %s.%s written", g, o), file.exists(tab_file(rdA, "rmats", g, o)))
bp <- rd_csv(tab_file(rdA, "rmats", "treatment_up", "go_bp"))
check("treatment_up (translation genes): at least one GO BP term is significant (p.adjust below 0.05)", nrow(bp) >= 1 && all(as.numeric(bp$p.adjust) < 0.05))
anc <- c(TERM_UP, unlist(as.list(GOBPANCESTOR[TERM_UP]), use.names = FALSE))
check("treatment_up: the terms contain the queried term GO:0006412 or a close parent", any(bp$ID %in% anc))
uni_ens <- unique(tbl$gene_id[!is.na(tbl$padj)])
p_tab <- as.numeric(bp$pvalue); names(p_tab) <- bp$ID
chk_ids <- unique(c(intersect(TERM_UP, bp$ID), head(bp$ID[order(p_tab)], 5)))
p_ind <- vapply(chk_ids, function(t) indep_p("human", "BP", t, S_h$up, uni_ens), numeric(1))
check(sprintf("treatment_up: the p-values of the queried term and of the 5 smallest equal an independent phyper to 1e-6 (relative; %s)", paste(signif(p_ind, 3), collapse = " ")),
      length(chk_ids) >= 1 && rel_eq(unname(p_tab[chk_ids]), unname(p_ind)))
check("treatment_up: the queried term GO:0006412 is in the table and its p-value equals phyper", TERM_UP %in% bp$ID && rel_eq(p_tab[[TERM_UP]], indep_p("human", "BP", TERM_UP, S_h$up, uni_ens)))
check("treatment_up: every term of the table has the BH adjusted p-value of the table at or above its p-value", all(as.numeric(bp$p.adjust) >= as.numeric(bp$pvalue)))
cu <- rd_csv(tab_file(rdA, "rmats", "control_up", "go_bp"))
check("control_up (12 immune-response genes) is a separate enrichment: its table does not equal the treatment_up table", !identical(cu$ID, bp$ID))
alls <- rd_csv(tab_file(rdA, "rmats", "all", "go_bp"))
p_all <- as.numeric(alls$pvalue); names(p_all) <- alls$ID
all_ens <- c(S_h$up, S_h$down)
chk_all <- head(alls$ID[order(p_all)], 3)
check("the group all (42 genes) matches the independent phyper for its 3 smallest p-values", length(chk_all) >= 1 && rel_eq(unname(p_all[chk_all]), unname(vapply(chk_all, function(t) indep_p("human", "BP", t, all_ens, uni_ens), numeric(1)))))
check("the dotplot figure of the contrast is written", file.exists(file.path(rdA, "figures", paste0(C_SYN, ".enrich.dotplot.png"))) && file.size(file.path(rdA, "figures", paste0(C_SYN, ".enrich.dotplot.png"))) > 0)
check("the default render prints none of the skip, organism and rebuild messages", !has(hA, "skipped: fewer than") && !has(hA, "organism not supported") && !has(hA, "require an image rebuild"))
check("the report states the universe (testable genes) and the thresholds", has(hA, "testable") && has(hA, "treatment minus control"))
check("the unit functions: group_genes splits by the direction of the top feature, to_entrez counts the unmapped",
      { cl <- list(genes = c("a", "b", "c"), tested = c("a", "b", "c", "d"), dir = data.frame(gene_id = c("a", "b"), direction = c("treatment_up", "control_up"), stringsAsFactors = FALSE))
        g <- E6$group_genes(cl); identical(g$treatment_up, "a") && identical(g$control_up, "b") && identical(sort(g$all), c("a", "b", "c")) } &&
      { m <- E6$to_entrez(c(S_h$up[1:3], "ENSG00099999999"), "org.Hs.eg.db"); m$n_ids == 4 && m$n_mapped == 3 && m$n_unmapped == 1 && length(m$entrez) >= 3 })
check("requested_dbs reads GO, KEGG and Reactome from a comma or space separated value, any case",
      identical(E6$requested_dbs("GO,KEGG"), c("go", "kegg")) && identical(E6$requested_dbs("GO KEGG Reactome"), c("go", "kegg", "reactome")) && identical(E6$requested_dbs("go"), "go"))

# gene sets of Rmd 06: a flagged (low_support) rMATS event is not called, and versioned Ensembl identifiers equal their unversioned forms
S_v <- S_h
tv <- synth_table(S_v$up, S_v$down, S_v$bg, S_v$na)
flag_gene <- S_v$na[1]   # a real gene, testable here, with a significant event that is flagged low_support
tv$padj[tv$gene_id == flag_gene] <- 0.001; tv$effect[tv$gene_id == flag_gene] <- 0.5; tv$direction[tv$gene_id == flag_gene] <- "treatment_up"; tv$low_support[tv$gene_id == flag_gene] <- TRUE
ver <- S_v$up[2:6]; tv$gene_id[tv$gene_id %in% ver] <- paste0(tv$gene_id[tv$gene_id %in% ver], ".3")
cl_v <- E6$calls_from_events(tv, 0.05, 0.1)
check("calls_from_events (06): the flagged gene is not called, versions are removed (the 5 versioned genes are called as unversioned identifiers)",
      !(flag_gene %in% cl_v$genes) && flag_gene %in% cl_v$tested && all(ver %in% cl_v$genes) && !any(grepl("\\.3$", cl_v$genes)) && length(cl_v$genes) == 42)
check("calls_from_gene_table (06): versions are removed on the gene-level table too",
      identical(E6$calls_from_gene_table(data.frame(gene_id = c("ENSG00000000001.3", "ENSG00000000002"), padj = c(0.01, 0.5), stringsAsFactors = FALSE),
                data.frame(gene_id = "ENSG00000000001.3", feature_id = "E1", padj = 0.01, effect = 1, direction = "treatment_up", stringsAsFactors = FALSE), 0.05)$genes, "ENSG00000000001"))
rdV <- new_dir("synth_versioned"); dir.create(file.path(rdV, "rds")); saveRDS(tv, file.path(rdV, "rds", paste0(C_SYN, ".rmats.rds")))
writeLines(c("contrast,treatment,control", paste0(C_SYN, ",T1,C1")), dpath(rdV))
okV <- render(rdV, "06_enrichment.Rmd", contrasts = dpath(rdV)); sV <- if (okV) summ_of(rdV) else NULL
check("06 renders with a flagged significant event and versioned identifiers: group sizes 42/30/12 (the flagged gene is not counted), the versioned genes map (1 unmapped in the 30), the flagged gene is testable (universe 2001)",
      okV && identical(as.integer(srow(sV, "rmats", "all")$genes), 42L) && identical(as.integer(srow(sV, "rmats", "treatment_up")$genes), 30L) && identical(as.integer(srow(sV, "rmats", "control_up")$genes), 12L) &&
      identical(as.integer(srow(sV, "rmats", "treatment_up")$genes_unmapped), 1L) && all(as.integer(sV$universe_genes) == 2001L))
# B: control_up has 7 genes: skipped, the other groups run
S_b <- synth("human", 29, 7)
rdB <- synth_dir("synth_human_7", S_b)
okB <- render(rdB, "06_enrichment.Rmd", contrasts = dpath(rdB))
check("06 renders when one direction group has fewer than 10 genes", okB)
hB <- html_of(rdB, "06_enrichment.html"); sB <- summ_of(rdB)
check("the group of 7 genes is skipped with the message 'skipped: fewer than 10 genes (n=7)' and no table is written for it",
      has(hB, "skipped: fewer than 10 genes (n=7)") && identical(srow(sB, "rmats", "control_up")$status, "skipped") &&
      !any(file.exists(tab_file(rdB, "rmats", "control_up", c("go_bp", "go_mf", "go_cc")))))
check("the groups of 30 and 37 genes still run and write tables", all(file.exists(tab_file(rdB, "rmats", "treatment_up", "go_bp"), tab_file(rdB, "rmats", "all", "go_bp"))) && identical(as.integer(srow(sB, "rmats", "all")$genes), 37L))
tbB <- readRDS(file.path(rdB, "rds", paste0(C_SYN, ".rmats.rds")))
check("the universe of the second list is again the testable genes", all(as.integer(sB$universe_genes) == length(unique(tbB$gene_id[!is.na(tbB$padj)]))))
# C: all groups below 10 genes
S_c <- synth("human", 3, 5)   # 3 mapped genes and 1 unmapped identifier up, 5 down
rdC <- synth_dir("synth_human_few", S_c)
okC <- render(rdC, "06_enrichment.Rmd", contrasts = dpath(rdC))
hC <- html_of(rdC, "06_enrichment.html"); sC <- summ_of(rdC)
check("06 renders when every group has fewer than 10 genes; all skipped with their n (4, 5, 9) and no GO table or dotplot is written",
      okC && has(hC, "skipped: fewer than 10 genes (n=9)") && has(hC, "skipped: fewer than 10 genes (n=4)") && has(hC, "skipped: fewer than 10 genes (n=5)") &&
      all(sC$status == "skipped") && length(list.files(file.path(rdC, "tables"), pattern = "\\.go_")) == 0 && !file.exists(file.path(rdC, "figures", paste0(C_SYN, ".enrich.dotplot.png"))) && has(hC, "no dotplot"))
check("the skipped groups still report the universe size (testable genes)", all(as.integer(sC$universe_genes) == length(unique(readRDS(file.path(rdC, "rds", paste0(C_SYN, ".rmats.rds")))$gene_id[!is.na(readRDS(file.path(rdC, "rds", paste0(C_SYN, ".rmats.rds")))$padj)]))))

# D: mouse
S_m <- synth("mouse", 29, 12)
rdM <- synth_dir("synth_mouse", S_m)
okM <- render(rdM, "06_enrichment.Rmd", contrasts = dpath(rdM), organism = "mouse")
check("06 renders on the SYNTHETIC mouse list (org.Mm.eg.db)", okM)
if (okM) {
  bpm <- rd_csv(tab_file(rdM, "rmats", "treatment_up", "go_bp")); sM <- summ_of(rdM); tbM <- readRDS(file.path(rdM, "rds", paste0(C_SYN, ".rmats.rds")))
  pm <- as.numeric(bpm$pvalue); names(pm) <- bpm$ID; um <- unique(tbM$gene_id[!is.na(tbM$padj)])
  ids_m <- head(bpm$ID[order(pm)], 3)
  check("mouse: at least one GO BP term, the group sizes are 30 and 12, and the p-values of the 3 smallest equal an independent phyper (org.Mm.eg.db)",
        nrow(bpm) >= 1 && identical(as.integer(srow(sM, "rmats", "treatment_up")$genes), 30L) && identical(as.integer(srow(sM, "rmats", "control_up")$genes), 12L) &&
        rel_eq(unname(pm[ids_m]), unname(vapply(ids_m, function(t) indep_p("mouse", "BP", t, S_m$up, um), numeric(1)))))
}
# E: organism not supported, F: GO plus KEGG, G: enrichment switched off, H: KEGG alone
rdZ <- synth_dir("synth_zebrafish", S_h)
okZ <- render(rdZ, "06_enrichment.Rmd", contrasts = dpath(rdZ), organism = "zebrafish")
hZ <- html_of(rdZ, "06_enrichment.html")
check("organism zebrafish: the render succeeds, prints 'organism not supported' and writes no enrichment table",
      okZ && has(hZ, "organism not supported") && length(list.files(file.path(rdZ, "tables"), pattern = "enrich")) == 0)
rdK <- synth_dir("synth_kegg", S_h)
okK <- render(rdK, "06_enrichment.Rmd", contrasts = dpath(rdK), enrich_dbs = "GO,KEGG")
hK <- html_of(rdK, "06_enrichment.html")
check("enrich_dbs 'GO,KEGG': the render succeeds, prints 'require an image rebuild' and still writes the GO tables",
      okK && has(hK, "require an image rebuild") && file.exists(tab_file(rdK, "rmats", "treatment_up", "go_bp")))
rdK2 <- synth_dir("synth_kegg_only", S_h)
okK2 <- render(rdK2, "06_enrichment.Rmd", contrasts = dpath(rdK2), enrich_dbs = "kegg")
check("enrich_dbs 'kegg' alone: the render succeeds, the rebuild message is printed and no GO table is written",
      okK2 && has(html_of(rdK2, "06_enrichment.html"), "require an image rebuild") && length(list.files(file.path(rdK2, "tables"), pattern = "\\.go_")) == 0)
rdO <- synth_dir("synth_off", S_h)
okO <- render(rdO, "06_enrichment.Rmd", contrasts = dpath(rdO), enrich = FALSE)
check("enrich = FALSE: the render succeeds and writes no table", okO && length(list.files(file.path(rdO, "tables"), pattern = "enrich")) == 0)
rdN <- new_dir("no_rds"); writeLines(c("contrast,treatment,control", paste0(C_SYN, ",T1,C1")), dpath(rdN))
okN <- render(rdN, "06_enrichment.Rmd", contrasts = dpath(rdN))
check("a report directory without RDS: 06 renders, names the missing file and writes no table", okN && has(html_of(rdN, "06_enrichment.html"), "not found") && length(list.files(file.path(rdN, "tables"), pattern = "enrich")) == 0)

# --- Rmd 06 and 07 on the real fixtures (FDR 0.9, |dPSI| 0.05, min_reads 0: gene lists exist)
rt_contrasts <- file.path(rt, "261003_a_contrasts.csv"); C <- "GBR_vs_YRI"
mods <- "rmats suppa dexseq_exon edger_exon dtu"
rdR <- new_dir_all("real")
perm <- list(contrasts = rt_contrasts, modules = mods, fdr = 0.9, dpsi = 0.05, min_reads = 0)
ok_all <- TRUE
for (r in c("01_import.Rmd", "02_rmats.Rmd", "03_suppa.Rmd", "04_usage.Rmd", "05_cross_tool.Rmd", "06_enrichment.Rmd")) {
  ok1 <- do.call(render, c(list(rdR, r), perm)); ok_all <- ok_all && ok1
  check(sprintf("%s renders on the real fixtures at FDR 0.9, |dPSI| 0.05, min_reads 0", r), ok1)
}
if (ok_all) {
  sR <- rd_csv(file.path(rdR, "tables", paste0(C, ".enrich.summary.csv")))
  sets <- rd_csv(file.path(rdR, "tables", paste0(C, ".cross_tool.sets.csv")))
  allg <- sR[sR$group == "all", , drop = FALSE]
  check("real fixtures: for each tool the group 'all' has the genes of the cross-tool sets and the universe is its tested genes",
        identical(sort(allg$tool), sort(sets$tool)) && all(as.integer(allg$genes[match(sets$tool, allg$tool)]) == as.integer(sets$significant_genes)) &&
        all(as.integer(allg$universe_genes[match(sets$tool, allg$tool)]) == as.integer(sets$tested_genes)))
  big <- sR[as.integer(sR$genes) >= 10, , drop = FALSE]
  check(sprintf("real fixtures: %d groups have at least 10 genes; each reports a mapped and an unmapped count that add up to its genes", nrow(big)),
        nrow(big) >= 3 && all(!is.na(suppressWarnings(as.integer(big$genes_unmapped)))) && all(as.integer(big$genes_mapped) + as.integer(big$genes_unmapped) == as.integer(big$genes)) &&
        all(as.integer(big$universe_mapped) + as.integer(big$universe_unmapped) == as.integer(big$universe_genes)))
  check("real fixtures: a group with fewer than 10 genes is skipped, a group that ran has at least 10 genes, and at least one group ran",
        all(sR$status[as.integer(sR$genes) < 10] == "skipped") && all(as.integer(sR$genes[sR$status == "ran"]) >= 10) && any(sR$status == "ran"))
  h6 <- html_of(rdR, "06_enrichment.html")
  check("real fixtures: the report shows unmapped counts and the universe", has(h6, "unmapped") && has(h6, "testable"))
  ran <- sR[sR$status == "ran", , drop = FALSE]
  check("real fixtures: every group that ran has its three tables", all(file.exists(unlist(lapply(seq_len(nrow(ran)), function(i) file.path(rdR, "tables", paste0(C, ".enrich.", ran$tool[i], ".", ran$group[i], ".", c("go_bp", "go_mf", "go_cc"), ".csv")))))))
  check("real fixtures: no table is written for a skipped group", {
    sk <- sR[sR$status == "skipped", , drop = FALSE]
    all(!file.exists(unlist(lapply(seq_len(nrow(sk)), function(i) file.path(rdR, "tables", paste0(C, ".enrich.", sk$tool[i], ".", sk$group[i], ".go_bp.csv")))))) })
  # --- Rmd 07 on the complete report directory
  ok7 <- do.call(render, c(list(rdR, "07_summary.Rmd"), perm))
  check("07 renders on the complete report directory", ok7)
  if (ok7) {
    h7 <- paste(readLines(file.path(rdR, "report_summary.html"), warn = FALSE), collapse = "\n")
    dec <- function(x) { x <- gsub("&amp;", "&", x, fixed = TRUE); x <- gsub("&lt;", "<", x, fixed = TRUE); x <- gsub("&gt;", ">", x, fixed = TRUE); gsub("&quot;", "\"", x, fixed = TRUE) }
    rows <- regmatches(h7, gregexpr("(?s)<tr[^>]*>.*?</tr>", h7, perl = TRUE))[[1]]
    cells <- lapply(rows, function(r) trimws(dec(gsub("<[^>]+>", "", regmatches(r, gregexpr("(?s)<t[dh][^>]*>.*?</t[dh]>", r, perl = TRUE))[[1]]))))
    shown <- vapply(cells, paste, character(1), collapse = "\r")
    for (nm in c("rmats.counts", "suppa.counts", "usage.counts", "cross_tool.sets", "cross_tool.concordance", "enrich.summary")) {
      f <- file.path(rdR, "tables", paste0(C, ".", nm, ".csv")); d <- rd_csv(f)
      want <- vapply(seq_len(nrow(d)), function(i) paste(unlist(d[i, ]), collapse = "\r"), character(1))
      check(sprintf("07 shows the %d rows of %s.%s.csv exactly as the table has them", nrow(d), C, nm), nrow(d) > 0 && all(want %in% shown))
    }
    hrefs <- sub("^href=\"(.*)\"$", "\\1", regmatches(h7, gregexpr("href=\"[^\"]*\"", h7))[[1]])
    hrefs <- unique(dec(hrefs)); hrefs <- hrefs[grepl("^(0[1-9]_[^/]*\\.html|figures/|tables/)", hrefs)]
    local_files <- c(list.files(rdR, pattern = "^0[1-6]_.*\\.html$"), file.path("figures", list.files(file.path(rdR, "figures"))), file.path("tables", list.files(file.path(rdR, "tables"))))
    check("07 links every report, figure and table that exists in the report directory", all(local_files %in% hrefs))
    check("07 links only files that exist (every relative link resolves)", all(file.exists(file.path(rdR, hrefs))))
    check("07 lists the six reports as produced (no 'not produced' line)", !grepl("not produced", h7, fixed = TRUE))
    check("07 states the thresholds of the render (0.9, 0.05, 0, JC)", all(vapply(c("FDR", "0.9", "dPSI", "0.05", "JC", "minimum read support"), function(s) grepl(s, h7, fixed = TRUE), logical(1))))
    check("07 names the contrast with treatment GBR and control YRI and lists the modules", grepl("GBR_vs_YRI", h7, fixed = TRUE) && grepl("GBR", h7, fixed = TRUE) && grepl("YRI", h7, fixed = TRUE) && grepl("dexseq_exon", h7, fixed = TRUE))
  }
}

# --- Rmd 07 on a scratch report directory: two reports, one figure, one table; a third report absent; no counts table
rdS <- new_dir("scratch07")
writeLines("<html><body><p>import</p></body></html>", file.path(rdS, "01_import.html"))
writeLines(c("<html><body>", "<p>suppa</p>", "<p>WARNING: no rMATS sample list found for X &amp; Y</p>", "</body></html>"), file.path(rdS, "03_suppa.html"))
dir.create(file.path(rdS, "figures")); dir.create(file.path(rdS, "tables"))
file.create(file.path(rdS, "figures", "SYN.example.png")); writeLines("a,b\n1,2", file.path(rdS, "tables", "SYN.example.csv"))
writeLines(c("contrast,treatment,control", "SYN,T1,C1", "OTHER,T2,C2"), dpath(rdS))
ok7s <- render(rdS, "07_summary.Rmd", contrasts = dpath(rdS), modules = "rmats suppa", fdr = 0.037, dpsi = 0.123, min_reads = 17, rmats_type = "JCEC")
h7s <- paste(readLines(file.path(rdS, "report_summary.html"), warn = FALSE), collapse = "\n")
h7n <- gsub("[[:space:]]+", " ", h7s)
check("07 renders on a directory with two reports, one figure and one table, and without any counts table", ok7s)
check("07 links the two reports that exist", has(h7s, "href=\"01_import.html\"") && has(h7s, "href=\"03_suppa.html\""))
check("07 does not link the reports that are absent (02, 04, 05, 06, and itself is not listed) and says 'not produced'",
      !any(vapply(c("02_rmats.html", "04_usage.html", "05_cross_tool.html", "06_enrichment.html"), function(f) has(h7s, paste0("href=\"", f, "\"")), logical(1))) && has(h7n, "not produced"))
check("07 links the figure and the table that exist and no other", has(h7s, "href=\"figures/SYN.example.png\"") && has(h7s, "href=\"tables/SYN.example.csv\"") &&
      length(regmatches(h7s, gregexpr("href=\"(figures|tables)/[^\"]*\"", h7s))[[1]]) == 2)
check("07 reports the missing counts tables instead of failing (not found)", has(h7n, "not found"))
check("07 states the thresholds read from the parameters (0.037, 0.123, 17, JCEC)", all(vapply(c("0.037", "0.123", "17", "JCEC"), function(s) has(h7n, s), logical(1))))
check("07 lists the contrasts with treatment and control (SYN: T1, C1; OTHER: T2, C2) and the modules found", all(vapply(c("SYN", "T1", "C1", "OTHER", "T2", "C2", "rmats, suppa"), function(s) has(h7n, s), logical(1))))
check("07 prints the sign rules, with the same wording as the other reports", has(h7n, "treatment minus control") && has(h7n, "sign flipped") && has(h7n, "IncLevelDifference") && has(h7n, "kept"))
check("07 collects the WARNING line of a report (with the entity decoded) and names the report", has(h7n, "WARNING: no rMATS sample list found for X &amp; Y") || has(h7n, "WARNING: no rMATS sample list found for X & Y"))
check("07 says that the warnings are recovered from the HTML of the reports", has(h7n, "recovered from the HTML"))
# a directory with no WARNING line and one counts table that exists
rdW <- new_dir("scratch07b"); writeLines("<html><body><p>import</p></body></html>", file.path(rdW, "01_import.html"))
dir.create(file.path(rdW, "tables")); writeLines(c("type,events,significant", "SE,10,3"), file.path(rdW, "tables", "SYN.rmats.counts.csv"))
writeLines(c("contrast,treatment,control", "SYN,T1,C1"), dpath(rdW))
ok7w <- render(rdW, "07_summary.Rmd", contrasts = dpath(rdW), modules = "rmats")
h7w <- gsub("[[:space:]]+", " ", paste(readLines(file.path(rdW, "report_summary.html"), warn = FALSE), collapse = " "))
check("07 renders with one counts table, shows it (SE 10 3) and says that no WARNING line was found", ok7w && has(h7w, ">SE</td>") && has(h7w, ">10</td>") && has(h7w, "No WARNING line"))
# unit: warnings from html
wf <- file.path(work, "w.html"); writeLines(c("<p>a <code>b</code> WARNING: x <code>y</code> z</p>", "<p>warning: lower case is not a hit</p>", "<p>WARNING: x <code>y</code> z</p>"), wf)
check("html_warnings strips the tags, keeps each WARNING line once and ignores lower case", identical(E7$html_warnings(wf), c("WARNING: x y z", "WARNING: x y z")[1]) || identical(E7$html_warnings(wf), "WARNING: x y z"))
check("list_files lists only what exists, sorted", { d <- file.path(work, "lf"); dir.create(d); file.create(file.path(d, c("b.png", "a.png", "c.txt"))); identical(E7$list_files(d, "\\.png$"), c("a.png", "b.png")) && length(E7$list_files(file.path(work, "absent"), "\\.png$")) == 0 })
check("read_counts reads every field as text and returns NULL for a missing file", { f <- file.path(work, "c.csv"); writeLines(c("a,b", "1,NA", "2,0.50"), f); d <- E7$read_counts(f); identical(d$b, c("NA", "0.50")) && is.null(E7$read_counts(file.path(work, "absent.csv"))) })
finish()
