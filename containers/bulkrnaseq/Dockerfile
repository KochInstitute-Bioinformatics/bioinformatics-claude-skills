FROM rocker/verse:4.5.1
LABEL description="Docker Rstudio image KI IGB bulk rnaseq and crispr - February 2026"

RUN apt-get update && apt-get install -y libhdf5-dev build-essential libxt6 libxml2-dev libssl-dev \
libv8-dev libsodium-dev libglpk40 libgdal-dev libboost-dev libomp-dev libbamtools-dev \
libboost-iostreams-dev libboost-log-dev libbz2-dev libboost-system-dev libboost-test-dev \
libcurl4-openssl-dev libz-dev libarmadillo-dev

RUN R -e "install.packages(c('cowplot','knitr','msigdbr','remotes','metap','devtools','R.utils','ggalt','BiocManager','openxlsx','XML','gprofiler2','cluster'),repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('gplots','MetaCycle','ashr'),repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('markdown'),repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('ggpubr'),repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('effectsize'),repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('DESeq2','tximport','Rsamtools','Rhtslib','GenomicAlignments','ComplexHeatmap','rtracklayer','fgsea','apeglm','edgeR'))"
RUN R -e "BiocManager::install(c('MAGeCKFlute'))"
RUN R -e "BiocManager::install(c('org.Hs.eg.db','org.Mm.eg.db'))"
RUN R -e "install.packages(c('rebus'),repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('GGally'),repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('VennDiagram'),repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('ggVennDiagram'),repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('UpSetR'),repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('colorRamp2'),repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('esquisse'),repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('biomaRt'))"
RUN R -e "install.packages(c('gganimate'),repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('gifski'),repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('ggbreak'),repos='http://cran.rstudio.com/')"
RUN R -e "install.packages(c('tidyverse', 'reprex', 'matrixStats', 'ggrepel', 'stringr', 'ComplexUpset', 'scales', 'viridis', 'outliers', 'showtext'),repos='http://cran.rstudio.com/')"

# Added for the bulk-rnaseq-pipeline skill's new modules: interactive
# consolidated report (DT), ORA (clusterProfiler), mixed-effects DE
# (variancePartition/dream). crosstalk/plotly/shiny/rsconnect are already
# present via rocker/verse and existing package dependencies.
RUN R -e "install.packages(c('DT'),repos='http://cran.rstudio.com/')"
RUN R -e " \
  options(repos = c(CRAN = 'https://cloud.r-project.org')); \
  install.packages('lme4'); \
  BiocManager::install(c('ggtangle', 'clusterProfiler', 'DOSE', 'GO.db', 'variancePartition'), ask = FALSE, update = FALSE); \
  stopifnot(requireNamespace('clusterProfiler', quietly = TRUE), \
            requireNamespace('DOSE', quietly = TRUE), \
            requireNamespace('GO.db', quietly = TRUE), \
            requireNamespace('variancePartition', quietly = TRUE)) \
"

# GseaVis is GitHub-only (not on CRAN or Bioconductor) - a prior attempt to
# install it via install.packages(repos=CRAN) silently no-opped (warning, not
# a build failure), leaving it missing from the image. It also depends on
# clusterProfiler/DOSE/GO.db (installed just above), so it must come after
# that block, not before. Install via remotes (already installed above) and
# fail the build loudly if it doesn't resolve.
#
# 'enrichit' (CRAN, published 2026-07-01) is a newer C++ backend that
# clusterProfiler's plotting internals (used by GseaVis::gseaNb) lazy-load at
# CALL time via loadNamespace(), not at package-load time - so a
# requireNamespace('GseaVis') check alone passes even when enrichit is
# missing, and the failure only surfaces when gseaNb() actually runs. Install
# it explicitly and smoke-test gseaNb() with a tiny synthetic GSEA result so
# the build fails loudly instead of silently at render time.
RUN R -e " \
  options(repos = c(CRAN = 'https://cloud.r-project.org')); \
  install.packages('enrichit'); \
  remotes::install_github('junjunlab/GseaVis', dependencies = TRUE, upgrade = 'never'); \
  stopifnot(requireNamespace('GseaVis', quietly = TRUE), \
            requireNamespace('enrichit', quietly = TRUE)); \
  library(clusterProfiler); \
  set.seed(1); \
  gl <- sort(rnorm(200), decreasing = TRUE); \
  names(gl) <- paste0('GENE', seq_along(gl)); \
  t2g <- data.frame(term = rep('SET_A', 20), gene = names(gl)[1:20]); \
  gs <- GSEA(geneList = gl, TERM2GENE = t2g, pvalueCutoff = 1, minGSSize = 5, seed = TRUE); \
  stopifnot(nrow(gs@result) > 0); \
  p <- GseaVis::gseaNb(object = gs, geneSetID = gs@result\$ID[1], subPlot = 2, addPval = TRUE); \
  stopifnot(!is.null(p)) \
"
# Added for the bulk-rnaseq-pipeline skill's DTU module (Rmd 04):
# DRIMSeq (filtering) -> DEXSeq (testing) -> stageR (gene/transcript OFDR).
RUN R -e " \
  options(repos = c(CRAN = 'https://cloud.r-project.org')); \
  BiocManager::install(c('DRIMSeq', 'DEXSeq', 'stageR'), ask = FALSE, update = FALSE); \
  stopifnot(requireNamespace('DRIMSeq', quietly = TRUE), \
            requireNamespace('DEXSeq', quietly = TRUE), \
            requireNamespace('stageR', quietly = TRUE)) \
"
