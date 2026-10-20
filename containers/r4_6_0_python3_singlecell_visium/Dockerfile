# R 4.6.0 + Seurat 5 (stable) + single-cell stack
# Fixed: removed visium-hd branch override, pinned Seurat/SeuratObject,
#        SeuratWrappers from R-universe, scDblFinder consolidated.
FROM rocker/verse:4.6.0
LABEL maintainer="yann" \
      description="Rocker Verse 4.6.0 + Seurat 5 (stable), Signac, Azimuth, SeuratWrappers, harmony, CellBender, scVI"

ENV DEBIAN_FRONTEND=noninteractive \
    RENV_CONFIG_REPOS_OVERRIDE=https://cloud.r-project.org \
    R_REMOTES_NO_ERRORS_FROM_WARNINGS=true

# Build-time GitHub token to avoid the unauthenticated 60 req/hr API rate
# limit during the many remotes::install_github() calls below. Not an ENV,
# so it is not persisted into the final image's runtime environment or
# layer metadata. Pass at build time: docker build --build-arg GITHUB_PAT=<token> ...
ARG GITHUB_PAT

# ---------------------------------------------------------------------------
# System dependencies
# ---------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    pkg-config \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libcairo2-dev \
    libxt-dev \
    zlib1g-dev \
    libomp-dev \
    libgsl-dev \
    libhdf5-dev \
    libfftw3-dev \
    libsodium-dev \
    libv8-dev \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libglpk40 libglpk-dev \
    libarmadillo-dev \
    libboost-all-dev \
    python3 python3-venv python3-dev python3-pip \
 && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# Step 1: Pin SeuratObject FIRST, then Seurat — must be the same v5 release.
# Installing SeuratObject before Seurat ensures no downgrade happens later.
# DO NOT install the visium-hd GitHub branch — those features are in Seurat 5.x CRAN.
# ---------------------------------------------------------------------------
RUN R --vanilla -q -e " \
  options(Ncpus = max(1L, parallel::detectCores() - 1L), \
          repos = c(CRAN = 'https://cloud.r-project.org')); \
  install.packages('remotes'); \
  install.packages('BiocManager'); \
  remotes::install_version('SeuratObject', version = '5.0.2', \
    repos = 'https://cloud.r-project.org', upgrade = 'never'); \
  remotes::install_version('Seurat', version = '5.2.1', \
    repos = 'https://cloud.r-project.org', upgrade = 'never'); \
"

# ---------------------------------------------------------------------------
# Step 2: CRAN packages
# ---------------------------------------------------------------------------
RUN R --vanilla -q -e " \
  options(Ncpus = max(1L, parallel::detectCores() - 1L), \
          repos = c(CRAN = 'https://cloud.r-project.org')); \
  cran_pkgs <- c( \
    'Signac', \
    'Matrix', 'dplyr', 'tidyr', 'purrr', 'stringr', \
    'cowplot', 'ggplot2', 'ggrepel', 'patchwork', 'ggrastr', 'ggpubr', 'ggseqlogo', \
    'data.table', 'matrixStats', \
    'clustree', 'SoupX', 'harmony', \
    'openxlsx', 'writexl', 'kableExtra', \
    'gprofiler2', 'VAM', 'hexbin', 'ComplexUpset', \
    'HGNChelper', 'scales', \
    'rliger', 'terra' \
  ); \
  to_install <- setdiff(cran_pkgs, rownames(installed.packages())); \
  if (length(to_install)) install.packages(to_install); \
"

# ---------------------------------------------------------------------------
# Step 3: Bioconductor packages
# ---------------------------------------------------------------------------
RUN R --vanilla -q -e " \
  options(Ncpus = max(1L, parallel::detectCores() - 1L), \
          repos = c(CRAN = 'https://cloud.r-project.org')); \
  BiocManager::install(c( \
    'SingleR', 'celldex', \
    'scDblFinder', 'SingleCellExperiment', \
    'slingshot', 'tradeSeq', 'condiments', 'scRNAseq', \
    'edgeR', 'apeglm', 'DESeq2', \
    'clusterProfiler', \
    'dittoSeq', 'scDataviz', \
    'GSEABase', 'ComplexHeatmap', 'GSVA', 'UCell', \
    'TFBSTools', 'motifmatchr', 'chromVAR', \
    'BiocGenerics', 'DelayedArray', 'SummarizedExperiment', \
    'batchelor', 'HDF5Array', \
    'pcaMethods', \
    'BSgenome.Mmusculus.UCSC.mm10', 'EnsDb.Mmusculus.v79', 'JASPAR2024' \
  ), ask = FALSE, update = FALSE); \
"

# ---------------------------------------------------------------------------
# Step 4: GitHub packages — installed AFTER Seurat is finalised.
# No try()/silent=TRUE so failures are visible at build time.
# ---------------------------------------------------------------------------
RUN R --vanilla -q -e " \
  options(Ncpus = max(1L, parallel::detectCores() - 1L), \
          repos = c(CRAN = 'https://cloud.r-project.org')); \
  remotes::install_github('satijalab/azimuth', ref = 'master',  upgrade = 'never'); \
  remotes::install_github('chris-mcginnis-ucsf/DoubletFinder',            upgrade = 'never'); \
  remotes::install_github('mojaveazure/seurat-disk',                       upgrade = 'never'); \
  remotes::install_github('cole-trapnell-lab/monocle3',                    upgrade = 'never'); \
  remotes::install_github('samuel-marsh/scCustomize',                      upgrade = 'never'); \
  remotes::install_github('saezlab/liana',                                 upgrade = 'never'); \
  remotes::install_github('velocyto-team/velocyto.R',                      upgrade = 'never'); \
  remotes::install_github('hypercompetent/colorway',                       upgrade = 'never'); \
  remotes::install_github('kharchenkolab/pagoda2',                         upgrade = 'never'); \
  remotes::install_github('powellgenomicslab/DropletQC',                   upgrade = 'never'); \
"

# ---------------------------------------------------------------------------
# Python — CellBender (CPU build, separate venv)
# ---------------------------------------------------------------------------
RUN python3 -m venv /opt/py/cellbender \
 && /opt/py/cellbender/bin/pip install --upgrade --no-cache-dir pip \
 && /opt/py/cellbender/bin/pip install --no-cache-dir \
      --index-url https://download.pytorch.org/whl/cpu torch \
 && /opt/py/cellbender/bin/pip install --no-cache-dir cellbender \
 && ln -s /opt/py/cellbender/bin/cellbender /usr/local/bin/cellbender

# ---------------------------------------------------------------------------
# Python — scVI (separate venv)
# ---------------------------------------------------------------------------
RUN python3 -m venv /opt/py/scvi \
 && /opt/py/scvi/bin/pip install --upgrade --no-cache-dir pip setuptools wheel \
 && /opt/py/scvi/bin/pip install --no-cache-dir \
      --index-url https://download.pytorch.org/whl/cpu \
      torch torchvision torchaudio \
 && /opt/py/scvi/bin/pip install --no-cache-dir \
      scvi-tools scanpy numpy pandas \
 && ln -s /opt/py/scvi/bin/python3 /usr/local/bin/scvi-python

# ---------------------------------------------------------------------------
# SeuratWrappers — installed LAST so a failure doesn't break cached layers above.
# Tries R-universe first (pre-built binary), falls back to GitHub HEAD.
# ---------------------------------------------------------------------------
RUN R --vanilla -q -e " \
  options(repos = c(CRAN = 'https://cloud.r-project.org')); \
  remotes::install_github('satijalab/seurat-wrappers', ref = 'HEAD', upgrade = 'never'); \
  stopifnot(requireNamespace('SeuratWrappers', quietly = TRUE)); \
  message('SeuratWrappers installed: ', packageVersion('SeuratWrappers')); \
"

WORKDIR /work
CMD ["bash"]

RUN R -e "install.packages('msigdbr')"

# ---------------------------------------------------------------------------
# Integration-benchmark metrics: LISI (iLISI / cLISI) and kBET.
# Neither is on CRAN/Bioconductor — installed from GitHub. Used to score
# batch mixing vs biological conservation when comparing integration methods
# (Harmony / CCA / RPCA / FastMNN / scVI). No try()/silent so build failures
# are visible. Package names are 'lisi' and 'kBET'.
# ---------------------------------------------------------------------------
RUN R --vanilla -q -e " \
  options(Ncpus = max(1L, parallel::detectCores() - 1L), \
          repos = c(CRAN = 'https://cloud.r-project.org')); \
  remotes::install_github('immunogenomics/LISI', ref = 'HEAD', upgrade = 'never'); \
  remotes::install_github('theislab/kBET',       ref = 'HEAD', upgrade = 'never'); \
  stopifnot(requireNamespace('lisi', quietly = TRUE), \
            requireNamespace('kBET', quietly = TRUE)); \
  message('lisi installed: ', packageVersion('lisi'), \
          ' | kBET installed: ', packageVersion('kBET')); \
"

# ---------------------------------------------------------------------------
# scVI is already installed above in the /opt/py/scvi venv (scvi-tools 1.4.3).
# SeuratWrappers::scVIIntegration drives scVI through reticulate, which finds
# no Python by default. Point reticulate at the venv interpreter so the R-side
# scVI integration path works out of the box. Override at runtime with
# RETICULATE_PYTHON / reticulate::use_python() if a different env is needed.
# ---------------------------------------------------------------------------
ENV RETICULATE_PYTHON=/opt/py/scvi/bin/python3

