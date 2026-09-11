FROM rocker/r-ver:4.4.3

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    libz-dev \
    libgit2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libsodium-dev \
    libuv1-dev \
    && rm -rf /var/lib/apt/lists/*
# libsodium-dev is required to build the R `sodium` package, which `plumber` depends on.
# libuv1-dev is required to build `fs`, which cascades to sass -> bslib -> rmarkdown -> tidyverse.
# Without either, install.packages() silently skips the failed package and the build exits 0,
# then the container crash-loops at `library(plumber)` / `library(tidyverse)`.

# CRAN packages. WGCNA is deliberately NOT here: it depends on Bioconductor packages
# (GO.db, impute, preprocessCore, AnnotationDbi) that CRAN cannot supply, so it must be
# installed AFTER the Bioconductor step below.
RUN Rscript -e "\
  install.packages(c( \
    'BiocManager', 'remotes', 'pak', \
    'plumber', 'future', 'promises', 'uuid', 'jsonlite', \
    'tidyverse', 'patchwork', 'igraph', 'ggraph', 'cowplot', 'enrichR', 'corrplot' \
  ), repos='https://cloud.r-project.org')"

# Bioconductor packages. The first four are WGCNA's dependencies and MUST precede it —
# without them WGCNA fails to install, hdWGCNA then fails, yet remotes only WARNS, so the
# build would exit 0 and ship an image with no hdWGCNA. The last four are hdWGCNA's own deps.
RUN Rscript -e "\
  BiocManager::install(c('GO.db', 'impute', 'preprocessCore', 'AnnotationDbi', \
    'BiocGenerics', 'GenomicRanges', 'GeneOverlap', 'UCell'), \
    ask=FALSE, update=FALSE)"

# WGCNA now that its Bioconductor dependencies are present
RUN Rscript -e "install.packages('WGCNA', repos='https://cloud.r-project.org')"

# Seurat and SeuratDisk
# NOTE: SeuratDisk (LoadH5Seurat) is unmaintained as of 2022.
# If your .h5Seurat files were written by Seurat v5, switch to
# readRDS() + .rds format and remove this dependency.
RUN Rscript -e "\
  install.packages('Seurat', repos='https://cloud.r-project.org'); \
  remotes::install_github('mojaveazure/seurat-disk', ref='877d4e18ab38c686f5db54f8cd290274ccdbe295')"

# hdWGCNA pinned to the dev-branch commit that produced the published analysis
# (v0.4.12, 2026-07-29). Do not float this back to ref='dev': module assignments
# are version-sensitive and the paper's numbers came from this commit.
RUN Rscript -e "remotes::install_github('smorabit/hdWGCNA', ref='e3344d1f7bbac4264adf94f7aa31e0802fa8282d')"

# harmony -- used by the DATA_*/build_seurat_*.R scripts for animal/donor batch
# correction before saving a Seurat object (not called by the bridge itself).
# R.utils -- data.table::fread() needs it to read .gz files directly.
RUN Rscript -e "install.packages(c('harmony', 'R.utils'), repos='https://cloud.r-project.org')"

# ggforestplot -- hdWGCNA's PlotDMEsLollipop loads this at call time and, if
# missing, tries to auto-install it itself (which fails without devtools and
# crashes the caller). Install it ahead of time instead.
RUN Rscript -e "remotes::install_github('NightingaleHealth/ggforestplot', ref='547617e63fa481a5f28ffc56c07d46be4af688b2')"

# Build-time tripwire: install.packages() and remotes::install_github() only WARN on a failed
# package and let the build exit 0 — so any missing runtime dependency ships silently and the
# container crash-loops at startup. This loads EVERY package the service uses at runtime
# (entrypoint.R + plumber.R top-level + pipeline_libs()); a missing one fails the build loudly.
RUN Rscript -e "pkgs <- c('plumber','future','promises','uuid','jsonlite','Seurat','SeuratDisk','WGCNA','hdWGCNA','tidyverse','cowplot','patchwork','igraph','enrichR'); invisible(lapply(pkgs, function(p) suppressMessages(library(p, character.only=TRUE)))); cat('ALL', length(pkgs), 'runtime packages load OK; hdWGCNA', as.character(packageVersion('hdWGCNA')), '\n')"

WORKDIR /app
# seuratdisk_compat.R is sourced by plumber.R at startup — without it the container cannot boot.
COPY entrypoint.R plumber.R seuratdisk_compat.R ./

EXPOSE 8100

CMD ["Rscript", "entrypoint.R"]
