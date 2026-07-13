FROM rocker/r-ver:4.4

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
    && rm -rf /var/lib/apt/lists/*

# CRAN + Bioconductor packages
# Install in dependency order to minimise recompilation on layer cache miss
RUN Rscript -e "\
  install.packages(c( \
    'BiocManager', 'remotes', 'pak', \
    'plumber', 'future', 'promises', 'uuid', 'jsonlite', \
    'tidyverse', 'patchwork', 'igraph', 'ggraph', \
    'WGCNA' \
  ), repos='https://cloud.r-project.org')"

# Seurat and SeuratDisk
# NOTE: SeuratDisk (LoadH5Seurat) is unmaintained as of 2022.
# If your .h5Seurat files were written by Seurat v5, switch to
# readRDS() + .rds format and remove this dependency.
RUN Rscript -e "\
  install.packages('Seurat', repos='https://cloud.r-project.org'); \
  remotes::install_github('mojaveazure/seurat-disk')"

# Bioconductor packages required by hdWGCNA
RUN Rscript -e "\
  BiocManager::install(c('BiocGenerics', 'GenomicRanges', 'GeneOverlap', 'UCell'), \
    ask=FALSE, update=FALSE)"

# hdWGCNA from GitHub dev branch
RUN Rscript -e "remotes::install_github('smorabit/hdWGCNA', ref='dev')"

WORKDIR /app
COPY entrypoint.R plumber.R ./

EXPOSE 8100

CMD ["Rscript", "entrypoint.R"]
