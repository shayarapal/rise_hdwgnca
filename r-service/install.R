# Local R package installation for r-service
# Run once before using `make run-r`:
#   Rscript r-service/install.R
#
# WARNING: First run takes 1-3 hours depending on your machine and CRAN mirror.
# Docker users do NOT need this — the Dockerfile handles it.

options(repos = c(CRAN = "https://cloud.r-project.org"))
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")

cat("==> Installing CRAN packages...\n")
install.packages(c(
  "BiocManager",
  "remotes",
  "plumber",
  "future",
  "promises",
  "uuid",
  "jsonlite",
  "Seurat",
  "WGCNA",
  "igraph",
  "tidyverse",
  "ggraph",
  "patchwork"
), dependencies = TRUE)

cat("==> Installing SeuratDisk (LoadH5Seurat)...\n")
# NOTE: SeuratDisk is unmaintained since 2022. If your .h5Seurat files were
# written by Seurat v5, switch to readRDS() + .rds format instead.
remotes::install_github("mojaveazure/seurat-disk")

cat("==> Installing Bioconductor packages...\n")
BiocManager::install(
  c("BiocGenerics", "GenomicRanges", "GeneOverlap", "UCell"),
  ask    = FALSE,
  update = FALSE
)

cat("==> Installing hdWGCNA (dev branch)...\n")
remotes::install_github("smorabit/hdWGCNA", ref = "dev")

cat("==> All packages installed successfully.\n")
