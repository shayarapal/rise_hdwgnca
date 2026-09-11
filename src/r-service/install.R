# Local R package installation for r-service
# Run once before using `make run-r`:
#   Rscript r-service/install.R
#
# WARNING: First run takes 1-3 hours depending on your machine and CRAN mirror.
# Docker users do NOT need this — the Dockerfile handles it.

options(repos = c(CRAN = "https://cloud.r-project.org"))
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = "true")

cat("==> Installing BiocManager + remotes...\n")
install.packages(c("BiocManager", "remotes"))

cat("==> Installing Bioconductor packages...\n")
# GO.db/impute/preprocessCore/AnnotationDbi are WGCNA's dependencies and must be installed
# BEFORE WGCNA (they are not on CRAN). The rest are hdWGCNA's deps. Getting this order wrong
# is what silently ships an image with no WGCNA/hdWGCNA — see the Dockerfile comment.
BiocManager::install(
  c("GO.db", "impute", "preprocessCore", "AnnotationDbi",
    "BiocGenerics", "GenomicRanges", "GeneOverlap", "UCell"),
  ask    = FALSE,
  update = FALSE
)

cat("==> Installing CRAN packages (WGCNA now that its Bioc deps exist)...\n")
install.packages(c(
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
  "patchwork",
  "cowplot",
  "enrichR"
), dependencies = TRUE)

cat("==> Installing SeuratDisk (LoadH5Seurat)...\n")
# NOTE: SeuratDisk is unmaintained since 2022. If your .h5Seurat files were
# written by Seurat v5, switch to readRDS() + .rds format instead.
remotes::install_github("mojaveazure/seurat-disk", ref = "877d4e18ab38c686f5db54f8cd290274ccdbe295")

cat("==> Installing hdWGCNA (pinned commit)...\n")
# v0.4.12, 2026-07-29. Do not float back to ref = "dev" -- module assignments are
# version-sensitive and the paper's numbers came from this commit.
remotes::install_github("smorabit/hdWGCNA", ref = "e3344d1f7bbac4264adf94f7aa31e0802fa8282d")

# Fail loudly if the core stack did not actually install (remotes only warns on missing deps).
cat("==> Verifying core stack loads...\n")
suppressMessages({library(hdWGCNA); library(WGCNA)})
cat("==> All packages installed successfully. hdWGCNA",
    as.character(packageVersion("hdWGCNA")), "OK\n")
