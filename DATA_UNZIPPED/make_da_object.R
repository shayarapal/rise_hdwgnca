# Build a DA-only analysis object for the hdWGCNA bridge.
#
# The full object is 107,113 nuclei, but the CACNA1D network is built on the 8,216
# dopaminergic ones (7.7%). The r-service runs plan(multisession, workers = 2), so each
# worker holds its own copy -- loading the full object twice is what hits R's 32 GB
# vector limit. Two reductions, neither of which changes any hdWGCNA result:
#
#   1. subset to cell_type == "Dopaminergic Neurons"
#   2. drop the scale.data layer (2,000 x 107,113 doubles ~ 1.7 GB). ScaleData() output is
#      used for PCA/clustering, which are already done; hdWGCNA reads counts/data + the pca
#      reduction, never scale.data.
#
# Kept deliberately: the `pca` reduction (MetacellsByGroups hardcodes reduction = "pca" in
# r-service/plumber.R) and VariableFeatures (SelectNetworkGenes defaults to
# gene_select = "variable", i.e. the 2,000 HVGs, which were computed on the full dataset).
suppressPackageStartupMessages(library(Seurat))

# Data root. Resolution order:
#   1. SHARED_DIR env var (host runs)      -> $SHARED_DIR/DATA_UNZIPPED
#   2. /shared/DATA_UNZIPPED if it exists  -> inside the r-service container
#   3. ./DATA_UNZIPPED relative to cwd     -> repo checkout fallback
# The data itself lives outside the repo; see .env.example and the README.
root <- local({
  sd <- Sys.getenv("SHARED_DIR", "")
  if (nzchar(sd)) return(file.path(sd, "DATA_UNZIPPED"))
  if (dir.exists("/shared/DATA_UNZIPPED")) return("/shared/DATA_UNZIPPED")
  "DATA_UNZIPPED"
})
if (!dir.exists(root)) {
  stop("Data root not found: ", root,
       "\nSet SHARED_DIR to the directory holding DATA_UNZIPPED/ ",
       "(see .env.example), or download the raw data from GEO accession GSE243639.")
}
seu  <- readRDS(file.path(root, "seurat_GSE243639_SNc.rds"))

cat("=== full object\n")
cat(sprintf("  cells        %d\n  size in RAM  %.2f GB\n", ncol(seu),
            as.numeric(object.size(seu)) / 2^30))

seu <- subset(seu, subset = cell_type == "Dopaminergic Neurons")
stopifnot(ncol(seu) > 0)

# hdWGCNA never reads scale.data; PCA/clustering that needed it are already done.
if ("scale.data" %in% Layers(seu[["RNA"]])) {
  seu[["RNA"]]$scale.data <- NULL
}

cat("\n=== DA-only object\n")
cat(sprintf("  cells        %d\n  size in RAM  %.2f GB\n", ncol(seu),
            as.numeric(object.size(seu)) / 2^30))
cat("  reductions   ", paste(Reductions(seu), collapse = ", "), "\n")
cat("  layers       ", paste(Layers(seu[["RNA"]]), collapse = ", "), "\n")
cat("  var features ", length(VariableFeatures(seu)), "\n")
cat("  cell_type    ", paste(levels(factor(seu$cell_type)), collapse = ", "), "\n")

# MetacellsByGroups drops any cell_type x sample_id group below min_cells (default 100).
tb <- table(seu$sample_id)
cond <- unique(seu[[c("sample_id", "condition")]])
keep <- names(tb)[tb >= 100]
cat(sprintf("  donors >=100 DA cells: %d of %d (%d PD, %d control) -- these drive the network\n",
    length(keep), length(tb),
    sum(cond$condition[match(keep, cond$sample_id)] == "PD"),
    sum(cond$condition[match(keep, cond$sample_id)] == "control")))

out <- file.path(root, "seurat_GSE243639_SNc_DA.rds")
saveRDS(seu, out)
cat("\nSaved", out, sprintf("(%.2f GB on disk)\n", file.size(out) / 2^30))
