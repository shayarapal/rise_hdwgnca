log_mem <- function(label) {
  gc()
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), label))
  flush(stdout())
}

library(Seurat)
library(dplyr)
library(data.table)
library(Matrix)
library(harmony)

# Empirical re-derivation of the mouse nFeature_RNA QC floor, mirroring the
# diagnostic that justified the human data's 750-gene floor (a documented
# ~7,874-cell "shadow" cluster at floor=200, only ~55% TH+). The mouse build
# script's 500 floor cites "the same rationale" but was never independently
# checked against this dataset's own shadow-cluster evidence -- this does
# that check. Only the LOWER floor is being tested; upper bound (8000) and
# percent.mt (<5) are unchanged and not in question here.

root <- "/shared/DATA_GSE233866"
raw  <- file.path(root, "GSE233866_untreated_counts.csv.gz")

header_line <- readLines(gzfile(raw), n = 1)
cell_barcodes <- strsplit(header_line, ",", fixed = TRUE)[[1]]
dt <- fread(raw, header = FALSE, skip = 1)
stopifnot(ncol(dt) == length(cell_barcodes) + 1)
gene_names <- dt[[1]]
counts <- as.matrix(dt[, -1, with = FALSE])
rownames(counts) <- gene_names
colnames(counts) <- cell_barcodes
rm(dt); gc()
counts <- Matrix(counts, sparse = TRUE)
sample_id <- sub("_[ACGT]+-1$", "", colnames(counts))

seu <- CreateSeuratObject(counts = counts, project = "GSE233866_qcdiag",
                          min.cells = 3, min.features = 200)
seu$sample_id <- sample_id
rm(counts); gc()
log_mem("loaded, floor=200 only")

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
cat("Pre-QC (floor=200) cells:", ncol(seu), "\n")
cat("nFeature_RNA distribution at floor=200:\n")
print(summary(seu$nFeature_RNA))
print(quantile(seu$nFeature_RNA, probs = c(0.01, 0.05, 0.10, 0.25, 0.50)))

# Only apply the UPPER bound + mito filter here -- deliberately leave the
# lower floor at 200 so any shadow cluster stays in the object to be found.
seu <- subset(seu, subset = nFeature_RNA < 8000 & percent.mt < 5)
cat("Post upper-bound/mito filter cells:", ncol(seu), "\n")

seu <- NormalizeData(seu) |>
       FindVariableFeatures(nfeatures = 2000) |>
       ScaleData() |>
       RunPCA(npcs = 30)
seu <- harmony::RunHarmony(seu, group.by.vars = "sample_id")
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30) |>
       FindNeighbors(reduction = "harmony", dims = 1:30) |>
       FindClusters(resolution = 0.5)
log_mem("clustered at floor=200")

da_markers <- intersect(c("Th", "Slc6a3", "Ddc", "Slc18a2"), rownames(seu))
seu <- AddModuleScore(seu, features = list(da_markers), name = "da_score_")
da_col <- grep("^da_score_", colnames(seu[[]]), value = TRUE)[1]

th_expr <- FetchData(seu, vars = "Th")[, 1]
seu$th_positive <- th_expr > 0

diag_df <- data.frame(
  cluster       = Idents(seu),
  nFeature_RNA  = seu$nFeature_RNA,
  percent_mt    = seu$percent.mt,
  da_score      = seu[[da_col]][, 1],
  th_positive   = seu$th_positive
)

summary_by_cluster <- diag_df %>%
  group_by(cluster) %>%
  summarise(
    n_cells           = n(),
    median_nFeature   = median(nFeature_RNA),
    mean_da_score     = mean(da_score),
    pct_below_500     = mean(nFeature_RNA < 500) * 100,
    pct_th_positive   = mean(th_positive) * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(mean_da_score))

cat("\n--- Per-cluster diagnostic (floor=200), sorted by DA score ---\n")
print(as.data.frame(summary_by_cluster), row.names = FALSE)

cat("\n--- Clusters that would be called 'Dopaminergic Neurons' (mean_da_score > 0) ---\n")
da_like <- summary_by_cluster[summary_by_cluster$mean_da_score > 0, ]
print(as.data.frame(da_like), row.names = FALSE)

# Flag: among DA-like clusters, is there one with low median nFeature AND
# low TH-positivity fraction -- i.e. a shadow cluster the same shape as the
# human diagnostic found (high module score via dropout-driven correlation,
# but not actually mostly TH+).
cat("\n--- Shadow-cluster check: DA-like clusters with median_nFeature < 500 ---\n")
shadow_candidates <- da_like[da_like$median_nFeature < 500, ]
if (nrow(shadow_candidates) > 0) {
  print(as.data.frame(shadow_candidates), row.names = FALSE)
} else {
  cat("None. No DA-scoring cluster has a median nFeature_RNA below 500.\n")
}

write.csv(summary_by_cluster, file.path(root, "results", "qc_floor_diagnostic_untreated.csv"), row.names = FALSE)
cat("\nDIAGNOSTIC COMPLETE\n")
