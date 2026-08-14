setwd("/app")
source("/app/plumber.R")
pipeline_libs()

# Pulls per-cell UMAP coordinates + region/condition labels + CACNA1D expression
# from the two already-built joint Seurat objects (UMAP was computed once on the
# FULL population in build_seurat_GSE233866.R / build_seurat_lesion_intact.R,
# before any SNc/VTA or intact/lesioned subsetting -- so this is a real shared
# embedding per cohort, not six independently-computed ones stitched together).
# Output is a flat CSV for a combined 3D visualization; no network is rebuilt.

extract_one <- function(path, condition_label, has_condition_col) {
  obj <- load_seurat(path)
  emb <- Embeddings(obj, "umap")
  df <- data.frame(
    cell_id  = colnames(obj),
    umap_1   = emb[, 1],
    umap_2   = emb[, 2],
    region   = obj$region,
    sample_id = obj$sample_id,
    cacna1d  = as.numeric(GetAssayData(obj, layer = "data")["Cacna1d", ])
  )
  df$condition <- if (has_condition_col) obj$condition else condition_label
  df <- df[!is.na(df$region) & df$region %in% c("SNc", "VTA"), ]
  df
}

df_healthy <- extract_one(
  "/shared/DATA_GSE233866/seurat_GSE233866_SNc_VTA.rds",
  condition_label = "healthy",
  has_condition_col = FALSE
)
cat("healthy:", nrow(df_healthy), "cells\n")

df_lesion_arm <- extract_one(
  "/shared/DATA_GSE233866/seurat_lesion_intact_SNc_VTA.rds",
  condition_label = NA,
  has_condition_col = TRUE
)
cat("lesion arm (intact+lesioned):", nrow(df_lesion_arm), "cells\n")

combined <- rbind(df_healthy, df_lesion_arm)
combined$group <- paste(combined$region, combined$condition, sep = "_")
print(table(combined$group))

out_path <- "/shared/DATA_GSE233866/results/umap_combined_all_networks.csv"
write.csv(combined, out_path, row.names = FALSE)
cat("Wrote", nrow(combined), "rows to", out_path, "\n")
