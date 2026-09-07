# One-off relabel: cluster 7 -> "Unassigned".
#
# The marker-scoring step assigned cluster 7 to "Dopaminergic Neurons" because it won the
# TH/SLC6A3 panel -- but it won weakly (score 0.329 vs cluster 4's 1.961) and only ~40% of
# its nuclei express TH/SLC6A3, against 80/88% in cluster 4. It also carries ~5x cluster 4's
# GAD1/GAD2 expression and the highest library complexity in the dataset, i.e. it is a mixed
# neuronal population (DA + GABAergic, plus likely doublets), not a clean DA cluster.
# Leaving it in would build the CACNA1D network on ~30% non-dopaminergic cells.
#
# NOTE: rerunning build_seurat_GSE243639.R will re-introduce this -- the scoring rule only
# requires score > 0, which cluster 7 passes. Re-run this script after any rebuild.
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
f    <- file.path(root, "seurat_GSE243639_SNc.rds")
seu  <- readRDS(f)

stopifnot("cell_type" %in% colnames(seu[[]]))
target <- "7"
stopifnot(target %in% levels(factor(seu$seurat_clusters)))

before <- table(seu$cell_type)
ct <- as.character(seu$cell_type)
ct[as.character(seu$seurat_clusters) == target] <- "Unassigned"
seu$cell_type <- factor(ct)
Idents(seu) <- "cell_type"

cat("=== cell_type before -> after\n")
print(cbind(before = as.integer(before[names(before)]),
            after  = as.integer(table(seu$cell_type)[names(before)])) |>
      `rownames<-`(names(before)))

da <- sum(seu$cell_type == "Dopaminergic Neurons")
cat(sprintf("\nDopaminergic Neurons now: %d cells (%.1f%% of dataset)\n", da, 100 * da / ncol(seu)))
stopifnot(da > 0)

saveRDS(seu, f)
cat("Saved", f, "\n")
