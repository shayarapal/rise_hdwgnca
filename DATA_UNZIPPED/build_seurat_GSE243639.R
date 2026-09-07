# =============================================================
# GSE243639 — human substantia nigra pars compacta (SNc) snRNA-seq
# 29 donors: 15 Parkinson's + 14 control. 10x v3. Goal: CACNA1D
# co-expression network, PD vs control, in SNc dopaminergic neurons.
# Run in RStudio with Seurat installed.
# =============================================================
library(Seurat)
library(dplyr)
# later: library(harmony); library(WGCNA); library(hdWGCNA)

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
inp  <- file.path(root, "SEURAT_INPUT")                 # 29 per-sample folders (barcodes/features/matrix)
meta <- read.csv(file.path(root, "GSE243639_metadata_clean.csv"), stringsAsFactors = FALSE)

# -----------------------------------------------------------------
# 1. LOAD each sample, tag donor-level clinical metadata, merge
# -----------------------------------------------------------------
# Read10X() cannot read these dirs: it expects either v2 (genes.tsv) or gzipped v3
# (*.tsv.gz). These are uncompressed v3, so it appends .gz and fails. ReadMtx takes
# explicit paths. features.tsv is the standard 3 columns; column 2 = gene symbol.
obj_list <- lapply(seq_len(nrow(meta)), function(i) {
  d <- file.path(inp, meta$sample_dir[i])
  m <- ReadMtx(
    mtx            = file.path(d, "matrix.mtx"),
    cells          = file.path(d, "barcodes.tsv"),
    features       = file.path(d, "features.tsv"),
    feature.column = 2
  )
  s <- CreateSeuratObject(counts = m, project = meta$sample_id[i],
                          min.cells = 3, min.features = 200)
  s$sample_id <- meta$sample_id[i]
  s$condition <- meta$condition[i]                       # PD / control
  s$age <- meta$age[i]; s$sex <- meta$sex[i]
  s$pmi <- meta$pmi_hours[i]; s$rin <- meta$rin[i]
  s$braak <- meta$braak[i]; s$cerad <- meta$cerad[i]
  RenameCells(s, add.cell.id = meta$sample_id[i])
})
# merge(x, y = list) once; Reduce() pairwise re-copies the growing object 28 times.
seu <- merge(obj_list[[1]], y = obj_list[-1])
seu <- JoinLayers(seu)                                   # Seurat v5: unify layers after merge

# -----------------------------------------------------------------
# 2. QC (single-nuclei: low mito expected; drop ambient/doublets)
# -----------------------------------------------------------------
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
# Floor raised 200 -> 750. At 200, low-complexity nuclei (median ~680 genes) formed a
# "shadow" cluster of each cell type; the dopaminergic shadow was 7,874 cells at ~55% TH+,
# where the missing TH is dropout, not absence. WGCNA is unusually sensitive to this:
# correlations across mostly-zero cells track detection depth rather than biology.
seu <- subset(seu, subset = nFeature_RNA > 750 & nFeature_RNA < 8000 & percent.mt < 5)

# -----------------------------------------------------------------
# 3. Normalize, reduce, correct DONOR batch, cluster
#    (batch correction matters: 29 donors, PD vs control confound)
# -----------------------------------------------------------------
seu <- NormalizeData(seu) |>
       FindVariableFeatures(nfeatures = 2000) |>
       ScaleData() |>
       RunPCA(npcs = 30)
# Harmony on: six clusters in the previous run were 82-100% a single donor. Those cells are
# high quality, so the QC floor above will not touch them -- only batch correction will.
# RunPCA must stay: harmony builds on it, and r-service/plumber.R hardcodes reduction="pca"
# for MetacellsByGroups, so the pca reduction has to survive into the saved object.
seu <- harmony::RunHarmony(seu, group.by.vars = "sample_id")
red <- "harmony"
seu <- RunUMAP(seu, reduction = red, dims = 1:30) |>
       FindNeighbors(reduction = red, dims = 1:30) |>
       FindClusters(resolution = 0.5)

# checkpoint: everything above is ~1h of compute. Relabelling below is seconds.
saveRDS(seu, file.path(root, "seurat_GSE243639_SNc_clustered.rds"))

# -----------------------------------------------------------------
# 4. LABEL CELL TYPES  ->  creates the `cell_type` column the GUI needs
# -----------------------------------------------------------------
# Assign labels from marker expression rather than a hand-written cluster->label map:
# resolution 0.5 yields far more clusters than there are cell types (one sample alone
# gives 13), several clusters map to the same type, and cluster numbering is not stable.
# AddModuleScore gives each cluster a background-corrected score per marker panel; the
# winning panel names the cluster. Same markers as the FeaturePlot below.
panels <- list(
  "Dopaminergic Neurons" = c("TH", "SLC6A3"),
  "Oligodendrocytes"     = c("MOBP", "PLP1"),
  "Astrocytes"           = c("AQP4", "GFAP"),
  "Microglia"            = c("CSF1R", "P2RY12"),
  "OPCs"                 = c("PDGFRA", "VCAN"),
  "GABAergic Neurons"    = c("GAD1", "GAD2"),
  "Endothelial"          = c("CLDN5")
)
panels <- lapply(panels, intersect, rownames(seu))
stopifnot(lengths(panels) > 0)

seu <- AddModuleScore(seu, features = panels, name = "panel_")
score_cols <- paste0("panel_", seq_along(panels))

# mean score per cluster x panel
scores <- sapply(score_cols, function(c) tapply(seu[[c]][, 1], Idents(seu), mean))
colnames(scores) <- names(panels)
write.csv(scores, file.path(root, "cluster_marker_scores.csv"))   # audit trail
print(round(scores, 3))

# A cluster whose best panel still scores <= 0 is not enriched for ANY of these markers
# (AddModuleScore is control-gene corrected, so 0 == background). Calling it a cell type
# anyway would be an invention -- leave it Unassigned rather than inflate a real label.
best  <- colnames(scores)[max.col(scores, ties.method = "first")]
best[apply(scores, 1, max) <= 0] <- "Unassigned"
names(best) <- rownames(scores)

seu <- RenameIdents(seu, best)
seu$cell_type <- Idents(seu)          # <-- THIS is the "CELL TYPE COLUMN" = "cell_type"

# sanity check before exporting:
print(table(seu$cell_type))            # confirm "Dopaminergic Neurons" exists and has cells
stopifnot("Dopaminergic Neurons" %in% levels(factor(seu$cell_type)))

# visual confirmation of the assignment above
p <- FeaturePlot(seu, unlist(panels, use.names = FALSE))
ggplot2::ggsave(file.path(root, "markers_featureplot.png"), p, width = 16, height = 12, dpi = 120)

# -----------------------------------------------------------------
# 5. SAVE — .rds  (upload THIS to the form; it now accepts .rds)
# -----------------------------------------------------------------
# NOTE: this stays a Seurat v5 object (Assay5). Make sure the r-service reading it
# back uses layer= (not the defunct slot=) when pulling counts/data.
saveRDS(seu, file.path(root, "seurat_GSE243639_SNc.rds"))   # <-- point the form's file input here

# -----------------------------------------------------------------
# 6. NEXT — CACNA1D network, PD vs control (this is what the GUI's soft-power
#    step begins; GROUP NAME = "Dopaminergic Neurons", CELL TYPE COLUMN = "cell_type"):
#   hdWGCNA: SetupForWGCNA -> MetacellsByGroups(group.by=c("condition","sample_id"))
#   Build a network per condition (SetDatExpr group_name="PD" / "control"),
#   locate CACNA1D's module in each, then compare with ModulePreservation +
#   differential correlation of CACNA1D partners. Enrich modules vs MitoCarta /
#   oxidative-stress / ubiquitination / PD gene sets.
#   (This dataset is SNc only — it gives PD-vs-control, not SNc-vs-VTA.)
# -----------------------------------------------------------------
cat("Done. Saved seurat_GSE243639_SNc.rds\n")
