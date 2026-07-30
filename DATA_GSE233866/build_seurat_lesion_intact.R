# =============================================================
# GSE233866 -- mouse midbrain snRNA-seq, LESIONED vs INTACT hemisphere.
# 6 animals, each contributing both a 6-OHDA-lesioned ("L") and the
# contralateral non-lesioned ("C"/intact) hemisphere -- a within-animal
# paired design. This is the mouse disease-model equivalent of the human
# PD-vs-control axis: sample prefix encodes BOTH the animal (numeric id,
# shared between its L/C pair -- the true batch/donor variable) and the
# condition (L=lesioned / C=intact -- the biological variable of interest,
# which must NOT be harmony-corrected away, exactly as condition=PD/control
# was never batch-corrected in the human build).
# =============================================================
library(Seurat)
library(dplyr)
library(data.table)
library(Matrix)
library(harmony)

root <- "/shared/DATA_GSE233866"
raw  <- file.path(root, "GSE233866_lesion_intact_counts.csv.gz")

# -----------------------------------------------------------------
# 1. LOAD -- same header-off-by-one quirk as the untreated file (see
#    build_seurat_GSE233866.R): read the header line separately rather
#    than trusting fread's header=TRUE autodetection.
# -----------------------------------------------------------------
header_line <- readLines(gzfile(raw), n = 1)
cell_barcodes <- strsplit(header_line, ",", fixed = TRUE)[[1]]

dt <- fread(raw, header = FALSE, skip = 1)
stopifnot(ncol(dt) == length(cell_barcodes) + 1)
gene_names <- dt[[1]]
counts <- as.matrix(dt[, -1, with = FALSE])
rownames(counts) <- gene_names
colnames(counts) <- cell_barcodes
rm(dt); gc()

stopifnot(sum(duplicated(gene_names)) == 0, sum(duplicated(cell_barcodes)) == 0)
counts <- Matrix(counts, sparse = TRUE)

# barcodes are "<condition><animal>_<barcode>-1", e.g. "L735_..." / "C735_..."
prefix    <- sub("_[ACGTN]+-1$", "", colnames(counts))
condition <- ifelse(grepl("^L", prefix), "lesioned", "intact")
animal_id <- sub("^[LC]", "", prefix)

seu <- CreateSeuratObject(counts = counts, project = "GSE233866_lesion",
                          min.cells = 3, min.features = 200)
seu$sample_id <- animal_id       # batch variable: 6 animals
seu$condition <- condition       # biological variable: lesioned vs intact
rm(counts); gc()

cat("condition counts (pre-QC):\n"); print(table(seu$condition))
cat("animal x condition (pre-QC):\n"); print(table(seu$sample_id, seu$condition))

# -----------------------------------------------------------------
# 2. QC
# -----------------------------------------------------------------
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
cat("Pre-QC cells:", ncol(seu), "\n")
print(summary(seu$nFeature_RNA))
print(summary(seu$percent.mt))

seu <- subset(seu, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 5)
cat("Post-QC cells:", ncol(seu), "\n")

# -----------------------------------------------------------------
# 3. Normalize, reduce, correct ANIMAL batch (not condition), cluster
# -----------------------------------------------------------------
seu <- NormalizeData(seu) |>
       FindVariableFeatures(nfeatures = 2000) |>
       ScaleData() |>
       RunPCA(npcs = 30)
seu <- harmony::RunHarmony(seu, group.by.vars = "sample_id")
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30) |>
       FindNeighbors(reduction = "harmony", dims = 1:30) |>
       FindClusters(resolution = 0.5)

saveRDS(seu, file.path(root, "seurat_lesion_intact_clustered.rds"))

# -----------------------------------------------------------------
# 4. IDENTIFY DOPAMINERGIC NEURONS
# -----------------------------------------------------------------
da_markers <- intersect(c("Th", "Slc6a3", "Ddc", "Slc18a2"), rownames(seu))
stopifnot(length(da_markers) > 0)
seu <- AddModuleScore(seu, features = list(da_markers), name = "da_score_")
da_col <- grep("^da_score_", colnames(seu[[]]), value = TRUE)[1]

da_by_cluster <- tapply(seu[[da_col]][, 1], Idents(seu), mean)
print(round(sort(da_by_cluster, decreasing = TRUE), 3))
da_clusters <- names(da_by_cluster)[da_by_cluster > 0]
stopifnot(length(da_clusters) > 0)

seu$cell_type <- ifelse(as.character(Idents(seu)) %in% da_clusters,
                         "Dopaminergic Neurons", "Non-DA")
print(table(seu$cell_type, seu$condition))

# -----------------------------------------------------------------
# 5. SPLIT DA NEURONS INTO SNc vs VTA via Sox6 (SNc) / Calb1 (VTA)
#    -- same rule as the untreated build, applied independently here
#    since lesioning could plausibly shift these markers' expression.
# -----------------------------------------------------------------
stopifnot(all(c("Sox6", "Calb1") %in% rownames(seu)))
da_cells <- colnames(seu)[seu$cell_type == "Dopaminergic Neurons"]
sox6  <- FetchData(seu, vars = "Sox6",  cells = da_cells)[, 1]
calb1 <- FetchData(seu, vars = "Calb1", cells = da_cells)[, 1]

region <- rep(NA_character_, length(da_cells))
region[sox6 > 0  & calb1 == 0] <- "SNc"
region[calb1 > 0 & sox6  == 0] <- "VTA"
names(region) <- da_cells

seu$region <- NA_character_
seu@meta.data[da_cells, "region"] <- region
print(table(seu$region, seu$condition, useNA = "always"))
stopifnot(sum(!is.na(seu$region)) > 0)

saveRDS(seu, file.path(root, "seurat_lesion_intact_SNc_VTA.rds"))
cat("Done. Saved seurat_lesion_intact_SNc_VTA.rds\n")
