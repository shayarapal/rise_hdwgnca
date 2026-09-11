# =============================================================
# GSE233866 -- mouse midbrain snRNA-seq (untreated), Yaghmaeian Salmani et al.
# 6 animals (DatCre;TRAP, C57BL6, P90/P530). Goal: split into SNc vs VTA
# dopaminergic neurons (this dataset ships no per-cell region label -- the
# paper itself derives the split from near-mutually-exclusive Sox6/Calb1
# expression, which we replicate here) for a same-species, same-study
# CACNA1D network comparison against the human SNc data already analyzed.
# =============================================================
library(Seurat)
library(dplyr)
library(data.table)
library(Matrix)
library(harmony)

root <- "/shared/DATA_GSE233866"
raw  <- file.path(root, "GSE233866_untreated_counts.csv.gz")

# -----------------------------------------------------------------
# 1. LOAD -- dense gene x cell CSV. The header line has one FEWER field
#    than the data rows (no name for the gene-symbol column), which
#    fread's header=TRUE auto-detection mis-parses: it silently discards
#    the real header and uses the first DATA row (gene "Xkr4"'s counts)
#    as column names instead, losing that gene and every cell barcode.
#    Read the header explicitly and skip it when reading the data instead
#    of letting fread guess.
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

# cell barcodes are "s<sample>_<barcode>-1" -- the sample id IS the batch
# variable (6 animals), analogous to sample_id in the human build script
sample_id <- sub("_[ACGT]+-1$", "", colnames(counts))

seu <- CreateSeuratObject(counts = counts, project = "GSE233866",
                          min.cells = 3, min.features = 200)
seu$sample_id <- sample_id
rm(counts); gc()

# -----------------------------------------------------------------
# 2. QC -- mouse mitochondrial genes use the lowercase "mt-" prefix
# -----------------------------------------------------------------
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
cat("Pre-QC cells:", ncol(seu), "\n")
print(summary(seu$nFeature_RNA))
print(summary(seu$percent.mt))

# Same rationale as the human build: WGCNA is sensitive to detection-depth-
# driven correlations in low-complexity nuclei, so floor above the raw
# 200-gene CreateSeuratObject cutoff.
seu <- subset(seu, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 5)
cat("Post-QC cells:", ncol(seu), "\n")

# -----------------------------------------------------------------
# 3. Normalize, reduce, correct ANIMAL batch, cluster
# -----------------------------------------------------------------
seu <- NormalizeData(seu) |>
       FindVariableFeatures(nfeatures = 2000) |>
       ScaleData() |>
       RunPCA(npcs = 30)
seu <- harmony::RunHarmony(seu, group.by.vars = "sample_id")
red <- "harmony"
seu <- RunUMAP(seu, reduction = red, dims = 1:30) |>
       FindNeighbors(reduction = red, dims = 1:30) |>
       FindClusters(resolution = 0.5)

saveRDS(seu, file.path(root, "seurat_GSE233866_clustered.rds"))

# -----------------------------------------------------------------
# 4. IDENTIFY DOPAMINERGIC NEURONS (this dataset is mDA + non-mDA mixed,
#    per its own SRA metadata -- cell_type: "mDA, non-mDA")
# -----------------------------------------------------------------
da_markers <- intersect(c("Th", "Slc6a3", "Ddc", "Slc18a2"), rownames(seu))
stopifnot(length(da_markers) > 0)
seu <- AddModuleScore(seu, features = list(da_markers), name = "da_score_")
da_col <- grep("^da_score_", colnames(seu[[]]), value = TRUE)[1]

# mean DA-marker score per cluster; a cluster is "Dopaminergic Neurons" if its
# mean score clears background (0, since AddModuleScore is control-gene
# corrected) -- mirrors the human build's "don't invent a label below
# background" rule.
da_by_cluster <- tapply(seu[[da_col]][, 1], Idents(seu), mean)
print(round(sort(da_by_cluster, decreasing = TRUE), 3))
da_clusters <- names(da_by_cluster)[da_by_cluster > 0]
stopifnot(length(da_clusters) > 0)

seu$cell_type <- ifelse(as.character(Idents(seu)) %in% da_clusters,
                         "Dopaminergic Neurons", "Non-DA")
print(table(seu$cell_type))

# -----------------------------------------------------------------
# 5. SPLIT DA NEURONS INTO SNc vs VTA via Sox6 (SNc) / Calb1 (VTA)
#    -- replicates the paper's own "near-mutually-exclusive expression"
#    logic; this dataset provides no pre-computed region label.
# -----------------------------------------------------------------
stopifnot(all(c("Sox6", "Calb1") %in% rownames(seu)))
da_cells <- colnames(seu)[seu$cell_type == "Dopaminergic Neurons"]
sox6  <- FetchData(seu, vars = "Sox6",  cells = da_cells)[, 1]
calb1 <- FetchData(seu, vars = "Calb1", cells = da_cells)[, 1]

region <- rep(NA_character_, length(da_cells))
region[sox6 > 0  & calb1 == 0] <- "SNc"
region[calb1 > 0 & sox6  == 0] <- "VTA"
# cells expressing both or neither are genuinely ambiguous under this rule --
# leave as NA rather than force an assignment (matches the "don't invent a
# label" principle used for cell-type calling in the human build).
names(region) <- da_cells

seu$region <- NA_character_
seu@meta.data[da_cells, "region"] <- region
print(table(seu$region, useNA = "always"))
stopifnot(sum(!is.na(seu$region)) > 0)
stopifnot("SNc" %in% seu$region, "VTA" %in% seu$region)

saveRDS(seu, file.path(root, "seurat_GSE233866_SNc_VTA.rds"))
cat("Done. Saved seurat_GSE233866_SNc_VTA.rds\n")
