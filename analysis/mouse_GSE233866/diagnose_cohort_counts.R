library(Seurat)
library(data.table)
library(Matrix)

load_raw <- function(raw_path) {
  header_line <- readLines(gzfile(raw_path), n = 1)
  cell_barcodes <- strsplit(header_line, ",", fixed = TRUE)[[1]]
  dt <- fread(raw_path, header = FALSE, skip = 1)
  gene_names <- dt[[1]]
  counts <- as.matrix(dt[, -1, with = FALSE])
  rownames(counts) <- gene_names
  colnames(counts) <- cell_barcodes
  rm(dt); gc()
  Matrix(counts, sparse = TRUE)
}

cat("========== UNTREATED (healthy baseline) ==========\n")
counts_u <- load_raw("/shared/DATA_GSE233866/GSE233866_untreated_counts.csv.gz")
cat("Raw barcodes in deposited file:", ncol(counts_u), "\n")
sample_id_u <- sub("_[ACGT]+-1$", "", colnames(counts_u))
cat("Animals represented (raw):", length(unique(sample_id_u)), "\n")
print(table(sample_id_u))

seu_u <- CreateSeuratObject(counts = counts_u, project = "u", min.cells = 3, min.features = 200)
seu_u$sample_id <- sub("_[ACGT]+-1$", "", colnames(seu_u))
cat("After CreateSeuratObject(min.cells=3, min.features=200):", ncol(seu_u), "\n")
rm(counts_u); gc()

seu_u[["percent.mt"]] <- PercentageFeatureSet(seu_u, pattern = "^mt-")
seu_u_qc <- subset(seu_u, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 5)
cat("After real QC (nFeature 500-8000, mt<5):", ncol(seu_u_qc), "\n")
cat("Per-animal counts after QC:\n")
print(table(seu_u_qc$sample_id))
rm(seu_u, seu_u_qc); gc()

cat("\n========== LESION ARM (intact + lesioned) ==========\n")
counts_l <- load_raw("/shared/DATA_GSE233866/GSE233866_lesion_intact_counts.csv.gz")
cat("Raw barcodes in deposited file:", ncol(counts_l), "\n")
prefix_l <- sub("_[ACGTN]+-1$", "", colnames(counts_l))
condition_l <- ifelse(grepl("^L", prefix_l), "lesioned", "intact")
animal_l <- sub("^[LC]", "", prefix_l)
cat("Animals represented (raw):", length(unique(animal_l)), "\n")
cat("Raw condition counts:\n")
print(table(condition_l))
cat("Raw animal x condition:\n")
print(table(animal_l, condition_l))

seu_l <- CreateSeuratObject(counts = counts_l, project = "l", min.cells = 3, min.features = 200)
seu_l$sample_id <- sub("^[LC]", "", sub("_[ACGTN]+-1$", "", colnames(seu_l)))
seu_l$condition <- ifelse(grepl("^L", sub("_[ACGTN]+-1$", "", colnames(seu_l))), "lesioned", "intact")
cat("After CreateSeuratObject(min.cells=3, min.features=200):", ncol(seu_l), "\n")
rm(counts_l); gc()

seu_l[["percent.mt"]] <- PercentageFeatureSet(seu_l, pattern = "^mt-")
seu_l_qc <- subset(seu_l, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 5)
cat("After real QC (nFeature 500-8000, mt<5):", ncol(seu_l_qc), "\n")
cat("Post-QC condition counts:\n")
print(table(seu_l_qc$condition))
cat("Post-QC animal x condition:\n")
print(table(seu_l_qc$sample_id, seu_l_qc$condition))

cat("\nDIAGNOSTIC COMPLETE\n")
