log_mem <- function(label) {
  gc()
  vmrss_gb <- tryCatch({
    lines <- readLines("/proc/self/status")
    kb <- as.numeric(gsub("[^0-9]", "", grep("^VmRSS:", lines, value = TRUE)))
    kb / 1024 / 1024
  }, error = function(e) NA)
  cat(sprintf("[MEM %s] RSS=%.2f GB  label=%s\n", format(Sys.time(), "%H:%M:%S"), vmrss_gb, label))
  flush(stdout())
}

setwd("/app")
source("/app/plumber.R")
pipeline_libs()
log_mem("after pipeline_libs()")

# THE region x disease-state comparison: SNc vs VTA, restricted to LESIONED
# cells only from both regions -- "what's more active in the vulnerable
# region specifically during neurodegeneration," not baseline regional
# structure (that's snc_network/vta_network + the preservation test) and not
# disease-vs-healthy within one region (that's snc_dme). One combined
# network pooling both regions' lesioned cells, then FindDMEs by region.
p <- list(
  h5seurat_path = "/shared/DATA_GSE233866/seurat_lesion_intact_SNc_VTA.rds",
  out_dir       = "/shared/DATA_GSE233866/results/region_dme_lesioned",
  cell_type_col = "cell_type",
  group_by      = c("cell_type", "sample_id"),
  group_name    = "Dopaminergic Neurons",
  wgcna_name    = "cacna1d_mouse_region_lesioned"
)
p <- with_defaults(p)
dir.create(p$out_dir, recursive = TRUE, showWarnings = FALSE)

obj_full <- load_seurat(p$h5seurat_path)
obj_lesioned <- subset_condition(obj_full, "condition", "lesioned")
rm(obj_full); gc()
# subset_condition requires an exact single value; region has two real values
# (SNc/VTA) plus NA (non-DA or ambiguous) -- keep only the two real ones.
obj <- obj_lesioned[, obj_lesioned$region %in% c("SNc", "VTA")]
rm(obj_lesioned); gc()
log_mem("after load+subset(lesioned, SNc+VTA)")
cat("n_cells:", ncol(obj), "\n")
print(table(obj$region))

# group_name must match ONE cell_type value; here all cells are already
# Dopaminergic Neurons, so this is a no-op restriction as in the other
# GSE233866 scripts. group_by's second column (sample_id) is still the
# batch variable; region is the biological variable of interest.
obj <- setup_and_metacells(obj, p)
log_mem("after setup_and_metacells")

obj <- test_soft_powers(obj, p)
log_mem("after test_soft_powers")

rec <- jsonlite::fromJSON(file.path(p$out_dir, "soft_power_recommendation.json"))
soft_power <- if (!is.null(rec$recommended_power) && !is.na(rec$recommended_power)) {
  rec$recommended_power
} else {
  stop("No soft power reached the scale-free threshold -- inspect soft_power_table.csv")
}
cat("Using soft_power =", soft_power, "\n")
p$soft_power <- soft_power

obj <- ConstructNetwork(
  obj,
  soft_power    = p$soft_power,
  networkType   = p$network_type,
  tom_outdir    = p$out_dir,
  tom_name      = p$wgcna_name,
  overwrite_tom = TRUE,
  wgcna_name    = p$wgcna_name
)
log_mem("after ConstructNetwork")

png(file.path(p$out_dir, "dendrogram.png"), width = 2400, height = 1400, res = 200)
PlotDendrogram(obj, main = "Mouse SNc vs VTA, lesioned only — gene dendrogram", wgcna_name = p$wgcna_name)
dev.off()

obj <- ScaleData(obj, features = GetWGCNAGenes(obj, p$wgcna_name), verbose = FALSE)
log_mem("after ScaleData")

obj <- ModuleEigengenes(obj, group.by.vars = p$group_by[[2]], wgcna_name = p$wgcna_name)
obj <- ModuleConnectivity(obj, group.by = p$cell_type_col, group_name = p$group_name, wgcna_name = p$wgcna_name)
log_mem("after ModuleEigengenes+Connectivity")

write_module_outputs(obj, p)
log_mem("after write_module_outputs")

cat("--- CACNA1D module in combined lesioned-only (SNc+VTA) network ---\n")
mods <- GetModules(obj, wgcna_name = p$wgcna_name)
print(mods[mods$gene_name == "Cacna1d", ])

# -----------------------------------------------------------------
# DME: SNc vs VTA, within lesioned cells only -- what's more active in
# the vulnerable region specifically during neurodegeneration.
# -----------------------------------------------------------------
barcodes1 <- colnames(obj)[obj$region == "SNc"]
barcodes2 <- colnames(obj)[obj$region == "VTA"]
cat("barcodes1 (SNc, lesioned):", length(barcodes1), " barcodes2 (VTA, lesioned):", length(barcodes2), "\n")

DMEs <- FindDMEs(
  obj,
  barcodes1  = barcodes1,
  barcodes2  = barcodes2,
  features   = "MEs",
  harmonized = TRUE,
  wgcna_name = p$wgcna_name,
  test.use   = "wilcox"
)
write.csv(DMEs, file.path(p$out_dir, "dme_results.csv"), row.names = FALSE)
cat("--- DME results (all modules, SNc vs VTA within lesioned) ---\n")
print(DMEs)

png(file.path(p$out_dir, "dme_volcano.png"), width = 2000, height = 1600, res = 200)
print(PlotDMEsVolcano(obj, DMEs, wgcna_name = p$wgcna_name))
dev.off()

png(file.path(p$out_dir, "dme_lollipop.png"), width = 2000, height = 1600, res = 200)
print(PlotDMEsLollipop(obj, DMEs, wgcna_name = p$wgcna_name, pvalue = "p_val_adj"))
dev.off()
log_mem("after DME + plots")

cat("DIAGNOSTIC COMPLETE\n")
