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

# Mirrors run_dme_region_intact.R / run_dme_region_lesioned.R, but on the
# healthy-baseline (untreated) cohort instead of the lesion arm -- no
# "condition" column exists here, so there is no condition-subsetting step.
# This is the formal DME counterpart the project was missing for the healthy
# state: healthy_baseline/snc_vs_vta_preservation/ tests whether SNc's
# MODULES structurally reproduce in VTA (ModulePreservation); this tests
# whether a module's ACTIVITY LEVEL differs between SNc and VTA (FindDMEs)
# -- a different statistic, not previously computed for this cohort.
p <- list(
  h5seurat_path = "/shared/DATA_GSE233866/seurat_GSE233866_SNc_VTA.rds",
  out_dir       = "/shared/DATA_GSE233866/results/region_dme_healthy",
  cell_type_col = "cell_type",
  group_by      = c("cell_type", "sample_id"),
  group_name    = "Dopaminergic Neurons",
  wgcna_name    = "cacna1d_mouse_region_healthy"
)
p <- with_defaults(p)
dir.create(p$out_dir, recursive = TRUE, showWarnings = FALSE)

obj_full <- load_seurat(p$h5seurat_path)
# subset_condition requires an exact single value; region has two real values
# (SNc/VTA) plus NA (non-DA or ambiguous) -- keep only the two real ones.
obj <- obj_full[, obj_full$region %in% c("SNc", "VTA")]
rm(obj_full); gc()
log_mem("after load+subset(healthy, SNc+VTA)")
cat("n_cells:", ncol(obj), "\n")
print(table(obj$region))

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
PlotDendrogram(obj, main = "Mouse SNc vs VTA, healthy baseline — gene dendrogram", wgcna_name = p$wgcna_name)
dev.off()

obj <- ScaleData(obj, features = GetWGCNAGenes(obj, p$wgcna_name), verbose = FALSE)
log_mem("after ScaleData")

obj <- ModuleEigengenes(obj, group.by.vars = p$group_by[[2]], wgcna_name = p$wgcna_name)
obj <- ModuleConnectivity(obj, group.by = p$cell_type_col, group_name = p$group_name, wgcna_name = p$wgcna_name)
log_mem("after ModuleEigengenes+Connectivity")

write_module_outputs(obj, p)
log_mem("after write_module_outputs")

cat("--- CACNA1D module in combined healthy-baseline (SNc+VTA) network ---\n")
mods <- GetModules(obj, wgcna_name = p$wgcna_name)
print(mods[mods$gene_name == "Cacna1d", ])

# -----------------------------------------------------------------
# DME: SNc vs VTA, within healthy-baseline cells -- is there already a
# regional asymmetry with no disease model involved at all.
# -----------------------------------------------------------------
barcodes1 <- colnames(obj)[obj$region == "SNc"]
barcodes2 <- colnames(obj)[obj$region == "VTA"]
cat("barcodes1 (SNc, healthy):", length(barcodes1), " barcodes2 (VTA, healthy):", length(barcodes2), "\n")

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
cat("--- DME results (all modules, SNc vs VTA within healthy baseline) ---\n")
print(DMEs)

png(file.path(p$out_dir, "dme_volcano.png"), width = 2000, height = 1600, res = 200)
print(PlotDMEsVolcano(obj, DMEs, wgcna_name = p$wgcna_name))
dev.off()

png(file.path(p$out_dir, "dme_lollipop.png"), width = 2000, height = 1600, res = 200)
print(PlotDMEsLollipop(obj, DMEs, wgcna_name = p$wgcna_name, pvalue = "p_val_adj"))
dev.off()
log_mem("after DME + plots")

cat("DIAGNOSTIC COMPLETE\n")
