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

# Human mirror of run_dme_snc.R: ONE combined network across both conditions
# (PD + control), then FindDMEs compares hMEs between condition-defined cell
# subsets within that same network -- unlike the earlier separate
# build_control_network.R / build_pd_network.R runs (different module
# definitions each, not directly comparable gene-for-gene).
#
# Efficiency note: earlier human runs subset by condition only (~80-90k
# cells/condition, all cell types -- DA-neuron restriction happened only
# implicitly inside MetacellsByGroups), which is why they needed 25-34GB.
# Pre-filtering to cell_type="Dopaminergic Neurons" FIRST here cuts this to
# ~10,388 cells total (both conditions combined) before any heavy step.
p <- list(
  h5seurat_path = "/shared/seurat_GSE243639_SNc.rds",
  out_dir       = "/shared/results/human_snc_dme",
  cell_type_col = "cell_type",
  group_by      = c("cell_type", "sample_id"),
  group_name    = "Dopaminergic Neurons",
  wgcna_name    = "cacna1d_human_combined"
)
p <- with_defaults(p)
dir.create(p$out_dir, recursive = TRUE, showWarnings = FALSE)

obj <- subset_condition(load_seurat(p$h5seurat_path), "cell_type", "Dopaminergic Neurons")
log_mem("after load+subset(Dopaminergic Neurons, both conditions)")
cat("n_cells:", ncol(obj), "\n")
print(table(obj$condition))

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
PlotDendrogram(obj, main = "Human SNc (PD + control combined) — gene dendrogram", wgcna_name = p$wgcna_name)
dev.off()

obj <- ScaleData(obj, features = GetWGCNAGenes(obj, p$wgcna_name), verbose = FALSE)
log_mem("after ScaleData")

obj <- ModuleEigengenes(obj, group.by.vars = p$group_by[[2]], wgcna_name = p$wgcna_name)
obj <- ModuleConnectivity(obj, group.by = p$cell_type_col, group_name = p$group_name, wgcna_name = p$wgcna_name)
log_mem("after ModuleEigengenes+Connectivity")

write_module_outputs(obj, p)
log_mem("after write_module_outputs")

cat("--- CACNA1D module in combined human network ---\n")
mods <- GetModules(obj, wgcna_name = p$wgcna_name)
print(mods[mods$gene_name == "CACNA1D", ])

# -----------------------------------------------------------------
# DME: PD vs control, within this one combined human SNc network
# -----------------------------------------------------------------
barcodes1 <- colnames(obj)[obj$condition == "PD"]
barcodes2 <- colnames(obj)[obj$condition == "control"]
cat("barcodes1 (PD):", length(barcodes1), " barcodes2 (control):", length(barcodes2), "\n")

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
cat("--- DME results (all modules) ---\n")
print(DMEs)

png(file.path(p$out_dir, "dme_volcano.png"), width = 2000, height = 1600, res = 200)
print(PlotDMEsVolcano(obj, DMEs, wgcna_name = p$wgcna_name))
dev.off()

png(file.path(p$out_dir, "dme_lollipop.png"), width = 2000, height = 1600, res = 200)
print(PlotDMEsLollipop(obj, DMEs, wgcna_name = p$wgcna_name, pvalue = "p_val_adj"))
dev.off()
log_mem("after DME + plots")

cat("DIAGNOSTIC COMPLETE\n")
