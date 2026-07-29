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

p <- list(
  h5seurat_path = "/shared/seurat_GSE243639_SNc.rds",
  out_dir       = "/shared/results/pd_network_final",
  cell_type_col = "cell_type",
  group_by      = c("cell_type", "sample_id"),
  group_name    = "Dopaminergic Neurons",
  wgcna_name    = "cacna1d_pd",
  soft_power    = 5,
  condition_col = "condition",
  ref_group     = "PD"
)
p <- with_defaults(p)
dir.create(p$out_dir, recursive = TRUE, showWarnings = FALSE)

obj <- subset_condition(load_seurat(p$h5seurat_path), p$condition_col, p$ref_group)
log_mem("after load+subset(PD)")

obj <- setup_and_metacells(obj, p)
log_mem("after setup_and_metacells(PD)")

obj <- test_soft_powers(obj, p)
log_mem("after test_soft_powers(PD)")

obj <- ConstructNetwork(
  obj,
  soft_power    = p$soft_power,
  networkType   = p$network_type,
  tom_outdir    = p$out_dir,
  tom_name      = p$wgcna_name,
  overwrite_tom = TRUE,
  wgcna_name    = p$wgcna_name
)
log_mem("after ConstructNetwork(PD)")

png(file.path(p$out_dir, "dendrogram.png"), width = 1400, height = 800, res = 150)
PlotDendrogram(obj, main = "PD — Dopaminergic Neuron gene dendrogram", wgcna_name = p$wgcna_name)
dev.off()
log_mem("after PlotDendrogram(PD)")

obj <- ScaleData(obj, features = GetWGCNAGenes(obj, p$wgcna_name), verbose = FALSE)
log_mem("after ScaleData(PD, WGCNA genes)")

obj <- ModuleEigengenes(obj, group.by.vars = p$group_by[[2]], wgcna_name = p$wgcna_name)
log_mem("after ModuleEigengenes(PD)")

obj <- ModuleConnectivity(obj, group.by = p$cell_type_col, group_name = p$group_name, wgcna_name = p$wgcna_name)
log_mem("after ModuleConnectivity(PD)")

write_module_outputs(obj, p)
log_mem("after write_module_outputs(PD)")

plot_list_hme <- ModuleFeaturePlot(
  obj,
  wgcna_name   = p$wgcna_name,
  features     = "hMEs",
  order_points = TRUE,
  reduction    = "umap"
)
stopifnot(length(plot_list_hme) > 0)
p_hme <- wrap_plots(plot_list_hme, ncol = 4)
ggsave(file.path(p$out_dir, "module_eigengene_umap.png"), plot = p_hme, width = 16, height = 12, dpi = 150)
log_mem("after ModuleFeaturePlot(PD)")

# ModuleCorrelogram uses base R graphics -- pdf() device, not ggsave()
pdf(file.path(p$out_dir, "module_correlogram.pdf"), width = 8, height = 7)
ModuleCorrelogram(obj, wgcna_name = p$wgcna_name)
dev.off()
png(file.path(p$out_dir, "module_correlogram.png"), width = 1000, height = 875, res = 125)
ModuleCorrelogram(obj, wgcna_name = p$wgcna_name)
dev.off()
log_mem("after ModuleCorrelogram(PD)")

obj <- ModuleExprScore(
  obj,
  n_genes    = 25,
  method     = "UCell",
  wgcna_name = p$wgcna_name
)
log_mem("after ModuleExprScore(PD, UCell)")

plot_list_scores <- ModuleFeaturePlot(
  obj,
  wgcna_name   = p$wgcna_name,
  features     = "scores",
  order_points = TRUE,
  ucell        = TRUE,
  reduction    = "umap"
)
stopifnot(length(plot_list_scores) > 0)
p_scores <- wrap_plots(plot_list_scores, ncol = 4)
ggsave(file.path(p$out_dir, "module_scores_umap.png"), plot = p_scores, width = 16, height = 12, dpi = 150)
log_mem("after ModuleFeaturePlot(PD, scores)")

cat("DIAGNOSTIC COMPLETE\n")
