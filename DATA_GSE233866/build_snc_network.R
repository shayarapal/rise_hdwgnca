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
  h5seurat_path = "/shared/DATA_GSE233866/seurat_GSE233866_SNc_VTA.rds",
  out_dir       = "/shared/DATA_GSE233866/results/snc_network",
  cell_type_col = "cell_type",
  group_by      = c("cell_type", "sample_id"),
  group_name    = "Dopaminergic Neurons",
  wgcna_name    = "cacna1d_mouse_snc",
  condition_col = "region",
  ref_group     = "SNc"
)
p <- with_defaults(p)
dir.create(p$out_dir, recursive = TRUE, showWarnings = FALSE)

obj <- subset_condition(load_seurat(p$h5seurat_path), p$condition_col, p$ref_group)
log_mem("after load+subset(SNc)")
cat("n_cells:", ncol(obj), "\n")

obj <- setup_and_metacells(obj, p)
log_mem("after setup_and_metacells(SNc)")

obj <- test_soft_powers(obj, p)
log_mem("after test_soft_powers(SNc)")

# Pick soft power from THIS dataset's own recommendation -- the human
# dataset's soft_power=5 has no bearing on this mouse dataset.
rec <- jsonlite::fromJSON(file.path(p$out_dir, "soft_power_recommendation.json"))
soft_power <- if (!is.null(rec$recommended_power) && !is.na(rec$recommended_power)) {
  rec$recommended_power
} else {
  stop("No soft power reached the scale-free threshold for SNc -- inspect soft_power_table.csv")
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
log_mem("after ConstructNetwork(SNc)")

png(file.path(p$out_dir, "dendrogram.png"), width = 1400, height = 800, res = 150)
PlotDendrogram(obj, main = "Mouse SNc — Dopaminergic Neuron gene dendrogram", wgcna_name = p$wgcna_name)
dev.off()

obj <- ScaleData(obj, features = GetWGCNAGenes(obj, p$wgcna_name), verbose = FALSE)
log_mem("after ScaleData(SNc)")

obj <- ModuleEigengenes(obj, group.by.vars = p$group_by[[2]], wgcna_name = p$wgcna_name)
obj <- ModuleConnectivity(obj, group.by = p$cell_type_col, group_name = p$group_name, wgcna_name = p$wgcna_name)
log_mem("after ModuleEigengenes+Connectivity(SNc)")

write_module_outputs(obj, p)
log_mem("after write_module_outputs(SNc)")

plot_list_hme <- ModuleFeaturePlot(obj, wgcna_name = p$wgcna_name, features = "hMEs", order_points = TRUE, reduction = "umap")
if (length(plot_list_hme) > 0) {
  ggsave(file.path(p$out_dir, "module_eigengene_umap.png"), plot = wrap_plots(plot_list_hme, ncol = 4), width = 16, height = 12, dpi = 150)
}

pdf(file.path(p$out_dir, "module_correlogram.pdf"), width = 8, height = 7)
ModuleCorrelogram(obj, wgcna_name = p$wgcna_name)
dev.off()
png(file.path(p$out_dir, "module_correlogram.png"), width = 1000, height = 875, res = 125)
ModuleCorrelogram(obj, wgcna_name = p$wgcna_name)
dev.off()
log_mem("after ModuleCorrelogram(SNc)")

obj <- ModuleExprScore(obj, n_genes = 25, method = "UCell", wgcna_name = p$wgcna_name)
plot_list_scores <- ModuleFeaturePlot(obj, wgcna_name = p$wgcna_name, features = "scores", order_points = TRUE, ucell = TRUE, reduction = "umap")
if (length(plot_list_scores) > 0) {
  ggsave(file.path(p$out_dir, "module_scores_umap.png"), plot = wrap_plots(plot_list_scores, ncol = 4), width = 16, height = 12, dpi = 150)
}
log_mem("after ModuleExprScore+FeaturePlot(SNc)")

dbs <- c("GO_Biological_Process_2023", "KEGG_2021_Human")
obj <- RunEnrichr(obj, dbs = dbs, max_genes = 100, wgcna_name = p$wgcna_name)
enrich_df <- GetEnrichrTable(obj, wgcna_name = p$wgcna_name)
write.csv(enrich_df, file.path(p$out_dir, "enrichr_table.csv"), row.names = FALSE)
EnrichrBarPlot(obj, outdir = file.path(p$out_dir, "enrichr_plots"), n_terms = 10, plot_size = c(5, 7), logscale = TRUE, wgcna_name = p$wgcna_name)
log_mem("after RunEnrichr(SNc)")

cat("DIAGNOSTIC COMPLETE\n")
