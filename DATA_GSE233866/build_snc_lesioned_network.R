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
  h5seurat_path = "/shared/DATA_GSE233866/seurat_lesion_intact_SNc_VTA.rds",
  out_dir       = "/shared/DATA_GSE233866/results/snc_lesioned_network",
  cell_type_col = "cell_type",
  group_by      = c("cell_type", "sample_id"),
  group_name    = "Dopaminergic Neurons",
  wgcna_name    = "cacna1d_mouse_snc_lesioned",
  condition_col = "condition",
  ref_group     = "lesioned"
)
p <- with_defaults(p)
dir.create(p$out_dir, recursive = TRUE, showWarnings = FALSE)

# Restrict to SNc first, THEN to the lesioned condition within SNc.
obj_full <- load_seurat(p$h5seurat_path)
obj_snc  <- subset_condition(obj_full, "region", "SNc")
rm(obj_full); gc()
obj <- subset_condition(obj_snc, p$condition_col, p$ref_group)
rm(obj_snc); gc()
log_mem("after load+subset(SNc, lesioned)")
cat("n_cells:", ncol(obj), "\n")

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
PlotDendrogram(obj, main = "Mouse SNc (lesioned) — Dopaminergic Neuron gene dendrogram", wgcna_name = p$wgcna_name)
dev.off()

obj <- ScaleData(obj, features = GetWGCNAGenes(obj, p$wgcna_name), verbose = FALSE)
log_mem("after ScaleData")

obj <- ModuleEigengenes(obj, group.by.vars = p$group_by[[2]], wgcna_name = p$wgcna_name)
obj <- ModuleConnectivity(obj, group.by = p$cell_type_col, group_name = p$group_name, wgcna_name = p$wgcna_name)
log_mem("after ModuleEigengenes+Connectivity")

write_module_outputs(obj, p)
log_mem("after write_module_outputs")

panel_height <- function(n, ncol = 4, per_row_in = 3.6, min_in = 10) max(min_in, ceiling(n / ncol) * per_row_in)

plot_list_hme <- ModuleFeaturePlot(obj, wgcna_name = p$wgcna_name, features = "hMEs", order_points = TRUE, reduction = "umap")
if (length(plot_list_hme) > 0) {
  ggsave(file.path(p$out_dir, "module_eigengene_umap.png"), plot = wrap_plots(plot_list_hme, ncol = 4),
         width = 20, height = panel_height(length(plot_list_hme)), dpi = 250, limitsize = FALSE)
}

pdf(file.path(p$out_dir, "module_correlogram.pdf"), width = 10, height = 9)
ModuleCorrelogram(obj, wgcna_name = p$wgcna_name)
dev.off()
png(file.path(p$out_dir, "module_correlogram.png"), width = 2200, height = 1900, res = 200)
ModuleCorrelogram(obj, wgcna_name = p$wgcna_name)
dev.off()
log_mem("after ModuleCorrelogram")

obj <- ModuleExprScore(obj, n_genes = 25, method = "UCell", wgcna_name = p$wgcna_name)
plot_list_scores <- ModuleFeaturePlot(obj, wgcna_name = p$wgcna_name, features = "scores", order_points = TRUE, ucell = TRUE, reduction = "umap")
if (length(plot_list_scores) > 0) {
  ggsave(file.path(p$out_dir, "module_scores_umap.png"), plot = wrap_plots(plot_list_scores, ncol = 4),
         width = 20, height = panel_height(length(plot_list_scores)), dpi = 250, limitsize = FALSE)
}
log_mem("after ModuleExprScore+FeaturePlot")

hme_df <- as.data.frame(GetMEs(obj, harmonized = TRUE, wgcna_name = p$wgcna_name))
hme_df$sample_id <- obj$sample_id[rownames(hme_df)]
mod_cols <- setdiff(colnames(hme_df), "sample_id")
hme_by_sample <- aggregate(hme_df[mod_cols], by = list(sample_id = hme_df$sample_id), FUN = mean)
mat <- as.matrix(hme_by_sample[mod_cols])
rownames(mat) <- hme_by_sample$sample_id
mat_z <- scale(mat)
df_long <- as.data.frame(as.table(mat_z))
colnames(df_long) <- c("sample_id", "module", "z_hME")
p_heatmap <- ggplot(df_long, aes(x = module, y = sample_id, fill = z_hME)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "z(hME)") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Module", y = "Animal", title = "Mouse SNc (lesioned) — mean harmonized module eigengene per animal")
ggsave(file.path(p$out_dir, "module_by_animal_heatmap.png"), plot = p_heatmap, width = 14, height = 6, dpi = 250)
log_mem("after module_by_animal_heatmap")

dbs <- c("GO_Biological_Process_2023", "KEGG_2019_Mouse")
obj <- RunEnrichr(obj, dbs = dbs, max_genes = 100, wgcna_name = p$wgcna_name)
enrich_df <- GetEnrichrTable(obj, wgcna_name = p$wgcna_name)
write.csv(enrich_df, file.path(p$out_dir, "enrichr_table.csv"), row.names = FALSE)
EnrichrBarPlot(obj, outdir = file.path(p$out_dir, "enrichr_plots"), n_terms = 10, plot_size = c(5, 7), logscale = TRUE, wgcna_name = p$wgcna_name)
log_mem("after RunEnrichr")

cat("DIAGNOSTIC COMPLETE\n")
