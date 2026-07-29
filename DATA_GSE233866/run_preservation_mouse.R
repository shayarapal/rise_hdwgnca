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
  h5seurat_path  = "/shared/DATA_GSE233866/seurat_GSE233866_SNc_VTA.rds",
  out_dir        = "/shared/DATA_GSE233866/results/snc_vs_vta_preservation",
  cell_type_col  = "cell_type",
  group_by       = c("cell_type", "sample_id"),
  group_name     = "Dopaminergic Neurons",
  wgcna_name     = "cacna1d_mouse_preservation",
  soft_power     = 5,
  condition_col  = "region",
  ref_group      = "SNc",
  query_group    = "VTA",
  n_permutations = 200
)
p <- with_defaults(p)
p$preservation_name <- paste0(p$ref_group, "-vs-", p$query_group)
dir.create(p$out_dir, recursive = TRUE, showWarnings = FALSE)

run_module_preservation(p)
log_mem("after run_module_preservation")

cat("DIAGNOSTIC COMPLETE\n")
