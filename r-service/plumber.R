library(plumber)
library(future)
library(promises)
library(uuid)
library(jsonlite)

# Null-coalescing: not in base R, avoid rlang dependency
`%||%` <- function(a, b) if (!is.null(a)) a else b

# load_seurat(): accepts .rds (Seurat v5 native) or .h5Seurat (legacy, needs shims).
# Exported to the future workers along with the pipeline helpers that call it.
source("seuratdisk_compat.R")

# In-memory job store. Lives in the main plumber process.
# Keys are UUIDs; values are lists with $status, $out_dir, $error.
jobs_env <- new.env(parent = emptyenv())

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

new_job <- function(out_dir) {
  job_id <- UUIDgenerate()
  jobs_env[[job_id]] <- list(status = "running", out_dir = out_dir, error = NULL)
  job_id
}

finish_job <- function(job_id) {
  jobs_env[[job_id]]$status <- "done"
}

fail_job <- function(job_id, err) {
  jobs_env[[job_id]]$status <- "failed"
  jobs_env[[job_id]]$error  <- conditionMessage(err)
}

require_job <- function(job_id, res) {
  if (!exists(job_id, envir = jobs_env, inherits = FALSE)) {
    res$status <- 404L
    stop(paste("job not found:", job_id))
  }
  jobs_env[[job_id]]
}

save_png <- function(plot_obj, path) {
  png(path, width = 1200, height = 800, res = 150)
  print(plot_obj)
  dev.off()
}

# ---------------------------------------------------------------------------
# Pipeline helpers (called inside futures — no reference to jobs_env here)
# ---------------------------------------------------------------------------

run_setup_through_soft_powers <- function(p) {
  library(Seurat)
  library(SeuratDisk)
  library(hdWGCNA)
  library(patchwork)

  dir.create(p$out_dir, recursive = TRUE, showWarnings = FALSE)

  obj <- load_seurat(p$h5seurat_path)   # .rds or .h5Seurat

  # gene_select/fraction per the hdWGCNA tutorial. Without them SelectNetworkGenes falls back
  # to gene_select = "variable", i.e. VariableFeatures() — genes chosen for variance ACROSS all
  # cell types, not genes expressed within the one being analysed.
  obj <- SetupForWGCNA(
    obj,
    gene_select = "fraction",
    fraction    = 0.05,
    wgcna_name  = p$wgcna_name
  )

  # MetacellsByGroups is kNN-based and takes no seed argument. Seed it here so /analyze rebuilds
  # the same metacells whose soft-power curve was inspected in /test-soft-powers — otherwise the
  # chosen soft power is applied to a different metacell set.
  set.seed(42)
  obj <- MetacellsByGroups(
    seurat_obj  = obj,
    group.by    = p$group_by,
    reduction   = "pca",
    k           = p$k,
    max_shared  = p$max_shared,
    ident.group = p$cell_type_col,
    wgcna_name  = p$wgcna_name
  )

  obj <- NormalizeMetacells(obj, wgcna_name = p$wgcna_name)

  obj <- SetDatExpr(
    obj,
    group_name = p$group_name,
    group.by   = p$cell_type_col,
    assay      = "RNA",
    layer      = "data",
    wgcna_name = p$wgcna_name
  )

  obj <- TestSoftPowers(obj,
    networkType = p$network_type,
    wgcna_name  = p$wgcna_name
  )

  sp_plot <- wrap_plots(PlotSoftPowers(obj), ncol = 2)
  save_png(sp_plot, file.path(p$out_dir, "soft_power_plot.png"))
}

run_full_pipeline <- function(p) {
  library(Seurat)
  library(SeuratDisk)
  library(hdWGCNA)
  library(patchwork)
  library(jsonlite)

  dir.create(p$out_dir, recursive = TRUE, showWarnings = FALSE)

  obj <- load_seurat(p$h5seurat_path)   # .rds or .h5Seurat

  # gene_select/fraction per the hdWGCNA tutorial. Without them SelectNetworkGenes falls back
  # to gene_select = "variable", i.e. VariableFeatures() — genes chosen for variance ACROSS all
  # cell types, not genes expressed within the one being analysed.
  obj <- SetupForWGCNA(
    obj,
    gene_select = "fraction",
    fraction    = 0.05,
    wgcna_name  = p$wgcna_name
  )

  # MetacellsByGroups is kNN-based and takes no seed argument. Seed it here so /analyze rebuilds
  # the same metacells whose soft-power curve was inspected in /test-soft-powers — otherwise the
  # chosen soft power is applied to a different metacell set.
  set.seed(42)
  obj <- MetacellsByGroups(
    seurat_obj  = obj,
    group.by    = p$group_by,
    reduction   = "pca",
    k           = p$k,
    max_shared  = p$max_shared,
    ident.group = p$cell_type_col,
    wgcna_name  = p$wgcna_name
  )

  obj <- NormalizeMetacells(obj, wgcna_name = p$wgcna_name)

  obj <- SetDatExpr(
    obj,
    group_name = p$group_name,
    group.by   = p$cell_type_col,
    assay      = "RNA",
    layer      = "data",
    wgcna_name = p$wgcna_name
  )

  obj <- TestSoftPowers(obj,
    networkType = p$network_type,
    wgcna_name  = p$wgcna_name
  )

  sp_plot <- wrap_plots(PlotSoftPowers(obj), ncol = 2)
  save_png(sp_plot, file.path(p$out_dir, "soft_power_plot.png"))

  # tom_outdir defaults to the relative path "TOM", which lands in the r-service working
  # directory rather than the job's out_dir — so concurrent jobs (workers = 2) sharing a
  # wgcna_name collide, and overwrite_tom = FALSE can trip over a stale TOM on a rerun.
  obj <- ConstructNetwork(
    obj,
    soft_power  = p$soft_power,
    networkType = p$network_type,
    tom_outdir  = p$out_dir,
    tom_name    = p$wgcna_name,
    wgcna_name  = p$wgcna_name
  )

  obj <- ModuleEigengenes(
    obj,
    group.by.vars = p$group_by[[2]],
    wgcna_name    = p$wgcna_name
  )

  obj <- ModuleConnectivity(
    obj,
    group.by   = p$cell_type_col,
    group_name = p$group_name,
    wgcna_name = p$wgcna_name
  )

  modules <- GetModules(obj, wgcna_name = p$wgcna_name)
  write.csv(modules, file.path(p$out_dir, "modules.csv"), row.names = FALSE)

  # PlotKMEs returns a ggplot of per-gene kME values coloured by module.
  # OPEN: if you prefer a connectivity graph, replace with ModuleNetworkPlot()
  # which writes per-module PNGs to a directory but has no single combined return.
  net_plot <- PlotKMEs(obj, ncol = 5)
  save_png(net_plot, file.path(p$out_dir, "network_plot.png"))
}

# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------

#* @post /test-soft-powers
#* @serializer unboxedJSON
function(req, res) {
  p <- req$body
  required <- c("h5seurat_path", "out_dir", "cell_type_col",
                 "group_by", "group_name", "wgcna_name")
  missing <- setdiff(required, names(p))
  if (length(missing)) {
    res$status <- 400L
    return(list(error = paste("missing fields:", paste(missing, collapse = ", "))))
  }

  p$k            <- p$k %||% 25L
  p$max_shared   <- p$max_shared %||% 15L
  p$network_type <- p$network_type %||% "signed"

  job_id <- new_job(p$out_dir)

  future_promise(run_setup_through_soft_powers(p)) %...>%
    (function(...) finish_job(job_id)) %...!%
    (function(err) fail_job(job_id, err))

  list(job_id = job_id)
}

#* @post /analyze
#* @serializer unboxedJSON
function(req, res) {
  p <- req$body
  required <- c("h5seurat_path", "out_dir", "cell_type_col",
                 "group_by", "group_name", "wgcna_name", "soft_power")
  missing <- setdiff(required, names(p))
  if (length(missing)) {
    res$status <- 400L
    return(list(error = paste("missing fields:", paste(missing, collapse = ", "))))
  }

  p$k            <- p$k %||% 25L
  p$max_shared   <- p$max_shared %||% 15L
  p$network_type <- p$network_type %||% "signed"

  job_id <- new_job(p$out_dir)

  future_promise(run_full_pipeline(p)) %...>%
    (function(...) finish_job(job_id)) %...!%
    (function(err) fail_job(job_id, err))

  list(job_id = job_id)
}

#* @get /status/<job_id>
#* @serializer unboxedJSON
function(job_id, res) {
  tryCatch({
    job <- require_job(job_id, res)
    out <- list(status = job$status)
    if (!is.null(job$error)) out$error <- job$error
    out
  }, error = function(e) {
    list(error = conditionMessage(e))
  })
}

#* @get /results/<job_id>/modules
#* @serializer unboxedJSON
function(job_id, res) {
  tryCatch({
    job  <- require_job(job_id, res)
    path <- file.path(job$out_dir, "modules.csv")
    if (!file.exists(path)) {
      res$status <- 404L
      return(list(error = "modules.csv not found; job may still be running"))
    }
    df <- read.csv(path, stringsAsFactors = FALSE)
    toJSON(df, dataframe = "rows", auto_unbox = TRUE)
  }, error = function(e) {
    res$status <- 500L
    list(error = conditionMessage(e))
  })
}

#* @get /results/<job_id>/plot
#* @serializer contentType list(type="image/png")
function(job_id, res) {
  tryCatch({
    job  <- require_job(job_id, res)
    path <- file.path(job$out_dir, "network_plot.png")
    if (!file.exists(path)) {
      res$status <- 404L
      return(charToRaw("network_plot.png not found"))
    }
    readBin(path, "raw", n = file.info(path)$size)
  }, error = function(e) {
    res$status <- 500L
    charToRaw(conditionMessage(e))
  })
}

#* @get /results/<job_id>/soft-power-plot
#* @serializer contentType list(type="image/png")
function(job_id, res) {
  tryCatch({
    job  <- require_job(job_id, res)
    path <- file.path(job$out_dir, "soft_power_plot.png")
    if (!file.exists(path)) {
      res$status <- 404L
      return(charToRaw("soft_power_plot.png not found"))
    }
    readBin(path, "raw", n = file.info(path)$size)
  }, error = function(e) {
    res$status <- 500L
    charToRaw(conditionMessage(e))
  })
}
