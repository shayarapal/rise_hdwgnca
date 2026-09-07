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
  png(path, width = 2400, height = 1600, res = 200)
  print(plot_obj)
  dev.off()
}

# Serve a file from a job's out_dir as JSON or PNG. Both result serializers 404
# identically when the job is still running and the file does not exist yet.
serve_csv <- function(job_id, res, filename) {
  tryCatch({
    job  <- require_job(job_id, res)
    path <- file.path(job$out_dir, filename)
    if (!file.exists(path)) {
      res$status <- 404L
      return(list(error = paste(filename, "not found; job may still be running")))
    }
    # Return the data.frame and let the serializer encode it. Calling toJSON() here
    # instead double-encodes: the serializer JSON-encodes the already-JSON string, and
    # the client receives a quoted string rather than an array of row objects.
    read.csv(path, stringsAsFactors = FALSE)
  }, error = function(e) {
    res$status <- 500L
    list(error = conditionMessage(e))
  })
}

# The power table plus both recommendations, as one object.
#
# The recommendations are recomputed from the CSV rather than read back from the JSON the
# pipeline wrote: fromJSON() turns a JSON null into an R NULL, which the serializer then emits
# as {} — so a round-trip would silently destroy exactly the "no power qualifies" signal the
# caller needs. The CSV is the single source of truth and recommend_power() is pure.
#
# read.csv(check.names = FALSE): the response keys must be exactly GetPowerTable()'s column
# names — SFT.R.sq, mean.k., median.k., max.k. — because the frontend indexes them by string.
# R's default name mangling would rewrite those dots and silently break the client.
serve_soft_powers <- function(job_id, res) {
  tryCatch({
    job   <- require_job(job_id, res)
    table <- file.path(job$out_dir, "soft_power_table.csv")
    if (!file.exists(table)) {
      res$status <- 404L
      return(list(error = "soft power table not found; job may still be running"))
    }
    pt  <- read.csv(table, stringsAsFactors = FALSE, check.names = FALSE)
    out <- recommend_power(pt)
    out$table <- pt
    out
  }, error = function(e) {
    res$status <- 500L
    list(error = conditionMessage(e))
  })
}

serve_png <- function(job_id, res, filename) {
  tryCatch({
    job  <- require_job(job_id, res)
    path <- file.path(job$out_dir, filename)
    if (!file.exists(path)) {
      res$status <- 404L
      return(charToRaw(paste(filename, "not found")))
    }
    readBin(path, "raw", n = file.info(path)$size)
  }, error = function(e) {
    res$status <- 500L
    charToRaw(conditionMessage(e))
  })
}

# Shared request validation + defaults for all POST endpoints.
BASE_FIELDS <- c("h5seurat_path", "out_dir", "cell_type_col",
                 "group_by", "group_name", "wgcna_name")

# The only values SelectNetworkGenes accepts; it stop()s on anything else.
GENE_SELECT_VALUES <- c("variable", "fraction", "all", "custom")
DEFAULT_FRACTIONS  <- c(0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25, 0.30)
SFT_THRESHOLD      <- 0.8   # WGCNA scale-free topology fit convention
SFT_MIN_POWER      <- 3     # ConstructNetwork()'s own min_power default

with_defaults <- function(p) {
  p$k            <- p$k %||% 25L
  p$max_shared   <- p$max_shared %||% 15L
  p$network_type <- p$network_type %||% "signed"
  p$gene_select  <- p$gene_select %||% "fraction"
  p$fraction     <- p$fraction %||% 0.05
  p
}

# Reject a bad gene_select at the HTTP boundary. Otherwise SelectNetworkGenes stop()s deep
# inside a future worker and the caller gets it as an opaque job failure minutes later.
# Returns NULL when valid, an error string otherwise.
validate_gene_select <- function(p) {
  gs <- p$gene_select %||% "fraction"
  if (!(is.character(gs) && length(gs) == 1L && gs %in% GENE_SELECT_VALUES)) {
    return(paste0("gene_select must be one of: ", paste(GENE_SELECT_VALUES, collapse = ", ")))
  }
  # "custom" needs a gene_list, which this bridge takes no parameter for.
  if (identical(gs, "custom")) {
    return("gene_select 'custom' is not supported: this service accepts no gene_list parameter")
  }
  if (identical(gs, "fraction")) {
    f <- p$fraction %||% 0.05
    if (!(is.numeric(f) && length(f) == 1L && is.finite(f) && f > 0 && f <= 1)) {
      return("fraction must be a number in (0, 1]")
    }
  }
  NULL
}

# ConstructNetwork(soft_power = NULL) silently resolves to Inf when no power clears the
# scale-free fit (min(numeric(0))). Refuse a non-finite power at the boundary instead.
validate_soft_power <- function(p) {
  sp <- p$soft_power
  if (!(is.numeric(sp) && length(sp) == 1L && is.finite(sp) && sp >= 1)) {
    return("soft_power must be a finite number >= 1")
  }
  NULL
}

# ---------------------------------------------------------------------------
# Pipeline helpers (called inside futures — no reference to jobs_env here)
# ---------------------------------------------------------------------------

# enableWGCNAThreads() spawns a worker cluster inside the calling process. Each job already
# runs in a future multisession worker (entrypoint.R: workers = 2), so this nests a second
# parallel layer: up to 2 x WGCNA_THREADS extra R processes, each able to receive a copy of
# the exported data. allowWGCNAThreads() sets the thread count WGCNA reads without spawning a
# cluster — safer, marginally slower — so it is the fallback. WGCNA_THREADS=1 disables both.
wgcna_threads <- function() {
  n <- suppressWarnings(as.integer(Sys.getenv("WGCNA_THREADS", "4")))
  if (is.na(n) || n < 1L) n <- 1L
  if (n == 1L) {
    try(WGCNA::disableWGCNAThreads(), silent = TRUE)
    return(invisible(1L))
  }
  tryCatch(
    WGCNA::enableWGCNAThreads(nThreads = n),
    error = function(e) try(WGCNA::allowWGCNAThreads(nThreads = n), silent = TRUE)
  )
  invisible(n)
}

# enrichR is attached to match the reference analysis environment; nothing here calls it.
# Load order matters: WGCNA masks stats::cor/dist, and igraph/tidyverse mask base and dplyr
# verbs, so the later attach wins for anything they collide on.
pipeline_libs <- function() {
  library(Seurat)
  library(SeuratDisk)
  library(WGCNA)
  library(hdWGCNA)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(igraph)
  library(enrichR)
  patch_hdwgcna_select_network_genes()
  wgcna_threads()
}

# Check metadata columns up front. Without this a missing column surfaces much later as
# whatever the first NULL-indexing consumer happens to throw — table(md[[NULL]]) reports
# "all arguments must have the same length", which says nothing about the real problem.
require_columns <- function(obj, cols) {
  have <- colnames(obj[[]])
  missing <- setdiff(unique(cols), have)
  if (length(missing)) {
    stop("column(s) not found in object metadata: ", paste(missing, collapse = ", "),
         ". Available: ", paste(have, collapse = ", "))
  }
}

# subset() with a column name held in a variable fights Seurat's NSE, so select
# cells by name instead.
subset_condition <- function(obj, col, value) {
  require_columns(obj, col)
  cells <- colnames(obj)[as.character(obj[[col]][, 1]) == value]
  if (!length(cells)) {
    stop("no cells with ", col, " == '", value, "'; values present: ",
         paste(unique(as.character(obj[[col]][, 1])), collapse = ", "))
  }
  subset(obj, cells = cells)
}

# gene_select/fraction default to the hdWGCNA tutorial's "fraction"/0.05. The default matters:
# SelectNetworkGenes otherwise falls back to gene_select = "variable", i.e. VariableFeatures() —
# genes chosen for variance ACROSS all cell types, not genes expressed within the one being
# analysed. Use POST /gene-selection to see how many genes each fraction actually keeps.
#
# Neither argument is a formal of SetupForWGCNA; both reach SelectNetworkGenes through `...`.
# If a future hdWGCNA release stops forwarding, they would be silently ignored and the gene set
# would silently change — hence the length check below, which also turns SelectNetworkGenes'
# stop()-at-zero-genes into an error naming the parameter that caused it.
#
# MetacellsByGroups is kNN-based and takes no seed argument. Seed it here so /analyze rebuilds
# the same metacells whose soft-power curve was inspected in /test-soft-powers — otherwise the
# chosen soft power is applied to a different metacell set.
setup_and_metacells <- function(obj, p) {
  obj <- SetupForWGCNA(
    obj,
    gene_select = p$gene_select,
    fraction    = p$fraction,
    wgcna_name  = p$wgcna_name
  )

  n_genes <- length(GetWGCNAGenes(obj, wgcna_name = p$wgcna_name))
  if (n_genes == 0L) {
    stop("gene_select='", p$gene_select, "' (fraction=", p$fraction,
         ") selected 0 genes; try a smaller fraction")
  }

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

  SetDatExpr(
    obj,
    group_name = p$group_name,
    group.by   = p$cell_type_col,
    assay      = "RNA",
    layer      = "data",
    wgcna_name = p$wgcna_name
  )
}

# Pick the soft power from the numeric power table rather than by eye.
#
# ConstructNetwork(soft_power = NULL) does, verbatim:
#   GetPowerTable() %>% subset(SFT.R.sq >= 0.8 & Power > min_power) %>% .$Power %>% min
# When nothing qualifies that is min(numeric(0)) == Inf, which it passes into WGCNA without
# complaint. recommended_power reproduces the qualifying rule exactly — so the number we
# recommend is the network ConstructNetwork would actually build — but returns NA, never Inf.
#
# The Power > min_power(=3) floor means hdWGCNA can never auto-select powers 1-3, so
# smallest_power (no floor) is reported alongside it and can legitimately be lower.
recommend_power <- function(pt) {
  ok <- !is.na(pt$SFT.R.sq) & pt$SFT.R.sq >= SFT_THRESHOLD

  pick <- function(idx) if (any(idx)) as.integer(min(pt$Power[idx])) else NA_integer_
  recommended <- pick(ok & pt$Power > SFT_MIN_POWER)
  smallest    <- pick(ok)

  max_r2 <- if (all(is.na(pt$SFT.R.sq))) NA_real_ else max(pt$SFT.R.sq, na.rm = TRUE)

  warning_msg <- NA_character_
  if (is.na(recommended)) {
    warning_msg <- paste0(
      "no power in the tested grid reaches SFT.R.sq >= ", SFT_THRESHOLD,
      " above min_power=", SFT_MIN_POWER, " (max observed ", round(max_r2, 3),
      "); a network built at any of these powers is not scale-free"
    )
  }

  list(
    recommended_power = recommended,
    smallest_power    = smallest,
    differ            = !is.na(recommended) && !is.na(smallest) && recommended != smallest,
    sft_threshold     = SFT_THRESHOLD,
    min_power         = SFT_MIN_POWER,
    max_sft_r_sq      = max_r2,
    warning           = warning_msg
  )
}

test_soft_powers <- function(obj, p) {
  obj <- TestSoftPowers(obj,
    networkType = p$network_type,
    wgcna_name  = p$wgcna_name
  )

  pt <- GetPowerTable(obj, wgcna_name = p$wgcna_name)
  write.csv(pt, file.path(p$out_dir, "soft_power_table.csv"), row.names = FALSE)

  rec <- recommend_power(pt)
  # na = "null" so an absent recommendation serializes to JSON null. A bare R NULL would
  # encode as {} and the frontend could not distinguish it from a missing field.
  writeLines(
    jsonlite::toJSON(rec, auto_unbox = TRUE, na = "null", null = "null"),
    file.path(p$out_dir, "soft_power_recommendation.json")
  )

  # PlotSoftPowers marks selected_power on the curve; NA would draw a line at NA, so pass NULL.
  sel <- if (is.na(rec$recommended_power)) NULL else rec$recommended_power
  gg  <- wrap_plots(
    PlotSoftPowers(obj, selected_power = sel, wgcna_name = p$wgcna_name),
    ncol = 2
  )

  save_png(gg, file.path(p$out_dir, "soft_power_plot.png"))
  # PDF is a convenience export. ggsave on a patchwork object goes through grid.draw.patchwork,
  # an untested path against ggplot2 4.x — a failure here must not lose the PNG and the table.
  try(ggsave(file.path(p$out_dir, "soft_power_plot.pdf"), gg, width = 10, height = 5),
      silent = TRUE)

  obj
}

# tom_outdir defaults to the relative path "TOM", which lands in the r-service working
# directory rather than the job's out_dir — so concurrent jobs (workers = 2) sharing a
# wgcna_name collide.
#
# overwrite_tom defaults to FALSE, which reuses any TOM already sitting in out_dir. Now that
# gene_select/fraction are request parameters, that stale TOM is not merely old — it was built
# on a DIFFERENT GENE SET, so reusing it silently produces modules for genes the caller did not
# ask for. Always rebuild.
build_network <- function(obj, p) {
  obj <- ConstructNetwork(
    obj,
    soft_power    = p$soft_power,
    networkType   = p$network_type,
    tom_outdir    = p$out_dir,
    tom_name      = p$wgcna_name,
    overwrite_tom = TRUE,
    wgcna_name    = p$wgcna_name
  )

  # Harmonized ModuleEigengenes (group.by.vars set) requires a populated scale.data layer, but
  # hdWGCNA only guards it by checking @commands for a ScaleData *record* — which an input object
  # can carry while the layer itself is empty. The DA object drops scale.data to save memory
  # (CHANGELOG 0.3.0) yet keeps the command, so the guard passes and the harmonization then
  # multiplies against an empty matrix → "non-conformable arguments", failing every /analyze.
  # Re-run ScaleData on the network genes so the layer actually exists — this is exactly what
  # hdWGCNA's own error message ("Need to run ScaleData before ... group.by.vars") asks for.
  obj <- ScaleData(obj, features = GetWGCNAGenes(obj, p$wgcna_name), verbose = FALSE)

  obj <- ModuleEigengenes(
    obj,
    group.by.vars = p$group_by[[2]],
    wgcna_name    = p$wgcna_name
  )

  ModuleConnectivity(
    obj,
    group.by   = p$cell_type_col,
    group_name = p$group_name,
    wgcna_name = p$wgcna_name
  )
}

# PlotKMEs returns a ggplot of per-gene kME values coloured by module.
# OPEN: if you prefer a connectivity graph, replace with ModuleNetworkPlot()
# which writes per-module PNGs to a directory but has no single combined return.
write_module_outputs <- function(obj, p) {
  modules <- GetModules(obj, wgcna_name = p$wgcna_name)
  write.csv(modules, file.path(p$out_dir, "modules.csv"), row.names = FALSE)
  save_png(PlotKMEs(obj, ncol = 5), file.path(p$out_dir, "network_plot.png"))
}

# Cells per donor per condition, and whether each clears MetacellsByGroups(min_cells = 100).
# Written next to the result because it is the n behind it: PD brains have lost the very
# dopaminergic cells being counted, so far fewer PD donors clear the threshold than control.
write_donor_counts <- function(obj, p) {
  md          <- obj[[]]
  sample_col  <- p$group_by[[2]]
  counts      <- as.data.frame(table(md[[sample_col]], md[[p$condition_col]]),
                               stringsAsFactors = FALSE)
  colnames(counts) <- c("sample_id", "condition", "n_cells")
  counts <- counts[counts$n_cells > 0, ]
  counts$clears_min_cells <- counts$n_cells >= 100
  counts <- counts[order(counts$condition, -counts$n_cells), ]
  write.csv(counts, file.path(p$out_dir, "donor_counts.csv"), row.names = FALSE)
}

# ---------------------------------------------------------------------------
# Job entrypoints
# ---------------------------------------------------------------------------

# hdWGCNA ships no tuner for gene_select/fraction — there is no TestSoftPowers analogue for it.
# What it does offer is a cheap sweep: gene selection runs on the raw object BEFORE metacells,
# so every candidate fraction can be evaluated from a single load by calling the real
# SelectNetworkGenes (via SetupForWGCNA) once per fraction and counting what survives.
#
# Deliberately NOT optimised into one pass that computes each gene's expressed-cell fraction and
# derives all the counts arithmetically: that would reimplement SelectNetworkGenes, and this
# service never reimplements hdWGCNA. N passes over the counts matrix is the price.
run_gene_selection <- function(p) {
  pipeline_libs()
  dir.create(p$out_dir, recursive = TRUE, showWarnings = FALSE)

  obj <- load_seurat(p$h5seurat_path)   # loaded ONCE for the whole sweep
  if (!is.null(p$condition_col)) {
    obj <- subset_condition(obj, p$condition_col, p$ref_group)
  }
  require_columns(obj, c(p$cell_type_col, unlist(p$group_by)))

  n_all <- nrow(obj)

  # SelectNetworkGenes stop()s at 0 genes and warns at <= 100, so a fraction set too high is an
  # error, not an empty row. Isolate each candidate: one unusable fraction must not lose the
  # whole sweep, it must just report why it failed.
  rows <- lapply(p$fractions, function(f) {
    warn <- NA_character_
    err  <- NA_character_
    n <- withCallingHandlers(
      tryCatch({
        tmp <- SetupForWGCNA(obj, wgcna_name = "sweep",
                             gene_select = "fraction", fraction = f)
        length(GetWGCNAGenes(tmp, wgcna_name = "sweep"))
      }, error = function(e) {
        err <<- conditionMessage(e)
        NA_integer_
      }),
      warning = function(w) {
        warn <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    )
    data.frame(
      fraction         = f,
      n_genes          = n,
      n_all_genes      = n_all,
      pct_of_all_genes = round(100 * n / n_all, 2),
      warning          = warn,
      error            = err,
      stringsAsFactors = FALSE
    )
  })

  df <- do.call(rbind, rows)
  write.csv(df, file.path(p$out_dir, "gene_selection.csv"), row.names = FALSE)

  ok <- df[!is.na(df$n_genes), ]
  if (nrow(ok)) {
    gg <- ggplot(ok, aes(x = fraction, y = n_genes)) +
      geom_line() +
      geom_point(size = 2) +
      # SelectNetworkGenes warns below this; a network on <= 100 genes is not worth building.
      geom_hline(yintercept = 100, linetype = "dashed") +
      labs(x = "fraction (min. share of cells expressing the gene)",
           y = "genes selected") +
      theme_cowplot()
    save_png(gg, file.path(p$out_dir, "gene_selection_plot.png"))
  }
  invisible(NULL)
}

run_setup_through_soft_powers <- function(p) {
  pipeline_libs()
  dir.create(p$out_dir, recursive = TRUE, showWarnings = FALSE)

  obj <- load_seurat(p$h5seurat_path)   # .rds or .h5Seurat

  # When comparing conditions, the soft power must be chosen from the curve of the
  # network /module-preservation will actually build — the reference condition alone,
  # not the pooled object.
  if (!is.null(p$condition_col)) {
    obj <- subset_condition(obj, p$condition_col, p$ref_group)
  }

  obj <- setup_and_metacells(obj, p)
  test_soft_powers(obj, p)
  invisible(NULL)
}

run_full_pipeline <- function(p) {
  pipeline_libs()
  dir.create(p$out_dir, recursive = TRUE, showWarnings = FALSE)

  obj <- load_seurat(p$h5seurat_path)   # .rds or .h5Seurat
  obj <- setup_and_metacells(obj, p)
  obj <- test_soft_powers(obj, p)
  obj <- build_network(obj, p)
  write_module_outputs(obj, p)
  invisible(NULL)
}

# Reference/query module preservation, following the call order in the hdWGCNA tutorial
# (https://smorabit.github.io/hdWGCNA/articles/module_preservation.html): the reference
# network is built first, ProjectModules runs BEFORE the query's metacells, and SetDatExpr
# runs on the reference before the query.
run_module_preservation <- function(p) {
  pipeline_libs()
  dir.create(p$out_dir, recursive = TRUE, showWarnings = FALSE)

  obj <- load_seurat(p$h5seurat_path)
  require_columns(obj, c(p$condition_col, p$cell_type_col, unlist(p$group_by)))
  write_donor_counts(obj, p)

  ref <- subset_condition(obj, p$condition_col, p$ref_group)
  rm(obj); gc()

  # Reference: the same network /analyze builds, restricted to one condition.
  ref <- setup_and_metacells(ref, p)
  ref <- test_soft_powers(ref, p)
  ref <- build_network(ref, p)
  write_module_outputs(ref, p)

  # query isn't touched until here, so it is loaded fresh now rather than kept in memory
  # (as a second full Seurat object) throughout ref's setup/network-construction above.
  query <- subset_condition(load_seurat(p$h5seurat_path), p$condition_col, p$query_group)

  # ProjectModules runs SetupForWGCNA and ModuleEigengenes on the query internally, and
  # leaves "projected" as its active hdWGCNA experiment — hence no wgcna_name on the query
  # calls below, matching the tutorial.
  query <- ProjectModules(
    seurat_obj      = query,
    seurat_ref      = ref,
    wgcna_name      = p$wgcna_name,
    wgcna_name_proj = "projected",
    assay           = "RNA"
  )

  set.seed(42)
  query <- MetacellsByGroups(
    seurat_obj  = query,
    group.by    = p$group_by,
    reduction   = "pca",
    k           = p$k,
    max_shared  = p$max_shared,
    ident.group = p$cell_type_col
  )
  query <- NormalizeMetacells(query)

  # ModuleConnectivity leaves its own datExpr on the reference; the tutorial re-sets both
  # here so the two matrices being compared are built the same way.
  ref <- SetDatExpr(
    ref,
    group_name = p$group_name,
    group.by   = p$cell_type_col,
    assay      = "RNA",
    layer      = "data",
    wgcna_name = p$wgcna_name
  )
  query <- SetDatExpr(
    query,
    group_name = p$group_name,
    group.by   = p$cell_type_col,
    assay      = "RNA",
    layer      = "data"
  )

  query <- ModulePreservation(
    query,
    seurat_ref     = ref,
    name           = p$preservation_name,
    verbose        = 3,
    n_permutations = p$n_permutations
  )

  pres <- GetModulePreservation(query, p$preservation_name)
  z    <- pres$Z
  obs  <- pres$obs
  z$module   <- rownames(z)
  obs$module <- rownames(obs)
  pres_df <- merge(z, obs, by = "module", suffixes = c(".Z", ".obs"))

  # Least-preserved module first — that is the end of the table the comparison is about.
  if ("Zsummary.pres" %in% colnames(pres_df)) {
    pres_df <- pres_df[order(pres_df$Zsummary.pres), ]
  }
  write.csv(pres_df, file.path(p$out_dir, "preservation.csv"), row.names = FALSE)

  save_png(
    wrap_plots(
      PlotModulePreservation(query, name = p$preservation_name, statistics = "summary"),
      ncol = 2
    ),
    file.path(p$out_dir, "preservation_plot.png")
  )
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------

#* @post /test-soft-powers
#* @serializer unboxedJSON
function(req, res) {
  p <- req$body
  missing <- setdiff(BASE_FIELDS, names(p))
  if (length(missing)) {
    res$status <- 400L
    return(list(error = paste("missing fields:", paste(missing, collapse = ", "))))
  }
  # condition_col is optional, but is meaningless without the condition to subset to.
  if (!is.null(p$condition_col) && is.null(p$ref_group)) {
    res$status <- 400L
    return(list(error = "ref_group is required when condition_col is given"))
  }
  bad <- validate_gene_select(p)
  if (!is.null(bad)) {
    res$status <- 400L
    return(list(error = bad))
  }

  p <- with_defaults(p)
  job_id <- new_job(p$out_dir)

  future_promise(run_setup_through_soft_powers(p)) %...>%
    (function(...) finish_job(job_id)) %...!%
    (function(err) fail_job(job_id, err))

  list(job_id = job_id)
}

#* @post /gene-selection
#* @serializer unboxedJSON
function(req, res) {
  p <- req$body
  missing <- setdiff(BASE_FIELDS, names(p))
  if (length(missing)) {
    res$status <- 400L
    return(list(error = paste("missing fields:", paste(missing, collapse = ", "))))
  }
  if (!is.null(p$condition_col) && is.null(p$ref_group)) {
    res$status <- 400L
    return(list(error = "ref_group is required when condition_col is given"))
  }

  fractions <- suppressWarnings(as.numeric(p$fractions %||% DEFAULT_FRACTIONS))
  if (!length(fractions) || any(is.na(fractions)) ||
      any(!is.finite(fractions)) || any(fractions <= 0 | fractions > 1)) {
    res$status <- 400L
    return(list(error = "fractions must be numbers in (0, 1]"))
  }

  p <- with_defaults(p)
  p$fractions <- sort(unique(fractions))
  job_id <- new_job(p$out_dir)

  future_promise(run_gene_selection(p)) %...>%
    (function(...) finish_job(job_id)) %...!%
    (function(err) fail_job(job_id, err))

  list(job_id = job_id)
}

#* @post /analyze
#* @serializer unboxedJSON
function(req, res) {
  p <- req$body
  missing <- setdiff(c(BASE_FIELDS, "soft_power"), names(p))
  if (length(missing)) {
    res$status <- 400L
    return(list(error = paste("missing fields:", paste(missing, collapse = ", "))))
  }
  bad <- validate_gene_select(p) %||% validate_soft_power(p)
  if (!is.null(bad)) {
    res$status <- 400L
    return(list(error = bad))
  }

  p <- with_defaults(p)
  job_id <- new_job(p$out_dir)

  future_promise(run_full_pipeline(p)) %...>%
    (function(...) finish_job(job_id)) %...!%
    (function(err) fail_job(job_id, err))

  list(job_id = job_id)
}

#* @post /module-preservation
#* @serializer unboxedJSON
function(req, res) {
  p <- req$body
  required <- c(BASE_FIELDS, "soft_power", "condition_col", "ref_group", "query_group")
  missing <- setdiff(required, names(p))
  if (length(missing)) {
    res$status <- 400L
    return(list(error = paste("missing fields:", paste(missing, collapse = ", "))))
  }
  if (identical(p$ref_group, p$query_group)) {
    res$status <- 400L
    return(list(error = "ref_group and query_group must differ"))
  }
  bad <- validate_gene_select(p) %||% validate_soft_power(p)
  if (!is.null(bad)) {
    res$status <- 400L
    return(list(error = bad))
  }

  p <- with_defaults(p)
  # 250 is the value used in the hdWGCNA module preservation tutorial; the
  # ModulePreservation() default is 500.
  p$n_permutations    <- p$n_permutations %||% 250L
  p$preservation_name <- p$preservation_name %||% paste0(p$ref_group, "-vs-", p$query_group)

  job_id <- new_job(p$out_dir)

  future_promise(run_module_preservation(p)) %...>%
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
function(job_id, res) serve_csv(job_id, res, "modules.csv")

# Not unboxedJSON: its toJSON() defaults encode NA_integer_ as the STRING "NA" and NULL as {},
# so a missing recommendation would reach the client as a truthy value instead of null —
# precisely wrong in the one case (no power reaches 0.8) the caller most needs to detect.
#* @get /results/<job_id>/soft-powers
#* @serializer json list(auto_unbox = TRUE, na = "null", null = "null")
function(job_id, res) serve_soft_powers(job_id, res)

#* @get /results/<job_id>/gene-selection
#* @serializer unboxedJSON
function(job_id, res) serve_csv(job_id, res, "gene_selection.csv")

#* @get /results/<job_id>/gene-selection-plot
#* @serializer contentType list(type="image/png")
function(job_id, res) serve_png(job_id, res, "gene_selection_plot.png")

#* @get /results/<job_id>/preservation
#* @serializer unboxedJSON
function(job_id, res) serve_csv(job_id, res, "preservation.csv")

#* @get /results/<job_id>/donor-counts
#* @serializer unboxedJSON
function(job_id, res) serve_csv(job_id, res, "donor_counts.csv")

#* @get /results/<job_id>/plot
#* @serializer contentType list(type="image/png")
function(job_id, res) serve_png(job_id, res, "network_plot.png")

#* @get /results/<job_id>/soft-power-plot
#* @serializer contentType list(type="image/png")
function(job_id, res) serve_png(job_id, res, "soft_power_plot.png")

#* @get /results/<job_id>/preservation-plot
#* @serializer contentType list(type="image/png")
function(job_id, res) serve_png(job_id, res, "preservation_plot.png")
