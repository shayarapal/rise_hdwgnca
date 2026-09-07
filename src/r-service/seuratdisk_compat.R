# Input loading for the bridge: .rds (Seurat v5 native) and .h5Seurat (legacy).
#
# Why this file exists: SeuratDisk has been unmaintained since 2022 and does not work
# against the installed SeuratObject 5.x / R 4.5. It breaks in three places, and it
# cannot serialise a v5 Assay5 at all -- there is no maintained h5Seurat writer for
# Seurat v5. .rds is the supported format (see install.R). h5Seurat is kept working
# for files written by older Seurat, hence the shims.

# 1+2. GetAssayData()/SetAssayData(slot=) became *defunct* in SeuratObject 5.0.0
#      (renamed layer=); SeuratDisk still calls them with slot=.
# 3.   GetImages() calls unlist(x = assays.images, index$global$images) -- the second
#      argument binds to unlist's `recursive` parameter instead of being concatenated.
#      R >= 4.5 rejects a zero-length `recursive`, so LoadH5Seurat() errors on every
#      object with no spatial images. Restores the intended c(unlist(...), ...).
# Must run inside the process that calls LoadH5Seurat (i.e. the future worker).
patch_seuratdisk <- function() {
  put <- function(env, nm, val) {
    if (!nm %in% ls(env, all.names = TRUE)) return(invisible(FALSE))
    unlockBinding(nm, env); assign(nm, val, envir = env); lockBinding(nm, env)
  }
  imp <- parent.env(asNamespace("SeuratDisk"))
  put(imp, "GetAssayData", function(object, slot, layer, ...) {
    if (!missing(slot) && missing(layer)) layer <- slot
    if (missing(layer)) layer <- "data"
    SeuratObject::GetAssayData(object = object, layer = layer, ...)
  })
  put(imp, "SetAssayData", function(object, slot, new.data, layer, ...) {
    if (!missing(slot) && missing(layer)) layer <- slot
    if (missing(layer)) layer <- "data"
    SeuratObject::SetAssayData(object = object, layer = layer, new.data = new.data, ...)
  })
  put(asNamespace("SeuratDisk"), "GetImages", function(images, index, assays = NULL) {
    if (isFALSE(images)) return(NULL)
    if (!is.null(images) && all(is.na(images))) return(index$global$images)
    if (is.null(images)) images <- unique(unlist(lapply(names(index), function(x) index[[x]]$images)))
    assays <- SeuratDisk:::GetAssays(assays = assays, index = index)
    assays.images <- lapply(names(assays), function(x) index[[x]]$images)
    intersect(images, unique(c(unlist(assays.images), index$global$images)))
  })
  invisible(TRUE)
}

# hdWGCNA::SelectNetworkGenes (gene_select="fraction" path) binarizes the counts matrix
# with `cur[cur > 0] <- 1` on a dgCMatrix chunk. Boolean-matrix-indexed assignment on a
# sparse Matrix is a well-documented Matrix-package memory trap: it does not reliably stay
# sparse internally, and on a chunk sized by ncol(expr_mat) (cells) rather than genes this
# balloons to tens of GB for cell counts in the tens of thousands, killing the process
# before SetupForWGCNA's caller ever gets a chance to run. `cur@x[cur@x > 0] <- 1` sets the
# same stored (nonzero) values to 1 directly, producing an identical binarized matrix
# without ever leaving sparse representation. Must run before hdWGCNA's namespace is used.
patch_hdwgcna_select_network_genes <- function() {
  fixed <- function(seurat_obj, gene_select = "variable", fraction = 0.05,
                     group.by = NULL, gene_list = NULL, assay = NULL, wgcna_name = NULL) {
    if (is.null(wgcna_name)) wgcna_name <- seurat_obj@misc$active_wgcna
    if (!(gene_select %in% c("variable", "fraction", "all", "custom"))) {
      stop(paste0("Invalid selection gene_select: ", gene_select,
                  ". Valid gene_selects are variable, fraction, all, or custom."))
    }
    if (is.null(assay)) assay <- Seurat::DefaultAssay(seurat_obj)
    expr_mat <- if (hdWGCNA:::CheckSeurat5()) {
      SeuratObject::LayerData(seurat_obj, layer = "counts", assay = assay)
    } else {
      Seurat::GetAssayData(seurat_obj, slot = "counts", assay = assay)
    }

    if (gene_select == "fraction") {
      n_chunks <- ceiling(ncol(expr_mat) / 10000)
      chunks <- if (n_chunks == 1) factor(rep(1), levels = 1) else cut(1:nrow(expr_mat), n_chunks)

      expr_mat <- do.call(rbind, lapply(levels(chunks), function(x) {
        cur <- expr_mat[chunks == x, ]
        cur@x[cur@x > 0] <- 1   # fixed: stays sparse, same result as cur[cur > 0] <- 1
        cur
      }))

      if (!is.null(group.by)) {
        groups <- unique(seurat_obj@meta.data[[group.by]])
        group_gene_list <- list()
        for (cur_group in groups) {
          cur_expr <- expr_mat[, seurat_obj@meta.data[[group.by]] == cur_group]
          gene_filter <- Matrix::rowSums(cur_expr) >= round(fraction * ncol(cur_expr))
          group_gene_list[[cur_group]] <- rownames(seurat_obj)[gene_filter]
        }
        gene_list <- unique(unlist(group_gene_list))
      } else {
        gene_filter <- Matrix::rowSums(expr_mat) >= round(fraction * ncol(seurat_obj))
        gene_list <- rownames(seurat_obj)[gene_filter]
      }
    } else if (gene_select == "variable") {
      gene_list <- Seurat::VariableFeatures(seurat_obj)
    } else if (gene_select == "all") {
      gene_list <- rownames(seurat_obj)
    } else if (gene_select == "custom") {
      gene_list <- unique(gene_list)
      if (!all(gene_list %in% rownames(seurat_obj))) {
        stop("Some selected features are not found in rownames(seurat_obj).")
      }
      if (!is.null(gene_list) && !is.character(gene_list)) {
        stop("Invalid type for gene_list, must be a character vector.")
      }
    }

    if (length(gene_list) == 0) stop("No genes found")
    if (length(gene_list) <= 100) {
      warning(paste0("Very few genes selected (", length(gene_list),
                     "), perhaps use a different method to select genes."))
    }
    hdWGCNA:::SetWGCNAGenes(seurat_obj, gene_list, wgcna_name)
  }
  assignInNamespace("SelectNetworkGenes", fixed, ns = "hdWGCNA")
  invisible(TRUE)
}

# Dispatch on extension. Without the guards below, handing this a non-HDF5 file (an .R
# script, an .rds) drops straight into the HDF5 C library, which reports "unable to read
# superblock" across ~20 lines of H5F.c stack instead of saying "that isn't an h5Seurat".
load_seurat <- function(path) {
  if (!nzchar(path %||% "")) stop("no input path given")
  if (!file.exists(path)) stop("file not found: ", path)

  ext <- tolower(tools::file_ext(path))
  obj <- switch(ext,
    rds = readRDS(path),
    h5seurat = { patch_seuratdisk(); SeuratDisk::LoadH5Seurat(path, verbose = FALSE) },
    stop("unsupported input format '.", ext, "': expected a .rds or .h5Seurat file, got ", path)
  )

  if (!inherits(obj, "Seurat")) {
    stop("file is not a Seurat object (got ", class(obj)[1], "): ", path)
  }
  obj
}
