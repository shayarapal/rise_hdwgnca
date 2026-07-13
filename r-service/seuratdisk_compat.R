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
