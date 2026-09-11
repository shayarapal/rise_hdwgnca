load("/shared/DATA_GSE233866/results/snc_intact_network/cacna1d_mouse_snc_intact_TOM.rda")
cat("Objects loaded:\n")
print(ls())
for (nm in ls()) {
  obj <- get(nm)
  cat(sprintf("\n--- %s ---\n", nm))
  cat("class:", class(obj), "\n")
  if (is.matrix(obj) || inherits(obj, "dist")) {
    cat("dim:", paste(dim(as.matrix(obj)), collapse=" x "), "\n")
    rn <- rownames(as.matrix(obj))
    cat("has rownames:", !is.null(rn), "\n")
    if (!is.null(rn)) cat("first 5 rownames:", paste(head(rn,5), collapse=", "), "\n")
  }
}
