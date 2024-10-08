minmax <- function(x) {
  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)

  if (min_val == max_val) {
    stop("min = max")
  }

  (x - min_val) / (max_val - min_val)
}

mf_minmax <- function(x) {
  t(apply(x, 1, minmax))
}
