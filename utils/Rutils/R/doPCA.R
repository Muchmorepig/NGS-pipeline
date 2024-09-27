#' @export
doPCA <- function(indata,
                  coldata,
                  intgroup = "groups",
                  ntop = ifelse(0.1 * nrow(indata) > 2000, 0.1 * nrow(indata), 2000),
                  pcsToUse = 1:3) {
  indata <- as.matrix(indata)
  rv <- apply(indata, 1, var)

  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(indata[select, ]), rank. = length(pcsToUse))

  percentVar <- pca$sdev^2 / sum(pca$sdev^2)

  if (!all(intgroup %in% names(coldata))) {
    stop("the specified 'intgroup' does not match any column in 'coldata'")
  }

  group <- if (is.factor(coldata[[intgroup]]) || is.character(coldata[[intgroup]])) {
    factor(coldata[[intgroup]])
  } else {
    stop("'intgroup' must be a factor or character vector")
  }

  d <- setNames(as.data.frame(pca$x[, pcsToUse, drop = FALSE]), paste0("PC", pcsToUse))
  d$group <- group
  d$samples <- colnames(indata)
  attr(d, "percentVar") <- percentVar[pcsToUse]
  return(d)
}
