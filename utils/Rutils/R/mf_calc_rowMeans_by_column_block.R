#' @export
calc_RowMeansByColumnBlocks <- function(data, block_size = 300) {
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("Input data must be a data frame or matrix.")
  }
  
  num_cols <- ncol(data)
  num_rows <- nrow(data)
  
  new_num_cols <- ceiling(num_cols / block_size)
  
  res <- matrix(NA, nrow = num_rows, ncol = new_num_cols)
  
  for (i in 1:new_num_cols) {
    start_col <- (i - 1) * block_size + 1
    end_col <- min(i * block_size, num_cols)
    
    res[, i] <- rowMeans(data[, start_col:end_col, drop = FALSE], na.rm = TRUE)
  }
  
  res <- as.data.frame(res)
  rownames(res) <- rownames(data)
  return(res)
}

