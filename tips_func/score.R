min_max_scale <- function(x) {
  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)
  
  if (min_val == max_val) {
    stop("min = max")
  }
  
  (x - min_val) / (max_val - min_val)
}

"%!in%" <- function(x, y) !("%in%"(x, y))

sce.calc.score <- function(sce, use_genes) {
  stopifnot(all(use_genes %in% rownames(sce)))
  
  x <- assay(sce, "counts")
  # total UMI per cell
  sf <- colSums(x)
  
  # dividing UMIs of each gene by total UMIs per cell and scaling/normalizing
  score <- t(x[use_genes, , drop = F]) / sf * 1e6
  score <- scale(score)
  
  score <- score[, rowSums(apply(score, 1, is.nan)) == 0]
  
  # taking the mean in each cell
  score <- rowMeans(score)
  
  q01 <- quantile(score, .01)
  q99 <- quantile(score, .99)
  score[score < q01] <- q01
  score[score > q99] <- q99
  score <- min_max_scale(score)
  
  return(score)
}

# Function for calculating the tau score for each cell type. x= vector of numbers for which we wanna calculate the tau score.
tau <- function(x) {
  if (!is.numeric(x)) {
    stop("Input must be a numeric vector")
  }
  if (length(x) <= 1) {
    warning("Input vector length must be greater than 1")
    return(NA)
  }
  
  max_x <- max(x, na.rm = TRUE)
  
  if (max_x > 0) {
    relative_expression <- 1 - (x / max_x)
    tau_value <- sum(relative_expression, na.rm = TRUE) / (length(x) - 1)
  } else {
    tau_value <- NA
  }
  
  return(tau_value)
}
