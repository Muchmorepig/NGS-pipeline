#' Smooth Data by Row
#'
#' This function smooths the data by row using either moving average or loess smoothing.
#' @param data A matrix or data frame where smoothing is to be applied.
#' @param method Method of smoothing: either "moving_average" or "loess".
#' @param window_size Size of the moving average window (only applicable if method is "moving_average").
#' @param span Span parameter for loess (only applicable if method is "loess").
#' @return A matrix of smoothed data.
#' @import zoo
#' @export
smooth_byRow <- function(
    data,
    method = c("moving_average", "loess"),
    window_size = 3,
    span = 0.2) {
  # Match the method argument to ensure it is valid
  method <- match.arg(method)

  if (method == "moving_average") {
    # Apply the moving average smoothing
    smoothed <- t(apply(data, 1, function(x) rollmean(x, k = window_size, fill = NA, align = "center")))
  } else if (method == "loess") {
    # Apply LOESS smoothing
    smoothed <- t(apply(data, 1, function(x) predict(loess(x ~ seq_along(x), span = span))))
  }

  # Return the smoothed data
  return(smoothed)
}
