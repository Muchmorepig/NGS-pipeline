library(zoo)

smooth_byRow <- function(data, method = c("moving_average", "loess"), 
                                  window_size = 3, span = 0.2) {
  method <- match.arg(method)
  
  if (method == "moving_average") {
    smoothed <- t(apply(data, 1, function(x) rollmean(x, k = window_size, fill = NA, align = "center")))
  } else if (method == "loess") {
    smoothed <- t(apply(data, 1, function(x) predict(loess(x ~ seq_along(x), span = span))))
  }
  
  return(smoothed)
}
