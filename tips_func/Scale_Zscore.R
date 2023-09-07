#Scale_Zscore_1 <- function(x) {
#  rm <- rowMeans(x)
#  x <- sweep(x, 1, rm)
#  sx <- apply(x, 1, sd)
#  x <- sweep(x, 1, sx, "/")
#  
#}

Scale_Zscore <- function(data){
  scale_data <- apply(data, 1, function(y) {
    (y - mean(y)) / (sd(y)^as.logical(mean(y)))
  })
  
  scale_data_t <- t(scale_data)
  
  return(scale_data_t)
}
