rToP <- function(rvalue, n = dim(ATAC_CPM_mean)[2]){
  # Define a function to get pvalue from correlation given a sample size
  # 这里n是自由度，等于样本数目
  # 用来计算相关性的显著性
  t <- abs(rvalue)*sqrt(n-2) / sqrt(1-rvalue^2)
  2*pt(t, n-2, lower=FALSE)
}