#' @export
CompInfo <- function(colData = colData,condition = "group"){
  suppressMessages(
    require(magrittr)
  )
  con <- colData[,condition] %>% unique()
  l <- length(con) - 1
  c <- c()
  for (i in 0:(length(con)-2)) {
    c[seq(l*(l-i)-i , l*(l-i))] <- seq(l*(l-i)-i , l*(l-i))
  }
  c <- na.omit(c) %>% as.numeric()
  comp <- expand.grid(con, con, stringsAsFactors = F) %>%
    subset(Var1 != Var2) %>%
    mutate(key = paste0(Var1,"_VS_",Var2)) %>% 
    .[c,]
  rownames(comp) <- NULL
  return(comp)
}