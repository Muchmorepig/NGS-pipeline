mf_mean1 <- function(mat, info, name) {
  suppressMessages(require(magrittr))
  info <- as.data.frame(info)
  nv <- info[, which(colnames(info) != name)]
  names(nv) <- info[, name]
  n <- names(nv) %>% unique()
  li <- list()
  for (i in n) {
    print(i)
    li[i] <- mat[, c(which(names(nv) == i))] %>%
      apply(., 1, mean) %>%
      as.data.frame()
  }
  dat <- as.data.frame(li, row.names = rownames(mat))
  return(dat)
}

mf_mean2 <- function(mat, samples = colnames(mat), group) {
  suppressMessages(require(magrittr))
  if(length(samples) != length(group))
    stop("Check：group's length = samples ?")
  
  names(samples) <- group
  nv <- unique(group)
  li <- lapply(nv, FUN = function(x){
    cat("Group:", x, "\n")
    dat <- mat[, names(samples) == x] %>%
      apply(., 1, mean) %>%
      as.data.frame()
    
    colnames(dat) <- x
    return(dat)
  })
  return(do.call(cbind, li))
}
