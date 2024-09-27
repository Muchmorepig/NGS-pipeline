#' @export

mf_mean <- function(mat, group, samples = colnames(mat)) {
  if (length(samples) != length(group)) {
    stop("Check: group's length = samples ?")
  }

  names(samples) <- group
  nv <- unique(group)
  li <- lapply(nv, FUN = function(x) {
    message("Group: ", x)
    dat <- rowMeans(mat[, names(samples) == x])
    return(setNames(data.frame(dat), x))
  })

  df <- do.call(cbind, li)
  rownames(df) <- rownames(mat)
  return(df)
}

# mf_mean1 <- function(mat, info, name) {
#   suppressMessages(require(magrittr))
#   info <- as.data.frame(info)
#   nv <- info[, which(colnames(info) != name)]
#   names(nv) <- info[, name]
#   n <- names(nv) %>% unique()
#   li <- list()
#   for (i in n) {
#     print(i)
#     li[i] <- mat[, c(which(names(nv) == i))] %>%
#       apply(., 1, mean) %>%
#       as.data.frame()
#   }
#   dat <- as.data.frame(li, row.names = rownames(mat))
#   return(dat)
# }
