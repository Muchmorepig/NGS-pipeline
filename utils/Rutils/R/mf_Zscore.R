#' @export

mf_zscore <- function(dat) {
    scale_data <- apply(dat, 1, function(y) {
        if (sd(y) == 0) {
            rep(0, length(y))
        } else {
            (y - mean(y)) / sd(y)
        }
    })
    return(t(scale_data))
}

# mf_zscore <- function(dat) {
#     scale_data <- apply(dat, 1, function(y) {
#         if (sd(y) == 0) {
#             (y - mean(y)) / (sd(y)^0)
#         } else {
#             (y - mean(y)) / sd(y)
#         }
#     })
#     scale_data_t <- t(scale_data)
#     return(scale_data_t)
# }
