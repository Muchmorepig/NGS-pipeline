counts2TPM <- function(count = count, efflength = efflen) {
    RPK <- count / (efflength / 1000) # reads per kilobase --长度标准化
    PMSC_rpk <- sum(RPK) / 1e6 # RPK的 per million scaling factor --深度标准化
    RPK / PMSC_rpk
}