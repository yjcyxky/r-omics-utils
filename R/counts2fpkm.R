# More details on https://cloud.tencent.com/developer/article/1669450

#' Convert counts to FPKM
#'
#' @param counts A data frame which rownames are ensembl id and colnames are sample name.
#' @return a data frame which contains FPKM values.
#' @export
counts2fpkm <- function(counts) {
    n <- sum(counts)
    exp(log(counts) + log(1e9) - log(length) - log(n))
}

#' Convert FPKM to TPM
#'
#' @param fpkmExp A data frame with FPKM values (rows are genes and columns are samples)
#' @return a data frame with TPM values
#' @export
fpkm2tpm <- function(fpkmExp) {
    tpm <- apply(fpkmExp, 2, function(fpkm) {
        exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
    })
    as.data.frame(tpm)
}
