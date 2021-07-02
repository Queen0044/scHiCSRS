#' log_normalization
#' Preprocess the data by normalizing the data so that each single cell have the same sequencing depth as the median sequencing depth and take a log transformation with pseudocount 1.
#'
#' Take log normalization of the data.
#'
#' @param x Combination of berved single cell HiC matrix, with each column being the upper triangular of a single cell HiC matrix.
#'
#' @return A matrix after normalization and log transformation.
#' @export
#'
#' @examples
#' log_normalization(X_count)
log_normalization = function(x){

    x <- as.matrix(x)

    n <- dim(x)[2]

    gene.exprs.count <- rowSums(x != 0)

    sf <- colSums(x)/median(colSums(x))

    return(log(sweep(x, 2, sf, '/')+1))

}
