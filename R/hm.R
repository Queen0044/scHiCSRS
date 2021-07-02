#' hm
#'
#' This function draws heatmap of HiC data so that we can visually compares the imputation results.
#'
#' @param datvec A vector of upper triangular mamtrix.
#' @param n Number of bins.
#'
#' @return Heatmap of the matrix.
#' @export
#'
#' @examples
#' # Heatmap of the first single cell in the simulated data.
#' hm(simudat[,1], 61)
#'
hm <- function(datvec, n) {
  library(plsgenomics)
  normmatrix <- function(matr) {
    maxvalue <- max(matr[upper.tri(matr)])
    minvalue <- min(matr[upper.tri(matr)])
    normmatr <- (matr-minvalue)/(maxvalue-minvalue)
    return(normmatr)
  }

  mat <- matrix(0, n, n)
  mat[upper.tri(mat, diag=FALSE)] <- datvec
  return(matrix.heatmap(normmatrix(mat)))
}

