#' calc.size.factor
#' Calculate the size factor
#'
#' @param x Combination of berved single cell HiC matrix, with each column being the upper triangular of a single cell HiC matrix.
#'
#' @return A list of size factor, defined as the relative sequencing depth to the median sequencing depth among all single cells.
#' @export
#'
#' @examples
#' calc.size.factor(X_count)

calc.size.factor <- function(x) {

  sf <- Matrix::colSums(x)/median(Matrix::colSums(x))

  scale.sf <- 1

  list(unname(sf), scale.sf)

}

