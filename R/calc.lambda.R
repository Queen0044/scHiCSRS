#' calc.lambda
#' Calculate lambda that is used as the tuning parameter in the objective function.
#'
#' @param X_count Combination of berved single cell HiC matrix, with each column being the upper triangular of a single cell HiC matrix.
#'
#' @return A number that is used as the tuning parameter in the objective function.
#' @export
#'
#' @examples
#' calc.lambda(X_count, percent)
calc.lambda <- function(X_count){

  X <- log_normalization(X_count)

  lambda <- sd(X)

}
