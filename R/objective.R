#' objective
#' Define the objective function.
#'
#' @param X Combination of berved single cell HiC matrix, with each column being the upper triangular of a single cell HiC matrix.
#' @param lambda The tunning parameter of objective function.
#' @param A Coefficient matrix A
#' @param B Coefficient matrix B
#'
#' @return
#' @export
#'
#' @examples
#' objective(X, lambda1, A, B)
objective <- function(X, lambda, A, B){

  objval <- (0.5*norm(X - A%*%X - X%*%B, 'F')^2 + lambda*sum(abs(A)) + lambda*sum(abs(B)))

  return(objval)
}
