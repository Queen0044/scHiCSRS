#' regularizer_define
#' Define regularizer
#'
#' @param weight_matrix Weighte matrix
#' @param lambda1 Lambda1
#' @param lambda2 Lambda2
#'
#' @return
#' @export
#'
#' @examples
#' regularizer_define(weight_matrix = weight_matrix ,lambda1 = lambda1,  lambda2 = lambda2)
regularizer_define <- function(weight_matrix ,lambda1 = 1.0, lambda2 = 1e10){
  lambda1 * k_sum(k_abs(weight_matrix), axis = c(1,2)) + lambda2 * tf$linalg$trace(tf$square(weight_matrix))
}

