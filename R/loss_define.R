#' loss_define
#' Calculate the loss value.
#'
#' @param y_true The true matrix.
#' @param y_pred The predicted matrix.
#'
#' @return The sum of squared distance between the predicted and the true matrix.
#' @export
#'
#' @examples
#' loss_define(y_true, y_pred)
loss_define <- function(y_true, y_pred){
  0.5 * k_sum(k_square(y_true - y_pred), axis = c(1,2))
}
