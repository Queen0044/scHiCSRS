#' lnorm
#' Normalize and log transform data with 1.01 pseudocount. Used as the input of Gaussian mixture model.
#'
#' @param raw_count
#'
#' @return
#' @export
#'
#' @examples
lnorm <- function(raw_count) {
  totalCounts_by_cell = colSums(raw_count)
  totalCounts_by_cell[totalCounts_by_cell == 0] = 1
  raw_count = sweep(raw_count, MARGIN = 2, median(totalCounts_by_cell)/totalCounts_by_cell, FUN = "*")
  if (min(raw_count) < 0) {
    stop("smallest read count cannot be negative!")
  }
  count_lnorm = log10(raw_count + 1.01)
  return(count_lnorm)
}
