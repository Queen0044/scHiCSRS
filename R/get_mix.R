###

#' get_mix
#' Estimate parameters in the Gaussian mixture distribution.
#'
#' @param xdata A vector of data.
#' @param nmix Number of Gaussian mixtures.
#'
#' @return
#' @export
#'
#' @import mclust
#' @examples
#' get_mix(rnorm(100), nmix=3)
#'
get_mix <- function(xdata, nmix) {
  fit = Mclust(xdata, G=nmix, model="V")
  paramt = c(fit$parameters$pro, fit$parameters$mean,fit$parameters$variance$sigmasq)
  names(paramt)=NULL
  return(paramt)
}

