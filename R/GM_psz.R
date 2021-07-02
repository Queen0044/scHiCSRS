#' GM_psz
#' Calculate Psz of a imputed matrix through a Gaussian mixture mdoel.
#'
#' @param imputed A imputed matrix.
#' @param observed The observed matrix.
#' @param nmix Number of Gaussian mixtures.
#'
#' @return
#' @export
#'
#' @examples
#' GM_Psz(simudat_res)
GM_Psz <- function(observed, imputed, nmix){
   prob = matrix(0, nrow = nrow(imputed), ncol = ncol(imputed))
   prob[which(apply(observed,1,sum)==0),]=1
   dat1=observed
   dat1[which(apply(observed,1,sum)==0),]=1
   index=(dat1==0)
   dat2=imputed[index]
   paramt=get_mix(dat2, nmix)
   pp=dmix(dat2,paramt)
   prob[index]=pp
   return(prob)
}
