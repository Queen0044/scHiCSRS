#' dmix
#' Estimate probability of SZ.
#'
#' @param x A vector of data.
#' @param pars The estimated Gaussian mixture parameters.
#'
#' @return
#' @export
#'
#' @examples
#' dat=rnorm(100)
#' paramt=get_mix(dat, nmix=3)
#' dmix(dat,paramt)

dmix <- function(x,pars) {
  output=c()
  G=length(pars)/3
  d=rep(0,(G-1))
  for (gg in 1:(G-1)) {
    d[gg]=pars[G+gg+1]-pars[G+gg]
  }
  if((10*d[1])<=d[2]){
    ind=c(1,2)
  }else{
    if((10*d[2])<=d[3]){
      ind=c(1,2,3)
    }else{
      ind=c(1)
    }
  }
  for (j in 1:length(x)) {
    comp=c()
    for (i in 1:G) {
      comp[i]=pars[i]*dnorm(x[j], mean = pars[i+G], sd=sqrt(pars[i+2*G]))
    }
    output[j]=sum(comp[ind])/sum(comp)
  }
  return(output)
}
