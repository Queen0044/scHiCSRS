####################################################
###             Set PTSZ=0.95                     ##
####################################################
#' summa2
#'
#'This function calculates PTDO when fix PTSZ=0.95.
#' @param single Observed single cells matrix with each column being the upper triangular of a single cell.
#' @param truecount Underline true counts from simulation.
#' @param result Result form SRS function.
#'
#' @return A vector of PTDO and its SD when fixing PTSZ to be 0.95, and the threshold used in that case.
#' @export
#'
#' @examples
#' summa2(data$singledat, data$truecount, res)
summa2 <- function(single, truecount, result){

  single_sum <- apply(single, 2, sum)
  max_single <- max(single_sum)
  singlelam <- single_sum/max_single

  ## true 0 positions figured out by MCMC
  PTSZfun <- function(thresh) {
    posi <- which(result$pii>=thresh,TRUE)
    posirows <- NULL
    for (i in 1:(dim(posi)[1])) {
      posirows <- c(posirows,matrow(posi[i,1],posi[i,2]))
    }
    IMP <- result$IMP1
    IMP[posirows,]<-0
    PTSZ=list()
    PTDO=list()

    for(j in 1:ncol(single)){
      indexobserved0=(single[,j]==0)
      predictv=IMP[,j][indexobserved0]
      Truevalue=truecount[,j][indexobserved0]

      PTSZ[[j]]=sum(Truevalue==0 & predictv==0)/sum(Truevalue==0)

      PTDO[[j]]=sum(Truevalue>0 & predictv>0)/sum(Truevalue>0)
    }
    PTSZ=unlist(PTSZ)
    PTDO=unlist(PTDO)
    return(mean(PTSZ))
  }

  threshold <- seq(0.001, max(result$pii), length.out = 200)
  PTSZout <- unlist(lapply(threshold,PTSZfun))
  tt <- rev(threshold)[which.min(abs(rev(PTSZout)-0.95))]

  posi <- which(result$pii>=tt,TRUE)
  posirows <- NULL
  for (i in 1:(dim(posi)[1])) {
    posirows <- c(posirows,matrow(posi[i,1],posi[i,2]))
  }

  IMP <- result$IMP1
  IMP[posirows,]<-0

  PTSZ=list()
  PTDO=list()
  for(j in 1:ncol(single)){
    indexobserved0=(single[,j]==0)
    predictv=IMP[,j][indexobserved0]
    Truevalue=truecount[,j][indexobserved0]

    PTSZ[[j]]=sum(Truevalue==0 & predictv==0)/sum(Truevalue==0)
    PTDO[[j]]=sum(Truevalue>0 & predictv>0)/sum(Truevalue>0)
  }

  PTSZ=unlist(PTSZ)
  PTDO=unlist(PTDO)
  summa_mean=data.frame(PTDO=mean(PTDO), SD2=sd(PTDO), thresh=tt)
  rownames(summa_mean)=NULL
  return(summa_mean)
}
