#' SRS_ROC
#'
#'This package draws ROC (Receiver operating characteristic) curve to visually demonstrate ability
#'to tell SZ from DO.
#'
#' @param single Observed sngle cell with each column being the upper triangular of single cell.
#' @param truecount Underline true count of simulated data.
#' @param result Result from SRS funciton.
#'
#' @return A plot of ROC curve.
#' @export
#'
#' @examples
#' SRS_ROC(simudat,simudat_true,simudat_res)
SRS_ROC <- function(single, truecount, result){

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
    IMP <-result$IMP1
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

  ptsz_range <- seq(0,1,length.out = 100)
  output=NULL
  for (k in 1:length(ptsz_range)) {
    tt <- rev(threshold)[which.min(abs(rev(PTSZout)-ptsz_range[k]))]

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
    summa_mean=data.frame(PTSZ=mean(PTSZ),PTDO=mean(PTDO),thresh=tt)
    output <- rbind(output,summa_mean)
  }
  #return(output)
  plot(1-output$PTSZ,output$PTDO, type = "b", xlab = "1-PTSZ", ylab = "PTDO",
       pch=3,col="darkorchid", lwd=2, ylim = c(0,1))
}

