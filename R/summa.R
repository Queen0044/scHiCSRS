#' summa
#'
#' This function summarizes the accuracy of simulated data.
#'
#' @param single Observed single cells matrix, with each column being the upper triangular of a single cell.
#' @param truecount Underline true count of the simulated data.
#' @param imputeduted Imputed matrix from SRS.
#'
#' @return A list of accuracy measurements with mean and standard error.
#' @export
#'
#' @examples
#' summa(data$singledat, data$truecount,res$imputed)
summa <- function(truecount,predicted,observed,prob,cut=0.5)
{
  PTSZ=list()
  PTDO=list()
  MSE=list()
  cor_0=list()
  cor_sampling=list()
  cor_all=list()
  RE_sampling=data.frame()
  RE_all=data.frame()
  AE_0=data.frame()
  AE_sampling=data.frame()
  AE_all=data.frame()

  predicted[prob>cut] = 0

  for(j in 1:ncol(observed))
  {
    indexobserved0=(observed[,j]==0)
    indexnotrue0=(truecount[,j]>0)

    predictv=predicted[,j][indexobserved0]
    Truevalue=truecount[,j][indexobserved0]

    if(sum(Truevalue==0)>0){
      PTSZ[[j]]=sum(Truevalue==0 & predictv==0)/sum(Truevalue==0)
    }else
    {
      PTSZ[[j]]=NA
    }

    PTDO[[j]]=sum(Truevalue>0 & predictv>0)/sum(Truevalue>0)

    MSE[[j]]=mean((predicted[,j]-truecount[,j])^2)

    RE_sampling=rbind(RE_sampling, data.frame(RE=(abs(predictv-Truevalue)[Truevalue>0]/Truevalue[Truevalue>0]),Sample=j))
    RE_all=rbind(RE_all, data.frame(RE=(abs(predicted[,j]-truecount[,j])[indexnotrue0]/truecount[,j][indexnotrue0]),Sample=j))

    AE_0=rbind(AE_0, data.frame(AE=abs(predictv-Truevalue),Sample=j))
    AE_sampling=rbind(AE_sampling, data.frame(AE=abs(predictv[Truevalue>0]-Truevalue[Truevalue>0]),Sample=j))
    AE_all=rbind(AE_all, data.frame(AE=abs(predicted[,j]-truecount[,j]),Sample=j))

    cor_0[[j]]=cor(Truevalue,predictv)
    cor_sampling[[j]]=cor(Truevalue[Truevalue>0],predictv[Truevalue>0])
    cor_all[[j]]=cor(predicted[,j],observed[,j])
  }

  PTSZ=unlist(PTSZ)
  PTDO=unlist(PTDO)
  MSE=unlist(MSE)
  cor_0=unlist(cor_0)
  cor_sampling=unlist(cor_sampling)
  cor_all=unlist(cor_all)

  out1=data.frame(PTSZ=mean(PTSZ),PTDO=mean(PTDO),MSE=mean(MSE),
                  RE_sampling=mean(RE_sampling$RE),RE_all=mean(RE_all$RE),
                  AE_0=mean(AE_0$AE),AE_sampling=mean(AE_sampling$AE),AE_all=mean(AE_all$AE),
                  cor_0=mean(cor_0),cor_sampling=mean(cor_sampling),cor_all=mean(cor_all))
  out2=data.frame(PTSZ=sd(PTSZ),PTDO=sd(PTDO),MSE=sd(MSE),
                  RE_sampling=sd(RE_sampling$RE),RE_all=sd(RE_all$RE),
                  AE_0=sd(AE_0$AE),AE_sampling=sd(AE_sampling$AE),AE_all=sd(AE_all$AE),
                  cor_0=sd(cor_0),cor_sampling=sd(cor_sampling),cor_all=sd(cor_all))
  rownames(out1)=NULL
  rownames(out2)=NULL
  output = list(out1,out2)
  names(output) = c("mean", "sd")
  return(output)
}
