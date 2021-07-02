## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(SRS)
library(mclust)

## -----------------------------------------------------------------------------
options(digits = 2)
library(SRS)
data(simudat)
head(simudat)

## -----------------------------------------------------------------------------
data("simudat")
#simudat_res=SRS(simudat, windowsize=2, nbins=61, epochs = 100, estimates.only = FALSE)

## -----------------------------------------------------------------------------
data("simudat_res")
head(simudat_res)

## ---- warning = FALSE, message = FALSE----------------------------------------
  Psz=GM_Psz(simudat, simudat_res, nmix=4)
  head(Psz)

## ---- warning = FALSE, message = FALSE----------------------------------------
par(mar = c(0.4,0.4,0.4,0.4))
par(mfrow=c(1,2))
hm(simudat[,1], 61)
hm(simudat_res[,1], 61)

## ---- warning = FALSE, message = FALSE----------------------------------------
data("str1")
head(str1)

## ---- warning = FALSE, message = FALSE----------------------------------------
set.seed(1234)
#Generate 100 random type1 single cells
data <- generate_single(data=str1, alpha_0=5.6,alpha_1=-1, beta_l=0.9,beta_g=0.9,beta_m=0.9, alpha=0.2, n_single=10) 

## ---- warning = FALSE, message = FALSE----------------------------------------
data("simudat_true")
options(digits = 2)
summa(simudat_true, simudat_res,simudat, Psz)

