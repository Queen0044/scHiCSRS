% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summa3.R
\name{summa3}
\alias{summa3}
\title{summa3}
\usage{
summa3(single, truecount, result)
}
\arguments{
\item{single}{Observed single cells matrix with each column being the upper triangular of a single cell.}

\item{truecount}{Underline true counts from simulation.}

\item{result}{Result from SRS funtion.}
}
\value{
A vector of PTSZ and its SD when fixing PTDO to be 0.80, and the threshold used in that case.
}
\description{
This function calculates PTDO when fix PTDO=0.80.
}
\examples{
summa3(singledat, truecount, res)
}
