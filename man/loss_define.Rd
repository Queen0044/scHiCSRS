% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loss_define.R
\name{loss_define}
\alias{loss_define}
\title{loss_define
Calculate the loss value.}
\usage{
loss_define(y_true, y_pred)
}
\arguments{
\item{y_true}{The true matrix.}

\item{y_pred}{The predicted matrix.}
}
\value{
The sum of squared distance between the predicted and the true matrix.
}
\description{
loss_define
Calculate the loss value.
}
\examples{
loss_define(y_true, y_pred)
}
