\name{degreesOfFreedom,AnnualAggLossDevModelOutput-method}
\alias{degreesOfFreedom,AnnualAggLossDevModelOutput-method}
\title{A method to plot and/or return the posterior of the degrees of freedom for the Student-t in  models in lossDev.}
\description{A method to plot and/or return the posterior of the degrees of freedom for the Student-\eqn{t} in  models in \pkg{lossDev}.}
\value{Mainly called for the side effect of plotting.  Also returns a named array with some select quantiles of the posterior for the degrees of freedom.  Returned invisibly.}
\docType{methods}
\seealso{\code{\link{degreesOfFreedom}}}
\arguments{\item{object}{The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the degrees of freedom.}
\item{plotDensity}{A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}
\item{plotTrace}{A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}}
