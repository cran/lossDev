\name{stochasticInflationStationaryMean,AnnualAggLossDevModelOutput-method}
\alias{stochasticInflationStationaryMean,AnnualAggLossDevModelOutput-method}
\title{A method to plot and/or return the posterior of the stochastic inflation stationary mean for models in lossDev.}
\description{A method to plot and/or return the posterior of the stochastic inflation stationary mean for models in \pkg{lossDev}.}
\value{Mainly called for the side effect of plotting.  Also returns a named array with select quantiles of the stochastic inflation stationary mean.  Returned invisibly.}
\docType{methods}
\seealso{\code{\link{stochasticInflationStationaryMean}}
\code{\link{stochasticInflationRhoParameter}}
\code{\link{stochasticInflation}}}
\arguments{\item{object}{The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the stochastic inflation stationary mean.}
\item{plotDensity}{A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}
\item{plotTrace}{A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}}
