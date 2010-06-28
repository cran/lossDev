\name{stochasticInflationStationaryMean}
\alias{stochasticInflationStationaryMean}
\title{A generic function to plot and/or return the posterior of the stochastic inflation stationary mean for models in lossDev.}
\description{A generic function to plot and/or return the posterior of the stochastic inflation stationary mean for models in \pkg{lossDev}.}
\details{If the model incorporates a stochastic rate of inflation, then that rate is assumed to follow (on the log scale) an autoregressive process of order 1.
(The autoregressive process of order 1 is the discrete equivalent to an Ornstein-Uhlenbeck process.)
This function plots the posterior for the stationary mean (on the log scale), assuming such a mean was estimated.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=stochasticInflationStationaryMean,AnnualAggLossDevModelOutput-method]{stochasticInflationStationaryMean("AnnualAggLossDevModelOutput")}}
\code{\link{stochasticInflationRhoParameter}}
\code{\link{stochasticInflation}}}
\arguments{\item{object}{The object from which to plot and/or return the stochastic inflation stationary mean.}
\item{plotDensity}{A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}
\item{plotTrace}{A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}}
