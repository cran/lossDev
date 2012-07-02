\name{stochasticInflationRhoParameter}
\alias{stochasticInflationRhoParameter}
\title{A generic function to plot and/or return the posterior of the stochastic inflation rho parameter for models in lossDev.}
\description{A generic function to plot and/or return the posterior of the stochastic inflation rho parameter for models in \pkg{lossDev}.}
\details{If the model incorporates a stochastic rate of inflation, then that rate is assumed to follow (on the log scale) an autoregressive process of order 1.
(The autoregressive process of order 1 is the discrete equivalent to an Ornstein-Uhlenbeck process.)
This function plots the posterior for the \eqn{rho} parameter, assuming one was estimated.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=stochasticInflationRhoParameter,AnnualAggLossDevModelOutput-method]{stochasticInflationRhoParameter("AnnualAggLossDevModelOutput")}}
\code{\link{stochasticInflationStationaryMean}}
\code{\link{stochasticInflation}}}
\arguments{\item{object}{The object from which to plot and/or return the stochastic inflation \eqn{rho} parameter.}
\item{plotDensity}{A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.}
\item{plotTrace}{A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.}}
