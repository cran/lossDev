\name{gompertzParameters,AnnualAggLossDevModelOutputWithZeros-method}
\alias{gompertzParameters,AnnualAggLossDevModelOutputWithZeros-method}
\title{A method to plot and/or return the posterior of the parameters for the gompertz curve which describes the probability of payment.}
\description{A method to plot and/or return the posterior of the parameters for the gompertz curve which describes the probability of payment.}
\details{The scale parameter describes how steep the curve is.
Larger values are steeper.
Positive values indicate that the probability of a positive payment should decrease with development time.
(The scale is restricted to be positive.)

The fifty.fifty parameter gives the point (in development time) when the gompertz curve gives a probability of fifty percent.}
\value{Mainly called for the side effect of plotting.}
\docType{genericFunction}
\seealso{\code{\link{gompertzParameters}}}
\arguments{\item{object}{The object from which to plot and/or return the parameters.}
\item{parameter}{A character describing which parameter to plot. \dQuote{scale} for the scale parameter. \dQuote{fifty.fifty} for the point at which the gompertz give a probably of fifty percent.}
\item{plotDensity}{A logical value. If \code{TRUE}, then the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}
\item{plotTrace}{A logical value. If \code{TRUE}, then the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}}
