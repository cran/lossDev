\name{standardDeviationForScaleInnovation}
\alias{standardDeviationForScaleInnovation}
\title{A generic function to plot and/or return the posterior of the standard deviation for the innovation in the scale parameter for models in lossDev.}
\description{A generic function to plot and/or return the posterior of the standard deviation for the innovation in the scale parameter for models in \pkg{lossDev}.}
\details{Changes in the scale parameter (see \code{\link{scaleParameter}}) are assumed to follow a second-order random walk on the log scale.
This function plots the posterior standard deviation for this random walk.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=standardDeviationForScaleInnovation,AnnualAggLossDevModelOutput-method]{standardDeviationForScaleInnovation("AnnualAggLossDevModelOutput")}}
\code{\link{standardDeviationVsDevelopmentTime}}}
\arguments{\item{object}{The object from which to plot and/or return the standard deviation for the innovation in the log of the scale parameter.}
\item{plotDensity}{A logical value. If \code{TRUE}, then the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}
\item{plotTrace}{A logical value. If \code{TRUE}, then the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.}}
