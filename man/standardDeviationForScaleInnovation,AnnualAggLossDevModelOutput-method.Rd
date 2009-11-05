\name{standardDeviationForScaleInnovation,AnnualAggLossDevModelOutput-method}
\alias{standardDeviationForScaleInnovation,AnnualAggLossDevModelOutput-method}
\title{A method to plot and/or return the posterior of the standard deviation for the innovation in the scale parameter for models in lossDev.}
\description{A method to plot and/or return the posterior of the standard deviation for the innovation in the scale parameter for models in \pkg{lossDev}.}
\value{Mainly called for the side effect of plotting.  Also returns a named array with some select quantiles of the posterior for the standard deviation in question.  Returned invisibly.}
\docType{methods}
\seealso{\code{\link{standardDeviationForScaleInnovation}}}
\seealso{\code{\link{scaleParameter}}}
\seealso{\code{\link{standardDeviationVsDevelopmentTime}}}
\arguments{\item{object}{The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the standard deviation for the innovation in the log of the scale parameter.}
\item{plotDensity}{A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.}
\item{plotTrace}{A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.}}
