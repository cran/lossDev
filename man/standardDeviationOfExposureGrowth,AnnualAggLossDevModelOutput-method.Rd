\name{standardDeviationOfExposureGrowth,AnnualAggLossDevModelOutput-method}
\alias{standardDeviationOfExposureGrowth,AnnualAggLossDevModelOutput-method}
\title{A method to plot and/or return the posterior of the standard deviation of rhe exposure growth rate for models in lossDev.}
\description{A method to plot and/or return the posterior of the standard deviation of rhe exposure growth rate for models in \pkg{lossDev}.}
\value{Mainly called for the side effect of plotting.  Also returns a named array with select quantiles of the posterior for the standard deviation of exposure growth.  Returned invisibly.}
\docType{methods}
\seealso{\code{\link{standardDeviationOfExposureGrowth}}}
\seealso{\code{\link{exposureGrowth}}}
\seealso{\code{\link{meanExposureGrowth}}}
\arguments{\item{object}{The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the standard deviation of the exposure growth rate.}
\item{plotDensity}{A logical value. If \code{TRUE}, then the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.}
\item{plotTrace}{A logical value. If \code{TRUE}, then the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.}}
