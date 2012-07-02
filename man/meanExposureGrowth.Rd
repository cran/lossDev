\name{meanExposureGrowth}
\alias{meanExposureGrowth}
\title{A generic function to plot and/or return the posterior of the mean exposure growth for models in lossDev.}
\description{A generic function to plot and/or return the posterior of the mean exposure growth for models in \pkg{lossDev}.}
\details{(Optionally) exposure growth is modeled as an ar1 process.  This inherently assumes that periods of high exposure growth are (or at least have the possibility of being) followed by continued high periods.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=meanExposureGrowth,AnnualAggLossDevModelOutput-method]{meanExposureGrowth("AnnualAggLossDevModelOutput")}}}
\arguments{\item{object}{The object from which to plot and/or return the mean exposure growth.}
\item{plotDensity}{A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}
\item{plotTrace}{A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}}
