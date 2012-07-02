\name{skewnessParameter,AnnualAggLossDevModelOutput-method}
\alias{skewnessParameter,AnnualAggLossDevModelOutput-method}
\title{A method to plot and/or return the posterior of the skewness parameter for models in lossDev.}
\description{A method to plot and/or return the posterior of the skewness parameter for models in \pkg{lossDev}.}
\details{The skewness parameter does not directly correspond to the degree of skewness. However, all else being equal, a larger (in magnitude) skewness parameter indicates a higher degree of skewness,
and a skewness parameter of zero equates to zero skew.}
\references{Kim, Y., and J. McCulloch (2007) \dQuote{The Skew-Student Distribution with Application to U.S. Stock Market Returns and the Equity Premium,} Department of Economics, Ohio State University, October 2007}
\value{Mainly called for the side effect of plotting.  But also returns a named array with some select quantiles of the posterior for the skewness parameter.  Returned invisibly.}
\docType{methods}
\seealso{\code{\link{skewnessParameter}}}
\arguments{\item{object}{The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the skewness parameter.}
\item{plotDensity}{A logical value. If \code{TRUE}, then the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}
\item{plotTrace}{A logical value. If \code{TRUE}, then the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}}
