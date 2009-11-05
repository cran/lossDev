\name{degreesOfFreedom}
\alias{degreesOfFreedom}
\title{A generic function to plot and/or return the posterior of the degrees of freedom for the Student-t in  models in lossDev.}
\description{A generic function to plot and/or return the posterior of the degrees of freedom for the Student-\eqn{t} in  models in \pkg{lossDev}.}
\details{When there is zero skew, the degrees of freedom are the degrees of freedom for the non-skewed \eqn{t}.}
\references{Kim, Y., and J. McCulloch (2007) \dQuote{The Skew-Student Distribution with Application to U.S. Stock Market Returns and the Equity Premium,} Department of Economics, Ohio State University, October 2007.}
\value{Mainly called for the side effect of plotting.}
\docType{genericFunction}
\seealso{\code{\link[=degreesOfFreedom,AnnualAggLossDevModelOutput-method]{degreesOfFreedom("AnnualAggLossDevModelOutput")}}}
\arguments{\item{object}{The object from which to plot and/or return the degrees of freedom.}
\item{plotDensity}{A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}
\item{plotTrace}{A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}}
