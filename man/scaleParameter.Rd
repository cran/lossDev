\name{scaleParameter}
\alias{scaleParameter}
\title{A generic function to plot and/or return the posterior of the scale parameter for the Student-t measurement equation for models in lossDev.}
\description{A generic function to plot and/or return the posterior of the scale parameter for the Student-\eqn{t} measurement equation for models in \pkg{lossDev}.}
\details{As the degrees of freedom of the \eqn{t} goes to infinity, the scale parameter is the standard deviation of the resulting normal distribution (assuming zero skew).}
\value{Mainly called for the side effect of plotting.}
\docType{genericFunction}
\seealso{\code{\link[=scaleParameter,AnnualAggLossDevModelOutput-method]{scaleParameter("AnnualAggLossDevModelOutput")}}}
\arguments{\item{object}{The object from which to plot and/or return the scale parameter.}
\item{column}{The scale parameter is allowed to vary with development time. Setting \code{column} results in the plotting and returning of the scale parameter corresponding to that column. Default value is \code{1}.}
\item{plotDensity}{A logical value. If \code{TRUE}, then the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}
\item{plotTrace}{A logical value. If \code{TRUE}, then the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}}
