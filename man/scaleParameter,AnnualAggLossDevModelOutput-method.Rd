\name{scaleParameter,AnnualAggLossDevModelOutput-method}
\alias{scaleParameter,AnnualAggLossDevModelOutput-method}
\title{A method to plot and/or return the posterior of the scale parameter for the Student-t measurement equation for models in lossDev.}
\description{A method to plot and/or return the posterior of the scale parameter for the Student-\eqn{t} measurement equation for models in \pkg{lossDev}.}
\value{Mainly called for the side effect of plotting.  Also returns a named array with select quantiles of the posterior for the scale parameter.  Returned invisibly.}
\docType{methods}
\seealso{\code{\link{scaleParameter}}}
\arguments{\item{object}{The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the scale parameter.}
\item{column}{The scale parameter is allowed to vary with development time. Setting \code{column} results in the plotting and returning of the scale parameter corresponding to that column. Default value is \code{1}.}
\item{plotDensity}{A logical value. If \code{TRUE}, then the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}
\item{plotTrace}{A logical value. If \code{TRUE}, then the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}}
