\name{triResi,AnnualAggLossDevModelOutput-method}
\alias{triResi,AnnualAggLossDevModelOutput-method}
\title{A method to plot and/or return residuals for models in the lossDev package.}
\description{A method to plot and/or return residuals for models in the \pkg{lossDev} package.}
\details{Because the model is Bayesian, each residual comes as a distribution.  To ease graphical interpretation, only the median for each residual is plotted/returned.
The residual is defined as the observed value minus the posterior mean; if standardized, it is also divided by its posterior standard deviation.}
\value{Mainly called for the side effect of plotting.  Also returns a named array with the same structure as the input triangle.  Returned invisibly.}
\docType{methods}
\seealso{\code{\link{triResi}}
\code{\link{QQPlot}}}
\arguments{\item{object}{The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the residuals.}
\item{timeAxis}{A character value describing along which of the three (3) time axis to plot the residuals. \sQuote{dy} for development year time, \sQuote{cy} for calendar year time, \sQuote{ey} for exposure year time.}
\item{standardize}{A logical value.  If \code{TRUE}, the plotted and returned residuals are normalized to their respective standard deviation.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
