\name{triResi}
\alias{triResi}
\title{A generic function to plot and/or return residuals for models in the lossDev package.}
\description{A generic function to plot and/or return residuals for models in the \pkg{lossDev} package.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=triResi,AnnualAggLossDevModelOutput-method]{triResi("AnnualAggLossDevModelOutput")}}
\code{\link{QQPlot}}}
\arguments{\item{object}{The object from which to plot and/or return the residuals.}
\item{standardize}{A logical value.  If \code{TRUE}, the plotted and returned residuals are normalized to their respective standard deviation.}
\item{timeAxis}{A character value describing along which of the three time axes to plot the residuals: \sQuote{dy} for development year time, \sQuote{cy} for calendar year time, \sQuote{ey} for exposure year time.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
