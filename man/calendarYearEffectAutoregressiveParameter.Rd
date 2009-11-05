\name{calendarYearEffectAutoregressiveParameter}
\alias{calendarYearEffectAutoregressiveParameter}
\title{A generic function to plot and/or return the posterior of the autoregressive parameter for the calendar year effect for models in lossDev.}
\description{A generic function to plot and/or return the posterior of the autoregressive parameter for the calendar year effect for models in \pkg{lossDev}.}
\value{Mainly called for the side effect of plotting.}
\docType{genericFunction}
\seealso{\code{\link[=calendarYearEffectAutoregressiveParameter,AnnualAggLossDevModelOutput-method]{calendarYearEffectAutoregressiveParameter("AnnualAggLossDevModelOutput")}}}
\seealso{\code{\link{standardDeviationOfCalendarYearEffect}}}
\seealso{\code{\link{calendarYearEffect}}}
\seealso{\code{\link{calendarYearEffectErrors}}}
\arguments{\item{object}{The object from which to plot and/or return the autoregressive parameter which is associated with the calendar year effect.}
\item{plotDensity}{A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}
\item{plotTrace}{A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.}}
