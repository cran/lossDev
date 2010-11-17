\name{calendarYearEffectErrorTracePlot}
\alias{calendarYearEffectErrorTracePlot}
\title{A generic function to generate the trace plots for select calendar year effect errors.}
\description{A generic function to generate the trace plots for select calendar year effect errors.}
\details{The calendar year effect is comprised of two components: 1) a prior expected value that may be unique to every cell and 2) a diagonal-specific error term.
This function generates trace plots for the diagonal specific error terms only.}
\value{NULL invisibly.  Only called for the side effect of plotting.}
\seealso{\code{\link[=calendarYearEffectErrorTracePlot,AnnualAggLossDevModelOutput-method]{calendarYearEffectErrorTracePlot("AnnualAggLossDevModelOutput")}}
\code{\link{calendarYearEffectErrors}}}
\arguments{\item{object}{The object from which to generate the trace plots.}
\item{elements}{A numeric vector indicating the elements for which to plot the trace.  Valid values are 2 through the total number of exposure years(observed and forecast).  If NULL, values are selected automatically.}}
