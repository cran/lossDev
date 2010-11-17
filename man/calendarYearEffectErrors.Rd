\name{calendarYearEffectErrors}
\alias{calendarYearEffectErrors}
\title{A generic function to plot and/or return predicted and forecast calendar year effect errors for models in lossDev.}
\description{A generic function to plot and/or return predicted and forecast calendar year effect errors for models in \pkg{lossDev}.}
\details{The calendar year effect is comprised of two components: 1) a prior expected value which may be unique to every cell (subject to weights and bounds) and 2) a diagonal-specific error term.
This function only plots and returns the error term, which includes an autoregressive component if the model is estimated with such a feature.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=calendarYearEffectErrors,AnnualAggLossDevModelOutput-method]{calendarYearEffectErrors("AnnualAggLossDevModelOutput")}}
\code{\link{calendarYearEffect}}
\code{\link{autoregressiveParameter}}
\code{\link{standardDeviationOfCalendarYearEffect}}
\code{\link{calendarYearEffectErrorTracePlot}}}
\arguments{\item{object}{The object from which to plot and/or return the calendar year effect errors.}
\item{extraYears}{An integer expressing the (maximum) number of years to plot (beyond the final observed calendar year).  Must be greater than or equal to zero.  Default is 15.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
