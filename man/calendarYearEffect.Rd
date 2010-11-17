\name{calendarYearEffect}
\alias{calendarYearEffect}
\title{A generic function to plot and/or return the predicted and forecast calendar year effects for models in lossDev.}
\description{A generic function to plot and/or return the predicted and forecast calendar year effects for models in \pkg{lossDev}.}
\details{The calendar year effect is comprised of two components: 1) a prior expected value that may be unique to every cell (subject to weights and bounds) and 2) a diagonal-specific error term.
This function plots and returns the factor resulting from the combined effect of these two, which includes an autoregressive component if the model is estimated with such a feature.

The first cell is \code{NA}. Values in the first column represent the rate of inflation/escalation to the corresponding cell from the cell in the same column but previous row.
Values in the 2nd column and beyond represent the rate of inflation/escalation to the corresponding cell from the cell in the same row but previous column.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=calendarYearEffect,AnnualAggLossDevModelOutput-method]{calendarYearEffect("AnnualAggLossDevModelOutput")}}
\code{\link{calendarYearEffectErrors}}
\code{\link{autoregressiveParameter}}
\code{\link{standardDeviationOfCalendarYearEffect}}
\code{\link{calendarYearEffectErrorTracePlot}}}
\arguments{\item{object}{The object from which to plot and/or return the calendar year effect.}
\item{restrictedSize}{A logical value.  If \code{TRUE}, the plotted calendar year effect is restricted to the square of dimension equal to the observed triangle with which the model was estimated.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
