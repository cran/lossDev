\name{calendarYearEffect,AnnualAggLossDevModelOutput-method}
\alias{calendarYearEffect,AnnualAggLossDevModelOutput-method}
\title{A method to plot and/or return predicted and forecast calendar year effects for models in lossDev.}
\description{A method to plot and/or return predicted and forecast calendar year effects for models in \pkg{lossDev}.}
\details{The first cell is \code{NA}. Values in the first column represent the rate of inflation/escalation to the corresponding cell from the cell in the same column but previous row.
Values in the 2nd column and beyond represent the rate of inflation/escalation to the corresponding cell from the cell in the same row but previous column.}
\value{Mainly called for the side effect of plotting.  Also returns a named array with the median predicted values (not on the log scale).  Returned invisibly.}
\docType{methods}
\seealso{\code{\link{calendarYearEffect}}
\code{\link{calendarYearEffectErrors}}
\code{\link{autoregressiveParameter}}
\code{\link{standardDeviationOfCalendarYearEffect}}
\code{\link{calendarYearEffectErrorTracePlot}}}
\arguments{\item{object}{The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the calendar year effect.}
\item{restrictedSize}{A logical value.  If \code{TRUE}, the plotted calendar year effect is restricted to the square of dimension equal to the observed triangle with which the model was estimated.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
