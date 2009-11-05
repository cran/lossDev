\name{calendarYearEffectErrors,AnnualAggLossDevModelOutput-method}
\alias{calendarYearEffectErrors,AnnualAggLossDevModelOutput-method}
\title{A method to plot and/or return predicted and forecast calendar year effect errors for models in lossDev.}
\description{A method to plot and/or return predicted and forecast calendar year effect errors for models in \pkg{lossDev}.}
\value{Mainly called for the side effect of plotting.  Also returns a named array with the median predicted errors (not on the log scale).  Returned invisibly.}
\docType{methods}
\seealso{\code{\link{calendarYearEffectErrors}}}
\seealso{\code{\link{calendarYearEffect}}}
\seealso{\code{\link{autoregressiveParameter}}}
\seealso{\code{\link{standardDeviationOfCalendarYearEffect}}}
\seealso{\code{\link{calendarYearEffectErrorTracePlot}}}
\arguments{\item{object}{The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the calendar year effect errors.}
\item{extraYears}{An integer expressing the (maximum) number of years to plot (beyond the final observed calendar year).  Must greater than or equal to zero.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
