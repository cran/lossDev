\name{finalCumulativeDiff,AnnualAggLossDevModelOutputWithZeros-method}
\alias{finalCumulativeDiff,AnnualAggLossDevModelOutputWithZeros-method}
\title{A method to plot and/or return the difference between final actual and predicted cumulative payments.}
\description{A method to plot and/or return the difference between final actual and predicted cumulative payments.}
\details{This method accounts for zero payments. By weighting estimated predicted payments by the probably that the payment is greater than zero.

The relative difference (x/y - 1) between the final observed cumulative payment and the corresponding predicted cumulative payment is plotted for each exposure year.
The horizontal lines of each box represent (starting from the top) the 90th, 75th, 50th, 20th, and 10th percentiles.  Exposure years in which all cumulative payments are \code{NA} are omitted.}
\value{Mainly called for the side effect of plotting the difference between final actual and predicted cumulative payments by exposure year.  Also returns a named array for the percentiles in the plot.  Returned invisibly.}
\docType{methods}
\seealso{\code{\link{accountForZeroPayments}}}
\seealso{\code{\link{finalCumulativeDiff}}}
\seealso{\code{\link{finalCumulativeDiff,AnnualAggLossDevModelOutput-method}}}
\arguments{\item{object}{The object of type \code{AnnualAggLossDevModelOuputWithZeros} from which to plot and/or return the difference between final actual and predicted cumulative payments.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}
\item{expYearRange}{Either a range of years (for example c(1995, 2006)) or one of the keywords \dQuote{all} or \dQuote{fullyObs}.}}
