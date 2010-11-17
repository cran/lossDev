\name{predictedPayments}
\alias{predictedPayments}
\title{A generic function to plot predicted vs actual payments for models from the lossDev package.}
\description{A generic function to plot predicted vs actual payments for models from the \pkg{lossDev} package.}
\details{Because the model is Bayesian, each estimated payment comes as a distribution.
The median of this distribution is used as a point estimate when plotting and/or returning values.
Note: One cannot calculate the estimated incremental payments from the estimated cumulative payments (and vice versa) since the median of sums need not be equal to the sum of medians.

If \code{mergePredictedWithObserved=TRUE} and \code{type="incremental"}, then any observed incremental payment will be used in place of its corresponding incremental payment.
If \code{mergePredictedWithObserved=TRUE} and \code{type="cumulative"}, then only predicted incremental payments (by row) to the right of the last observed cumulative value will enter the calculation.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=predictedPayments,AnnualAggLossDevModelOutput-method]{predictedPayments("AnnualAggLossDevModelOutput")}}}
\arguments{\item{object}{The object from which to plot predicted vs actual payments and from which to return predicted payments.}
\item{type}{A single character value specifying whether to plot/return the predicted incremental or cumulative payments. Valid values are \dQuote{incremental} or \dQuote{cumulative.}  See details as to why these may not match up.}
\item{logScale}{A logical value.  If \code{TRUE}, then values are plotted on a log scale.}
\item{mergePredictedWithObserved}{A logical value.  See details.}
\item{plotObservedValues}{A logical value.  If \code{FALSE}, then only the predicted values are plotted.}
\item{plotPredictedOnlyWhereObserved}{A logical value.  If \code{TRUE}, then only the predicted incremental payments with valid corresponding observed (log) incremental payment are plotted. Ignored for \code{type="cumulative"}.}
\item{quantiles}{A vector of quantiles for the predicted payments to return.  Useful for constructing credible intervals.}
\item{plot}{A logical value. If \code{TRUE}, then the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
