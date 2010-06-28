\name{predictedPayments,AnnualAggLossDevModelOutput-method}
\alias{predictedPayments,AnnualAggLossDevModelOutput-method}
\title{A method to plot predicted vs actual payments for models from the lossDev package.}
\description{A method to plot predicted vs actual payments for models from the \pkg{lossDev} package.}
\details{Because the model is Bayesian, each estimated payment comes as a distribution.
The median of this distribution is used as a point estimate when plotting and/or returning values.
Note: One cannot calculate the estimated incremental payments from the estimated cumulative payments (and vice versa) since the median of sums need not be equal to the sum of medians.

If \code{mergePredictedWithObserved=TRUE} and \code{type="incremental"}, then any observed incremental payment will be used in place of its corresponding incremental payment.
If \code{mergePredictedWithObserved=TRUE} and \code{type="cumulative"}, then only predicted incremental payments (by row) to the right of the last observed cumulative value will enter the calculation.}
\value{Mainly called for the side effect of plotting.  Also returns a named array (with the same structure as the input triangle) containing the predicted log incremental payments.  Returned invisibly.}
\docType{methods}
\seealso{\code{\link{predictedPayments}}}
\arguments{\item{object}{The object of type \code{AnnualAggLossDevModelOutput} from which to plot predicted vs actual payments and return predicted payments.}
\item{type}{A singe character value specifying whether to plot/return the predicted incremental or cumulative payments. Valid values are "incremental" or "cumulative."  See details as to why these may not match up.}
\item{logScale}{A logical value.  If \code{TRUE}, then values are plotted on a log scale.}
\item{mergePredictedWithObserved}{A logical value.  If \code{TRUE}, then the returned values treat observed incremental payments at "face value"; otherwise predicted values are used in place of observed values.}
\item{plotObservedValues}{A logical value.  If \code{FALSE}, then only the predicted values are plotted.}
\item{plotPredictedOnlyWhereObserved}{A logical value.  See details.}
\item{quantiles}{A vector of quantiles for the predicted payments to return.  Usefull for constructing credible intervals.}
\item{plot}{A logical value. If \code{TRUE}, then the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
