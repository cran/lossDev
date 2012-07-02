\name{consumptionPath,StandardAnnualAggLossDevModelOutput-method}
\alias{consumptionPath,StandardAnnualAggLossDevModelOutput-method}
\title{A method to plot and/or return the estimated consumption path vs development year time for standard models.}
\description{A method to plot and/or return the estimated consumption path vs development year time for standard models.}
\details{At the heart of aggregate loss development models in \pkg{lossDev} is the consumption path.
The consumption path is (on a log scale) the trajectory of incremental payments absent any calendar year effects and wisht exposure normalized to the first row.
Note that the measurement error term is (possibly) a skewed \eqn{t} and as such (possibly) has a non zero mean.   The consumption path is absent any such shifts due to skewness.
This is a method that allows for the retrieval and illustration of this consumption path.

The standard model has a common consumption path for all exposure years.

Because the model is Bayesian, the estimated consumption path comes as a distribution; only the median is plotted and/or returned.}
\value{Mainly called for the side effect of plotting.  Also returns the plotted statistics.  Returned invisibly.}
\docType{methods}
\seealso{\code{\link{consumptionPath}}
\code{\link[=consumptionPath,BreakAnnualAggLossDevModelOutput-method]{consumptionPath("BreakAnnualAggLossDevModelOutput")}}}
\arguments{\item{object}{The object from which to plot and/or return the estimated consumption path.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
