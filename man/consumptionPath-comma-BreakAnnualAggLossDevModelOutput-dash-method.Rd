\name{consumptionPath,BreakAnnualAggLossDevModelOutput-method}
\alias{consumptionPath,BreakAnnualAggLossDevModelOutput-method}
\title{A method to plot and/or return the estimated consumption path vs development year time for break models.}
\description{A method to plot and/or return the estimated consumption path vs development year time for break models.}
\details{At the heart of aggregate loss development models in \pkg{lossDev} is the consumption path.
The consumption path is (on a log scale) the trajectory of incremental payments absent any calendar year effects and normalized to the first row.
Note that the measurement error term is (possibly) a skewed \eqn{t} and as such (possibly) has a non zero mean.   The consumption path is absent any such shifts in due to skewness.
This is a method to allow for the retrieval and illustration of this consumption path.

The break model has a two consumption paths: exposure years prior to the structural break follow one path and exposure years after the break follow another.
The slope of the consumption path for exposure years prior to the break is used to extend the consuption path for exposure years post break to the end of the triangle.

Because the model is Bayesian, the estimated consumption paths come as distributions; only the medians are plotted and/or returned.}
\value{Mainly called for the side effect of plotting.  Also returns the plotted statistics.  Returned invisibly.}
\docType{methods}
\seealso{\code{\link{consumptionPath}}
\code{\link[=consumptionPath,StandardAnnualAggLossDevModelOutput-method]{consumptionPath("StandardAnnualAggLossDevModelOutput")}}}
\arguments{\item{object}{The object from which to plot and/or return the estimated consumption paths.}
\item{plot}{A logical value. If \code{TRUE}, then the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
