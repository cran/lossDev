\name{consumptionPath}
\alias{consumptionPath}
\title{A generic function to plot and/or return the estimated consumption path vs development year time.}
\description{A generic function to plot and/or return the estimated consumption path vs development year time.}
\details{At the heart of aggregate loss development models in \pkg{lossDev} is the consumption path.
The consumption path is (on a log scale) the trajectory of incremental payments absent calendar year effects and with exposure normalized to the first row.
Note that the measurement error term is (possibly) a skewed \eqn{t} and as such (possibly) has a non zero mean.   The consumption path is absent any such shifts due to skewness.
This is a generic function that allows for the retrieval and illustration of this consumption path.}
\value{Mainly called for the side effect of plotting.  Also returns the plotted statistics.  Returned invisibly.}
\seealso{\code{\link[=consumptionPath,StandardAnnualAggLossDevModelOutput-method]{consumptionPath("StandardAnnualAggLossDevModelOutput")}}
\code{\link[=consumptionPath,BreakAnnualAggLossDevModelOutput-method]{consumptionPath("BreakAnnualAggLossDevModelOutput")}}
\code{\link{consumptionPathTracePlot}}}
\arguments{\item{object}{The object from which to plot and/or return the estimated consumption path.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
