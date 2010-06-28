\name{lossDevelopmentFactors}
\alias{lossDevelopmentFactors}
\title{A generic function to plot and/or return a table of predicted age-to-age loss development factors (or link ratios).}
\description{A generic function to plot and/or return a table of predicted age-to-age loss development factors (or link ratios).}
\details{While the model estimates ultimate losses directly, comparisons of predicted to observed development factors can give the user a better feel for the model's adequacy.
Since the model is Bayesian, each development factor comes as a distribution.  Only the median, as a point estimate, are plotted/returned.

The age-to-age factors are the ratios of the cumulative paid values at one period to the previous period.
Note that the median of products is not the product of medians, and thus it is possible (or rather likely) that age-to-age factors will not line up with age-to-ultimate factors (see \code{\link{tailFactor}}).}
\value{Mainly called for the side effect of plotting.  Also returns a numeric matrix of plotted statistics.}
\seealso{\code{\link[=lossDevelopmentFactors,AnnualAggLossDevModelOutput-method]{lossDevelopmentFactors("AnnualAggLossDevModelOutput")}}
\code{\link{tailFactor}}}
\arguments{\item{object}{The object from which to plot and/or return loss development factors.}
\item{cex.text}{The \code{cex} value supplied to \code{text}. Adjusts the relative size of text.}
\item{linespace}{Adjusts the spacing between observed and predicted values.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
