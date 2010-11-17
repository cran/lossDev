\name{firstYearInNewRegime}
\alias{firstYearInNewRegime}
\title{A generic function to plot and/or return the posterior change point.}
\description{A generic function to plot and/or return the posterior change point.}
\details{When incorporating a structural break, the user has the option of specifying either 1) the first year in which the new regime applies or 2) a (inclusive) range in which the first year in the new regime applies.
If the user specifies a range, the actual year is estimated as a model parameter.
This function allows for the retrieval/illustration of the posterior for this estimate.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=firstYearInNewRegime,BreakAnnualAggLossDevModelOutput-method]{firstYearInNewRegime("BreakAnnualAggLossDevModelOutput")}}
\code{\link{firstYearInNewRegimeTracePlot}}}
\arguments{\item{object}{The object from which to plot and/or return the posterior change point estimate.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
