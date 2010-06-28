\name{exposureGrowthTracePlot}
\alias{exposureGrowthTracePlot}
\title{A generic function to generate the trace plots for select exposure growth rates.}
\description{A generic function to generate the trace plots for select exposure growth rates.}
\value{NULL invisibly.  Only called for the side effect of plotting.}
\seealso{\code{\link[=exposureGrowthTracePlot,AnnualAggLossDevModelOutput-method]{exposureGrowthTracePlot("AnnualAggLossDevModelOutput")}}
\code{\link{exposureGrowth}}}
\arguments{\item{object}{The object from which to generate the trace plots.}
\item{elements}{A numeric vector indicating the elements for which to plot the trace.  Valid values are 2 through the total number of exposure years.  If NULL, values are selected automatically.}}
