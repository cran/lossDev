\name{rateOfDecayTracePlot}
\alias{rateOfDecayTracePlot}
\title{A generic function to plot the trace plots for select rate of decay values.}
\description{A generic function to plot the trace plots for select rate of decay values.}
\value{NULL invisibly.  Only called for the side effect of plotting.}
\docType{genericFunction}
\seealso{\code{\link[=rateOfDecayTracePlot,StandardAnnualAggLossDevModelOutput-method]{rateOfDecayTracePlot("StandardAnnualAggLossDevModelOutput")}}}
\seealso{\code{\link{rateOfDecay}}}
\arguments{\item{object}{The object from which to generate the trace plots.}
\item{elements}{A numeric vector indicating for which elements to plot the trace.  Valid values are 2 through the number of columns in the observed triangle.  If NULL, values are selected automatically.}
\item{\dots}{Additional arguments used by methods.}}
