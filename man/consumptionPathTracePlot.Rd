\name{consumptionPathTracePlot}
\alias{consumptionPathTracePlot}
\title{A generic function to generate the trace plots for select consumption path values.}
\description{A generic function to generate the trace plots for select consumption path values.}
\value{NULL invisibly.  Only called for the side effect of plotting.}
\seealso{\code{\link[=consumptionPathTracePlot,StandardAnnualAggLossDevModelOutput-method]{consumptionTracePlot("StandardAnnualAggLossDevModelOutput")}}
\code{\link{consumptionPath}}}
\arguments{\item{object}{The object from which to generate the trace plots.}
\item{elements}{A numeric vector indicating the elements for which to plot the trace.  Valid values are 1 through the number of development years (columns) in the observed triangle.  If NULL, values are selected automatically.}
\item{\dots}{Additional arguments used by methods.}}
