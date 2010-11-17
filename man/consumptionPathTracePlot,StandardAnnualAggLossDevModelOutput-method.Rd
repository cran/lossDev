\name{consumptionPathTracePlot,StandardAnnualAggLossDevModelOutput-method}
\alias{consumptionPathTracePlot,StandardAnnualAggLossDevModelOutput-method}
\title{A method to generate the trace plots for select consumption path values.}
\description{A method to generate the trace plots for select consumption path values.}
\value{NULL invisibly.  Only called for the side effect of plotting.}
\docType{methods}
\seealso{\code{\link{consumptionPathTracePlot}}
\code{\link{consumptionPath}}}
\arguments{\item{object}{The object of type \code{StandardAnnualAggLossDevModelOutput} from which to generate the trace plots.}
\item{elements}{A numeric vector indicating the elements for which to plot the trace.  Valid values are 1 through the number of development years (columns) in the observed triangle.  If NULL, values are selected automatically.}
\item{\dots}{Additional arguments used by other methods.  Not utilized by this method.}}
