\name{rateOfDecayTracePlot,BreakAnnualAggLossDevModelOutput-method}
\alias{rateOfDecayTracePlot,BreakAnnualAggLossDevModelOutput-method}
\title{A method to generate the trace plots for select rate of decay values.}
\description{A method to generate the trace plots for select rate of decay values.}
\value{NULL invisibly.  Only called for the side effect of plotting.}
\docType{methods}
\seealso{\code{\link{rateOfDecayTracePlot}}
\code{\link{rateOfDecay}}}
\arguments{\item{object}{The object of type \code{BreakAnnualAggLossDevModelOutput} from which to generate the trace plots.}
\item{elements}{A numeric vector indicating for which elements to plot the trace.  Valid values are 2 through the number of columns in the observed triangle.  If NULL, values are selected automatically.}
\item{preBreak}{A logical value indicating whether to plot the trace for the pre-break rate(s) of decay or the post-break rate(s) of decay.}
\item{\dots}{Additional arguments used by other methods.  Not utilized by this method.}}
