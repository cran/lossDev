\name{rateOfDecay}
\alias{rateOfDecay}
\title{A generic function to plot and/or return the esimtated rate of decay vs development year time.}
\description{A generic function to plot and/or return the esimtated rate of decay vs development year time.}
\details{The simplest definition of the rate of decay is the exponentiated first difference of the \link[=consumptionPath]{consumption path}.
This is a generic function to allow for the retrieval and illustration of the rate of decay.}
\value{Mainly called for the side effect of plotting.  Also returns the plotted statistics.  Returned invisibly.}
\docType{genericFunction}
\seealso{\code{\link[=rateOfDecay,StandardAnnualAggLossDevModelOutput-method]{rateOfDecay("StandardAnnualAggLossDevModelOutput")}}}
\seealso{\code{\link[=rateOfDecay,BreakAnnualAggLossDevModelOutput-method]{rateOfDecay("BreakAnnualAggLossDevModelOutput")}}}
\seealso{\code{\link{consumptionPath}}}
\seealso{\code{\link{rateOfDecayTracePlot}}}
\arguments{\item{object}{The object from which to plot and/or return the estimated rate of decay.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
