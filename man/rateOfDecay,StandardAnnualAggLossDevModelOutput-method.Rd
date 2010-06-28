\name{rateOfDecay,StandardAnnualAggLossDevModelOutput-method}
\alias{rateOfDecay,StandardAnnualAggLossDevModelOutput-method}
\title{A method to plot and/or return the esimtated rate of decay vs development year time for standard models.}
\description{A method to plot and/or return the esimtated rate of decay vs development year time for standard models.}
\details{The simplest definition of the rate of decay is the exponentiated first difference of the \link[=consumptionPath]{consumption path}.
The standard model has a common rate of decay for all exposure years.
This is a method to allow for the retrieval and illustration of the rate of decay.

Because the model is Bayesian, the estimated rate of decay comes as a distribution; only the median is plotted and/or returned.}
\value{Mainly called for the side effect of plotting.  Also returns the plotted statistics.  Returned invisibly.}
\docType{methods}
\seealso{\code{\link{rateOfDecay}}
\code{\link[=rateOfDecay,BreakAnnualAggLossDevModelOutput-method]{rateOfDecay("BreakAnnualAggLossDevModelOutput")}}
\code{\link{consumptionPath}}}
\arguments{\item{object}{The object from which to plot and/or return the estimated rate of decay.}
\item{plot}{A logical value. If \code{TRUE}, then the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
