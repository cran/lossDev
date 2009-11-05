\name{rateOfDecay,BreakAnnualAggLossDevModelOutput-method}
\alias{rateOfDecay,BreakAnnualAggLossDevModelOutput-method}
\title{A method to plot and/or return the esimtated rate of decay vs development year time for break models.}
\description{A method to plot and/or return the esimtated rate of decay vs development year time for break models.}
\details{The simplest definition of the rate of decay is the exponentiated first difference of the \link[=consumptionPath]{consumption path}.
The break model has two rates of decay.  One which applies to exposure years prior to a structural break.  And another which applies after the break.
This is a method to allow for the retrieval and illustration of these rates of decay.

Because the model is Bayesian, the estimated rates of decay come as distributions; only the medians are plotted and/or returned.}
\value{Mainly called for the side effect of plotting.  Also returns the plotted statistics.  Returned invisibly.}
\docType{methods}
\seealso{\code{\link{rateOfDecay}}}
\seealso{\code{\link[=rateOfDecay,StandardAnnualAggLossDevModelOutput-method]{rateOfDecay("StandardAnnualAggLossDevModelOutput")}}}
\seealso{\code{\link{consumptionPath}}}
\arguments{\item{object}{The object from which to plot and/or return the estimated rate of decay.}
\item{plot}{A logical value. If \code{TRUE}, then the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
