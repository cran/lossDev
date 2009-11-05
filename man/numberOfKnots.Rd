\name{numberOfKnots}
\alias{numberOfKnots}
\title{A generic function to plot and/or return the posterior number of knots.}
\description{A generic function to plot and/or return the posterior number of knots.}
\details{The \link[=consumptionPath]{consumption path} (or calendar year effect and exposure growth adjusted log incremental payments) is modeled as a linear spline.
The number of knots (or places where the spline changes slope) in this spline is endogenous to the model and estimated by way of Reversible Jump Markov Chain Monte Carlo simulation.}
\value{Mainly called for the side effect of plotting.  Also returns statics on the number of knots.  Returned invisibly.}
\docType{genericFunction}
\seealso{\code{\link{consumptionPath}}}
\seealso{\code{\link[=numberOfKnots,StandardAnnualAggLossDevModelOutput-method]{numberOfKnots("StandardAnnualAggLossDevModelOutput")}}}
\seealso{\code{\link[=numberOfKnots,BreakAnnualAggLossDevModelOutput-method]{numberOfKnots("BreakAnnualAggLossDevModelOutput")}}}
\arguments{\item{object}{The object from which to plot the number of knots.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
