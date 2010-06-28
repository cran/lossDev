\name{numberOfKnots,BreakAnnualAggLossDevModelOutput-method}
\alias{numberOfKnots,BreakAnnualAggLossDevModelOutput-method}
\title{A method to plot and/or return the posterior number of knots.}
\description{A method to plot and/or return the posterior number of knots.}
\details{The break model has to consumption paths.  This method will plot the number of knots for each path.}
\value{Mainly called for the side effect of plotting.  Also returns a list of length 2 with each element containing a named vector with names equal to the number of knots and values equal to the density.  Returned invisibly.}
\docType{methods}
\seealso{\code{\link{consumptionPath}}
\code{\link{numberOfKnots}}
\code{\link[=numberOfKnots,StandardAnnualAggLossDevModelOutput-method]{numberOfKnots("StandardAnnualAggLossDevModelOutput")}}}
\arguments{\item{object}{The object from which to plot the number of knots.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
