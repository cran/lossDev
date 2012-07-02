\name{tailFactor,BreakAnnualAggLossDevModelOutput-method}
\alias{tailFactor,BreakAnnualAggLossDevModelOutput-method}
\title{A method to plot and/or return the predicted tail factors for a specific attachment point.}
\description{A method to plot and/or return the predicted tail factors for a specific attachment point.}
\details{The tail factor is the ratio of the estimated ultimate loss to cumulative loss at some point in development time.
This is a method to allow for the retrieval and illustration of the tail factor by exposure year.

Because the model is Bayesian, each tail factor comes as a distribution.  To ease graphical interpretation, only the median for each factor is plotted/returned.
See for more details \code{\link{tailFactor}}.

For comparison purposes, the function returns three separated tail factors for three scenarios.  Theses three tail factors are returned as a list with the following names and meanings:
\describe{
\item{\dQuote{Actual}}{
These are the tail factors estimated when taking the break into consideration.
}
\item{\dQuote{AsIfPostBreak}}{
These are the tail factors estimated when assuming all years where in the post-break regime.
}
\item{\dQuote{AsIfPreBreak}}{
These are the tail factors estimated when assuming all years where in the pre-break regime.
}
}}
\value{Mainly called for the side effect of plotting.  Also returns tail factors for \emph{all} attachment points through \code{finalAttachment}.  See Details. Returned invisibly.}
\docType{methods}
\seealso{\code{\link{tailFactor}}
\code{\link[=tailFactor,StandardAnnualAggLossDevModelOutput-method]{tailFactor("StandardAnnualAggLossDevModelOutput")}}}
\arguments{\item{object}{The object from which to plot the predicted tail factors and return tail factors for \emph{all} attachment points.}
\item{attachment}{An integer value specifying the attachment point for the tail.  Must be at least 1. See Details for more info.}
\item{useObservedValues}{A logical value.  If \code{TRUE}, observed values are substituted for predicted values whenever possible in the calculation.  If \code{FALSE}, only predicted values are used.}
\item{firstIsHalfReport}{A logical value or \code{NA}.  See Details for more information.}
\item{finalAttachment}{An integer value must be at least 1 default value is \code{attachment}.  A call to \code{tailFactor} returns (invisibly) a matrix of tail factors through this value.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}
\item{expYearRange}{Either a range of years (for example c(1995, 2006)) or one of the keywords \dQuote{all} or \dQuote{fullyObs}.}}
