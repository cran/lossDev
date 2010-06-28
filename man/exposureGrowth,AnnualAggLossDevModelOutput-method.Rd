\name{exposureGrowth,AnnualAggLossDevModelOutput-method}
\alias{exposureGrowth,AnnualAggLossDevModelOutput-method}
\title{A method to plot and/or return the posterior predicted exposure growth (corresponding to eta in the model).}
\description{A method to plot and/or return the posterior predicted exposure growth (corresponding to \emph{eta} in the model).}
\value{Mainly called for the side effect of plotting the exposure growth.  Also returns a named numeric vector for the median of the posterior for the exposure growth on the real (not log) scale.  Returned invisibly.}
\docType{methods}
\seealso{\code{\link{exposureGrowth}}
\code{\link{exposureGrowthTracePlot}}}
\arguments{\item{object}{The object from which to plot and/or return the posterior predicted exposure growth.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
