\name{stochasticInflation,AnnualAggLossDevModelOutput-method}
\alias{stochasticInflation,AnnualAggLossDevModelOutput-method}
\title{A method to plot and/or return predicted and forecast stochastic inflation rates for models in lossDev.}
\description{A method to plot and/or return predicted and forecast stochastic inflation rates for models in \pkg{lossDev}.}
\value{Mainly called for the side effect of plotting.  Also returns a named array of the median predicted inflation rate (not on the log scale).  Returned invisibly.}
\docType{methods}
\seealso{\code{\link{stochasticInflation}}
\code{\link{stochasticInflationRhoParameter}}
\code{\link{stochasticInflationStationaryMean}}}
\arguments{\item{object}{The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return predicted and forecast stochastic inflation rates.}
\item{extraYears}{An integer expressing the (maximum) number of years to plot (beyond the final observed year).  Must be at least zero.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
