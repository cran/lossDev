\name{stochasticInflation}
\alias{stochasticInflation}
\title{A generic function to plot and/or return predicted and forecast stochastic inflation rates for models in lossDev.}
\description{A generic function to plot and/or return predicted and forecast stochastic inflation rates for models in \pkg{lossDev}.}
\details{If the model incorporates a stochastic rate of inflation, then that rate is assumed to follow (on the log scale) an autoregressive process of order 1.
(The autoregressive process of order 1 is the discrete equivalent to an Ornstein-Uhlenbeck process.)
This function plots the median of the posterior predictive distribution for stochastic inflation (not on the log scale) rates by year.
Values are returned prior to the application of any limits or weights.
Note that for years where observed values are supplied, the model takes those values at face value.}
\value{Mainly called for the side effect of plotting.}
\seealso{\code{\link[=stochasticInflation,AnnualAggLossDevModelOutput-method]{stochasticInflation("AnnualAggLossDevModelOutput")}}
\code{\link{stochasticInflationRhoParameter}}
\code{\link{stochasticInflationStationaryMean}}}
\arguments{\item{object}{The object from which to plot and/or return the stochastic inflation rates.}
\item{extraYears}{An integer expressing the (maximum) number of years to plot (beyond the final observed year).  Must be at least zero.  Default is 15.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
