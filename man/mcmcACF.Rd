\name{mcmcACF}
\alias{mcmcACF}
\title{A generic function to plot autocorrelations found in the MCMC samples for select parameters.}
\description{A generic function to plot autocorrelations found in the \acronym{MCMC} samples for select parameters.}
\details{Chains with high autocorrelation require a longer burnin and more samples to fully explore the parameter space.}
\value{Called for the side effect of plotting.}
\seealso{\code{\link[=mcmcACF,StandardAnnualAggLossDevModelOutput-method]{mcmcACF("StandardAnnualAggLossDevModelOutput")}}
\code{\link[=mcmcACF,BreakAnnualAggLossDevModelOutput-method]{mcmcACF("BreakAnnualAggLossDevModelOutput")}}}
\arguments{\item{object}{The object from which to plot autocorrelations.}}
