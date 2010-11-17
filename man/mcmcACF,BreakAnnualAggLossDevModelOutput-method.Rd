\name{mcmcACF,BreakAnnualAggLossDevModelOutput-method}
\alias{mcmcACF,BreakAnnualAggLossDevModelOutput-method}
\title{A method to plot autocorrelations found in the MCMC samples for select parameters.}
\description{A method to plot autocorrelations found in the \acronym{MCMC} samples for select parameters.}
\details{Chains with high autocorrelation require a longer burnin and more samples to fully explore the parameter space.}
\value{Called for the side effect of plotting.}
\docType{methods}
\seealso{\code{\link{mcmcACF}}
\code{\link[=mcmcACF,StandardAnnualAggLossDevModelOutput-method]{mcmcACF("StandardAnnualAggLossDevModelOutput")}}}
\section{TODO}{Add option to plot other values.  Currently only plots \dQuote{First Rate Of Decay} (for the pre and post break), \dQuote{First Calendar Year Effect Error,} and \dQuote{First Exposure Year Growth,}}
\arguments{\item{object}{The object from which to plot autocorrelations.}}
