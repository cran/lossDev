\name{estimate.priors}
\alias{estimate.priors}
\title{A function to estimate priors for the gompertz curve.}
\usage{estimate.priors(p)}
\description{A function to estimate priors for the gompertz curve. Intended for internal use only.}
\details{The function uses \code{nlm} to minimize the squared error.}
\value{A vector equal in length to the number of columns in x representing the empirical probably of payment.}
\arguments{\item{p}{A vector of the form returned by \code{\link{calculateProbOfPayment}}. \code{NA}s are allowed.}}
