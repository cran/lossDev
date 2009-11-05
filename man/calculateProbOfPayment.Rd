\name{calculateProbOfPayment}
\alias{calculateProbOfPayment}
\title{A function to calculate an empirical vector of the probability of payment.}
\usage{calculateProbOfPayment(x)}
\description{A function to calculate an empirical vector of the probability of payment. Intended for internal use only.}
\value{A vector equal in length to the number of columns in x representing the empirical probably of payment.}
\arguments{\item{x}{The matrix of the form returned by \code{\link{getPaymentNoPaymentMatrix}}.}}
