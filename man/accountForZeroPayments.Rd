\name{accountForZeroPayments}
\alias{accountForZeroPayments}
\title{A function to take a triangle estimated without considering zero payments, and account for the possibility of zero payments.}
\usage{accountForZeroPayments(object, burnIn=1000, nAddapt=1000)}
\description{A function to take a triangle estimated without considering zero payments, and account for the possibility of zero payments.}
\details{As incremental payments are modeled on the log scale, zero payments (and negative payments) are treated as missing values.
So, without somehow accounting for zero payments, the estimated payments would be overstated.
Zero payments are accounted for by weighting the predicted payment (given that the payment is greater than zero) with the probability that this payment is zero.
(Negative payments are not (currently) accounted for.)
Currently the trajectory for this probably follows a gompertz curve and is constant across exposure years.
This is currently implemented as a function but may be switched to a method.}
\arguments{\item{object}{The object containing the triangle estimated without accounting for zero payments.}
\item{burnIn}{An integer to represent the number of initial \acronym{MCMC} iterations to be discarded. (The adaptive phase (\code{nAddapt}) is not considered part of \code{burnIn}.)}
\item{nAddapt}{The length of the adaptive phase for the \acronym{MCMC} algorithm. (Default is \code{trunc(burnIn/4)+1}.)}}
