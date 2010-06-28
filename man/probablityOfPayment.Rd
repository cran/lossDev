\name{probablityOfPayment}
\alias{probablityOfPayment}
\title{A generic function to plot the probability of a payment.}
\description{A generic function to plot the probability of a payment.}
\details{Because the model is Bayesian, each estimated payment comes as a distribution.
The median of this distribution is used as a point estimate when plotting and/or returning values.
Note: Negative payments are treated as missing and are not accounted for.}
\value{Mainly called for the side effect of plotting.  Also returns a matrix containing the (median) probably of payment.  Returned invisibly.}
\seealso{\code{\link{accountForZeroPayments}}}
\arguments{\item{object}{The object from which to plot the probability of a payment.}
\item{plot}{A logical value. If \code{TRUE}, then the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
