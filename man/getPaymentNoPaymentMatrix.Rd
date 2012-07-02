\name{getPaymentNoPaymentMatrix}
\alias{getPaymentNoPaymentMatrix}
\title{A function to turn a matrix of incremental payments into zero or ones depending upon whether a payment is positive.}
\usage{getPaymentNoPaymentMatrix(object)}
\description{A function to turn a matrix of incremental payments into zero or ones depending upon whether a payment is positive. Intended for internal use only.}
\details{The conversion rule is as follows.  If \code{NA}, then \code{NA}. Else if greater than zero, then 1.  Else if equal to zero, then zero. Else \code{NA}.}
\value{A matrix of zero or one (or \code{NA}) matching the structure of in input matrix.}
\arguments{\item{object}{The matrix of incremental payments.}}
