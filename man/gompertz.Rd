\name{gompertz}
\alias{gompertz}
\title{The gompertz function.}
\usage{gompertz(x, scale, fifty.fifty)}
\description{The gompertz function. Intended for internal use only.}
\details{This function is used as the probably of observing a zero payment.
(Or one minus the probably of observing a positive payment.)
(Note that negative payments are assumed to be missing values.)}
\value{An object of type \code{AnnualAggLossDevModelOutputWithZeros} and \code{AnnualAggLossDevModelOutput}.}
\arguments{\item{x}{The value(s) at which to evaluate the gompertz function.}
\item{scale}{The scale parameter should always (or at least for how it is used in lossDev) be positive as this indicates an increasing probably of a zero payment.}
\item{fifty.fifty}{The value at which the gompertz function returns 0.5.}}
