\name{getJagsData}
\alias{getJagsData}
\title{A method to retrieve data for a JAGS model.}
\description{A method to retrieve data for a \acronym{JAGS} model. Intended for internal use only.}
\details{Any classes derived from class \code{LossDevModelInput} must provide a method to override this generic function.
This overriding method should return a named list of data values.}
\value{This method returns a named list of data values.}
\arguments{\item{object}{The object from which to retrieve initial values.}}
