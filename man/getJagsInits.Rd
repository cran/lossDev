\name{getJagsInits}
\alias{getJagsInits}
\title{A method to retrieve initial values for a JAGS model.}
\description{A method to retrieve initial values for a \acronym{JAGS} model. Intended for internal use only.}
\details{All classes derived from class \code{LossDevModelInput} should override this method.
This the overriding method should return a parameterless function that, when evaluated, returns a named list of initial values.}
\value{A function. Overriding methods should return a parameterless function that, when evaluated, returns a named list of initial values.}
\docType{genericFunction}
\arguments{\item{object}{The object from which to retrieve initial values.}}
