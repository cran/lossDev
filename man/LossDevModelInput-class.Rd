\name{LossDevModelInput-class}
\alias{LossDevModelInput-class}
\title{The base input class for all models in lossDev.}
\description{The base input class for all models in \pkg{lossDev}.}
\details{\code{LossDevModelInput} is the base input class for all model objects.
Derived classes must contain all needed data to construct input and initials for the \acronym{JAGS} model.
These are accessed via the S4 methods \dQuote{\code{\link{getJagsData}}} and \dQuote{\code{\link{getJagsInits}}.}}
\docType{class}
\seealso{\code{\linkS4class{LossDevModelOutput}}}
