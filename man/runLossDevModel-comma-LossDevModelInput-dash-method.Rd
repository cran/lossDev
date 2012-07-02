\name{runLossDevModel,LossDevModelInput-method}
\alias{runLossDevModel,LossDevModelInput-method}
\title{A method to run models in lossDev.}
\description{A method to run models in \pkg{lossDev}.}
\details{This method returns a valid output object or flags an error.
This method is suitable for classes properly derived from class \code{LossDevModelInput} that properly override \dQuote{\code{\link{getJagsData}}} and \dQuote{\code{\link{getJagsInits}}}
and whose output type has a valid \code{getModelOutputNodes} method.

\pkg{lossDev} sets the seed in each chain from a random number generated inside of \code{R}.
So to make a run reproducible, all one must do is set the seed (using \code{set.seed}) in \code{R} prior to the execution of this method.}
\value{An object of class \code{LossDevModelOutput}.}
\docType{methods}
\seealso{\code{\link{runLossDevModel}} \code{\link{set.seed}}}
\arguments{\item{object}{The object of type \code{LossDevModelInput} containing the model to estimate.}
\item{burnIn}{An integer to represent the number of initial \acronym{MCMC} iterations to be discarded. (The adaptive phase (\code{nAddapt}) is not considered part of \code{burnIn}.)}
\item{sampleSize}{An integer to represent the number of \acronym{MCMC} iterations to execute following the \code{burnIn}. (Actual number of iterations kept approximately \code{sampleSize / thin}.)}
\item{thin}{Keep every \code{thin}'th value of \code{sampleSize}.}
\item{nChains}{The number of \acronym{MCMC} chains to run.}
\item{nAddapt}{The length of the adaptive phase for the \acronym{MCMC} algorithm. (Default is \code{trunc(burnIn/4)+1}.)}}
