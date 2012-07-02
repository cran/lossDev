\name{lossDevOptions}
\alias{lossDevOptions}
\title{Options for lossDev.}
\usage{lossDevOptions(...)}
\description{Options for \pkg{lossDev}.}
\details{Currently the only options are \code{keepCodaOnDisk} and \code{logsplinePenaltyFunction}.

\describe{
\item{\code{logsplinePenaltyFunction}}{
When drawing kernal density plots using the \pkg{logspline}, it maybe desirable to specify a penalty to smooth the density (See \code{?logspline}).
This value must be a function which takes one paramter (a vector of the sampled data points) and returns one value -- the penalty.
The default returns the the log of the number of draws.
}

}}
\value{The current (or altered) list of option settings is returned.}
\arguments{\item{\dots}{named values to set.  If empty, only the current list of option settings is returned.}}
