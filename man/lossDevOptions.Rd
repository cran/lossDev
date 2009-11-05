\name{lossDevOptions}
\alias{lossDevOptions}
\title{Options for lossDev.}
\usage{lossDevOptions(...)}
\description{Options for \pkg{lossDev}.}
\details{Currently the only options are \code{keepCodaOnDisk} and \code{logsplinePenaltyFunction}.

\describe{
\item{\code{keepCodaOnDisk}}{
If \code{TRUE} (the default is FALSE), then \pkg{filehash} will be used to store the coda for every node in a temporary file.
This reduces the required memory and can improve copying performance.
Since copied objects refer to the same file, copy time can be greatly reduced.
If \code{FALSE}, then coda's will be kept in memory.
Changing the value will only be taking into account on a "going forward" basis.
}
\item{\code{logsplinePenaltyFunction}}{
When drawing kernal density plots using the \pkg{logspline}, it maybe desirable to specify a penalty to smooth the density (See \code{?logspline}).
This value must be a function which takes one paramter (a vector of the sampled data points) and returns one value -- the penalty.
The default returns the the log of the number of draws.
}

}}
\value{The current (or altered) list of option settings is returned.}
\arguments{\item{\dots}{named values to set.  If empty, only the current list of option settings is returned.}}
