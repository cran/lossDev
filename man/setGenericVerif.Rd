\name{setGenericVerif}
\alias{setGenericVerif}
\title{A Safe Version of setGeneric.}
\usage{setGenericVerif(name, ...)}
\description{A Safe Version of \code{setGeneric}. Intended for internal use only.}
\details{\code{setGeneric} will overwrite existing generic functions.  This will result in the loss of all methods already associated with that generic.
\code{setGenericVerif} only sets the generic if it is not already a generic.
If a generic by the name of \code{name} already exists, a warning is issued and NULL is returned.  Otherwise \code{setGeneric} is called and its value returned.}
\value{\code{setGenericVerif} really exists for its side effect; but returns the value returned by \code{setGeneric} or NULL.}
\seealso{\code{\link{setGeneric}}}
\arguments{\item{name}{The character string name of the generic function.}
\item{\dots}{Additional arguments to pass to \code{setGeneric}.}}
