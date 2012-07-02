\name{getNonMissingIndexes}
\alias{getNonMissingIndexes}
\title{A function to return the index non-missing values of a 2d container.}
\usage{getNonMissingIndexes(m)}
\description{A function to return the index non-missing values of a 2d container.}
\value{A 2d array of integers.  Each rows corresponds to a non-missing value in \code{m}.  First column is the row.  Second column is the column.}
\keyword{internal}
\arguments{\item{m}{A 2d container on which \code{dim} can be called and is subsetable.}}
