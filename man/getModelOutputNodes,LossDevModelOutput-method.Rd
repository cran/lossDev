\name{getModelOutputNodes,LossDevModelOutput-method}
\alias{getModelOutputNodes,LossDevModelOutput-method}
\title{A method to determine which JAGS nodes the output object is expecting.}
\description{A method to determine which \acronym{JAGS} nodes the output object is expecting.  Intended for internal use only.}
\details{This method works by examining all the slot names in \code{object} and returning all of them except for a slot named \dQuote{input.}}
\value{A character vector with each element corresponding to the name of an expected node.}
\docType{methods}
\seealso{\code{\link{getModelOutputNodes}}}
\arguments{\item{object}{The object of type \code{LossDevModelOutput} from which to get the expected nodes.}}
