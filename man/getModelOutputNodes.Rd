\name{getModelOutputNodes}
\alias{getModelOutputNodes}
\title{A method to determine which JAGS nodes the output object is expecting.}
\description{A method to determine which \acronym{JAGS} nodes the output object is expecting.  Intended for internal use only.}
\details{This method is used to query the nodes a particular \code{LossDevModelOutput} is expecting.}
\value{A character vector with each element corresponding to the name of an expected node.}
\docType{genericFunction}
\seealso{\code{\link[=getModelOutputNodes,LossDevModelOutput-method]{getModelOutputNodes("LossDevModelOutput")}}}
\arguments{\item{object}{The object from which to get the expected nodes.}}
