\name{AnnualAggLossDevModelOutput-class}
\alias{AnnualAggLossDevModelOutput-class}
\title{The base output class for all aggregate annual models.}
\description{The base output class for all aggregate annual models.}
\details{\code{AnnualAggLossDevModelOutput} is the base output class for all  aggregate annual model objects.
Derived classes should contain all output from a \acronym{JAGS} run of the input object in the slot \dQuote{input}.
Currenly only the slot \dQuote{input} is allowed to be a non-model node.  All other nodes should be the exact name of some settable node in the model.
This is because \code{getModelOutputNodes} currently looks at the slot names to determine what values to set; only slot \dQuote{input} is known to be a slot other than a settable node.
This class is derived from \code{LossDevModelOutput}}
\docType{class}
\seealso{\code{\linkS4class{LossDevModelOutput}}}
