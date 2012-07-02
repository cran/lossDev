\name{LossDevModelOutput-class}
\alias{LossDevModelOutput-class}
\title{The base output class for all models in lossDev.}
\description{The base output class for all models in \pkg{lossDev}.}
\details{\code{LossDevModelOutput} is the base class for all output model objects.
Derived classes should contain all needed output from a \acronym{JAGS} run and the input object associated with that run in the slot \dQuote{input.}
Currenly, only the slot \dQuote{input} is allowed to be a non-model node.  All other nodes must be the exact name of some settable node in the model.
This is because \code{getModelOutputNodes} currently looks at the slot names to determine what values to set; only slot \dQuote{input} is known to be a slot other than a settable node.
Any class adding extra non-model node slots must also override the method \code{getModelOutputNodes}.}
\docType{class}
\seealso{\code{\linkS4class{AnnualAggLossDevModelInput}}}
