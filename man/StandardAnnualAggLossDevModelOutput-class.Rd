\name{StandardAnnualAggLossDevModelOutput-class}
\alias{StandardAnnualAggLossDevModelOutput-class}
\title{The final output class for all standard aggregate annual models.}
\description{The final output class for all standard aggregate annual models.}
\details{\code{StandardAnnualAggLossDevModelOutput} is the final output class for all  standard aggregate annual model objects.
Currenly, only the slot \dQuote{input} is allowed to be a non-model node.  All other nodes should be the exact name of some settable node in the model.
This is because \code{getModelOutputNodes} currently looks at the slot names to determine what values to set; only slot \dQuote{input} is known to be a slot other than a settable node.
This class is derived from \code{AnnualAggLossDevModelOutput}.}
\docType{class}
\seealso{\code{\linkS4class{AnnualAggLossDevModelOutput}}}
