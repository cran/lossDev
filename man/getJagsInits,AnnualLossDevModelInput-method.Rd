\name{getJagsInits,AnnualLossDevModelInput-method}
\alias{getJagsInits,AnnualLossDevModelInput-method}
\title{A method to collect all the needed initial values common to both the standard model and the break model.}
\description{A method to collect all the needed initial values common to both the standard model and the break model. Intended for internal use only.}
\details{There are currenly two types of \code{AnnualAggLossDevModel}s (break and standard).  These models have many features in common, and code to create initial values for these common features is placed in this method.
The derived classes \code{StandardAnnualAggLossDevModelInput} and \code{BreakAnnualAggLossDevModelInput} call this method via \code{NextMethod()} and then return a new function.}
\value{A named list of the specific model elements.  See details for more information.}
\docType{methods}
\arguments{\item{object}{An object of type \code{AnnualAggLossDevModelInput} from which to collect the needed initial values for the model.}}
