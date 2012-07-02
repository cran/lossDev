\name{getJagsInits,StandardAnnualLossDevModelInput-method}
\alias{getJagsInits,StandardAnnualLossDevModelInput-method}
\title{A method to collect all the needed initial values unique to the standard model.}
\description{A method to collect all the needed initial values unique to the standard model. Intended for internal use only.}
\details{There are currently two types of \code{AnnualAggLossDevModel}s (break and standard). Code needed to create initial values specific the standard model is placed in this method.
This method returns a parameterless function which, when called, first calls the function returned by \code{NextMethod()} and then appends the following initial values onto the list returned by that function:

\describe{
\item{\code{S.}}{The initial values for the spline node.  Needed because \emph{dspline} cannot create initial values.}
}}
\value{A named list of the specific model elements.  See details for more information.}
\docType{methods}
\seealso{\code{\link{getJagsInits}}}
\arguments{\item{object}{An object of type \code{StandardAnnualAggLossDevModelInput} from which to collect the needed initial values for the model.}}
