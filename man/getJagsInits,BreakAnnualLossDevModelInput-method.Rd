\name{getJagsInits,BreakAnnualLossDevModelInput-method}
\alias{getJagsInits,BreakAnnualLossDevModelInput-method}
\title{A method to collect all the needed model initial values unique to the break model.}
\description{A method to collect all the needed model initial values unique to the break model.}
\details{There are currently two types of \code{AnnualAggLossDevModel}s (break and standard). Code needed to create initial values specific the break model is placed in this method.
This method returns a parameterless function which when called first calls the function returned by \code{NextMethod()} and then to the list returned by that function appends the following initial values.

\describe{
\item{\code{R.}}{The initial values for the spline node.  Needed because dspline cannot create initial values.}
}}
\value{A named list of the specific model elements.  See details for more info.}
\docType{methods}
\seealso{\code{\link{getJagsInits}}}
\arguments{\item{object}{An object of type \code{BreakAnnualAggLossDevModelInput} from which to collect the needed model initial values.}}
