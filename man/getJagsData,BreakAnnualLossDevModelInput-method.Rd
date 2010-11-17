\name{getJagsData,BreakAnnualLossDevModelInput-method}
\alias{getJagsData,BreakAnnualLossDevModelInput-method}
\title{A method to collect all the needed model input specific to the break model.}
\description{A method to collect all the needed model input specific to the break model.}
\details{There are currently two types of \code{AnnualAggLossDevModel}s (break and standard).  These models have many data elements in common.
This method appends only the elements specific to the break model onto the list created by a call to \code{NextMethod()}.

The following elements are appended:
\describe{
\item{\code{x.0}}{Two values. The lower bound for the number of knots. First is for first spline.}
\item{\code{x.r}}{Two values. The upper bound for the number of knots. First is for first spline.}
\item{\code{break.row}}{Two integer values given the (inclusive) range of rows in which the first year in the new regime could occur.}
\item{\code{break.row.priors}}{The parameters for the beta distribution which serves as the prior for the location of the structural break.}
\item{\code{K.trim}}{The maximum number of columns in the post-break triangle.}
\item{\code{beta.prior}}{A matrix giving the prior for the location of knots.  First column is for the pre-break spline.  Second is for the post-break spline.}
\item{\code{mu.number.of.knots.prior}}{A matrix giving the prior for the mean of the number of knots.  First column is for the pre-break spline.  Second is for the post-break spline}
\item{\code{number.of.knots.ubound}}{A vector giving the upper bound for the number of knots.  First is for the pre-break spline.  Second is for the post-break spline}
}}
\value{A named list of the specific model elements.  See details for more info.}
\docType{methods}
\arguments{\item{object}{An object of type \code{BreakAnnualAggLossDevModelInput} from which to collect the needed model input.}}
