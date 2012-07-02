\name{slot,NodeOutput,character-method}
\alias{slot,NodeOutput,character-method}
\title{A method to override the behavoir of the function slot.}
\description{A method to override the behavoir of the function \code{slot}.}
\details{In order to enhance the memory management, coda files are optionally stored on the harddrive in temporary files and loaded on an as needed basis.
By overriding this function, we are able to make this seamless.
Overriding the function \code{slot} is a slight abuse and, as such, may in the future be replaced by an accessor function.}
\value{Only if name is exactly \dQuote{value} will the method return the \code{mcarray} containing the coda.  Otherwise, it returns the result of \code{callNextMethod()}.}
\docType{methods}
\seealso{\code{\link{slot}}}
\arguments{\item{object}{The object of type \code{NodeOutput} with the slot to look up.}
\item{name}{A character value giving the name of the slot to look up.}}
