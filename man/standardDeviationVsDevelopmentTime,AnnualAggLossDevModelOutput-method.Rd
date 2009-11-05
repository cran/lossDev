\name{standardDeviationVsDevelopmentTime,AnnualAggLossDevModelOutput-method}
\alias{standardDeviationVsDevelopmentTime,AnnualAggLossDevModelOutput-method}
\title{A method to plot and/or return the posterior estimated standard deviation by development year.}
\description{A method to plot and/or return the posterior estimated standard deviation by development year.}
\details{Aggregate loss development models in \pkg{lossDev} allow for changes (by development year) in the measurement error around the log incremental payments.
This is a method that allows for the retrieval and illustration of this standard deviation.}
\value{Mainly called for the side effect of plotting.  Also returns a numeric vector of the plotted statistics.  Returned invisibly.}
\docType{methods}
\seealso{\code{\link{standardDeviationVsDevelopmentTime}}}
\arguments{\item{object}{The object of type \code{AnnualAggLossDevModelOutput} from which to plot and/or return the estimated standard deviation by development year.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
