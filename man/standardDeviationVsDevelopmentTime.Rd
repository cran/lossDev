\name{standardDeviationVsDevelopmentTime}
\alias{standardDeviationVsDevelopmentTime}
\title{A generic function to plot and/or return the posterior estimated standard deviation by development year.}
\description{A generic function to plot and/or return the posterior estimated standard deviation by development year.}
\details{Aggregate loss development models in \pkg{lossDev} allow for changes (by development year) in the measurement error around the log incremental payments.
This is a generic function that allows for the retrieval and illustration of this standard deviation.}
\value{Mainly called for the side effect of plotting.}
\docType{genericFunction}
\seealso{\code{\link[=standardDeviationVsDevelopmentTime,AnnualAggLossDevModelOutput-method]{standardDeviationVsDevelopmentTime("AnnualAggLossDevModelOutput")}}}
\arguments{\item{object}{The object from which to plot and/or return the estimated standard deviation by development year.}
\item{plot}{A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.}}
