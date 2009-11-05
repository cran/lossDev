\name{makeBreakAnnualInput}
\alias{makeBreakAnnualInput}
\title{Create an Object of class BreakAnnualAggLossDevModelInput.}
\usage{makeBreakAnnualInput(
incremental.payments=decumulate(cumulative.payments),
first.year.in.new.regime=trunc(median(as.integer(dimnames(incremental.payments)[[1]]))),
prior.for.first.year.in.new.regime=c(2,2),
extra.dev.years=1,
extra.exp.years=1,
non.stoch.inflation.rate=0,
non.stoch.inflation.weight=1,
stoch.inflation.rate=0,
stoch.inflation.weight=1-non.stoch.inflation.weight,
stoch.inflation.lower.bound=-1,
stoch.inflation.upper.bound=Inf,
known.stoch.inflation.mean=NA,
known.stoch.inflation.persistence=NA,
total.dev.years=extra.dev.years+dim(incremental.payments)[2],
total.exp.years=extra.exp.years+dim(incremental.payments)[1],
cumulative.payments=cumulate(incremental.payments),
exp.year.type=c('ambiguous', 'py', 'ay'),
prior.for.knot.locations.pre.break=NA,
prior.for.knot.locations.post.break=NA,
use.skew.t=FALSE,
bound.for.skewness.parameter=10,
last.column.with.scale.innovation=dim(incremental.payments)[2],
use.ar1.in.calendar.year=FALSE,
use.ar1.in.exposure.growth=TRUE,
projected.rate.of.decay=NA)}
\description{Create an Object of class \code{BreakAnnualAggLossDevModelInput}.}
\details{The loss development models require a lot of input.  Much of the input is directly dependent on the values of other input.
As such, this function facilitates much of the work of setting model parameters and determining which output to collect.

The function creates an object of class \code{BreakAnnualAggLossDevModelInput}.

Many arguments and much functionality is described in \code{\link{makeStandardAnnualInput}}.

\describe{
\item{\code{first.year.in.new.regime}}{
The break model allows for a structural break along the exposure year axis in the consumption path.
The slope of the pre-break consumption path is used to extend the post-break consumption path.
The time when the break occurs can either be specified with 100% certainty by the user, or it can be estimated by the model.
To specify \emph{the} first exposure year in which a new consumption path applies, the user should supply a single value to the parameter \code{first.year.in.new.regime}.
To have the model estimate the point in time when the break occurs, the user should supply a range (min and max) for the possible values to \code{first.year.in.new.regime}.
Note also that the number of rows above and below the break must be at least 4 and as such a triangle of size 8 is the smallest triangle which can be estimated with a break.
}
\item{\code{prior.for.first.year.in.new.regime}}{
The prior for the \code{first.year.in.new.regime} is a (scaled and discretized) beta distribtion.
A value of \code{c(1,1)} would indicate a uniform distribution.
See \code{\link{firstYearInNewRegime}}.
}
\item{\code{prior.for.knot.locations.pre.break} and \code{prior.for.knot.locations.post.break}}{
If these values are \code{NA} (the default), then \code{prior.for.knot.locations.pre.break} will be assigned a value of 2.
And \code{prior.for.knot.locations.post.break} wil be assigned a value of \code{1 + (num.years.in.post.break.period + 0.5 * num.years.in.break.period)/num.years.in.triangle.}.
These values must either both be \code{NA} or both be set to numbers.
}
}}
\references{Kim, Y., and J. McCulloch (2007) \dQuote{The Skew-Student Distribution with Application to U.S. Stock Market Returns and the Equity Premium,} Department of Economics, Ohio State University, October 2007}
\value{An object of class \code{AggModelInput}.  This the model specified by the returned object must then be estimated using the function \code{runLossDevModel}.}
\arguments{\item{incremental.payments}{A square matrix of incremental payments.  Row names should correspond to the exposure year. Only upper-left (including the diagonal) of Triangle may have non-missing values.  Lower-right must be \code{NA}.}
\item{first.year.in.new.regime}{May be one of two types.  1) A single numeric value.  2) A numeric vector of length 2. See Details.}
\item{prior.for.first.year.in.new.regime}{A numeric vector of length 2. See Details.}
\item{extra.exp.years}{A single integer value (\code{total.exp.years} overrides) greater than or equal to 1 (default is 1) specifying the number of additional exposure years (or rows in the triangle) to project.}
\item{extra.dev.years}{A single integer value (\code{total.dev.years} overrides) greater than or equal to 1 (default is 1) specifying the additional number of development years (or columns in the triangle) to project.}
\item{non.stoch.inflation.rate}{May be one of three types (See \emph{Inflation Rate} in Details). 1) A single numeric value. 2) A vector of numerics (of specific length). 3) A matrix of numerics (of specific dim).}
\item{non.stoch.inflation.weight}{May be one of three types (See \emph{Inflation Rate} in Details). 1) A single numeric value. 2) A vector of numerics (of specific length). 3) A matrix of numerics (of specific dim).}
\item{stoch.inflation.rate}{May be one of two types (See \emph{Inflation Rate} in Details). 1) A single numeric value of \emph{zero}. 2) A vector of numerics (of specific length).}
\item{stoch.inflation.weight}{May be one of three types (See \emph{Inflation Rate} in Details). 1) A single numerical value. 2) A vector of numerics (of specific length). 3) A matrix of numerics (of specific dim).}
\item{stoch.inflation.lower.bound}{May be one of three types (See \emph{Inflation Rate} in Details). 1) A single numeric value. 2) A vector of numerics (of specific length). 3) A matrix of numerics (of specific dim).}
\item{stoch.inflation.upper.bound}{May be one of three types (See \emph{Inflation Rate} in Details). 1) A single numeric value. 2) A vector of numerics (of specific length). 3) A matrix of numerics (of specific dim).}
\item{known.stoch.inflation.mean}{May be one of two types (See \emph{Inflation Rate} in Details). 1) A single numeric value. 2) \code{NA}.}
\item{known.stoch.inflation.persistence}{May be one of two types (See \emph{Inflation Rate} in Details). 1) A single numeric value. 2) \code{NA}.}
\item{total.exp.years}{A single integer value (overrides \code{extra.exp.years}) specifying the last exposure year to project.  Must be at least the number of rows in \code{incremental.payments} + 1.}
\item{total.dev.years}{A single integer value (overrides \code{extra.dev.years}) specifying the last development year to project.  Must be at least the number of columns in \code{incremental.payments} + 1 .}
\item{cumulative.payments}{A numeric matrix with the same dim and dimnames as \code{incremental.payments}.  Must be a possible cumulative payment triangle of \code{incremental.payments}.  (See \emph{Cumulative Payments} Section).}
\item{exp.year.type}{A single character value indicating the type of exposure years:  \sQuote{ambiguous}, \sQuote{py}, and \sQuote{ay} mean \sQuote{Exposure Year}, \sQuote{Policy Year}, and \sQuote{Accident Year}; respectively.}
\item{prior.for.knot.locations.pre.break}{A single numeric value of at least 1.  The prior for the location of knots is a scaled beta with parameters \code{c(1,prior.for.knot.locations.pre.break)}.}
\item{prior.for.knot.locations.post.break}{See \code{prior.for.knot.locations.pre.break}. Large values produce stable consumption paths at high development years.}
\item{use.skew.t}{A logical value.  If \code{TRUE} the model assumes the observed and estimated log incremental payments are realizations from a skewed \eqn{t} distribution; if \code{FALSE} it will assume zero skewness. (See Reference.)}
\item{bound.for.skewness.parameter}{A positive numerical value representing the symetric boundaries for the skewness parameter.  In most cases, the default should be sufficient. Ignored if \code{use.skew.t=FALSE}.}
\item{last.column.with.scale.innovation}{A single integer must be at least 1 and at most the number of columns in \code{incremental.payments}.  See \emph{Measurment Error-Second Order Random Walk} in Details.}
\item{use.ar1.in.calendar.year}{A logical value.  The calendar year effect errors may (at the users digression) include an autoregressive process of order 1.  \code{TRUE} turns on the ar1 process, \code{FALSE} turns it off.}
\item{use.ar1.in.exposure.growth}{A logical value.  The exposure growth errors may (at the users discretion) include an autoregressive process of order 1.  \code{TRUE} (the Default) turns on the ar1 process, \code{FALSE} turns it off.}
\item{projected.rate.of.decay}{May be one of three types (See \emph{Projected Rate of Decay} in Details). 1) \code{NA}. 2) A matrix of numerics (of specific dim). 3) A named list.}}
