\name{makeStandardAnnualInput}
\alias{makeStandardAnnualInput}
\title{Create an Object of class StandardAnnualAggLossDevModelInput.}
\usage{makeStandardAnnualInput(
incremental.payments=decumulate(cumulative.payments),
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
prior.for.knot.locations=2,
use.skew.t=FALSE,
bound.for.skewness.parameter=10,
last.column.with.scale.innovation=dim(incremental.payments)[2],
use.ar1.in.calendar.year=FALSE,
use.ar1.in.exposure.growth=TRUE,
projected.rate.of.decay=NA)}
\description{Create an Object of class \code{StandardAnnualAggLossDevModelInput}.}
\details{The loss development models require extensive input.  Much of the input is directly dependent on the values of other input.
As such, this function facilitates much of the work of setting model parameters and determining which output to collect.

The function creates an object of class \code{StandardAnnualAggLossDevModelInput}.

\describe{
\item{\bold{Inflation Rate}}{
Workers Compensation indemnity statutes can be extremely complex.
In order to provide a sufficient level of flexibility, the model allows for a combination for 3 inflation/escalation rates: one non-stochastic rate, one stochastic rate, and one zero rate.
All inflation rates are on the real scale, that is, a 5 percent increase is entered as \dQuote{0.05.}
Inflation rates are only used to create an expected inflation rate.  Actual calendar year effects deviate from this prior along the diagional.
\describe{
\item{Non-Stochastic Rate of Inflation/Escalation}{
The user must specify a non-stochastic rate of inflation for all observed and future time periods.
Assuming a non-stochastic inflation rate of zero turns this feature off.
\describe{
\item{\code{non.stoch.inflation.rate}}{
The non-stochastic inflation rate is specified by the parameter \code{non.stoch.inflation.rate}.
The default value is set to \dQuote{0.}
There are three possible ways to specify \code{non.stoch.inflation.rate}:
\describe{
\item{If specified as a single value}{
The inflation in every cell is assumed to be this value.
For instance, choosing zero results in no non-stochastic inflation rate.
}
\item{If specified as a vector}{
Its length must be two less than \code{total.exp.years + total.dev.years} and must have \code{names} starting at \code{dimnames(incremental.payments) [[1]][2]} and increasing by one.
The first value corresponds to the inflation rate from the first to second diagonal (or cell 1,1 to cells 1,2 and 2,1).
The second value corresponds to the inflation rate from the second to third diagonal. etc.
}
\item{If specified as a matrix}{
It must have number of rows equal to \code{total.exp.years} and number of columns equal to \code{total.dev.years}.
The row names must be equal to the row names of \code{incremental.payments}.
Each cell represents the inflation to that cell.
Cells in the first column represent inflation down rows.
Cells in other columns represent inflation from left to right.
Cell 1,1 is ignored.
}
}
}
\item{\code{non.stoch.inflation.weight}}{
The user has the option of allowing the percentage of dollars that inflate/escalate at the non-stochastic rate to vary.
This is done by properly setting the value of \code{non.stoch.inflation.weight}.
The default value is set to \dQuote{1.}
\code{non.stoch.inflation.weight} must be at least zero and at most one.  Also \code{non.stoch.inflation.weight + stoch.inflation.weight} must be at most one.
There are three possible ways to specify \code{non.stoch.inflation.weight}:
\describe{
\item{If specified as a single value}{
The proportion of dollars that inflate at the rate of \code{non.stoch.inflatoin.rate} in every cell is assumed to be this value.
For instance, choosing zero results in the model assuming that no non-stochastic inflation applies.
}
\item{If specified as a vector}{
Its length must be \code{total.dev.years}.
This vector is repeated to produce a matrix such that each row equals this vector.
This matrix is then assumed to be the supplied for this parameter (See below).
See next item.
}
\item{If specified as a matrix}{
It must have number of rows equal to \code{total.exp.years} and number of columns equal to \code{total.dev.years}.
The row names must be equal to the row names of \code{incremental.payments}.
Each cell represents the proportion of dollars in the previous cell inflating at the \code{non.stoch.inflation.rate} rate to that cell.
Cells in the first column represent inflation down rows.
Cells in other columns represent inflation from left to right.
Cell 1,1 is ignored.
}
}
}
}
}
\item{Stochastic Rate of Inflation/Escalation}{
The stochastic inflation rate is modeled as an AR1 process (on the log scale).
The user must supply a stochastic rate of inflation for at least all observed diagonals and has the option of supplying values for diagonals beyond the final diagonal in the observed triangle.
As such, the stochasticity of the stochastic rate of inflation only applies to future rates of inflation.
Values not supplied for future diagonals are predicted using the estimated parameters for the supplied inflation series.
If the user does not wish to supply/assume a stochastic rate of inflation, he may set the stochastic inflation rate to zero for all periods
(\code{stoch.inflation.rate=0}).
\describe{
\item{\code{stoch.inflation.rate}}{
The stochastic inflation rate is specified by the parameter \code{stoch.inflation.rate}.
The default value is set to \dQuote{0.}
There are two possible ways to specify \code{stoch.inflation.rate}:
\describe{
\item{If specified as a single 0 (zero)}{
It is assumed that no stochastic inflation applies.
}
\item{If specified as a vector}{
Its length must be at least one less than the number of rows in \code{incremental.payments}.
Future values are simulated as needed.
(No element in \code{stoch.inflation.rate} may be \code{NA}; only supply values for known rates of inflation.)
The vector must be named and the names must be consecutive integers.
The names must contain the second and last values of \code{dimnames(incremental.payments)[[1]]}.
The names represent the calendar year of inflation.
}
}
}
\item{\code{stoch.inflation.weight}}{
The default value is set to \dQuote{\code{1 - non.stoch.inflation.weight}}.
See \code{non.stoch.inflation.weight}.
}
\item{\code{stoch.inflation.lower.bound}}{
The user has the option of putting a lower bound on the stochastic rate of inflation/escalation.
\code{stoch.inflation.rate} should not account for \code{stoch.inflation.lower.bound} as this would result in incorrectly simulated future inflation rates.
The model properly applies the \code{stoch.inflation.lower.bound} to both observed and future rates of inflation prior to \dQuote{inflating} dollars in the triangle.
The user accounts for lower bounds by properly setting the value of \code{stoch.inflation.lower.bound}.
\code{stoch.inflation.lower.bound} should be on the real scale (\code{x[i] / x[i-1] - 1}).
The default value is set to \code{-1}, which is no bound.
Bounds are applied prior to weights.
There are three possible ways to specify \code{stoch.inflation.weight}:
\describe{
\item{If specified as a single value}{
All stochastic rates of inflation in the triangle are bound by this value.
}
\item{If specified as a vector}{
Its length must be \code{total.dev.years}.
This vector is repeated to produce a matrix with each rows equal to it.
It will then be as if this matrix was supplied for this parameter.
See next item.
}
\item{If specified as a matrix}{
It must have number of rows equal to \code{total.exp.years} and number of columns equal to \code{total.dev.years}.
The row names should be equal to the row names of \code{incremental.payments}.
Each cell represents the bound on the stochastic rate of inflation from the previous cell to that cell.
Cells in the first column represent inflation down rows.
Cells in other columns represent inflation from left to right.
Cell 1,1 is ignored.
}
}
}
\item{\code{stoch.inflation.upper.bound}}{
Default value is \code{Inf}, which is no bound.
See \code{stoch.inflation.lower.bound}.
}
\item{\code{known.stoch.inflation.mean}}{
The AR1 process used to simulate future stochastic rates of inflation has a non-zero mean which is (at least approximately) equal to the historical mean.
If for some reason (i.e., the series of inflation rates is too short) the user believes this historical mean poorly represents future rates of inflation, the user can override the estimation process.
If \code{known.stoch.inflation.mean} is set to \code{NA}, the mean is estimated.  Otherwise the mean is assumed (with certainty) to be the specified value.
(Note that since the estimation takes place on the log scale, the logarithm of \code{known.stoch.inflation.mean} plus 1 is used as the mean on the log scale;
this results in \code{known.stoch.inflation.mean} being the geometric mean.
}
\item{\code{known.stoch.inflation.persistence}}{
The user has the option of overriding the estimation of the \eqn{rho} coefficient in the AR1 process used to simulate future stochastic rates of inflation.
If \code{known.stoch.inflation.persistence} is set to \code{NA}, the \eqn{rho} coefficient is estimated; otherwise it is taken as this value.
}
}
}
}
}
\item{\bold{Projected Rate of Decay}}{
The model estimates an inflation-free rate of decay up to the size of the triangle (\code{dim(incremental.payments)[2]}).
Since no data is available to the model beyond this point, the model must be supplied with an assumption for future values.
This assumption is supplied via the parameter \code{projected.rate.of.decay} in one of two ways:
\describe{
\item{A single value of \code{NA}}{
The inflation-free rate of decay for the final development year in the triangle (\code{dim(incremental.payments)[2]}) is used for all future development years through development year \code{total.dev.years}.
}
\item{A named list}{
This is currently a rather low level interface allowing for much flexibility at the expense of some ease of usability.
(In the (possibly near) future, it is likely that the burden placed on the user will be eased.)
The named list must contain three elements (with an optional fourth element):
\describe{
\item{\code{last.non.zero.payment}}{
The named element \code{last.non.zero.payment} refers to the last development year in which a non-zero payment may occur.
It is assumed that all subsequent incremental payments are exactly zero.
Must be at most \code{total.dev.years} and at least \code{dim(incrementals)[2]+2}.
\code{last.non.zero.payment} must be specified as a named numeric vector of length \code{total.exp.years}.
Element one then refers to the first exposure year; and element two to the second.
The names must start with \code{dimnames(incrementals)[[1]][1]} and increment by 1.
\code{last.non.zero.payment} can be omitted.  If omitted, then the \code{last.non.zero.payment} for each exposure year is \code{total.dev.years}.
}
\item{\code{final.rate.of.decay.mean}}{
The parameter \code{final.rate.of.decay.mean} refers to the rate of decay in \code{last.non.zero.payment}.
As the rate of decay is estimated on the log scale, \code{final.rate.of.decay.mean} is first increased by 1 and then the log of the resulting value is taken; the result is then treated as the mean on the log scale.
This results in \code{final.rate.of.decay.mean} being the geometric mean.
(\code{final.rate.of.decay.mean=0.05} means that the mean of the final rate of decay (on the log scale) is \code{log(1.05)}.)
This may be supplied as a single numeric value, in which case all exposure years are assumed to have the same \code{final.rate.of.decay.mean}.
Or this may be supplied as a named numeric vector of length \code{total.exp.years}.
Element one then refers to the first exposure year, and element two refers to the second, etc.
The names of the vector must start with \code{dimnames(incrementals)[[1]][1]} and increment by \code{1}.
(A Technical Note: \code{last.non.zero.payment} actually defines the final rate of decay since the rate of decay from it to the next development year is \code{-1.00}.
Thus a better name for \code{final.rate.of.decay.mean} does not refer to the \emph{final} rate of decay but rather the \emph{penultimate} rate of decay.)
}
\item{\code{final.rate.of.decay.sd}}{
The parameter \code{final.rate.of.decay.sd} refers to the standard deviation around \code{final.rate.of.decay.mean}.
This may be supplied as a single numeric value, in which case all exposure years are assumed to have the same \code{final.rate.of.decay.sd}.
Or this may be supplied as a named numeric vector of length \code{total.exp.years}.
Element one then refers to the first exposure year, and element two refers to the second, etc.
The names of the vector must start with \code{dimnames(incrementals)[[1]][1]} and increment by \code{1}.
}
\item{\code{final.rate.of.decay.weight}}{
The basic concept behind a user-supplied rate of decay is that of a weighted average between the final estimated rate of decay (the one associated with \code{dim(incrementals)[2]}) and
the final supplied value \code{final.rate.of.decay.mean}.
Thus the user must supply a vector of weights to describe the path from the estimated value to the projected value.
The \dQuote{weight} is the weight (on the log scale) of the user-supplied \code{final.rate.of.decay}.
So, a value of \code{1} means that the projected rate of decay should take on \code{final.rate.of.decay} and a value of \code{0} means that the projected rate of decay should equal the final estimated rate of decay.
Note that while these weights are allowed to be outside the interval \code{[0,1]}, convention should dictate they remain within this interval.
These weights are supplied by specifying a value for \code{final.rate.of.decay.weight}, which may be specified in two ways.
It may be specified as a numeric matrix with number of columns equal to \code{total.dev.years - dim(incremental.payments)[2]} and number of rows equal to \code{total.exp.years},
in which case the first row corresponds to the first exposure year in \code{incremental.payments} and the first column corresponds to the development year immediately following the last one in \code{incremental.payments}.
It may be specified as a numeric vector of length \code{total.dev.years - dim(incremental.payments)[2]}, in which case it is assumed to be the same for all exposure years.
}
}
}
}
}
\item{\bold{Cumulative Triangle}}{
As a safeguard, the cumulative and incremental triangles must match up.  They are tested by taking the cumulative triangle and then decumulating it with the function \code{decumulate}.
The decumulated triangle is then compared to the incremental triangle.  Except for \code{NA}'s in either triangle, the values in the triangles must match.  (Errors due to numerical precision are accounted for.)
}
\item{\bold{Measurment Error Second-Order Random Walk}}{
The model allows for changes in the scale of the measurement error with development time.
Changes in the scale parameter (see \code{\link{scaleParameter}}) are assumed to follow a second-order random walk on the log scale.
The model allows for imposed stationarity in the scale parameter by setting of the value \code{last.column.with.scale.innovation}.
\code{last.column.with.scale.innovation} must be coercible to an integer, be least 1, and at most the number of columns in \code{incremental.payments}.
A value of 1 constrains the scale parameter in each column to the same value.
A value of \code{dim(incremental.payments)[2]} allows for all columns to have their own scale parameter (but smoothed by the second-order random walk).
The effective value used by the model for this argument is the minimum of the supplied value and the last column in \code{incremental.payments} with a non-missing value (on the log scale).
Note that since the scale parameter is assumed to follow a second-order random walk, a value of \code{2} results in an (effectively) unconstrained estimation for the scale parameter in the first column.
}
}}
\references{Kim, Y., and J. McCulloch (2007) \dQuote{The Skew-Student Distribution with Application to U.S. Stock Market Returns and the Equity Premium,} Department of Economics, Ohio State University, October 2007}
\value{An object of class \code{AggModelInput}.  The model specified by the returned object must then be estimated using the function \code{runLossDevModel}.}
\arguments{\item{incremental.payments}{A square matrix of incremental payments.  Row names must correspond to the exposure year. Only the upper-left (including the diagonal) of this matrix may have non-missing values.  Lower-right must be \code{NA}.}
\item{extra.exp.years}{A single integer value (\code{total.exp.years} overrides) greater than or equal to 1 (default is 1) specifying the number of additional exposure years (or rows in the triangle) to project forward.}
\item{extra.dev.years}{A single integer value (\code{total.dev.years} overrides) greater than or equal to 1 (default is 1) specifying the additional number of development years (or columns in the triangle) to project forward.}
\item{non.stoch.inflation.rate}{May be one of three types (See \emph{Inflation Rate} in Details): 1) A single numeric value; 2) a vector of numerics (of specific length); 3) a matrix of numerics (of specific dim).}
\item{non.stoch.inflation.weight}{May be one of three types (See \emph{Inflation Rate} in Details): 1) A single numeric value; 2) a vector of numerics (of specific length); 3) a matrix of numerics (of specific dim).}
\item{stoch.inflation.rate}{May be one of two types (See \emph{Inflation Rate} in Details): 1) A single numeric value of \emph{zero}; 2) a vector of numerics (of specific length).}
\item{stoch.inflation.weight}{May be one of three types (See \emph{Inflation Rate} in Details): 1) A single numerical value; 2) a vector of numerics (of specific length); 3) a matrix of numerics (of specific dim).}
\item{stoch.inflation.lower.bound}{May be one of three types (See \emph{Inflation Rate} in Details): 1) A single numeric value; 2) a vector of numerics (of specific length); 3) a matrix of numerics (of specific dim).}
\item{stoch.inflation.upper.bound}{May be one of three types (See \emph{Inflation Rate} in Details): 1) A single numeric value; 2) a vector of numerics (of specific length); 3) a matrix of numerics (of specific dim).}
\item{known.stoch.inflation.mean}{May be one of two types (See \emph{Inflation Rate} in Details): 1) A single numeric value; 2) \code{NA}.}
\item{known.stoch.inflation.persistence}{May be one of two types (See \emph{Inflation Rate} in Details): 1) A single numeric value; 2) \code{NA}.}
\item{total.exp.years}{A single integer value (overrides \code{extra.exp.years}) specifying the last exposure year to project forward.  Must be at least the number of rows in \code{incremental.payments} + 1.}
\item{total.dev.years}{A single integer value (overrides \code{extra.dev.years}) specifying the last development year to project forward.  Must be at least the number of columns in \code{incremental.payments} + 1 .}
\item{cumulative.payments}{A numeric matrix with the same dim and dimnames as \code{incremental.payments}.  Must be a possible cumulative payment triangle of \code{incremental.payments}.  (See \emph{Cumulative Payments} Section.)}
\item{exp.year.type}{A single character value indicating the type of exposure years:  \sQuote{ambiguous}, \sQuote{py}, and \sQuote{ay} mean \sQuote{Exposure Year}, \sQuote{Policy Year}, and \sQuote{Accident Year}; respectively.}
\item{prior.for.knot.locations}{A single numeric value of at least 1.  The prior for the location of knots is a scaled beta with parameters \code{c(1,prior.for.knot.locations)}.  Large values produce stable consumption paths at high development years.}
\item{use.skew.t}{A logical value.  If \code{TRUE}, the model assumes that the observed and estimated log incremental payments are realizations from a skewed \eqn{t} distribution; if \code{FALSE} it assumes zero skewness. (See Reference.)}
\item{bound.for.skewness.parameter}{A positive numerical value representing the symetric boundaries for the skewness parameter.  In most cases, the default should be sufficient. Ignored if \code{use.skew.t=FALSE}.}
\item{last.column.with.scale.innovation}{A single integer. Must be at least 1 and at most the number of columns in \code{incremental.payments}.  See \emph{Measurment Error-Second Order Random Walk} in Details.}
\item{use.ar1.in.calendar.year}{A logical value.  The calendar year effect errors may (at the users discretion) include an autoregressive process of order 1.  \code{TRUE} turns on the ar1 process, \code{FALSE} (the Default) turns it off.}
\item{use.ar1.in.exposure.growth}{A logical value.  The exposure growth errors may (at the users discretion) include an autoregressive process of order 1.  \code{TRUE} (the Default) turns on the ar1 process, \code{FALSE} turns it off.}
\item{projected.rate.of.decay}{May be one of three types (See \emph{Projected Rate of Decay} in Details): 1) \code{NA}; 2) a matrix of numerics (of specific dim); 3) a named list.}}
