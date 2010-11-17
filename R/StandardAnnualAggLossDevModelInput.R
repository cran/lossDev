##################################################################################################
##                                                                                              ##
##    lossDev is an R-package.                                                                  ##
##    It is a Bayesian time series model of loss development.                                   ##
##    Features include skewed Student-t distribution with time-varying scale parameters,        ##
##    an expert prior for the calendar year effect,                                             ##
##    and accommodation for structural breaks in the consumption path of services.              ##
##                                                                                              ##
##    Copyright © 2009, National Council On Compensation Insurance Inc.,                        ##
##                                                                                              ##
##    This file is part of lossDev.                                                             ##
##                                                                                              ##
##    lossDev is free software: you can redistribute it and/or modify                           ##
##    it under the terms of the GNU General Public License as published by                      ##
##    the Free Software Foundation, either version 3 of the License, or                         ##
##    (at your option) any later version.                                                       ##
##                                                                                              ##
##    This program is distributed in the hope that it will be useful,                           ##
##    but WITHOUT ANY WARRANTY; without even the implied warranty of                            ##
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                             ##
##    GNU General Public License for more details.                                              ##
##                                                                                              ##
##    You should have received a copy of the GNU General Public License                         ##
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.                     ##
##                                                                                              ##
##################################################################################################


##This is intended to be a final class for models without a break.
##This class is derived from the AnnualAggLosDevModelInput class.
##' @include zzz.R
##' @include AnnualAggLossDevModelInput.R
NULL

##' The final input class for models without a break.
##'
##' \code{StandardAnnualAggLossDevModelInput} is the final input class for models without a break.
##' @name StandardAnnualAggLossDevModelInput-class
##' @docType class
##' @seealso \code{\linkS4class{LossDevModelInput}}
##'  \code{\linkS4class{StandardAnnualAggLossDevModelInput}}
setClass(
         'StandardAnnualAggLossDevModelInput',
         representation(priorForKnotPositions='numeric',
                        priorForNumberOfKnots='numeric'),
         contains='AnnualAggLossDevModelInput')


##' Create an Object of class \code{StandardAnnualAggLossDevModelInput}.
##'
##' The loss development models require extensive input.  Much of the input is directly dependent on the values of other input.
##' As such, this function facilitates much of the work of setting model parameters and determining which output to collect.
##'
##' The function creates an object of class \code{StandardAnnualAggLossDevModelInput}.
##'
##' \describe{
##'   \item{\bold{Inflation Rate}}{
##'     Workers Compensation indemnity statutes can be extremely complex.
##'     In order to provide a sufficient level of flexibility, the model allows for a combination for 3 inflation/escalation rates: one non-stochastic rate, one stochastic rate, and one zero rate.
##'     All inflation rates are on the real scale, that is, a 5 percent increase is entered as \dQuote{0.05.}
##'     Inflation rates are only used to create an expected inflation rate.  Actual calendar year effects deviate from this prior along the diagional.
##'     \describe{
##'       \item{Non-Stochastic Rate of Inflation/Escalation}{
##'         The user must specify a non-stochastic rate of inflation for all observed and future time periods.
##'         Assuming a non-stochastic inflation rate of zero turns this feature off.
##'         \describe{
##'           \item{\code{non.stoch.inflation.rate}}{
##'             The non-stochastic inflation rate is specified by the parameter \code{non.stoch.inflation.rate}.
##'             The default value is set to \dQuote{0.}
##'             There are three possible ways to specify \code{non.stoch.inflation.rate}:
##'             \describe{
##'               \item{If specified as a single value}{
##'               The inflation in every cell is assumed to be this value.
##'               For instance, choosing zero results in no non-stochastic inflation rate.
##'               }
##'               \item{If specified as a vector}{
##'               Its length must be two less than \code{total.exp.years + total.dev.years} and must have \code{names} starting at \code{dimnames(incremental.payments) [[1]][2]} and increasing by one.
##'               The first value corresponds to the inflation rate from the first to second diagonal (or cell 1,1 to cells 1,2 and 2,1).
##'               The second value corresponds to the inflation rate from the second to third diagonal. etc.
##'               }
##'               \item{If specified as a matrix}{
##'               It must have number of rows equal to \code{total.exp.years} and number of columns equal to \code{total.dev.years}.
##'               The row names must be equal to the row names of \code{incremental.payments}.
##'               Each cell represents the inflation to that cell.
##'               Cells in the first column represent inflation down rows.
##'               Cells in other columns represent inflation from left to right.
##'               Cell 1,1 is ignored.
##'               }
##'             }
##'           }
##'           \item{\code{non.stoch.inflation.weight}}{
##'             The user has the option of allowing the percentage of dollars that inflate/escalate at the non-stochastic rate to vary.
##'             This is done by properly setting the value of \code{non.stoch.inflation.weight}.
##'             The default value is set to \dQuote{1.}
##'             \code{non.stoch.inflation.weight} must be at least zero and at most one.  Also \code{non.stoch.inflation.weight + stoch.inflation.weight} must be at most one.
##'             There are three possible ways to specify \code{non.stoch.inflation.weight}:
##'             \describe{
##'               \item{If specified as a single value}{
##'                 The proportion of dollars that inflate at the rate of \code{non.stoch.inflatoin.rate} in every cell is assumed to be this value.
##'                 For instance, choosing zero results in the model assuming that no non-stochastic inflation applies.
##'               }
##'               \item{If specified as a vector}{
##'                 Its length must be \code{total.dev.years}.
##'                 This vector is repeated to produce a matrix such that each row equals this vector.
##'                 This matrix is then assumed to be the supplied for this parameter (See below).
##'                 See next item.
##'               }
##'               \item{If specified as a matrix}{
##'                 It must have number of rows equal to \code{total.exp.years} and number of columns equal to \code{total.dev.years}.
##'                 The row names must be equal to the row names of \code{incremental.payments}.
##'                 Each cell represents the proportion of dollars in the previous cell inflating at the \code{non.stoch.inflation.rate} rate to that cell.
##'                 Cells in the first column represent inflation down rows.
##'                 Cells in other columns represent inflation from left to right.
##'                 Cell 1,1 is ignored.
##'               }
##'             }
##'           }
##'         }
##'       }
##'       \item{Stochastic Rate of Inflation/Escalation}{
##'         The stochastic inflation rate is modeled as an AR1 process (on the log scale).
##'         The user must supply a stochastic rate of inflation for at least all observed diagonals and has the option of supplying values for diagonals beyond the final diagonal in the observed triangle.
##'         As such, the stochasticity of the stochastic rate of inflation only applies to future rates of inflation.
##'         Values not supplied for future diagonals are predicted using the estimated parameters for the supplied inflation series.
##'         If the user does not wish to supply/assume a stochastic rate of inflation, he may set the stochastic inflation rate to zero for all periods
##'         (\code{stoch.inflation.rate=0}).
##'         \describe{
##'           \item{\code{stoch.inflation.rate}}{
##'             The stochastic inflation rate is specified by the parameter \code{stoch.inflation.rate}.
##'             The default value is set to \dQuote{0.}
##'             There are two possible ways to specify \code{stoch.inflation.rate}:
##'             \describe{
##'               \item{If specified as a single 0 (zero)}{
##'                 It is assumed that no stochastic inflation applies.
##'               }
##'               \item{If specified as a vector}{
##'                 Its length must be at least one less than the number of rows in \code{incremental.payments}.
##'                 Future values are simulated as needed.
##'                 (No element in \code{stoch.inflation.rate} may be \code{NA}; only supply values for known rates of inflation.)
##'                 The vector must be named and the names must be consecutive integers.
##'                 The names must contain the second and last values of \code{dimnames(incremental.payments)[[1]]}.
##'                 The names represent the calendar year of inflation.
##'               }
##'             }
##'           }
##'           \item{\code{stoch.inflation.weight}}{
##'             The default value is set to \dQuote{\code{1 - non.stoch.inflation.weight}}.
##'             See \code{non.stoch.inflation.weight}.
##'           }
##'           \item{\code{stoch.inflation.lower.bound}}{
##'             The user has the option of putting a lower bound on the stochastic rate of inflation/escalation.
##'             \code{stoch.inflation.rate} should not account for \code{stoch.inflation.lower.bound} as this would result in incorrectly simulated future inflation rates.
##'             The model properly applies the \code{stoch.inflation.lower.bound} to both observed and future rates of inflation prior to \dQuote{inflating} dollars in the triangle.
##'             The user accounts for lower bounds by properly setting the value of \code{stoch.inflation.lower.bound}.
##'             \code{stoch.inflation.lower.bound} should be on the real scale (\code{x[i] / x[i-1] - 1}).
##'             The default value is set to \code{-1}, which is no bound.
##'             Bounds are applied prior to weights.
##'             There are three possible ways to specify \code{stoch.inflation.weight}:
##'             \describe{
##'               \item{If specified as a single value}{
##'               All stochastic rates of inflation in the triangle are bound by this value.
##'               }
##'               \item{If specified as a vector}{
##'                 Its length must be \code{total.dev.years}.
##'                 This vector is repeated to produce a matrix with each rows equal to it.
##'                 It will then be as if this matrix was supplied for this parameter.
##'                 See next item.
##'               }
##'               \item{If specified as a matrix}{
##'                 It must have number of rows equal to \code{total.exp.years} and number of columns equal to \code{total.dev.years}.
##'                 The row names should be equal to the row names of \code{incremental.payments}.
##'                 Each cell represents the bound on the stochastic rate of inflation from the previous cell to that cell.
##'                 Cells in the first column represent inflation down rows.
##'                 Cells in other columns represent inflation from left to right.
##'                 Cell 1,1 is ignored.
##'               }
##'             }
##'           }
##'           \item{\code{stoch.inflation.upper.bound}}{
##'           Default value is \code{Inf}, which is no bound.
##'           See \code{stoch.inflation.lower.bound}.
##'           }
##'           \item{\code{known.stoch.inflation.mean}}{
##'           The AR1 process used to simulate future stochastic rates of inflation has a non-zero mean which is (at least approximately) equal to the historical mean.
##'           If for some reason (i.e., the series of inflation rates is too short) the user believes this historical mean poorly represents future rates of inflation, the user can override the estimation process.
##'           If \code{known.stoch.inflation.mean} is set to \code{NA}, the mean is estimated.  Otherwise the mean is assumed (with certainty) to be the specified value.
##'           (Note that since the estimation takes place on the log scale, the logarithm of \code{known.stoch.inflation.mean} plus 1 is used as the mean on the log scale;
##'           this results in \code{known.stoch.inflation.mean} being the geometric mean.
##'           }
##'           \item{\code{known.stoch.inflation.persistence}}{
##'           The user has the option of overriding the estimation of the \eqn{rho} coefficient in the AR1 process used to simulate future stochastic rates of inflation.
##'           If \code{known.stoch.inflation.persistence} is set to \code{NA}, the \eqn{rho} coefficient is estimated; otherwise it is taken as this value.
##'           }
##'         }
##'       }
##'     }
##'   }
##'   \item{\bold{Projected Rate of Decay}}{
##'     The model estimates an inflation-free rate of decay up to the size of the triangle (\code{dim(incremental.payments)[2]}).
##'     Since no data is available to the model beyond this point, the model must be supplied with an assumption for future values.
##'     This assumption is supplied via the parameter \code{projected.rate.of.decay} in one of two ways:
##'     \describe{
##'       \item{A single value of \code{NA}}{
##'         The inflation-free rate of decay for the final development year in the triangle (\code{dim(incremental.payments)[2]}) is used for all future development years through development year \code{total.dev.years}.
##'       }
##'       \item{A named list}{
##'         This is currently a rather low level interface allowing for much flexibility at the expense of some ease of usability.
##'         (In the (possibly near) future, it is likely that the burden placed on the user will be eased.)
##'         The named list must contain three elements (with an optional fourth element):
##'         \describe{
##'           \item{\code{last.non.zero.payment}}{
##'             The named element \code{last.non.zero.payment} refers to the last development year in which a non-zero payment may occur.
##'             It is assumed that all subsequent incremental payments are exactly zero.
##'             Must be at most \code{total.dev.years} and at least \code{dim(incrementals)[2]+2}.
##'             \code{last.non.zero.payment} must be specified as a named numeric vector of length \code{total.exp.years}.
##'             Element one then refers to the first exposure year; and element two to the second.
##'             The names must start with \code{dimnames(incrementals)[[1]][1]} and increment by 1.
##'             \code{last.non.zero.payment} can be omitted.  If omitted, then the \code{last.non.zero.payment} for each exposure year is \code{total.dev.years}.
##'           }
##'           \item{\code{final.rate.of.decay.mean}}{
##'             The parameter \code{final.rate.of.decay.mean} refers to the rate of decay in \code{last.non.zero.payment}.
##'             As the rate of decay is estimated on the log scale, \code{final.rate.of.decay.mean} is first increased by 1 and then the log of the resulting value is taken; the result is then treated as the mean on the log scale.
##'             This results in \code{final.rate.of.decay.mean} being the geometric mean.
##'             (\code{final.rate.of.decay.mean=0.05} means that the mean of the final rate of decay (on the log scale) is \code{log(1.05)}.)
##'             This may be supplied as a single numeric value, in which case all exposure years are assumed to have the same \code{final.rate.of.decay.mean}.
##'             Or this may be supplied as a named numeric vector of length \code{total.exp.years}.
##'             Element one then refers to the first exposure year, and element two refers to the second, etc.
##'             The names of the vector must start with \code{dimnames(incrementals)[[1]][1]} and increment by \code{1}.
##'             (A Technical Note: \code{last.non.zero.payment} actually defines the final rate of decay since the rate of decay from it to the next development year is \code{-1.00}.
##'             Thus a better name for \code{final.rate.of.decay.mean} does not refer to the \emph{final} rate of decay but rather the \emph{penultimate} rate of decay.)
##'           }
##'           \item{\code{final.rate.of.decay.sd}}{
##'             The parameter \code{final.rate.of.decay.sd} refers to the standard deviation around \code{final.rate.of.decay.mean}.
##'             This may be supplied as a single numeric value, in which case all exposure years are assumed to have the same \code{final.rate.of.decay.sd}.
##'             Or this may be supplied as a named numeric vector of length \code{total.exp.years}.
##'             Element one then refers to the first exposure year, and element two refers to the second, etc.
##'             The names of the vector must start with \code{dimnames(incrementals)[[1]][1]} and increment by \code{1}.
##'           }
##'           \item{\code{final.rate.of.decay.weight}}{
##'             The basic concept behind a user-supplied rate of decay is that of a weighted average between the final estimated rate of decay (the one associated with \code{dim(incrementals)[2]}) and
##'              the final supplied value \code{final.rate.of.decay.mean}.
##'             Thus the user must supply a vector of weights to describe the path from the estimated value to the projected value.
##'             The \dQuote{weight} is the weight (on the log scale) of the user-supplied \code{final.rate.of.decay}.
##'             So, a value of \code{1} means that the projected rate of decay should take on \code{final.rate.of.decay} and a value of \code{0} means that the projected rate of decay should equal the final estimated rate of decay.
##'             Note that while these weights are allowed to be outside the interval \code{[0,1]}, convention should dictate they remain within this interval.
##'             These weights are supplied by specifying a value for \code{final.rate.of.decay.weight}, which may be specified in two ways.
##'             It may be specified as a numeric matrix with number of columns equal to \code{total.dev.years - dim(incremental.payments)[2]} and number of rows equal to \code{total.exp.years},
##'             in which case the first row corresponds to the first exposure year in \code{incremental.payments} and the first column corresponds to the development year immediately following the last one in \code{incremental.payments}.
##'             It may be specified as a numeric vector of length \code{total.dev.years - dim(incremental.payments)[2]}, in which case it is assumed to be the same for all exposure years.
##'           }
##'         }
##'       }
##'     }
##'   }
##'   \item{\bold{Cumulative Triangle}}{
##'     As a safeguard, the cumulative and incremental triangles must match up.  They are tested by taking the cumulative triangle and then decumulating it with the function \code{decumulate}.
##'     The decumulated triangle is then compared to the incremental triangle.  Except for \code{NA}'s in either triangle, the values in the triangles must match.  (Errors due to numerical precision are accounted for.)
##'   }
##'   \item{\bold{Measurment Error Second-Order Random Walk}}{
##'     The model allows for changes in the scale of the measurement error with development time.
##'     Changes in the scale parameter (see \code{\link{scaleParameter}}) are assumed to follow a second-order random walk on the log scale.
##'     The model allows for imposed stationarity in the scale parameter by setting of the value \code{last.column.with.scale.innovation}.
##'     \code{last.column.with.scale.innovation} must be coercible to an integer, be least 1, and at most the number of columns in \code{incremental.payments}.
##'     A value of 1 constrains the scale parameter in each column to the same value.
##'     A value of \code{dim(incremental.payments)[2]} allows for all columns to have their own scale parameter (but smoothed by the second-order random walk).
##'     The effective value used by the model for this argument is the minimum of the supplied value and the last column in \code{incremental.payments} with a non-missing value (on the log scale).
##'     Note that since the scale parameter is assumed to follow a second-order random walk, a value of \code{2} results in an (effectively) unconstrained estimation for the scale parameter in the first column.
##'   }
##'   \item{\bold{Number of Knots}}{
##'     The consumption path is modeled (on the log scale) as a linear spline.
##'     The number of knots in the spline is endogenous to the model estimation.
##'     This allows the model to adapt to loss triangles of varying complexity.
##'     In order to allow the user to give more credibility to models of a lower complexity (and possibly avoid overfitting), a truncated negative binomial prior is placed on the number of knots.
##'     The truncation point is adjusted with the size of the triangle.
##'     The paramters for this negative binomial can be specified with the argument \code{prior.for.number.of.knots}, which should be a vector of two elements.
##'     The first element should be any realy number greater than zero.
##'     If this number is an integer, it can be interpreted as the number of failures until the experiment is stopped.
##'     The second parameter should be a real number greater than 0 but less than 1 and represents the success probability.
##'     The mean of the (un-truncated) negative binomial is then \code{prior.for.number.of.knots[1] * prior.for.number.of.knots[2]/ (1 - prior.for.number.of.knots[2])}.
##'     And the variance is \code{prior.for.number.of.knots[1] * prior.for.number.of.knots[2]/ (1 - prior.for.number.of.knots[2]) ^ 2}.
##'   }
##' }
##'
##' @param incremental.payments A square matrix of incremental payments.  Row names must correspond to the exposure year. Only the upper-left (including the diagonal) of this matrix may have non-missing values.  Lower-right must be \code{NA}.
##' @param extra.exp.years A single integer value (\code{total.exp.years} overrides) greater than or equal to 1 (default is 1) specifying the number of additional exposure years (or rows in the triangle) to project forward.
##' @param extra.dev.years A single integer value (\code{total.dev.years} overrides) greater than or equal to 1 (default is 1) specifying the additional number of development years (or columns in the triangle) to project forward.
##'
##' @param non.stoch.inflation.rate May be one of three types (See \emph{Inflation Rate} in Details): 1) A single numeric value; 2) a vector of numerics (of specific length); 3) a matrix of numerics (of specific dim).
##' @param non.stoch.inflation.weight May be one of three types (See \emph{Inflation Rate} in Details): 1) A single numeric value; 2) a vector of numerics (of specific length); 3) a matrix of numerics (of specific dim).
##'
##' @param stoch.inflation.rate May be one of two types (See \emph{Inflation Rate} in Details): 1) A single numeric value of \emph{zero}; 2) a vector of numerics (of specific length).
##' @param stoch.inflation.weight May be one of three types (See \emph{Inflation Rate} in Details): 1) A single numerical value; 2) a vector of numerics (of specific length); 3) a matrix of numerics (of specific dim).
##' @param stoch.inflation.lower.bound  May be one of three types (See \emph{Inflation Rate} in Details): 1) A single numeric value; 2) a vector of numerics (of specific length); 3) a matrix of numerics (of specific dim).
##' @param stoch.inflation.upper.bound  May be one of three types (See \emph{Inflation Rate} in Details): 1) A single numeric value; 2) a vector of numerics (of specific length); 3) a matrix of numerics (of specific dim).
##' @param known.stoch.inflation.mean   May be one of two types (See \emph{Inflation Rate} in Details): 1) A single numeric value; 2) \code{NA}.
##' @param known.stoch.inflation.persistence  May be one of two types (See \emph{Inflation Rate} in Details): 1) A single numeric value; 2) \code{NA}.
##'
##' @param total.exp.years A single integer value (overrides \code{extra.exp.years}) specifying the last exposure year to project forward.  Must be at least the number of rows in \code{incremental.payments} + 1.
##' @param total.dev.years A single integer value (overrides \code{extra.dev.years}) specifying the last development year to project forward.  Must be at least the number of columns in \code{incremental.payments} + 1 .
##'
##' @param cumulative.payments A numeric matrix with the same dim and dimnames as \code{incremental.payments}.  Must be a possible cumulative payment triangle of \code{incremental.payments}.  (See \emph{Cumulative Payments} Section.)
##'
##' @param exp.year.type A single character value indicating the type of exposure years:  \sQuote{ambiguous}, \sQuote{py}, and \sQuote{ay} mean \sQuote{Exposure Year}, \sQuote{Policy Year}, and \sQuote{Accident Year}; respectively.
##' @param prior.for.knot.locations A single numeric value of at least 1.  The prior for the location of knots is a scaled beta with parameters \code{c(1,prior.for.knot.locations)}.  Large values produce stable consumption paths at high development years.
##' @param prior.for.number.of.knots A two element vector giving the paramters for the prior number of knots. (See \emph{Number of Knots} in Details)
##' @param use.skew.t A logical value.  If \code{TRUE}, the model assumes that the observed and estimated log incremental payments are realizations from a skewed \eqn{t} distribution; if \code{FALSE} it assumes zero skewness. (See Reference.)
##' @param bound.for.skewness.parameter A positive numerical value representing the symetric boundaries for the skewness parameter.  In most cases, the default should be sufficient. Ignored if \code{use.skew.t=FALSE}.
##' @param last.column.with.scale.innovation A single integer. Must be at least 1 and at most the number of columns in \code{incremental.payments}.  See \emph{Measurment Error-Second Order Random Walk} in Details.
##' @param use.ar1.in.calendar.year A logical value.  The calendar year effect errors may (at the users discretion) include an autoregressive process of order 1.  \code{TRUE} turns on the ar1 process, \code{FALSE} (the Default) turns it off.
##' @param use.ar1.in.exposure.growth A logical value.  The exposure growth errors may (at the users discretion) include an autoregressive process of order 1.  \code{TRUE} (the Default) turns on the ar1 process, \code{FALSE} turns it off.
##'
##' @param projected.rate.of.decay  May be one of three types (See \emph{Projected Rate of Decay} in Details): 1) \code{NA}; 2) a matrix of numerics (of specific dim); 3) a named list.

##'

##'
##' @references
##' Kim, Y., and J. McCulloch (2007) \dQuote{The Skew-Student Distribution with Application to U.S. Stock Market Returns and the Equity Premium,} Department of Economics, Ohio State University, October 2007
##'
##'
##' @return An object of class \code{AggModelInput}.  The model specified by the returned object must then be estimated using the function \code{runLossDevModel}.
##' @export
##' @usage
##' makeStandardAnnualInput(
##'   incremental.payments=decumulate(cumulative.payments),
##'   extra.dev.years=1,
##'   extra.exp.years=1,
##'   non.stoch.inflation.rate=0,
##'   non.stoch.inflation.weight=1,
##'   stoch.inflation.rate=0,
##'   stoch.inflation.weight=1-non.stoch.inflation.weight,
##'   stoch.inflation.lower.bound=-1,
##'   stoch.inflation.upper.bound=Inf,
##'   known.stoch.inflation.mean=NA,
##'   known.stoch.inflation.persistence=NA,
##'   total.dev.years=extra.dev.years+dim(incremental.payments)[2],
##'   total.exp.years=extra.exp.years+dim(incremental.payments)[1],
##'   cumulative.payments=cumulate(incremental.payments),
##'   exp.year.type=c('ambiguous', 'py', 'ay'),
##'   prior.for.knot.locations=2,
##'   prior.for.number.of.knots=c(3, 1/7),
##'   use.skew.t=FALSE,
##'   bound.for.skewness.parameter=10,
##'   last.column.with.scale.innovation=dim(incremental.payments)[2],
##'   use.ar1.in.calendar.year=FALSE,
##'   use.ar1.in.exposure.growth=TRUE,
##'   projected.rate.of.decay=NA)
##'
##'
makeStandardAnnualInput <- function(incremental.payments=decumulate(cumulative.payments),
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
                                    prior.for.number.of.knots=c(3, 1/7),
                                    use.skew.t=FALSE,
                                    bound.for.skewness.parameter=10,
                                    last.column.with.scale.innovation=dim(incremental.payments)[2],
                                    use.ar1.in.calendar.year=FALSE,
                                    use.ar1.in.exposure.growth=TRUE,
                                    projected.rate.of.decay=NA)
{

    force(incremental.payments)
    force(cumulative.payments)

    ans <- new('StandardAnnualAggLossDevModelInput')

    ##set the model file and output type
    ans@modelFile <- 'standard.model.txt'
    ans@outputType <- 'StandardAnnualAggLossDevModelOutput'

    ##we shouldn't have to check incrementals.payments because validObject will do that for us
    ##but we do some tests so we can make sure to supply meaningful error messages

    if(!is.matrix(incremental.payments))
        stop('"is.matrix(incremental.payments)" did not return true. "incremental.payments" must be a square numeric MATRIX')
    if(!is.numeric(incremental.payments))
        stop('"is.numeric(incremental.payments)" did not return true. "incremental.payments" must be a square NUMERIC matrix')
    if(dim(incremental.payments)[1] != dim(incremental.payments)[2])
         stop('"dim(incremental.payments)[1] != dim(incremental.payments)[2]".  "incremental.payments" must be a SQUARE numeric matrix')
    if(dim(incremental.payments)[1] <= 5)
        stop('"dim(incremental.payments)[1] <= 5". "incremental.payments" must be of dim >= 6, please supply a larger triangle.')
    if(any(is.infinite(incremental.payments)))
        stop('"incremental.payments" must be finite')

    if(any(na.omit(as.vector(incremental.payments <= 0))))
        warning('"incremental.payments" contains non-positive values.  These will be treated as missing an can thus result it an over estimated ultimate loss.  Be sure to check how predicted and actual cumulatives line up.')
     if(any(na.omit(as.vector(incremental.payments == 0))))
         message('"incremental.payments" contains incremental payments of zero.  You may wish to use the function "acccountForZeroPayments".')

    i <- rep(1:dim(incremental.payments)[1], dim(incremental.payments)[1])
    j <- rep(1:dim(incremental.payments)[1], rep(dim(incremental.payments)[1],dim(incremental.payments)[1]))
    if(any(!is.na(incremental.payments[i + j - 1 > dim(incremental.payments)[1]])))
        stop('The lower right half of "incremental.payments" must ALL be "NA".')
    rm(i,j)

    ##check to make sure dimnames(incremental.payments) is correct
    if(is.null(dimnames(incremental.payments)) || is.null(dimnames(incremental.payments)[[1]]))
        stop('dimnames(incremental.payments)[[1]]" cannot be "NULL"')
    tmp.exp.years <- dimnames(incremental.payments)[[1]]
    suppressWarnings(tmp.exp.years <- as.numeric(tmp.exp.years))
    if(any(is.na(tmp.exp.years)) || any(as.integer(tmp.exp.years) != tmp.exp.years))
        stop('"dimnames(incremental.payments)[[1]]" must contain integer values stored as charcters. c(1990, 1991, ...)')
    if(is.unsorted(tmp.exp.years) || any(tmp.exp.years[-1] - tmp.exp.years[-length(tmp.exp.years)] != 1))
        stop('"dimnames(incremental.payments)[[1]]" must contain consecutive integer values stored as charcters. c(1990, 1991, ...)')
    ans@exposureYears <- as.integer(tmp.exp.years)
    rm(tmp.exp.years)

    ##warning('need to check "incremental.payments" to ensure it has "enough" values which will be non-missing after taking the log')
    ##for now we just let it crash if the user supplies poorly populated triangles

    ans@incrementals <- incremental.payments
    dimnames(ans@incrementals) <- NULL

    ##check consistency of total.dev.years and  total.exp.year
    ##first total.dev.years
    if(!is.vector(total.dev.years, mode='numeric') || length(total.dev.years) != 1)
        stop('"total.dev.years" must be a numeric vector and have length 1')
    if(is.infinite(total.dev.years))
        stop('"total.dev.years" must be finite')
    if(is.na(total.dev.years))
        stop('"total.dev.years" cannot be NA')
    if(!is.numeric(total.dev.years) || as.integer(total.dev.years)!=total.dev.years)
        stop('"total.dev.years" must be numeric and coercible to an integer')
    if(total.dev.years < dim(incremental.payments)[2] + 1)
        stop('"total.dev.years < dim(incremental.payments)[2] + 1".  "total.dev.years" must be at least dim(incremental.payments)[2] + 1')
    ans@totalDevYears <- as.integer(total.dev.years)

    ##next total.exp.years
    if(!is.vector(total.exp.years, mode='numeric') || length(total.exp.years) != 1)
        stop('"total.exp.years" must be a vector and have length 1')
    if(is.infinite(total.exp.years))
        stop('"total.exp.years" must be finite')
    if(is.na(total.exp.years))
        stop('"total.exp.years" cannot be NA')
    if(!is.numeric(total.dev.years) || as.integer(total.exp.years) != total.exp.years)
        stop('"total.exp.years" must be numeric and coercible to an integer')
    if(total.exp.years < dim(incremental.payments)[1])
        stop('"total.exp.years < dim(incremental.payments)[1] + 1".  "total.exp.years" must be at least dim(incremental.payments)[1] + 1')
    ans@totalExpYears <- as.integer(total.exp.years)


    extract.matrix.from.1.vector.or.matrix <- function(arg, horz.copy)#if horz.copy and get(arg) is a vector, a matrix with rows equal to get(arg) will be created, otherwise it will be assumed to be a diagonal effect
    {

        arg.v <- get(arg, parent.frame(), inherits=FALSE)

         ##check somethings for arg.v regardless of how it is supplied
        if(is.null(arg.v))
            stop(paste('"', arg, '" cannot be NULL', sep=''))
        if(!is.numeric(arg.v))
            stop(paste('"', arg,'" must be numeric', sep=''))
        if(length(arg.v) < 1)
            stop(paste('"', arg, '" must have length of at least 1', sep=''))
        if(any(is.na(arg.v)))
            stop(paste('"', arg, '" cannot be NA', sep=''))

        ##now we account for three possible ways to enter arg.
        if(is.vector(arg.v, mode='numeric') && length(arg.v) == 1) {
            ## arg.v is a single value and should be used for all inflation rates
            ans <- array(arg.v, c(total.exp.years, total.dev.years))
            return(ans)

        } else if(is.vector(arg.v, mode='numeric')) {
            ## arg.v is a vector,
            if(horz.copy)
            {
                ##arg.v should be copied down columns
                ##first element is for first column, second is for second, etc
                if(length(arg.v) != total.dev.years)
                    stop(paste('If "',arg,'" is supplied as a vector, it must be of length "total.dev.years".  See documentation for this function.', sep=''))

                ans <- rep(arg.v, rep(total.exp.years, total.dev.years))
                dim(ans) <- c(total.exp.years, total.dev.years)
                return(ans)
            } else {
                ##arg.v applies to diagonals
                ##first element is for first diagional, second is for second, etc
                if(length(arg.v) != total.dev.years + total.exp.years - 2)
                    stop(paste('If "',arg,'" is supplied as a vector, it must be of length "total.dev.years + total.exp.years - 2".  See documentation for this function.', sep=''))
                if(is.null(names(arg.v)))
                    stop(paste('"names(', arg, '" cannot be NULL if supplied as a vector. See documentation of this function.', sep=''))
                tmp.years <- names(arg.v)
                suppressWarnings(tmp.years <- as.numeric(tmp.years))
                if(any(is.na(tmp.years)) || any(as.integer(tmp.years) != tmp.years))
                    stop(paste('"names(', arg,')" must contain integer values stored as charcters. c(1990, 1991, ...)', sep=''))
                if(is.unsorted(tmp.years) || any(tmp.years[-1] - tmp.years[-length(tmp.years)] != 1) || tmp.years[1] != ans@exposureYears[2])
                    stop(paste('"names(', arg,')" must contain consecutive integer values stored as charcters. c(1990, 1991, ...) and start at dimnames(incremental.payments)[[1]][2]', sep=''))
                rm(tmp.years)

                i <- rep(1:total.exp.years, total.dev.years)
                j <- rep(1:total.dev.years, rep(total.exp.years,total.dev.years))
                arg.v.0 <- c(0, arg.v)
                ans <- arg.v.0[i + j - 1]
                dim(ans) <- c(total.exp.years, total.dev.years)
                return(ans)
            }
        } else if(is.matrix(arg.v)) {
            ## arg.v is a matrix and should simply be copied
            if(dim(arg.v)[1] != total.exp.years || dim(arg.v)[2] != total.dev.years)
                stop(paste('If "', arg,'" is supplied as a matrix, it must be of dim "total.exp.years" and "total.dev.years".  See documentation for this function.', sep=''))
            if(is.null(dimnames(arg.v)) || is.null(dimnames(arg.v)[[1]]))
                stop(paste('"dimnames(', arg,')[[1]]" cannot be NULL if it is supplied as a matrix.', sep=''))
            tmp.years <- dimnames(arg.v)[[1]]
            suppressWarnings(tmp.years <- as.numeric(tmp.years))
            if(any(is.na(tmp.years)) || any(as.integer(tmp.years) != tmp.years))
                stop(paste('"dimnames(', arg, ')[[1]]" must contain integer values stored as charcters. c(1990, 1991, ...)', sep=''))
            if(is.unsorted(tmp.years) || any(tmp.years[-1] - tmp.years[-length(tmp.years)] != 1) || tmp.years[1] != ans@exposureYears[1])
                stop(paste('"dimnames(', arg, ')[[1]]" must contain consecutive integer values stored as charcters. c(1990, 1991, ...) and start at dimnames(incremental.payments)[[1]][1]', sep=''))
            rm(tmp.years)
            ans <- arg.v
            dimnames(ans) <- NULL
            return(ans)

        } else {
            ## don't know what non.stoch.inflation.rate is
                stop(paste('"', arg, '" must either be a vector of lenght 1,',
                           'a vector of length "', ifelse(horz.copy, 'total.dev.years', 'total.dev.years + total.exp.years - 2'),'",',
                           'or a matrix of dim "total.exp.years" and "total.dev.years".  See documentation for this function.', sep=''))
        }
    }

    ##check somethings for non.stoch.inflation.rate regardless of how it is supplied
    if(any(non.stoch.inflation.rate <= -1))
        stop('"non.stoch.inflation.rate" must be > -1')
    if(any(is.infinite(non.stoch.inflation.rate)))
        stop('"non.stoch.inflation.rate" must be finite')
    if(any(abs(non.stoch.inflation.rate) > 0.5))
        warning('You have specified a non-stochastic rate of inflation of more than 50% in magnitude.  Are you sure this is correct?')
    ans@nonStochInflationRate <- extract.matrix.from.1.vector.or.matrix(arg='non.stoch.inflation.rate', horz.copy=FALSE)

    ##check somethings for non.stoch.inflation.weight regardless of how it is supplied
    if(any(non.stoch.inflation.weight > 1) || any(non.stoch.inflation.weight < 0))
        stop('"non.stoch.inflation.weight" must be in [0,1]')
    ans@nonStochInflationWeight <- extract.matrix.from.1.vector.or.matrix(arg='non.stoch.inflation.weight', horz.copy=TRUE)


    ##pull out stoch.inflation.rate, It can be zero or a vector
    if(!is.numeric(stoch.inflation.rate))
        stop('"stoch.inflation.rate" must be numeric')
    if(any(is.na(stoch.inflation.rate)))
        stop('"stoch.inflation.rate" cannot be NA')
    if(any(stoch.inflation.rate <= -1))
        stop('"stoch.inflation.rate" must be > -1')
    if(any(is.infinite(stoch.inflation.rate)))
        stop('"stoch.inflation.rate" must be  finite')
    if(any(abs(stoch.inflation.rate) > 0.5))
        warning('You have specified a stochastic rate of inflation of more than 50% in magnitude.  Are you sure this is correct?')

    if(identical(stoch.inflation.rate,0))
    {
        ##by convention, stochInflationRate == 0 will mean to turn of stoch inflation
        ans@stochInflationRate <- 0
        ans@stochInflationYears <- integer()
        ans@stochInflationWeight <- array(NA, c(0,0))
        ans@stochInflationLowerBound <- array(NA, c(0,0))
        ans@stochInflationUpperBound <- array(NA, c(0,0))
        ans@knownStochInflationMean <- numeric()
        ans@knownStochInflationPersistence <- numeric()

    } else if (is.vector(stoch.inflation.rate, mode='numeric')) {
        if(all(stoch.inflation.rate == 0))
            warning('"stoch.inflation.rate" contains more than 1 element and all of these elements are (0) zero.  If you wish not to specify a stochastic inflation rate, you should set this argument to a single zero.')

        tmp.years <- names(stoch.inflation.rate)
        suppressWarnings(tmp.years <- as.numeric(tmp.years))
        if(any(is.na(tmp.years)) || any(as.integer(tmp.years) != tmp.years))
            stop('"names(stoch.inflation.rate)" must contain integer values stored as charcters. c(1990, 1991, ...)')
        if(is.unsorted(tmp.years) || any(tmp.years[-1] - tmp.years[-length(tmp.years)] != 1) || any(!(ans@exposureYears %in% tmp.years)))
            stop('"names(stoch.inflation.rate)" must contain consecutive integer values stored as charcters. c(1990, 1991, ...) and contain dimnames(incremental.payments)[[1]]')
        ans@stochInflationYears <- as.integer(tmp.years)
        rm(tmp.years)

        ans@stochInflationRate <- stoch.inflation.rate
        names(ans@stochInflationRate) <- NULL

        ##check somethings for stoch.inflation.weight regardless of how it is supplied
        if(any(stoch.inflation.weight > 1) || any(stoch.inflation.weight < 0))
            stop('"stoch.inflation.weight" must be in [0,1]')
        ans@stochInflationWeight <- extract.matrix.from.1.vector.or.matrix(arg='stoch.inflation.weight', horz.copy=TRUE)

        ##check to make sure that the two weights don't exceed 1
        if(any(ans@stochInflationWeight + ans@nonStochInflationWeight > 1))
            stop('"non.stoch.inflation.weight + stoch.inflation.weight" must be at most 1')

        ##check somethings for stoch.inflation.lower.bound regardless of how it is supplied
        if(any(stoch.inflation.lower.bound < -1))
            stop('"stoch.inflation.lower.bound" must be at least -1')
        if(any(!is.finite(stoch.inflation.lower.bound)))
            stop('"stoch.inflation.lower.bound" must be finite')
        ans@stochInflationLowerBound <- extract.matrix.from.1.vector.or.matrix(arg='stoch.inflation.lower.bound', horz.copy=TRUE)
        if(any(stoch.inflation.upper.bound <= -1))
            stop('"stoch.inflation.upper.bound" must be greater than -1')
        ans@stochInflationUpperBound <- extract.matrix.from.1.vector.or.matrix(arg='stoch.inflation.upper.bound', horz.copy=TRUE)

        ##check to make sure the bounds make sense
        if(any(ans@stochInflationLowerBound > ans@stochInflationUpperBound))
            stop('"stoch.inflation.lower.bound" must be less than or equal to "stoch.inflation.upper.bound"')


        ## known.stoch.inflation.mean
        if(!is.na(known.stoch.inflation.mean))
        {
            if(!is.vector(known.stoch.inflation.mean, 'numeric') || length(known.stoch.inflation.mean) != 1)
                stop('"known.stoch.inflation.mean" must be a vector and have length 1')
            if(is.infinite(known.stoch.inflation.mean))
                stop('"known.stoch.inflation.mean" must be finite or NA')
            if(!is.numeric(known.stoch.inflation.mean))
                stop('"known.stoch.inflation.mean" must be numeric or NA')
            if(known.stoch.inflation.mean <= -1)
                stop('"known.stoch.inflation.mean <= -1".  "known.stoch.inflation.mean" must be > -1')
        }
        ans@knownStochInflationMean <- as.numeric(known.stoch.inflation.mean)

        ## known.stoch.inflation.persistence
        if(!is.na(known.stoch.inflation.persistence))
        {
            if(!is.vector(known.stoch.inflation.persistence, mode='numeric') || length(known.stoch.inflation.persistence) != 1)
                stop('"known.stoch.inflation.persistence" must be a vector and have length 1')
            if(is.infinite(known.stoch.inflation.persistence))
                stop('"known.stoch.inflation.persistence" must be finite or NA')
            if(!is.numeric(known.stoch.inflation.persistence))
                stop('"known.stoch.inflation.persistence" must be numeric or NA')
            if(known.stoch.inflation.persistence < 0 || known.stoch.inflation.persistence >= 1)
                stop('"known.stoch.inflation.persistence" must be in [0,1)')
        }
        ans@knownStochInflationPersistence <- as.numeric(known.stoch.inflation.persistence)



    } else {
        stop('"stoch.inflation.rate" must either be zero or a vector with a length of at least "dim(incremental.payments)[1] - 1"')
    }


    ##check cumulative.payments we will do this by creating an incremental triangle from the cumulative triangle and then comparing the two
    if(!is.matrix(cumulative.payments) || !is.numeric(cumulative.payments))
        stop('"cumulative.payments" must be a numeric matrix')
    if(!identical(dim(cumulative.payments), dim(incremental.payments)))
        stop('"dim(cumulative.payments)" must be the same as "dim(incremental.payments)"')
    if(!identical(dimnames(cumulative.payments), dimnames(incremental.payments)))
        stop('"dimnames(cumulative.payments)" must be the same as "dimnames(incremental.payments"')
    if(any(is.infinite(cumulative.payments)))
        stop('"cumulative.payments" must be finite')
    i <- rep(1:dim(cumulative.payments)[1], dim(cumulative.payments)[1])
    j <- rep(1:dim(cumulative.payments)[1], rep(dim(cumulative.payments)[1],dim(cumulative.payments)[1]))
    if(any(!is.na(cumulative.payments[i + j - 1 > dim(cumulative.payments)[1]])))
        stop('The lower right half of "cumulative.payments" must ALL be "NA".')
    rm(i,j)

    inc.from.cumulative <- cumulative.payments[,1]
    for(j in 2:dim(cumulative.payments)[2])
        inc.from.cumulative <- cbind(inc.from.cumulative, cumulative.payments[,j] - cumulative.payments[,j-1])
    if(!all.equal(inc.from.cumulative[!is.na(inc.from.cumulative) & !is.na(incremental.payments)], incremental.payments[!is.na(inc.from.cumulative) & !is.na(incremental.payments)]))
        stop('"cumulative.payments" is not a possible cumulative triangle of "incremental.payments".')

    ans@cumulatives <- cumulative.payments
    dimnames(ans@cumulatives) <- NULL

    ##Rate of Decay
    if(identical(projected.rate.of.decay, NA)) {
        ans@lastNonZeroPayment <- as.integer(rep(total.dev.years, total.exp.years))
        ans@finalRateOfDecayMean <- rep(0, total.exp.years)
        ans@finalRateOfDecaySD <- rep(1, total.exp.years)
        ans@finalRateOfDecayWeight <- array(0,c(total.exp.years, total.dev.years - dim(incremental.payments)[2]))
    } else if (is.list(projected.rate.of.decay)) {

        if(is.null(names(projected.rate.of.decay)) || any(!(c('final.rate.of.decay.mean', 'final.rate.of.decay.sd', 'final.rate.of.decay.weight') %in% names(projected.rate.of.decay))))
            stop('If "projected.rate.of.decay" is specified as a list, it must be named and those names must contain, "final.rate.of.decay.mean", "final.rate.of.decay.sd", and "final.rate.of.decay.weight"')

        if(length(projected.rate.of.decay) != 3 && (length(projected.rate.of.decay) != 4 && !('last.non.zero.payment' %in% names(projected.rate.of.decay))))
            stop('If "projected.rate.of.decay" is specified as a list, its names must only contain, "final.rate.of.decay.mean", "final.rate.of.decay.sd", and "final.rate.of.decay.weight" (and optionally "last.non.zero.payment")')

        ##last.non.zero.payment
        if('last.non.zero.payment' %in% names(projected.rate.of.decay))
        {
            last.non.zero.payment <- projected.rate.of.decay[['last.non.zero.payment']]
            if(!is.vector(last.non.zero.payment, mode='numeric') || length(last.non.zero.payment) != total.exp.years)
                stop('"projected.rate.of.decay[[\'last.non.zero.payment\']]" must be a numeric vector of length "total.exp.years"')
            if(any(is.infinite(last.non.zero.payment)))
                stop('"projected.rate.of.decay[[\'last.non.zero.payment\']]" must be finite')
            if(any(is.na(last.non.zero.payment)))
                stop('"projected.rate.of.decay[[\'last.non.zero.payment\']]" must not be NA')
            if(!is.numeric(last.non.zero.payment) || any(as.integer(last.non.zero.payment) != last.non.zero.payment))
                stop('"projected.rate.of.decay[[\'last.non.zero.payment\']]" must be numeric and coercible to an integer')
            if(any(last.non.zero.payment > total.dev.years) || any(last.non.zero.payment < dim(incremental.payments)[2] + 1))
                stop('"projected.rate.of.decay[[\'last.non.zero.payment\']]" must be at least dim(incremental.payments)[2] + 1 and at most "total.dev.years"')

            if(is.null(names(last.non.zero.payment)))
                stop('"names(projected.rate.of.decay[[\'last.non.zero.payment\']])" cannot be NULL if supplied as a vector. See documentation of this function.')
            tmp.years <- names(last.non.zero.payment)
            suppressWarnings(tmp.years <- as.numeric(tmp.years))
            if(any(is.na(tmp.years)) || any(as.integer(tmp.years) != tmp.years))
                stop('"names(projected.rate.of.decay[[\'last.non.zero.payment\']])" must contain integer values stored as charcters. c(1990, 1991, ...)')
            if(is.unsorted(tmp.years) || any(tmp.years[-1] - tmp.years[-length(tmp.years)] != 1) || tmp.years[1] != ans@exposureYears[1])
                stop('"names(projected.rate.of.decay[[\'last.non.zero.payment\']])" must contain consecutive integer values stored as charcters. c(1990, 1991, ...) and start at dimnames(incremental.payments)[[1]][1]')
            rm(tmp.years)


        } else {
            last.non.zero.payment <- rep(total.dev.years, total.exp.years)
        }
        ans@lastNonZeroPayment <- as.integer(last.non.zero.payment)

        ##final.rate.of.decay.mean
        final.rate.of.decay.mean <- projected.rate.of.decay[['final.rate.of.decay.mean']]
        if(!is.vector(final.rate.of.decay.mean, mode='numeric') || (length(final.rate.of.decay.mean)!= 1 && length(final.rate.of.decay.mean) != total.exp.years))
            stop('"projected.rate.of.decay[[\'final.rate.of.decay\']]" must be a numeric vector of length "total.exp.years" or 1')
        if(any(is.infinite(final.rate.of.decay.mean)))
            stop('"projected.rate.of.decay[[\'final.rate.of.decay.mean\']]" must be finite')
        if(any(is.na(final.rate.of.decay.mean)))
            stop('"projected.rate.of.decay[[\'final.rate.of.decay.mean\']]" must not be NA')
        if(!is.numeric(final.rate.of.decay.mean))
            stop('"projected.rate.of.decay[[\final.rate.of.decay.mean\']]" must be numeric')
        if(any(final.rate.of.decay.mean <= -1))
            stop('"projected.rate.of.decay[[\'final.rate.of.decay.mean\']]" must be at > -1')

        if(length(final.rate.of.decay.mean) == 1)
        {
            final.rate.of.decay.mean <- rep(final.rate.of.decay.mean, total.exp.years)
        } else {
            if(is.null(names(final.rate.of.decay.mean)))
                stop('"names(projected.rate.of.decay[[\'final.rate.of.decay.mean\']])" cannot be NULL if supplied as a vector. See documentation of this function.')
            tmp.years <- names(final.rate.of.decay.mean)
            suppressWarnings(tmp.years <- as.numeric(tmp.years))
            if(any(is.na(tmp.years)) || any(as.integer(tmp.years) != tmp.years))
                stop('"names(projected.rate.of.decay[[\'final.rate.of.decay.mean\']])" must contain integer values stored as charcters. c(1990, 1991, ...)')
            if(is.unsorted(tmp.years) || any(tmp.years[-1] - tmp.years[-length(tmp.years)] != 1) || tmp.years[1] != ans@exposureYears[1])
                stop('"names(projected.rate.of.decay[[\'final.rate.of.decay.mean\']])" must contain consecutive integer values stored as charcters. c(1990, 1991, ...) and start at dimnames(incremental.payments)[[1]][1]')
            rm(tmp.years)

        }
        ans@finalRateOfDecayMean <- final.rate.of.decay.mean

        ##final.rate.of.decay.sd
        final.rate.of.decay.sd <- projected.rate.of.decay[['final.rate.of.decay.sd']]
        if(!is.vector(final.rate.of.decay.sd, mode='numeric') || (length(final.rate.of.decay.sd)!= 1 && length(final.rate.of.decay.sd) != total.exp.years))
            stop('"projected.rate.of.decay[[\'final.rate.of.decay\']]" must be a numeric vector of length "total.exp.years" or 1')
        if(any(is.infinite(final.rate.of.decay.sd)))
            stop('"projected.rate.of.decay[[\'final.rate.of.decay.sd\']]" must be finite')
        if(any(is.na(final.rate.of.decay.sd)))
            stop('"projected.rate.of.decay[[\'final.rate.of.decay.sd\']]" must not be NA')
        if(!is.numeric(final.rate.of.decay.sd))
            stop('"projected.rate.of.decay[[\final.rate.of.decay.sd\']]" must be numeric')
        if(any(final.rate.of.decay.sd <= 0))
            stop('"projected.rate.of.decay[[\'final.rate.of.decay.sd\']]" must be a > 0')

        if(length(final.rate.of.decay.sd) == 1)
        {
            final.rate.of.decay.sd <- rep(final.rate.of.decay.sd, total.exp.years)
        } else {
            if(is.null(names(final.rate.of.decay.sd)))
                stop('"names(projected.rate.of.decay[[\'final.rate.of.decay.sd\']])" cannot be NULL if supplied as a vector. See documentation of this function.')
            tmp.years <- names(final.rate.of.decay.sd)
            suppressWarnings(tmp.years <- as.numeric(tmp.years))
            if(any(is.na(tmp.years)) || any(as.integer(tmp.years) != tmp.years))
                stop('"names(projected.rate.of.decay[[\'final.rate.of.decay.sd\']])" must contain integer values stored as charcters. c(1990, 1991, ...)')
            if(is.unsorted(tmp.years) || any(tmp.years[-1] - tmp.years[-length(tmp.years)] != 1) || tmp.years[1] != ans@exposureYears[1])
                stop('"names(projected.rate.of.decay[[\'final.rate.of.decay.sd\']])" must contain consecutive integer values stored as charcters. c(1990, 1991, ...) and start at dimnames(incremental.payments)[[1]][1]')
            rm(tmp.years)

        }
        ans@finalRateOfDecaySD <- final.rate.of.decay.sd

        ##final.rate.of.decay.weight
        final.rate.of.decay.weight <- projected.rate.of.decay[['final.rate.of.decay.weight']]
        if((!is.vector(final.rate.of.decay.weight, mode='numeric') || length(final.rate.of.decay.weight) != total.dev.years - dim(incremental.payments)[2]) &&
           (!is.matrix(final.rate.of.decay.weight) || !all(dim(final.rate.of.decay.weight) == c(total.exp.years, total.dev.years - dim(incremental.payments)[2]))))
            stop('"projected.rate.of.decay[[\'final.rate.of.decay.weight\']]" must be a numeric vector of length "total.dev.years - dim(incremental.payments)[2]" or a matrix of dim total.exp.years, total.dev.years - dim(incremental.payments)[2]')
        if(any(is.infinite(final.rate.of.decay.weight)))
            stop('"projected.rate.of.decay[[\'final.rate.of.decay.weight\']]" must be finite')
        if(any(is.na(final.rate.of.decay.weight)))
            stop('"projected.rate.of.decay[[\'final.rate.of.decay.weight\']]" must not be NA')
        if(!is.numeric(final.rate.of.decay.weight))
            stop('"projected.rate.of.decay[[\final.rate.of.decay.weight\']]" must be numeric')


        if(is.vector(final.rate.of.decay.weight, mode='numeric'))
        {
            final.rate.of.decay.weight <- rep(final.rate.of.decay.weight, rep(total.exp.years, total.dev.years - dim(incremental.payments)[2]))
            dim(final.rate.of.decay.weight) <- c(total.exp.years, total.dev.years - dim(incremental.payments)[2])
        } else {
            if(is.null(dimnames(final.rate.of.decay.weight)) || is.null(dimnames(final.rate.of.decay.weight)[[1]]))
                stop('"dimnames(projected.rate.of.decay[[\'final.rate.of.decay.weight\']])[[1]]" cannot be NULL if supplied as a matrix. See documentation of this function.')
            tmp.years <- dimnames(final.rate.of.decay.weight)[[1]]
            suppressWarnings(tmp.years <- as.numeric(tmp.years))
            if(any(is.na(tmp.years)) || any(as.integer(tmp.years) != tmp.years))
                stop('"dimnames(projected.rate.of.decay[[\'final.rate.of.decay.weight\']])[[1]]" must contain integer values stored as charcters. c(1990, 1991, ...)')
            if(is.unsorted(tmp.years) || any(tmp.years[-1] - tmp.years[-length(tmp.years)] != 1) || tmp.years[1] != ans@exposureYears[1])
                stop('"dimnames(projected.rate.of.decay[[\'final.rate.of.decay.sd\']])[[1]]" must contain consecutive integer values stored as charcters. c(1990, 1991, ...) and start at dimnames(incremental.payments)[[1]][1]')
            rm(tmp.years)

        }
        ans@finalRateOfDecayWeight <- final.rate.of.decay.weight


    } else {
        stop('"projected.rate.of.decay" must either be a single NA or a list.  See documentation')
    }

    ans@triangleType <- match.arg(exp.year.type) #The type of triangle loaded.

    if(!is.numeric(prior.for.knot.locations) || length(prior.for.knot.locations) != 1)
        stop('"prior.for.knot.locations" must be a numeric of length 1')
    if(prior.for.knot.locations < 1)
        stop('"prior.for.knot.locations" must be at least 1')
    ans@priorForKnotPositions <- prior.for.knot.locations

    if(!is.numeric(prior.for.number.of.knots) || length(prior.for.number.of.knots) != 2)
          stop('"prior.for.number.of.knots" must be a numeric of length 2')
    if(prior.for.number.of.knots[1] <= 0)
        stop('"prior.for.number.of.knots[1]" must be greater than zero')
    if(prior.for.number.of.knots[2] <= 0 ||prior.for.number.of.knots[2] >= 1 )
        stop('"prior.for.number.of.knots[2]" must be greater than zero but less than one')
    ans@priorForNumberOfKnots <- prior.for.number.of.knots




    if(use.skew.t)
        ans@allowForSkew <- TRUE
    else
        ans@allowForSkew <- FALSE

    if(!is.numeric(bound.for.skewness.parameter) || length(bound.for.skewness.parameter) != 1 || !is.finite(bound.for.skewness.parameter) || bound.for.skewness.parameter == 0)
        stop('"bound.for.skewness.parameter" must be a numeric of length 1 and be finite and non-zero')

    if(bound.for.skewness.parameter < 0)
    {
        warning('"bound.for.skewness.parameter" was supplied as a negative value. The absolute value will be used.')
        bound.for.skewness.parameter <- abs(bound.for.skewness.parameter)
    }

    ans@skewnessParameterBounds <- bound.for.skewness.parameter * c(-1, 1)

    if(!is.numeric(last.column.with.scale.innovation) || !identical(length(last.column.with.scale.innovation), as.integer(1)))
        stop('"last.column.with.scale.innovation" must be a single numeric')
    if(as.integer(last.column.with.scale.innovation) != last.column.with.scale.innovation)
        stop('"last.column.with.scale.innovation" must be coercible to an integer')
    if(last.column.with.scale.innovation < 1 || last.column.with.scale.innovation > dim(incremental.payments)[2])
        stop('"last.column.with.scale.innovation" must be at least 1 and at most  dim(incremental.payments)[2]')

    ans@noChangeInScaleParameterAfterColumn <- as.integer(last.column.with.scale.innovation)


    if(use.ar1.in.calendar.year)
        ans@ar1InCalendarYearEffect <- TRUE
    else
        ans@ar1InCalendarYearEffect <- FALSE

    if(use.ar1.in.exposure.growth)
        ans@ar1InExposureGrowth <- TRUE
    else
        ans@ar1InExposureGrowth <- FALSE


    if(!validObject(ans))
        stop('Unable to create a valid input object.')

    return(invisible(ans))



}



##' A method to collect all the needed model input specific to the standard model. Intended for internal use only.
##'
##' There are currently two types of \code{AnnualAggLossDevModel}s (break and standard).  These models have many data elements in common.
##' This method appends only the elements specific to the standard model onto the list created by a call to \code{NextMethod()}.
##'
##' The following elements are appended:
##' \describe{
##'   \item{\code{x.0}}{A single value. The lower bound for the location of knots.}
##'   \item{\code{x.r}}{A single value. The upper bound for the location of knots.}
##'   \item{\code{beta.prior}}{A vector giving the prior for the location of knots.}
##'   \item{\code{mu.number.of.knots.prior}}{The priors for the mean of the number of knots.}
##'   \item{\code{number.of.knots.ubound}}{The upper bound for the number of knots.}
##' }
##' @name getJagsData,StandardAnnualLossDevModelInput-method
##' @param object An object of type \code{StandardAnnualAggLossDevModelInput} from which to collect the needed model input.
##' @return A named list of the specific model elements.  See details for more info.
##' @docType methods
setMethod(
          'getJagsData',
          signature(object='StandardAnnualAggLossDevModelInput'),
          function(object)
      {
          ans <- callNextMethod()
          ##return(ans)
          log.inc.index <- ans$log.inc.index



          ##calculate the last column with two non-missing values
          x.r <- function()
          {
              log.inc <- ans$log.inc

              for(i in dim(log.inc)[2]:1)
              {
                  tmp <- sum(!is.na(log.inc[,i]))

                  if(tmp >= 2)
                      return(i)

              }

              stop('"log.inc" does not have a column with at least two non missing values')

          }


          ans$x.0 <- 1
          ans$x.r <- x.r()
          ans$beta.prior <- c(1,object@priorForKnotPositions)

          ans$mu.number.of.knots.prior <- object@priorForNumberOfKnots
          ans$mu.number.of.knots.prior[2] <- (1 - ans$mu.number.of.knots.prior[2]) / ans$mu.number.of.knots.prior[2]
          ans$number.of.knots.ubound <- trunc(ans$K/2) + 1
          return(ans)



      })




##' A method to collect all the needed initial values unique to the standard model. Intended for internal use only.
##'
##' There are currently two types of \code{AnnualAggLossDevModel}s (break and standard). Code needed to create initial values specific the standard model is placed in this method.
##' This method returns a parameterless function which, when called, first calls the function returned by \code{NextMethod()} and then appends the following initial values onto the list returned by that function:
##'
##' \describe{
##'  \item{\code{S.}}{The initial values for the spline node.  Needed because \emph{dspline} cannot create initial values.}
##' }
##' @name getJagsInits,StandardAnnualLossDevModelInput-method
##' @param object An object of type \code{StandardAnnualAggLossDevModelInput} from which to collect the needed initial values for the model.
##' @return A named list of the specific model elements.  See details for more information.
##' @docType methods
##' @seealso \code{\link{getJagsInits}}
setMethod(
          'getJagsInits',
          signature(object='StandardAnnualAggLossDevModelInput'),
          function(object)
      {
          super.f <- callNextMethod()
          K <- getTriDim(object)[1]
          S. <- array(0,c(K,4))
          S.[,2] <- 1
          function()
          {
              ans <- super.f()
              ans$S. <- S.
              return(ans)
          }
      })
