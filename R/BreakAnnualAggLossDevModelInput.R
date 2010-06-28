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

##This is intended to be a final class for models with a break.
##This class is derived from the AnnualAggLosDevModelInput class.
##' @include zzz.R
##' @include AnnualAggLossDevModelInput.R
NULL

##' The final input class for models with a break.
##'
##' \code{BreakAnnualAggLossDevModelInput} is the final input class for models with a break.
##' @name BreakAnnualAggLossDevModelInput-class
##' @docType class
##' @seealso \code{\linkS4class{LossDevModelInput}}
##'  \code{\linkS4class{StandardAnnualAggLossDevModelInput}}
setClass(
         'BreakAnnualAggLossDevModelInput',
         representation(
                        rangeForFirstYearInNewRegime='integer',            #The range of exposure years describing the break.
                        priorsForFirstYearInNewRegime='numeric',           #The prior for the break.
                        priorForKnotPositionsPreBreak='numeric',
                        priorForKnotPositionsPostBreak='numeric',
                        priorForNumberOfKnotsPreBreak='numeric',
                        priorForNumberOfKnotsPostBreak='numeric'),
         contains='AnnualAggLossDevModelInput')



##' Create an Object of class \code{BreakAnnualAggLossDevModelInput}.
##'
##' The loss development models require a lot of input.  Much of the input is directly dependent on the values of other input.
##' As such, this function facilitates much of the work of setting model parameters and determining which output to collect.
##'
##' The function creates an object of class \code{BreakAnnualAggLossDevModelInput}.
##'
##' Many arguments and much functionality is described in \code{\link{makeStandardAnnualInput}}.
##'
##' \describe{
##'   \item{\code{first.year.in.new.regime}}{
##'     The break model allows for a structural break along the exposure year axis in the consumption path.
##'     The slope of the pre-break consumption path is used to extend the post-break consumption path.
##'     The time when the break occurs can either be specified with 100% certainty by the user, or it can be estimated by the model.
##'     To specify \emph{the} first exposure year in which a new consumption path applies, the user should supply a single value to the parameter \code{first.year.in.new.regime}.
##'     To have the model estimate the point in time when the break occurs, the user should supply a range (min and max) for the possible values to \code{first.year.in.new.regime}.
##'     Note also that the number of rows above and below the break must be at least 4 and as such a triangle of size 8 is the smallest triangle which can be estimated with a break.
##'   }
##'   \item{\code{prior.for.first.year.in.new.regime}}{
##'     The prior for the \code{first.year.in.new.regime} is a (scaled and discretized) beta distribtion.
##'     A value of \code{c(1,1)} would indicate a uniform distribution.
##'     See \code{\link{firstYearInNewRegime}}.
##'   }
##'   \item{\code{prior.for.knot.locations.pre.break} and \code{prior.for.knot.locations.post.break}}{
##'     If these values are \code{NA} (the default), then \code{prior.for.knot.locations.pre.break} will be assigned a value of 2.
##'     And \code{prior.for.knot.locations.post.break} wil be assigned a value of \code{1 + (num.years.in.post.break.period + 0.5 * num.years.in.break.period)/num.years.in.triangle.}.
##'     These values must either both be \code{NA} or both be set to numbers.
##'   }
##' }
##'
##' @param incremental.payments A square matrix of incremental payments.  Row names should correspond to the exposure year. Only upper-left (including the diagonal) of Triangle may have non-missing values.  Lower-right must be \code{NA}.
##' @param first.year.in.new.regime May be one of two types.  1) A single numeric value.  2) A numeric vector of length 2. See Details.
##' @param prior.for.first.year.in.new.regime A numeric vector of length 2. See Details.
##'
##' @param extra.exp.years A single integer value (\code{total.exp.years} overrides) greater than or equal to 1 (default is 1) specifying the number of additional exposure years (or rows in the triangle) to project.
##' @param extra.dev.years A single integer value (\code{total.dev.years} overrides) greater than or equal to 1 (default is 1) specifying the additional number of development years (or columns in the triangle) to project.
##'
##' @param non.stoch.inflation.rate May be one of three types (See \emph{Inflation Rate} in Details). 1) A single numeric value. 2) A vector of numerics (of specific length). 3) A matrix of numerics (of specific dim).
##' @param non.stoch.inflation.weight May be one of three types (See \emph{Inflation Rate} in Details). 1) A single numeric value. 2) A vector of numerics (of specific length). 3) A matrix of numerics (of specific dim).
##'
##' @param stoch.inflation.rate May be one of two types (See \emph{Inflation Rate} in Details). 1) A single numeric value of \emph{zero}. 2) A vector of numerics (of specific length).
##' @param stoch.inflation.weight May be one of three types (See \emph{Inflation Rate} in Details). 1) A single numerical value. 2) A vector of numerics (of specific length). 3) A matrix of numerics (of specific dim).
##' @param stoch.inflation.lower.bound  May be one of three types (See \emph{Inflation Rate} in Details). 1) A single numeric value. 2) A vector of numerics (of specific length). 3) A matrix of numerics (of specific dim).
##' @param stoch.inflation.upper.bound  May be one of three types (See \emph{Inflation Rate} in Details). 1) A single numeric value. 2) A vector of numerics (of specific length). 3) A matrix of numerics (of specific dim).
##' @param known.stoch.inflation.mean   May be one of two types (See \emph{Inflation Rate} in Details). 1) A single numeric value. 2) \code{NA}.
##' @param known.stoch.inflation.persistence  May be one of two types (See \emph{Inflation Rate} in Details). 1) A single numeric value. 2) \code{NA}.
##'
##' @param total.exp.years A single integer value (overrides \code{extra.exp.years}) specifying the last exposure year to project.  Must be at least the number of rows in \code{incremental.payments} + 1.
##' @param total.dev.years A single integer value (overrides \code{extra.dev.years}) specifying the last development year to project.  Must be at least the number of columns in \code{incremental.payments} + 1 .
##'
##' @param cumulative.payments A numeric matrix with the same dim and dimnames as \code{incremental.payments}.  Must be a possible cumulative payment triangle of \code{incremental.payments}.  (See \emph{Cumulative Payments} Section).
##'
##' @param exp.year.type A single character value indicating the type of exposure years:  \sQuote{ambiguous}, \sQuote{py}, and \sQuote{ay} mean \sQuote{Exposure Year}, \sQuote{Policy Year}, and \sQuote{Accident Year}; respectively.
##' @param prior.for.knot.locations.pre.break A single numeric value of at least 1.  The prior for the location of knots is a scaled beta with parameters \code{c(1,prior.for.knot.locations.pre.break)}.
##' @param prior.for.number.of.knots.pre.break A two element vector giving the paramters for the prior number of knots.
##' @param prior.for.knot.locations.post.break See \code{prior.for.knot.locations.pre.break}. Large values produce stable consumption paths at high development years.
##' @param prior.for.number.of.knots.post.break A two element vector giving the paramters for the prior number of knots.
##'
##' @param use.skew.t A logical value.  If \code{TRUE} the model assumes the observed and estimated log incremental payments are realizations from a skewed \eqn{t} distribution; if \code{FALSE} it will assume zero skewness. (See Reference.)
##' @param bound.for.skewness.parameter A positive numerical value representing the symetric boundaries for the skewness parameter.  In most cases, the default should be sufficient. Ignored if \code{use.skew.t=FALSE}.
##' @param last.column.with.scale.innovation A single integer must be at least 1 and at most the number of columns in \code{incremental.payments}.  See \emph{Measurment Error-Second Order Random Walk} in Details.
##' @param use.ar1.in.calendar.year A logical value.  The calendar year effect errors may (at the users digression) include an autoregressive process of order 1.  \code{TRUE} turns on the ar1 process, \code{FALSE} turns it off.
##' @param use.ar1.in.exposure.growth A logical value.  The exposure growth errors may (at the users discretion) include an autoregressive process of order 1.  \code{TRUE} (the Default) turns on the ar1 process, \code{FALSE} turns it off.
##'
##' @param projected.rate.of.decay  May be one of three types (See \emph{Projected Rate of Decay} in Details). 1) \code{NA}. 2) A matrix of numerics (of specific dim). 3) A named list.

##'

##'
##' @references
##' Kim, Y., and J. McCulloch (2007) \dQuote{The Skew-Student Distribution with Application to U.S. Stock Market Returns and the Equity Premium,} Department of Economics, Ohio State University, October 2007
##'
##'
##' @return An object of class \code{AggModelInput}.  This the model specified by the returned object must then be estimated using the function \code{runLossDevModel}.
##' @export
##' @usage
##' makeBreakAnnualInput(
##'   incremental.payments=decumulate(cumulative.payments),
##'   first.year.in.new.regime=trunc(median(as.integer(dimnames(incremental.payments)[[1]]))),
##'   prior.for.first.year.in.new.regime=c(2,2),
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
##'   prior.for.knot.locations.pre.break=NA,
##'   prior.for.number.of.knots.pre.break=c(3, 1/7),
##'   prior.for.knot.locations.post.break=NA,
##'   prior.for.number.of.knots.post.break=c(3, 1/7),
##'   use.skew.t=FALSE,
##'   bound.for.skewness.parameter=10,
##'   last.column.with.scale.innovation=dim(incremental.payments)[2],
##'   use.ar1.in.calendar.year=FALSE,
##'   use.ar1.in.exposure.growth=TRUE,
##'   projected.rate.of.decay=NA)
makeBreakAnnualInput <- function(incremental.payments=decumulate(cumulative.payments),
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
                                 prior.for.number.of.knots.pre.break=c(3,1/7),
                                 prior.for.knot.locations.post.break=NA,
                                 prior.for.number.of.knots.post.break=c(3,1/7),
                                 use.skew.t=FALSE,
                                 bound.for.skewness.parameter=10,
                                 last.column.with.scale.innovation=dim(incremental.payments)[2],
                                 use.ar1.in.calendar.year=FALSE,
                                 use.ar1.in.exposure.growth=TRUE,
                                 projected.rate.of.decay=NA)
{

    standard.ans <- makeStandardAnnualInput(incremental.payments=incremental.payments,
                                            extra.dev.years=extra.dev.years,
                                            extra.exp.years=extra.exp.years,
                                            non.stoch.inflation.rate=non.stoch.inflation.rate,
                                            non.stoch.inflation.weight=non.stoch.inflation.weight,
                                            stoch.inflation.rate=stoch.inflation.rate,
                                            stoch.inflation.weight=stoch.inflation.weight,
                                            stoch.inflation.lower.bound=stoch.inflation.lower.bound,
                                            stoch.inflation.upper.bound=stoch.inflation.upper.bound,
                                            known.stoch.inflation.mean=known.stoch.inflation.mean,
                                            known.stoch.inflation.persistence=known.stoch.inflation.persistence,
                                            total.dev.years=total.dev.years,
                                            total.exp.years=total.exp.years,
                                            cumulative.payments=cumulative.payments,
                                            exp.year.type=exp.year.type,

                                            use.skew.t=use.skew.t,
                                            bound.for.skewness.parameter=bound.for.skewness.parameter,
                                            last.column.with.scale.innovation=last.column.with.scale.innovation,
                                            use.ar1.in.calendar.year=use.ar1.in.calendar.year,
                                            use.ar1.in.exposure.growth=use.ar1.in.exposure.growth,
                                            projected.rate.of.decay=projected.rate.of.decay)

    ans <- new('BreakAnnualAggLossDevModelInput')

    for(s in slotNames(standard.ans))
    {
        if(s %in% slotNames(ans))
        {
            slot(ans, s) <-  slot(standard.ans, s)
        }
    }

    rm(standard.ans)




    ##perform some constancy checks
    if(getTriDim(ans)[1] < 8)
        stop('The triangle supplied to "incremental.payments" is too small to run with a break.  Try running the standard model.')

    suppressWarnings(int.first.year.in.new.regime <- as.integer(first.year.in.new.regime))

    if(any(int.first.year.in.new.regime != first.year.in.new.regime))
        stop('"first.year.in.new.regime" must be coercible to an integer')

    if(length(int.first.year.in.new.regime) > 2 || length(int.first.year.in.new.regime) == 0  || !all(int.first.year.in.new.regime %in% ans@exposureYears))
        stop('"first.year.in.new.regime" must have a length of either 1 or 2 and must only contain values in the exposure years of incremental.payments."')

    if(length(int.first.year.in.new.regime) == 1)
    {
        ans@rangeForFirstYearInNewRegime <- rep(int.first.year.in.new.regime, 2)
    } else {
        if(int.first.year.in.new.regime[1] > int.first.year.in.new.regime[2])
        {
            warning('"int.first.year.in.new.regime[1] > int.first.year.in.new.regime[2]" these values will be resorted')
            int.first.year.in.new.regime <- sort(int.first.year.in.new.regime)
        }
        ans@rangeForFirstYearInNewRegime <- int.first.year.in.new.regime
    }

    if(length(prior.for.knot.locations.pre.break) == 1 && is.na(prior.for.knot.locations.pre.break) &&
       length(prior.for.knot.locations.post.break) == 1 && is.na(prior.for.knot.locations.post.break))
    {
        ans@priorForKnotPositionsPreBreak <- 2
        ans@priorForKnotPositionsPostBreak <- 1 +
            (sum(ans@exposureYears > max(ans@rangeForFirstYearInNewRegime)) + .5 * (1 + diff(ans@rangeForFirstYearInNewRegime))) /
                length(ans@exposureYears)

    } else {

        if(!is.numeric(prior.for.knot.locations.pre.break) || length(prior.for.knot.locations.pre.break) != 1)
            stop('"prior.for.knot.locations.pre.break" and "prior.for.knot.locations.post.break" must both be numeric of length 1.')
        if(prior.for.knot.locations.pre.break < 1)
            stop('"prior.for.knot.locations.pre.break" must be at least 1')
        ans@priorForKnotPositionsPreBreak <- prior.for.knot.locations.pre.break



        if(!is.numeric(prior.for.knot.locations.post.break) || length(prior.for.knot.locations.post.break) != 1)
            stop('"prior.for.knot.locations.pre.break" and "prior.for.knot.locations.post.break" must both be numeric of length 1')
        if(prior.for.knot.locations.post.break < 1)
            stop('"prior.for.knot.locations.post.break" must be at least 1')

         ans@priorForKnotPositionsPostBreak <- prior.for.knot.locations.post.break
    }


   if(!is.numeric(prior.for.number.of.knots.pre.break) || length(prior.for.number.of.knots.pre.break) != 2)
          stop('"prior.for.number.of.knots.pre.break" must be a numeric of length 2')
    if(prior.for.number.of.knots.pre.break[1] <= 0)
        stop('"prior.for.number.of.knots.pre.break[1]" must be greater than zero')
    if(prior.for.number.of.knots.pre.break[2] <= 0 ||prior.for.number.of.knots.pre.break[2] >= 1 )
        stop('"prior.for.number.of.knots.pre.break[2]" must be greater than zero but less than one')
    ans@priorForNumberOfKnotsPreBreak <- prior.for.number.of.knots.pre.break

    if(!is.numeric(prior.for.number.of.knots.post.break) || length(prior.for.number.of.knots.post.break) != 2)
          stop('"prior.for.number.of.knots.post.break" must be a numeric of length 2')
    if(prior.for.number.of.knots.post.break[1] <= 0)
        stop('"prior.for.number.of.knots.post.break[1]" must be greater than zero')
    if(prior.for.number.of.knots.post.break[2] <= 0 ||prior.for.number.of.knots.post.break[2] >= 1 )
        stop('"prior.for.number.of.knots.post.break[2]" must be greater than zero but less than one')
    ans@priorForNumberOfKnotsPostBreak <- prior.for.number.of.knots.post.break

    if(sum(ans@exposureYears < ans@rangeForFirstYearInNewRegime[1]) < 4)
        stop('The minimum value for "first.year.in.new.regime" is too small.  There must be at least 4 years in the pre-break period.')

    if(sum(ans@exposureYears >= ans@rangeForFirstYearInNewRegime[2]) < 4)
        stop('The maximum value for "first.year.in.new.regime" is too large.  There must be at least 4 years in the post-break period.')

    if(!is.numeric(prior.for.first.year.in.new.regime) || length(prior.for.first.year.in.new.regime) != 2)
        stop('The value supplied for "prior.for.first.year.in.new.regime" must be a numeric vector of length 2.')

    if(any(prior.for.first.year.in.new.regime <=0))
        stop('"prior.for.first.year.in.new.regime" must be greater than 0')

    ans@priorsForFirstYearInNewRegime <- prior.for.first.year.in.new.regime


    ans@outputType <- 'BreakAnnualAggLossDevModelOutput'
    ans@modelFile <- 'break.model.txt'
    return(invisible(ans))
}



##' A method to collect all the needed model input specific to the break model.
##'
##' There are currently two types of \code{AnnualAggLossDevModel}s (break and standard).  These models have many data elements in common.
##' This method appends only the elements specific to the break model onto the list created by a call to \code{NextMethod()}.
##'
##' The following elements are appended:
##' \describe{
##'   \item{\code{x.0}}{Two values. The lower bound for the number of knots. First is for first spline.}
##'   \item{\code{x.r}}{Two values. The upper bound for the number of knots. First is for first spline.}
##'   \item{\code{break.row}}{Two integer values given the (inclusive) range of rows in which the first year in the new regime could occur.}
##'   \item{\code{break.row.priors}}{The parameters for the beta distribution which serves as the prior for the location of the structural break.}
##'   \item{\code{K.trim}}{The maximum number of columns in the post-break triangle.}
##'   \item{\code{beta.prior}}{A matrix giving the prior for the location of knots.  First column is for the pre-break spline.  Second is for the post-break spline.}
##'   \item{\code{mu.number.of.knots.prior}}{A matrix giving the prior for the mean of the number of knots.  First column is for the pre-break spline.  Second is for the post-break spline}
##'   \item{\code{number.of.knots.ubound}}{A vector giving the upper bound for the number of knots.  First is for the pre-break spline.  Second is for the post-break spline}
##' }
##' @name getJagsData,BreakAnnualLossDevModelInput-method
##' @param object An object of type \code{BreakAnnualAggLossDevModelInput} from which to collect the needed model input.
##' @return A named list of the specific model elements.  See details for more info.
##' @docType methods
setMethod(
          'getJagsData',
          signature(object='BreakAnnualAggLossDevModelInput'),
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

          x.r.post.break <- function()
          {
              log.inc <- ans$log.inc
              log.inc <- log.inc[object@exposureYears >= max(object@rangeForFirstYearInNewRegime),]
              for(i in dim(log.inc)[2]:1)
              {
                  tmp <- sum(!is.na(log.inc[,i]))

                  if(tmp >= 2)
                      return(i)

              }
              stop('"log.inc" for the post period does not have a column with at least two non missing values')
          }

          ans$x.0 <- c(1, 1)
          ans$x.r <- c(x.r(), x.r.post.break())
          ans$break.row <- match(object@rangeForFirstYearInNewRegime, object@exposureYears)
          ans$break.row.priors <- object@priorsForFirstYearInNewRegime
          ans$K.trim <- ans$K - ans$break.row[1] + 1

          ans$beta.prior <- array(NA, c(2,2))
          ans$beta.prior[,1] <- c(1,object@priorForKnotPositionsPreBreak)
          ans$beta.prior[,2] <- c(1,object@priorForKnotPositionsPostBreak)

          ans$mu.number.of.knots.prior <- array(NA, c(2, 2))
          ans$mu.number.of.knots.prior[,1] <- object@priorForNumberOfKnotsPreBreak
          ans$mu.number.of.knots.prior[2,1] <- (1 - ans$mu.number.of.knots.prior[2,1]) / ans$mu.number.of.knots.prior[2,1]
          ans$mu.number.of.knots.prior[,2] <- object@priorForNumberOfKnotsPostBreak
          ans$mu.number.of.knots.prior[2,2] <- (1 - ans$mu.number.of.knots.prior[2,2]) / ans$mu.number.of.knots.prior[2,2]
          ans$number.of.knots.ubound <- numeric(2)
          ans$number.of.knots.ubound[1] <- trunc(ans$K/2) + 1
          ans$number.of.knots.ubound[2] <- trunc(ans$K.trim/2) + 1

          return(ans)
      })




##' A method to collect all the needed model initial values unique to the break model.
##'
##' There are currently two types of \code{AnnualAggLossDevModel}s (break and standard). Code needed to create initial values specific the break model is placed in this method.
##' This method returns a parameterless function which when called first calls the function returned by \code{NextMethod()} and then to the list returned by that function appends the following initial values.
##'
##' \describe{
##'  \item{\code{R.}}{The initial values for the spline node.  Needed because dspline cannot create initial values.}
##' }
##' @name getJagsInits,BreakAnnualLossDevModelInput-method
##' @param object An object of type \code{BreakAnnualAggLossDevModelInput} from which to collect the needed model initial values.
##' @return A named list of the specific model elements.  See details for more info.
##' @docType methods
##' @seealso \code{\link{getJagsInits}}
setMethod(
          'getJagsInits',
          signature(object='BreakAnnualAggLossDevModelInput'),
          function(object)
      {
          super.f <- callNextMethod()
          K <- getTriDim(object)[1]
          R. <- array(0, c(K, 6))
          R.[,2] <- 1
          function()
          {
              ans <- super.f()
              ans$R. <- R.
              return(ans)
          }
      })
