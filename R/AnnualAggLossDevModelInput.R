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


##All Input classes for Annual Aggregate models should be derived from this base class.
##This class is derived from the LossDevModelInput class.
##' @include zzz.R
##' @include LossDevModelInput.R
NULL

##' The base class for all aggregate models.
##'
##' \code{AnnualAggLossDevModelInput} is the base class for all aggregate model input objects.  It is derived from \code{LossDevModelInput}.
##' @name AnnualAggLossDevModelInput-class
##' @docType class
##' @seealso \code{\linkS4class{LossDevModelInput}}, \code{\linkS4class{StandardAnnualAggLossDevModelInput}}
setClass(
         'AnnualAggLossDevModelInput',
         representation(
                        triangleType='character',                 #The type of triangle loaded.
                        totalDevYears='integer',                  #The total number of columns including the projection period.
                        totalExpYears='integer',                  #The total number of rows including the projection period.
                        incrementals='matrix',                    #the incremental triangle, must be square
                        cumulatives='matrix',                     #the cumulative triangle, must be square
                        exposureYears='integer',                  #a vector of years which correspond to the rows in the triangel
                        stochInflationRate='numeric',             #The normal inflation rate.  (x[i] / x[i-1] - 1)
                        knownStochInflationMean='numeric',        #The normal inflation rate.  (x[i] / x[i-1] - 1)
                        knownStochInflationPersistence='numeric',
                        stochInflationWeight='matrix',            #The normal inflation rate.  (x[i] / x[i-1] - 1)
                        stochInflationUpperBound='matrix',        #The normal inflation rate.  (x[i] / x[i-1] - 1)
                        stochInflationLowerBound='matrix',        #The normal inflation rate.  (x[i] / x[i-1] - 1)
                        stochInflationYears='integer',            #The corresponding years of elements in inflationRate and logInflationRate.
                        nonStochInflationRate='matrix',           #The normal inflation rate.  (x[i] / x[i-1] - 1)
                        nonStochInflationWeight='matrix',         #The normal inflation rate.  (x[i] / x[i-1] - 1)
                        lastNonZeroPayment='integer',
                        finalRateOfDecayMean='numeric',
                        finalRateOfDecaySD='numeric',
                        finalRateOfDecayWeight='matrix',
                        allowForSkew='logical',                   #should the skewed-t be turned on?
                        skewnessParameterBounds='numeric',
                        ar1InCalendarYearEffect='logical',        #should the calendar year effect include an ar1 or white-noise error term?
                        ar1InExposureGrowth='logical',            #should the exposure growth include an ar1 or white-noise error term?
                        noChangeInScaleParameterAfterColumn='integer', #valid values are 1 through K, 1 means all columns have same scale, K means all have different, acutal value will be truncated to last observed column
                        'VIRTUAL'),
         contains='LossDevModelInput')

##' A function to return a "nice" character string for the exposure year label in charts.
##'
##' @param object An object of type \code{AnnualAggLossDevModelInput}.
##' @return A character.
getExposureYearLabel <- function(object)
{
    stopifnot(is(object, 'AnnualAggLossDevModelInput'))
    if(length(object@triangleType) != 1)
        stop('length(object@triangleType) must be 1')

    i <- match(object@triangleType, c('ambiguous', 'py', 'ay'))
    if(is.na(i))
        stop('object@triangleType must be one of "ambiguous", "py", or "ay"')

    return(c('Exposure Year', 'Policy Year', 'Accident Year')[i])
}

##' A function to return the index non-missing values of a 2d container.
##'
##' @param m A 2d container on which \code{dim} can be called and is subsetable.
##' @return A 2d array of integers.  Each rows corresponds to a non-missing value in \code{m}.  First column is the row.  Second column is the column.
##' @keywords internal
getNonMissingIndexes <- function(m)
{
    stopifnot(
              !is.null(dim(m)),
              length(dim(m)) == 2)

    first.non.missing <- TRUE
    for(i in 1:dim(m)[1])
    {
        for(j in 1:dim(m)[2])
        {
            if(!is.na(m[i,j]))
            {
                if(first.non.missing)
                {
                    ans <- c(i,j)
                    first.non.missing <- FALSE
                } else {

                    ans <- rbind(ans, c(i,j))
                }
            }
        }
    }
    return(ans)
}


##' A generic function to get the dimension of the observed triangle.  Intended for internal use only.
##'
##' The dimension of a triangle is the number of rows and number of columns.
##' @name getTriDim
##' @param object The object from which to get the dim of the observed triangle.
##' @return An integer vector of the dimension of the observed triangle.
##' @docType genericFunction
setGenericVerif('getTriDim',
                 function(object)
                 standardGeneric('getTriDim'))

##' A method to get the dimension of the observed triangle.  Intended for internal use only.
##'
##' @name getTriDim,AnnualAggLossDevModelInput-method
##' @param object An object of type \code{AnnualLossDevModelInput}.
##' @return An integer vector of the dimension of the observed triangle.
##' @docType methods
##' @seealso \code{\link{getTriDim}}
setMethod('getTriDim',
          signature(object='AnnualAggLossDevModelInput'),
          function(object)
          {
            return(dim(object@incrementals))
          })

##' A method to create all the needed JAGS input common to both the standard model and the break model. Intended for internal use only.
##'
##' There are currenly two types of \code{AnnualAggLossDevModel}s (break and standard).  These models have many data elements in common and the code to create these common elements is placed in this method.
##' The derived classes \code{StandardAnnualAggLossDevModelInput} and \code{BreakAnnualAggLossDevModelInput} should call this method via \code{NextMethod()} and then should append any needed data.
##'
##' The returned value is a list containing the following elements with the following meanings:
##' \describe{
##'   \item{\code{log.inc}}{A square matrix with only the upper right containing non-missing values of the log of the incremental payments.}
##'   \item{\code{log.inc.index}}{A 2-column matrix.  Each row corresponds to a non-missing value in \code{log.inc}.  The first column gives the row of the non-missing value, the second the column.}
##'   \item{\code{log.inc.index.length}}{The number of rows in \code{log.inc.index}.}
##'   \item{\code{smooth.tau.h.2.log}}{A 2-valued vector giving the smoothing parameters for the changing (in development time) variance and skewness.}
##'   \item{\code{h.same.as.last}}{A vector of 1's and 0's equal in length to the columns in \code{log.inc}.  A value of one means the variance and skewness in the corresponding column should be the same as in the previous column.}
##'   \item{\code{K}}{Dimension of \code{log.inc}.}
##'   \item{\code{H}}{Number of additional rows to estimate.}
##'   \item{\code{L.vec}}{Vector equal in length to the number of total rows (observed plus forecast).  Each value is the column number of  the final non-zero incremental payment.}
##'   \item{\code{delta.p}}{Vector equal in length to the number of total rows (observed plus forecast).  First column is (on the log scale) the mean of the final decay rate.  Second is standard deviation.}
##'   \item{\code{v}}{Matrix with number of rows equal to total rows (observed plus forecast) and number of columns equal to total development years (observed plus forecast) minus the number of observed development years.}
##'   \item{\code{ar1}}{Zero if calendar year effect does not have an ar1 component; 1 if it does.}
##'   \item{\code{rho.prior}}{The parameters for the beta prior of the autoregressive parameter in calendar year effect.  Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
##'   \item{\code{ar1.eta}}{Zero if exposure growth does not have an ar1 component; 1 if it does.}
##'   \item{\code{rho.eta.prior}}{The parameters for the beta prior of the autoregressive parameter in exposure growth.  Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
##'   \item{\code{allow.for.skew}}{Zero tells the model to assume zero skewness, one to use the skewed \eqn{t}.}
##'   \item{\code{precision.for.skewness}}{The prior precision for the skewness parameter.  Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
##'   \item{\code{df.for.skewness}}{The prior df for the skewness parameter.  Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
##'   \item{\code{bounds.for.skewness}}{The the skewness parameter is bounted to prevent chains from getting stuck.  Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
##'   \item{\code{precision.for.eta.mu}}{The prior precision for \code{eta.mu}.  Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
##'   \item{\code{df.k}}{The parameter for the chi-sqr prior on the degrees of freedom.  Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
##'   \item{\code{df.bounds}}{The lower and upper bounds on the degrees of freedom.  Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
##'   \item{\code{a.ou.prior}}{The parameters for the beta prior of the autoregressive parameter for the stochastic inflation. Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
##'   \item{\code{precision.for.b.ou}}{The prior precision for the intercept in the stochastic inflation ar1 process.  Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
##'   \item{\code{sigma.eta.bounds}}{The lower and upper bounds on the standard deviation of the exposure year growth rate.  Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
##'   \item{\code{sigma.kappa.bounds}}{The lower and upper bounds on the standard deviation of the calendar year growth rate. Stored in R rather than in the model so that modifications are accurately reflected in charts for the prior.}
##'   \item{\code{stoch.log.inf}}{A vector giving the log inflation rate; missing future values are coded as NA.}
##'   \item{\code{w.stoch.pct.inf.}}{Matrix with number of rows equal to total (observed and forecast) exposure years and number of columns equal to total development years. The portion of dollars inflating.}
##'   \item{\code{non.stoch.log.inf}}{Matrix with number of rows equal to total (observed and forecast) exposure years and number of columns equal total development years. Non-stochastic inflation rate.  Cannot be \code{NA}.}
##'   \item{\code{w.non.stoch.pct.inf}}{Matrix with number of rows equal to total (observed and forecast) exposure years and number of columns equal to total development years. The percent of dollars inflating.}
##'   \item{\code{P}}{Single value.  The position in \code{stoch.log.inf.c} corresponding to the final observed diagonal in \code{log.inc}.}
##'   \item{\code{stoch.log.inf.lower.bound}}{Matrix with number of rows equal to total (observed and forecast) exposure years and number of columns equal total development years. Floor for inflation rate.}
##'   \item{\code{stoch.log.inf.upper.bound}}{Matrix with number of rows equal to total (observed and forecast) exposure years and number of columns equal total development years. Ceiling for inflation rate.}
##'   \item{\code{estimate.a.ou}}{Single value (zero or one). Should the auto correlation coefficient for the Ornstein-Uhlenbeck process be estimated or should \code{fixed.a.ou} be used?}
##'   \item{\code{fixed.a.ou}}{Single value. Possible non-stochastic  auto-correlation coefficient for the Ornstein-Uhlenbeck process.}
##'   \item{\code{estimate.b.ou}}{Single value (zero or one). Should the stochastic rate of inflation have an estimated constant term or should \code{fixed.b.ou} be used?}
##'   \item{\code{fixed.b.ou}}{Single value. Possible non-stochastic constant term for stochastic inflation rate.}
##'   \item{\code{stoch.log.inf.known.mu}}{Single value.  Added to the log stochastic inflation rate after the ar1 estimation process.}
##' }
##' @name getJagsData,AnnualLossDevModelInput-method
##' @param object An object of type \code{AnnualAggLossDevModelInput} from which to collect the needed model input.
##' @return A named list of the specific model elements.  See Details for more information.
##' @docType methods
setMethod(
          'getJagsData',
          signature(object='AnnualAggLossDevModelInput'),
          function(object)
      {
          ans <- list()

          ans$log.inc <- object@incrementals
          ans$log.inc[ans$log.inc <= 0] <- NA
          ans$log.inc <- log(ans$log.inc)

          ans$log.inc.index <- getNonMissingIndexes(ans$log.inc)
          ans$log.inc.index.length <- dim(ans$log.inc.index)[1]

          ans$K <- getTriDim(object)[1]
          ans$H <- object@totalExpYears - ans$K
          ans$L.vec <- object@lastNonZeroPayment

          ans$smooth.tau.h.2.log.innov <- c(5, 0.5)

          supplied.last.column.with.different.scale <- object@noChangeInScaleParameterAfterColumn
          max.last.column.with.different.scale <- max(ans$log.inc.index[,2])
          last.column.with.different.scale <- max(1,
                                                  min(supplied.last.column.with.different.scale,
                                                      max.last.column.with.different.scale))


          ans$h.same.as.last <- c(rep(0, last.column.with.different.scale),
                                                rep(1, ans$K - last.column.with.different.scale))

          ans$delta.p <- cbind(log(1 + object@finalRateOfDecayMean), object@finalRateOfDecaySD)
          ans$v <- object@finalRateOfDecayWeight

          ans$ar1 <- ifelse(object@ar1InCalendarYearEffect, 1, 0)
          ans$rho.prior <- c(1, 1.5)

          ans$ar1.eta <- ifelse(object@ar1InExposureGrowth, 1, 0)
          ans$rho.eta.prior <- c(1, 1.5)


          ans$allow.for.skew <- ifelse(object@allowForSkew, 1, 0)
          ans$precision.for.skewness <- 0.5
          ans$df.for.skewness <- 2.1
          ans$bounds.for.skewness <- object@skewnessParameterBounds

          ans$precision.for.eta.mu <- 1

          ans$df.k <- 8

          if(object@allowForSkew)
              ans$df.bounds <- c(4,50)
          else
              ans$df.bounds <- c(2,50)


          ans$a.ou.prior <- c(1.1, 1.1)
          ans$precision.for.b.ou <- 1.0E-3

          ans$sigma.eta.bounds <- c(0,10)
          ans$sigma.kappa.bounds <- c(0,10)

          ans$non.stoch.log.inf <- log(1 + object@nonStochInflationRate)
          ans$w.non.stoch.pct.inf <-  object@nonStochInflationWeight

          if(identical(object@stochInflationRate, 0))
          {
              ans$P <- getTriDim(object)[1] - 1
              N <- ans$H + max(ans$L.vec) - 1
              ans$stoch.log.inf.c <- c(0, rep(NA, ans$P - 1), rep(NA, N))

              ans$w.stoch.pct.inf <- array(0, c(object@totalExpYears, object@totalDevYears))
              ans$stoch.log.inf.upper.bound <- array(3000, c(object@totalExpYears, object@totalDevYears))
              ans$stoch.log.inf.lower.bound <- array(-3000, c(object@totalExpYears, object@totalDevYears))

              ans$estimate.a.ou <- 1
              ans$fixed.a.ou <- 0

              ans$estimate.b.ou <- 1
              ans$fixed.b.ou <- 0
              ans$stoch.log.inf.known.mu <- 0

          } else {
              ans$stoch.log.inf.c <- log(1 + object@stochInflationRate)
              ##we are guaranteed that stochInflationRate will cover exposureYears
              ##we just have to caclulate P
              ##we don't have to worry if stochInflationRate is too long because the model file will take care of that
              ans$P <- match(max(object@exposureYears), object@stochInflationYears)
              ##number of diagionals beyond the last observed diagional (K+1) for which to simulate rates of inflation
              N <- ans$H + max(ans$L.vec) - 1
              times <- N - (length(ans$stoch.log.inf.c) - ans$P)
              if(times > 0)
                  ans$stoch.log.inf.c <- c(ans$stoch.log.inf.c, rep(NA, times))

              ans$w.stoch.pct.inf <- object@stochInflationWeight

              j.inf <- 3000
              suppressWarnings(ans$stoch.log.inf.upper.bound  <- log(1 + object@stochInflationUpperBound))
              ans$stoch.log.inf.upper.bound[!is.finite(ans$stoch.log.inf.upper.bound)] <- j.inf
              ans$stoch.log.inf.upper.bound[is.na(ans$stoch.log.inf.upper.bound)] <- j.inf
              ans$stoch.log.inf.upper.bound[ ans$stoch.log.inf.upper.bound > j.inf] <- j.inf

              suppressWarnings(ans$stoch.log.inf.lower.bound  <- log(1 + object@stochInflationLowerBound))
              ans$stoch.log.inf.lower.bound[!is.finite(ans$stoch.log.inf.lower.bound)] <- -j.inf
              ans$stoch.log.inf.lower.bound[is.na(ans$stoch.log.inf.lower.bound)] <- -j.inf
              ans$stoch.log.inf.lower.bound[ans$stoch.log.inf.lower.bound < -j.inf] <- -j.inf

              if(is.na(object@knownStochInflationPersistence))
              {
                  ans$estimate.a.ou <- 1
                  ans$fixed.a.ou <- 0
              } else {
                  ans$estimate.a.ou <- 0
                  ans$fixed.a.ou <- object@knownStochInflationPersistence
              }

              if(is.na(object@knownStochInflationMean))
              {
                  ans$estimate.b.ou <- 1
                  ans$fixed.b.ou <- 0
                  ans$stoch.log.inf.known.mu <- 0
              } else {
                  ans$estimate.b.ou <- 0
                  ans$fixed.b.ou <- 0
                  ans$stoch.log.inf.known.mu <- log(1 + object@knownStochInflationMean)

                  ans$stoch.log.inf.c <- ans$stoch.log.inf.c -  ans$stoch.log.inf.known.mu
              }
          }

          return(ans)

      })

##' A method to collect all the needed initial values common to both the standard model and the break model. Intended for internal use only.
##'
##' There are currenly two types of \code{AnnualAggLossDevModel}s (break and standard).  These models have many features in common, and code to create initial values for these common features is placed in this method.
##' The derived classes \code{StandardAnnualAggLossDevModelInput} and \code{BreakAnnualAggLossDevModelInput} call this method via \code{NextMethod()} and then return a new function.
##'
##' @name getJagsInits,AnnualLossDevModelInput-method
##' @param object An object of type \code{AnnualAggLossDevModelInput} from which to collect the needed initial values for the model.
##' @return A named list of the specific model elements.  See details for more information.
##' @docType methods
setMethod(
          'getJagsInits',
          signature(object='AnnualAggLossDevModelInput'),
          function(object)
      {
          K <- getTriDim(object)[1]
          jd <- getJagsData(object)

          stoch.log.inf.c.obs <- jd$stoch.log.inf.c
          stoch.log.inf.c.init <- rep(NA, length(stoch.log.inf.c.obs))

          center.h <- 0.5
          function()
          {


              stoch.log.inf.c.init[is.na(stoch.log.inf.c.obs)] <- rnorm(length(stoch.log.inf.c.obs[is.na(stoch.log.inf.c.obs)]), mean(stoch.log.inf.c.obs, na.rm=TRUE), 0.01)

              beta.stoch <- c(rbeta(1,100, 100) * range(jd$bounds.for.skewness) + jd$bouds.for.skewness[1], rnorm(1, 0, 0.001))
              beta.stoch.abs <- abs(beta.stoch)
              beta.stoch <- beta.stoch[beta.stoch.abs == min(beta.stoch.abs)][1]
              list(
                   h.stoch=runif(2, center.h / 1.01, center.h * 1.01),
                   h.2.log.stoch=c(rep(NA,2), rnorm(K-2, log(center.h), 0.01)), #very tight to avoid divergent chains
                   tau.h.2.log.innov=rgamma(1,shape=jd$smooth.tau.h.2.log.innov[1], rate=jd$smooth.tau.h.2.log.innov[2]),
                   sigma.kappa=runif(1,0.01,0.1),
                   sigma.eta=runif(1,0.01,0.1),
                   sigma.ou.=runif(1, 0.1 / 1.05, 0.1 * 1.05),
                   a.ou.stoch= runif(1, 0.20, .3),
                   b.ou.stoch=rnorm(1, mean=0, sd=0.1),
                   beta.stoch=beta.stoch,
                   stoch.log.inf.c=stoch.log.inf.c.init
                   )
          }
      })
