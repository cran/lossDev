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


##' @include zzz.R
##' @include AnnualAggLossDevModelOutput.R
##' @include StandardAnnualAggLossDevModelOutput.R
##' @include BreakAnnualAggLossDevModelOutput.R
NULL

##' The class to handle incremental payments of zero for the standard model.
##'
##' \code{StandardAnnualAggLossDevModelOutputWithZeros} is a special class designed to be merged with aggregate annual model objects.
##' It adds one extra node \code{prob.of.non.zero.payment} to the list of slots.
##'
##' @name StandardAnnualAggLossDevModelOutputWithZeros-class
##' @docType class
##' @seealso \code{\linkS4class{LossDevModelOutput}}
##' @seealso \code{\linkS4class{AnnualAggLossDevModelOutputWithZeros}}
setClass(
         'StandardAnnualAggLossDevModelOutputWithZeros',
         representation(
                        prob.of.non.zero.payment='NodeOutput',
                        scale.log='NodeOutput',
                        fifty.fifty.log='NodeOutput'),
         contains=c('StandardAnnualAggLossDevModelOutput'))

##' The class to handle incremental payments of zero for the break model.
##'
##' \code{BreakAnnualAggLossDevModelOutputWithZeros} is a special class designed to be merged with aggregate annual model objects.
##' It adds one extra node \code{prob.of.non.zero.payment} to the list of slots.
##'
##' @name BreakAnnualAggLossDevModelOutputWithZeros-class
##' @docType class
##' @seealso \code{\linkS4class{LossDevModelOutput}}
##' @seealso \code{\linkS4class{AnnualAggLossDevModelOutputWithZeros}}
setClass(
         'BreakAnnualAggLossDevModelOutputWithZeros',
         representation(
                        prob.of.non.zero.payment='NodeOutput',
                        scale.log='NodeOutput',
                        fifty.fifty.log='NodeOutput'),
         contains=c('BreakAnnualAggLossDevModelOutput'))

##' The parent of \code{StandardAnnualAggLossDevModelOutputWithZeros} and \code{BreakAnnualAggLossDevModelOutputWithZeros}.
##'
##' To avoid creating multiple inheritance (directly), this class is created using \code{setClassUnion}.
##' The union consists of classes \code{StandardAnnualAggLossDevModelOutputWithZeros} and \code{StandardAnnualAggLossDevModelOutputWithZeros}.
##'
##' @name AnnualAggLossDevModelOutputWithZeros-class
##' @docType class
##' @seealso \code{\linkS4class{LossDevModelOutput}}
##' @seealso \code{\linkS4class{BreakAnnualAggLossDevModelOutputWithZeros}}
##' @seealso \code{\linkS4class{StandardAnnualAggLossDevModelOutputWithZeros}}
setClassUnion('AnnualAggLossDevModelOutputWithZeros', c('StandardAnnualAggLossDevModelOutputWithZeros', 'BreakAnnualAggLossDevModelOutputWithZeros'))

##' The gompertz function. Intended for internal use only.
##'
##' This function is used as the probably of observing a zero payment.
##' (Or one minus the probably of observing a positive payment.)
##' (Note that negative payments are assumed to be missing values.)
##' @param x The value(s) at which to evaluate the gompertz function.
##' @param scale The scale parameter should always (or at least for how it is used in lossDev) be positive as this indicates an increasing probably of a zero payment.
##' @param fifty.fifty The value at which the gompertz function returns 0.5.
##' @return An object of type \code{AnnualAggLossDevModelOutputWithZeros} and \code{AnnualAggLossDevModelOutput}.
gompertz <- function(x, scale, fifty.fifty)
{
    b <- scale
    c <- -log(-log(0.5)) / b - fifty.fifty

    return(exp(-exp(-b * (x + c))))
}

##' A function to take a triangle estimated without considering zero payments, and account for the possibility of zero payments.
##'
##' As incremental payments are modeled on the log scale, zero payments (and negative payments) are treated as missing values.
##' So, without somehow accounting for zero payments, the estimated payments would be overstated.
##' Zero payments are accounted for by weighting the predicted payment (given that the payment is greater than zero) with the probability that this payment is zero.
##' (Negative payments are not (currently) accounted for.)
##' Currently the trajectory for this probably follows a gompertz curve and is constant across exposure years.
##' This is currently implemented as a function but may be switched to a method.
##'
##' @param object The object containing the triangle estimated without accounting for zero payments.
##' @param burnIn An integer to represent the number of initial \acronym{MCMC} iterations to be discarded. (The adaptive phase (\code{nAddapt}) is not considered part of \code{burnIn}.)
##' @param nAddapt The length of the adaptive phase for the \acronym{MCMC} algorithm. (Default is \code{trunc(burnIn/4)+1}.)
##' @export
##  #import rjags only do this in zzz.R
accountForZeroPayments <- function(object, burnIn=1000, nAddapt=1000)
{
    time.begin <- Sys.time()

    if(is(object, 'StandardAnnualAggLossDevModelOutput'))
    {
        ans <- new('StandardAnnualAggLossDevModelOutputWithZeros')
    } else if(is(object, 'BreakAnnualAggLossDevModelOutput'))
    {
         ans <- new('BreakAnnualAggLossDevModelOutputWithZeros')
    } else {
        stop('Current type of "object" is unsupported')
    }

    for(s in slotNames(object))
        slot(ans, s) <- slot(object, s)

    rm(object)

    u <- getPaymentNoPaymentMatrix(ans@input)

    p.emp <- calculateProbOfPayment(u)

    if(all(p.emp == 1))
        stop('This function can only be called on triangles with "zero" payments.')

    priorsForProbOfPayment <- estimate.priors(p.emp)

    jags.data <- getJagsData(ans@input)
    if(dim(u)[1] != jags.data$K)
        stop('error "dim(u)[1] != jags.data$K"')

    jags.data.new <- list()
    jags.data.new$u <- u
    jags.data.new$scale.prior <- priorsForProbOfPayment['scale']
    jags.data.new$fifty.fifty.prior <- priorsForProbOfPayment['fifty.fifty']
    jags.data.new$K <- jags.data$K
    jags.data.new$H <- jags.data$H
    jags.data.new$L.vec <- jags.data$L.vec



    parameters.to.save. <- c('prob.of.non.zero.payment', 'scale.log', 'fifty.fifty.log')
    ##print(parameters.to.save.)


    ##We will NOT count the addaptive phase as part of the burnin.
    ##ans@burnIn <- as.integer(burnIn)
    ##ans@sampleSize <- as.integer(sampleSize)
    ##ans@thin <- as.integer(thin)
    ##ans@nChains <- as.integer(nChains)

    ##warning('figure out how to set "nAddapt", "burnIn," and "thin"')
    ##nAddapt <- 1000
    eta.mu <- slot(ans@inc.pred, 'value')
    nChains <- dim(eta.mu)['chain']
    ##burnIn <- 1000
    sampleSize <- dim(eta.mu)['iteration']
    thin <- 1

    rngs <- rep(paste('base', c('Wichmann-Hill', 'Marsaglia-Multicarry', 'Super-Duper', 'Mersenne-Twister'), sep='::'), length.out=nChains)

    gen.seed <- function() ceiling(runif(1, 1, 10000))

    rng.seeds <- gen.seed()
    for(i in seq(1, nChains)[-1])
    {
        prop.seed <- gen.seed()
        while(prop.seed %in%  rng.seeds)
            prop.seed <- gen.seed()

        rng.seeds[i] <- prop.seed
    }

    #master.inits.f <- getJagsInits(object)

    inits.f <- function(chain)
    {
        ans <- list()
        ans[['.RNG.name']] <- rngs[chain]
        ans[['.RNG.seed']] <- rng.seeds[chain]

        return(ans)

    }





    message(paste('Preparing Jags Model\nadapting for', nAddapt, 'iterations\n\n'))
    jm <- jags.model(file=file.path(myLibPath(), myPkgName(), 'models', 'probOfPayment.model.txt'),
                     data=jags.data.new,
                     inits=inits.f,
                     n.chains=nChains,
                     n.adapt=nAddapt)


    message(paste('Burning-In Jags Model for', burnIn, 'iterations\n', 'Total Burn-In = ', burnIn))
    update(jm, burnIn)

    message(paste('Sampling Jags Model for', sampleSize, 'iterations Thin =', thin,'\n', 'This will result in ~', sampleSize / thin, 'Samples'))
    output <- jags.samples(jm, parameters.to.save., sampleSize, thin)



    for(i in parameters.to.save.)
        slot(ans,i) <- newNodeOutput(output[[i]])
    ##slot(ans,i) <- new('NodeOutput', value=new('safe.mcmc', value=output[[i]]))

    if(!validObject(ans))
        stop('A valid output could not be created')

          print(paste('Update took', format( Sys.time() - time.begin)))

    return(invisible(ans))


}

##' A function to turn a matrix of incremental payments into zero or ones depending upon whether a payment is positive. Intended for internal use only.
##'
##' The conversion rule is as follows.  If \code{NA}, then \code{NA}. Else if greater than zero, then 1.  Else if equal to zero, then zero. Else \code{NA}.
##'
##' @param object The matrix of incremental payments.
##' @return A matrix of zero or one (or \code{NA}) matching the structure of in input matrix.
getPaymentNoPaymentMatrix <- function(object)
{
    if(!is(object, 'AnnualAggLossDevModelInput'))
        stop('this function currently only supports objects of type "AnnualAggLossDevModelInput"')

    inc <- object@incrementals

    ans <- array(NA, dim(inc))

    f <- function(x)
    {
        if(is.na(x)) {
            NA
        } else if(x > 0) {
            1
        } else if (x == 0) {
            0
        } else {
            NA
        }
    }
    ans <- apply(inc, c(1,2), f)
}

##' A function to calculate an empirical vector of the probability of payment. Intended for internal use only.
##'
##'
##' @param x The matrix of the form returned by \code{\link{getPaymentNoPaymentMatrix}}.
##' @return A vector equal in length to the number of columns in x representing the empirical probably of payment.
calculateProbOfPayment <- function(x)
{
    K <- dim(x)[2]

    ans <- numeric(K)
    for(i in 1:K)
        ans[i] <- sum(x[,i] == 1, na.rm=TRUE) / sum(x[,i] %in% c(0,1), na.rm=TRUE)

    ans[which(is.na(ans))] <- NA

    return(ans)

}

##' A function to estimate priors for the gompertz curve. Intended for internal use only.
##'
##' The function uses \code{nlm} to minimize the squared error.
##'
##' @param p A vector of the form returned by \code{\link{calculateProbOfPayment}}. \code{NA}s are allowed.
##' @return A vector equal in length to the number of columns in x representing the empirical probably of payment.
##' @importFrom stats nlm
estimate.priors <- function(p)
{
    if(all(p == 1))
        stop('This function can only be called on p vectors with some probably of a "zero" payment.')

    K <- length(p)

    scale <- 1
    fifty.fifty <- mean(min(which(!is.na(p) & p != 1)), K)

    f <- function(x)
    {
        p.hat <- 1 - gompertz(1:K, x[1], x[2])
        sse <- sum((p.hat - p)^2, na.rm=TRUE)
        return(sse)
    }

    lse <- nlm(f, c(scale, fifty.fifty))

    if(! lse$code %in% c(1, 2, 3))
        stop('unable to find proper priors')
    ans <- lse$estimate

    names(ans) <- c('scale', 'fifty.fifty')

    return(ans)


}


##' A method to plot and/or return the difference between final actual and predicted cumulative payments.
##'
##' This method accounts for zero payments. By weighting estimated predicted payments by the probably that the payment is greater than zero.
##'
##' The relative difference (x/y - 1) between the final observed cumulative payment and the corresponding predicted cumulative payment is plotted for each exposure year.
##' The horizontal lines of each box represent (starting from the top) the 90th, 75th, 50th, 20th, and 10th percentiles.  Exposure years in which all cumulative payments are \code{NA} are omitted.
##'
##' @name finalCumulativeDiff,AnnualAggLossDevModelOutputWithZeros-method
##' @param object The object of type \code{AnnualAggLossDevModelOuputWithZeros} from which to plot and/or return the difference between final actual and predicted cumulative payments.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @param expYearRange Either a range of years (for example c(1995, 2006)) or one of the keywords \dQuote{all} or \dQuote{fullyObs}.
##' @return Mainly called for the side effect of plotting the difference between final actual and predicted cumulative payments by exposure year.  Also returns a named array for the percentiles in the plot.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{accountForZeroPayments}}
##' @seealso \code{\link{finalCumulativeDiff}}
##' @seealso \code{\link{finalCumulativeDiff,AnnualAggLossDevModelOutput-method}}
setMethod('finalCumulativeDiff',
          signature(object='AnnualAggLossDevModelOutputWithZeros'),
          function(object, plot, expYearRange)
      {

          tmp <- slot(object@inc.pred, 'value') * slot(object@prob.of.non.zero.payment, 'value')
          object@inc.pred <- newNodeOutput(tmp)

          current.class <- class(object)
          i <- match(current.class, is(object))
          f <- selectMethod('finalCumulativeDiff', is(object)[i+1])
          return(f(object, plot))

      })

##' A method to plot and/or return the predicted tail factors for a specific attachment point.
##'
##' This method accounts for zero payments. By weighting estimated predicted payments by the probably that the payment is greater than zero.
##'
##' The tail factor is the ratio of the estimated ultimate loss to cumulative loss at some point in development time.
##' This is a method to allow for the retrieval and illustration of the tail factor by exposure year.
##'
##' Because the model is Bayesian, each tail factor comes as a distribution.  To ease graphical interpretation, only the median for each factor is plotted/returned.
##' See for more details \code{\link{tailFactor}}.
##'
##' For comparison purposes, the function returns three separated tail factors for three scenarios.  Theses three tail factors are returned as a list with the following names and meanings:
##' \describe{
##'   \item{\dQuote{Actual}}{
##'     These are the tail factors estimated when taking the break into consideration.
##'   }
##'   \item{\dQuote{AsIfPostBreak}}{
##'     These are the tail factors estimated when assuming all years where in the post-break regime.
##'   }
##'   \item{\dQuote{AsIfPreBreak}}{
##'     These are the tail factors estimated when assuming all years where in the pre-break regime.
##'   }
##' }
##'
##' @name tailFactor,BreakAnnualAggLossDevModelOutputWithZeros-method
##' @param object The object from which to plot the predicted tail factors and return tail factors for \emph{all} attachment points.
##' @param attachment An integer value specifying the attachment point for the tail.  Must be at least 1. See Details for more info.
##' @param useObservedValues A logical value.  If \code{TRUE}, observed values are substituted for predicted values whenever possible in the calculation.  If \code{FALSE}, only predicted values are used.
##' @param firstIsHalfReport A logical value or \code{NA}.  See Details for more information.
##' @param finalAttachment An integer value must be at least 1 default value is \code{attachment}.  A call to \code{tailFactor} returns (invisibly) a matrix of tail factors through this value.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @param expYearRange Either a range of years (for example c(1995, 2006)) or one of the keywords \dQuote{all} or \dQuote{fullyObs}.
##' @return Mainly called for the side effect of plotting.  Also returns tail factors for \emph{all} attachment points through \code{finalAttachment}.  See Details. Returned invisibly.
##' @docType methods
##' @seealso \code{\link{accountForZeroPayments}}
##' @seealso \code{\link{tailFactor}}
##' @seealso \code{\link[=tailFactor,BreakAnnualAggLossDevModelOutput-method]{tailFactor("BreakAnnualAggLossDevModelOutput")}}
##' @seealso \code{\link[=tailFactor,StandardAnnualAggLossDevModelOutputWithZeros-method]{tailFactor("StandardAnnualAggLossDevModelOutputWithZeros")}}
##' @seealso \code{\link[=tailFactor,StandardAnnualAggLossDevModelOutput-method]{tailFactor("StandardAnnualAggLossDevModelOutput")}}
setMethod('tailFactor',
          signature(object='BreakAnnualAggLossDevModelOutputWithZeros'),
          function(object, attachment, useObservedValues, firstIsHalfReport, finalAttachment, plot, expYearRange)
      {

          prob.of.non.zero.payment.coda <-  slot(object@prob.of.non.zero.payment, 'value')

          tmp.pre <- slot(object@inc.brk.pre, 'value')
          tmp.post <- slot(object@inc.brk.post, 'value')

          tmp.pre <- tmp.pre * prob.of.non.zero.payment.coda
          tmp.post <- tmp.post * prob.of.non.zero.payment.coda


          object@inc.brk.pre <- newNodeOutput(tmp.pre)
          object@inc.brk.post <- newNodeOutput(tmp.post)
          rm(tmp.pre)
          rm(tmp.post)

          tmp <- slot(object@inc.pred, 'value') *  prob.of.non.zero.payment.coda
          object@inc.pred <- newNodeOutput(tmp)
          rm(tmp)



          current.class <- class(object)
          i <- match(current.class, is(object))
          f <- selectMethod('tailFactor', is(object)[i+1])
          return(f(object, attachment, useObservedValues, firstIsHalfReport, finalAttachment, plot, expYearRange))

      })

##' A method to plot and/or return the predicted tail factors for a specific attachment point.
##'
##' This method accounts for zero payments. By weighting estimated predicted payments by the probably that the payment is greater than zero.
##'
##' The tail factor is the ratio of the estimated ultimate loss to cumulative loss at some point in development time.
##' This is a method to allow for the retrieval and illustration of the tail factor by exposure year.
##'
##' Because the model is Bayesian, each tail factor comes as a distribution.  To ease graphical interpretation, only the median for each factor is plotted/returned.
##' See for more details \code{\link{tailFactor}}.
##'
##' For comparison purposes, the function returns three separated tail factors for three scenarios.  Theses three tail factors are returned as a list with the following names and meanings:
##' \describe{
##'   \item{\dQuote{Actual}}{
##'     These are the tail factors estimated when taking the break into consideration.
##'   }
##'   \item{\dQuote{AsIfPostBreak}}{
##'     These are the tail factors estimated when assuming all years where in the post-break regime.
##'   }
##'   \item{\dQuote{AsIfPreBreak}}{
##'     These are the tail factors estimated when assuming all years where in the pre-break regime.
##'   }
##' }
##'
##' @name tailFactor,StandardAnnualAggLossDevModelOutputWithZeros-method
##' @param object The object from which to plot the predicted tail factors and return tail factors for \emph{all} attachment points.
##' @param attachment An integer value specifying the attachment point for the tail.  Must be at least 1. See Details for more info.
##' @param useObservedValues A logical value.  If \code{TRUE}, observed values are substituted for predicted values whenever possible in the calculation.  If \code{FALSE}, only predicted values are used.
##' @param firstIsHalfReport A logical value or \code{NA}.  See Details for more information.
##' @param finalAttachment An integer value must be at least 1 default value is \code{attachment}.  A call to \code{tailFactor} returns (invisibly) a matrix of tail factors through this value.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @param expYearRange Either a range of years (for example c(1995, 2006)) or one of the keywords \dQuote{all} or \dQuote{fullyObs}.
##' @return Mainly called for the side effect of plotting.  Also returns tail factors for \emph{all} attachment points through \code{finalAttachment}.  See Details. Returned invisibly.
##' @docType methods
##' @seealso \code{\link{accountForZeroPayments}}
##' @seealso \code{\link{tailFactor}}
##' @seealso \code{\link[=tailFactor,BreakAnnualAggLossDevModelOutput-method]{tailFactor("BreakAnnualAggLossDevModelOutput")}}
##' @seealso \code{\link[=tailFactor,BreakAnnualAggLossDevModelOutputWithZeros-method]{tailFactor("BreakAnnualAggLossDevModelOutputWithZeros")}}
##' @seealso \code{\link[=tailFactor,StandardAnnualAggLossDevModelOutput-method]{tailFactor("StandardAnnualAggLossDevModelOutput")}}
setMethod('tailFactor',
          signature(object='StandardAnnualAggLossDevModelOutputWithZeros'),
          function(object, attachment, useObservedValues, firstIsHalfReport, finalAttachment, plot, expYearRange)
      {

          tmp <- slot(object@inc.pred, 'value') * slot(object@prob.of.non.zero.payment, 'value')
          object@inc.pred <- newNodeOutput(tmp)

          current.class <- class(object)
          i <- match(current.class, is(object))
          f <- selectMethod('tailFactor', is(object)[i+1])
          return(f(object, attachment, useObservedValues, firstIsHalfReport, finalAttachment, plot, expYearRange))

      })

##' A method to plot predicted vs actual payments for models from the \pkg{lossDev} package.
##'
##' This method accounts for zero payments. By weighting estimated predicted payments by the probably that the payment is greater than zero.
##'
##' Because the model is Bayesian, each estimated payment comes as a distribution.
##' The median of this distribution is used as a point estimate when plotting and/or returning values.
##' Note: One cannot calculate the estimated incremental payments from the estimated cumulative payments (and vice versa) since the median of sums need not be equal to the sum of medians.
##'
##' @name predictedPayments,AnnualAggLossDevModelOutputWithZeros-method
##' @param object The object of type \code{AnnualAggLossDevModelOutputWithZeros} from which to plot predicted vs actual payments and return predicted payments.
##' @param type A singe character value specifying whether to plot/return the predicted incremental or cumulative payments. Valid values are "incremental" or "cumulative."  See details as to why these may not match up.
##' @param logScale A logical value.  If \code{TRUE}, then values are plotted on a log scale.
##' @param mergePredictedWithObserved A logical value.  If \code{TRUE}, then the returned values treat observed incremental payments at "face value"; otherwise predicted values are used in place of observed values.
##' @param plotObservedValues A logical value.  If \code{FALSE}, then only the predicted values are plotted.
##' @param plotPredictedOnlyWhereObserved A logical value.  If \code{TRUE}, then only the predicted incremental payments with valid corresponding observed (log) incremental payment are plotted. Ignored for \code{type="cumulative"}.
##' @param quantiles A vector of quantiles for the predicted payments to return.  Usefull for constructing credible intervals.
##' @param plot A logical value. If \code{TRUE}, then the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a named array (with the same structure as the input triangle) containing the predicted log incremental payments.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{accountForZeroPayments}}
##' @seealso \code{\link{predictedPayments}}
setMethod('predictedPayments',
          signature(object='AnnualAggLossDevModelOutputWithZeros'),
          f <- function(object, type, logScale, mergePredictedWithObserved, plotObservedValues, plotPredictedOnlyWhereObserved, quantiles, plot)
      {


          tmp <- slot(object@inc.pred, 'value') * slot(object@prob.of.non.zero.payment, 'value')
          object@inc.pred <- newNodeOutput(tmp)

          current.class <- class(object)
          i <- match(current.class, is(object))
          f <- selectMethod('predictedPayments', is(object)[i+1])

          return(f(object, type, logScale, mergePredictedWithObserved, plotObservedValues, plotPredictedOnlyWhereObserved, quantiles, plot))

      })

##' A generic function to plot the probability of a payment.
##'
##' Because the model is Bayesian, each estimated payment comes as a distribution.
##' The median of this distribution is used as a point estimate when plotting and/or returning values.
##' Note: Negative payments are treated as missing and are not accounted for.
##'
##' @param object The object from which to plot the probability of a payment.
##' @param plot A logical value. If \code{TRUE}, then the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a matrix containing the (median) probably of payment.  Returned invisibly.
##' @name probablityOfPayment
##' @seealso \code{\link{accountForZeroPayments}}
##' @exportMethod probablityOfPayment
##' @docType genericFunction
setGenericVerif('probablityOfPayment',
                function(object, plot=TRUE)
            {
                standardGeneric('probablityOfPayment')
            })


##' A method to plot the probability of a payment.
##'
##' Because the model is Bayesian, each estimated payment comes as a distribution.
##' The median of this distribution is used as a point estimate when plotting and/or returning values.
##' Note: Negative payments are treated as missing and are not accounted for.
##'
##' @param object The object from which to plot the probability of a payment.
##' @param plot A logical value. If \code{TRUE}, then the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a matrix containing the (median) probably of payment.  Returned invisibly.
##' @name probablityOfPayment,AnnualAggLossDevModelOutputWithZeros-method
##' @seealso \code{\link{accountForZeroPayments}}
##' @docType methods
setMethod('probablityOfPayment',
          signature(object='AnnualAggLossDevModelOutputWithZeros'),
          f <- function(object,  plot)
      {

          if(plot)
          {


              f.plot <- function()
              {
                  u <- getPaymentNoPaymentMatrix(object@input)

                  p.emp <- calculateProbOfPayment(u)


                  matplot(t(object@prob.of.non.zero.payment@median),
                          type='l',
                          lty=1,
                          col=1,
                          xlab='Development Year',
                          ylab='Probability of Payment',
                          lwd=2,
                          cex.axis=1.25,
                          cex.lab=1.25)

                  points(p.emp)
              }
              f.legend <- function()
              {

                  legend('center',
                         c("Fitted","Empirical"),
                         col = c('black','black'),
                         lwd=c(1),
                         lty=c(1, NA),
                         pch=c(NA, 1),
                         horiz=TRUE,
                         bty='n',
                         xpd=NA)
              }

              plot.top.bottom(f.plot, f.legend)
          }

          ans <- object@prob.of.non.zero.payment@median

          dimnames(ans)[[1]] <- object@input@exposureYears[1] - 1 + 1:dim(ans)[1]

          return(invisible(ans))

      })


##' A generic function to plot and/or return the posterior of the parameters for the gompertz curve which describes the probability of payment.
##'
##' The scale parameter describes how steep the curve is.
##' Larger values are steeper.
##' Positive values indicate that the probability of a positive payment should decrease with development time.
##' (The scale is restricted to be positive.)
##'
##' The fifty.fifty parameter gives the point (in development time) when the gompertz curve gives a probability of fifty percent.
##'
##' @name gompertzParameters
##' @param object The object from which to plot and/or return the parameters.
##' @param parameter A character describing which parameter to plot. \dQuote{scale} for the scale parameter. \dQuote{fifty.fifty} for the point at which the gompertz give a probably of fifty percent.
##' @param plotDensity A logical value. If \code{TRUE}, then the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, then the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @return Mainly called for the side effect of plotting.
##' @docType genericFunction
##' @seealso \code{\link{gompertzParameters,AnnualAggLossDevModelOutputWithZeros-method}}
##' @exportMethod gompertzParameters
setGenericVerif('gompertzParameters',
                function(object,  parameter=c('scale', 'fifty.fifty'), plotDensity=TRUE, plotTrace=TRUE)
            {
                standardGeneric('gompertzParameters')
            })


##' A method to plot and/or return the posterior of the parameters for the gompertz curve which describes the probability of payment.
##'
##' The scale parameter describes how steep the curve is.
##' Larger values are steeper.
##' Positive values indicate that the probability of a positive payment should decrease with development time.
##' (The scale is restricted to be positive.)
##'
##' The fifty.fifty parameter gives the point (in development time) when the gompertz curve gives a probability of fifty percent.
##'
##' @name gompertzParameters,AnnualAggLossDevModelOutputWithZeros-method
##' @param object The object from which to plot and/or return the parameters.
##' @param parameter A character describing which parameter to plot. \dQuote{scale} for the scale parameter. \dQuote{fifty.fifty} for the point at which the gompertz give a probably of fifty percent.
##' @param plotDensity A logical value. If \code{TRUE}, then the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, then the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @return Mainly called for the side effect of plotting.
##' @docType genericFunction
##' @seealso \code{\link{gompertzParameters}}
setMethod('gompertzParameters',
          signature(object='AnnualAggLossDevModelOutputWithZeros'),
          function(object,  parameter, plotDensity, plotTrace)
      {


          parameter <- match.arg(parameter)

          if(parameter == 'scale')
              coda <- exp(slot(object@scale.log, 'value')[1,,])
          else
              coda <- exp(slot(object@fifty.fifty.log, 'value')[1,,])


          ans <- plot.density.and.or.trace(coda=coda,
                                           plotDensity = plotDensity ,
                                           plotTrace =   plotTrace,
                                           draw.prior=FALSE,
                                           nice.parameter.name=paste('Gompertz Parameter:', parameter),
                                           zero.line=FALSE)

          return(invisible(ans))
      })



