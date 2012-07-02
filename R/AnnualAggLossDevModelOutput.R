##################################################################################################
##                                                                                              ##
##    lossDev is an R-package.                                                                  ##
##    It is a Bayesian time series model of loss development.                                   ##
##    Features include skewed Student-t distribution with time-varying scale parameters,        ##
##    an expert prior for the calendar year effect,                                             ##
##    and accommodation for structural breaks in the consumption path of services.              ##
##                                                                                              ##
##    Copyright © 2009, 2010, 2011, 2012 National Council On Compensation Insurance Inc.,       ##
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
##    You should have receed a copy of the GNU General Public License                           ##
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.                     ##
##                                                                                              ##
##################################################################################################

##' @include zzz.R
##' @include NodeOutput.R
##' @include LossDevModelOutput.R
NULL

##' The base output class for all aggregate annual models.
##'
##' \code{AnnualAggLossDevModelOutput} is the base output class for all  aggregate annual model objects.
##' Derived classes should contain all output from a \acronym{JAGS} run of the input object in the slot \dQuote{input}.
##' Currenly only the slot \dQuote{input} is allowed to be a non-model node.  All other nodes should be the exact name of some settable node in the model.
##' This is because \code{getModelOutputNodes} currently looks at the slot names to determine what values to set; only slot \dQuote{input} is known to be a slot other than a settable node.
##' This class is derived from \code{LossDevModelOutput}
##' @name AnnualAggLossDevModelOutput-class
##' @docType class
##' @seealso \code{\linkS4class{LossDevModelOutput}}
setClass(
         'AnnualAggLossDevModelOutput',
         representation(inc.pred='NodeOutput',
                        eta='NodeOutput',
                        eta.mu='NodeOutput',
                        sigma.eta='NodeOutput',
                        sigma.kappa='NodeOutput',
                        kappa.log.error='NodeOutput',
                        rho='NodeOutput',
                        rho.eta='NodeOutput',
                        h='NodeOutput',
                        sigma.h.2.log.innov='NodeOutput',
                        beta='NodeOutput',
                        df='NodeOutput',
                        k='NodeOutput',
                        mu.upper.left='NodeOutput',
                        a.ou='NodeOutput',
                        b.ou='NodeOutput',
                        stoch.log.inf.pred='NodeOutput',
                        kappa='NodeOutput',
                        delta.tail='NodeOutput',
                        #omega.obs='NodeOutput',
                        'VIRTUAL'),
         contains='LossDevModelOutput')

##' A generic function to plot and/or return the posterior predicted exposure growth (corresponding to \emph{eta} in the model).
##'
##' @name exposureGrowth
##' @param object The object from which to plot and/or return the posterior predicted exposure growth.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting the exposure growth.  Also returns a named numeric vector for the median of the posterior for the exposure growth on the real (not log) scale.  Returned invisibly.
##' @seealso \code{\link[=exposureGrowth,AnnualAggLossDevModelOutput-method]{exposureGrowth("AnnualAggLossDevModelOutput")}}
##'  \code{\link{exposureGrowthTracePlot}}
##' @exportMethod exposureGrowth
setGenericVerif('exposureGrowth',
                function(object, plot=TRUE)
                standardGeneric('exposureGrowth'))

##' A method to plot and/or return the posterior predicted exposure growth (corresponding to \emph{eta} in the model).
##'
##' @name exposureGrowth,AnnualAggLossDevModelOutput-method
##' @param object The object from which to plot and/or return the posterior predicted exposure growth.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting the exposure growth.  Also returns a named numeric vector for the median of the posterior for the exposure growth on the real (not log) scale.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{exposureGrowth}}
##'  \code{\link{exposureGrowthTracePlot}}
setMethod('exposureGrowth',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object, plot)
      {

          K <- getTriDim(object@input)[1]
          eta <- object@eta@median[-1]

          obs.years <- object@input@exposureYears[-1]
          pred.years <- 1:(length(eta) - (K-1)) + max(obs.years)

          eta.obs <- eta[1:length(obs.years)]
          eta.pred <- eta[1:length(pred.years) + length(obs.years)]

          ans <- c(eta.obs, eta.pred)
          names(ans) <- c(obs.years, pred.years)

          if(plot)
          {

              f.plot <- function()
              {

                  plot(
                       x=range(obs.years, pred.years),
                       y=range(eta),
                       xlab=getExposureYearLabel(object@input),
                       ylab="Rate of Exposure Growth (Net of Calendar Year Effect)",
                       type='n',
                       cex.axis=1.25,
                       cex.lab=1.25)

                  lines(
                        x=obs.years,
                        y=eta.obs,
                        type='o',
                        lty=1,
                        pch=1,
                        lwd=1)

                  lines(
                        x=pred.years,
                        y=eta.pred,
                        type='o',
                        lty=3,
                        pch=20,
                        lwd=2)

                  abline(h=median(exp(slot(object@eta.mu, 'value')) - 1),
                         col='black',
                         lwd=2,
                         lty=2)
              }
              f.legend <- function()
              {

                  legend('center',
                         c('Rate of Exposure Growth','Future Rate of Growth','Stationary Mean'),
                         col=c('black','black','black'),
                         lwd=c(1,2,2),
                         pch=c(1,20,NA),
                         lty=c(1,3,2),
                         horiz=TRUE,
                         bty='n',
                         xpd=NA)
              }

              plot.top.bottom(f.plot, f.legend)
          }


          return(invisible(ans))

      })



##' A generic function to plot and/or return the difference between final actual and predicted cumulative payments.
##'
##' @name finalCumulativeDiff
##' @param object The object from which to plot and/or return the difference.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @param expYearRange Either a range of years (for example c(1995, 2006)) or one of the keywords \dQuote{all} or \dQuote{fullyObs}.
##' @return Mainly called for the side effect of plotting.
##' @seealso \code{\link[=finalCumulativeDiff,AnnualAggLossDevModelOutput-method]{finalCumulativeDiff("AnnualAggLossDevModelOutput")}}
##' @exportMethod finalCumulativeDiff
setGenericVerif('finalCumulativeDiff',
                function(object, plot=TRUE, expYearRange='all')
                standardGeneric('finalCumulativeDiff'))

##' A method to plot and/or return the difference between final actual and predicted cumulative payments.
##'
##' The relative difference (x/y - 1) between the final observed cumulative payment and the corresponding predicted cumulative payment is plotted for each exposure year.
##' The horizontal lines of each box represent (starting from the top) the 90th, 75th, 50th, 20th, and 10th percentiles.  Exposure years in which all cumulative payments are \code{NA} are omitted.
##'
##' If \code{expYearRange} is \dQuote{fullyObs}, then only exposure years with a non missing value in the first column will be plotted.
##'
##' @name finalCumulativeDiff,AnnualAggLossDevModelOutput-method
##' @param object The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the difference between final actual and predicted cumulative payments.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @param expYearRange Either a range of years (for example c(1995, 2006)) or one of the keywords \dQuote{all} or \dQuote{fullyObs}.
##' @return Mainly called for the side effect of plotting the difference between final actual and predicted cumulative payments by exposure year.  Also returns a named array for the percentiles in the plot.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{finalCumulativeDiff}}
setMethod('finalCumulativeDiff',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object, plot, expYearRange)
      {

          K <- getTriDim(object@input)[1]
          inc.pred.coda <- slot(object@inc.pred, 'value')[1:K, 1:K,,]
          cumulatives <- object@input@cumulatives
          exp.years <- object@input@exposureYears

          if(is.character(expYearRange))
          {
              if(length(expYearRange) != 1)
                  stop('"expYearRange" must be of length one if it is a character')
              if(expYearRange != 'all' && expYearRange != 'fullyObs')
                  stop('"expYearRange" must be one of "all" or "fullyObs" if it is supplied as a character')
              if(expYearRange == 'all')
                  expYearRange <- range(exp.years)
              else
                  expYearRange <- range(exp.years[which(!is.na(cumulatives[,1]))])
          } else {

              if(!all(as.integer(expYearRange) == expYearRange))
                  stop('"expYearRange" must be supplied as an integer')
              if(length(expYearRange) != 2)
                    stop('"expYearRange" must have length 2')
              if(max(exp.years) < max(expYearRange) || min(exp.years) > min(expYearRange))
                  stop('"expYearRange" must be a subset of the actual exposure years')
          }


          cumulative.resi.stats <- array(NA, c(5, K), dimnames=list(c('10%', '25%', '50%', '75%', '90%'), exp.years))

          for(i in 1:K)
          {
              tmp <- which(!is.na(cumulatives[i,]))
              if(length(tmp) == 0)
              {
                  cumulative.resi.stats[,i] <- NA
                  next
              }else{
                  last.obs.cumulative.column <- max(tmp)

                  if(length(tmp) == 1 && tmp[1] == 1)
                      diff <-   cumulatives[i,last.obs.cumulative.column] /apply(inc.pred.coda[i,1,,], c(1,2), sum) - 1
                  else
                      diff <-   cumulatives[i,last.obs.cumulative.column] / apply(inc.pred.coda[i,1:last.obs.cumulative.column,,], c(2,3), sum) - 1

                  stats <- quantile(diff, c(.1, .25, .5, .75, .9))

                  cumulative.resi.stats[names(stats),i] <- stats

              }

          }

          if(plot)
          {

              expYearRange.seq <- seq(expYearRange[1], expYearRange[2])
              plot(
                   x=range(exp.years) + c(-1, +1),
                   y=range(as.vector(cumulative.resi.stats[,as.character(expYearRange.seq) ]), na.rm=TRUE),
                   type='n',
                   xlab=getExposureYearLabel(object@input),
                   ylab="Relative Difference Between Actual and Estimated Cumulatives",
                   cex.axis=1.25,
                   cex.lab=1.25)

              abline(h=0,col='gray23',lwd=2,lty='dashed')

              for(i in seq_along(expYearRange.seq))
              {
                  year.i <- expYearRange.seq[i]
                  i. <- match(year.i, object@input@exposureYears)

                  ##draw median to make it thick
                  off.set <- .45
                  lines(x=c(year.i-off.set, year.i+off.set),
                        y=rep(cumulative.resi.stats['50%',i.],2),
                        lwd=2)

                  ##upper 25%
                  off.set <- .45
                  upper.lower <- c('75%','50%')
                  lines(x=c(year.i-off.set, year.i+off.set, year.i+off.set, year.i-off.set, year.i-off.set),
                        y=cumulative.resi.stats[upper.lower[c(1,1,2,2,1)],i.])

                  ##lower 25%
                  off.set <- .45
                  upper.lower <- c('50%','25%')
                  lines(x=c(year.i-off.set, year.i+off.set, year.i+off.set, year.i-off.set, year.i-off.set),
                        y=cumulative.resi.stats[upper.lower[c(1,1,2,2,1)],i.])


                  ##lower%
                  off.set <- .25
                  upper.lower <- c('25%','10%')
                  lines(x=c(year.i-off.set, year.i+off.set, year.i+off.set, year.i-off.set, year.i-off.set),
                        y=cumulative.resi.stats[upper.lower[c(1,1,2,2,1)],i.])


                  ##upper%
                  off.set <- .25
                  upper.lower <- c('90%','75%')
                  lines(x=c(year.i-off.set, year.i+off.set, year.i+off.set, year.i-off.set, year.i-off.set),
                        y=cumulative.resi.stats[upper.lower[c(1,1,2,2,1)],i.])


              }
          }


          return(invisible(cumulative.resi.stats))


      }
          )

##' A generic function to plot and/or return residuals for models in the \pkg{lossDev} package.
##'
##' @name triResi
##' @param object The object from which to plot and/or return the residuals.
##' @param standardize A logical value.  If \code{TRUE}, the plotted and returned residuals are normalized to their respective standard deviation.
##' @param timeAxis A character value describing along which of the three time axes to plot the residuals: \sQuote{dy} for development year time, \sQuote{cy} for calendar year time, \sQuote{ey} for exposure year time.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.
##' @seealso \code{\link[=triResi,AnnualAggLossDevModelOutput-method]{triResi("AnnualAggLossDevModelOutput")}}
##'  \code{\link{QQPlot}}
##' @exportMethod triResi
setGenericVerif('triResi',
                function(object, timeAxis=c('dy', 'cy', 'ey'), standardize=TRUE, plot=TRUE)
                standardGeneric('triResi'))

##' A method to plot and/or return residuals for models in the \pkg{lossDev} package.
##'
##' Because the model is Bayesian, each residual comes as a distribution.  To ease graphical interpretation, only the median for each residual is plotted/returned.
##' The residual is defined as the observed value minus the posterior mean; if standardized, it is also divided by its posterior standard deviation.
##'
##' @name triResi,AnnualAggLossDevModelOutput-method
##' @param object The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the residuals.
##' @param timeAxis A character value describing along which of the three (3) time axis to plot the residuals. \sQuote{dy} for development year time, \sQuote{cy} for calendar year time, \sQuote{ey} for exposure year time.
##' @param standardize A logical value.  If \code{TRUE}, the plotted and returned residuals are normalized to their respective standard deviation.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a named array with the same structure as the input triangle.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{triResi}}
##'  \code{\link{QQPlot}}
setMethod('triResi',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object, timeAxis, standardize, plot)
      {
          timeAxis <- match.arg(timeAxis)

          K <- getTriDim(object@input)[1]

          log.inc <- object@input@incrementals
          log.inc[log.inc <= 0] <- NA
          log.inc <- log(log.inc)

          mu <- slot(object@mu.upper.left, 'value')
          beta <- slot(object@beta, 'value')[1,,]
          v <- slot(object@df, 'value')[1,,]
          h <- slot(object@h, 'value')


          resi <- array(NA, c(K, K), list(object@input@exposureYears, NULL))

          v.factor <-  v / (v - 2)
          h.squared <- h ^ 2
          if(standardize)
          {

              h.to.the.forth <- h.squared ^ 2
              var.second.factor <- 2 * beta ^ 2 *v.factor ^ 2 / (v - 4)
              for(i in 1:K)
                  for(j in 1:K)
                  {
                      if(is.na(log.inc[i,j]))
                          next
                      resi[i,j] <- median((log.inc[i,j] -
                                           (mu[i,j,,] + beta * h.squared[j,,] * v.factor)) /
                                          sqrt(h.squared[j,,] * v.factor + h.to.the.forth[j,,] * var.second.factor))
                  }
          } else {

              for(i in 1:K)
                  for(j in 1:K)
                  {
                      if(is.na(log.inc[i,j]))
                          next
                      resi[i,j] <- median(log.inc[i,j] -
                                          (mu[i,j,,] + beta * h.squared[j,,] * v.factor))
                  }
          }

          if(identical(timeAxis, 'dy'))
          {
              f.plot <- function()
              {
                  plot(x=c(1,K),
                       y=range(resi, na.rm=TRUE),
                       ylab=ifelse(standardize, 'Standardized Residuals', 'Residuals'),
                       xlab='Development Year',
                       type="n",
                       cex.axis=1.25,
                       cex.lab=1.25)
                  abline(a=0,b=0,lwd=2,col='black',lty='dashed')

                  for(i in 1:K)
                      points(x=rep(i,K),
                             y=resi[,i])

                  points(x=1:K, y=apply(resi, 2, median, na.rm=TRUE), lwd=3, type="h", col='red') #bars
                  points(x=1:K, y=apply(resi, 2, median, na.rm=TRUE), pch=20, type="p", col='red') #pinheads
              }

          } else if(identical(timeAxis, 'ey')) {

              exp.years <- object@input@exposureYears
              f.plot <- function()
              {
                  plot(x=range(exp.years),
                       y=range(resi, na.rm=TRUE),
                       ylab=ifelse(standardize, 'Standardized Residuals', 'Residuals'),
                       xlab=getExposureYearLabel(object@input),
                       type="n",
                       cex.axis=1.25,
                       cex.lab=1.25)
                  abline(a=0,b=0,lwd=2,col='black',lty='dashed')

                  for(i in 1:K)
                      points(x=rep(exp.years[i],K),
                             y=resi[i,])

                  points(x=exp.years, y=apply(resi, 1, median, na.rm=TRUE), lwd=3, type="h", col='red') #bars
                  points(x=exp.years, y=apply(resi, 1, median, na.rm=TRUE), pch=20, type="p", col='red') #pinheads
              }

          } else if(identical(timeAxis, 'cy')) {

              cal.years <- object@input@exposureYears
              f.plot <- function()
              {
                  plot(x=range(cal.years),
                       y=range(resi, na.rm=TRUE),
                       ylab=ifelse(standardize, 'Standardized Residuals', 'Residuals'),
                       xlab='Calendar Year',
                       type="n",
                       cex.axis=1.25,
                       cex.lab=1.25)
                  abline(a=0,b=0,lwd=2,col='black',lty='dashed')

                  i <- rep(1:K, K)
                  j <- rep(1:K, rep(K, K))

                  for(k in 1:K)
                  {
                      sub <- resi[i+j-1 == k]
                      l <- length(sub)
                      points(x=rep(cal.years[k],l),
                             y=sub)

                      points(x=cal.years[k],
                             y=median(sub, na.rm=TRUE),
                             lwd=3,
                             type='h',
                             col='red')

                      points(x=cal.years[k],
                             y=median(sub, na.rm=TRUE),
                             lwd=3,
                             type='p',
                             col='red')


                  }

              }
          }

          f.legend <- function()
          {

              legend('center',
                     'Median of Residuals',
                     col=c('red'),
                     pch=c(20),
                     horiz=TRUE,
                     bty='n',
                     xpd=NA)
          }


          if(plot)
              plot.top.bottom(f.plot, f.legend)

          return(invisible(resi))
      })

##' A generic function to plot a Q-Q plot for models in the \pkg{lossDev} package.
##'
##' This function plots sorted observed log incremental payments vs sorted predicted log incremental payments.
##' Credible intervals are also plotted.
##'
##' @name QQPlot
##' @param object The object from which to plot the values.
##' @return Mainly called for the side effect of plotting.
##' @seealso \code{\link[=QQPlot,AnnualAggLossDevModelOutput-method]{QQPlot("AnnualAggLossDevModelOutput")}}
##'  \code{\link{triResi}}
##' @exportMethod QQPlot
setGenericVerif('QQPlot',
                function(object)
                standardGeneric('QQPlot'))

##' A method to plot a Q-Q plot for models in the \pkg{lossDev} package.
##'
##' This function plots sorted observed log incremental payments vs sorted predicted log incremental payments.
##' Credible intervals are also plotted.
##'
##' @name QQPlot,AnnualAggLossDevModelOutput-method
##' @param object The object of type \code{AnnualAggLossDevModelOuput} from which to plot the values.
##' @return NULL. Called for the side effect of plotting.
##' @docType methods
##' @seealso \code{\link{QQPlot}}
##'  \code{\link{triResi}}
setMethod('QQPlot',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object)
      {
          ##plot obs log.inc vs log.inc.pred
          ##by first sorting for every draw the log.inc.pred and then taking quantiles over the sorted values

          f.plot <- function()
          {
              K <- getTriDim(object@input)[1]
              inc.obs <- object@input@incrementals
              inc.pred <- slot(object@inc.pred, 'value')[1:K, 1:K, , ]

              log.inc.obs <- inc.obs
              log.inc.obs[log.inc.obs <= 0] <- NA
              log.inc.obs <- log(log.inc.obs)
              log.inc.obs.not.na <- !is.na(log.inc.obs)

              obs.s <- sort(as.vector(inc.obs[log.inc.obs.not.na]))

              pred.s <- apply(inc.pred, c(3,4), function(x) sort(as.vector(x[log.inc.obs.not.na])))

              pred.s.q <- apply(pred.s, 1, quantile, c(0.05, 0.5, 0.95))



              plot(x=range(obs.s),
                   xlab='Sorted Observed Incrementals (Log Scale)',
                   y=range(pred.s.q),
                   ylab='Sorted Predicted Incrementals (Log Scale)',
                   type='n',
                   cex.axis=1.25,
                   cex.lab=1.25,
                   log='xy')

              lines(x=obs.s,
                    y=pred.s.q[1,])

              points(x=obs.s,
                     y=pred.s.q[2,],
                     cex=1.3)


              lines(x=obs.s,
                    y=pred.s.q[3,])

              abline(a=0,b=1,col='red',lty=2)
          }

          f.legend <- function()
          {

              legend('center',
                     c('Median', '90 Percent\nCredible Intervals', '45 Degree Line'),
                     col=c('black', 'black', 'red'),
                     lty=c(NA, 1, 2),
                     pch=c(1, NA, NA),
                     horiz=TRUE,
                     bty='n',
                     xpd=NA)
          }

          plot.top.bottom(f.plot, f.legend)


      })

##' A generic function to plot and/or return the posterior of the skewness parameter for models in \pkg{lossDev}.
##'
##' The skewness parameter does not directly correspond to the degree of skewness.  However, all else being equal, a larger (in magnitude) skewness parameter indicates a higher degree of skewness,
##' and a skewness parameter of zero equates to zero skew.
##'
##' @references
##' Kim, Y., and J. McCulloch (2007) \dQuote{The Skew-Student Distribution with Application to U.S. Stock Market Returns and the Equity Premium,} Department of Economics, Ohio State University, October 2007
##'
##' @name skewnessParameter
##' @param object The object from which to plot and/or return the skewness parameter.
##' @param plotDensity A logical value. If \code{TRUE}, then the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, then the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @return Mainly called for the side effect of plotting.
##' @seealso \code{\link[=skewnessParameter,AnnualAggLossDevModelOutput-method]{skewnessParameter("AnnualAggLossDevModelOutput")}}
##' @exportMethod skewnessParameter
setGenericVerif('skewnessParameter',
                function(object, plotDensity=TRUE, plotTrace=TRUE)
                standardGeneric('skewnessParameter'))

##' A method to plot and/or return the posterior of the skewness parameter for models in \pkg{lossDev}.
##'
##' The skewness parameter does not directly correspond to the degree of skewness. However, all else being equal, a larger (in magnitude) skewness parameter indicates a higher degree of skewness,
##' and a skewness parameter of zero equates to zero skew.
##'
##' @references
##' Kim, Y., and J. McCulloch (2007) \dQuote{The Skew-Student Distribution with Application to U.S. Stock Market Returns and the Equity Premium,} Department of Economics, Ohio State University, October 2007
##'
##' @name skewnessParameter,AnnualAggLossDevModelOutput-method
##' @param object The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the skewness parameter.
##' @param plotDensity A logical value. If \code{TRUE}, then the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, then the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  But also returns a named array with some select quantiles of the posterior for the skewness parameter.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{skewnessParameter}}
##' @importFrom stats integrate
setMethod('skewnessParameter',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object,  plotDensity, plotTrace)
      {

          if(!object@input@allowForSkew)
          {
              warning('Cannot call "skewnessParameter" unless the model was estimated with a skewed-t. Returning "NULL" invisibly.')
              return(invisible(NULL))
          }

          jd <- getJagsData(object@input)
          precision.for.skewness <- jd$precision.for.skewness
          df.for.skewness <- jd$precision.for.skewness
          mu <- 0
          d.un <- function(x)
          {


              gamma((df.for.skewness+1)/2) / gamma(df.for.skewness / 2) * (precision.for.skewness/ df.for.skewness / pi) ^ 0.5 * (1 + precision.for.skewness / df.for.skewness * (x - mu)^2) ^ (-(df.for.skewness + 1) / 2)

          }

          l <- integrate (d.un, lower = -Inf, upper = jd$bounds.for.skewness[1])$value
          u <- integrate (d.un, lower = -Inf, upper = jd$bounds.for.skewness[2])$value
          d <- function(x)
          {
              if(x < jd$bounds.for.skewness[1] || x > jd$bounds.for.skewness[2])
                  return(0)

              return(d.un(x) / (u - l))

          }
          ans <- plot.density.and.or.trace(coda=slot(object@beta, 'value')[1,,],
                                           plotDensity = plotDensity ,
                                           plotTrace =   plotTrace,
                                           d.prior=d,
                                           nice.parameter.name='Skewness Parameter',
                                           zero.line=TRUE,
                                           lower.bound=jd$bounds.for.skewness[1],
                                           upper.bound=jd$bounds.for.skewness[2])

          return(invisible(ans))

      })

##' A generic function to plot and/or return the posterior of the autoregressive parameter for models in \pkg{lossDev}.
##'
##'
##'
##' @name autoregressiveParameter
##' @param object The object from which to plot and/or return the autoregressive parameter.
##' @param plotDensity A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @return Mainly called for the side effect of plotting.
##' @seealso \code{\link[=autoregressiveParameter,AnnualAggLossDevModelOutput-method]{autoregressiveParameter("AnnualAggLossDevModelOutput")}}
##'  \code{\link{standardDeviationOfCalendarYearEffect}}
##'  \code{\link{calendarYearEffect}}
##'  \code{\link{calendarYearEffectErrors}}
##' @exportMethod autoregressiveParameter
setGenericVerif('autoregressiveParameter',
                function(object, plotDensity=TRUE, plotTrace=TRUE)
            {
                .Deprecated('calendarYearEffectAutoregressiveParameter')
                standardGeneric('autoregressiveParameter')
            })

##' A method to plot and/or return the posterior of the autoregressive parameter for models in \pkg{lossDev}.
##'
##' @name autoregressiveParameter,AnnualAggLossDevModelOutput-method
##' @param object The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the autoregressive parameter.
##' @param plotDensity A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a named array with select quantiles of the posterior for the autoregressive parameter.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{autoregressiveParameter}}
##'  \code{\link{standardDeviationOfCalendarYearEffect}}
##'  \code{\link{calendarYearEffect}}
##'  \code{\link{calendarYearEffectErrors}}
setMethod('autoregressiveParameter',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object,  plotDensity, plotTrace)
      {
          calendarYearEffectAutoregressiveParameter(object, plotDensity, plotTrace)
      })


##' A generic function to plot and/or return the posterior of the autoregressive parameter for the calendar year effect for models in \pkg{lossDev}.
##'
##'
##'
##' @name calendarYearEffectAutoregressiveParameter
##' @param object The object from which to plot and/or return the autoregressive parameter which is associated with the calendar year effect.
##' @param plotDensity A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @return Mainly called for the side effect of plotting.
##' @seealso \code{\link[=calendarYearEffectAutoregressiveParameter,AnnualAggLossDevModelOutput-method]{calendarYearEffectAutoregressiveParameter("AnnualAggLossDevModelOutput")}}
##'  \code{\link{standardDeviationOfCalendarYearEffect}}
##'  \code{\link{calendarYearEffect}}
##'  \code{\link{calendarYearEffectErrors}}
##' @exportMethod calendarYearEffectAutoregressiveParameter
setGenericVerif('calendarYearEffectAutoregressiveParameter',
                function(object, plotDensity=TRUE, plotTrace=TRUE)
            {
                standardGeneric('calendarYearEffectAutoregressiveParameter')
            })

##' A method to plot and/or return the posterior of the autoregressive parameter for the calendar year effect for models in \pkg{lossDev}.
##'
##' @name calendarYearEffectAutoregressiveParameter,AnnualAggLossDevModelOutput-method
##' @param object The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the autoregressive parameter which is associated with the calendar year effect.
##' @param plotDensity A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a named array with select quantiles of the posterior for the autoregressive parameter.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{calendarYearEffectAutoregressiveParameter}}
##'  \code{\link{standardDeviationOfCalendarYearEffect}}
##'  \code{\link{calendarYearEffect}}
##'  \code{\link{calendarYearEffectErrors}}
setMethod('calendarYearEffectAutoregressiveParameter',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object,  plotDensity, plotTrace)
      {

          if(!object@input@ar1InCalendarYearEffect)
          {
              warning('Cannot call "calendarYearEffectAutoregressiveParameter" unless the model was estimated with a autoregressive error term in the calendar year effect. Returning "NULL" invisibly.')
              return(invisible(NULL))
          }

          jd <- getJagsData(object@input)
          ans <- plot.density.and.or.trace(coda=slot(object@rho, 'value')[1,,],
                                           plotDensity = plotDensity ,
                                           plotTrace =   plotTrace,
                                           d.prior=function(x) dbeta(x, jd$rho.prior[1], jd$rho.prior[2]),
                                           nice.parameter.name='Calendar Year AR Parameter',
                                           zero.line=FALSE,
                                           lower.bound=0,
                                           upper.bound=1)

          return(invisible(ans))

      })

##' A generic function to plot and/or return the posterior of the autoregressive parameter for the exposure growth for models in \pkg{lossDev}.
##'
##'
##'
##' @name exposureGrowthAutoregressiveParameter
##' @param object The object from which to plot and/or return the autoregressive parameter which is associated with exposure growth.
##' @param plotDensity A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @return Mainly called for the side effect of plotting.
##' @seealso \code{\link[=exposureGrowthAutoregressiveParameter,AnnualAggLossDevModelOutput-method]{exposureGrowthAutoregressiveParameter("AnnualAggLossDevModelOutput")}}
##' @exportMethod exposureGrowthAutoregressiveParameter
setGenericVerif('exposureGrowthAutoregressiveParameter',
                function(object, plotDensity=TRUE, plotTrace=TRUE)
            {
                standardGeneric('exposureGrowthAutoregressiveParameter')
            })

##' A method to plot and/or return the posterior of the autoregressive parameter for the exposure growth for models in \pkg{lossDev}.
##'
##' @name exposureGrowthAutoregressiveParameter,AnnualAggLossDevModelOutput-method
##' @param object The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the autoregressive parameter which is associated with exposure growth.
##' @param plotDensity A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a named array with select quantiles of the posterior for the autoregressive parameter.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{exposureGrowthAutoregressiveParameter}}
setMethod('exposureGrowthAutoregressiveParameter',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object,  plotDensity, plotTrace)
      {

          if(!object@input@ar1InExposureGrowth)
          {
              warning('Cannot call "exposureGrowthAutoregressiveParameter" unless the model was estimated with a autoregressive error term in the calendar year effect. Returning "NULL" invisibly.')
              return(invisible(NULL))
          }

          jd <- getJagsData(object@input)
          ans <- plot.density.and.or.trace(coda=slot(object@rho.eta, 'value')[1,,],
                                           plotDensity = plotDensity ,
                                           plotTrace =   plotTrace,
                                           d.prior=function(x) dbeta(x, jd$rho.eta.prior[1], jd$rho.eta.prior[2]),
                                           nice.parameter.name='Exposure Growth AR Parameter',
                                           zero.line=FALSE,
                                           lower.bound=0,
                                           upper.bound=1)

          return(invisible(ans))

      })




##' A generic function to plot and/or return the posterior of the mean exposure growth for models in \pkg{lossDev}.
##'
##' (Optionally) exposure growth is modeled as an ar1 process.  This inherently assumes that periods of high exposure growth are (or at least have the possibility of being) followed by continued high periods.
##'
##' @name meanExposureGrowth
##' @param object The object from which to plot and/or return the mean exposure growth.
##' @param plotDensity A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @return Mainly called for the side effect of plotting.
##' @seealso \code{\link[=meanExposureGrowth,AnnualAggLossDevModelOutput-method]{meanExposureGrowth("AnnualAggLossDevModelOutput")}}
##' @exportMethod meanExposureGrowth
setGenericVerif('meanExposureGrowth',
                function(object, plotDensity=TRUE, plotTrace=TRUE)
                standardGeneric('meanExposureGrowth'))

##' A method to plot and/or return the posterior of the mean exposure growth for models in \pkg{lossDev}.
##'
##'
##' @name meanExposureGrowth,AnnualAggLossDevModelOutput-method
##' @param object The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the mean exposure growth.
##' @param plotDensity A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a named array with select quantiles of the posterior for the mean exposure growth.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{meanExposureGrowth}}
setMethod('meanExposureGrowth',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object,  plotDensity, plotTrace)
      {
          ans <- plot.density.and.or.trace(coda=slot(object@eta.mu, 'value')[1,,],
                                           plotDensity = plotDensity ,
                                           plotTrace =   plotTrace,
                                           d.prior=function(x) dnorm(x, 0, sqrt(1/getJagsData(object@input)$precision.for.eta.mu)),
                                           nice.parameter.name='Mean Exposure Growth',
                                           zero.line=TRUE)

          return(invisible(ans))
      })


##' A generic function to plot and/or return the posterior of the degrees of freedom for the Student-\eqn{t} in  models in \pkg{lossDev}.
##'
##' When there is zero skew, the degrees of freedom are the degrees of freedom for the non-skewed \eqn{t}.
##'
##' @references
##' Kim, Y., and J. McCulloch (2007) \dQuote{The Skew-Student Distribution with Application to U.S. Stock Market Returns and the Equity Premium,} Department of Economics, Ohio State University, October 2007.
##'
##' @name degreesOfFreedom
##' @param object The object from which to plot and/or return the degrees of freedom.
##' @param plotDensity A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @return Mainly called for the side effect of plotting.
##' @seealso \code{\link[=degreesOfFreedom,AnnualAggLossDevModelOutput-method]{degreesOfFreedom("AnnualAggLossDevModelOutput")}}
##' @exportMethod degreesOfFreedom
setGenericVerif('degreesOfFreedom',
                function(object, plotDensity=TRUE, plotTrace=TRUE)
                standardGeneric('degreesOfFreedom'))

##' A method to plot and/or return the posterior of the degrees of freedom for the Student-\eqn{t} in  models in \pkg{lossDev}.
##'
##'
##' @name degreesOfFreedom,AnnualAggLossDevModelOutput-method
##' @param object The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the degrees of freedom.
##' @param plotDensity A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a named array with some select quantiles of the posterior for the degrees of freedom.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{degreesOfFreedom}}
setMethod('degreesOfFreedom',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object,  plotDensity, plotTrace)
      {
          jd <- getJagsData(object@input)
          ans <- plot.density.and.or.trace(coda=slot(object@df, 'value')[1,,],
                                           plotDensity = plotDensity ,
                                           plotTrace =   plotTrace,
                                           d.prior=function(x) dchisq(x, df=jd$df.k) / (pchisq(jd$df.bounds[2], jd$df.k) - pchisq(jd$df.bounds[1], jd$df.k)),
                                           nice.parameter.name='Degrees of Freedom',
                                           zero.line=FALSE,
                                           lower.bound=jd$df.bounds[1],
                                           upper.bound=jd$df.bounds[2])

          return(invisible(ans))
      })

##' A generic function to plot and/or return the posterior of the standard deviation of the exposure growth rate for models in \pkg{lossDev}.
##'
##' @name standardDeviationOfExposureGrowth
##' @param object The object from which to plot and/or return the standard deviation of the exposure growth rate.
##' @param plotDensity A logical value. If \code{TRUE}, then the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, then the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @return Mainly called for the side effect of plotting.
##' @seealso \code{\link[=standardDeviationOfExposureGrowth,AnnualAggLossDevModelOutput-method]{standardDeviationOfExposureGrowth("AnnualAggLossDevModelOutput")}}
##' @exportMethod standardDeviationOfExposureGrowth
setGenericVerif('standardDeviationOfExposureGrowth',
                function(object, plotDensity=TRUE, plotTrace=TRUE)
                standardGeneric('standardDeviationOfExposureGrowth'))

##' A method to plot and/or return the posterior of the standard deviation of rhe exposure growth rate for models in \pkg{lossDev}.
##'
##' @name standardDeviationOfExposureGrowth,AnnualAggLossDevModelOutput-method
##' @param object The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the standard deviation of the exposure growth rate.
##' @param plotDensity A logical value. If \code{TRUE}, then the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, then the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a named array with select quantiles of the posterior for the standard deviation of exposure growth.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{standardDeviationOfExposureGrowth}}
##'  \code{\link{exposureGrowth}}
##'  \code{\link{meanExposureGrowth}}
setMethod('standardDeviationOfExposureGrowth',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object,  plotDensity, plotTrace)
      {

          jd <- getJagsData(object@input)
          ans <- plot.density.and.or.trace(coda=slot(object@sigma.eta, 'value')[1,,],
                                           plotDensity = plotDensity ,
                                           plotTrace =   plotTrace,
                                           d.prior=function(x) dunif(x, jd$sigma.eta.bounds[1], jd$sigma.eta.bounds[2]),
                                           nice.parameter.name='Exposure Growth Standard Deviation',
                                           zero.line=FALSE,
                                           lower.bound=jd$sigma.eta.bounds[1])

          return(invisible(ans))
      })

##' A generic function to plot and/or return the posterior of the standard deviation of the calendar year effect for models in \pkg{lossDev}.
##'
##' @name standardDeviationOfCalendarYearEffect
##' @param object The object from which to plot and/or return the standard deviation of the calendar year effect.
##' @param plotDensity A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.
##' @return Mainly called for the side effect of plotting.
##' @exportMethod standardDeviationOfCalendarYearEffect
##' @seealso \code{\link[=standardDeviationOfCalendarYearEffect,AnnualAggLossDevModelOutput-method]{standardDeviationOfCalendarYearEffect("AnnualAggLossDevModelOutput")}}
##'  \code{\link{calendarYearEffect}}
##'  \code{\link{calendarYearEffectErrors}}
##'  \code{\link{autoregressiveParameter}}
setGenericVerif('standardDeviationOfCalendarYearEffect',
                function(object, plotDensity=TRUE, plotTrace=TRUE)
                standardGeneric('standardDeviationOfCalendarYearEffect'))

##' A method to plot and/or return the posterior of the standard deviation of the calendar year effect for models in \pkg{lossDev}.
##'
##' @name standardDeviationOfCalendarYearEffect,AnnualAggLossDevModelOutput-method
##' @param object The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the standard deviation of the calendar year effect.
##' @param plotDensity A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a named array with select quantiles of the posterior for the standard deviation of the calendar year effect.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{standardDeviationOfCalendarYearEffect}}
##'  \code{\link{calendarYearEffect}}
##'  \code{\link{calendarYearEffectErrors}}
##'  \code{\link{autoregressiveParameter}}
setMethod('standardDeviationOfCalendarYearEffect',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object,  plotDensity, plotTrace)
      {

          jd <- getJagsData(object@input)
          ans <- plot.density.and.or.trace(coda=slot(object@sigma.kappa, 'value')[1,,],
                                           plotDensity = plotDensity ,
                                           plotTrace =   plotTrace,
                                           d.prior=function(x) dunif(x, jd$sigma.kappa.bounds[1], jd$sigma.kappa.bounds[2]),
                                           nice.parameter.name='Calendar Effect Standard Deviation',
                                           zero.line=FALSE,
                                           lower.bound=jd$sigma.kappa.bounds[1])

          return(invisible(ans))
      })

##' A generic function to plot and/or return the posterior of the standard deviation for the innovation in the scale parameter for models in \pkg{lossDev}.
##'
##' Changes in the scale parameter (see \code{\link{scaleParameter}}) are assumed to follow a second-order random walk on the log scale.
##' This function plots the posterior standard deviation for this random walk.
##'
##' @name standardDeviationForScaleInnovation
##' @param object The object from which to plot and/or return the standard deviation for the innovation in the log of the scale parameter.
##' @param plotDensity A logical value. If \code{TRUE}, then the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, then the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.
##' @return Mainly called for the side effect of plotting.
##' @seealso \code{\link[=standardDeviationForScaleInnovation,AnnualAggLossDevModelOutput-method]{standardDeviationForScaleInnovation("AnnualAggLossDevModelOutput")}}
##'  \code{\link{standardDeviationVsDevelopmentTime}}
##' @exportMethod standardDeviationForScaleInnovation
setGenericVerif('standardDeviationForScaleInnovation',
                function(object, plotDensity=TRUE, plotTrace=TRUE)
                standardGeneric('standardDeviationForScaleInnovation'))

##' A method to plot and/or return the posterior of the standard deviation for the innovation in the scale parameter for models in \pkg{lossDev}.
##'
##' @name standardDeviationForScaleInnovation,AnnualAggLossDevModelOutput-method
##' @param object The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the standard deviation for the innovation in the log of the scale parameter.
##' @param plotDensity A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a named array with some select quantiles of the posterior for the standard deviation in question.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{standardDeviationForScaleInnovation}}
##'  \code{\link{scaleParameter}}
##'  \code{\link{standardDeviationVsDevelopmentTime}}
setMethod('standardDeviationForScaleInnovation',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object,  plotDensity, plotTrace)
      {

          if(object@input@noChangeInScaleParameterAfterColumn <= 2)
          {
              warning('Cannot call "standardDeviationForScaleInnovation" unless the model was estimated with at least three columns with different scales. Returning "NULL" invisibly.')
              return(invisible(NULL))
          }

          jd <- getJagsData(object@input)
          ans <- plot.density.and.or.trace(coda=slot(object@sigma.h.2.log.innov, 'value')[1,,],
                                           plotDensity = plotDensity,
                                           plotTrace =   plotTrace,
                                           draw.prior = FALSE,
                                           nice.parameter.name='Scale Innovation Standard Deviation',
                                           zero.line=FALSE,
                                           lower.bound=0)

          return(invisible(ans))
      })



##' A generic function to plot and/or return the posterior of the scale parameter for the Student-\eqn{t} measurement equation for models in \pkg{lossDev}.
##'
##' As the degrees of freedom of the \eqn{t} goes to infinity, the scale parameter is the standard deviation of the resulting normal distribution (assuming zero skew).
##'
##' @name scaleParameter
##' @param object The object from which to plot and/or return the scale parameter.
##' @param column The scale parameter is allowed to vary with development time. Setting \code{column} results in the plotting and returning of the scale parameter corresponding to that column. Default value is \code{1}.
##' @param plotDensity A logical value. If \code{TRUE}, then the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, then the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @return Mainly called for the side effect of plotting.
##' @seealso \code{\link[=scaleParameter,AnnualAggLossDevModelOutput-method]{scaleParameter("AnnualAggLossDevModelOutput")}}
##' @exportMethod scaleParameter
setGenericVerif('scaleParameter',
                function(object, column=1, plotDensity=TRUE, plotTrace=TRUE)
            {
                standardGeneric('scaleParameter')
            })

##' A method to plot and/or return the posterior of the scale parameter for the Student-\eqn{t} measurement equation for models in \pkg{lossDev}.
##'
##' @name scaleParameter,AnnualAggLossDevModelOutput-method
##' @param object The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the scale parameter.
##' @param column The scale parameter is allowed to vary with development time. Setting \code{column} results in the plotting and returning of the scale parameter corresponding to that column. Default value is \code{1}.
##' @param plotDensity A logical value. If \code{TRUE}, then the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, then the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a named array with select quantiles of the posterior for the scale parameter.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{scaleParameter}}
setMethod('scaleParameter',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object,  column, plotDensity, plotTrace)
      {

          if(!is.numeric(column))
              stop('"column" must be numeric')
          if(!identical(length(column), as.integer(1)))
              stop('"column" must be of length 1')
          if(column < 1 || column > getTriDim(object@input)[1])
              stop('"column" must be greater than 0 and less than the number of columns in the supplied incremental triangle.')

          jd <- getJagsData(object@input)
          ans <- plot.density.and.or.trace(coda=slot(object@h, 'value')[column,,],
                                           plotDensity = plotDensity ,
                                           plotTrace =   plotTrace,
                                           draw.prior=FALSE,
                                           nice.parameter.name=paste('Scale Parameter:', column),
                                           zero.line=FALSE,
                                           lower.bound=0)

          return(invisible(ans))
      })

##' A generic function to plot and/or return the posterior of the stochastic inflation rho parameter for models in \pkg{lossDev}.
##'
##' If the model incorporates a stochastic rate of inflation, then that rate is assumed to follow (on the log scale) an autoregressive process of order 1.
##' (The autoregressive process of order 1 is the discrete equivalent to an Ornstein-Uhlenbeck process.)
##' This function plots the posterior for the \eqn{rho} parameter, assuming one was estimated.
##'
##' @name stochasticInflationRhoParameter
##' @param object The object from which to plot and/or return the stochastic inflation \eqn{rho} parameter.
##' @param plotDensity A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.
##' @return Mainly called for the side effect of plotting.
##' @seealso \code{\link[=stochasticInflationRhoParameter,AnnualAggLossDevModelOutput-method]{stochasticInflationRhoParameter("AnnualAggLossDevModelOutput")}}
##'  \code{\link{stochasticInflationStationaryMean}}
##'  \code{\link{stochasticInflation}}
##' @exportMethod stochasticInflationRhoParameter
setGenericVerif('stochasticInflationRhoParameter',
                function(object, plotDensity=TRUE, plotTrace=TRUE)
            {
                standardGeneric('stochasticInflationRhoParameter')
            })

##' A method to plot and/or return the posterior of the stochastic inflation \eqn{rho} parameter for models in \pkg{lossDev}.
##'
##' @name stochasticInflationRhoParameter,AnnualAggLossDevModelOutput-method
##' @param object The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the stochastic inflation \eqn{rho} parameter.
##' @param plotDensity A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a named array with select quantiles of the \eqn{rho} parameter.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{stochasticInflationRhoParameter}}
##'  \code{\link{stochasticInflationStationaryMean}}
##'  \code{\link{stochasticInflation}}
setMethod('stochasticInflationRhoParameter',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object, plotDensity, plotTrace)
      {
          if(identical(object@input@stochInflationRate,0) || !is.na(object@input@knownStochInflationPersistence))
          {
              warning('Cannot call "stochasticInflationRhoParameter" unless 1) there is a stochastic rate of inflation and 2) the rho is not known. Returning "NULL" invisibly.')
              return(invisible(NULL))
          }

          jd <- getJagsData(object@input)
          ans <- plot.density.and.or.trace(coda=slot(object@a.ou, 'value')[1,,],
                                           plotDensity = plotDensity ,
                                           plotTrace =   plotTrace,
                                           d.prior = function(x) dbeta(x, jd$a.ou.prior[1], jd$a.ou.prior[2]),
                                           nice.parameter.name='Inflation Autoregressive Parameter',
                                           zero.line=FALSE,
                                           lower.bound=0,
                                           upper.bound=1)

          return(invisible(ans))
      })

##' A generic function to plot and/or return the posterior of the stochastic inflation stationary mean for models in \pkg{lossDev}.
##'
##' If the model incorporates a stochastic rate of inflation, then that rate is assumed to follow (on the log scale) an autoregressive process of order 1.
##' (The autoregressive process of order 1 is the discrete equivalent to an Ornstein-Uhlenbeck process.)
##' This function plots the posterior for the stationary mean (on the log scale), assuming such a mean was estimated.
##'
##' @name stochasticInflationStationaryMean
##' @param object The object from which to plot and/or return the stochastic inflation stationary mean.
##' @param plotDensity A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @return Mainly called for the side effect of plotting.
##' @seealso \code{\link[=stochasticInflationStationaryMean,AnnualAggLossDevModelOutput-method]{stochasticInflationStationaryMean("AnnualAggLossDevModelOutput")}}
##'  \code{\link{stochasticInflationRhoParameter}}
##'  \code{\link{stochasticInflation}}
##' @exportMethod stochasticInflationStationaryMean
setGenericVerif('stochasticInflationStationaryMean',
                function(object, plotDensity=TRUE, plotTrace=TRUE)
            {
                standardGeneric('stochasticInflationStationaryMean')
            })

##' A method to plot and/or return the posterior of the stochastic inflation stationary mean for models in \pkg{lossDev}.
##'
##' @name stochasticInflationStationaryMean,AnnualAggLossDevModelOutput-method
##' @param object The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the stochastic inflation stationary mean.
##' @param plotDensity A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, then only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a named array with select quantiles of the stochastic inflation stationary mean.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{stochasticInflationStationaryMean}}
##'  \code{\link{stochasticInflationRhoParameter}}
##'  \code{\link{stochasticInflation}}
setMethod('stochasticInflationStationaryMean',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object, plotDensity, plotTrace)
      {

          if(identical(object@input@stochInflationRate,0) || !is.na(object@input@knownStochInflationMean))
          {
              warning('Cannot call "stochasticInflationStationaryMean" unless 1) there is a stochastic rate of inflation and 2) the mean is not known. Returning "NULL" invisibly.')
              return(invisible(NULL))
          }


          coda <- slot(object@b.ou, 'value')[1,,] / (1 - slot(object@a.ou, 'value')[1,,])
          ans <- plot.density.and.or.trace(coda=coda,
                                           plotDensity = plotDensity ,
                                           plotTrace =   plotTrace,
                                           draw.prior = FALSE,
                                           nice.parameter.name='(Log) Inflation Stationary Mean',
                                           zero.line=TRUE)

          return(invisible(ans))
      })

##' A generic function to plot and/or return predicted and forecast stochastic inflation rates for models in \pkg{lossDev}.
##'
##' If the model incorporates a stochastic rate of inflation, then that rate is assumed to follow (on the log scale) an autoregressive process of order 1.
##' (The autoregressive process of order 1 is the discrete equivalent to an Ornstein-Uhlenbeck process.)
##' This function plots the median of the posterior predictive distribution for stochastic inflation (not on the log scale) rates by year.
##' Values are returned prior to the application of any limits or weights.
##' Note that for years where observed values are supplied, the model takes those values at face value.
##'
##' @name stochasticInflation
##' @param object The object from which to plot and/or return the stochastic inflation rates.
##' @param extraYears An integer expressing the (maximum) number of years to plot (beyond the final observed year).  Must be at least zero.  Default is 15.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.
##' @seealso \code{\link[=stochasticInflation,AnnualAggLossDevModelOutput-method]{stochasticInflation("AnnualAggLossDevModelOutput")}}
##'  \code{\link{stochasticInflationRhoParameter}}
##'  \code{\link{stochasticInflationStationaryMean}}
##' @exportMethod stochasticInflation
setGenericVerif('stochasticInflation',
                function(object, extraYears=15, plot=TRUE)
            {
                if(!is.numeric(extraYears))
                    stop('"extraYears" must be numeric')
                if(!identical(length(extraYears), as.integer(1)))
                    stop('"extraYears" must be of length 1')
                if(extraYears < 0)
                    stop('"extraYears" must be at least zero.')

                standardGeneric('stochasticInflation')
            })

##' A method to plot and/or return predicted and forecast stochastic inflation rates for models in \pkg{lossDev}.
##'
##' @name stochasticInflation,AnnualAggLossDevModelOutput-method
##' @param object The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return predicted and forecast stochastic inflation rates.
##' @param extraYears An integer expressing the (maximum) number of years to plot (beyond the final observed year).  Must be at least zero.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a named array of the median predicted inflation rate (not on the log scale).  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{stochasticInflation}}
##'  \code{\link{stochasticInflationRhoParameter}}
##'  \code{\link{stochasticInflationStationaryMean}}
setMethod('stochasticInflation',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object, extraYears, plot)
      {

          if(identical(object@input@stochInflationRate,0))
          {
              warning('Cannot call "stochasticInflation" unless there is a stochastic rate of inflation. Returning "NULL" invisibly.')
              return(invisible(NULL))
          }

          #don't have to check extraYears because the generic does that for us

          #the median holds up under log and exp so we don't have to calculate this draw by draw
          simulated.inflation.rate <- exp(object@stoch.log.inf.pred@median) - 1
          simulated.inflation.rate.ci <- apply(exp(slot(object@stoch.log.inf.pred, 'value')) - 1, 1, quantile, c(0.05, 0.95))

          total.simulated.years <- length(simulated.inflation.rate)
          simulated.years <- min(object@input@stochInflationYears) - 1 + 1:total.simulated.years
          names(simulated.inflation.rate) <- simulated.years

          observed.years <- object@input@stochInflationYears
          total.observed.years <- length(observed.years)
          observed.inflation.rate <- object@input@stochInflationRate

          if(total.observed.years >= total.simulated.years)
              extra.years <- 0
          else
              extra.years <- min(extraYears,
                                 total.simulated.years - total.observed.years)



          suppressWarnings( stat.mean <- stochasticInflationStationaryMean(object, plotDensity=FALSE, plotTrace=FALSE))
          if(is.null(stat.mean))
          {
              stat.mean <- object@input@knownStochInflationMean
          } else {
              stat.mean <- stat.mean['50%']
          }

          f.plot <- function()
          {
              plot(x=range(observed.years) + c(0, extra.years),
                   y=range(observed.inflation.rate, simulated.inflation.rate[1:(total.observed.years + extra.years)], stat.mean),
                   xlab="Calendar Year",
                   ylab="Rate of Inflation (Actual and Predicted)",
                   type='n',
                   cex.axis=1.25,
                   cex.lab=1.25)

              abline(h=stat.mean,
                     lwd=2,
                     col='gray',
                     lty=3)

              lines(
                    x=observed.years,
                    y=observed.inflation.rate,
                    lwd=3,
                    col='gray')

              lines(
                    x=observed.years,
                    y=simulated.inflation.rate[1:total.observed.years],
                    lwd=2,
                    col='black')

              if(extra.years > 0 )
              {
                  lines(
                        x=max(observed.years) + 1:extra.years,
                        y=simulated.inflation.rate[total.observed.years + 1:extra.years],
                        lwd=2,
                        col='black',
                        lty=1)

                  for(ind in c('5%', '95%'))
                  {
                      lines(
                            x=max(observed.years) + 1:extra.years,
                            y=simulated.inflation.rate.ci[ind, total.observed.years + 1:extra.years],
                            lwd=2,
                            col='gray',
                            lty=2)
                  }

              }
          }

          f.legend <- function()
          {
              if(extra.years == 0)
                  legend('center',
                         c('Actual','Predicted', 'Stationary\nMean'),
                         col = c('gray','black', 'gray'),
                         lwd=c(3,2,2),
                         lty=c(1,1,3),
                         horiz=TRUE,
                         xpd=NA,
                         bty='n')
              else
                  legend('center',
                         c('Actual', 'Predicted/\nForecast', '90 Percent\nCredible Interval', 'Stationary\nMean'),
                         col = c('gray', 'black', 'gray', 'gray'),
                         lwd=c(3, 2, 2, 2),
                         lty=c(1, 1, 2, 3),
                         horiz=TRUE,
                         xpd=NA,
                         bty='n')
          }


          if(plot)
              plot.top.bottom(f.plot, f.legend)



          return(invisible(simulated.inflation.rate))
      })

##' A generic function to plot and/or return predicted and forecast calendar year effect errors for models in \pkg{lossDev}.
##'
##' The calendar year effect is comprised of two components: 1) a prior expected value which may be unique to every cell (subject to weights and bounds) and 2) a diagonal-specific error term.
##' This function only plots and returns the error term, which includes an autoregressive component if the model is estimated with such a feature.
##'
##' @name calendarYearEffectErrors
##' @param object The object from which to plot and/or return the calendar year effect errors.
##' @param extraYears An integer expressing the (maximum) number of years to plot (beyond the final observed calendar year).  Must be greater than or equal to zero.  Default is 15.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.
##' @seealso \code{\link[=calendarYearEffectErrors,AnnualAggLossDevModelOutput-method]{calendarYearEffectErrors("AnnualAggLossDevModelOutput")}}
##'  \code{\link{calendarYearEffect}}
##'  \code{\link{autoregressiveParameter}}
##'  \code{\link{standardDeviationOfCalendarYearEffect}}
##'  \code{\link{calendarYearEffectErrorTracePlot}}
##' @exportMethod calendarYearEffectErrors
setGenericVerif('calendarYearEffectErrors',
                function(object, extraYears=15, plot=TRUE)
            {
                if(!is.numeric(extraYears))
                    stop('"extraYears" must be numeric')
                if(!identical(length(extraYears), as.integer(1)))
                    stop('"extraYears" must be of length 1')
                if(extraYears < 0)
                    stop('"extraYears" must be at least zero.')

                standardGeneric('calendarYearEffectErrors')
            })

##' A method to plot and/or return predicted and forecast calendar year effect errors for models in \pkg{lossDev}.
##'
##' @name calendarYearEffectErrors,AnnualAggLossDevModelOutput-method
##' @param object The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the calendar year effect errors.
##' @param extraYears An integer expressing the (maximum) number of years to plot (beyond the final observed calendar year).  Must greater than or equal to zero.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a named array with the median predicted errors (not on the log scale).  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{calendarYearEffectErrors}}
##'  \code{\link{calendarYearEffect}}
##'  \code{\link{autoregressiveParameter}}
##'  \code{\link{standardDeviationOfCalendarYearEffect}}
##'  \code{\link{calendarYearEffectErrorTracePlot}}
setMethod('calendarYearEffectErrors',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object, extraYears, plot)
      {

          ##don't have to check extraYears because the generic does that for us

          ##the median holds up under log and exp so we don't have to calculate this draw by draw
          kappa.error <- exp(object@kappa.log.error@median[-(1:2)]) - 1 #first value is for diagonal not in the triangle, second is for the first diagonal in the triangle which has no identifiable effect
          total.years <- length(kappa.error)
          years <- min(object@input@exposureYears) + 1 - 1 + 1:total.years
          names(kappa.error) <- years

          observed.years <- object@input@exposureYears[-1]
          total.observed.years <- length(observed.years)

          if(total.observed.years >= total.years)
              extra.years <- 0
          else
              extra.years <- min(extraYears,
                                 total.years - total.observed.years)



          f.plot <- function()
          {
              plot(x=range(observed.years) + c(0, extra.years),
                   y=range(kappa.error[1:(total.observed.years + extra.years)]),
                   xlab="Calendar Year",
                   ylab="Calendar Effect Error",
                   type='n',
                   cex.axis=1.25,
                   cex.lab=1.25)

              abline(h=0,
                     lty=2)

              lines(
                    x=observed.years,
                    y=kappa.error[1:total.observed.years],
                    lwd=2,
                    col='black')

              if(extra.years > 0 )
                  lines(
                        x=max(observed.years) + 1:extra.years,
                        y=kappa.error[total.observed.years + 1:extra.years],
                        lwd=2,
                        col='gray')


          }

          f.legend <- function()
          {
              if(extra.years == 0)
                  legend('center',
                         c('Estimated'),
                         col = c('black'),
                         lwd=c(2),
                         horiz=TRUE,
                         xpd=NA,
                         bty='n')
              else
                  legend('center',
                         c('Estimated', 'Predicted'),
                         col = c('black', 'gray'),
                         lwd=c(2, 2),
                         horiz=TRUE,
                         xpd=NA,
                         bty='n')
          }


          if(plot)
              plot.top.bottom(f.plot, f.legend)



          return(invisible(kappa.error))
      })


##' A generic function to plot and/or return the predicted and forecast calendar year effects for models in \pkg{lossDev}.
##'
##' The calendar year effect is comprised of two components: 1) a prior expected value that may be unique to every cell (subject to weights and bounds) and 2) a diagonal-specific error term.
##' This function plots and returns the factor resulting from the combined effect of these two, which includes an autoregressive component if the model is estimated with such a feature.
##'
##' The first cell is \code{NA}. Values in the first column represent the rate of inflation/escalation to the corresponding cell from the cell in the same column but previous row.
##' Values in the 2nd column and beyond represent the rate of inflation/escalation to the corresponding cell from the cell in the same row but previous column.
##'
##' @name calendarYearEffect
##' @param object The object from which to plot and/or return the calendar year effect.
##' @param restrictedSize A logical value.  If \code{TRUE}, the plotted calendar year effect is restricted to the square of dimension equal to the observed triangle with which the model was estimated.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.
##' @seealso \code{\link[=calendarYearEffect,AnnualAggLossDevModelOutput-method]{calendarYearEffect("AnnualAggLossDevModelOutput")}}
##'  \code{\link{calendarYearEffectErrors}}
##'  \code{\link{autoregressiveParameter}}
##'  \code{\link{standardDeviationOfCalendarYearEffect}}
##'  \code{\link{calendarYearEffectErrorTracePlot}}
##' @exportMethod calendarYearEffect
setGenericVerif('calendarYearEffect',
                function(object, restrictedSize=FALSE, plot=TRUE)
            {
                standardGeneric('calendarYearEffect')
            })

##' A method to plot and/or return predicted and forecast calendar year effects for models in \pkg{lossDev}.
##'
##' The first cell is \code{NA}. Values in the first column represent the rate of inflation/escalation to the corresponding cell from the cell in the same column but previous row.
##' Values in the 2nd column and beyond represent the rate of inflation/escalation to the corresponding cell from the cell in the same row but previous column.
##'
##' @name calendarYearEffect,AnnualAggLossDevModelOutput-method
##' @param object The object of type \code{AnnualAggLossDevModelOuput} from which to plot and/or return the calendar year effect.
##' @param restrictedSize A logical value.  If \code{TRUE}, the plotted calendar year effect is restricted to the square of dimension equal to the observed triangle with which the model was estimated.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a named array with the median predicted values (not on the log scale).  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{calendarYearEffect}}
##'  \code{\link{calendarYearEffectErrors}}
##'  \code{\link{autoregressiveParameter}}
##'  \code{\link{standardDeviationOfCalendarYearEffect}}
##'  \code{\link{calendarYearEffectErrorTracePlot}}
setMethod('calendarYearEffect',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object, restrictedSize, plot)
      {
          if(plot)
          {
              K <- getTriDim(object@input)[1]
              kappa <- object@kappa@median
              kappa[1,1] <- NA

              if(restrictedSize)
                  kappa <- kappa[1:K,1:K]

              all.exp.years <- min(object@input@exposureYears) - 1 + 1:dim(kappa)[1]

              data <- expand.grid(row=all.exp.years, column=1:ncol(kappa))
              data$z <- as.vector(as.numeric(kappa))

              print(
                    levelplot(z ~  column * row,
                              data,
                              aspect='iso',
                              ylim=(range(all.exp.years) + c(-0.5, 0.5))[c(2,1)],
                              col.regions=grey(seq(from=1, to=0, length.out=100)),
                              xlab='Development Year',
                              ylab=getExposureYearLabel(object@input)))
          }

          ans <- object@kappa@median
          ans[1,1] <- NA

          dimnames(ans)[[1]] <- min(object@input@exposureYears) - 1 + 1:dim(ans)[1]

          return(invisible(ans))
      })


##' A generic function to plot predicted vs actual payments for models from the \pkg{lossDev} package.
##'
##' Because the model is Bayesian, each estimated payment comes as a distribution.
##' The median of this distribution is used as a point estimate when plotting and/or returning values.
##' Note: One cannot calculate the estimated incremental payments from the estimated cumulative payments (and vice versa) since the median of sums need not be equal to the sum of medians.
##'
##' If \code{mergePredictedWithObserved=TRUE} and \code{type="incremental"}, then any observed incremental payment will be used in place of its corresponding incremental payment.
##' If \code{mergePredictedWithObserved=TRUE} and \code{type="cumulative"}, then only predicted incremental payments (by row) to the right of the last observed cumulative value will enter the calculation.
##'
##' @name predictedPayments
##' @param object The object from which to plot predicted vs actual payments and from which to return predicted payments.
##' @param type A single character value specifying whether to plot/return the predicted incremental or cumulative payments. Valid values are \dQuote{incremental} or \dQuote{cumulative.}  See details as to why these may not match up.
##' @param logScale A logical value.  If \code{TRUE}, then values are plotted on a log scale.
##' @param mergePredictedWithObserved A logical value.  See details.
##' @param plotObservedValues A logical value.  If \code{FALSE}, then only the predicted values are plotted.
##' @param plotPredictedOnlyWhereObserved A logical value.  If \code{TRUE}, then only the predicted incremental payments with valid corresponding observed (log) incremental payment are plotted. Ignored for \code{type="cumulative"}.
##' @param quantiles A vector of quantiles for the predicted payments to return.  Useful for constructing credible intervals.
##' @param plot A logical value. If \code{TRUE}, then the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.
##' @seealso \code{\link[=predictedPayments,AnnualAggLossDevModelOutput-method]{predictedPayments("AnnualAggLossDevModelOutput")}}
##' @exportMethod predictedPayments
setGenericVerif('predictedPayments',
                function(object, type=c('incremental', 'cumulative'), logScale=TRUE, mergePredictedWithObserved=FALSE, plotObservedValues=logScale, plotPredictedOnlyWhereObserved=FALSE, quantiles=c(0.05, .5, 0.95), plot=TRUE)
                standardGeneric('predictedPayments'))

##' A method to plot predicted vs actual payments for models from the \pkg{lossDev} package.
##'
##' Because the model is Bayesian, each estimated payment comes as a distribution.
##' The median of this distribution is used as a point estimate when plotting and/or returning values.
##' Note: One cannot calculate the estimated incremental payments from the estimated cumulative payments (and vice versa) since the median of sums need not be equal to the sum of medians.
##'
##' If \code{mergePredictedWithObserved=TRUE} and \code{type="incremental"}, then any observed incremental payment will be used in place of its corresponding incremental payment.
##' If \code{mergePredictedWithObserved=TRUE} and \code{type="cumulative"}, then only predicted incremental payments (by row) to the right of the last observed cumulative value will enter the calculation.
##'
##' @name predictedPayments,AnnualAggLossDevModelOutput-method
##' @param object The object of type \code{AnnualAggLossDevModelOutput} from which to plot predicted vs actual payments and return predicted payments.
##' @param type A singe character value specifying whether to plot/return the predicted incremental or cumulative payments. Valid values are "incremental" or "cumulative."  See details as to why these may not match up.
##' @param logScale A logical value.  If \code{TRUE}, then values are plotted on a log scale.
##' @param mergePredictedWithObserved A logical value.  If \code{TRUE}, then the returned values treat observed incremental payments at "face value"; otherwise predicted values are used in place of observed values.
##' @param plotObservedValues A logical value.  If \code{FALSE}, then only the predicted values are plotted.
##' @param plotPredictedOnlyWhereObserved A logical value.  See details.
##' @param quantiles A vector of quantiles for the predicted payments to return.  Usefull for constructing credible intervals.
##' @param plot A logical value. If \code{TRUE}, then the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a named array (with the same structure as the input triangle) containing the predicted log incremental payments.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{predictedPayments}}
setMethod('predictedPayments',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object, type, logScale, mergePredictedWithObserved, plotObservedValues, plotPredictedOnlyWhereObserved, quantiles, plot)
      {
          type <- match.arg(type)

          K <- getTriDim(object@input)[1]

          inc.obs <- object@input@incrementals
          cumul.obs <- object@input@cumulatives

          inc.pred.coda <- slot(object@inc.pred, 'value')
          inc.pred.median <- object@inc.pred@median

          if(type == 'cumulative')
          {
              cumul.pred.coda <- array(NA, dim(inc.pred.coda))
              cumul.pred.coda[,1,,] <-  inc.pred.coda[,1,,]

              for(i in 2:dim(inc.pred.coda)[2])
                  cumul.pred.coda[,i,,] <- cumul.pred.coda[,i-1,,] + inc.pred.coda[,i,,]

              cumul.pred.median <- apply(cumul.pred.coda, c(1,2), median)

          }

          ##trim down the predicted values to the size of the original triangle
          inc.pred.median.trim <- inc.pred.median[1:K, 1:K]
          ##get rid of predicted values which do not correspond to observed values
          inc.pred.median.trim.only.where.obs <- inc.pred.median.trim
          inc.pred.median.trim.only.where.obs[inc.obs <= 0 | is.na(inc.obs)] <- NA




          if(type == 'incremental')
          {
              ans <- inc.pred.coda
              dimnames(ans)[[1]] <-  min(object@input@exposureYears) - 1 + 1:dim(ans)[1]

              if(mergePredictedWithObserved)
                  ans[1:K, 1:K, , ][!is.na(inc.obs)] <- inc.obs[!is.na(inc.obs)]

              ans <- apply(ans, c(1, 2), quantile, quantiles)


          } else {

              if(!mergePredictedWithObserved)
              {
                  ans <- cumul.pred.coda
                  dimnames(ans)[[1]] <-  min(object@input@exposureYears) - 1 + 1:dim(ans)[1]
                  ans <- apply(ans, c(1, 2), quantile, quantiles)

              } else {

                  tmp <- array(NA,
                               dim(inc.pred.coda),
                               dimnames=list(min(object@input@exposureYears) - 1 + 1:dim(inc.pred.coda)[1], NULL))

                  tmp[1:K, 1:K, , ] <- cumul.obs



                  for(i in 1:K)
                  {
                      j.lower <- which(!is.na(cumul.obs[i,]))
                      if(length(j.lower) == 0)
                      {
                          j.lower <- 1
                      } else {
                          j.lower <- max(j.lower)
                      }

                      for(j in (j.lower+1):(dim(inc.pred.coda)[2]))
                          tmp[i,j,,] <- tmp[i,j-1,,] +  inc.pred.coda[i,j,,]
                  }

                  tmp[(K+1):(dim(cumul.pred.coda)[1]), , , ] <- cumul.pred.coda[(K+1):(dim(cumul.pred.coda)[1]), , , ]
                  ans <- apply(tmp, c(1,2), quantile, probs = quantiles, na.rm = TRUE)

              }

          }





          plot.f <- function()
          {
              if(type == 'incremental')
              {
                  if(plotPredictedOnlyWhereObserved)
                  {
                      inc.pred.for.plot <- inc.pred.median.trim.only.where.obs
                  } else {
                      inc.pred.for.plot <- inc.pred.median
                  }

                  if(logScale)
                  {
                      inc.pred.for.plot[inc.pred.for.plot <= 0] <- NA
                      inc.obs.for.plot <- inc.obs
                      inc.obs.for.plot[inc.obs.for.plot <= 0] <- NA
                  } else {
                      inc.obs.for.plot <- inc.obs
                  }

                  x.range <- c(1, dim(inc.pred.for.plot)[2])

                  if(plotObservedValues)
                      y.range <- range(inc.pred.for.plot, inc.obs.for.plot, na.rm=TRUE)
                  else
                      y.range <- range(inc.pred.for.plot, na.rm=TRUE)

                  plot(x=x.range,
                       y=y.range,
                       ylab='Incremental Payments',
                       xlab='Development Year',
                       type='n',
                       log=ifelse(logScale, 'y', ''),
                       cex.axis=1.25,
                       cex.lab=1.25)

                  for(i in 1:dim(inc.pred.for.plot)[1])
                  {
                      tmp <- inc.pred.for.plot[i,]
                      lines(x=1:length(tmp),
                            y=tmp,
                            col=get.color(i),
                            type=ifelse(plotPredictedOnlyWhereObserved, 'o', 'l'),
                            pch=20)
                  }

                  if(plotObservedValues)
                  {
                      for(i in 1:dim(inc.obs.for.plot)[1])
                      {
                          tmp <- inc.obs.for.plot[i,]
                          lines(x=1:length(tmp),
                                y=tmp,
                                type='p',
                                col=get.color(i))
                      }
                  }
              } else {

                  cumul.pred.for.plot <- cumul.pred.median

                  if(logScale)
                  {
                      cumul.pred.for.plot[cumul.pred.for.plot <= 0] <- NA
                      cumul.obs.for.plot <- cumul.obs
                      cumul.obs.for.plot[cumul.obs.for.plot <= 0] <- NA
                  } else {
                      cumul.obs.for.plot <- cumul.obs
                  }

                  x.range <- c(1, dim(cumul.pred.for.plot)[2])

                  if(plotObservedValues)
                      y.range <- range(cumul.pred.for.plot, cumul.obs.for.plot, na.rm=TRUE)
                  else
                      y.range <- range(cumul.pred.for.plot, na.rm=TRUE)

                  plot(x=x.range,
                       y=y.range,
                       ylab='Cumulative Payments',
                       xlab='Development Year',
                       type='n',
                       log=ifelse(logScale, 'y', ''),
                       cex.axis=1.25,
                       cex.lab=1.25)

                  for(i in 1:dim(cumul.pred.for.plot)[1])
                  {
                      tmp <- cumul.pred.for.plot[i,]
                      lines(x=1:length(tmp),
                            y=tmp,
                            col=get.color(i))
                  }

                  if(plotObservedValues)
                  {
                      for(i in 1:dim(cumul.obs.for.plot)[1])
                      {
                          tmp <- cumul.obs.for.plot[i,]
                          lines(x=1:length(tmp),
                                y=tmp,
                                type='p',
                                col=get.color(i))
                      }
                  }
              }
          }

          legend.f <- function()
          {
              if(plotObservedValues)
              {
                  legend('center',
                         legend=c('Predicted','Observed'),
                         col='black',
                         lwd=2,
                         lty=c(1, NA),
                         pch=c(ifelse(type=='incremental' && plotPredictedOnlyWhereObserved, 20, NA), 1),
                         horiz=TRUE,
                         bty='n',
                         xpd=NA)
              }

          }

          if(plot)
              plot.top.bottom(plot.f, legend.f)


          return(invisible(ans))
      })

##' A generic function to plot and/or return the posterior estimated standard deviation by development year.
##'
##' Aggregate loss development models in \pkg{lossDev} allow for changes (by development year) in the measurement error around the log incremental payments.
##' This is a generic function that allows for the retrieval and illustration of this standard deviation.
##'
##' @name standardDeviationVsDevelopmentTime
##' @param object The object from which to plot and/or return the estimated standard deviation by development year.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.
##' @seealso \code{\link[=standardDeviationVsDevelopmentTime,AnnualAggLossDevModelOutput-method]{standardDeviationVsDevelopmentTime("AnnualAggLossDevModelOutput")}}
##' @exportMethod standardDeviationVsDevelopmentTime
setGenericVerif('standardDeviationVsDevelopmentTime',
                function(object, plot=TRUE)
                standardGeneric('standardDeviationVsDevelopmentTime'))

##' A method to plot and/or return the posterior estimated standard deviation by development year.
##'
##' Aggregate loss development models in \pkg{lossDev} allow for changes (by development year) in the measurement error around the log incremental payments.
##' This is a method that allows for the retrieval and illustration of this standard deviation.
##'
##' @name standardDeviationVsDevelopmentTime,AnnualAggLossDevModelOutput-method
##' @param object The object of type \code{AnnualAggLossDevModelOutput} from which to plot and/or return the estimated standard deviation by development year.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a numeric vector of the plotted statistics.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{standardDeviationVsDevelopmentTime}}
setMethod('standardDeviationVsDevelopmentTime',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object, plot)
      {
          K <- getTriDim(object@input)[1]

          beta <- slot(object@beta, 'value')[1,,]
          v <- slot(object@df, 'value')[1,,]
          h <- slot(object@h, 'value')

          st.d <- array(NA, c(3, K), list(c('95%', '50%', '5%'), NULL))

          for(i in 1:K)
          {
              tmp <- quantile(
                              sqrt(h[i,,] ^ 2 * v / (v - 2) + 2 * beta ^ 2 * h[i,,] ^ 4 * v ^ 2 / (v - 2) ^ 2 / (v - 4)),
                              c(.95, .5, .05))

              st.d[names(tmp),i] <- tmp
          }

          f.plot <- function()
          {
               matplot(y=t(st.d),
                      x=1:K,
                      col=c('grey', 'black', 'grey'),
                      lwd=c(3,2,3),
                      lty=c(2,1,2),
                      xlab='Development Year',
                      ylab="Standard Deviation in Measurement Equation",
                      cex.axis=1.25,
                      cex.lab=1.25,
                       type='l')
           }

          f.legend <- function()
          {
              legend('center',
                     legend=c('Median','90 Percent Credible Interval'),
                     col=c('black', 'grey'),
                     lwd=c(2,3),
                     lty=c(1,2),
                     horiz=TRUE,
                     bty='n',
                     xpd=NA)
          }

          if(plot)
              plot.top.bottom(f.plot, f.legend)


          return(invisible(st.d))

      })


##' A generic function to generate the trace plots for select exposure growth rates.
##'
##'
##' @name exposureGrowthTracePlot
##' @param object The object from which to generate the trace plots.
##' @param elements A numeric vector indicating the elements for which to plot the trace.  Valid values are 2 through the total number of exposure years.  If NULL, values are selected automatically.
##' @return NULL invisibly.  Only called for the side effect of plotting.
##' @seealso \code{\link[=exposureGrowthTracePlot,AnnualAggLossDevModelOutput-method]{exposureGrowthTracePlot("AnnualAggLossDevModelOutput")}}
##'  \code{\link{exposureGrowth}}
##' @exportMethod exposureGrowthTracePlot
setGenericVerif('exposureGrowthTracePlot',
                function(object, elements=NULL)
                standardGeneric('exposureGrowthTracePlot'))

##' A method to generate the trace plots for select exposure growth rates.
##'
##'
##' @name exposureGrowthTracePlot,AnnualAggLossDevModelOutput-method
##' @param object The object of type \code{AnnualAggLossDevModelOutput} from which to generate the trace plots.
##' @param elements A numeric vector indicating the elements for which to plot the trace. If NULL, values are selected automatically.
##' @return NULL invisibly.  Only called for the side effect of plotting.
##' @docType methods
##' @seealso \code{\link{exposureGrowthTracePlot}}
##'  \code{\link{exposureGrowth}}
setMethod('exposureGrowthTracePlot',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object, elements)
      {
          lb <- 2
          ub <- object@input@totalExpYears

          if(is.null(elements))
          {
              elements <- c(lb, floor((lb + ub) / 2), ub)
          } else {

              if(!is.numeric(elements) || !is.vector(elements))
                  stop('"elements" must either be "NULL" or a numeric vector')

              if(any(elements != as.integer(elements)))
                  stop('"elements" must be coercible to integers')

              if(any(elements > ub) || any(elements < lb))
                  stop(paste('"elements" must be at most the total number of exposure years (', ub,') and at least 2', sep=''))
          }

          plot.trace.plots(slot(object@eta, 'value')[elements,,], paste('Exposure Growth :', elements, sep=''))


      })

##' A generic function to generate the trace plots for select calendar year effect errors.
##'
##' The calendar year effect is comprised of two components: 1) a prior expected value that may be unique to every cell and 2) a diagonal-specific error term.
##' This function generates trace plots for the diagonal specific error terms only.
##'
##' @name calendarYearEffectErrorTracePlot
##' @param object The object from which to generate the trace plots.
##' @param elements A numeric vector indicating the elements for which to plot the trace.  Valid values are 2 through the total number of exposure years(observed and forecast).  If NULL, values are selected automatically.
##' @return NULL invisibly.  Only called for the side effect of plotting.
##' @seealso \code{\link[=calendarYearEffectErrorTracePlot,AnnualAggLossDevModelOutput-method]{calendarYearEffectErrorTracePlot("AnnualAggLossDevModelOutput")}}
##'  \code{\link{calendarYearEffectErrors}}
##' @exportMethod calendarYearEffectErrorTracePlot
setGenericVerif('calendarYearEffectErrorTracePlot',
                function(object, elements=NULL)
                standardGeneric('calendarYearEffectErrorTracePlot'))

##' A method to generate the trace plots for select calendar year effect errors.
##'
##' @name calendarYearEffectErrorTracePlot,AnnualAggLossDevModelOutput-method
##' @param object The object of type \code{AnnualAggLossDevModelOutput} from which to generate the trace plots.
##' @param elements A numeric vector indicating the elements for which to plot the trace.  Valid values are 2 through the total number of exposure years. If NULL, values are selected automatically.
##' @return NULL invisibly.  Only called for the side effect of plotting.
##' @docType methods
##' @seealso \code{\link{calendarYearEffectErrorTracePlot}}
##'  \code{\link{calendarYearEffectErrors}}
setMethod('calendarYearEffectErrorTracePlot',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object, elements)
      {
          lb <- 2
          ub <- object@input@totalExpYears

          if(is.null(elements))
          {
              elements <- c(lb, floor((lb + ub) / 2), ub)
          } else {

              if(!is.numeric(elements) || !is.vector(elements))
                  stop('"elements" must either be "NULL" or a numeric vector')

              if(any(elements != as.integer(elements)))
                  stop('"elements" must be coercible to integers')

              if(any(elements > ub) || any(elements < lb))
                  stop(paste('"elements" must be at most the total number of exposure years (', ub,') and at least 3', sep=''))
          }


          elements <- elements + 1 #in the model file, the first identifiable calendar year effect is numbered 3
          plot.trace.plots(exp(slot(object@kappa.log.error, 'value')[elements,,]) - 1, paste('Calendar Year Effect Error :', elements, sep=''))


      })


##' A generic function to plot and/or return a table of predicted age-to-age loss development factors (or link ratios).
##'
##' While the model estimates ultimate losses directly, comparisons of predicted to observed development factors can give the user a better feel for the model's adequacy.
##' Since the model is Bayesian, each development factor comes as a distribution.  Only the median, as a point estimate, are plotted/returned.
##'
##' The age-to-age factors are the ratios of the cumulative paid values at one period to the previous period.
##' Note that the median of products is not the product of medians, and thus it is possible (or rather likely) that age-to-age factors will not line up with age-to-ultimate factors (see \code{\link{tailFactor}}).
##'
##' @name lossDevelopmentFactors
##' @param object The object from which to plot and/or return loss development factors.
##' @param cex.text The \code{cex} value supplied to \code{text}. Adjusts the relative size of text.
##' @param linespace Adjusts the spacing between observed and predicted values.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a numeric matrix of plotted statistics.
##' @seealso \code{\link[=lossDevelopmentFactors,AnnualAggLossDevModelOutput-method]{lossDevelopmentFactors("AnnualAggLossDevModelOutput")}}
##'  \code{\link{tailFactor}}
##' @exportMethod lossDevelopmentFactors
setGenericVerif('lossDevelopmentFactors',
                function(object, cex.text=.77, linespace=0.5,  plot=TRUE)
            {
                if(!is.numeric(cex.text) || length(cex.text) != 1)
                    stop('"cex.text" must be numeric of length 1')

                if(!is.numeric(linespace) || length(linespace) != 1)
                    stop('"linespace" must be numeric of length 1')

                standardGeneric('lossDevelopmentFactors')
            })

##' A method to plot and/or return a table of predicted age-to-age loss development factors (or link ratios).
##'
##' @name  lossDevelopmentFactors,AnnualAggLossDevModelOutput-method
##' @param object The object of type /code{AnnualAggLossDevModelOutput} from which to plot and/or return the factors.
##' @param cex.text The \code{cex} value supplied to \code{text}. Adjusts the relative size of text.
##' @param linespace Adjusts the spacing between observed and predicted values.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting, but also returns a numeric matrix of plotted statistics.
##' @docType methods
##' @seealso \code{\link{lossDevelopmentFactors}}
##'  \code{\link{tailFactor}}
setMethod('lossDevelopmentFactors',
          signature(object='AnnualAggLossDevModelOutput'),
          function(object, cex.text, linespace, plot)
      {
          col.grey <- grey(0.45)

          K <- getTriDim(object@input)[1]
          obs.exp.years <- object@input@exposureYears
          pred.exp.years <- max(obs.exp.years) + 1:(object@input@totalExpYears - K)
          all.exp.years <- c(obs.exp.years, pred.exp.years)

          ey.type <- object@input@triangleType

          obs.ldf <- object@input@cumulatives[,-1] / object@input@cumulatives[,-K]
          obs.ldf <- rbind(obs.ldf, array(NA, c(length(pred.exp.years), dim(obs.ldf)[2])))

          ldf.pred <- array(NA, dim(obs.ldf))


          inc.pred <- slot(object@inc.pred, 'value')
          current.loss <- inc.pred[, 1, , ]
          for(j in 1:dim(ldf.pred)[2])
          {
              next.loss <- current.loss + inc.pred[, j + 1, , ]

              ldf.pred[,j] <- apply(next.loss / current.loss, 1, median, na.rm = TRUE)

              current.loss <- next.loss
          }

          if(ey.type=='py')
          {
              sub <- 'First Column is Link From Half to First'
          } else {
              sub <- 'First Column is Link From First to Second'
          }

          plot(x=range(1:(K-1)) + c(-0.5, 0.5),
               y=(range(obs.exp.years, pred.exp.years) + c(-0.5, 0.5))[c(2,1)],
               ylim=(range(obs.exp.years, pred.exp.years) + c(-0.5, 0.5))[c(2,1)],
               ylab=getExposureYearLabel(object@input),
               xlab='Development Year',
               cex.axis=1.25,
               cex.lab=1.25,
               type="n",
               font.main=1,
               cex.main=1.5,
               sub=sub)

          for (i in 1:length(all.exp.years))
          {
              for (j in 1:(K-1))
              {
                  text(x=j,
                       y=all.exp.years[i],
                       as.character(format(ldf.pred[i,j], digits=4, nsmall=3)),
                       font=1,
                       cex=cex.text,
                       col='black')

                  text(x=j,
                       y=all.exp.years[i]+linespace,
                       ifelse(is.na(obs.ldf[i,j]), '-', as.character(format(obs.ldf[i,j], digits=4, nsmall=3))),
                       font=1,
                       cex=cex.text,
                       col=col.grey)
              }
          }

      })


##' A generic function to plot the trace plots for select rate of decay values.
##'
##'
##' @name rateOfDecayTracePlot
##' @param object The object from which to generate the trace plots.
##' @param elements A numeric vector indicating for which elements to plot the trace.  Valid values are 2 through the number of columns in the observed triangle.  If NULL, values are selected automatically.
##' @param \dots Additional arguments used by methods.
##' @return NULL invisibly.  Only called for the side effect of plotting.
##' @seealso \code{\link[=rateOfDecayTracePlot,StandardAnnualAggLossDevModelOutput-method]{rateOfDecayTracePlot("StandardAnnualAggLossDevModelOutput")}}
##'  \code{\link{rateOfDecay}}
##' @exportMethod rateOfDecayTracePlot
setGenericVerif('rateOfDecayTracePlot',
                function(object, elements=NULL, ...)
                standardGeneric('rateOfDecayTracePlot'))

##' A generic function to generate the trace plots for select consumption path values.
##'
##' @name consumptionPathTracePlot
##' @param object The object from which to generate the trace plots.
##' @param elements A numeric vector indicating the elements for which to plot the trace.  Valid values are 1 through the number of development years (columns) in the observed triangle.  If NULL, values are selected automatically.
##' @param \dots Additional arguments used by methods.
##' @return NULL invisibly.  Only called for the side effect of plotting.
##' @seealso \code{\link[=consumptionPathTracePlot,StandardAnnualAggLossDevModelOutput-method]{consumptionTracePlot("StandardAnnualAggLossDevModelOutput")}}
##'  \code{\link{consumptionPath}}
##' @exportMethod consumptionPathTracePlot
setGenericVerif('consumptionPathTracePlot',
                function(object, elements=NULL, ...)
                standardGeneric('consumptionPathTracePlot'))


##' A generic function to plot and/or return the estimated consumption path vs development year time.
##'
##' At the heart of aggregate loss development models in \pkg{lossDev} is the consumption path.
##' The consumption path is (on a log scale) the trajectory of incremental payments absent calendar year effects and with exposure normalized to the first row.
##' Note that the measurement error term is (possibly) a skewed \eqn{t} and as such (possibly) has a non zero mean.   The consumption path is absent any such shifts due to skewness.
##' This is a generic function that allows for the retrieval and illustration of this consumption path.
##'
##' @name consumptionPath
##' @param object The object from which to plot and/or return the estimated consumption path.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns the plotted statistics.  Returned invisibly.
##' @seealso \code{\link[=consumptionPath,StandardAnnualAggLossDevModelOutput-method]{consumptionPath("StandardAnnualAggLossDevModelOutput")}}
##' \code{\link[=consumptionPath,BreakAnnualAggLossDevModelOutput-method]{consumptionPath("BreakAnnualAggLossDevModelOutput")}}
##' \code{\link{consumptionPathTracePlot}}
##' @exportMethod consumptionPath
setGenericVerif('consumptionPath',
                function(object, plot=TRUE)
                standardGeneric('consumptionPath'))

##' A generic function to plot and/or return the esimtated rate of decay vs development year time.
##'
##' The simplest definition of the rate of decay is the exponentiated first difference of the \link[=consumptionPath]{consumption path}.
##' This is a generic function to allow for the retrieval and illustration of the rate of decay.
##'
##' @name rateOfDecay
##' @param object The object from which to plot and/or return the estimated rate of decay.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns the plotted statistics.  Returned invisibly.
##' @seealso \code{\link[=rateOfDecay,StandardAnnualAggLossDevModelOutput-method]{rateOfDecay("StandardAnnualAggLossDevModelOutput")}}
##'  \code{\link[=rateOfDecay,BreakAnnualAggLossDevModelOutput-method]{rateOfDecay("BreakAnnualAggLossDevModelOutput")}}
##'  \code{\link{consumptionPath}}
##'  \code{\link{rateOfDecayTracePlot}}
##' @exportMethod rateOfDecay
setGenericVerif('rateOfDecay',
                function(object, plot=TRUE)
                standardGeneric('rateOfDecay'))

##' A generic function to plot and/or return the predicted tail factors for a specific attachment point.
##'
##' The tail factor is the ratio of the estimated ultimate loss to cumulative loss at some point in development time.
##' This is a generic function to allow for the retrieval and illustration the tail factor by exposure year.
##'
##' \bold{Note on \code{firstIsHalfReport} and \code{attachment}:} \code{firstIsHalfReport} refers to the first column of the triangle.
##' For policy year triangles, the first column is often referred to as a \dQuote{half-report}, the second column is called \dQuote{first-report}, the third column is called \dQuote{second-report}, etc.
##' If \code{firstIsHalfReport=TRUE}, then \code{tailFactor} will assume the triangle is arranged in such a way that the first column is the \dQuote{half-report}
##' and \code{attachment=1} indicates that the charted tail factor attaches at the cumulative loss through the second column.  If \code{firstIsHalfReport=FALSE},
##' then \code{attachment=1} indicates that the charted tail factor attaches at the cumulative loss through the first column.  Since \code{attachment} must be coercible to an integer,
##' it is impossible to plot half-to-ultimate tail factors; however, they are the first column in the returned matrix.
##'
##' \code{firstIsHalfReport} can be \code{NA} (the default)
##' if the exposure year type was specified to be one of \dQuote{policy year} or \dQuote{accident year} at the time the input object was constructed (see \code{\link{makeStandardAnnualInput}}
##' or \code{\link{makeBreakAnnualInput}}).  An exposure year type of \dQuote{policy year} corresponds to \code{firstIsHalfReport=TRUE},
##' and an exposure year type of \dQuote{accident year} corresponds to \code{firstIsHalfReport=FALSE}.  Setting \code{firstIsHalfReport} to a non-missing value will override this default.
##'
##' If \code{expYearRange} is \dQuote{fullyObs}, then only exposure years with a non missing value in the first column will be plotted.
##'
##' @name tailFactor
##' @param object The object from which to plot the predicted tail factors and return tail factors for \emph{all} attachment points.
##' @param attachment An integer value specifying the attachment point for the tail.  Must be at least 1. See Details for more information.
##' @param useObservedValues A logical value.  If \code{TRUE}, observed values are substituted for predicted values whenever possible in the calculation.  If \code{FALSE}, only predicted values are used.
##' @param firstIsHalfReport A logical value or \code{NA}.  See Details for more info.
##' @param finalAttachment An integer value must be at least 1. Default value is \code{attachment}.  A call to \code{tailFactor} will return (invisibly) a matrix of tail factors through this value.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @param expYearRange Either a range of years (for example c(1995, 2006)) or one of the keywords \dQuote{all} or \dQuote{fullyObs}.
##' @return Mainly called for the side effect of plotting.
##' @seealso \code{\link[=tailFactor,StandardAnnualAggLossDevModelOutput-method]{tailFactor("StandardAnnualAggLossDevModelOutput")}}
##'  \code{\link[=tailFactor,BreakAnnualAggLossDevModelOutput-method]{tailFactor("BreakAnnualAggLossDevModelOutput")}}
##' @exportMethod tailFactor
setGenericVerif('tailFactor',
                function(object, attachment, useObservedValues=FALSE, firstIsHalfReport=NA, finalAttachment=attachment, plot=TRUE, expYearRange='all')
            {
                if(!is.numeric(attachment))
                    stop('"attachment" must be numeric')

                if(!identical(length(attachment), as.integer(1)))
                    stop('"attachment" must be of length 1')

                if(is.na(attachment))
                    stop('"attachment" cannot be NA')

                if(as.integer(attachment) != attachment)
                    stop('"attachment" must be coercible to an integer')

                if(attachment <= 0)
                    stop('"attachment" must be at least 1')

                if(!is.numeric(finalAttachment))
                    stop('"finalAttachment" must be numeric')

                if(!identical(length(finalAttachment), as.integer(1)))
                    stop('"finalAttachment" must be of length 1')

                if(is.na(finalAttachment))
                    stop('"finalAttachment" cannot be NA')

                if(as.integer(finalAttachment) != finalAttachment)
                    stop('"finalAttachment" must be coercible to an integer')

                if(finalAttachment <= 0)
                    stop('"finalAttachment" must be at least 1')

                if(finalAttachment < attachment)
                    stop('"finalAttachment" must be at least equal to "attachment"')


                if(is.character(expYearRange))
                {
                    if(length(expYearRange) != 1)
                        stop('"expYearRange" must be of length one if it is a character')
                    if(expYearRange != 'all' && expYearRange != 'fullyObs')
                        stop('"expYearRange" must be one of "all" or "fullyObs" if it is supplied as a character')
                    ##if(expYearRange == 'all')
                        ##expYearRange <- range(exp.years)
                    ##else
                        ##expYearRange <- range(exp.years[which(!is.na(cumulatives[,1]))])
                } else {

                    if(!all(as.integer(expYearRange) == expYearRange))
                        stop('"expYearRange" must be supplied as an integer')
                    if(length(expYearRange) != 2)
                        stop('"expYearRange" must have length 2')
                    ##if(max(exp.years) < max(expYearRange) || min(exp.years) > min(expYearRange))
                        ##stop('"expYearRange" must be a subset of the actual exposure years')
                }



                standardGeneric('tailFactor')
            })

##' A generic function to plot and/or return the posterior number of knots.
##'
##' The \link[=consumptionPath]{consumption path} (or calendar year effect and exposure growth adjusted log incremental payments) is modeled as a linear spline.
##' The number of knots (or places where the spline changes slope) in this spline is endogenous to the model and estimated by way of Reversible Jump Markov Chain Monte Carlo simulation.
##'
##' @name numberOfKnots
##' @param object The object from which to plot the number of knots.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns statics on the number of knots.  Returned invisibly.
##' @seealso \code{\link{consumptionPath}}
##'  \code{\link[=numberOfKnots,StandardAnnualAggLossDevModelOutput-method]{numberOfKnots("StandardAnnualAggLossDevModelOutput")}}
##'  \code{\link[=numberOfKnots,BreakAnnualAggLossDevModelOutput-method]{numberOfKnots("BreakAnnualAggLossDevModelOutput")}}
##' @exportMethod numberOfKnots
setGenericVerif('numberOfKnots',
                function(object, plot=TRUE)
            {
                standardGeneric('numberOfKnots')
            })

##' A generic function to plot autocorrelations found in the \acronym{MCMC} samples for select parameters.
##'
##' Chains with high autocorrelation require a longer burnin and more samples to fully explore the parameter space.
##'
##' @name mcmcACF
##' @param object The object from which to plot autocorrelations.
##' @return Called for the side effect of plotting.
##' @seealso \code{\link[=mcmcACF,StandardAnnualAggLossDevModelOutput-method]{mcmcACF("StandardAnnualAggLossDevModelOutput")}}
##'  \code{\link[=mcmcACF,BreakAnnualAggLossDevModelOutput-method]{mcmcACF("BreakAnnualAggLossDevModelOutput")}}
##' @exportMethod mcmcACF
setGenericVerif('mcmcACF',
                function(object)
            {
                standardGeneric('mcmcACF')
            })
