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
##' @include NodeOutput.R
##' @include AnnualAggLossDevModelOutput.R
NULL

##' The final output class for all standard aggregate annual models.
##'
##' \code{StandardAnnualAggLossDevModelOutput} is the final output class for all  standard aggregate annual model objects.
##' Currenly, only the slot \dQuote{input} is allowed to be a non-model node.  All other nodes should be the exact name of some settable node in the model.
##' This is because \code{getModelOutputNodes} currently looks at the slot names to determine what values to set; only slot \dQuote{input} is known to be a slot other than a settable node.
##' This class is derived from \code{AnnualAggLossDevModelOutput}.
##' @name StandardAnnualAggLossDevModelOutput-class
##' @docType class
##' @seealso \code{\linkS4class{AnnualAggLossDevModelOutput}}
setClass(
         'StandardAnnualAggLossDevModelOutput',
         representation(S='NodeOutput'),
         contains='AnnualAggLossDevModelOutput')



##' A method to plot and/or return the estimated consumption path vs development year time for standard models.
##'
##' At the heart of aggregate loss development models in \pkg{lossDev} is the consumption path.
##' The consumption path is (on a log scale) the trajectory of incremental payments absent any calendar year effects and wisht exposure normalized to the first row.
##' Note that the measurement error term is (possibly) a skewed \eqn{t} and as such (possibly) has a non zero mean.   The consumption path is absent any such shifts due to skewness.
##' This is a method that allows for the retrieval and illustration of this consumption path.
##'
##' The standard model has a common consumption path for all exposure years.
##'
##' Because the model is Bayesian, the estimated consumption path comes as a distribution; only the median is plotted and/or returned.
##'
##' @name consumptionPath,StandardAnnualAggLossDevModelOutput-method
##' @param object The object from which to plot and/or return the estimated consumption path.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns the plotted statistics.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{consumptionPath}}
##'  \code{\link[=consumptionPath,BreakAnnualAggLossDevModelOutput-method]{consumptionPath("BreakAnnualAggLossDevModelOutput")}}
setMethod('consumptionPath',
          signature(object='StandardAnnualAggLossDevModelOutput'),
          function(object, plot)
      {

          cp <- object@S@median

          if(plot)
              plot(cp,
                   xlab='Development Year',
                   ylab="Calendar Year-Effect Adjusted Log Incremenal Payments",
                   type="b",
                   cex.axis=1.25,
                   cex.lab=1.25)

          return(invisible(cp))


      })

##' A method to plot and/or return the esimtated rate of decay vs development year time for standard models.
##'
##' The simplest definition of the rate of decay is the exponentiated first difference of the \link[=consumptionPath]{consumption path}.
##' The standard model has a common rate of decay for all exposure years.
##' This is a method to allow for the retrieval and illustration of the rate of decay.
##'
##' Because the model is Bayesian, the estimated rate of decay comes as a distribution; only the median is plotted and/or returned.
##'
##' @name rateOfDecay,StandardAnnualAggLossDevModelOutput-method
##' @param object The object from which to plot and/or return the estimated rate of decay.
##' @param plot A logical value. If \code{TRUE}, then the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns the plotted statistics.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{rateOfDecay}}
##'  \code{\link[=rateOfDecay,BreakAnnualAggLossDevModelOutput-method]{rateOfDecay("BreakAnnualAggLossDevModelOutput")}}
##'  \code{\link{consumptionPath}}
setMethod('rateOfDecay',
          signature(object='StandardAnnualAggLossDevModelOutput'),
          function(object, plot)
      {

          K <- getTriDim(object@input)[1]
          S <- slot(object@S, 'value')
          delta.obs <- apply(exp(S[-1,,] - S[-K,,]) - 1, 1, median)
          delta.tail <- object@delta.tail@median

          last.non.zero.payment <- object@input@lastNonZeroPayment

          ans <- array(NA,
                       c(dim(delta.tail)[1], length(delta.obs)))

          for(i in 1:(K-1))
              ans[,i] <- delta.obs[i]

          ans <- cbind(ans, delta.tail)
          for(i in 1:length(last.non.zero.payment))
          {
              ans[,1:dim(ans)[2] + 1 == last.non.zero.payment[i] + 1] <- -1
              ans[,1:dim(ans)[2] + 1 > last.non.zero.payment[i] + 1] <- NA
          }

          dimnames(ans) <- list(1:dim(ans)[1] + object@input@exposureYears[1] - 1,NULL)



          f.plot <- function()
          {
              plot(
                   x=c(2, max(last.non.zero.payment)),
                   y=range(delta.obs, delta.tail, 0),
                   xlab='Development Year',
                   ylab="Rate of Decay",
                   type="n",
                   cex.axis=1.25,
                   cex.lab=1.25
                   )

              lines(
                    x=2:K,
                    y=delta.obs,
                    lwd=2,
                    type='b')

              for(i in 1:length(last.non.zero.payment))
              {
                  lines(
                        x=(K + 1):last.non.zero.payment[i],
                        y=delta.tail[i,(K + 1):last.non.zero.payment[i] - K],
                        lwd=2,
                        col='black',
                        type="b",
                        lty=3,
                        pch=3)
              }

              abline(a=0,b=0,col='black',lwd=2,lty='dashed')
          }

          f.legend <- function()
          {

              legend('center',
                     c("Estimated","Projected"),
                     col=c('black','black'),
                     lwd=c(2,2),
                     lty=c(1,3),
                     pch=c(1,3),
                     horiz=TRUE,
                     bty='n',
                     xpd=NA)
          }

          if(plot)
              plot.top.bottom(f.plot, f.legend)

          return(invisible(ans))

      })

##' A method to plot and/or return the predicted tail factors for a specific attachment point.
##'
##' The tail factor is the ratio of the estimated ultimate loss to cumulative loss at some point in development time.
##' This is a method to allow for the retrieval and illustration of the tail factor by exposure year.
##'
##' Because the model is Bayesian, each tail factor comes as a distribution.  To ease graphical interpretation, only the median for each factor is plotted/returned.
##' See for more details \code{\link{tailFactor}}.
##'
##'
##' @name tailFactor,StandardAnnualAggLossDevModelOutput-method
##' @param object The object from which to plot the predicted tail factors and return tail factors for \emph{all} attachment points.
##' @param attachment An integer value specifying the attachment point for the tail.  Must be at least 1. See Details for more info.
##' @param useObservedValues A logical value.  If \code{TRUE}, observed values are substituted for predicted values whenever possible in the calculation.  If \code{FALSE}, only predicted values are used.
##' @param firstIsHalfReport A logical value or \code{NA}.  See Details for more information.
##' @param finalAttachment An integer value must be at least 1 default value is \code{attachment}.  A call to \code{tailFactor} returns (invisibly) a matrix of tail factors through this value.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @param expYearRange Either a range of years (for example c(1995, 2006)) or one of the keywords \dQuote{all} or \dQuote{fullyObs}.
##' @return Mainly called for the side effect of plotting.  Also returns tail factors for \emph{all} attachment points through \code{finalAttachment}.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{tailFactor}}
##'  \code{\link[=tailFactor,BreakAnnualAggLossDevModelOutput-method]{tailFactor("BreakAnnualAggLossDevModelOutput")}}
setMethod('tailFactor',
          signature(object='StandardAnnualAggLossDevModelOutput'),
          function(object, attachment, useObservedValues, firstIsHalfReport, finalAttachment, plot, expYearRange)
      {

          ey.type <- object@input@triangleType
          if(is.na(firstIsHalfReport))
          {
              if(ey.type == 'py')
                  firstIsHalfReport <- TRUE
              else if(ey.type == 'ay')
                  firstIsHalfReport <- FALSE
              else
                  stop('"firstIshalfReport" is missing and the exposure year type is not one of "policy year" or "accident year." Please set a value for "firstIsHalfReport".')
          }

          attachment.adj <- ifelse(firstIsHalfReport, attachment + 1, attachment)

          ##we don't have to test attachment nor finalAttachment because the generic function does that for us
          n.columns <- max(attachment.adj, finalAttachment)

          inc.pred <- slot(object@inc.pred, 'value')
          total.col <- dim(inc.pred)[2]
          total.rows <- dim(inc.pred)[1]
          total.exp.years <- object@input@exposureYears[1] + 1:total.rows - 1


          if(is.character(expYearRange))
          {
              ##if(length(expYearRange) != 1)
              ##stop('"expYearRange" must be of length one if it is a character')
              ##if(expYearRange != 'all' && expYearRange != 'fullyObs')
              ##stop('"expYearRange" must be one of "all" or "fullyObs" if it is supplied as a character')
              if(expYearRange == 'all')
                  expYearRange <- range(total.exp.years)
              else
                  expYearRange <- range(object@input@exposureYears[which(!is.na(object@input@cumulatives[,1]))],
                                        total.exp.years[total.exp.years > max(object@input@exposureYears)])
          } else {

              ##if(!all(as.integer(expYearRange) == expYearRange))
              ##stop('"expYearRange" must be supplied as an integer')
              ##if(length(expYearRange) != 2)
              ##stop('"expYearRange" must have length 2')
              if(max(total.exp.years) < max(expYearRange) || min(total.exp.years) > min(expYearRange))
                  stop('"expYearRange" must be a subset of the actual exposure years')
          }

          if(n.columns > total.col)
              stop('either "attachment" or "finalAttachment" is set larger than the number of total development years.  Set these to a smaller value or re-estiamte the triangle with additional predicted development years.')

          if(useObservedValues)
          {
              K <- getTriDim(object@input)[1]
              inc <- object@input@incrementals
              inc.pred[1:K,1:K,,][!is.na(inc)] <- inc[!is.na(inc)]
          }

          if(firstIsHalfReport)
              col.names <- paste(c('half', 1:(n.columns-1)), 'ult', sep='-to-')
          else
              col.names <- paste(1:n.columns, 'ult', sep='-to-')

          tail.matrix <- array(NA, c(total.rows, n.columns), list(total.exp.years, col.names))
          ult <- apply(inc.pred, c(1,3,4), sum)

          current.loss <- array(0, dim(inc.pred)[c(1,3,4)])
          for(j in 1:n.columns)
          {
              current.loss <- current.loss + inc.pred[,j,,]

              tail.matrix[,j] <- apply(ult / current.loss, c(1), median, na.rm=TRUE)

          }

          if(plot)
          {
              expYearRange.seq <- seq(expYearRange[1], expYearRange[2])

              plot(x=range(total.exp.years),
                   y=range(tail.matrix[as.character(expYearRange.seq), attachment.adj]),
                   xlab=getExposureYearLabel(object@input),
                   ylab=paste(attachment,"th to Ultimate Tail Factor",sep=""),
                   type='n',
                   cex.axis=1.25,
                   cex.lab=1.25)

              lines(x=expYearRange.seq,
                    y=tail.matrix[as.character(expYearRange.seq),attachment.adj],
                    type='b')


          }



          return(invisible(tail.matrix))
      })


##' A method to plot and/or return the posterior number of knots.
##'
##'
##' @name numberOfKnots,StandardAnnualAggLossDevModelOutput-method
##' @param object The object from which to plot the number of knots.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a named vector with names equal to the number of knots and values equal to the density.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{consumptionPath}}
##'  \code{\link{numberOfKnots}}
##'  \code{\link[=numberOfKnots,BreakAnnualAggLossDevModelOutput-method]{numberOfKnots("BreakAnnualAggLossDevModelOutput")}}
setMethod('numberOfKnots',
          signature(object='StandardAnnualAggLossDevModelOutput'),
          function(object, plot=TRUE)
      {

         k <- slot(object@k, 'value')
         k.vector <- as.vector(k)

         jd <- getJagsData(object@input)

         x.number.of.knots <- seq(from=0, to=jd$number.of.knots.ubound)
         r <- jd$mu.number.of.knots.prior[1]
         p <- 1 / (1 + jd$mu.number.of.knots.prior[2])


         plot.data <- array(NA,
                            c(length(x.number.of.knots), 3),
                            list(NULL, c('NumberOfKnots', 'Prior', 'Posterior')))

         k <- object@k@get.value.env$get.value()
         N <- prod(dim(k))
         for(i in seq_along(x.number.of.knots))
         {
             x.i <- x.number.of.knots[i]
             plot.data[i, 'NumberOfKnots'] <- x.i
             plot.data[i, 'Prior'] <- gamma(r + x.i) / factorial(x.i) / gamma(r) * (1 - p) ^ r * p ^ x.i
             plot.data[i, 'Posterior'] <- sum( k == x.i) / N
         }

         plot.data[,'Posterior'] <-   plot.data[,'Posterior'] / sum(plot.data[,'Posterior'])

         bar <- function(x, y, offset=c('L', 'R'), col=c('grey', 'black'))
         {
             epsilon <- 0.4

             offset <- match.arg(offset)
             col <- match.arg(col)

             if(offset == 'L')
                 epsilon <- epsilon * -1

                 rect(xleft = x + epsilon,
                      ybottom = 0,
                      xright = x,
                      ytop = y,
                      density = -1,
                      col=col,
                      border=NA)

         }


         f.plot <- function()
         {

             plot(x=range(plot.data[,'NumberOfKnots']) + c(-0.75, 0.75),
                  y=range(plot.data[,c('Prior', 'Posterior')]),
                  type='n',
                  cex.axis=1.25,
                  cex.lab=1.25,
                  xlab='Number of Knots',
                  ylab='Relative Frequency')

             for(i in 1:dim(plot.data)[1])
             {
                 bar(x=plot.data[i,'NumberOfKnots'],
                     y=plot.data[i,'Prior'],
                     col='grey',
                     offset='L')

                 bar(x=plot.data[i,'NumberOfKnots'],
                     y=plot.data[i,'Posterior'],
                     col='black',
                     offset='R')
             }

         }

        f.legend <- function()
          {

              legend('center',
                     c("Prior","Posterior"),
                     col=c('grey','black'),
                     pch=15,
                     pt.cex=3,
                     horiz=TRUE,
                     bty='n',
                     xpd=NA)
          }

          if(plot)
              plot.top.bottom(f.plot, f.legend)


         ans <- plot.data[, 'Posterior']
         names(ans) <- plot.data[, 'NumberOfKnots']

         return(invisible(ans))

      })


##' A method to generate the trace plots for select rate of decay values.
##'
##'
##' @name rateOfDecayTracePlot,StandardAnnualAggLossDevModelOutput-method
##' @param object The object of type \code{StandardAnnualAggLossDevModelOutput} from which to generate the trace plots.
##' @param elements A numeric vector indicating for which elements to plot the trace.  Valid values are 2 through the number of columns in the observed triangle.  If NULL, values are selected automatically.
##' @param \dots Additional arguments used by other methods.  Not utilized by this method.
##' @return NULL invisibly.  Only called for the side effect of plotting.
##' @docType methods
##' @seealso \code{\link{rateOfDecayTracePlot}}
##'  \code{\link{rateOfDecay}}
setMethod('rateOfDecayTracePlot',
          signature(object='StandardAnnualAggLossDevModelOutput'),
          function(object, elements=NULL, ...)
      {
          lb <- 2
          ub <- getTriDim(object@input)[2]

          if(is.null(elements))
          {
              elements <- c(lb, floor((lb + ub) / 2), ub)
          } else {

              if(!is.numeric(elements) || !is.vector(elements))
                  stop('"elements" must either be "NULL" or a numeric vector')

              if(any(elements != as.integer(elements)))
                  stop('"elements" must be coercible to integers')

              if(any(elements > ub) || any(elements < lb))
                  stop(paste('"elements" must be at most the number of columns in the observed triangle (', ub,') and at least 2', sep=''))
          }

          K <- ub
          S <- slot(object@S, 'value')
          delta.obs <- exp(S[-1,,] - S[-K,,]) - 1
          plot.trace.plots(delta.obs[elements-1,,], paste('Rate of Decay :', elements, sep=''))
      })


##' A method to generate the trace plots for select consumption path values.
##'
##'
##' @name consumptionPathTracePlot,StandardAnnualAggLossDevModelOutput-method
##' @param object The object of type \code{StandardAnnualAggLossDevModelOutput} from which to generate the trace plots.
##' @param elements A numeric vector indicating the elements for which to plot the trace.  Valid values are 1 through the number of development years (columns) in the observed triangle.  If NULL, values are selected automatically.
##' @param \dots Additional arguments used by other methods.  Not utilized by this method.
##' @return NULL invisibly.  Only called for the side effect of plotting.
##' @docType methods
##' @seealso \code{\link{consumptionPathTracePlot}}
##'  \code{\link{consumptionPath}}
setMethod('consumptionPathTracePlot',
          signature(object='StandardAnnualAggLossDevModelOutput'),
          function(object, elements=NULL, ...)
      {
          lb <- 1
          ub <- getTriDim(object@input)[2]

          if(is.null(elements))
          {
              elements <- c(lb, floor((lb + ub) / 2), ub)
          } else {

              if(!is.numeric(elements) || !is.vector(elements))
                  stop('"elements" must either be "NULL" or a numeric vector')

              if(any(elements != as.integer(elements)))
                  stop('"elements" must be coercible to integers')

              if(any(elements > ub) || any(elements < lb))
                  stop(paste('"elements" must be at most the number of columns in the observed triangle (', ub,') and at least 1', sep=''))
          }

          plot.trace.plots(slot(object@S, 'value')[elements,,], paste('Consumption Path :', elements, sep=''))
      })

##' A method to plot autocorrelations found in the \acronym{MCMC} samples for select parameters.
##'
##' Chains with high autocorrelation require a longer burnin and more samples to fully explore the parameter space.
##'
##' @name mcmcACF,StandardAnnualAggLossDevModelOutput-method
##' @param object The object from which to plot autocorrelations.
##' @return Called for the side effect of plotting.
##' @docType methods
##' @seealso \code{\link{mcmcACF}}
##'  \code{\link[=mcmcACF,BreakAnnualAggLossDevModelOutput-method]{mcmcACF("BreakAnnualAggLossDevModelOutput")}}
##' @TODO Add option to plot other values.  Currently only plots \dQuote{First Rate Of Decay,} \dQuote{First Calendar Year Effect Error,} and \dQuote{First Exposure Year Growth.}
setMethod('mcmcACF',
          signature(object='StandardAnnualAggLossDevModelOutput'),
          function(object)
      {
          ##Autocorrelation

          n.chain <- dim(slot(object@df, 'value'))['chain']
          K <- getTriDim(object@input)[1]

          op <- par(no.readonly=TRUE)
          on.exit(par(op))

          par(mfcol=c(3,n.chain))

          S <- slot(object@S, 'value')
          delta <- S[-1, , ] - S[-K, , ]




          for(i in 1:n.chain)
          {
              acf(delta[1,,i],
                  cex.axis=1.25,
                  cex.lab=1.25,
                  main=paste('First Rate Of Decay\nChain:', i))

              acf(slot(object@kappa.log.error, 'value')[3,,i],
                  cex.axis=1.25,
                  cex.lab=1.25,
                  main=paste('First Calendar Year Effect Error\nChain:', i))

              acf(slot(object@eta, 'value')[2,,i],
                  cex.axis=1.25,
                  cex.lab=1.25,
                  main=paste('First Exposure Year Growth\nChain:', i))
          }





      })
