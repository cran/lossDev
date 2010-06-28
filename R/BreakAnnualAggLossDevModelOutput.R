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
##' Currenly only the slot \dQuote{input} is allowed to be a non-model node.  All other nodes should be the exact name of some settable node in the model.
##' This is because \code{getModelOutputNodes} currently looks the slot names to determine what values to set; only slot \dQuote{input} is known to be a slot other than a settable node.
##' This class is derived from \code{AnnualAggLossDevModelOutput}
##' @name BreakAnnualAggLossDevModelOutput-class
##' @docType class
##' @seealso \code{\linkS4class{AnnualAggLossDevModelOutput}}
setClass(
         'BreakAnnualAggLossDevModelOutput',
         representation(
                        R='NodeOutput',
                        inc.brk.pre='NodeOutput',
                        inc.brk.post='NodeOutput',
                        first.row.in.post='NodeOutput'),
         contains='AnnualAggLossDevModelOutput')

##' A method to plot and/or return the estimated consumption path vs development year time for break models.
##'
##' At the heart of aggregate loss development models in \pkg{lossDev} is the consumption path.
##' The consumption path is (on a log scale) the trajectory of incremental payments absent any calendar year effects and normalized to the first row.
##' Note that the measurement error term is (possibly) a skewed \eqn{t} and as such (possibly) has a non zero mean.   The consumption path is absent any such shifts in due to skewness.
##' This is a method to allow for the retrieval and illustration of this consumption path.
##'
##' The break model has a two consumption paths: exposure years prior to the structural break follow one path and exposure years after the break follow another.
##' The slope of the consumption path for exposure years prior to the break is used to extend the consuption path for exposure years post break to the end of the triangle.
##'
##' Because the model is Bayesian, the estimated consumption paths come as distributions; only the medians are plotted and/or returned.
##'
##' @name consumptionPath,BreakAnnualAggLossDevModelOutput-method
##' @param object The object from which to plot and/or return the estimated consumption paths.
##' @param plot A logical value. If \code{TRUE}, then the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns the plotted statistics.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{consumptionPath}}
##' @seealso \code{\link[=consumptionPath,StandardAnnualAggLossDevModelOutput-method]{consumptionPath("StandardAnnualAggLossDevModelOutput")}}
setMethod('consumptionPath',
          signature(object='BreakAnnualAggLossDevModelOutput'),
          function(object, plot)
      {


          jd <- getJagsData(object@input)
          K <- getTriDim(object@input)[1]
          K.trim <- jd$K.trim

          cp <- object@R@median

          dimnames(cp) <- list(c('Pre-Break', 'Post-Break'), NULL)
          #cp[2,1:K.trim] <- NA

          if(plot)
          {
              f.plot <- function()
              {
                  plot(x=c(1,K),
                       y=range(cp, na.rm=TRUE),
                       xlab='Development Year',
                       ylab='Calendar Year-Effect Adjusted Log Incremenal Payments',
                       type='n',
                       cex.axis=1.25,
                       cex.lab=1.25)

                  lines(
                        x=1:K,
                        y=cp[1,],
                        type='b',
                        col='gray',
                        lwd=3)


                  lines(
                        x=1:K.trim,
                        y=cp[2,1:K.trim],
                        type='b',
                        col='black',
                        lwd=2)
              }

              f.legend <- function()
              {

                  legend('center',
                         c("Pre-Structrural Break","Post-Structural Break"),
                         col = c('gray','black'),
                         lwd=c(3,2),
                         horiz=TRUE,
                         bty='n',
                         xpd=NA)
              }

              plot.top.bottom(f.plot, f.legend)

          }

          return(invisible(cp))


      })



##' A method to plot and/or return the esimtated rate of decay vs development year time for break models.
##'
##' The simplest definition of the rate of decay is the exponentiated first difference of the \link[=consumptionPath]{consumption path}.
##' The break model has two rates of decay.  One which applies to exposure years prior to a structural break.  And another which applies after the break.
##' This is a method to allow for the retrieval and illustration of these rates of decay.
##'
##' Because the model is Bayesian, the estimated rates of decay come as distributions; only the medians are plotted and/or returned.
##'
##'
##' @name rateOfDecay,BreakAnnualAggLossDevModelOutput-method
##' @param object The object from which to plot and/or return the estimated rate of decay.
##' @param plot A logical value. If \code{TRUE}, then the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns the plotted statistics.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{rateOfDecay}}
##' @seealso \code{\link[=rateOfDecay,StandardAnnualAggLossDevModelOutput-method]{rateOfDecay("StandardAnnualAggLossDevModelOutput")}}
##' @seealso \code{\link{consumptionPath}}
setMethod('rateOfDecay',
          signature(object='BreakAnnualAggLossDevModelOutput'),
          function(object, plot)
      {
          jd <- getJagsData(object@input)
          K <- getTriDim(object@input)[1]
          K.trim <- jd$K.trim


          R <- slot(object@R, 'value')
          pre.delta.obs <- apply(exp(R[1,-1,,] - R[1,-K,,]) - 1, 1, median)
          post.delta.obs <- apply(exp(R[2,-1,,] - R[2,-K,,]) - 1, 1, median)
          delta.tail <- object@delta.tail@median
          last.non.zero.payment <- object@input@lastNonZeroPayment

          ans <- rbind(pre.delta.obs,
                       post.delta.obs)

          dimnames(ans) <- list(c('Pre-Break', 'Post-Break'), NULL)



          f.plot <- function()
          {


              plot(
                   x=c(2,max(last.non.zero.payment)),
                   y=range(pre.delta.obs, post.delta.obs, delta.tail, na.rm=TRUE),
                   xlab='Development Year',
                   ylab='Rate of Decay',
                   type='n',
                   cex.axis=1.25,
                   cex.lab=1.25)

              lines(
                    x=2:K,
                    y=pre.delta.obs,
                    type='b',
                    pch=20,
                    col='gray',
                    lwd=2)

              lines(
                    x=2:K,
                    y=post.delta.obs,
                    pch=1,
                    type='b',
                    col='black',
                    lwd=2)

              for(i in 1:(object@input@totalExpYears))
              {
                  lines(
                        x=(K+1):last.non.zero.payment[i],
                        delta.tail[i,1:(last.non.zero.payment[i] - K)],
                        pch=1,
                        lwd=2,
                        col='grey',
                        type='b',
                        lty=3)
              }

              abline(a=0,b=0,lwd=2,col='black',lty='dashed')

          }

          f.legend <- function()
          {
              legend('center',
                     c(
                       'Pre-Structrural\nBreak',
                       'Post-Structural\nBreak',
                       'Pre/Post-Structrural\nBreak'),
                     col=c('gray','black','grey'),
                     pch=c(20,1,1),
                     lwd=c(2,2,2),
                     lty=c(1,1,3),
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
##' @name tailFactor,BreakAnnualAggLossDevModelOutput-method
##' @param object The object from which to plot the predicted tail factors and return tail factors for \emph{all} attachment points.
##' @param attachment An integer value specifying the attachment point for the tail.  Must be at least 1. See Details for more info.
##' @param useObservedValues A logical value.  If \code{TRUE}, observed values are substituted for predicted values whenever possible in the calculation.  If \code{FALSE}, only predicted values are used.
##' @param firstIsHalfReport A logical value or \code{NA}.  See Details for more information.
##' @param finalAttachment An integer value must be at least 1 default value is \code{attachment}.  A call to \code{tailFactor} returns (invisibly) a matrix of tail factors through this value.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @param expYearRange Either a range of years (for example c(1995, 2006)) or one of the keywords \dQuote{all} or \dQuote{fullyObs}.
##' @return Mainly called for the side effect of plotting.  Also returns tail factors for \emph{all} attachment points through \code{finalAttachment}.  See Details. Returned invisibly.
##' @docType methods
##' @seealso \code{\link{tailFactor}}
##' @seealso \code{\link[=tailFactor,StandardAnnualAggLossDevModelOutput-method]{tailFactor("StandardAnnualAggLossDevModelOutput")}}
setMethod('tailFactor',
          signature(object='BreakAnnualAggLossDevModelOutput'),
          function(object, attachment, useObservedValues, firstIsHalfReport, finalAttachment, plot, expYearRange)
      {

          ey.type <- object@input@triangleType
          if(is.na(firstIsHalfReport))
          {
              if(ey.type == 'py')
                  firstIsHalfReport <- TRUE
              else if(ey.type == 'ay')
                  firstIsHalfReport <- TRUE
              else
                  stop('"firstIsHalfReport" is missing and the exposure year type is not one of "policy year" or "accident year." Please set a value for "firstIsHalfReport".')
          }

          attachment.adj <- ifelse(firstIsHalfReport, attachment + 1, attachment)

          ##we don't have to test attachment nor finalAttachment because the generic function does that for us
          n.columns <- max(attachment.adj, finalAttachment)

          inc.pred <- list()
          inc.pred[['Actual']] <- slot(object@inc.pred, 'value')
          inc.pred[['AsIfPreBreak']] <-  slot(object@inc.brk.pre, 'value')
          inc.pred[['AsIfPostBreak']] <- slot(object@inc.brk.post, 'value')

          total.col <- dim(inc.pred[[1]])[2]
          total.rows <- dim(inc.pred[[1]])[1]
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
              for(i in 1:length(inc.pred))
                  inc.pred[[i]][1:K,1:K,,][!is.na(inc)] <- inc[!is.na(inc)]
          }

          if(firstIsHalfReport)
              col.names <- paste(c('half', 1:(n.columns-1)), 'ult', sep='-to-')
          else
              col.names <- paste(1:n.columns, 'ult', sep='-to-')

          tail.list <- list()
          for(i in names(inc.pred))
          {
              tail.list[[i]] <- array(NA, c(total.rows, n.columns), list(total.exp.years, col.names))
              ult <- apply(inc.pred[[i]], c(1,3,4), sum)

              current.loss <- array(0, dim(inc.pred[[i]])[c(1,3,4)])
              for(j in 1:n.columns)
              {
                  current.loss <- current.loss + inc.pred[[i]][,j,,]

                  tail.list[[i]][,j] <- apply(ult / current.loss, c(1), median, na.rm=TRUE)

              }
          }

          if(plot)
          {

              plot.f <- function()
              {

                  expYearRange.seq <- seq(expYearRange[1], expYearRange[2])

                  y.range <- NA
                  for(i in names(tail.list))
                      y.range <- range(y.range, tail.list[[i]][as.character(expYearRange.seq), attachment.adj], na.rm=TRUE)

                  plot(x=range(total.exp.years),
                       y=y.range,
                       xlab=getExposureYearLabel(object@input),
                       ylab=paste(attachment,"th to Ultimate Tail Factor",sep=""),
                       type="n",
                       cex.axis=1.25,
                       cex.lab=1.25)





                  lines(x=expYearRange.seq,
                        y=tail.list$AsIfPreBreak[as.character(expYearRange.seq),attachment.adj],
                        type="b",
                        col="gray",
                        lty=2,
                        lwd=3)

                  lines(x=expYearRange.seq,
                        y=tail.list$AsIfPostBreak[as.character(expYearRange.seq),attachment.adj],
                        type="b",
                        col="gray",
                        lty=1,
                        lwd=3)

                  lines(x=expYearRange.seq,
                        y=tail.list$Actual[as.character(expYearRange.seq),attachment.adj],
                        type="b",
                        col='black',
                        lty=1,
                        lwd=2)



              }
              legend.f <- function()
              {
                  legend('center',
                         c("Actual","As if Post-Break","As if Pre-Break"),
                         col = c('black','gray','gray'),
                         lty=c(1,1,2),
                         lwd=c(2,2,2),
                         horiz=TRUE,
                         bty='n',
                         xpd=NA)
              }

              plot.top.bottom(plot.f, legend.f)
          }

          return(invisible(tail.list))
      })

##' A method to plot and/or return the posterior number of knots.
##'
##' The break model has to consumption paths.  This method will plot the number of knots for each path.
##'
##' @name numberOfKnots,BreakAnnualAggLossDevModelOutput-method
##' @param object The object from which to plot the number of knots.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a list of length 2 with each element containing a named vector with names equal to the number of knots and values equal to the density.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{consumptionPath}}
##' @seealso \code{\link{numberOfKnots}}
##' @seealso \code{\link[=numberOfKnots,StandardAnnualAggLossDevModelOutput-method]{numberOfKnots("StandardAnnualAggLossDevModelOutput")}}
setMethod('numberOfKnots',
          signature(object='BreakAnnualAggLossDevModelOutput'),
          function(object, plot=TRUE)
      {

          jd <- getJagsData(object@input)

          k <- slot(object@k, 'value')

          plot.data <- list()

          k.nice.names <- c('Pre-Structural Break', 'Post-Structural Break')
          k.short.names <- c('PreStructuralBreak', 'PostStructuralBreak')

         for(i in 1:2)
         {
             x.number.of.knots <- seq(from=0, to=jd$number.of.knots.ubound[i])
             r <- jd$mu.number.of.knots.prior[1,i]
             p <- 1 / (1 + jd$mu.number.of.knots.prior[2,i])
             N <- prod(dim(k[1,,]))

             plot.data[[i]] <- array(NA,
                                     c(length(x.number.of.knots), 3),
                                     list(NULL, c('NumberOfKnots', 'Prior', 'Posterior')))


             for(j in seq_along(x.number.of.knots))
             {
                 x.j <- x.number.of.knots[j]
                 plot.data[[i]][j, 'NumberOfKnots'] <- x.j
                 plot.data[[i]][j, 'Prior'] <- gamma(r + x.j) / factorial(x.j) / gamma(r) * (1 - p) ^ r * p ^ x.j
                 plot.data[[i]][j, 'Posterior'] <- sum( k == x.j) / N
             }

             plot.data[[i]][,'Posterior'] <-   plot.data[[i]][,'Posterior'] / sum(plot.data[[i]][,'Posterior'])
         }

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



          if(plot)
          {
              op <- par(no.readonly = TRUE)
              on.exit(par(op))

              par(mfrow=c(1,2), oma = c(3, 0, 0, 0))  #several plots in one row

              for(i in 1:2)
              {
                  plot(x=range(plot.data[[i]][,'NumberOfKnots']) + c(-0.75, 0.75),
                       y=range(plot.data[[i]][,c('Prior', 'Posterior')]),
                       type='n',
                       cex.axis=1.25,
                       cex.lab=1.25,
                       xlab='Number of Knots',
                       ylab='Relative Frequency',
                       main=k.nice.names[i],
                       cex.main=1.25,
                       font.main=1)

                  for(j in  1:dim(plot.data[[i]])[1])
                  {

                      bar(x=plot.data[[i]][j,'NumberOfKnots'],
                          y=plot.data[[i]][j,'Prior'],
                          col='grey',
                          offset='L')

                      bar(x=plot.data[[i]][j,'NumberOfKnots'],
                          y=plot.data[[i]][j,'Posterior'],
                          col='black',
                          offset='R')
                  }




                  if(i == 1)
                  {
                      legend(x=grconvertX( mean(par('omd')[1:2]), from='ndc', to='user'),
                             y=grconvertY( par('omd')[3], from='ndc', to='user'),
                             c("Prior","Posterior"),
                             col=c('grey','black'),
                             pch=15,
                             pt.cex=3,
                             horiz=TRUE,
                             bty='n',
                             xpd=NA,
                             xjust=0.5)
                  }

              }
          }



          ans <- list()

          for(i in 1:2)
          {

              ans[[k.short.names[i]]] <- plot.data[[i]][,'Posterior']
              names(ans[[k.short.names[i]]]) <- plot.data[[i]][,'NumberOfKnots']
          }
          return(invisible(ans))

      })

##' A method to generate the trace plots for select rate of decay values.
##'
##'
##' @name rateOfDecayTracePlot,BreakAnnualAggLossDevModelOutput-method
##' @param object The object of type \code{BreakAnnualAggLossDevModelOutput} from which to generate the trace plots.
##' @param elements A numeric vector indicating for which elements to plot the trace.  Valid values are 2 through the number of columns in the observed triangle.  If NULL, values are selected automatically.
##' @param preBreak A logical value indicating whether to plot the trace for the pre-break rate(s) of decay or the post-break rate(s) of decay.
##' @param \dots Additional arguments used by other methods.  Not utilized by this method.
##' @return NULL invisibly.  Only called for the side effect of plotting.
##' @docType methods
##' @seealso \code{\link{rateOfDecayTracePlot}}
##' @seealso \code{\link{rateOfDecay}}
setMethod('rateOfDecayTracePlot',
          signature(object='BreakAnnualAggLossDevModelOutput'),
          function(object, elements=NULL, preBreak=TRUE, ...)
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
          R.i <- ifelse(preBreak, 1, 2)
          R <- slot(object@R, 'value')
          delta.obs <- exp(R[R.i,-1,,] - R[R.i,-K,,]) - 1
          plot.trace.plots(delta.obs[elements-1,,], paste(ifelse(preBreak, 'Pre-Break', 'Post-Break'), ' Rate of Decay :', elements, sep=''))
      })


##' A method to generate the trace plots for select consumption path values.
##'
##'
##' @name consumptionPathTracePlot,BreakAnnualAggLossDevModelOutput-method
##' @param object The object of type \code{BreakAnnualAggLossDevModelOutput} from which to generate the trace plots.
##' @param elements A numeric vector indicating the elements for which to plot the trace.  Valid values are 1 through the number of development years (columns) in the observed triangle.  If NULL, values are selected automatically.
##' @param preBreak A logical value indicating whether to plot the trace for the pre-break consumption path or the post-break consumption path.
##' @param \dots Additional arguments used by other methods.  Not utilized by this method.
##' @return NULL invisibly.  Only called for the side effect of plotting.
##' @docType methods
##' @seealso \code{\link{consumptionPathTracePlot}}
##' @seealso \code{\link{consumptionPath}}
setMethod('consumptionPathTracePlot',
          signature(object='BreakAnnualAggLossDevModelOutput'),
          function(object, elements=NULL, preBreak=TRUE, ...)
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

          R.i <- ifelse(preBreak, 1, 2)
          plot.trace.plots(slot(object@R, 'value')[R.i,elements,,], paste(ifelse(preBreak, 'Pre-Break', 'Post-Break'), ' Consumption Path :', elements, sep=''))
      })


##' A method to plot autocorrelations found in the \acronym{MCMC} samples for select parameters.
##'
##' Chains with high autocorrelation require a longer burnin and more samples to fully explore the parameter space.
##'
##' @name mcmcACF,BreakAnnualAggLossDevModelOutput-method
##' @param object The object from which to plot autocorrelations.
##' @return Called for the side effect of plotting.
##' @docType methods
##' @seealso \code{\link{mcmcACF}}
##' @seealso \code{\link[=mcmcACF,StandardAnnualAggLossDevModelOutput-method]{mcmcACF("StandardAnnualAggLossDevModelOutput")}}
##' @TODO Add option to plot other values.  Currently only plots \dQuote{First Rate Of Decay} (for the pre and post break), \dQuote{First Calendar Year Effect Error,} and \dQuote{First Exposure Year Growth,}
setMethod('mcmcACF',
          signature(object='BreakAnnualAggLossDevModelOutput'),
          function(object)
      {
          ##Autocorrelation

          n.chain <- dim(slot(object@df, 'value'))['chain']
          K <- getTriDim(object@input)[1]

          op <- par(no.readonly=TRUE)
          on.exit(par(op))

          par(mfcol=c(4,n.chain))

          R <- slot(object@R, 'value')
          delta.1 <- R[1,-1, , ] - R[1,-K, , ]
          delta.2 <- R[2,-1, , ] - R[2,-K, , ]



          for(i in 1:n.chain)
          {
              acf(delta.1[1,,i],
                  cex.axis=1.25,
                  cex.lab=1.25,
                  main=paste('First Pre-Break Rate Of Decay\nChain:', i))

              acf(delta.2[1,,i],
                  cex.axis=1.25,
                  cex.lab=1.25,
                  main=paste('First Post-Break Rate Of Decay\nChain:', i))


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

##' A generic function to plot and/or return the posterior change point.
##'
##' When incorporating a structural break, the user has the option of specifying either 1) the first year in which the new regime applies or 2) a (inclusive) range in which the first year in the new regime applies.
##' If the user specifies a range, the actual year is estimated as a model parameter.
##' This function allows for the retrieval/illustration of the posterior for this estimate.
##'
##' @name firstYearInNewRegime
##' @param object The object from which to plot and/or return the posterior change point estimate.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.
##' @docType genericFunction
##' @seealso \code{\link[=firstYearInNewRegime,BreakAnnualAggLossDevModelOutput-method]{firstYearInNewRegime("BreakAnnualAggLossDevModelOutput")}}
##' @seealso \code{\link{firstYearInNewRegimeTracePlot}}
##' @exportMethod  firstYearInNewRegime
setGenericVerif('firstYearInNewRegime',
                function(object, plot=TRUE)
                standardGeneric('firstYearInNewRegime'))

##' A method to plot and/or return the posterior change point.
##'
##' @name firstYearInNewRegime,BreakAnnualAggLossDevModelOutput-method
##' @param object The object from which to plot and/or return the posterior change point estimate.
##' @param plot A logical value. If \code{TRUE}, the plot is generated and the statistics are returned; otherwise only the statistics are returned.
##' @return Mainly called for the side effect of plotting.  Also returns a named numeric vector with names corresponding to the first year in the new regime and with values corresponding to the density.  Returned invisibly.
##' @docType methods
##' @seealso \code{\link{firstYearInNewRegime}}
##' @seealso \code{\link{firstYearInNewRegimeTracePlot}}
setMethod('firstYearInNewRegime',
          signature(object='BreakAnnualAggLossDevModelOutput'),
          function(object, plot)
      {


          exposure.years <- object@input@exposureYears

          breaks <- seq(from=exposure.years[1] - 0.5,
                        to=max(exposure.years) + 0.5,
                        by=1)
          first.year.in.post.vec <- as.vector(slot(object@first.row.in.post, 'value')) + exposure.years[1] - 1

          if(plot)
          {
              plot.f <- function()
              {

                  hist(first.year.in.post.vec,
                       breaks=breaks,
                       include.lowest=TRUE,
                       plot=TRUE,
                       freq=FALSE,
                       all(diff(breaks)==1),
                       col=gray(0.8),
                       ylim=c(0,1),
                       xlim=range(breaks),
                       cex.axis=1.25,
                       cex.lab=1.25,
                       cex.main=1.25,
                       font.main=1,
                       xlab='First Year in Post-Structural Break Regime',
                       ylab='Relative Frequency',
                       main=NULL)


                  jd <- getJagsData(object@input)

                  first.year.in.new.regime.range <- object@input@rangeForFirstYearInNewRegime
                  first.year.in.new.regime.seq <- seq(from=first.year.in.new.regime.range[1], to=first.year.in.new.regime.range[2], by=1)
                  first.year.in.new.regime.breaks <- seq(from=first.year.in.new.regime.range[1] - 0.5,
                                                         to=max(first.year.in.new.regime.range) + 0.5,
                                                         by=1)


                  break.row.priors <- jd$break.row.priors

                  len <- length(first.year.in.new.regime.seq)
                  prior.prob <- numeric(len)

                  zero.one.endpoints <- seq(from=0, to=1, length.out=len + 1)

                  for(i in 1:len)
                  {

                      prior.prob[i] <- pbeta(zero.one.endpoints[i+1], break.row.priors[1], break.row.priors[2]) - pbeta(zero.one.endpoints[i], break.row.priors[1], break.row.priors[2])
                      lines(x=c(first.year.in.new.regime.breaks[i], first.year.in.new.regime.breaks[i], first.year.in.new.regime.breaks[i+1], first.year.in.new.regime.breaks[i+1], first.year.in.new.regime.breaks[i]),
                            y=c(0, prior.prob[i], prior.prob[i], 0, 0))
                  }
              }

              legend.f <- function()
              {
                  legend('center',
                         c("Posterior (gray) and Prior"),
                         bty='n',
                         xpd=NA,
                         horiz=TRUE)
              }

              plot.top.bottom(plot.f, legend.f)
          }

          h <- hist(first.year.in.post.vec,
                    breaks=breaks,
                    include.lowest=TRUE,
                    plot=FALSE)


              ans <- h$density
              names(ans) <- h$mids

          return(invisible(ans))


      })

##' A generic function to generate the trace plot for the posterior change point.
##'
##'
##' @name firstYearInNewRegimeTracePlot
##' @param object The object from which to generate the trace plot for the change point estimate.
##' @return Only called for the side effect of plotting.
##' @docType genericFunction
##' @seealso \code{\link[=firstYearInNewRegimeTracePlot,BreakAnnualAggLossDevModelOutput-method]{firstYearInNewRegimeTracePlot("BreakAnnualAggLossDevModelOutput")}}
##' @seealso \code{\link{firstYearInNewRegimeTracePlot}}
##' @exportMethod  firstYearInNewRegimeTracePlot
setGenericVerif('firstYearInNewRegimeTracePlot',
                function(object)
                standardGeneric('firstYearInNewRegimeTracePlot'))

##' A method to generate the trace plot for the posterior change point.
##'
##' @name firstYearInNewRegimeTracePlot,BreakAnnualAggLossDevModelOutput-method
##' @param object The object from which to generate the trace plot for the change point estimate.
##' @return Only called for the side effect of plotting.
##' @docType methods
##' @seealso \code{\link{firstYearInNewRegimeTracePlot}}
##' @seealso \code{\link{firstYearInNewRegime}}
setMethod('firstYearInNewRegimeTracePlot',
          signature(object='BreakAnnualAggLossDevModelOutput'),
          function(object)
      {
          exposure.years <- object@input@exposureYears
          first.year.in.post <- slot(object@first.row.in.post, 'value') + exposure.years[1] - 1
          plot.trace.plots(first.year.in.post, 'First Year in Post-Structural Break Regime')

      })

