##################################################################################################
##                                                                                              ##
##    lossDev is an R-package.                                                                  ##
##    It is a Bayesian time series model of loss development.                                   ##
##    Features include skewed Student-t distribution with time-varying scale parameters,        ##
##    an expert prior for the calendar year effect,                                             ##
##    and accommodation for structural breaks in the consumption path of services.              ##
##                                                                                              ##
##    Copyright © 2009, 2010, 2011 National Council On Compensation Insurance Inc.,             ##
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

##' A function to cumulate a triangle.
##'
##' @param triangle A matrix of incremental payments.
##' @return A matrix resulting from cumulating the input triangle.
##' @export
cumulate <- function(triangle)
{
    ans <- triangle[,1]
    for(i in 2:dim(triangle)[2])
        ans <- cbind(ans, apply(triangle[,1:i], 1, sum))

    dimnames(ans) <- dimnames(triangle)
    return(ans)
}

##' A function to decumulate a triangle.
##'
##' @param triangle A matrix of cumulative payments.
##' @return A matrix resulting from decumulating the input triangle.
##' @export
decumulate <- function(triangle)
{
    ans <- triangle[,1]
    for(i in 2:dim(triangle)[2])
        ans <- cbind(ans, triangle[,i] - triangle[,i-1])

    dimnames(ans) <- dimnames(triangle)
    return(ans)
}

##' A function to plot a top and bottom graph on the same chart. Intended for internal use only.
##'
##' Main use is to plot the legend of a graph in the \dQuote{bottom.}
##'
##' @param f.top The function to call for ploting the top graph.
##' @param f.bottom The function to call for ploting the bottom graph.
##' @param top.scale A number between zero and 1 indicating the bottom of the top graph.
##' @param bottom.scale A number between zero and 1 indicating the top of the bottom graph.
##' @return NULL invisibly.  This function is called for its side effects.
plot.top.bottom <- function(f.top, f.bottom, top.scale=.95, bottom.scale=.1)
{
    op <- par(no.readonly=TRUE)
    on.exit(par(op))

    plot.new()

    par(fig=c(0,1,1 - top.scale,1),
        new=TRUE)
    f.top()

    par(fig=c(0,1,0,bottom.scale),
        new=TRUE)
    f.bottom()
    return(invisible(NULL))
}

##' A function to plot a top, middle, and bottom graph on the same chart. Intended for internal use only.
##'
##' @param f.top The function to call for plotting the top graph.
##' @param f.middle The function to call for plotting the middle graph.
##' @param f.bottom The function to call for plotting the bottom graph.
##' @param top.scale A number between zero and 1 indicating the bottom of the top graph.
##' @param middle.top A number between zero and 1 indicating the top of the middle graph.
##' @param middle.scale A numer between zero and 1 indicating the size of the middle graph.
##' @param bottom.scale A number between zero and 1 indicating the top of the bottom graph.
##' @return NULL invisibly.  This function is called for its side effects.
plot.top.middle.bottom <- function(f.top, f.middle, f.bottom, top.scale=.525, middle.top=.525, middle.scale=.1, bottom.scale=.525)
{
    op <- par(no.readonly=TRUE)
    on.exit(par(op))

    plot.new()

    par(fig=c(0,1,1 - top.scale,1),
        new=TRUE)
    f.top()

    par(fig=c(0,1,middle.top - middle.scale, middle.top),
        new=TRUE)
    f.middle()

    par(fig=c(0,1,0,bottom.scale),
        new=TRUE)
    f.bottom()
    return(invisible(NULL))
}

##' A function to get color values.  Intended for internal use only.
##'
##' Currenly, the main feature is that it never returns yellow.
##' May, in the future, be expanded.
##' @param i An integer value of 1 or greater.
##' @return A color value suitable for ploting commands.
get.color <- function(i)
{
    color.vec <- c('black', 'red', 'green', 'blue', 'lightblue', 'purple', 'grey')
    i <- (i-1) %% length(color.vec) + 1
    return(color.vec[i])
}



##' A rather generic function to plot diagnostics for a single node (a one-dimensional node or a single slot from a multi-dimensional node). Intended for internal use only.
##'
##'
##' @param coda The code for the node.  Rows are iterations. Columns are chains.
##' @param plotDensity A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.
##' @param plotTrace A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.
##' @param d.prior A function that takes an array of values and returns the prior density evaluated at those values.
##' @param nice.parameter.name A character value to use as labels in plots.
##' @param zero.line A logical value. Should a verical zero line be drawn on the density plot?
##' @param lower.bound Can be missing, used by \code{density}: (\code{from}).
##' @param upper.bound Can be missing, used by \code{density}: (\code{to}).
##' @param draw.prior Should the prior be drawn?
##' @return Mainly called for the side effect of plotting. Also returns a vector of quantiles.
##' @usage
##'    plot.density.and.or.trace( coda,
##'                               plotDensity,
##'                               plotTrace,
##'                               d.prior,
##'                               nice.parameter.name,
##'                               zero.line=FALSE,
##'                               lower.bound=NA,
##'                               upper.bound=NA,
##'                               draw.prior=TRUE)

plot.density.and.or.trace <- function(coda,  plotDensity, plotTrace, d.prior, nice.parameter.name, zero.line=FALSE, lower.bound=NA, upper.bound=NA, draw.prior=TRUE)
{



    if(!is.logical(plotDensity) || length(plotDensity) != 1 || is.na(plotDensity))
        stop('"plotDensity" must either be "TRUE" or "FALSE"')

    if(!is.logical(plotTrace) || length(plotTrace) != 1 || is.na(plotTrace))
        stop('"plotTrace" must either be "TRUE" or "FALSE"')

    f.xx <- function(fit)
    {
        n <- 200
        l <- qlogspline(0.005, fit)
        u <- qlogspline(0.995, fit)
        xx <- seq(from=l, to=u, length.out=n)
    }


    plot.d <- function()
    {

        penalty <-  lossDevOptions()[['logsplinePenaltyFunction']](coda)
        if(is.na(lower.bound) && is.na(upper.bound))
        {
            ##post <- density(coda)
            fit <- logspline(coda, penalty=penalty)


        } else if(!is.na(lower.bound) && is.na(upper.bound))
        {
            ##post <- density(coda, from=lower.bound)
            fit <- logspline(coda, lbound=min(coda), penalty=penalty)


        } else if(!is.na(upper.bound) && is.na(lower.bound))
        {
            ##post <- density(coda, to=upper.bound)
            fit <- logspline(coda, ubound=max(coda), penalty=penalty)


        } else {
            ##post <- density(coda, from=lower.bound, to=upper.bound)
             fit <- logspline(coda, lbound=min(coda), ubound=max(coda), penalty=penalty)

        }

        xx <- f.xx(fit)
        yy <- dlogspline(xx, fit)
        post <- list(x=xx, y=yy)


        if(draw.prior)
            prior <- list(x=post$x,
                          y= d.prior(post$x))
        else
            prior <- list(x=post$x,
                          y= rep(NA, length(post$x)))

        plot(
             x=range(post$x, prior$x),
             y=range(post$y, prior$y, na.rm=TRUE),
             type='n',
             xlab=nice.parameter.name,
             ylab="Density",
             cex.axis=1.25,
             cex.lab=1.25)

        if(zero.line)
            abline(v=0)

        lines(
              x=prior$x,
              y=prior$y,
              col='grey',
              lwd=3)

        lines(
              x=post$x,
              y=post$y,
              col='black',
              lwd=2)

    }

    d.legend <- function()
    {
        if(draw.prior)
            legend('center',
                   c('Prior', 'Posterior'),
                   lwd=c(3,2),
                   col=c('grey', 'black'),
                   horiz=TRUE,
                   bty='n',
                   xpd=NA)
        else
            legend('center',
                   c('Posterior'),
                   lwd=c(2),
                   col=c('black'),
                   horiz=TRUE,
                   bty='n',
                   xpd=NA)

    }

    plot.t <- function()
    {
        matplot(
                coda,
                cex.axis=1.25,
                cex.lab=1.25,
                xlab="Sample",
                ylab=nice.parameter.name,
                type='l',
                lty=1)
    }

    if(plotDensity && plotTrace)
        plot.top.middle.bottom(plot.d, d.legend, plot.t)
    else if( plotDensity)
        plot.top.bottom(plot.d, d.legend)
    else if(plotTrace)
        plot.t()


    return(invisible(quantile(coda)))
}

##' A rather generic function to plot (multiple) trace plots in one call on one graph. Intended for internal use only.
##'
##' Plots a trace plot for each of the first dimensions in coda.
##'
##' @param coda The coda for the node(s):  first dimension indicates the node;  second is iterations; third is chains.
##' @param names A character vector equal in length to the first dim of coda representing the names of the nodes (these are used to label the trace plots).
##' @return \code{NULL} invisibly.  Only called for the side effect of plotting.
plot.trace.plots <- function(coda, names)
{

    if(length(dim(coda)) == 2)
    {
        is.single <- TRUE
        elements <- 1
    }
    else
    {
        is.single <- FALSE
        elements <- 1:dim(coda)[1]
    }

    op <- par(no.readonly=TRUE)
    on.exit(par(op))

    par(mfrow=c(length(elements), 1))


    if(is.single)
    {

            matplot(coda,
                    xlab='Sample',
                    ylab='Parameter Value',
                    type='l',
                    cex.axis=1.25,
                    cex.lab=1.25,
                    font.main=1,
                    main=names[1])
    } else {
        for(i in elements)
        {
            matplot(coda[i,,],
                    xlab='Sample',
                    ylab='Parameter Value',
                    type='l',
                    cex.axis=1.25,
                    cex.lab=1.25,
                    font.main=1,
                    main=names[i])
        }
    }

    return(invisible(NULL))
}

