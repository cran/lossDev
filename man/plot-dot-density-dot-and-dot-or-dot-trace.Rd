\name{plot.density.and.or.trace}
\alias{plot.density.and.or.trace}
\title{A rather generic function to plot diagnostics for a single node (a one-dimensional node or a single slot from a multi-dimensional node).}
\usage{plot.density.and.or.trace( coda,
plotDensity,
plotTrace,
d.prior,
nice.parameter.name,
zero.line=FALSE,
lower.bound=NA,
upper.bound=NA,
draw.prior=TRUE)}
\description{A rather generic function to plot diagnostics for a single node (a one-dimensional node or a single slot from a multi-dimensional node). Intended for internal use only.}
\value{Mainly called for the side effect of plotting. Also returns a vector of quantiles.}
\arguments{\item{coda}{The code for the node.  Rows are iterations. Columns are chains.}
\item{plotDensity}{A logical value. If \code{TRUE}, the density is plotted. If \code{plotTrace} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.}
\item{plotTrace}{A logical value. If \code{TRUE}, the trace is plotted. If \code{plotDensity} is also \code{TRUE}, then two plots are generated.  If they are both \code{FALSE}, only the statistics are returned.}
\item{d.prior}{A function that takes an array of values and returns the prior density evaluated at those values.}
\item{nice.parameter.name}{A character value to use as labels in plots.}
\item{zero.line}{A logical value. Should a verical zero line be drawn on the density plot?}
\item{lower.bound}{Can be missing, used by \code{density}: (\code{from}).}
\item{upper.bound}{Can be missing, used by \code{density}: (\code{to}).}
\item{draw.prior}{Should the prior be drawn?}}
