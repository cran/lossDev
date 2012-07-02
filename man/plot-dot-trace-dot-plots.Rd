\name{plot.trace.plots}
\alias{plot.trace.plots}
\title{A rather generic function to plot (multiple) trace plots in one call on one graph.}
\usage{plot.trace.plots(coda, names)}
\description{A rather generic function to plot (multiple) trace plots in one call on one graph. Intended for internal use only.}
\details{Plots a trace plot for each of the first dimensions in coda.}
\value{\code{NULL} invisibly.  Only called for the side effect of plotting.}
\arguments{\item{coda}{The coda for the node(s):  first dimension indicates the node;  second is iterations; third is chains.}
\item{names}{A character vector equal in length to the first dim of coda representing the names of the nodes (these are used to label the trace plots).}}
