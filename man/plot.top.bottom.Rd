\name{plot.top.bottom}
\alias{plot.top.bottom}
\title{A function to plot a top and bottom graph on the same chart.}
\usage{plot.top.bottom(f.top, f.bottom, top.scale=0.95, bottom.scale=0.1)}
\description{A function to plot a top and bottom graph on the same chart. Intended for internal use only.}
\details{Main use is to plot the legend of a graph in the \dQuote{bottom.}}
\value{NULL invisibly.  This function is called for its side effects.}
\arguments{\item{f.top}{The function to call for ploting the top graph.}
\item{f.bottom}{The function to call for ploting the bottom graph.}
\item{top.scale}{A number between zero and 1 indicating the bottom of the top graph.}
\item{bottom.scale}{A number between zero and 1 indicating the top of the bottom graph.}}
