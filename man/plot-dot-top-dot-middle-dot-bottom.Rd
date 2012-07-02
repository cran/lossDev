\name{plot.top.middle.bottom}
\alias{plot.top.middle.bottom}
\title{A function to plot a top, middle, and bottom graph on the same chart.}
\usage{plot.top.middle.bottom(f.top, f.middle, f.bottom, top.scale=0.525,
    middle.top=0.525, middle.scale=0.1, bottom.scale=0.525)}
\description{A function to plot a top, middle, and bottom graph on the same chart. Intended for internal use only.}
\value{NULL invisibly.  This function is called for its side effects.}
\arguments{\item{f.top}{The function to call for plotting the top graph.}
\item{f.middle}{The function to call for plotting the middle graph.}
\item{f.bottom}{The function to call for plotting the bottom graph.}
\item{top.scale}{A number between zero and 1 indicating the bottom of the top graph.}
\item{middle.top}{A number between zero and 1 indicating the top of the middle graph.}
\item{middle.scale}{A numer between zero and 1 indicating the size of the middle graph.}
\item{bottom.scale}{A number between zero and 1 indicating the top of the bottom graph.}}
