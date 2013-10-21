\name{plotdata}
\alias{plotdata}
\title{Prepare Quantities for Diagnostic Plots (Histogram and Uniform Q-Q Plots) Using ggplot2}
\usage{
  plotdata(gofm.obj, envelope = 0.95, binwidth = 0.1)
}
\arguments{
  \item{gofm.obj}{an object of class "gofm" from the
  \code{\link{nb.gof.m}} output.}

  \item{envelope}{the level of confidence band when drawing
  a uniform Q-Q plot.}

  \item{binwidth}{the bin width used when drawing a
  histogram.}
}
\value{
  a list of three components.
}
\description{
  This function prepares quantities for making diagnostic
  plots of negative binomial dispersion models. The result
  (a list) will be passed to the \code{\link{diagPlot}}
  function to make histogram of p-values and uniform Q-Q
  plots.
}
\keyword{internal}

