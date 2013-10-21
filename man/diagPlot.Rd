\name{diagPlot}
\alias{diagPlot}
\title{Diagnostic Plots for Graphical Checks of Negative Binomial Dispersion Models}
\usage{
  diagPlot(gofm.obj, type="qq", envelope=0.95,
    binwidth=0.1)
}
\arguments{
  \item{gofm.obj}{an object of class "gofm" from the
  \code{\link{nb.gof.m}} output.}

  \item{type}{which diagnostic plot to make. Either
  \code{qq} for a uniform Q-Q plot, or \code{hist} for a
  histogram of p-values.}

  \item{envelope}{the level of confidence band when drawing
  a uniform Q-Q plot.}

  \item{binwidth}{the bin width used when drawing a
  histogram.}
}
\value{
  A ggplot object, either a histogram (\code{type="hist"})
  or a Q-Q plot (\code{type="qq"}).
}
\description{
  This function makes diagnostic plots based on the
  goodness-of-fit test results (an "gofm" object) from
  \code{\link{nb.gof.m}}: a histogram of p-values or a
  uniform quantile-quantile (Q-Q) plot is available.
}
\author{
  Gu Mi <mig@stat.oregonstate.edu>, Yanming Di, Daniel
  Schafer
}
\references{
  See \url{https://github.com/gu-mi/NBGOF/wiki/} for more
  details.
}
\seealso{
  \code{\link{nb.gof.m}} for simulated data examples.
}

