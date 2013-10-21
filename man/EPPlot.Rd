\name{EPPlot}
\alias{EPPlot}
\title{Empirical Probability Plot for Univariate Negative Binomial Regression Model Diagnostics}
\usage{
  EPPlot(x, conf.env=0.95, data.note=NULL, leg.cex=1, ...)
}
\arguments{
  \item{x}{an object of class "gofv" from the
  \code{\link{nb.gof.v}} output.}

  \item{conf.env}{confidence level for the envelope.}

  \item{data.note}{a note on how data are simulated or the
  source of data.}

  \item{leg.cex}{size of the legend fonts.}

  \item{...}{for future use.}
}
\value{
  An empirical probability plot with points outside a
  pre-specified confidence envelope colored in red. GOF
  test p-values are also reported on the plot.
}
\description{
  Empirical probability plot based on a "gofv" object
  obtained from the \code{\link{nb.gof.v}} outputs.
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
  \code{\link{nb.gof.v}} for simulated data examples, and
  \code{\link{earthquake}} for a real data example.
}

