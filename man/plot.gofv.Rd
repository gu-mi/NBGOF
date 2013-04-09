\name{plot.gofv}
\alias{plot.gofv}
\title{Plot Methods for "gofv" Objects}
\usage{
  plot(x, conf.env=0.95, data.note=NULL, leg.cex=1, ...)
}
\arguments{
  \item{x}{an object of class "gofv" (from \link{nb_gof_v}
  output)}

  \item{conf.env}{confidence level for the envelope}

  \item{data.note}{a note on how data are simulated or the
  source of data}

  \item{leg.cex}{size of the legend fonts}

  \item{...}{for future use}
}
\value{
  An empirical probability plot with points outside a
  pre-specified confidence envelope colored in red.
}
\description{
  A plot method for empirical probability plot based on a
  "gofv" object obtained from the \code{\link{nb_gof_v}}
  outputs
}
\details{
  This is a generic function for plotting goodness-of-fit
  test results.
}
\author{
  Gu Mi, Yanming Di, Daniel Schafer
}
\seealso{
  See \code{\link{nb_gof_v}} for simulated data examples
}

