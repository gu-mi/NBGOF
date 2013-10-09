\name{plot.gofm}
\alias{plot.gofm}
\title{Plot Methods for "gofm" Objects (Histograms and Q-Q Uniform Plots of P-values)}
\usage{
  plot(x, type="qq", logscale=TRUE, ...)
}
\arguments{
  \item{x}{an object of class "gofm" from the
  \code{\link{nb.gof.m}} output}

  \item{type}{either a quantile-quantile uniform plot
  (\code{qq}), or a histogram (\code{hist})}

  \item{logscale}{whether use logscale for the uniform Q-Q
  plot (default is TRUE) as recommended in the \code{gap}
  package}

  \item{...}{for future use}
}
\value{
  A quantile-quantile uniform plot or a histogram of the
  GOF p-values.
}
\description{
  A plot method for evaluating the p-values based on a
  "gofm" object obtained from the \code{\link{nb.gof.m}}
  outputs.
}
\details{
  This is a generic function for plotting goodness-of-fit
  test results. Users can specify either a
  quantile-quantile (\code{qq}) uniform plot, or a
  histogram (\code{hist}) of the GOF p-values.
}
\author{
  Gu Mi <mig@stat.oregonstate.edu>, Yanming Di, Daniel
  Schafer
}
\seealso{
  \code{\link{nb.gof.m}} for simulated data examples.
}

