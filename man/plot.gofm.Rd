\name{plot.gofm}
\alias{plot.gofm}
\title{Plot Methods for "gofm" Objects}
\usage{
  plot(x, model=NULL, type="qq", logscale=TRUE, ...)
}
\arguments{
  \item{x}{an object of class "gofm" (from \link{nb_gof_m}
  output)}

  \item{model}{the model used to fit the data}

  \item{type}{either a quantile-quantile uniform plot, or a
  histogram}

  \item{logscale}{whether use logscale for the uniform Q-Q
  plot (default is TRUE)}

  \item{...}{for future use}
}
\value{
  A quantile-quantile uniform plot or a histogram
}
\description{
  A plot method for evaluating the p-values based on a
  "gofm" object obtained from the \code{\link{nb_gof_m}}
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
  See \code{\link{nb_gof_m}} for simulated data examples
}

