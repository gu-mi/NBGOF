\name{plot.gofm}
\alias{plot.gofm}
\title{Plot Methods for "gofm" Objects}
\usage{
  plot(x, conf.env=0.95, data.note=NULL, leg.cex=1, put.lab
    = TRUE, ...)
}
\arguments{
  \item{x}{an object of class "gofm" (from \link{nb_gof_m}
  output)}

  \item{conf.env}{confidence level for the envelope}

  \item{data.note}{a note on how data are simulated or the
  source of data}

  \item{leg.cex}{size of the legend fonts}

  \item{put.lab}{logical value: whether to put gene index
  labels on the plot}

  \item{...}{for future use}
}
\value{
  A diagnostic plot for evaluating model adequacy. Each
  gene index (row index) that falls outside a pre-specified
  confidence envelope is in red color and labels (if
  \code{put.lab=TRUE}).
}
\description{
  A plot method for empirical probability plot based on a
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

