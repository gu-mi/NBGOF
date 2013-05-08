\name{plot.md}
\alias{plot.md}
\title{mean-dispersion plot (log-log) scale for paper use}
\usage{
  plot.md(y, x, model = NULL, scatter = FALSE, legend =
    FALSE, data.type = NULL, ...)
}
\arguments{
  \item{y}{}

  \item{x}{an n-by-p design matrix.}

  \item{model}{a string of characters specifying the
  negative binomial model used to fit the data. Currently
  the following dispersion models are available to be
  checked for goodness-of-fit}

  \item{scatter}{}

  \item{legend}{}

  \item{data.type}{}
}
\value{
  An object
}
\description{
  This function is designed to
}
\details{
  When the response is a vector, we can use this function
  to test the goodness-of-fit of a specified negative
  binomial regression model.
}
\author{
  Gu Mi <mig@stat.oregonstate.edu>, Yanming Di, Daniel
  Schafer
}
\references{
  \url{https://github.com/gu-mi/NBGOF/wiki/}
}

