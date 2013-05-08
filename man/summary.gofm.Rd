\name{summary.gofm}
\alias{summary.gofm}
\title{Report summaries of the GOF test results for multivariate response}
\usage{
  summary(x, conf.env=0.95, data.note=NULL, ...)
}
\arguments{
  \item{x}{an object of class "gofm" (from \link{nb_gof_m}
  output)}

  \item{conf.env}{confidence level for the envelope}

  \item{data.note}{a note on how data are simulated or the
  source of data}

  \item{bin}{the number of bins used for calculating a
  standard chi-square p-value}

  \item{...}{for future use}
}
\value{
  Information of the dataset, simulation parameter
  specifications, exact binomial test results, and Monte
  Carlo GOF p-values
}
\description{
  A summary method of GOF test results based on a "gofm"
  object obtained from the \code{\link{nb_gof_m}} outputs
}
\details{
  This is a generic function for summarizing
  goodness-of-fit test results.
}
\author{
  Gu Mi, Yanming Di, Daniel Schafer
}
\seealso{
  See \code{\link{nb_gof_m}} for simulated data examples
}

