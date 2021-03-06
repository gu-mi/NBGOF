% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.gofm.R
\name{summary.gofm}
\alias{summary.gofm}
\alias{summary}
\title{Report Summaries of the GOF Test Results for NB Dispersion Models}
\usage{
summary(x, conf.env=0.95, data.note=NULL, ...)
}
\arguments{
\item{x}{an object of class "gofm" from the \code{\link{nb.gof.m}} output}

\item{conf.env}{confidence level for the envelope}

\item{data.note}{a note on how data are simulated or the source of data}

\item{...}{(for future use)}
}
\value{
Information of the dataset, simulation parameter specifications, exact binomial test
results, and Monte Carlo GOF p-values.
}
\description{
A summary method of GOF test results based on a "gofm" 
object obtained from the \code{\link{nb.gof.m}} outputs
}
\details{
This is a generic function for summarizing goodness-of-fit test results of testing
negative binomial dispersion models.
}
\references{
See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
}
\seealso{
\code{\link{nb.gof.m}} for simulated data examples, and \code{\link{arab}} for 
a real RNA-Seq data example.
}
\author{
Gu Mi <neo.migu@gmail.com>, Yanming Di, Daniel Schafer
}
