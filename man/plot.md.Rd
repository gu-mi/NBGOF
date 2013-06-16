\name{plot.md}
\alias{plot.md}
\title{The Mean-Dispersion Plot (on Log-log Scale) with Fitted Curves}
\usage{
  plot.md(y, x, model = NULL, scatter = FALSE, legend =
    FALSE, ...)
}
\arguments{
  \item{y}{an m-by-n count matrix of non-negative integers.
  For a typical RNA-Seq experiment, this is the read counts
  with m genes and n samples.}

  \item{x}{an n-by-p design matrix.}

  \item{model}{a string of characters specifying the
  negative binomial model used to fit the data. Currently
  the following dispersion models are implemented: common,
  NBP, trended (non-parametric), tagwise-common and
  tagwise-trended.}

  \item{scatter}{logical: plot the points (\code{TRUE} if
  first invoking the plot). If \code{TRUE}, then a scatter
  plot is drawn, and the fitted NBP curve is superimposed.}

  \item{legend}{logical: put the legend of dispersion
  models on plot (\code{TRUE} when plotting the last
  curve).}

  \item{...}{for future use.}
}
\value{
  A mean-dispersion plot (log-log scale) with fitted curves
  from the following NB dispersion models: common, NBP,
  trended (non-parametric), tagwise-common and
  tagwise-trended.
}
\description{
  This function provides a quick-and-dirty mean-dispersion
  plot (log-log scale) with relative mean frequencies on
  the x-axis and estimated NB dispersions on the y-axis.
  Several fitted curves from NB dispersion models are
  superimposed on the plot.
}
\details{
  \strong{This function was originally used solely for
  paper figure productions. In the final release of the R
  package, it may be made invisible to end-users. Also,
  current implementations are rather inflexible.}
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
  The Examples section of the \code{\link{arab}} dataset.
}

