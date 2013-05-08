\name{nb_gof_v}
\alias{nb_gof_v}
\title{Main Function of Implementing Simulation-based Goodness-of-Fit Tests on
the regression of a univariate count response on one or more explanatory variables}
\usage{
  nb_gof_v(y, x, lib.sizes=NULL, sim=999, model.fit =
    "NB2", est.method="MLE")
}
\arguments{
  \item{y}{an n-by-1 vector of non-negative integers. For a
  typical RNA-Seq experiment, this may represent the read
  counts for a single gene}

  \item{x}{an n-by-p design matrix.}

  \item{lib.sizes}{library sizes of a RNA-Seq experiment.
  Default is 1 for all samples}

  \item{sim}{number of simulations performed.}

  \item{model.fit}{a string of characters specifying the
  negative binomial model used to fit the data. Currently
  the following dispersion models are available to be
  checked for goodness-of-fit: \itemize{ \item NB2:
  conventional NB2 regression model (\code{NB2}) \item NBP:
  NBP regression model (see \code{\link{NBPSeq}}) package
  (\code{NBP}) } Users are recommended to specify
  \strong{exactly} the same characters as the ones in
  paratheses above for each NB regression model.}

  \item{est.method}{either "MLE" or "APLE" for the
  estimation methods used}
}
\value{
  An object of class "gofv" to which other methods (plot,
  summary, etc.) can be applied.
}
\description{
  This function is designed to test goodness-of-fit of an
  NB2 or NBP negative binomial regression model with a
  known design matrix. Estimation methods for NBP model
  fitting include MLE and APLE
}
\details{
  When the response is a vector, we can use this function
  to test the goodness-of-fit of a specified negative
  binomial regression model.
}
\examples{
## load package into session
library(NBGOF)

## basic set-up of the model
seed = 31513
n = 100
beta.v = c(1, -3)

## specify a design matrix
X = cbind(rep(1,n),seq(1.5, 3.5, length=n))
s = rep(1e6, n)
mu = s * exp(X \%*\% beta.v)
pi = mu/s  # relative frequency
sim = 499

## simulate an n-by-1 vector for the read counts of a single gene
set.seed(seed)
alpha1 = -1  # NB1 data
phi0 = 0.0001
alpha0 = log(phi0);
phi.nb1 = phi0 * pi^alpha1
y.nb1 = rnbinom(n, size=1/phi.nb1, mu=mu)

## implement the GOF test for testing an NB regression model adequacy of NB2 and NBP models
# pdf("gofv-result.pdf", width=14, height=7)
par(mfrow=c(1,2))
# NB2 model fit
gf.nb1.nb2 = nb_gof_v(y.nb1, X, s, sim=sim, model.fit="NB2")
plot(gf.nb1.nb2, conf.env=0.95, data.note = "NB1", pch=".", cex=5)
# NBP model fit
gf.nb1.nbp = nb_gof_v(y.nb1, X, s, sim=sim, model.fit="NBP")
plot(gf.nb1.nbp, conf.env=0.95, data.note = "NB1", pch=".", cex=5)
# dev.off()
}
\author{
  Gu Mi <mig@stat.oregonstate.edu>, Yanming Di, Daniel
  Schafer
}
\references{
  \url{https://github.com/gu-mi/NBGOF/wiki/}
}

