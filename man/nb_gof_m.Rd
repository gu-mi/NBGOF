\name{nb_gof_m}
\alias{nb_gof_m}
\title{Main Function of Implementing Simulation-based Goodness-of-Fit Tests on a
RNA-Seq Dataset}
\usage{
  nb_gof_m(counts, x, lib.sizes=colSums(counts), sim=199,
    model.fit="NB")
}
\arguments{
  \item{counts}{an m-by-n count matrix of non-negative
  integers. For a typical RNA-Seq experiment, this is the
  read counts with m genes and n samples.}

  \item{x}{an n-by-p design matrix.}

  \item{lib.sizes}{library sizes of a RNA-Seq experiment.
  Default is the column sums of the \code{counts} matrix.}

  \item{sim}{number of simulations performed.}

  \item{model.fit}{a string of characters specifying the
  negative binomial dispersion model used to fit the data.
  Currently the following dispersion models are available
  to be checked for goodness-of-fit: \itemize{ \item NB
  genewise dispersion model in the \code{\link{NBPSeq}}
  package (\code{NB}) \item NBP dispersion model in the
  \code{\link{NBPSeq}} package (\code{NBP}) \item NB common
  dispersion model in the \code{\link{edgeR}} package
  (\code{edgeR-common}) \item NB tagwise dispersion model
  in the \code{\link{edgeR}} package (\code{edgeR-tagwise})
  \item NB trended dispersion model in the
  \code{\link{edgeR}} package (\code{edgeR-trended}) }
  Users are recommended to specify \strong{exactly} the
  same characters as the ones in paratheses above for each
  dispersion model.}
}
\value{
  An object of class "gofm" to which other methods (plot,
  summary, etc.) can be applied.
}
\description{
  This function is designed to test goodness-of-fit of
  different negative binomial dispersion models for a
  RNA-Seq dataset with a known design matrix.
}
\details{
  When the response is a count matrix, we can use this
  function to test the goodness-of-fit of a specified
  negative binomial dispersion model.
}
\examples{
## load package into session
library(NBGOF)

## basic set-up of the model
seed = 31513
sim = 999
conf.env = 0.95
n = 100;
r = 14;
mu = matrix(seq(10, 1000, length=n), n, r);
lib.sizes = rep(1e6, r);
pi = mu/lib.sizes

## simulate an m-by-n count matrix mimicking a RNA-Seq dataset
set.seed(seed)
alpha1 = 0   # NB2 regression model
phi0 <- 0.8  # == phi
alpha0 = log(phi0);
phi.nb2c = phi0 * pi^alpha1;  # == phi0
cbind(mu[,1], phi.nb2c[,1])   # make sure phi's are in reasonable range
y = rnbinom(n * r, mu=mu, size=1/phi.nb2c)  # response matrix
dim(y) = dim(mu);
rownames(y) = paste("g", seq(1,n), sep="")
colnames(y) = paste("s", seq(1,r), sep="")

## specify a design matrix
grp.ids = as.factor(c(rep(1,8), rep(2,4), rep(3,2)))
x = model.matrix(~grp.ids)

## implement the GOF test for testing an NB dispersion model adequacy
## CAUTION: may be time-consuming depending on the size of data and simulations

## consider all dispersion estimation methods:
# pdf(file=file.path(path1,"gof-nb2comphi-95.pdf"), width=12, height=8)
par(mfrow=c(2,3))
#
fnb2.nb2c = nb_gof_m(counts=y,x=x,lib.sizes=lib.sizes, sim=sim,model.fit="NB")
plot(fnb2.nb2c, conf.env=conf.env, data.note="NB2", col="azure4", pch=".", cex=3)
#
fnbp.nb2c = nb_gof_m(counts=y,x=x,lib.sizes=lib.sizes, sim=sim,model.fit="NBP")
plot(fnbp.nb2c, conf.env=conf.env, data.note="NB2", col="azure4", pch=".", cex=3)
#
fcom.nb2c = nb_gof_m(counts=y,x=x,lib.sizes=lib.sizes, sim=sim,model.fit="edgeR-common")
plot(fcom.nb2c, conf.env=conf.env, data.note="NB2", col="azure4", pch=".", cex=3)
#
fgen.nb2c = nb_gof_m(counts=y,x=x,lib.sizes=lib.sizes, sim=sim,model.fit="edgeR-genewise")
plot(fgen.nb2c, conf.env=conf.env, data.note="NB2", col="azure4", pch=".", cex=3)
#
ftag.nb2c = nb_gof_m(counts=y,x=x,lib.sizes=lib.sizes, sim=sim,model.fit="edgeR-tagwise")
plot(ftag.nb2c, conf.env=conf.env, data.note="NB2", col="azure4", pch=".", cex=3)
#
ftrd.nb2c = nb_gof_m(counts=y,x=x,lib.sizes=lib.sizes, sim=sim,model.fit="edgeR-trended")
plot(ftrd.nb2c, conf.env=conf.env, data.note="NB2", col="azure4", pch=".", cex=3)
# dev.off()

## summarize the GOF test results:
summary(fnb2.nb2c, conf.env=0.95, data.note="NB2 Common Dispersion Data")
summary(fnbp.nb2c, conf.env=0.95, data.note="NB2 Common Dispersion Data")
summary(fcom.nb2c, conf.env=0.95, data.note="NB2 Common Dispersion Data")
summary(fgen.nb2c, conf.env=0.95, data.note="NB2 Common Dispersion Data")
summary(ftag.nb2c, conf.env=0.95, data.note="NB2 Common Dispersion Data")
summary(ftrd.nb2c, conf.env=0.95, data.note="NB2 Common Dispersion Data")
}
\author{
  Gu Mi <mig@stat.oregonstate.edu>, Yanming Di, Daniel
  Schafer
}
\references{
  \url{https://github.com/gu-mi/NBGOF/wiki/}
}

