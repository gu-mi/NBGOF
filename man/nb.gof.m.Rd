\name{nb.gof.m}
\alias{nb.gof.m}
\title{Main Function of Implementing Simulation-based Goodness-of-Fit Tests on an
RNA-Seq Dataset}
\usage{
  nb.gof.m(counts, x, lib.sizes=colSums(counts), sim=999,
    model="NBP", method="MAPL", min.n=100, prior.df = 10,
    design = "single", ncores = NULL)
}
\arguments{
  \item{counts}{an m-by-n count matrix of non-negative
  integers. For a typical RNA-Seq experiment, this is the
  read counts with m genes and n samples.}

  \item{x}{an n-by-p design matrix.}

  \item{lib.sizes}{library sizes of an RNA-Seq experiment.
  Default is the column sums of the \code{counts} matrix.}

  \item{sim}{number of simulations performed.}

  \item{model}{a string of characters specifying the
  negative binomial dispersion model used to fit the data.
  Currently the following dispersion models are available
  to be checked for goodness-of-fit: \itemize{ \item NBP
  dispersion model in the \code{NBPSeq} package
  (\code{NBP}) \item NBQ dispersion model in the
  \code{NBPSeq} pacakge (\code{NBQ}) \item NB common
  dispersion model in the \code{\link{edgeR}} package
  (\code{edgeR-common}) \item NB genewise dispersion model
  in the \code{\link{edgeR}} package
  (\code{edgeR-genewise}) \item NB trended (non-parametric)
  dispersion model in the \code{\link{edgeR}} package
  (\code{edgeR-trended}) \item NB tagwise-common dispersion
  model in the \code{\link{edgeR}} package
  (\code{edgeR-tagcom}) \item NB tagwise-trended dispersion
  model in the \code{\link{edgeR}} package
  (\code{edgeR-tagtrd}) } Users should specify
  \strong{exactly} the same characters as the ones in
  paratheses above for each dispersion model.}

  \item{method}{method for estimating dispersions. MAPL:
  maximum adjusted profile likelihood}

  \item{min.n}{for \code{edgeR-trended} model only: specify
  the minimim number of genes in a bin (default is 100).}

  \item{prior.df}{control of the shrinkage applied for
  \code{edgeR} genewise, tagwise (and its variants)
  dispersion models. Setting \code{prior.df=0} gives the
  genewise model without any shrinkage. Setting
  \code{prior.df=10} (default in \code{edgeR}) gives the
  tagwise model with shrinkage, either towards the common
  dispersion (\code{edgeR-tagcom}) or the global dispersion
  trend (\code{edgeR-tagtrd}).}

  \item{design}{specifications of testing (1) a
  single-group model (\code{single}); (2) a multiple-group
  model (\code{multiple}); (3) complex design with
  interactions (\code{complex}). \strong{The codes are only
  tested for single- and multiple-group cases, and the
  modeling of dispersions in the multiple-group case is
  still under consideration! For the complex design with
  interactions, though we can pass the design matrix
  directly in edgeR, the results may not make any sense!}}

  \item{ncores}{number of CPU cores to use. If unspecified,
  use the total number of CPUs minus 1.}
}
\value{
  An object of class "gofm" to which other methods can be
  applied.
}
\description{
  This function performs simulation-based goodness-of-fit
  tests of different negative binomial dispersion models
  for an RNA-Seq dataset with a known design matrix.
  Currently supported models are implemented in the
  \code{NBPSeq} and \code{edgeR} packages.
}
\details{
  When the response is a count matrix, we can use this
  function to test the goodness-of-fit of a specified
  negative binomial dispersion model. Specifically, this
  function calls \code{\link{model.nbp.m}} to fit the NBP
  dispersion model, calls \code{\link{model.nbq.m}} to fit
  the NBQ dispersion model, calls
  \code{\link{model.edgeR.common}} to fit the common
  dispersion model, calls
  \code{\link{model.edgeR.genewise}} to fit the pure
  genewise dispersion model (without shrinkage), calls
  \code{\link{model.edgeR.trended}} to fit the
  non-parametric trended dispersion model in \code{edgeR},
  calls \code{\link{model.edgeR.tagcom}} to fit the tagwise
  model with the dispersions shrinking towards a common
  dispersion, and calls \code{\link{model.edgeR.tagtrd}} to
  fit the tagwise model with the dispersions shrinking
  towards the global dispersion trend.
}
\examples{
## load package into session:
library(NBGOF)

## basic set-up of the model:
seed = 315825
sim = 999
conf.env = 0.95

## basic parameter specifications:
m = 1000; # 1000 genes
n = 6;    # 6 samples
s = 1e6;  # lib.sizes
offset = log(s);  # make sure offset for NBP functions should take the log
offset    # 13.82

## simulate coefficients:
beta = matrix(0, m, 1);  # beta0 only
set.seed(seed)
beta[,1] = rnorm(m, 5, 2) - offset;   # beta0 ~ normal (mean and a based on real RNA-Seq data)
beta  # m-by-1 matrix

## design matrix (single group only):
x = matrix(rep(1,n))
x

## specify mean levels:
mu = round(t(s * exp(x \%*\% t(beta))));
sum(rowSums(mu) == 0)   # 4
mu[rowSums(mu) == 0, ] = 1
sum(rowSums(mu) == 0)   # 0
pi = mu/s

## simulate an m-by-n count matrix mimicking a RNA-Seq dataset:
## simulate NB1.8 data with random uniform noise added to the dispersion
alpha1 = -0.2  # NB1.8
phi0 = 0.02
alpha0 = log(phi0)
phi.nbp = phi0 * pi^alpha1
range(phi.nbp)
cbind(mu[,1], phi.nbp[,1]);

## add noise:
a = 0.1
set.seed(seed)
phi.noi.vec = phi.nbp[ ,1] * exp(matrix(runif(m, -a, a), nr=m, nc=1))
phi.noi = matrix(phi.noi.vec, nr=m, nc=n)
range(phi.noi)
cbind(mu[,1], phi.noi[,1])  # make sure phi's are in reasonable range

## generate NBP response with added noise:
set.seed(seed)
y = rnbinom(m * n, mu=mu, size=1/phi.noi)
dim(y) = dim(mu)
rownames(y) = paste("g", seq(1,m), sep="")
colnames(y) = paste("s", seq(1,n), sep="")
plot(mu, phi.noi, log="xy")

## GOF tests for different dispersion models, using parallel computing:
## CAUTION: may be time-consuming depending on the size of data and simulations\\
nc = detectCores() - 1
fnbp.noip = nb.gof.m(counts=y, x=x, sim=sim, model="NBP", ncores=nc)
fnbq.noip = nb.gof.m(counts=y, x=x, sim=sim, model="NBQ", ncores=nc)
fcom.noip = nb.gof.m(counts=y, x=x, sim=sim, model="edgeR-common", ncores=nc)
fgen.noip = nb.gof.m(counts=y, x=x, sim=sim, model="edgeR-genewise", ncores=nc)
ftgc.noip = nb.gof.m(counts=y, x=x, sim=sim, model="edgeR-tagcom", ncores=nc)
ftgt.noip = nb.gof.m(counts=y, x=x, sim=sim, model="edgeR-tagtrd", ncores=nc)
ftrd.noip = nb.gof.m(counts=y, x=x, sim=sim, model="edgeR-trended", ncores=nc)

## summarize the GOF test results:
summary(fnbp.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
summary(fnbq.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
summary(fcom.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
summary(fgen.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
summary(ftgc.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
summary(ftgt.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")
summary(ftrd.noip, conf.env=conf.env, data.note="NB1.8-noise(a=0.1)")

## evaluate GOF p-values using histograms and quantile-quantile uniform plots:
# pdf(file="plots-NBP-01-1grp.pdf", width=8, height=24)
par(mfrow=c(6,2))
plot(fnbp.noip, type="hist")
plot(fnbp.noip, type="qq", pch=21, bg="blue", bty="n")
plot(fcom.noip, type="hist")
plot(fcom.noip, type="qq", pch=21, bg="blue", bty="n")
plot(fgen.noip, type="hist")
plot(fgen.noip, type="qq", pch=21, bg="blue", bty="n")
plot(ftgc.noip, type="hist")
plot(ftgc.noip, type="qq", pch=21, bg="blue", bty="n")
plot(ftgt.noip, type="hist")
plot(ftgt.noip, type="qq", pch=21, bg="blue", bty="n")
plot(ftrd.noip, type="hist")
plot(ftrd.noip, type="qq", pch=21, bg="blue", bty="n")
# dev.off()
}
\author{
  Gu Mi <mig@stat.oregonstate.edu>, Yanming Di, Daniel
  Schafer
}
\references{
  Mi G, Di Y and Schafer DW (2013). Goodness-of-Fit Tests
  and Model Diagnostics for Negative Binomial Regression of
  RNA sequencing Data. \emph{Biometrics} (submitted).

  Di Y, Schafer DW, Cumbie JS, and Chang JH (2011): "The
  NBP Negative Binomial Model for Assessing Differential
  Gene Expression from RNA-Seq", \emph{Statistical
  Applications in Genetics and Molecular Biology}, 10 (1).

  McCarthy DJ, Chen Y and Smyth GK (2012). Differential
  expression analysis of multifactor RNA-Seq experiments
  with respect to biological variation. \emph{Nucleic Acids
  Research} 40, 4288-4297.

  Cox, DR, and Reid, N (1987). Parameter orthogonality and
  approximate conditional inference. \emph{Journal of the
  Royal Statistical Society Series B} 49, 1-39.

  See \url{https://github.com/gu-mi/NBGOF/wiki/} for more
  details.
}

