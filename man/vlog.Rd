\name{vlog}
\alias{vlog}
\title{User-defined link function for gamma regression}
\usage{
  family = Gamma(link = vlog())
}
\description{
  User-defined link function for gamma regression. This is
  made available to be passed to the \code{family=}
  argument in the \code{glm} function
}
\examples{
library(NBGOF)
set.seed(123)
n = 1000
x = runif(n)
sh = 2
vv = vlog()
y = rgamma(n, scale=vv$linkinv(2+3*x)/sh, shape=sh)
fit1 = glm(y~x, family=Gamma(link=vv))
fit1
sum.fit1 = summary(fit1)
1/sum.fit1$dispersion   # compare with sh: close to each other as n increases
}
\author{
  Gu Mi <mig@stat.oregonstate.edu>, Yanming Di, Daniel
  Schafer
}
\keyword{internal}

