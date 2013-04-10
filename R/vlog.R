

#' @title User-defined link function for gamma regression
#' 
#' @description User-defined link function for gamma regression. This is made available to be passed
#' to the \code{family=} argument in the \code{glm} function
#' 
#' @usage
#' family = Gamma(link = vlog())
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' @export
#' 
#' @examples
#' 
#' library(NBGOF)
#' set.seed(123)
#' n = 1000                     
#' x = runif(n)
#' sh = 2   
#' vv = vlog()                     
#' y = rgamma(n, scale=vv$linkinv(2+3*x)/sh, shape=sh)
#' fit1 = glm(y~x, family=Gamma(link=vv))                       
#' fit1
#' sum.fit1 = summary(fit1)
#' 1/sum.fit1$dispersion   # compare with sh: close to each other as n increases
#' 
vlog <- function() {
  ## link
  linkfun <- function(y) log(exp(y)-1)
  ## inverse link
  linkinv <- function(eta)  log(exp(eta)+1)
  ## derivative of invlink wrt eta
  mu.eta <- function(eta) { 1/(exp(-eta) + 1) }
  valideta <- function(eta) TRUE
  link <- "log(exp(y)-1)"
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta, 
                 name = link),
            class = "link-glm")
}