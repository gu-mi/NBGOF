

##' @title Implement genewise dispersion approach (report dispersion and fitted means)
##' 
##' @param y
##' @param s
##' @param x
##' @param beta0
##' @param kappa
##' 
##' @return phi.hat and mu.hat
##' 
##' @export
##' 
##' @author Yanming Di Gu Mi
##' 
fit.nb.reg.1 = function(y, s, x, beta0=rep(NA, dim(x)[2]), kappa=NA) {
  
  ## Find preliminary estiamtes of mu assuming phi=0.1. Will serve as
  ## initial values for the later irls algorithm.
  mustart = irls.nb.1(y, s, x, 0.1, beta0)$mu;
  
  if (is.na(kappa)) {
    
    ## Log likelihood of log(kappa)
    ll =  function(lkappa) {
      kappa = exp(lkappa);
      res = irls.nb.1(y, s, x, 1/kappa, beta0, mustart);
      sum(dnbinom(y, size =kappa, mu=res$mu, log=TRUE)); 
      ## log.likelihood.nb(kappa, res$mu, y);
    }
    
    res = optimize(ll, c(log(1e-20), log(1e20)), maximum=TRUE);
    kappa.hat = exp(res$maximum);
  } else {
    kappa.hat = kappa
  }
  
  phi.hat = 1/kappa.hat;
  res = irls.nb.1(y, s, x, phi.hat, beta0, mustart);
  beta.hat = res$beta;
  mu.hat = res$mu;
  
  result = list(phi.hat = phi.hat, mu.hat = mu.hat)
  
  return(result)

}