
## code from Yanming on 2012/11/08 (concise version of estimate.disp.1.R)

## Estimate dispersion for each gene separately

library(NBPSeq);
library(numDeriv);

source("irls.nb.R")
source("irls.nb.1.R")

log.likelihood.nb= function(kappa, mu, y) {
  ## p = mu/(mu+kappa);
  ## l0= sum(log(gamma(kappa+y)) - log(gamma(kappa)) - log(gamma(1+y))
  ##       + y * log(mu) + kappa * log(kappa) - (y+kappa) * log(mu+kappa));
  
  sum(dnbinom(y, kappa, mu=mu, log=TRUE));
}

pl.phi.1 = function(phi, y, lib.sizes, xx,
                    beta0=rep(NA, dim(xx)[2]),
                    log.dispersion=TRUE,
                    print.level=0) {
  
  if (log.dispersion) {
    phi = exp(phi);
  }
  fit = irls.nb.1(y, lib.sizes, xx, phi, beta0);
  log.likelihood.nb(1/phi, fit$mu, y);
  
}

mle.phi.1 = function(likelihood, interval, y, s, x, beta0, 
                     information=TRUE,
                     log.dispersion=TRUE,
                     print.level=1) {
  
  ## Maximize the profile likelihood of phi
  obj = optimize(likelihood, interval=interval,
                 y=y, lib.sizes=s, xx=x, beta0=beta0,
                 log.dispersion=log.dispersion,
                 print.level=print.level-1, 
                 maximum=TRUE);
  
  if (information) {
    h = -hessian(likelihood, obj$maximum,
                 y=y, lib.sizes=s, xx=x, beta0=beta0,
                 log.dispersion=log.dispersion); 
  } else {
    h = NA;
  }
  
  ## Return the MLE
  list(maximum=obj$maximum, information=h);
}

estimate.disp.1 = function(counts, lib.sizes, x, beta0, 
                           likelihood, interval = c(-20, 7),
                           log.dispersion=TRUE,
                           information=TRUE,
                           print.level=0) {
  n = dim(counts)[1];
  
  phi = numeric(n);
  i.hat = numeric(n);
  
  ## Set up progress bar
  #pb=txtProgressBar(style=3);
  
  ## Estimate the disperison parameter for each gene separately
  for (i in 1:n) {
    #setTxtProgressBar(pb, i/n);
    res =  mle.phi.1(likelihood, interval, counts[i,], lib.sizes, x, beta0,
                     log.dispersion=log.dispersion,
                     information=information);
    phi[i]= res$maximum;
    i.hat[i] = res$information;
  }
  #close(pb);
  
  list(estimate=phi, information=i.hat);
}


nb.regression.1 = function(y, s, x, beta, likelihood=pl.phi.1) {
  
  if (!is.matrix(y)) {
    dim(y) = c(1, length(y));
  }
  
  ## By default, we estimate the MLE of log(phi)
  log.phi.hat = estimate.disp.1(y, s, x, beta,
                                information=FALSE,
                                likelihood=likelihood)$estimate;
  phi.hat = exp(log.phi.hat);
  
  ## Fit NB regression
  fit = irls.nb(y, s, x, phi.hat, beta0=beta);
  
  fit$phi = matrix(phi.hat, dim(y)[1], dim(y)[2]); 
  
  ## If you need variance, uncomment the following line
  fit$v = fit$mu + fit$phi * fit$mu^2;
  
  fit
}