irls.nb = function(y, s, x, phi, beta0, ..., print.level=0) {
  m = dim(y)[1];
  n = dim(y)[2];
  p = dim(x)[2];
  
  phi = matrix(phi, m, n, byrow=TRUE);
  
  if (print.level > 0)
    print("Estimating NB regression coefficients using IRLS.");
  
  res=list(mu=matrix(NA, m, n), beta=matrix(NA, m, p), phi=matrix(NA, m, n),
           conv=logical(m), iter=numeric(m));
  
  if (print.level > 1) {
    ## Set up progress bar
    pb=txtProgressBar(style=3);
  }
  
  for (i in 1:m) {
    if (print.level > 1) {
      setTxtProgressBar(pb, i/m);
    }   
    res0 = irls.nb.1(y[i,], s, x, phi[i,], beta0,
                     ...,
                     print.level=print.level-1);
    res$mu[i,] =  res0$mu;
    res$beta[i,] = res0$beta;
    res$conv[i] = res0$conv;
    res$iter[i] = res0$iter;
  }
  
  if (print.level>1) close(pb);
  
  res;
}