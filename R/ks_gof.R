
ks_gof = function(y0, x, model = "NB2", sim=1){
  
  if (model == "NB2"){
    
    n = length(y0)
    res.h.mat = matrix(0, nrow = sim, ncol = n)
    
    # get Pearson residuals on observed response y0:
    fit.nb2 = model_nb2_v(y=y0, x=x)
    res0 = fit.nb2$res.vec
    # ordered Pearson residuals:
    ord.res0 = sort(res0)
    
    # extract fitted quantities:
    mu.hat.v0 = fit.nb2$mu.hat.v
    phi0 = fit.nb2$phi
    
    pb = txtProgressBar(style=3)
    for (i in 1:sim){
      setTxtProgressBar(pb, i/sim)
      # re-simulate new responses y:
      y.h = rnbinom(n, mu = mu.hat.v0, size = 1/phi0)
      fit.nb2.h = model_nb2_v(y=y.h, x=x)
      res.h.mat[i, ] = fit.nb2.h$res.vec
    }
    close(pb)
    
    res.sort = t(apply(res.h.mat, 1, sort))
    ord.res.h = apply(res.sort, 2, median)
    #
    ks_test_res = ks.test(x=ord.res0, y=ord.res.h, alternative="two.sided")  
  }
  return(list(pval = ks_test_res$p.value,
              ord.res0 = ord.res0,
              ord.res.h = ord.res.h)
  )
  
}