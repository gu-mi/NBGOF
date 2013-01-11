
## For NB2 model fitting on the original & simulated datasets
## dispersion estimation: common dispersion in edgeR

# update: 2013/01/08

model_edgeR_common <- function(counts, x, lib.sizes = colSums(counts)){
  
  ## edgeR common dispersion:
  
  # convert model matrix into group index
  grp.ids = factor(apply(x, 1, function(x){paste(rev(x), collapse = ".")}), 
                   labels = seq(ncol(x)))
  d = DGEList(counts=counts, lib.size=lib.sizes, group = grp.ids)
  design = model.matrix(~grp.ids, data=d$samples)
  e.com = estimateGLMCommonDisp(d, design, verbose=TRUE)
  com.fit = glmFit(d, design, dispersion=e.com$common.dispersion)
  
  # extract quantities:
  mu.hat.m = com.fit$fitted.values  # mu may be close to 0
  phi.hat.m = com.fit$dispersion    # there may be NA's
  v = mu.hat.m + phi.hat.m * mu.hat.m^2 + 1e-14  # add epsilon to avoid NaN in v
  res.m = (counts - mu.hat.m) / sqrt(v) 
  
  res.om = t(apply(res.m, 1, sort))  # order each row first! (a matrix still)
  ord.res.v = as.vector(t(res.om))
  
  # save as a list
  model_com_m_obj = list(mu.hat.mat = mu.hat.m,
                         res.mat = res.m,
                         res.omat = res.om,
                         ord.res.vec = ord.res.v,
                         phi.hat.mat = phi.hat.m,
                         fit.com = com.fit
  )
  return(model_com_m_obj)
}