
## For NB2 model fitting on the original & simulated datasets
# dispersion estimation: tagwise dispersion in edgeR (default with shrinkage)

# update: 2012/12/04

model_edgeR_tagwise <- function(counts, x, lib.sizes = colSums(counts)){
  
#   nr = dim(counts)[1]
#   nc = dim(counts)[2]
#   n = nr * nc
  
  # initializations
#   mu.hat.m = matrix(0, nr = nr, nc = nc)
#   phi.hat.m = matrix(0, nr = nr, nc = nc)
#   res.m = matrix(0, nr = nr, nc = nc)
  
  ## edgeR tagwise dispersion:
  # convert model matrix into group index
  grp.ids = factor(apply(x, 1, function(x){paste(rev(x), collapse = ".")}), 
                   labels = seq(ncol(x)))
  d = DGEList(counts=counts, lib.size=lib.sizes, group = grp.ids)
  design = model.matrix(~grp.ids, data=d$samples)
  e.com = estimateGLMCommonDisp(d, design, verbose=TRUE)
  e.tag = estimateGLMTagwiseDisp(e.com, design)
  # e.tag: trend is FALSE since we didn't use estimateGLMTrendedDisp() beforehand
  tag.fit = glmFit(d, design, dispersion=e.tag$tagwise.dispersion)
  
  # extract quantities:
  mu.hat.m = tag.fit$fitted.values   # mu may be close to 0
  phi.hat.m = tag.fit$dispersion     # there may be NA's
  v = mu.hat.m + phi.hat.m * mu.hat.m^2 + 1e-14
  res.m = (counts - mu.hat.m) / sqrt(v) 
  res.om = t(apply(res.m, 1, sort))  # order each row first! (a matrix still)
  ord.res.v = as.vector(t(res.om))
  
  # save as a list
  model_tag_m_obj = list(mu.hat.mat = mu.hat.m,
                         res.mat = res.m,
                         res.omat = res.om,
                         ord.res.vec = ord.res.v,
                         phi.hat.mat = phi.hat.m
                         #fit.tag = tag.fit
  )
  return(model_tag_m_obj)
}