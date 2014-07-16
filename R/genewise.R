

##' @title Implement genewise dispersion approach (but without any tests)
##' 
##' @param nb.data
##' @param grp.ids
##' @param grp1
##' @param grp2
##' @param R
##' 
##' @return a data frame where each row is from hoa.1u() function output
##' 
##' @export
##' 
##' @author Yanming Di Gu Mi
##' 
genewise = function(nb.data, grp.ids, grp1, grp2, R = 100) {
  
  counts = nb.data$counts;
  eff.lib.sizes = nb.data$eff.lib.sizes;
  
  ## Extract relevant columns from the counts matrix
  id1 = grp.ids %in% grp1;
  id2 = grp.ids %in% grp2;
  
  y = cbind(counts[, id1, drop=FALSE], counts[, id2, drop=FALSE]);
  
  s = c(eff.lib.sizes[id1], eff.lib.sizes[id2]);
  
  ## Construct the model matrix
  treatment = unlist(list(grp.ids[id1], grp.ids[id2]));
  x = model.matrix(~factor(treatment, levels = c(grp1, grp2)));
  colnames(x) = c(grp1, grp2);
  
  ## The null hypothesis
  beta0 = c(NA, 0);
  n.pars = length(beta0) + 1;
  
  ## Test DE using HOA
  ## set.seed(999);
  
  m = nrow(y);
  n = ncol(y);
  ## m = 100;
  
#   res = data.frame(phi.hat = numeric(m), 
#                    mu.hat = I(matrix(NA, m,n)),
#                    stringsAsFactors=FALSE);
#   
#   for (i in 1:m) {
#     tmp = fit.nb.reg.1(y[i, ], s = s, x = x)
#     res[i,1] = tmp$phi.hat
#     res[i,-1] = tmp$mu.hat
#   }
  
  phi.hat = numeric(m)
  mu.hat = matrix(0, nr=m, nc=n)

  for (i in 1:m) {
    tmp = fit.nb.reg.1(y[i, ], s=s, x=x)
    phi.hat[i] = tmp$phi.hat
    mu.hat[i, ] = tmp$mu.hat
  }
  
  #return(list(phi.hat = res[,1], mu.hat = res[,-1]))
  return(list(phi.hat = phi.hat, mu.hat = mu.hat))
      
}



