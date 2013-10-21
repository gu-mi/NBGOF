
#' @title Calculations of coefficient of variations (CV) for subsetting read counts 
#'
#' @description This function returns a list of (1) matrix of CV for each treatment group; (2) gene index vector
#' 
#' @param counts raw RNA-Seq read counts matrix.
#' @param grp.ids a vector (converted to factor) indicating group memebership of samples.
#' @param cutoff the threshold of CV used to subset genes. Default is no subset.
#' 
#' @export
#' 
#' @usage CV(counts, grp.ids=NULL, cutoff=Inf)
#' 
#' @return A list of (1) matrix of CV for each treatment group; (2) gene index vector
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 

CV = function(counts, grp.ids=NULL, cutoff=Inf){
  m = dim(counts)[1]
  n.grp = length(unique(grp.ids))
  cv.mat = matrix(0, nr=m, nc=n.grp)
  gene.idx = numeric(m) + 1
  cv.tmp = t(apply(counts, 1, cv.vec, grp.ids))
  for (i in 1:m){
    cv.mat[i, ] = unique(cv.tmp[i, ])
    # if specify a cutoff of CV, then return gene indices of TRUE or FALSE (if any one of the treatment groups does)
    if (any(cv.mat[i, ] > cutoff)){
      gene.idx[i] = 0
    }
  }
  return(list(cv.mat = cv.mat, gene.idx = as.logical(gene.idx)))
}