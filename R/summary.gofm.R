
#' @title Report summaries of the GOF test results for multivariate response
#'
#' @description A summary method of GOF test results based on a "gofm" 
#' object obtained from the \code{\link{nb_gof_m}} outputs
#' 
#' @param x an object of class "gofm" (from \link{nb_gof_m} output)
#' @param conf.env confidence level for the envelope
#' @param data.note a note on how data are simulated or the source of data
#' @param bin the number of bins used for calculating a standard chi-square p-value
#' @param ... for future use
#' 
#' @method summary gofm
#' @rdname summary.gofm
#' @export
#' 
#' @details This is a generic function for summarizing goodness-of-fit test results.
#' @usage summary(x, conf.env=0.95, data.note=NULL, ...)
#' 
#' @return Information of the dataset, simulation parameter specifications, exact binomial test
#' results, and Monte Carlo GOF p-values
#' 
#' @seealso See \code{\link{nb_gof_m}} for simulated data examples
#' 
#' @author Gu Mi, Yanming Di, Daniel Schafer
#' 
summary.gofm = function(x, 
                        conf.env=0.95,
                        data.note=NULL,
                        ...){
  
  ## quantities from a "gofm" object:
  model.fit = x$model.fit
  counts.dim = x$counts.dim
  xx = x$design.mat
  pv.Vert = x$pv.Vert
  pv.Pear = x$pv.Pear
  pval.pear.median = x$pval.pear.median
  sim.size = x$sim
  dist.mat = x$dist.mat
  pear.mat = x$pear.mat
  
  #### -----------------------------------------------------------------------------------------
  ## envelope method: Vertical distance
  alpha = 1-conf.env
  n.pts = dim(dist.mat)[2]  # i.e. number of genes, m
  quant.95.vert = apply(dist.mat[1:sim.size, ], 2, quantile, probs=conf.env)  # on sim. data only!
  out.ind.vert = dist.mat[(sim.size+1), ] > quant.95.vert
  # indicator if d(obs.) > d(sim.95th), j=1,...,m
  nout.vert = sum(out.ind.vert)
  p.hat.vert = nout.vert/n.pts
  # binom.test results: one-sided test
  bt.vert = binom.test(x=nout.vert, n=n.pts, p=alpha, alternative="greater",
                  conf.level=0.95)
  bt.pval.vert = round(bt.vert$p.value, 6)
  
  ## envelope method: Pearson statistics
  alpha = 1-conf.env
  n.pts = dim(dist.mat)[2]  # i.e. number of genes, m
  quant.95.pear = apply(pear.mat[1:sim.size, ], 2, quantile, probs=conf.env)  # on sim. data only!
  out.ind.pear = pear.mat[(sim.size+1), ] > quant.95.pear
  # indicator if d(obs.) > d(sim.95th), j=1,...,m
  nout.pear = sum(out.ind.pear)
  p.hat.pear = nout.pear/n.pts
  # binom.test results: one-sided test
  bt.pear = binom.test(x=nout.pear, n=n.pts, p=alpha, alternative="greater",
                       conf.level=0.95)
  bt.pval.pear = round(bt.pear$p.value, 6)
  
  #### -----------------------------------------------------------------------------------------
  ## Fisher's method:
  pv.fisher.vert = chisq_gof(x)$fisher.vert
  pv.fisher.pear = chisq_gof(x)$fisher.pear
  
  
  #### -----------------------------------------------------------------------------------------
  # summaries
  cat("--------------------------------------------------------------- \n")
  cat("| Data simulated: ", data.note, "\n")
  cat("| NB model used: ", model.fit, "\n")
  cat("| Simulation size: ", sim.size, "\n")
  cat("| Count matrix dimension = ", counts.dim, "\n")
  cat("| # pts. outside envelope (Pear.) = ", nout.pear, " (",round(p.hat.pear*100,2),"%)", 
      sep="", "\n")
  cat("| # pts. outside envelope (Vert.) = ", nout.vert, " (",round(p.hat.vert*100,2),"%)", 
      sep="", "\n")
  cat("| GOF exact binom.test p-value (Pear.Stat.) = ", bt.pval.pear, "\n")
  cat("| GOF exact binom.test p-value (Vert.Dist.) = ", bt.pval.vert, "\n")
  cat("| Monte Carlo GOF p-value (Pear.Stat.) = ", pv.Pear, "\n")
  cat("| Monte Carlo GOF p-value (Vert.Dist.) = ", pv.Vert, "\n")
  cat("| Monte Carlo GOF p-value (P.S.median) = ", pval.pear.median, "\n")
  cat("| Fisher Method GOF p-value (Pear.Stat.) = ", pv.fisher.pear, "\n")
  cat("| Fisher Method GOF p-value (Vert.Dist.) = ", pv.fisher.vert, "\n")
  cat("-------------------------------------------------------------- \n")  

}