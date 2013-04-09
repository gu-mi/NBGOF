
#' @title Report summaries of the GOF test results for multivariate response
#'
#' @description A summary method of GOF test results based on a "gofm" 
#' object obtained from the \code{\link{nb_gof_m}} outputs
#' 
#' @param x an object of class "gofm" (from \link{nb_gof_m} output)
#' @param conf.env confidence level for the envelope
#' @param data.note a note on how data are simulated or the source of data
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
summary.gofm <- function(x, 
                      conf.env=0.95,
                      data.note=NULL,
                      ...){
  
  ## quantities from a "gofm" object:
  model.fit <- x$model.fit
  counts.dim = x$counts.dim
  xx = x$design.mat
  pv.G <- x$orth.pval
  sim.size <- x$sim
  dist.mat = x$dist.mat
  
  alpha = 1-conf.env
  n.pts = dim(dist.mat)[2]
  quant.95 = apply(dist.mat, 2, quantile, probs=conf.env)
  out.ind = dist.mat[(sim.size+1), ] > quant.95
  nout = sum(out.ind)
  p.hat = nout/n.pts
  
  # binom.test results: one-sided test
  bt = binom.test(x=nout, n=n.pts, p=alpha, alternative="greater",
                  conf.level=0.95)
  bt.pval = bt$p.value
  bt.ci = as.vector(round(bt$conf.int,6))
  
  # summaries
  cat("--------------------------------------------------------------- \n")
  cat("| Data simulated: ", data.note, "\n")
  cat("| NB model used: ", model.fit, "\n")
  cat("| Simulation size: ", sim.size, "\n")
  cat("| Count matrix dimension = ", counts.dim, "\n")
  cat("| # pts. outside envelope = ", nout, " (",round(p.hat*100,2),"%)", sep="", "\n")
  cat("| 95% binom.test conf.int = ", "(",bt.ci[1]*100,"%,", bt.ci[2]*100,"%)", sep="", "\n")
  cat("| GOF exact binom.test p-value = ", round(bt.pval,6), "\n")
  cat("| Monte Carlo GOF p-value = ", round(pv.G, 6), "\n")
  cat("-------------------------------------------------------------- \n")

}