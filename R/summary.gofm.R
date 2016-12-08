
#' @title Report Summaries of the GOF Test Results for NB Dispersion Models
#'
#' @description A summary method of GOF test results based on a "gofm" 
#' object obtained from the \code{\link{nb.gof.m}} outputs
#' 
#' @param x an object of class "gofm" from the \code{\link{nb.gof.m}} output
#' @param conf.env confidence level for the envelope
#' @param data.note a note on how data are simulated or the source of data
#' @param ... (for future use)
#' 
#' @method summary gofm
#' @rdname summary.gofm
#' @aliases summary
#' @export
#' 
#' @details This is a generic function for summarizing goodness-of-fit test results of testing
#' negative binomial dispersion models.
#' 
#' @usage summary(x, conf.env=0.95, data.note=NULL, ...)
#' 
#' @return Information of the dataset, simulation parameter specifications, exact binomial test
#' results, and Monte Carlo GOF p-values.
#' 
#' @seealso \code{\link{nb.gof.m}} for simulated data examples, and \code{\link{arab}} for 
#' a real RNA-Seq data example.
#' 
#' @author Gu Mi <neo.migu@gmail.com>, Yanming Di, Daniel Schafer
#' 
#' @references
#' See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
#'
summary.gofm = function(x, conf.env=0.95, data.note=NULL, ...){
  
  ## quantities from a "gofm" object:
  model = x$model
  counts.dim = x$counts.dim
  xx = x$design.mat
  sim.size = x$sim
  
  #### -----------------------------------------------------------------------------------------
  ## Fisher's method:
  pv.fisher.vert = chisq_gof(x)$fisher.vert
  pv.fisher.pear.1s = chisq_gof(x)$fisher.pear.1s
  
  #### -----------------------------------------------------------------------------------------
  # summaries
  cat("--------------------------------------------------------------- \n")
  cat("| Data simulated/used:", data.note, "\n")
  cat("| NB model used:", model, "\n")
  cat("| Simulation size:", sim.size, "\n")
  cat("| Count matrix dimension:", counts.dim, "\n")
  cat("| Fisher Method GOF p-value (Vert.Dist.(1s)) =", pv.fisher.vert, "\n")
  cat("| Fisher Method GOF p-value (Pear.Stat.(1s)) =", pv.fisher.pear.1s, "\n")
  cat("-------------------------------------------------------------- \n")  

}