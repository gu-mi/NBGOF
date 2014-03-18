
#' @title Prepare Quantities for Diagnostic Plots (Histogram and Uniform Q-Q Plots) Using ggplot2
#' 
#' @description This function prepares quantities for making diagnostic plots of negative binomial dispersion models. 
#' The result (a list) will be passed to the \code{\link{diagPlot}} function to make histogram of p-values and uniform Q-Q plots.
#' 
#' @param gofm.obj an object of class "gofm" from the \code{\link{nb.gof.m}} output.
#' 
#' @return a list of two components.
#' @keywords internal
#' 
diagdata = function(gofm.obj){
  
  # Produces x and y co-ordinates for uniform Q-Q plot
  # Arguments: gofm.obj: object of class "gofm" from nb.gof.m output
  #           
  # Output: List with 2 components:
  #   qq = data.frame with x and y co-ordinates of plot
  #   info = data.frame with name of tested model
  
  if ( class(gofm.obj) != "gofm") {
    stop("You must pass an object of class 'gofm' as the first argument!")
  }
  
  ## quantities from a "gofm" object:
  model = gofm.obj$model
  counts.dim = gofm.obj$counts.dim
  sim.size = gofm.obj$sim
  x = gofm.obj$v.pvals
  
  # in case there are any NA values in x
  ord.x = x[!is.na(x)][order(x[!is.na(x)])]
  n = length(ord.x)
  # ppoints(): used in qq plots to generate the set of probabilities at which to evaluate the inverse distribution
  P = ppoints(n)
  z = qunif(P)
  
  # construct data frames for quantities to be used
  qq = data.frame(px = z, py = ord.x)
  info = data.frame(model = model)
  
  # return as a list
  return(list(qq = qq, info = info))
}

#' @title Diagnostic Plots for Graphical Checks of Negative Binomial Dispersion Model Fits Using ggplot2
#' 
#' @description This function makes diagnostic plots based on the goodness-of-fit test results (an "gofm" object) from
#' \code{\link{nb.gof.m}}: a histogram of p-values or a uniform quantile-quantile (Q-Q) plot is available.
#' 
#' @param gofm.obj an object of class "gofm" from the \code{\link{nb.gof.m}} output.
#' @param type which diagnostic plot to make. Either \code{qq} for a uniform Q-Q plot, or \code{hist} for a histogram of p-values.
#' @param binwidth the bin width used when drawing a histogram.
#' 
#' @export
#' 
#' @usage diagPlot(gofm.obj, type="qq", binwidth=0.1)
#' 
#' @return A ggplot object, either a Q-Q plot (\code{type="qq"}) or a histogram (\code{type="hist"}).
#' 
#' @seealso \code{\link{nb.gof.m}} for simulated data examples.
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' 
#' @references See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
#' 
diagPlot = function(gofm.obj, type="qq", binwidth=0.1){
  
  if ( class(gofm.obj) != "gofm") {
    stop("You must pass an object of class 'gofm' as the first argument!")
  }
  if ( !type %in% c("qq", "hist") ) {
    stop("You must pass the name of plotting type as 'qq' or 'hist'!")
  }
  
  diagdata = diagdata(gofm.obj=gofm.obj)
  pv.fisher.vert = signif(chisq_gof(gofm.obj)$fisher.vert, 4)
  
  title = with(diagdata$info, paste(model))
  pv.note = paste(" = ", signif(pv.fisher.vert, 2), sep="")
  if (pv.fisher.vert < 0.0001) pv.note = " < 0.0001"
  if (pv.fisher.vert > 0.9999) pv.note = " > 0.9999"
  annot = paste("GOF p-value", pv.note, sep="")
  
  if (type == "qq") {
    qq = ggplot(data=diagdata$qq, aes(x = px, y = py)) + 
      geom_point(data=diagdata$qq, aes(x = px, y = py), colour = "Blue", alpha = 0.5, size=0.8, shape=1) +
      geom_abline(data=diagdata$info, aes(intercept = 0, slope = 1), linetype = "dotted") +
      theme_bw() +
      scale_x_continuous("", limits=c(0,1)) +
      scale_y_continuous("", limits=c(0,1)) +
      annotate("text", x = 0.38, y = 0.95, label=annot, size = 3.5) +
      #scale_colour_manual(values = "Blue", labels = eval(annot)) +
      ggtitle(title) + 
      theme(plot.title = element_text(face="bold", size=10),
            axis.title.x = element_text(face="bold", size=8),
            axis.title.y = element_text(face="bold", size=8, angle=90),
            axis.text.x  = element_text(size=8),
            axis.text.y  = element_text(size=8),
            plot.margin=unit(x=c(0,0.2,0,0.2), units="cm"),  # top, right, bottom, and left
            legend.justification=c(0,1), 
            legend.position=c(0,1),
            legend.title = element_blank(),
            legend.key = element_blank()
      )
    
    return(qq)
  }
  
  if (type == "hist"){
    title = with(diagdata$info, paste("Histogram of P-values for ", model, " Model", sep=""))
    hist = ggplot(data=diagdata$histogram, aes(x=pvals)) +
      geom_histogram(colour = "black", fill = "white", binwidth = binwidth) + 
      ggtitle(title) + theme_bw() +
      scale_x_continuous("P-values", limits = c(0,1)) + 
      scale_y_continuous("Frequency") +
      theme(plot.title = element_text(face="bold", size=10),
            axis.title.x = element_text(face="bold", size=8),
            axis.title.y = element_text(face="bold", size=8, angle=90),
            axis.text.x  = element_text(size=8),
            axis.text.y  = element_text(size=8),
            plot.margin=unit(x=c(0,0.2,0,0.2), units="cm"),  # top, right, bottom, and left
            legend.justification=c(0,1), 
            legend.position=c(0,1),
            legend.title = element_blank(),
            legend.key = element_blank()
      )
    return(hist)
  }
}
