
#' @title Prepare Quantities for Diagnostic Plots (Histogram and Uniform Q-Q Plots) Using ggplot2
#' @description This function prepares quantities for making diagnostic plots of negative binomial dispersion models. 
#' The result (a list) will be passed to the \code{\link{diagPlot}} function to make histogram of p-values and uniform Q-Q plots.
#' @param gofm.obj an object of class "gofm" from the \code{\link{nb.gof.m}} output.
#' @param envelope the level of confidence band when drawing a uniform Q-Q plot.
#' @param binwidth the bin width used when drawing a histogram.
#' @return a list of three components.
#' @keywords internal
#' 
plotdata = function(gofm.obj, envelope=0.95, binwidth=0.1){
  
  # Produces x and y co-ordinates for uniform Q-Q plot
  # Arguments: gofm.obj: object of class "gofm" from nb.gof.m output
  #           
  # Output: List with 3 components:
  #   qq = data.frame with x and y co-ordinates of plot; lower confident limit (lcl) and upper confident limit (ucl)
  #   info = data.frame with name of tested model, intercept and slope of reference line, and percentage of points outside envelope
  #   histogram = data.frame with quantities for making a histogram
  
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
  
  # pointwise confidence envelope
  conf = ifelse(envelope == FALSE, 0.95, envelope)
  crit.val = qnorm(1 - (1 - conf)/2)   # ~1.96 if conf = 0.95
  #
  Q.x = quantile(ord.x, c(.25,.75))
  Q.z = qunif(c(.25,.75))
  b = (Q.x[2] - Q.x[1])/(Q.z[2] - Q.z[1])
  a = Q.x[1] - b*Q.z[1]
  v = (b/dunif(z)) * sqrt(P*(1 - P)/n)
  m = a + b*z
  #
  lcl = m - crit.val * v  # lower confidence limit (lcl)
  ucl = m + crit.val * v  # upper confidence limit (ucl)
  n.out = sum(ord.x >= ucl) + sum(ord.x<=lcl)
  perc = round(n.out/n * 100, 2)
  
  # construct data frames for quantities to be used
  qq = data.frame(px = z, py = ord.x, lcl = lcl, ucl = ucl)
  info = data.frame(model = model, a=a[[1]], b=b[[1]], perc = perc)
  histogram = data.frame(pvals = x)
  
  # return as a list
  return(list(qq = qq, info = info, histogram = histogram))
}

#' @title Diagnostic Plots for Graphical Checks of Negative Binomial Dispersion Models
#' 
#' @description This function makes diagnostic plots based on the goodness-of-fit test results (an "gofm" object) from
#' \code{\link{nb.gof.m}}: a histogram of p-values or a uniform quantile-quantile (Q-Q) plot is available.
#' 
#' @param gofm.obj an object of class "gofm" from the \code{\link{nb.gof.m}} output.
#' @param type which diagnostic plot to make. Either \code{qq} for a uniform Q-Q plot, or \code{hist} for a histogram of p-values.
#' @param envelope the level of confidence band when drawing a uniform Q-Q plot.
#' @param binwidth the bin width used when drawing a histogram.
#' 
#' @export
#' 
#' @usage diagPlot(gofm.obj, type="qq", envelope=0.95, binwidth=0.1)
#' 
#' @return A ggplot object, either a histogram (\code{type="hist"}) or a Q-Q plot (\code{type="qq"}).
#' 
#' @seealso \code{\link{nb.gof.m}} for simulated data examples.
#' 
#' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
#' @references See \url{https://github.com/gu-mi/NBGOF/wiki/} for more details.
#' 
diagPlot = function(gofm.obj, type="qq", envelope=0.95, binwidth=0.1){
  
  plotdata = plotdata(gofm.obj, envelope=envelope)
  
  if (type == "qq"){
    title = with(plotdata$info, paste("Uniform Q-Q Plot of ", model, " Model", sep=""))
    annot = with(plotdata$info, paste(perc, "% Outside ", eval(envelope*100), "% Confidence Band", sep=""))
    
    # begin plotting
    qq = ggplot(data=plotdata$qq, aes(x = px, y = py)) + 
      geom_point(data=plotdata$qq, aes(x = px, y = py, colour = "Blue"), alpha = 1, size=2) +
      geom_line(data=plotdata$qq, aes(x = px, y = lcl), linetype="dashed", colour="Black") +
      geom_line(data=plotdata$qq, aes(x = px, y = ucl), linetype="dashed", colour="Black") + 
      geom_abline(data=plotdata$info, aes(intercept = a, slope = b), linetype = "dotted") +
      theme_grey() +
      scale_x_continuous("Theoretical Uniform Quantiles") +
      scale_y_continuous("P-value Quantiles") +
      scale_colour_manual(values = "Blue", labels = eval(annot)) +
      ggtitle(title) + 
      theme(plot.title = element_text(face="bold", size=14),
            axis.title.x = element_text(face="bold", size=12),
            axis.title.y = element_text(face="bold", size=12, angle=90),
            # panel.grid.major = element_blank(),
            # panel.grid.minor = element_blank(),
            legend.justification=c(1,0), 
            legend.position=c(1,0),
            legend.title = element_blank(),
            legend.key = element_blank()
      )
  }
  else if (type == "hist"){
    title = with(plotdata$info, paste("Histogram of P-values for ", model, " Model", sep=""))
    hist = ggplot(data=plotdata$histogram, aes(x=pvals)) +
      geom_histogram(colour = "black", fill = "white", binwidth = binwidth) + 
      ggtitle(title) + theme_grey() +
      scale_x_continuous("P-values", limits = c(0,1)) + 
      scale_y_continuous("Frequency") +
      theme(plot.title = element_text(face="bold", size=14),
            axis.title.x = element_text(face="bold", size=12),
            axis.title.y = element_text(face="bold", size=12, angle=90)
            # panel.grid.major = element_blank(),
            # panel.grid.minor = element_blank()
            )
  }
}


# # test
# qp0 = diagPlot(fnbp.arab.1g, type="hist", binwidth=0.1)
# print(qp0)
# qp1 = diagPlot(fnbp.arab.1g, type="qq", envelope=0.95)
# print(qp1)
# qp2 = diagPlot(fcom.arab.1g)
# print(qp2)
# qp3 = diagPlot(ftrd.arab.1g, envelope=0.99)
# print(qp3)
# multiplot(qp0, qp1, qp2, qp3, cols=2, layout=matrix(c(1,2,3,4), nrow=2, byrow=TRUE))


# #' @title Uniform Quantile-Quantile (Q-Q) Plot of Individual P-values for Model Diagnostics
# #' 
# #' @description This function provides effective visual diagnostics of the fitted NB dispersion model, 
# #' by drawing a uniform quantile-quantile plot of p-values from individual GOF p-value of each gene.
# #' 
# #' @param gofm.obj an object of class "gofm" from the \code{\link{nb.gof.m}} output
# #' @param col color of the points on the Q-Q plot
# #' @param col.lines color of the pointwise confidence envelope
# #' @param lwd line width
# #' @param pch point style
# #' @param cex point size
# #' @param lty line style
# #' @param envelope confidence level of the envelope; if FASLE, then no drawing
# #' @param xlab x-axis label
# #' @param ylab y-axis label
# #' @param grid whether to draw grids on background
# #' @param ... (for future use)
# #' 
# #' @return A Q-Q plot, with percentage of points outside the confidence envelope shown in legend.
# #' 
# #' @usage QQPlot(gofm.obj, col=palette()[1], col.lines=palette()[2], lwd=2, pch=".", cex=5, lty="dashed", 
# #' envelope=0.95, xlab="Expected", ylab="Observed", grid=TRUE, ...)
# #' 
# #' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di, Daniel Schafer
# #' @export
# #' 
# #' @references \url{https://github.com/gu-mi/NBGOF/wiki/}
# #' 
# QQPlot = function(gofm.obj, col=palette()[1], col.lines=palette()[2], lwd=2, pch=".", cex=5, lty="dashed",
#                   envelope=0.95, xlab="Expected", ylab="Observed", grid=TRUE, ...) {
#   
#   ## quantities from a "gofm" object:
#   model = gofm.obj$model
#   counts.dim = gofm.obj$counts.dim
#   sim.size = gofm.obj$sim
#   x = gofm.obj$v.pvals
#   
#   # in case there are any NA values in x
#   ord.x = x[!is.na(x)][order(x[!is.na(x)])]
#   n = length(ord.x)
#   # ppoints(): used in qq plots to generate the set of probabilities at which to evaluate the inverse distribution
#   P = ppoints(n)
#   z = qunif(P)
#   #
#   plot(z, ord.x, type="n", xlab=xlab, ylab=ylab, main=paste("Q-Q Uniform Plot of ", model, " Model", sep=""))
#   # grid on background
#   if(grid){
#     grid(lty=1, equilogs=FALSE)
#     box()
#   }
#   # draw points on Q-Q plot
#   points(z, ord.x, col=col, pch=pch, cex=cex)
#   
#   # pointwise confidence envelope
#   conf = ifelse(envelope == FALSE, 0.95, envelope)
#   crit.val = qnorm(1 - (1 - conf)/2)   # ~1.96 if conf = 0.95
#   
#   Q.x = quantile(ord.x, c(.25,.75))
#   Q.z = qunif(c(.25,.75))
#   b = (Q.x[2] - Q.x[1])/(Q.z[2] - Q.z[1])
#   a = Q.x[1] - b*Q.z[1]
#   v = (b/dunif(z)) * sqrt(P*(1 - P)/n)
#   m = a + b*z
#   
#   #   SE = (b/dunif(z, ...)) * sqrt(P*(1 - P)/n)
#   #   fit.value = a + b*z
#   #   upper = fit.value + zz*SE
#   #   lower = fit.value - zz*SE
#   #   if (envelope != FALSE) {
#   #     lines(z, upper, lty=2, lwd=lwd, col=col.lines)
#   #     lines(z, lower, lty=2, lwd=lwd, col=col.lines)
#   #   }
#   
#   #m = (1:n)/(n+1)
#   #v = sqrt((1:n) * (n - (1:n) + 1)/(n + 1)^2/(n + 2))
#   
#   lcl = m - crit.val * v  # lower confidence limit (lcl)
#   ucl = m + crit.val * v  # upper confidence limit (ucl)
#   
#   #   lid = (lcl > 0)    # valid id for lower limit
#   #   uid = (ucl <= 1)   # valid id for upper limit
#   #   la = z[lid]
#   #   lb = lcl[lid]
#   #   uc = z[uid]
#   #   ud = ucl[uid]
#   
#   # draw lines on Q-Q plot
#   lines(z, lcl, lty=lty, col=col.lines)   # lower limit
#   lines(z, ucl, lty=lty, col=col.lines)   # upper limit
#   abline(a, b, col=col.lines, lwd=lwd)
#   
#   # legend: percentage of points outside the confidence band
#   n.out = sum(ord.x >= ucl) + sum(ord.x<=lcl)
#   perc = round(n.out/n * 100, 2)
#   legend("topleft", bty="n", legend=paste("Percentage Outside = ", perc, "%", sep=""))
#   
# }


