
#' @title Multiple plot function
#' 
#' @description This function is used to plot multiple ggplot objects in a single page. The layout can be flexibly specified.
#' 
#' @param ... ggplot objects can be passed in ....
#' @param plotlist ggplot objects can also be passed to plotlist (as a list of ggplot objects).
#' @param cols number of columns in the layout.
#' @param layout a matrix specifying the layout. If present, 'cols' is ignored.
#' @param labs  global labels for x- and y-axes.
#' @param labpos  positions for the global labels.
#' 
#' @details If the layout is something like \code{matrix(c(1,2,3,3), nrow=2, byrow=TRUE)}, then plot 1 will go in the upper left, 
#' plot 2 will go in the upper right, and plot 3 will go all the way across the bottom.
#' 
#' @usage multiplot(..., plotlist=NULL, cols=1, layout=NULL, labs=list(), labpos=list(c(0.5, 0.02), c(0.02, 0.5)))
#' 
#' @seealso \code{\link{nb.gof.m}} for simulated data examples. Modifications of global labels made by Scott Chamberlain (StackOverflow)
#' 
#' @source http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#' 
#' @export
#' 
multiplot = function (..., plotlist = NULL, cols = 1, layout = NULL, labs=list(),
                      labpos=list(c(0.5, 0.01), c(0.01, 0.5))) {
  
  # make a list from the ... arguments and plotlist
  plots = c(list(...), plotlist)
  
  numPlots = length(plots)
  
  if (is.null(layout)) {
    # make the panel
    # ncol: number of columns of plots
    # nrow: number of rows of plots
    layout = matrix(seq(1, cols * ceiling(numPlots/cols)), 
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
  }
  
  else {
    # set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # make each plot, in the correct location
    for (i in 1:numPlots) {
      # get the i,j matrix positions of the regions that contain this subplot
      matchidx = as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, 
                                      layout.pos.col = matchidx$col))
    }
    
    # add global axis labels
    if (!length(labs) == 0){
      grid.text(labs[1], x=labpos[[1]][1], y=labpos[[1]][2], just="center", gp=gpar(fontsize=20, fontface="bold"))
      grid.text(labs[2], x=labpos[[2]][1], y=labpos[[2]][2], just="center",rot=90, gp=gpar(fontsize=20, fontface="bold"))
    }
  }
}