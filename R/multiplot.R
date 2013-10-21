

#' @title Multiple plot function
#' 
#' @description This function is used to plot multiple ggplot objects in a single page. The layout can be flexibly specified.
#' 
#' @param ... ggplot objects can be passed in ....
#' @param plotlist ggplot objects can also be passed to plotlist (as a list of ggplot objects).
#' @param cols number of columns in the layout.
#' @param layout a matrix specifying the layout. If present, 'cols' is ignored.
#' @details If the layout is something like \code{matrix(c(1,2,3,3), nrow=2, byrow=TRUE)}, then plot 1 will go in the upper left, 
#' plot 2 will go in the upper right, and plot 3 will go all the way across the bottom.
#' 
#' @usage multiplot(..., plotlist=NULL, cols=1, layout=NULL)
#' 
#' @seealso \code{\link{nb.gof.m}} for simulated data examples.
#' @source http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#' 
#' @export
#' 
multiplot = function(..., plotlist=NULL, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots = c(list(...), plotlist)
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout = matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
  } 
  else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx = as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}