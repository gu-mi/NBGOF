\name{multiplot}
\alias{multiplot}
\title{Multiple plot function}
\source{
  http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
}
\usage{
  multiplot(..., plotlist=NULL, cols=1, layout=NULL)
}
\arguments{
  \item{...}{ggplot objects can be passed in ....}

  \item{plotlist}{ggplot objects can also be passed to
  plotlist (as a list of ggplot objects).}

  \item{cols}{number of columns in the layout.}

  \item{layout}{a matrix specifying the layout. If present,
  'cols' is ignored.}
}
\description{
  This function is used to plot multiple ggplot objects in
  a single page. The layout can be flexibly specified.
}
\details{
  If the layout is something like matrix(c(1,2,3,3),
  nrow=2, byrow=TRUE), then plot 1 will go in the upper
  left, 2 will go in the upper right, and 3 will go all the
  way across the bottom.
}
\seealso{
  \code{\link{nb.gof.m}} for simulated data examples.
}

