printf = function(...) cat(sprintf(...));

check.equal = function(x, y) {
  id.nax = is.na(x);
  id.nay = is.na(y);
  if (any(id.nax != id.nay)) return(FALSE);

  id = !id.nax;

  return(range(x[id] - y[id]));
}

## An alternative to plot.default() for plotting a large number of
## densely distributed points.  This function can produce a visually
## almost identical plot using only a subset of the points.  This is
## particular useful for reducing output file size when plots are
## written to eps files.

smart.plot = function(x, y=NULL, xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL,
  log="", resolution=500, plot=TRUE, col=NULL, clip=Inf, color.clipped=TRUE, ...) {

  ## Arguments:
  ##
  ##   resolution: determines the distance below which points will be
  ##   considered as overlapping.
  ##
  ##   clip
  ##
  ##   use.color
  ##
  ##   other arguments are the same as in plot.default().
  ##
  ##
  ## Values (if plot=FALSE):
  ##
  ##   x, y: the x, y-coordinates of the subset of representative points
  ## 
  ##   id: the indicies of these points in the original data set
  ##  
  ##   freqs: the numbers of points that overlap with each representative point
  ##
  ##   col: colors determined by the freqs

  ## Details:
  ##
  ##   Writing plots with a large number of points to eps files can
  ##   result in big files and lead to very slow rendering time.
  ##
  ##   Usually for a large number of points, a lot of them will
  ##   overlap with each other. Plotting only a subset of selected
  ##   non-overlapping points can give visually almost identical
  ##   plots. Further more, the plots can be enhanced if using gray
  ##   levels (the default setting) that are proportional to the
  ##   number points overlapping with each plotted point.
  ## 
  ##   This function scans the points sequentially. For each unmarked
  ##   point that will be plotted, all points that overlap with it
  ##   will be marked and not to plotted, and the number of
  ##   overlapping points will be recorded. This is essentially
  ##   producing a 2d histogram. The freqs of the points will be
  ##   converted to gray levels, darker colors correspond to higher
  ##   freqs.


  ## These lines are copied from plot.default.
  xlabel = if (!missing(x)) deparse(substitute(x));
  ylabel = if (!missing(y)) deparse(substitute(y));
  xy = xy.coords(x, y, xlabel, ylabel, log);
  xlab = if (is.null(xlab)) xy$xlab else xlab;
  ylab = if (is.null(ylab)) xy$ylab else ylab;
  xlim = if (is.null(xlim)) range(xy$x[is.finite(xy$x)]) else xlim;
  ylim = if (is.null(ylim)) range(xy$y[is.finite(xy$y)]) else ylim;
  
  x = xy$x;
  y = xy$y;
  n = length(x);
  id = is.finite(xy$x) & is.finite(xy$y);
  id[id & (x < xlim[1] | x > xlim[2] | y < ylim[1] | y > ylim[2])]=FALSE;

  logxy = strsplit(log, NULL)[[1L]];
  if ("x" %in% logxy) {
    x[id] = log(x[id]);
    epsx = diff(log(xlim)) / resolution;
  } else
    epsx = diff(xlim) / resolution;

  if ("y" %in% logxy) {
    y[id] = log(y[id]);
    epsy = diff(log(ylim)) / resolution; 
  } else
    epsy = diff(ylim) / resolution; 

  counts = rep(0, n);
  i = 1;

  ## Scan the points and select non-overlapping ones
  while (i < n) {
    if (id[i]) {
      ids = ((i+1):n)[id[(i+1):n]];
      overlap = ids[abs(x[ids] - x[i]) < epsx & abs(y[ids] - y[i]) < epsy];
      id[overlap] = FALSE;
      counts[i] = length(overlap) + 1;
    }
    i = i + 1;
  }

  ## Sort the data so that points representing more points will be plotted later.
  id = (1:n)[id];
  id = id[order(counts[id])];
  counts = counts[id];

  if (is.null(col)) {
    ## Convert counts of overlapping points to gray levels
    counts.clipped = counts;
    id.clipped = counts > clip;
    counts.clipped[id.clipped] = clip;
    col = gray((1 - counts.clipped / max(counts.clipped)) * 0.8);
    if (color.clipped) {
      col[id.clipped] = rgb(counts[id.clipped]/max(counts), 0, 0);
    }
  }

  if (plot) {
    plot(x=xy$x[id], y=xy$y[id], xlim=xlim, ylim=ylim, log=log, xlab=xlab, ylab=ylab, col=col, ...);
    invisible();
  } else {
    list(x=xy$x[id], y=xy$y[id], freqs=counts, col = col, id=id);
  }
}


smart.pairs = function(x, log.hist=FALSE, ...) {
  n = dim(x)[2];

  par(mfrow=c(n, n));

  for (i in 1:n) {
    for (j in 1:n) {
      if (i==j) {
        if (log.hist) {
          hist(log(x[,i]));
        } else {
          hist(x[,i]);
        }
      } else {
        smart.plot(x[,c(i,j)], ...);
        abline(0, 1);
      }
    }
  }

  invisible();
}

example.smart.filter = function() {
  x = rnorm(20000);
  y = rnorm(20000);

  ## debug(smart.plot);
  ## debug(plot.default);
  ## undebug(plot.default);
  par(mfrow=c(2,3));
  plot(x, y, pch=19);
  ## Plot with translucent color
  plot(x, y, col=rgb(0, 0, 1, 0.1), pch=19);
  smart.plot(x, y, resolution=100, pch=19);
  smart.plot(x, y, pch=19);
  smart.plot(x, y, res = 100, pch="+");
  smart.plot(x, y, res = 100, pch="+", col=1);

  smart.plot(x, y, resolution=100, clip=5, pch=19);
  smart.plot(x, y, resolution= 100, pch=19, log="xy");
  obj = smart.plot(x, y, resolution= 100, pch=19, log="xy", plot=FALSE);

}
