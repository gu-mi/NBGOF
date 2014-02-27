\name{MDPlot}
\alias{MDPlot}
\title{Mean-Dispersion Plots with Curves of Fitted Negative Binomial Dispersion Models}
\usage{
  MDPlot(model.vec, counts, x, title=NULL)
}
\arguments{
  \item{model.vec}{a character vector specifying the names
  of the NB dispersion models. Currently supported include
  "NBP", "NBQ", "Common", "Tagwise-Common", "Tagwise-Trend"
  and "Trended"}

  \item{counts}{an m-by-n count matrix of non-negative
  integers. For a typical RNA-Seq experiment, this is the
  read counts with m genes and n samples}

  \item{x}{an n-by-p design matrix}

  \item{title}{title of the plot}
}
\value{
  A ggplot object of mean-dispersion plot with fitted NB
  dispersion model curves.
}
\description{
  This function makes mean-dispersion plots with curves of
  fitted negative binomial dispersion models. We currently
  consider only one group and fit an intercept-only model.
}
\examples{
# load package
library(NBGOF)

# load data and subset for illustration
data(arab)
counts = arab[1:5000,1:3]
x = matrix(c(1,1,1), nr=3)

# specify what models to plot
model.vec=c("Common", "NBP", "NBQ", "Trended", "Tagwise-Common", "Tagwise-Trend")

# begin plotting
pl.arab = MDPlot(model.vec, counts, x, title="Mean-Dispersion Plot with Fitted Dispersion Models (Arabidopsis Data)")
print(pl.arab)
}
\author{
  Gu Mi <mig@stat.oregonstate.edu>, Yanming Di, Daniel
  Schafer
}
\references{
  See \url{https://github.com/gu-mi/NBGOF/wiki/} for more
  details.
}

