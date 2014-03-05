NBGOF
=====

Goodness-of-Fit Tests and Model Diagnostics for Negative Binomial Regression of RNA sequencing Data

----------------------------------------------------------------------------------------------------

The NBGOF package implements goodness-of-fit (GOF) tests for negative binomial (NB) distributions and NB dispersion models, with applications in RNA-Seq data analysis. This package can be used to test the GOF of the NB2, NBP or Poisson **regression models**. It can also be used to test GOF for a variety of **NB dispersion models** in popular R/Bioconductor packages, including

* NBP dispersion model in the NBPSeq package (NBP)

* NBQ dispersion model in the NBPSeq package (NBQ)

* NB common dispersion model in the edgeR package (Common)

* NB genewise dispersion model in the edgeR package (Genewise)

* NB trended (non-parametric) dispersion model in the edgeR package (Trended)

* NB tagwise-common dispersion model in the edgeR package (Tagwise-Common)

* NB tagwise-trended dispersion model in the edgeR package (Tagwise-Trend)

To load the package into current R session, run

```{r}
library(NBGOF)
```
Two main functions for testing the adequacy of regression models and NB dispersion models are `nb.gof.v` and `nb.gof.m`. 
