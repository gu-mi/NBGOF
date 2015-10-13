NBGOF
=====

### Goodness-of-Fit Tests and Model Diagnostics for Negative Binomial Regression of RNA sequencing Data

`NBGOF` is an R package for implementing goodness-of-fit (GOF) tests for negative binomial (NB) distributions and NB dispersion models, with applications in RNA-Seq data analysis. This package can be used to test the GOF of the NB2, NBP or Poisson regression models. It can also be used to test GOF for a variety of NB dispersion models in popular R/Bioconductor packages, including

* NBP dispersion model in the NBPSeq package (NBP)

* NBQ dispersion model in the NBPSeq package (NBQ)

* NB common dispersion model in the edgeR package (Common)

* NB genewise dispersion model in the edgeR package (Genewise)

* NB trended (non-parametric) dispersion model in the edgeR package (Trended)

* NB tagwise-common dispersion model in the edgeR package (Tagwise-Common)

* NB tagwise-trended dispersion model in the edgeR package (Tagwise-Trend)

The methodologies are discussed in the manuscript **Goodness-of-Fit Tests and Model Diagnostics for Negative Binomial Regression of RNA Sequencing Data** (by Gu Mi, Yanming Di, and Daniel W. Schafer, PLOS ONE, 10(3)). The paper is freely available from this [link](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0119254). Functions in `NBGOF` have been used to generate all figures and tables displayed in the manuscript, plus some additional analysis tools for further investigations.

******

**The `NBGOF` source package can be downloaded [here](https://www.dropbox.com/s/4pb3wfiee48o6v7/NBGOF_0.1.8.tar.gz?dl=0).**

### Install
Installation from this repository requires the `devtools` package pre-installed. Run the following R command to install `NBGOF`:

```S
devtools::install_github('NBGOF', 'gu-mi')
```

To load the package into current R session, run

```S
library(NBGOF)
```

Two main functions for testing the adequacy of regression models and NB dispersion models are `nb.gof.v` and `nb.gof.m`, repectively. We provide two real datasets (`arab` and `earthquake`) to illustrate the use of this package.

******

We provide (Dropbox) links below to download R source codes and related supporting files for preparing the datasets and reproducing figures/tables in the manuscript. Some intermediate key results are also provided when necessary.

### Figures

* Figure 1 -- The mean-dispersion plot with six fitted dispersion models (common, NBP, NBQ, trended, tagwise-common and tagwise-trend) for the Arabidopsis RNA-Seq dataset (19,623 genes from three biological samples in the mock treatment group): [Download R files and results](https://www.dropbox.com/sh/x5quzc102xnjhqe/AADCi_1nBqaYwDlV7OIV0Bi_a?dl=0)
* Figure 2 -- Empirical probability plots with GOF test p-values for evaluating NB2 and NBP model fits on the earthquake dataset: [Download R files and results](https://www.dropbox.com/sh/kw4u8i2d3k4ie2h/AAAFmDwivEHSAuCpSrX2eKS-a?dl=0)
* Figure 3 -- Empirical probability plots and GOF p-values for testing NB2 (top row) and NBP (bottom) on four simulated datasets with sample size = 45: [Download R files and results](https://www.dropbox.com/sh/zen2m2cjsg1zh7v/AAAAXmj605u8hO1-dKcPFovva?dl=0)
* Figure 4 -- Uniform QQ plots of individual GOF test p-values for the Arabidopsis dataset (based on a random sample of 1,000 genes from six experimental units in two experimental groups): [Download R files and results](https://www.dropbox.com/sh/nvrbltddxnmkpjt/AABFYEO1SaEZD7LIm5Bf3UJGa?dl=0)
* SI Figure 1 -- Uniform QQ plots of individual GOF test p-values for the simulated ``NB2+noise'' dataset: [Download R files and results](https://www.dropbox.com/sh/lbipomac0s9t7gv/AAC7OhlfxyxL8KZelgOPzaVpa?dl=0)


### Tables
* Table 1 -- Type I error rates for 0.05-level NB2 and NBP GOF tests: [Download R files](https://www.dropbox.com/sh/131q2kw0skoy3du/AABYvOBGEPodEjOqB9cUbNlia?dl=0)
* Table 2 -- Rejection rates for 0.05-level NB2 and NBP GOF tests: [Download R files](https://www.dropbox.com/sh/iidjc6x00ufzl0m/AADzJzZwneq203TFTFn9ipTUa?dl=0)

******

All results (figures/tables) in the submitted manuscript was produced using the following versions of R and its packages:

```{r}
R version 3.1.0 Patched (2014-05-22 r65728)
Platform: x86_64-apple-darwin10.8.0 (64-bit)

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] plyr_1.8.1      edgeR_3.6.2     limma_3.20.7    NBPSeq_0.3.0    NBGOF_0.1.6     ggplot2_1.0.0   doMC_1.3.3     
[8] iterators_1.0.7 foreach_1.4.2  

loaded via a namespace (and not attached):
 [1] codetools_0.2-9   colorspace_1.2-4  digest_0.6.4      gtable_0.1.2      MASS_7.3-35       munsell_0.4.2    
 [7] numDeriv_2012.9-1 proto_0.3-10      Rcpp_0.11.3       reshape2_1.4      scales_0.2.4      splines_3.1.0    
[13] stringr_0.6.2     tools_3.1.0      
```

******

If you have any questions, please do not hesitate to email the repository maintainer (Gu Mi) at neo.migu@gmail.com. Thank you for your interests in our research work.
