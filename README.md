# NBGOF
=============================================
*********************************************

The R package **NBGOF** implements the `NBGOF` approach discussed in ... 

## NBGOF

### Model Formulations

We assume that the responses $Y_{ij}$, such as read counts of genes/exons, are presented in the form of a $m$-by-$n$ matrix. In a typical RNA-Seq experiment, each row of such matrix represents a gene and each column a replicate of a (biological) sample.

\textbf{NB2 Model}\\
The \texttt{glm.nb} function in the \texttt{MASS} package takes a vector as the response and fits an NB2 regression model with a common dispersion $\phi$. This function has convergence issues (\texttt{theta.ml} $\rightarrow\infty$), likely to occur when the maximum likelihood function is not strictly concave or badly behaved. Therefore, we switch to the \texttt{nb.regression.1} function for getting the MLE of $\phi_i$ for each row ($i=1,\cdots,m$), but adjusted profile likelihood estimates can also be obtained. In that case, the results should be identical to \texttt{edgeR} tagwise dispersion model with the prior degree of freedom $G_0=0$.

\textbf{NBP Model}\\
In the \texttt{NBPSeq} package, the \texttt{estimate.dispersion} function pools information across genes and biological samples by modeling $\phi_{ij}$ as a function of relative mean frequencies. The ``log-linear-rel-mean'' method assumes a parametric dispersion model as discussed in Equation \ref{eq:disp-ll}. It estimates different dispersions for each gene/sample combination (the same dispersion for all replicates in a sample but differ across samples). The function also works when a sample contains read counts close to zero (or even all zeros), but the dispersion estimates would be rather high. A first-stage screening of genes with very few counts might be necessary for more stable results. Such screening procedure is not uncommon in RNA-Seq analysis. In \citet{lund2012detecting}, for example, the authors filtered out genes with average counts (1) less than or equal to 1 or (2) fewer than 2 samples with non-zero counts.

The \texttt{glm.nbp.1} function takes a vector as the response and simultaneously estimates $(\beta,\alpha_0,\alpha_1)$. In the case of moderate-sized replicates for each gene, this function should be used with caution: though the coefficients $\beta$ may be accurately estimated, the dispersion estimates can be far off.

\textbf{Dispersion Models in \texttt{edgeR}}\\
For the \texttt{edgeR} common dispersion model, there are three methods available for estimating the dispersion $\phi$. The default \texttt{dispCoxReid} uses the Cox-Reid adjusted profile likelihood estimator which leads to the least biased estimate. \texttt{dispPearson} uses the pseduo-likelihood (PL) estimator which tends to under-estimate. \texttt{dispDeviance} uses the quasi-likelihood (QL) estimator which tends to over-estimate. We use the default setting so $\phi$ is estimated by the APLE. In the paper, this approach is termed ``\texttt{edgeR}-tagwise'', to be distinguished from the APL estimator without shrinkage ($G_0=0$, ``\texttt{edgeR}-genewise'').

## References
