# Single-cell Regulatory-driven Clustering (scregclust)

<!-- badges: start -->

<!-- badges: end -->

![A diagram illustrating the *scregclust* algorithm.](man/figures/overview_fig1A_bg.png "Illustration of the scregclust algorithm")

The goal of *scregclust* is to cluster genes by regulatory programs. To do so, genes are clustered into modules which in turn are associated with regulators. The algorithm alternates between associating regulators to modules and reallocating target genes into modules.

- The documentation for this package can be found at [https://scmethods.github.io/scregclust/](https://scmethods.github.io/scregclust/)
- A detailed description of the algorithm and an in-depth evaluation of its properties can be found in our original research article [Larsson, Held, et al. (2024) Reconstructing the regulatory programs underlying the phenotypic plasticity of neural cancers. Nature Communications 15, 9699 DOI 10.1038/s41467-024-53954-3](https://doi.org/10.1038/s41467-024-53954-3)

## Installation

You can install the stable version of *scregclust* from [CRAN](https://cran.r-project.org/package=scregclust) with

```r
install.packages("scregclust")
```

You can install the current development version of *scregclust* from [GitHub](https://github.com/scmethods/scregclust) with:

```r
# install.packages("devtools")
devtools::install_github("scmethods/scregclust")
```
