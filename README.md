# Single-cell Regulatory-driven Clustering (scregclust)

<!-- badges: start -->

<!-- badges: end -->

![A diagram illustrating the *scregclust* algorithm.](man/figures/overview_fig1A_bg.png "Illustration of the scregclust algorithm")

The goal of *scregclust* is to cluster genes by regulatory programs. To do so, clusters are associated with regulatory programs and target genes are allocated to clusters with best fitting regulatory programs.

-   The documentation for this package can be found at [https://sven-nelander.github.io/scregclust](https://sven-nelander.github.io/scregclust)
-   A detailed description of the algorithm and an in-depth evaluation of its properties can be found in our [pre-print on bioRxiv](https://doi.org/10.1101/2023.03.10.532041 "Reconstructing the regulatory programs underlying the phenotypic plasticity of neural cancers"){.uri}.

## Installation

You can install the development version of *scregclust* from [GitHub](https://github.com/sven-nelander/scregclust) with:

```r
# install.packages("devtools")
devtools::install_github("sven-nelander/scregclust")
```
