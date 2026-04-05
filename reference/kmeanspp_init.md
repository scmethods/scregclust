# Determine initial centers for the kmeans++ algorithm

Determine initial centers for the kmeans++ algorithm

## Usage

``` r
kmeanspp_init(n_cluster, x = NULL, dm = NULL)
```

## Arguments

- x:

  data matrix to be clustered

- dm:

  distance matrix (between rows of x; of class "dist")

## Value

Row indices of initial cluster centers of x
