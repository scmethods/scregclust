# Perform the k-means++ algorithm

Performs the k-means++ algorithm to cluster the rows of the input
matrix.

## Usage

``` r
kmeanspp(x, n_cluster, n_init_clusterings = 10L, n_max_iter = 10L)
```

## Arguments

- x:

  Input matrix (n x p)

- n_cluster:

  Number of clusters

- n_init_clusterings:

  Number of repeated random initializations to perform

- n_max_iter:

  Number of maximum iterations to perform in the k-means algorithm

## Value

An object of class
[`stats::kmeans`](https://rdrr.io/r/stats/kmeans.html).

## Details

Estimation is repeated

## References

David Arthur and Sergei Vassilvitskii. K-Means++: The advantages of
careful seeding. In Proceedings of the Eighteenth Annual ACM-SIAM
Symposium on Discrete Algorithms, SODA '07, pages 1027––1035. Society
for Industrial and Applied Mathematics, 2007.
