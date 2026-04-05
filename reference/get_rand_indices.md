# Compute Rand indices

Compute Rand indices for fitted scregclust object

## Usage

``` r
get_rand_indices(fit, groundtruth, adjusted = TRUE)
```

## Arguments

- fit:

  An object of class `scregclust`

- groundtruth:

  A known clustering of the target genes (integer vector)

- adjusted:

  If TRUE, the Adjusted Rand index is computed. Otherwise the ordinary
  Rand index is computed.

## Value

A [`data.frame`](https://rdrr.io/r/base/data.frame.html) containing the
Rand indices. Since there can be more than one final configuration for
some penalization parameters, Rand indices are averaged for each fixed
penalization parameter. Returned are the mean, standard deviation and
number of final configurations that were averaged.

## References

W. M. Rand (1971). "Objective criteria for the evaluation of clustering
methods". Journal of the American Statistical Association 66 (336):
846–850. DOI:10.2307/2284239

Lawrence Hubert and Phipps Arabie (1985). "Comparing partitions".
Journal of Classification. 2 (1): 193–218. DOI:10.1007/BF01908075
