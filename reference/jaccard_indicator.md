# Compute indicator matrix of pairwise distances smaller than threshold

Computes the Jaccard distance between rows of a matrix and returns a
sparse symmetric indicator matrix containing the entries with a distance
of less than a given upper bound. Note that the diagonal is always 1.

## Usage

``` r
jaccard_indicator(x, upper_bnd = 0.8)
```

## Arguments

- x:

  the input matrix with vectors to be compared in the rows.

- upper_bnd:

  pairs with a Jaccard distance below this upper bound are returned as 1
  while all others receive the entry 0.

## Value

A list of vectors describing a sparse lower triangular pattern matrix

- i:

  Row indices

- j:

  Column indices
