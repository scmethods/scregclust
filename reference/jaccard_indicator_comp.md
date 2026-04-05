# Perform the computations for thresholded Jaccard distance

Perform the computations for thresholded Jaccard distance

## Usage

``` r
jaccard_indicator_comp(gs, eps)
```

## Arguments

- gs:

  a list of integer vectors, one for each row, giving the column indices
  of the non-zero elements of the row or `NULL` if the whole row is
  empty.

- eps:

  an upper bound on the Jaccard distance (`1 - eps` becomes a lower
  bound on the Jaccard similarity)

## Value

A list with row and column indices in the \#row x \#row indicator matrix
specifying which rows in the original matrix had a distance of at most
`eps`.

## Details

This function is optimized for sparse matrices and computes the pairwise
Jaccard distances between the rows of the input matrix. Note that the
actual distance is not saved. Instead, a threshold (`eps`) is supplied
and an indicator matrix is returned, with a one indicating that the
distance is smaller than `eps` (equivalently, the Jaccard similarity is
larger than `1 - eps`).
