# Create a table of module overlap for two clusterings

Compares two clusterings and creates a table of overlap between them.
Module labels do not have to match.

## Usage

``` r
cluster_overlap(k1, k2)
```

## Arguments

- k1:

  First clustering

- k2:

  Second clustering

## Value

A matrix showing the module overlap with the labels of `k1` in the
columns and the labels of `k2` in the rows.
