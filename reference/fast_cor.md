# Fast computation of correlation

This uses a more memory-intensive but much faster algorithm than the
built-in `cor` function.

## Usage

``` r
fast_cor(x, y)
```

## Arguments

- x:

  first input matrix

- y:

  second input matrix

## Value

Correlations matrix between the columns of `x` and `y`

## Details

Computes the correlation between the columns of `x` and `y`.
