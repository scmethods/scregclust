# Compute OLS coefficients

If the design matrix has full column-rank, then use the normal least
squares estimate. Otherwise, use the Moore-Penrose inverse to compute
the least squares estimate.

## Usage

``` r
coef_ols(y, x)
```

## Arguments

- y:

  Target vector (n x 1)/matrix (n x m)

- x:

  Design matrix (n x p)

## Value

Vector of OLS coefficients
