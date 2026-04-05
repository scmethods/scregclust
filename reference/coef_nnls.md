# Compute NNLS coefficients

Computes non-negative least squares coefficients with a matrix right
hand side.

## Usage

``` r
coef_nnls(x, y, eps = 1e-12, max_iter = 1000L)
```

## Arguments

- x:

  Coefficient matrix (p x n matrix)

- y:

  Right hand side (p x m matrix)

- eps:

  Convergence tolerance

- max_iter:

  Maximum number of iterations

## Value

A list containing

- beta:

  The estimated coefficient matrix

- iterations:

  A vector containing the number of iterations needed for the `i`-th
  column in `y` in the `i`-th entry.

## References

Duy Khuong Nguyen and Tu Bao Ho. Accelerated anti-lopsided algorithm for
nonnegative least squares. International Journal of Data Science and
Analytics, 3(1):23–34, 2017.

Adapted from <https://github.com/khuongnd/nnls_antilopsided>
