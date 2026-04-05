# ADMM algorithm for solving the group-penalized least squares problem

Implements estimation of the coop-lasso problem.

## Usage

``` r
coop_lasso(
  y,
  x,
  lambda,
  weights,
  beta_0 = NULL,
  rho_0 = 0.2,
  alpha_0 = 1.5,
  n_update = 2L,
  eps_corr = 0.2,
  max_iter = 1000L,
  eps_rel = 1e-08,
  eps_abs = 1e-12,
  verbose = FALSE
)
```

## Arguments

- y:

  Target (n x m)

- x:

  Design matrix (n x p)

- lambda:

  Penalization parameter

- weights:

  A specific weight for each group (typically this is
  `sqrt(group size)`).

- beta_0:

  Initial value for coefficients, allowing for warm start. Can be set to
  NULL, which results in the initial `beta` being a zero matrix.

- rho_0:

  Initial ADMM step-size

- alpha_0:

  Initial ADMM relaxation parameter

- n_update:

  Number of steps in-between updates of the step-size/adaptation
  parameters

- eps_corr:

  Lower bound for the correlation in the step-size update steps

- max_iter:

  Maximum number of iterations

- eps_rel:

  Relative tolerance for convergence check

- eps_abs:

  Absolute tolerance for convergence check

- verbose:

  Whether or not information about the optimization process should be
  printed to the terminal

## Value

A list containing

- beta:

  The coefficients at convergence

- iterations:

  Number of iterations

## References

Xu et al. (2017) Adaptive relaxed ADMM: Convergence theory and practical
implementation. DOI 10.1109/CVPR.2017.765
