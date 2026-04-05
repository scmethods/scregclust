# Uncover gene modules and their regulatory programs from single-cell data

Use the scRegClust algorithm to determine gene modules and their
regulatory programs from single-cell data.

## Usage

``` r
scregclust(
  expression,
  genesymbols,
  is_regulator,
  penalization,
  n_modules,
  initial_target_modules = NULL,
  sample_assignment = NULL,
  center = TRUE,
  split1_proportion = 0.5,
  total_proportion = 1,
  split_indices = NULL,
  prior_indicator = NULL,
  prior_genesymbols = NULL,
  prior_baseline = 1e-06,
  prior_weight = 0.5,
  min_module_size = 0L,
  allocate_per_obs = TRUE,
  noise_threshold = 0.025,
  n_cycles = 50L,
  use_kmeanspp_init = TRUE,
  n_initializations = 50L,
  max_optim_iter = 10000L,
  tol_coop_rel = 1e-08,
  tol_coop_abs = 1e-12,
  tol_nnls = 1e-04,
  compute_predictive_r2 = TRUE,
  compute_silhouette = FALSE,
  nowarnings = FALSE,
  verbose = TRUE,
  quick_mode = FALSE,
  quick_mode_percent = 0.1
)
```

## Arguments

- expression:

  `p x n` matrix of pre-processed single cell expression data with `p`
  rows of genes and `n` columns of cells.

- genesymbols:

  A vector of gene names corresponding to rows of `expression`. Has to
  be of length `p`.

- is_regulator:

  An indicator vector where `1` indicates that the corresponding row in
  `expression` is a candidate regulator. All other rows represent target
  genes. Has to be of length `p`.

- penalization:

  Sparsity penalty related to the amount of regulators associated with
  each module. Either a single positive number or a vector of positive
  numbers.

- n_modules:

  Requested number of modules (integer). If this is provided without
  specifying `initial_target_modules`, then an initial module allocation
  is performed on the cross-correlation matrix of targets and genes on
  the first dataset after data splitting.

- initial_target_modules:

  The initial assignment of target genes to modules of length
  `sum(is_regulator == 0L)`. If this is not specified, then see
  `n_modules` regarding module initialization. If provided,
  `use_kmeanspp_init` and `n_initializations` are ignored.

- sample_assignment:

  A vector of sample assignment for each cell, can be used to perform
  the data splitting with stratification. Has to be of length `n`. No
  stratification if `NULL` is supplied.

- center:

  Whether or not genes should be centered within each subgroup defined
  in `sample_assignment`.

- split1_proportion:

  The proportion to use for the first dataset during data splitting. The
  proportion for the second dataset is `1 - split1_proportion`. If
  stratification with `sample_assignment` is used, then the proportion
  of each strata is controlled.

- total_proportion:

  Can be used to only use a proportion of the supplied observations. The
  proportion of the first dataset during data splitting in relation to
  the full dataset will be `total_proportion * split1_proportion`.

- split_indices:

  Can be used to provide an explicit data split. If this is supplied
  then `split1_proportion`, and `total_proportion` are ignored. Note
  that if `sample_assigment` is provided and `center == TRUE`, then
  subgroup centering will be performed as in the case of random
  splitting. A vector of length `n` containing entries 1 for cells in
  the first data split, 2 for cells in the second data split and `NA`
  for cells that should be excluded from the computations.

- prior_indicator:

  An indicator matrix (sparse or dense) of size `q x q` that indicates
  whether there is a known functional relationship between two genes.
  Ideally, this is supplied as a sparse matrix (`sparseMatrix` in the
  `Matrix` package). If not, then the matrix is converted to one.

- prior_genesymbols:

  A vector of gene names of length q corresponding to the rows/columns
  in `prior_indicator`. Does not have to be the same as `genesymbols`,
  but only useful if there is overlap.

- prior_baseline:

  A positive baseline for the network prior. The larger this parameter
  is, the less impact the network prior will have.

- prior_weight:

  A number between 0 and 1 indicating the strength of the prior in
  relation to the data. 0 ignores the prior and makes the algorithm
  completely data-driven. 1 uses only the prior during module
  allocation.

- min_module_size:

  Minimum required size of target genes in a module. Smaller modules are
  emptied.

- allocate_per_obs:

  Whether module allocation should be performed for each observation in
  the second data split separately. If `FALSE`, target genes are
  allocated into modules on the aggregate sum of squares across all
  observations in the second data split.

- noise_threshold:

  Threshold for the best \\R^2\\ of a target gene before it gets
  identified as noise.

- n_cycles:

  Number of maximum algorithmic cycles.

- use_kmeanspp_init:

  Use kmeans++ for module initialization if `initial_target_modules` is
  a single integer; otherwise use kmeans with random initial cluster
  centers

- n_initializations:

  Number of kmeans(++) initialization runs.

- max_optim_iter:

  Maximum number of iterations during optimization in the coop-Lasso and
  NNLS steps.

- tol_coop_rel:

  Relative convergence tolerance during optimization in the coop-Lasso
  step.

- tol_coop_abs:

  Absolute convergence tolerance during optimization in the coop-Lasso
  step.

- tol_nnls:

  Convergence tolerance during optimization in the NNLS step.

- compute_predictive_r2:

  Whether to compute predictive \\R^2\\ per module as well as regulator
  importance.

- compute_silhouette:

  Whether to compute silhouette scores for each target gene.

- nowarnings:

  When turned on then no warning messages are shown.

- verbose:

  Whether to print progress.

- quick_mode:

  Whether to use a reduced number of noise targets to speed up
  computations.

- quick_mode_percent:

  A number in \[0, 1) indicating the amount of noise targets to use in
  the re-allocation process if `quick_mode = TRUE`.

## Value

A list with S3 class `scregclust` containing

- penalization:

  The supplied `penalization` parameters

- results:

  A list of result lists (each with S3 class `scregclust_result`), one
  for each supplied `penalization` parameter. See below.

- initial_target_modules:

  Initial allocation of target genes into modules.

- split_indices:

  either verbatim the vector given as input or a vector encoding the
  splits as NA = not included, 1 = split 1 or 2 = split 2. Allows
  reproducibility of data splits.

For each supplied penalization parameter, `results` contains a list with

- the current `penalization` parameter,

- the supplied `genesymbols` after filtering (as used during fitting),

- the supplied `is_regulator` vector after filtering (as used during
  fitting),

- the number of fitted modules `n_modules`,

- whether the current run `converged` to a single configuration (as a
  boolean),

- as well as an `output` object containing the numeric results for each
  final configuration.

It is possible that the algorithm ends in a finite cycle of
configurations instead of a unique final configuration. Therefore,
`output` is a list with each element itself being a list with the
following contents:

- `reg_table`:

  a regulator table, a matrix of weights for each regulator and module

- `module`:

  vector of same length as `genesymbols` containing the module
  assignments for all genes with regulators marked as `NA`. Genes
  considered noise are marked as `-1`.

- `module_all`:

  same as `module`, however, genes that were marked as noise (-1 in
  `module`) are assigned to the module in which it has the largest
  \\R^2\\, even if it is below `noise_threshold`.

- `r2`:

  matrix of predictive \\R^2\\ value for each target gene and module

- `best_r2`:

  vector of best predictive \\R^2\\ for each gene (regulators marked
  with NA)

- `best_r2_idx`:

  module index corresponding to best predictive \\R^2\\ for each gene
  (regulators marked with NA)

- `r2_module`:

  a vector of predictive \\R^2\\ values for each module (included if
  `compute_predictive_r2 == TRUE`)

- `importance`:

  a matrix of importance values for each regulator (rows) and module
  (columns) (included if `compute_predictive_r2 == TRUE`)

- `r2_cross_module_per_target`:

  a matrix of cross module \\R^2\\ values for each target gene (rows)
  and each module (columns) (included if `compute_silhouette == TRUE`)

- `silhouette`:

  a vector of silhouette scores for each target gene (included if
  `compute_silhouette == TRUE`)

- `models`:

  regulator selection for each module as a matrix with regulators in
  rows and modules in columns

- `signs`:

  regulator signs for each module as a matrix with regulators in rows
  and modules in columns

- `weights`:

  average regulator coefficient for each module

- `coeffs`:

  list of regulator coefficient matrices for each module for all target
  genes as re-estimated in the NNLS step

- `sigmas`:

  matrix of residual variances, one per target gene in each module;
  derived from the residuals in NNLS step
