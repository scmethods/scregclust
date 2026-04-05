# Split Sample

Splits sample in train and test set

## Usage

``` r
split_sample(
  z,
  stratification,
  is_regulator,
  split_indices,
  split1_proportion,
  total_proportion,
  center
)
```

## Arguments

- z:

  matrix of single cell data with rows as genes and columns as cells.

- stratification:

  a vector by which the sampling will be stratified of length `ncol(z)`

- is_regulator:

  an indicator vector, telling which rows in `z` are candidate
  regulators

- split_indices:

  a vector of given split indices. can be `NULL`

- split1_proportion:

  proportion to include in first data split

- total_proportion:

  proportion of data to include overall in splitting

- center:

  TRUE if data should be row-centered. Set to FALSE otherwise.

## Value

a list containing

- z1_reg:

  first data split, TF-part

- z2_reg:

  second data split, TF-part

- z1_target:

  first data split, non-TF part

- z2_target:

  second data split, non-TF part

- split_indices:

  either verbatim the vector given as input or a vector encoding the
  splits as NA = not included, 1 = split 1 or 2 = split 2. Allows
  reproducibility of data splits.
