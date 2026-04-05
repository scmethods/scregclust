# Package data before clustering

Package data before clustering

## Usage

``` r
scregclust_format(expression_matrix, mode = c("TF", "kinase"))
```

## Arguments

- expression_matrix:

  The p x n gene expression matrix with gene symbols as rownames.

- mode:

  Determines which genes are considered to be regulators.

## Value

A list with

- genesymbols:

  The gene symbols extracted from the expression matrix

- sample_assignment:

  A vector filled with `1`'s of the same length as there are columns in
  the gene expression matrix.

- is_regulator:

  Whether a gene is considered to be a regulator or not, determined
  dependent on `mode`.

## See also

[`get_regulator_list()`](https://scmethods.github.io/scregclust/reference/get_regulator_list.md)
