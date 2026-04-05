# Return list of regulator genes

Return list of regulator genes

## Usage

``` r
get_regulator_list(mode = c("TF", "kinase"))
```

## Arguments

- mode:

  Determines which genes are considered to be regulators. Currently
  supports TF=transcription factors and kinases.

## Value

a list of gene symbols

## See also

[`scregclust_format()`](https://scmethods.github.io/scregclust/reference/scregclust_format.md)
