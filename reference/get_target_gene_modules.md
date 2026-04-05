# Extract target gene modules for given penalization parameters

Extract target gene modules for given penalization parameters

## Usage

``` r
get_target_gene_modules(fit, penalization = NULL)
```

## Arguments

- fit:

  An object of class `scregclust`

- penalization:

  A numeric vector of penalization parameters. The penalization
  parameters specified here must have been used used during fitting of
  the `fit` object.

## Value

A list of lists of final target modules. One list for each parameter in
`penalization`. The lists contain the modules of target genes for each
final configuration.
