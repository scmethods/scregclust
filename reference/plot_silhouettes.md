# Plot individual silhouette scores

Plot individual silhouette scores

## Usage

``` r
plot_silhouettes(list_of_fits, penalization, final_config = 1L)
```

## Arguments

- list_of_fits:

  A list of `scregclust` objects each fit to the same dataset across a
  variety of module counts (varying `n_modules` when running
  [`scregclust`](https://scmethods.github.io/scregclust/reference/scregclust.md)).

- penalization:

  Either a single numeric value requesting the results for the same
  penalty parameter across all fits in `list_of_fits`, or one for each
  individual fit.

- final_config:

  The final configuration that should be visualized. Either a single
  number to be used for all fits in `list_of_fits`, or one for each
  individual fit.

## Value

A ggplot2 plot showing the the silhouette scores for each supplied fit.
