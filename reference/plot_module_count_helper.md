# Plot average silhouette scores and average predictive \\R^2\\

Plot average silhouette scores and average predictive \\R^2\\

## Usage

``` r
plot_module_count_helper(list_of_fits, penalization)
```

## Arguments

- list_of_fits:

  A list of `scregclust` objects each fit to the same dataset across a
  variety of module counts (varying `n_modules` while running
  [`scregclust`](https://scmethods.github.io/scregclust/reference/scregclust.md)).

- penalization:

  Either a single numeric value requesting the results for the same
  penalty parameter across all fits in `list_of_fits`, or one for each
  individual fit.

## Value

A ggplot2 plot showing the average silhouette score and the average
predictive \\R^2\\
