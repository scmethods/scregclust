# Package index

## Setting up and performing clustering

Functions to prepare the input data and to perform single-cell
regulatory-driven clustering.

- [`scregclust()`](https://scmethods.github.io/scregclust/reference/scregclust.md)
  : Uncover gene modules and their regulatory programs from single-cell
  data
- [`scregclust_format()`](https://scmethods.github.io/scregclust/reference/scregclust_format.md)
  : Package data before clustering

## Plotting and evaluation

Functions which help in plotting and evaluating results.

- [`plot_module_count_helper()`](https://scmethods.github.io/scregclust/reference/plot_module_count_helper.md)
  : Plot average silhouette scores and average predictive \\R^2\\
- [`plot_regulator_network()`](https://scmethods.github.io/scregclust/reference/plot_regulator_network.md)
  : Plotting the regulatory table from scregclust as a directed graph
- [`plot_silhouettes()`](https://scmethods.github.io/scregclust/reference/plot_silhouettes.md)
  : Plot individual silhouette scores

## Utility functions

Functions that make accessing aspects of the results easier.

- [`get_avg_num_regulators()`](https://scmethods.github.io/scregclust/reference/get_avg_num_regulators.md)
  : Get the average number of active regulators per module
- [`get_num_final_configs()`](https://scmethods.github.io/scregclust/reference/get_num_final_configs.md)
  : Return the number of final configurations
- [`get_rand_indices()`](https://scmethods.github.io/scregclust/reference/get_rand_indices.md)
  : Compute Rand indices
- [`get_regulator_list()`](https://scmethods.github.io/scregclust/reference/get_regulator_list.md)
  : Return list of regulator genes
- [`get_target_gene_modules()`](https://scmethods.github.io/scregclust/reference/get_target_gene_modules.md)
  : Extract target gene modules for given penalization parameters

## Other helpers

- [`available_results()`](https://scmethods.github.io/scregclust/reference/available_results.md)
  : Extract final configurations into a data frame
- [`cluster_overlap()`](https://scmethods.github.io/scregclust/reference/cluster_overlap.md)
  : Create a table of module overlap for two clusterings
- [`fast_cor()`](https://scmethods.github.io/scregclust/reference/fast_cor.md)
  : Fast computation of correlation
- [`find_module_sizes()`](https://scmethods.github.io/scregclust/reference/find_module_sizes.md)
  : Determine module sizes
- [`kmeanspp()`](https://scmethods.github.io/scregclust/reference/kmeanspp.md)
  : Perform the k-means++ algorithm
