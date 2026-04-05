# Format count table nicely

Format count table nicely

## Usage

``` r
count_table(counts, title, row_names, col_width = 5)
```

## Arguments

- counts:

  a list of count vectors with `1 + n_cl` entries each. `NA` values are
  replaced with `-`

- title:

  title above the table

- row_names:

  a vector of row names, one for each count vector

- col_width:

  minimum width for columns

## Value

A string formatted as a table
