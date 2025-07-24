# scregclust 0.2.1

- Added a vignette illustrating how to supply your own regulator list

# scregclust 0.2.0-1

- Regulator importance was computed unnecessarily if there is only one
  regulator. A single regulator is always the most important regulator
  for a cluster.

# scregclust 0.2.0

## New features

- Quick Mode: Instead of trying to re-allocate all target genes that were
  allocated into the noise cluster, only a certain (random) percentage of
  these target genes is attempted to be re-allocated.
  
  `quick_mode = TRUE` has to be supplied as an argument to `scregclust` to
  activate this feature (off by default) and the percentage of
  noise target genes to re-allocate is given by `quick_mode_percent`,
  a number in [0, 1).

## Minor changes

- Added CRAN install instructions to the README

# scregclust 0.1.0

- First release on CRAN
