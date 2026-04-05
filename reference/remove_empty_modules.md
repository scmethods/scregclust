# Remove empty modules

Remove empty modules

## Usage

``` r
remove_empty_modules(module)
```

## Arguments

- module:

  Vector of module indices

## Value

The updated vector of module indices with empty modules removed.

## Details

Only iterates through modules with positive index, leaving the noise
module untouched.
