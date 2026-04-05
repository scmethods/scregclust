# Quick'n'dirty progress bar

Creates a progress bar and returns it as a string.

## Usage

``` r
progstr(step, n_steps, name, finished = FALSE, progress_length = 20L)
```

## Arguments

- step:

  current step being worked on

- n_steps:

  total number of steps

- name:

  name of the process

- finished:

  whether the process is finished

- progress_length:

  length of the progress bar in ascii signs

## Value

A string formatted as a progress bar
