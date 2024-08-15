#' Uncover gene modules and their regulatory programs from single-cell data
#'
#' Use the scRegClust algorithm to determine gene modules and their
#' regulatory programs from single-cell data.
#'
#' @param expression `p x n` matrix of pre-processed single cell expression
#'                   data with `p` rows of genes and `n` columns of cells.
#' @param genesymbols A vector of gene names corresponding to rows of
#'                    `expression`. Has to be of length `p`.
#' @param is_regulator An indicator vector where `1` indicates that the
#'                     corresponding row in `expression` is a candidate
#'                     regulator. All other rows represent target genes.
#'                     Has to be of length `p`.
#' @param penalization Sparsity penalty related to the amount of regulators
#'                     associated with each module. Either a single positive
#'                     number or a vector of positive numbers.
#' @param n_modules Requested number of modules (integer).
#'             If this is provided without specifying `initial_target_modules`,
#'             then an initial module allocation is performed on the
#'             cross-correlation matrix of targets and genes on the first
#'             dataset after data splitting.
#' @param initial_target_modules The initial assignment of target genes to
#'                 modules of length `sum(is_regulator == 0L)`.
#'                 If this is not specified, then see `n_modules` regarding
#'                 module initialization. If provided, `use_kmeanspp_init`
#'                 and `n_initializations` are ignored.
#' @param sample_assignment A vector of sample assignment for each cell, can
#'                          be used to perform the data splitting with
#'                          stratification. Has to be of length `n`.
#'                          No stratification if `NULL` is supplied.
#' @param center Whether or not genes should be centered within each subgroup
#'               defined in `sample_assignment`.
#' @param split1_proportion The proportion to use for the first dataset during
#'                         data splitting. The proportion for the second
#'                         dataset is `1 - split1_proportion`. If stratification
#'                         with `sample_assignment` is used, then the proportion
#'                         of each strata is controlled.
#' @param total_proportion Can be used to only use a proportion of the supplied
#'                         observations. The proportion of the first dataset
#'                         during data splitting in relation to the full
#'                         dataset will be
#'                         `total_proportion * split1_proportion`.
#' @param split_indices Can be used to provide an explicit data split. If this
#'                      is supplied then `split1_proportion`, and
#'                      `total_proportion` are ignored.
#'                      Note that if `sample_assigment` is provided and
#'                      `center == TRUE`, then subgroup centering will be
#'                      performed as in the case of random splitting.
#'                      A vector of length `n` containing entries 1 for cells
#'                      in the first data split, 2 for cells in the second
#'                      data split and `NA` for cells that should be excluded
#'                      from the computations.
#' @param prior_indicator An indicator matrix (sparse or dense) of size `q x q`
#'                        that indicates whether there is a known functional
#'                        relationship between two genes. Ideally, this is
#'                        supplied as a sparse matrix (`sparseMatrix`
#'                        in the `Matrix` package). If not, then the matrix
#'                        is converted to one.
#' @param prior_genesymbols A vector of gene names of length q corresponding
#'                          to the rows/columns in `prior_indicator`. Does not
#'                          have to be the same as `genesymbols`, but only
#'                          useful if there is overlap.
#' @param prior_baseline A positive baseline for the network prior. The larger
#'                       this parameter is, the less impact the network prior
#'                       will have.
#' @param prior_weight A number between 0 and 1 indicating the strength of the
#'                     prior in relation to the data. 0 ignores the prior and
#'                     makes the algorithm completely data-driven. 1 uses only
#'                     the prior during module allocation.
#' @param min_module_size Minimum required size of target genes in a module.
#'                        Smaller modules are emptied.
#' @param allocate_per_obs Whether module allocation should be performed for
#'                         each observation in the second data split separately.
#'                         If `FALSE`, target genes are allocated into modules
#'                         on the aggregate sum of squares across all
#'                         observations in the second data split.
#' @param noise_threshold Threshold for the best \eqn{R^2} of a target gene
#'                        before it gets identified as noise.
#' @param n_cycles Number of maximum algorithmic cycles.
#' @param use_kmeanspp_init Use kmeans++ for module initialization if
#'                          `initial_target_modules` is a single integer;
#'                          otherwise use kmeans with random initial cluster
#'                          centers
#' @param n_initializations Number of kmeans(++) initialization runs.
#' @param max_optim_iter Maximum number of iterations during optimization
#'                       in the coop-Lasso and NNLS steps.
#' @param tol_coop_rel Relative convergence tolerance during optimization
#'                     in the coop-Lasso step.
#' @param tol_coop_abs Absolute convergence tolerance during optimization
#'                     in the coop-Lasso step.
#' @param tol_nnls Convergence tolerance during optimization in the NNLS step.
#' @param compute_predictive_r2 Whether to compute predictive \eqn{R^2} per
#'                              module as well as regulator importance.
#' @param compute_silhouette Whether to compute silhouette scores for each
#'                           target gene.
#' @param nowarnings When turned on then no warning messages are shown.
#' @param verbose Whether to print progress.
#'
#' @return A list with S3 class `scregclust` containing
#'   \item{penalization}{The supplied `penalization` parameters}
#'   \item{results}{A list of result lists (each with S3 class
#'                  `scregclust_result`), one for each supplied `penalization`
#'                  parameter. See below.}
#'   \item{initial_target_modules}{Initial allocation of target genes into
#'                                 modules.}
#'   \item{split_indices}{either verbatim the vector given as input or
#'                        a vector encoding the splits as NA = not included,
#'                        1 = split 1 or 2 = split 2. Allows reproduciblity
#'                        of data splits.}
#'
#' For each supplied penalization parameter, `results` contains a list with
#' * the current `penalization` parameter,
#' * the supplied `genesymbols` after filtering (as used during fitting),
#' * the supplied `is_regulator` vector after filtering (as used during
#'   fitting),
#' * the number of fitted modules `n_modules`,
#' * whether the current run `converged` to a single configuration (as a
#'   boolean),
#' * as well as an `output` object containing the numeric results for each
#'   final configuration.
#'
#' It is possible that the algorithm ends in a finite cycle of configurations
#' instead of a unique final configuration.
#' Therefore, `output` is a list with each element itself being a list
#' with the following contents:
#' \describe{
#'   \item{`reg_table`}{a regulator table, a matrix of weights for each
#'                      regulator and module}
#'   \item{`module`}{vector of same length as `genesymbols` containing the
#'                   module assignments for all genes with regulators
#'                   marked as `NA`. Genes considered noise are marked as `-1`.}
#'   \item{`module_all`}{same as `module`, however, genes that were marked as
#'                       noise (-1 in `module`) are assigned to the
#'                       module in which it has the largest \eqn{R^2},
#'                       even if it is below `noise_threshold`.}
#'   \item{`r2`}{matrix of predictive \eqn{R^2} value for each target gene and
#'              module}
#'   \item{`best_r2`}{vector of best predictive \eqn{R^2} for each gene
#'                   (regulators marked with NA)}
#'   \item{`best_r2_idx`}{module index corresponding to best predictive
#'                        \eqn{R^2} for each gene (regulators marked with NA)}
#'   \item{`r2_module`}{a vector of predictive \eqn{R^2} values for each
#'                      module (included if `compute_predictive_r2 == TRUE`)}
#'   \item{`importance`}{a matrix of importance values for each regulator (rows)
#'                       and module (columns) (included if
#'                       `compute_predictive_r2 == TRUE`)}
#'   \item{`r2_cross_module_per_target`}{a matrix of cross module \eqn{R^2}
#'                                        values for each target gene (rows)
#'                                        and each module (columns) (included
#'                                        if `compute_silhouette == TRUE`)}
#'   \item{`silhouette`}{a vector of silhouette scores for each target gene
#'                       (included if `compute_silhouette == TRUE`)}
#'   \item{`models`}{regulator selection for each module as a matrix with
#'                   regulators in rows and modules in columns}
#'   \item{`signs`}{regulator signs for each module as a matrix with
#'                  regulators in rows and modules in columns}
#'   \item{`weights`}{average regulator coefficient for each module}
#'   \item{`coeffs`}{list of regulator coefficient matrices for each module
#'                   for all target genes as re-estimated in the NNLS step}
#'   \item{`sigmas`}{matrix of residual variances, one per target gene
#'                   in each module; derived from the residuals in NNLS step}
#' }
#'
#' @concept main
#'
#' @export
scregclust <- function(expression,
                       genesymbols,
                       is_regulator,
                       penalization,
                       n_modules,
                       initial_target_modules = NULL,
                       sample_assignment = NULL,
                       center = TRUE,
                       split1_proportion = 0.5,
                       total_proportion = 1,
                       split_indices = NULL,
                       prior_indicator = NULL,
                       prior_genesymbols = NULL,
                       prior_baseline = 1e-6,
                       prior_weight = 0.5,
                       min_module_size = 0L,
                       allocate_per_obs = TRUE,
                       noise_threshold = 0.025,
                       n_cycles = 50L,
                       use_kmeanspp_init = TRUE,
                       n_initializations = 50L,
                       max_optim_iter = 10000L,
                       tol_coop_rel = 1e-8,
                       tol_coop_abs = 1e-12,
                       tol_nnls = 1e-4,
                       compute_predictive_r2 = TRUE,
                       compute_silhouette = FALSE,
                       nowarnings = FALSE,
                       verbose = TRUE) {
  ###############################
  # START input validation
  ###############################

  if (verbose) {
    cat(paste0(cli::symbol$arrow_right, " Validating input"))
    cl <- TRUE
    start_time <- Sys.time()
  }

  if (!(is.matrix(expression) && is.numeric(expression))) {
    if (verbose && cl) {
      cat("\n")
      cl <- FALSE
    }
    cli::cli_abort("{.var expression} needs to be a numeric matrix.")
  }
  expression <- as.matrix(expression)
  p <- nrow(expression)
  n <- ncol(expression)

  if (length(genesymbols) != p) {
    if (verbose && cl) {
      cat("\n")
      cl <- FALSE
    }
    cli::cli_abort(c(
      "{.var genesymbols} needs to be of length {p}.",
      "x" = "Supplied vector has length {length(genesymbols)}.",
      "i" =
        "There needs to be one gene symbol for each gene in {.var expression}."
    ))
  } else {
    genesymbols <- genesymbols
  }

  if (!(
    length(is_regulator) == p
    && (
      all(unique(is_regulator) %in% c(0, 1))
      || all(unique(is_regulator) %in% c(TRUE, FALSE))
    )
  )) {
    if (verbose && cl) {
      cat("\n")
      cl <- FALSE
    }
    cli::cli_abort(c(
      paste(
        "{.var is_regulator} needs to be 0/1 or TRUE/FALSE vector",
        "of length {p}."
      ),
      "i" = paste(
        "{.var is_regulator} indicates for each gene in",
        "{.var expression} whether it is a regulator (1) or not (0)."
      )
    ))
  } else {
    is_regulator <- as.logical(is_regulator)
  }

  if (sum(is_regulator) == 0) {
    if (verbose && cl) {
      cat("\n")
      cl <- FALSE
    }
    cli::cli_abort(c(
      paste(
        "At least one regulator needs to be present in the dataset."
      ),
      "i" = paste(
        "{.var is_regulator} is not equal to TRUE/1 for any gene.",
      )
    ))
  }

  if (!(
    is.numeric(penalization)
    && length(penalization) >= 1L
    && all(penalization > 0)
  )) {
    if (verbose && cl) {
      cat("\n")
      cl <- FALSE
    }
    cli::cli_abort(c(
      "{.var penalization} is not supplied correctly.",
      "x" = paste(
        "{.var penalization} needs to be a numeric vector of positive numbers",
        "of minimum length 1."
      )
    ))
  } else {
    penalization <- as.double(penalization)
  }

  if (!(
    is.numeric(n_modules)
    && length(n_modules) == 1L
    && as.integer(n_modules) == n_modules
    && n_modules >= 1
    && n_modules <= sum(is_regulator == 0)
  )) {
    if (verbose && cl) {
      cat("\n")
      cl <- FALSE
    }
    cli::cli_abort(c(
      "{.var n_modules} is not supplied correctly.",
      "x" = paste(
        "An integer between 1 and the number of target genes",
        "({sum(is_regulator == 0)}) specifying the maximum number of modules."
      )
    ))
  }

  if (!is.null(initial_target_modules)) {
    if (!(
      is.numeric(initial_target_modules)
      && all(as.integer(initial_target_modules) == initial_target_modules)
      && length(initial_target_modules) == sum(is_regulator == 0)
      && all(initial_target_modules >= 1)
      && all(initial_target_modules <= n_modules)
    )) {
      if (verbose && cl) {
        cat("\n")
        cl <- FALSE
      }
      cli::cli_abort(c(
        "{.var initial_target_modules} is not supplied correctly.",
        "x" = paste(
          "A vector of initial module indices for the target genes of",
          "length {sum(is_regulator == 0)}."
        ),
        "x" = paste(
          "Entries need to be between 1 and {.var n_modules} = {n_modules}."
        )
      ))
    } else {
      initial_target_modules <- as.integer(initial_target_modules)
    }
  }

  if (!(
    is.null(sample_assignment)
    || (
      is.vector(sample_assignment)
      && length(sample_assignment) == n
    )
  )) {
    if (verbose && cl) {
      cat("\n")
      cl <- FALSE
    }
    cli::cli_abort(c(
      paste(
        "{.var sample_assignment} needs to be a vector of length {n}",
        "or {.code NULL}"
      ),
      "i" = paste(
        "If not {.code NULL}, then a sample assignment needs to be supplied",
        "for each column of {.var expression}."
      )
    ))
  } else {
    sample_assignment <- sample_assignment
  }

  if (!(
    is.logical(center)
    && length(center) == 1
  )) {
    if (verbose && cl) {
      cat("\n")
      cl <- FALSE
    }
    cli::cli_abort(
      "{.var center} needs to be TRUE or FALSE."
    )
  } else {
    center <- center
  }

  # If no split indices are provided, then the data split is determined
  # randomly dependent on split1_proportion and total_proportion.
  # Otherwise, the provided indices are used.
  if (is.null(split_indices)) {
    if (!(
      is.numeric(split1_proportion)
      && length(split1_proportion) == 1
      && 0 < split1_proportion
      && split1_proportion < 1
    )) {
      if (verbose && cl) {
        cat("\n")
        cl <- FALSE
      }
      cli::cli_abort(
        "{.var split1_proportion} needs to be between 0 and 1."
      )
    } else {
      split1_proportion <- as.double(split1_proportion)
    }

    if (!(
      is.numeric(total_proportion)
      && length(total_proportion) == 1
      && 0 < total_proportion
      && total_proportion <= 1
    )) {
      if (verbose && cl) {
        cat("\n")
        cl <- FALSE
      }
      cli::cli_abort(paste(
        "{.var total_proportion} needs to be between 0 (exclusive)",
        "and 1 (inclusive)."
      ))
    } else {
      total_proportion <- as.double(total_proportion)
    }

    n_samples_split1 <- floor(total_proportion * n * split1_proportion)
    n_samples_split2 <- floor(total_proportion * n * (1 - split1_proportion))

    if (
      n_samples_split1 <= sum(is_regulator == 1)
      || n_samples_split2 <= sum(is_regulator == 1)
    ) {
      if (verbose && cl) {
        cat("\n")
        cl <- FALSE
      }
      if (!nowarnings) {
        cli::cli_inform(c(
          "Few cells compared with regulators.",
          "*" = paste(
            "For optimal performance, each data split should contain more",
            "cells than there are regulators."
          ),
          "*" = "Consider reducing the number of regulators.",
          "i" = "Number of regulators = {sum(is_regulator == 1)}",
          "i" = "Cells in split 1 = {n_samples_split1}",
          "i" = "Cells in split 2 = {n_samples_split2}"
        ))
      }
    }
  } else {
    if (!(
      is.numeric(split_indices)
      && length(split_indices) == n
      && all(split_indices %in% c(NA, 1, 2))
      && sum(split_indices == 1, na.rm = TRUE) >= 1
      && sum(split_indices == 2, na.rm = TRUE) >= 1
    )) {
      if (verbose && cl) {
        cat("\n")
        cl <- FALSE
      }
      cli::cli_abort(c(
        "Wrong format for {.var split_indices}.",
        "i" = "{.code length(split_indices)} = {length(split_indices)}",
        "i" = "Should be {n}",
        "i" = "Should only contain entries 1, 2, or {.code NA}.",
        "i" = "Needs to contain at least one 1 and one 2."
      ))
    } else {
      split_indices <- as.integer(split_indices)
    }

    if (any(table(split_indices) <= sum(is_regulator == 1))) {
      if (verbose && cl) {
        cat("\n")
        cl <- FALSE
      }
      if (!nowarnings) {
        cli::cli_warn(c(
          "Few cells compared with regulators.",
          "*" = paste(
            "For optimal performance, each data split should contain more",
            "cells than there are regulators."
          ),
          "*" = "Consider reducing the number of regulators.",
          "i" = "Number of regulators = {sum(is_regulator == 1)}",
          "i" = "Cells in split 1 = {sum(split_indices == 1, na.rm = TRUE)}",
          "i" = "Cells in split 2 = {sum(split_indices == 2, na.rm = TRUE)}"
        ))
      }
    }
  }

  if (!is.null(prior_indicator)) {
    if (!(
      (
        methods::is(prior_indicator, "Matrix")
        || is.matrix(prior_indicator)
      )
      && length(unique(dim(prior_indicator))) == 1L
      && Matrix::isSymmetric(prior_indicator)
    )) {
      if (verbose && cl) {
        cat("\n")
        cl <- FALSE
      }
      cli::cli_abort(c(
        paste(
          "{.var prior_indicator} needs to be a square symmetric",
          "indicator matrix."
        ),
        "i" = paste(
          "Can be supplied as a normal R matrix or as a matrix created",
          "with the {.pkg Matrix} package (e.g. a sparse matrix)."
        )
      ))
    }

    if (is.null(prior_genesymbols)) {
      if (verbose && cl) {
        cat("\n")
        cl <- FALSE
      }
      cli::cli_abort(
        paste(
          "{.var prior_genesymbols} missing, but {.var prior_indicator}",
          "is available."
        )
      )
    }

    if (length(prior_genesymbols) != nrow(prior_indicator)) {
      if (verbose && cl) {
        cat("\n")
        cl <- FALSE
      }
      cli::cli_abort(c(
        paste(
          "{.var prior_genesymbols} needs to have one element for each",
          "row/column of {.var prior_indicator}"
        ),
        "x" = "Length {.var prior_genesymbols}: {length(prior_genesymbols)}",
        "x" = "Dimension of {.var prior_indicator}: {nrow(prior_indicator)}"
      ))
    }

    if (!(
      is.numeric(prior_baseline)
      && length(prior_baseline) == 1L
      && prior_baseline > 0
    )) {
      if (verbose && cl) {
        cat("\n")
        cl <- FALSE
      }
      cli::cli_abort(
        "{.var prior_baseline} needs to be a positive number."
      )
    } else {
      prior_baseline <- as.double(prior_baseline)
    }

    if (!(
      is.numeric(prior_weight)
      && length(prior_weight) == 1L
      && prior_weight >= 0
      && prior_weight <= 1
    )) {
      if (verbose && cl) {
        cat("\n")
        cl <- FALSE
      }
      cli::cli_abort(
        "{.var prior_weight} needs to be a number between 0 and 1."
      )
    } else {
      prior_weight <- as.double(prior_weight)
    }

    # Remove diagonal elements and save as list of indices
    Matrix::diag(prior_indicator) <- 0
    nonzero <- Matrix::which(prior_indicator != 0, arr.ind = TRUE)
    o1 <- order(nonzero[, 1])
    nz1 <- nonzero[, 1][o1]
    nz2 <- nonzero[, 2][o1]
    nonzero <- NULL

    idx_rle <- rle(nz1)
    idx <- c(0L, cumsum(idx_rle$lengths))
    prior_indicator_list <- vector("list", nrow(prior_indicator))
    m <- 1L
    len <- length(idx_rle$values)
    for (l in seq_len(nrow(prior_indicator))) {
      if (m > len) {
        break
      }

      if (idx_rle$values[m] == l) {
        prior_indicator_list[[l]] <- nz2[(idx[m] + 1):idx[m + 1L]]
        m <- m + 1L
      } else {
        prior_indicator_list[[l]] <- vector("integer", 0)
      }
    }
  } else {
    # Necessary for technical reasons
    prior_indicator_list <- vector("list", 0)
    prior_genesymbols <- vector("character", 0)
    prior_baseline <- 1e-6
    prior_weight <- 0
  }

  if (!(
    length(min_module_size) == 1
    && as.integer(min_module_size) == min_module_size
    && min_module_size >= 0
  )) {
    if (verbose && cl) {
      cat("\n")
      cl <- FALSE
    }
    cli::cli_abort(
      "{.var min_module_size} needs to be a non-negative integer."
    )
  } else {
    min_module_size <- as.integer(min_module_size)
  }

  if (!(
    is.logical(allocate_per_obs)
    && length(allocate_per_obs) == 1
  )) {
    if (verbose && cl) {
      cat("\n")
      cl <- FALSE
    }
    cli::cli_abort(
      "{.var allocate_per_obs} needs to be TRUE or FALSE."
    )
  } else {
    allocate_per_obs <- allocate_per_obs
  }

  if (!(
    is.numeric(noise_threshold)
    && length(noise_threshold) == 1
    && noise_threshold >= 0
  )) {
    if (verbose && cl) {
      cat("\n")
      cl <- FALSE
    }
    cli::cli_abort(
      "{.var noise_threshold} needs to be a non-negative number."
    )
  } else {
    noise_threshold <- as.double(noise_threshold)
  }

  if (!(
    length(n_cycles) == 1
    && as.integer(n_cycles) == n_cycles
    && n_cycles >= 0
  )) {
    if (verbose && cl) {
      cat("\n")
      cl <- FALSE
    }
    cli::cli_abort(
      "{.var n_cycles} needs to be a non-negative integer."
    )
  } else {
    n_cycles <- as.integer(n_cycles)
  }

  if (is.null(initial_target_modules)) {
    if (!(
      is.logical(use_kmeanspp_init)
      && length(use_kmeanspp_init) == 1
    )) {
      if (verbose && cl) {
        cat("\n")
        cl <- FALSE
      }
      cli::cli_abort(
        "{.var use_kmeanspp_init} needs to be TRUE or FALSE."
      )
    } else {
      use_kmeanspp_init <- use_kmeanspp_init
    }

    if (!(
      length(n_initializations) == 1
      && as.integer(n_initializations) == n_initializations
      && n_initializations > 0
    )) {
      if (verbose && cl) {
        cat("\n")
        cl <- FALSE
      }
      cli::cli_abort(
        "{.var n_initializations} needs to be a positive integer."
      )
    } else {
      n_initializations <- as.integer(n_initializations)
    }
  }

  if (!(
    length(max_optim_iter) == 1
    && as.integer(max_optim_iter) == max_optim_iter
    && max_optim_iter > 0
  )) {
    if (verbose && cl) {
      cat("\n")
      cl <- FALSE
    }
    cli::cli_abort(
      "{.var max_optim_iter} needs to be a positive integer."
    )
  } else {
    max_optim_iter <- as.integer(max_optim_iter)
  }

  if (!(
    length(tol_coop_rel) == 1
    && tol_coop_rel > 0
  )) {
    if (verbose && cl) {
      cat("\n")
      cl <- FALSE
    }
    cli::cli_abort(
      "{.var tol_coop_rel} needs to be a positive scalar."
    )
  } else {
    tol_coop_rel <- as.double(tol_coop_rel)
  }

  if (!(
    length(tol_coop_abs) == 1
    && tol_coop_abs > 0
  )) {
    if (verbose && cl) {
      cat("\n")
      cl <- FALSE
    }
    cli::cli_abort(
      "{.var tol_coop_abs} needs to be a positive scalar."
    )
  } else {
    tol_coop_abs <- as.double(tol_coop_abs)
  }

  if (!(
    length(tol_nnls) == 1
    && tol_nnls > 0
  )) {
    if (verbose && cl) {
      cat("\n")
      cl <- FALSE
    }
    cli::cli_abort(
      "{.var tol_nnls} needs to be a positive scalar."
    )
  } else {
    tol_nnls <- as.double(tol_nnls)
  }

  if (!(
    is.logical(compute_predictive_r2)
    && length(compute_predictive_r2) == 1
  )) {
    if (verbose && cl) {
      cat("\n")
      cl <- FALSE
    }
    cli::cli_abort(
      "{.var compute_predictive_r2} needs to be TRUE or FALSE."
    )
  } else {
    compute_predictive_r2 <- compute_predictive_r2
  }

  if (!(
    is.logical(compute_silhouette)
    && length(compute_silhouette) == 1
  )) {
    if (verbose && cl) {
      cat("\n")
      cl <- FALSE
    }
    cli::cli_abort(
      "{.var compute_silhouette} needs to be TRUE or FALSE."
    )
  } else {
    compute_silhouette <- compute_silhouette
  }

  if (!(
    is.logical(verbose)
    && length(verbose) == 1
  )) {
    if (verbose && cl) {
      cat("\n")
      cl <- FALSE
    }
    cli::cli_abort(
      "{.var verbose} needs to be TRUE or FALSE."
    )
  } else {
    verbose <- verbose
  }

  if (verbose) {
    end_time <- Sys.time()
    if (cl) {
      # ANSI code for clearing the line and reset cursor to the beginning
      # of the line
      cat("\33[2K\r")
    }
    cli::cli_alert_success(paste0(
      "Input validated ",
      cli::col_blue(
        "(in ", prettyunits::pretty_dt(end_time - start_time), ")"
      )
    ))
  }

  ###############################
  # END input validation
  ###############################



  ###############################
  # START sample splitting
  ###############################

  if (verbose) {
    cat(paste0(cli::symbol$arrow_right, " Split samples"))
    cl <- TRUE
    start_time <- Sys.time()
  }

  if (is.null(sample_assignment)) {
    stratification <- rep.int(1L, n)
  } else {
    # Create stratification indices from 1..length(unique(sample_assignment))
    assignment_vals <- unique(sample_assignment)

    stratification <- rep.int(0L, n)
    for (i in seq_along(assignment_vals)) {
      stratification[sample_assignment == assignment_vals[i]] <- i
    }
  }

  split_z <- split_sample(
    expression,
    stratification,
    is_regulator,
    split_indices,
    split1_proportion,
    total_proportion,
    center
  )

  z1_reg <- split_z$z1_reg
  z2_reg <- split_z$z2_reg
  z1_target <- split_z$z1_target
  z2_target <- split_z$z2_target
  split_indices <- split_z$split_indices

  # Remove unnecessary variable to save memory
  split_z <- NULL

  # Check for constant genes after data splitting -- remove them if necessary
  constant_reg <- unique(c(
    which(apply(z1_reg, 2, sd) == 0),
    which(apply(z2_reg, 2, sd) == 0)
  ))
  constant_target <- unique(c(
    which(apply(z1_target, 2, sd) == 0),
    which(apply(z2_target, 2, sd) == 0)
  ))

  gs_constant_reg <- genesymbols[is_regulator == 1][constant_reg]
  gs_constant_target <- genesymbols[is_regulator == 0][constant_target]

  if (length(c(gs_constant_reg, gs_constant_target)) > 0) {
    if (verbose) {
      cl <- FALSE
      cat("\n")
      cli::cli_alert_info(
        paste(
          "The following genes are constant in at least one data split and",
          "will be removed."
        )
      )
      cli::cli_ul()
      if (length(gs_constant_reg) > 0) {
        cli::cli_li("Regulators: {gs_constant_reg}")
      }
      if (length(gs_constant_target) > 0) {
        cli::cli_li("Target genes: {gs_constant_target}")
      }
      cli::cli_end()
    }

    gs_remove <- (
      genesymbols %in% c(gs_constant_reg, gs_constant_target)
    )

    gs_constant_reg <- genesymbols[is_regulator == 1] %in% gs_constant_reg
    gs_constant_target <- (
      genesymbols[is_regulator == 0] %in% gs_constant_target
    )

    genesymbols <- genesymbols[!gs_remove]
    is_regulator <- is_regulator[!gs_remove]
    p <- length(genesymbols)

    z1_reg <- z1_reg[, !gs_constant_reg, drop = FALSE]
    z2_reg <- z2_reg[, !gs_constant_reg, drop = FALSE]
    z1_target <- z1_target[, !gs_constant_target, drop = FALSE]
    z2_target <- z2_target[, !gs_constant_target, drop = FALSE]

    if (!is.null(initial_target_modules)) {
      initial_target_modules <- initial_target_modules[!gs_constant_target]
    }

    if (verbose) {
      cli::cli_alert_info(paste(
        "{sum(gs_constant_reg)} regulator{?s} and",
        "{sum(gs_constant_target)} target gene{?s} removed."
      ))
      cli::cli_alert_info(paste(
        "{sum(is_regulator)} regulator{?s} and {sum(!is_regulator)} target",
        "gene{?s} remaining."
      ))
    }
  }

  if (verbose) {
    end_time <- Sys.time()
    if (cl) {
      # ANSI code for clearing the line and reset cursor to the beginning
      # of the line
      cat("\33[2K\r")
    }
    cli::cli_alert_success(paste0(
      "Samples split ",
      cli::col_blue(
        "(in ", prettyunits::pretty_dt(end_time - start_time), ")"
      )
    ))
  }

  ###############################
  # END sample splitting
  ###############################



  ###############################
  # START pre-computing
  ###############################

  if (verbose) {
    sb <- cli::cli_status("{cli::symbol$arrow_right} Precomputing")
    start_time <- Sys.time()
  }

  # Scale first data split for regression preprocessing
  # No scaling of the target genes. Their scale is implicitly dealt with below.
  z1_target_centered <- scale(z1_target, scale = FALSE)
  z1_target_sds <- apply(z1_target, 2, sd)
  z1_reg_scaled <- scale(z1_reg)

  # Scale second data split as first one since it will be used as a
  # test set in Step 2
  z2_target_centered <- scale(
    z2_target, colMeans(z1_target), scale = FALSE
  )
  z2_reg_scaled <- scale(
    z2_reg, colMeans(z1_reg), apply(z1_reg, 2, sd)
  )

  n_target <- ncol(z1_target_centered) # number of targets
  n_reg <- ncol(z1_reg_scaled) # number of predictors
  n1 <- nrow(z1_target_centered) # cells in first split
  n2 <- nrow(z2_target_centered) # cells in second split

  # Remove unnecessary variables to save memory
  z1_reg <- NULL
  z2_reg <- NULL
  z1_target <- NULL
  z2_target <- NULL

  # Extract genesymbols for TFs and non-TFs
  genesymbols_reg <- genesymbols[is_regulator == 1L]
  genesymbols_target <- genesymbols[is_regulator == 0L]

  # Extract subset of prior_indicator that is actually relevant for
  # target genesymbols
  genes_overlap <- intersect(genesymbols_target, prior_genesymbols)
  prior_idx_prior <- as.vector(stats::na.omit(
    match(genes_overlap, prior_genesymbols)
  ))
  prior_idx_target <- as.vector(stats::na.omit(
    match(genes_overlap, genesymbols_target)
  ))

  prior_indicator_tmp <- lapply(
    seq_along(genesymbols_target),
    function(i) vector("integer", 0)
  )
  prior_indicator_tmp[prior_idx_target] <- lapply(
    prior_indicator_list[prior_idx_prior],
    function(idx) {
      which(genesymbols_target %in% prior_genesymbols[idx]) - 1
    }
  )

  prior_indicator_list <- prior_indicator_tmp

  if (n1 >= n_reg) {
    # Default case: Enough cells to perform OLS

    # Pre-compute target gene standard deviation
    beta_init <- coef_ols(z1_target_centered, z1_reg_scaled)
    init_df <- n_reg
  } else {
    # High-dimensional case. More regulators than cells in the first data split
    # Use a residual variance estimator which is stable in p > n scenario
    #
    # Equations 15 and below in
    #
    # L. H. Dicker. Variance estimation in high-dimensional linear models.
    # Biometrika, 101(2):269–284, 2014.
    #
    # Seem like a good choice but cannot figure out how to make them
    # always positive. Use ridge regression for now.

    # m1 <- sum(z1_reg_scaled^2) / (n1 * n_reg)
    # ztz <- crossprod(z1_reg_scaled)
    # m2 <- (
    #   sum((ztz / n1)^2) / n_reg
    #   - sum(z1_reg_scaled^2)^2 / (n1 * n_reg * n1^2)
    # )

    # # Estimate residual standard deviation for each response
    # # Complex correlation case
    # z1_target_centered_res_sds <- sqrt(
    #   (1 + n_reg * m1^2 / ((n1 + 1) * m2)) / n1 * colSums(z1_target_centered^2)
    #   - m1 / (n1 * (n1 + 1) * m2)
    #   * colSums(crossprod(z1_reg_scaled, z1_target_centered)^2)
    # )

    # # Orthogonal predictors case
    # z1_target_centered_res_sds <- sqrt(
    #   (n_reg + n1 + 1) / (n1 * (n1 + 1)) * colSums(z1_target_centered^2)
    #   - 1 / (n1 * (n1 + 1))
    #   * colSums(crossprod(z1_reg_scaled, z1_target_centered)^2)
    # )

    ztz <- crossprod(z1_reg_scaled)
    es <- eigen(ztz, symmetric = TRUE, only.values = TRUE)$values

    lambda_minimal <- es[n1 - 1] + 1e-4 # Minimal eigenvalue + fudge factor
    # Effective degrees of freedom
    init_df <- sum(es / (es + lambda_minimal))
    while (init_df >= n1) {
      lambda_minimal <- 2 * lambda_minimal
      init_df <- sum(es / (es + lambda_minimal))
    }

    beta_init <- coef_ridge(z1_target_centered, z1_reg_scaled, lambda_minimal)
  }

  # Estimate residual standard deviation for each response
  z1_target_centered_res_sds <- sqrt(
    colSums((z1_target_centered - z1_reg_scaled %*% beta_init)^2)
    / (n1 - init_df)
  )

  beta_adapt_sq <- (
    beta_init %*% diag(
      1 / z1_target_centered_res_sds, nrow = n_target, ncol = n_target
    )
  )^2

  # Pre-compute cross-correlation matrices
  cross_corr1 <- fast_cor(z1_target_centered, z1_reg_scaled)
  # cross_corr2 <- fast_cor(z2_target_centered, z2_reg_scaled)

  # Initial module allocation
  if (is.null(initial_target_modules)) {
    # Find tentative modules with k-means (use k-means++ to choose
    # good initial clusterings)
    if (use_kmeanspp_init) {
      k_start_cl <- kmeanspp(
        cross_corr1,
        n_cluster = n_modules,
        n_init_clusterings = n_initializations,
        n_max_iter = 100L
      )
    } else {
      k_start_cl <- stats::kmeans(
        cross_corr1,
        n_modules,
        iter.max = 100L,
        nstart = n_initializations
      )$cluster
    }
  } else {
    k_start_cl <- initial_target_modules
  }

  if (verbose) {
    end_time <- Sys.time()
    cli::cli_alert_success(paste0(
      "Precomputing done ",
      cli::col_blue(
        "(in ", prettyunits::pretty_dt(end_time - start_time), ")"
      )
    ))
  }

  ###############################
  # END pre-computing
  ###############################

  ###############################
  # START additional setup
  ###############################

  if (verbose) {
    cli::cli_status_update(
      id = sb,
      "{cli::symbol$arrow_right} Additional setup"
    )
    start_time <- Sys.time()
  }

  if (allocate_per_obs) {
    # Pre-allocate large SSQ array to speed up computations below
    sq_residuals_test <- alloc_array(z2_target_centered^2, n_modules)
    attr(sq_residuals_test, "dim") <- c(
      n_modules, n2, n_target
    )
  } else {
    sq_residuals_test <- NULL
  }

  if (verbose) {
    end_time <- Sys.time()
    cli::cli_alert_success(paste0(
      "Additional setup done ",
      cli::col_blue(
        "(in ", prettyunits::pretty_dt(end_time - start_time), ")"
      )
    ))
  }

  ###############################
  # END additional setup
  ###############################

  if (verbose) {
    cli::cli_status_clear(id = sb)
    cat("\n")
    cli::cli_text(
      cli::col_green(cli::symbol$square_small_filled),
      " All ready for clustering."
    )

    target_counts <- c(
      NA_integer_, sapply(seq_len(n_modules), function(s) sum(k_start_cl == s))
    )
    cat("\n")
    cat(count_table(
      list(target_counts),
      title = "Initial counts",
      row_names = "Targets"
    ))
  }

  ###############################
  # START algorithm
  ###############################

  results <- vector("list", length(penalization))
  for (l in seq_along(penalization)) {
    if (verbose) {
      cat("\n")
      cli::cli_text(
        cli::col_yellow("{.strong #}"),
        " {.strong Clustering} with penalization = ",
        cli::col_blue("{penalization[l]}")
      )
    }

    # On the last cycle we only re-estimate the regulator coefficients in Step 1
    last_cycle <- FALSE

    # Track whether the algorithm converged (or ended up in a cycle)
    converged <- FALSE

    k <- k_start_cl

    # Module history to determine convergence
    n_k <- 2
    k_history <- list()
    k_history[[1]] <- k

    best_r2 <- list()

    coeffs_final <- list()
    sigmas_final <- list()
    r2_final <- list()
    best_r2_final <- list()
    best_r2_idx_final <- list()
    models_final <- list()
    signs_final <- list()
    weights_final <- list()

    final_cycle_length <- 1L

    ############################################################################
    ## Start cycles ############################################################
    ############################################################################

    for (cycle in seq_len(n_cycles + 1L)) {
      # Check if we are on the last cycle
      if (cycle == n_cycles + 1L) {
        if (n_cycles > 0L) {
          # If we arrive here then no convergence occured. Inform the
          # user about this.
          cli::cli_alert_info(paste(
            "Reached maximum number of iterations without convergence."
          ))
          cli::cli_text(
            cli::col_yellow(cli::symbol$square_small_filled),
            " Stopping"
          )
        }
        last_cycle <- TRUE
      }

      if (verbose) {
        cat("\n")
        if (last_cycle) {
          cli::cli_text(
            cli::symbol$pointer,
            cli::col_blue(" {.strong Last cycle}"),
            " No module allocation"
          )
        } else {
          cli::cli_text(
            cli::symbol$pointer,
            cli::col_blue(" {.strong Cycle} "),
            "{cycle} / {n_cycles}"
          )
        }
      }

      ###############################################
      ## START Step 1: Finding the best regulators ##
      ###############################################

      for (m in seq_len(final_cycle_length)) {
        # If the final cycle is longer than 1, then repeat the estimation
        # of coefficients and therefore models, weights, signs for each of
        # the final cycle configurations.
        k <- k_history[[length(k_history) - m + 1L]]

        models <- matrix(FALSE, nrow = n_reg, ncol = n_modules)
        weights <- matrix(0, nrow = n_reg, ncol = n_modules)
        signs <- matrix(NA_real_, nrow = n_reg, ncol = n_modules)

        if (last_cycle) {
          coeffs_final[[m]] <- vector("list", n_modules)
        }

        if (verbose) {
          start_time <- Sys.time()
          msg <- "{cli::symbol$arrow_right} {.strong Step 1:} Select regulators"
          if (last_cycle && final_cycle_length > 1) {
            msg <- paste(
              msg, cli::col_yellow(sprintf("[final configuration #%d]", m))
            )
          }
          sb <- cli::cli_status(msg)
        }

        for (j in seq_len(n_modules)) {
          if (verbose) {
            module_progstr <- progstr(j, n_modules, "modules")
            cli::cli_status_update(
              id = sb,
              paste(msg, module_progstr)
            )
          }

          z1_target_centered_cl <- z1_target_centered[, k == j, drop = FALSE]

          n_target_cl <- ncol(z1_target_centered_cl) # number of genes in module

          if (n_target_cl > 0) {
            # Compute weights based on initial (OLS or ridge) coefficient
            # estimates and residual standard deviation
            ws <- (
              rep.int(sqrt(n_target_cl), n_reg)
              / sqrt(sqrt(rowSums(beta_adapt_sq[, k == j, drop = FALSE])))
            )

            admm_fit <- coop_lasso(
              (
                (
                  z1_target_centered_cl
                  %*% diag(
                    1 / z1_target_centered_res_sds[k == j],
                    nrow = n_target_cl,
                    ncol = n_target_cl
                  )
                ) / sqrt(nrow(z1_target_centered_cl))
              ),
              z1_reg_scaled / sqrt(nrow(z1_reg_scaled)),
              penalization[l], ws,
              eps_rel = tol_coop_rel,
              eps_abs = tol_coop_abs,
              max_iter = max_optim_iter,
              verbose = FALSE
            )

            beta <- (
              admm_fit$beta %*% diag(
                z1_target_centered_res_sds[k == j],
                nrow = n_target_cl,
                ncol = n_target_cl
              )
            )
            # ADMM only soft-thresholds the other variable (zeta), not beta.
            # At convergence, both should be very close but ensure that
            # almost-zeros are actual zeros.
            beta[abs(beta) < 1e-6] <- 0

            models[, j] <- rowSums(abs(beta) > 0) > 0
            weights[, j] <- rowMeans(beta)
            signs[models[, j], j] <- sign(weights[models[, j], j])
          }
        }

        if (last_cycle) {
          models_final[[m]] <- models
          signs_final[[m]] <- signs
          weights_final[[m]] <- weights
        }

        if (verbose) {
          end_time <- Sys.time()
          cli::cli_status_clear(id = sb)
          msg <- "{.strong Step 1:} Regulators selected"
          if (last_cycle && final_cycle_length > 1) {
            msg <- paste0(
              msg, " ", cli::col_yellow(sprintf("[final configuration #%d]", m))
            )
          }
          cli::cli_alert_success(paste0(
            msg,
            " ",
            cli::col_blue(
              "(in ", prettyunits::pretty_dt(end_time - start_time), ")"
            )
          ))
          if (last_cycle) {
            target_counts <- sapply(
              c(-1, seq_len(n_modules)), function(s) sum(k == s)
            )
            reg_counts <- c(NA_integer_, colSums(models))
            tbl_title <- if (final_cycle_length > 1) {
              sprintf("Final counts for configuration #%d", m)
            } else {
              "Final counts"
            }
            cat(count_table(
              list(target_counts, reg_counts),
              title = tbl_title,
              row_names = c("Targets", "Regulators")
            ))
            cat("\n")
          }
        }


        ########################################################
        ## START Step 2: Determining target gene coefficients ##
        ########################################################

        # If there are no predictors in a model, then the best prediction
        # is constant zero, since there is no intercept. The sum-of-squares
        # for each gene is then just the squared response.
        sum_squares_train <- matrix(
          rep.int(colSums(z1_target_centered^2), n_modules),
          nrow = n_target,
          ncol = n_modules
        )

        sum_squares_test <- matrix(
          rep.int(colSums(z2_target_centered^2), n_modules),
          nrow = n_target,
          ncol = n_modules
        )

        if (allocate_per_obs && !last_cycle) {
          # Reset values in SSQ array
          reset_array(sq_residuals_test, z2_target_centered^2)
        }

        if (verbose) {
          start_time <- Sys.time()
          sb <- cli::cli_status(paste(
            "{cli::symbol$arrow_right} {.strong Step 2a:}",
            "Re-estimate coefficients"
          ))
        }

        for (j in seq_len(n_modules)) {
          if (verbose) {
            module_progstr <- progstr(j, n_modules, "modules")
            cli::cli_status_update(
              id = sb,
              paste(
                "{cli::symbol$arrow_right} {.strong Step 2a:}",
                "Re-estimate coefficients",
                module_progstr
              )
            )
          }

          # Get regulators used in model for module j
          reg_cl <- which(models[, j] == TRUE)
          if (length(reg_cl) > 0L) {
            if (length(reg_cl) > n1 && !nowarnings) {
              cli::cli_inform(c(
                paste(
                  "More selected regulators in module {j} than cells.",
                  "Results may be instable."
                ),
                "i" = paste(
                  "Consider reducing the number of regulators or increasing",
                  "the penalty parameter."
                ),
                "i" = "Number of regulators = {length(reg_cl)}",
                "i" = "Number of cells = {n1}"
              ))
            }

            z1_reg_scaled_cl <- z1_reg_scaled[, reg_cl, drop = FALSE]
            z2_reg_scaled_cl <- z2_reg_scaled[, reg_cl, drop = FALSE]

            signs_cl <- signs[reg_cl, j]
            # Adjust the sign of the predicting regulators and estimate
            # coefficients using non-negative least squares for
            # sign-consistent estimates (following Meinshausen, 2013)
            z1_reg_scaled_cl_sign_corrected <- (
              z1_reg_scaled_cl
              %*% diag(
                signs_cl,
                nrow = length(signs_cl),
                ncol = length(signs_cl)
              )
            )

            beta_hat_nnls <- coef_nnls(
              z1_reg_scaled_cl_sign_corrected,
              z1_target_centered %*% diag(
                1 / z1_target_sds,
                nrow = n_target,
                ncol = n_target
              ),
              eps = tol_nnls, max_iter = max_optim_iter
            )$beta * signs_cl * z1_target_sds

            residuals_train_nnls <- (
              z1_target_centered - z1_reg_scaled_cl %*% beta_hat_nnls
            )
            residuals_test_nnls <- (
              z2_target_centered - z2_reg_scaled_cl %*% beta_hat_nnls
            )

            # NNLS squared residuals, SSQ and R-square
            if (allocate_per_obs && !last_cycle) {
              sq_residuals_test[j, , ] <- residuals_test_nnls^2
            }
            sum_squares_test[, j] <- colSums(residuals_test_nnls^2)
            sum_squares_train[, j] <- colSums(residuals_train_nnls^2)

            if (last_cycle) {
              coeffs_final[[m]][[j]] <- beta_hat_nnls
            }
          }
        }

        if (verbose) {
          end_time <- Sys.time()
          cli::cli_status_clear(id = sb)
          cli::cli_alert_success(paste(
            "{.strong Step 2a:} Coefficients re-estimated",
            cli::col_blue(
              "(in ", prettyunits::pretty_dt(end_time - start_time), ")"
            )
          ))
        }

        r2_test <- 1 - (
          sum_squares_test / colSums(z2_target_centered^2)
        )
        best_r2[[cycle]] <- apply(r2_test, 1, max)

        resid_var <- t(
          t(sum_squares_train)
          / (n1 - colSums(models))
        )

        if (last_cycle) {
          sigmas_final[[m]] <- sqrt(resid_var)
          r2_final[[m]] <- r2_test
          best_r2_final[[m]] <- best_r2[[cycle]]
          best_r2_idx_final[[m]] <- apply(r2_test, 1, which.max)
        }
      }

      # Do not perform allocation on last cycle
      if (last_cycle) {
        break
      }

      if (verbose) {
        start_time <- Sys.time()
        sb <- cli::cli_status(paste(
          "{cli::symbol$arrow_right} {.strong Step 2b:}",
          "Allocating modules"
        ))
      }

      if (allocate_per_obs) {
        update_order <- sample(length(k), length(k)) - 1L
        idx_new <- allocate_into_modules(
          sq_residuals_test,
          resid_var,
          prior_indicator_list,
          k,
          update_order,
          prior_baseline,
          prior_weight
        )

        idx <- idx_new
      } else {
        # We could choose modules on aggregated results instead of
        # per observation vote
        idx <- k
        module_totals <- sapply(seq_len(n_modules), function(i) sum(idx == i))

        # Model log-likelihood
        model_log_likelihood <- (
          n2 * log(2 * pi * resid_var)
          + sum_squares_test / resid_var
        )

        # Convert to probabilities across modules
        max_model_log_likelihood <- apply(model_log_likelihood, 1, max)
        model_log_likelihood_m_max <- (
          model_log_likelihood - max_model_log_likelihood
        )
        model_log_prob <- (
          model_log_likelihood_m_max
          - log(rowSums(exp(model_log_likelihood_m_max)))
        )

        # Update one gene at a time to ensure correctness of the prior
        # After changing module assignment of a single gene, the prior needs
        # to be recalculated, which is what we are doing here.
        update_order <- sample(length(idx), length(idx))
        for (gene in update_order) {
          prior_overlap_modules <- sapply(
            # indices are 0-based for C++ code
            prior_indicator_list[[gene]] + 1,
            function(i) idx[i]
          )
          prior_frac <- sapply(
            seq_len(n_modules),
            function(i) sum(prior_overlap_modules == i)
          )

          prior_frac[module_totals > 0] <- (
            prior_frac[module_totals > 0] / module_totals[module_totals > 0]
          )
          prior_frac <- prior_frac + prior_baseline
          prior_log_prob <- log(prior_frac) - log(sum(prior_frac))

          total_model_log_scores <- (
            (1 - prior_weight) * model_log_prob[gene, ]
            + prior_weight * prior_log_prob
          )

          idx[gene] <- which.min(total_model_log_scores)
        }
      }

      # Put all genes hard to predict in noise module, marked by -1
      idx[which(best_r2[[cycle]] < noise_threshold)] <- -1L

      # Enforce minimum module size by emptying small modules.
      # Observations in those modules get assigned to the noise module
      # but can be re-assigned in the next cycle.
      module_size <- sapply(seq_len(n_modules), function(i) sum(idx == i))
      idx[idx %in% which(module_size < min_module_size)] <- -1

      k <- as.integer(idx)

      if (verbose) {
        end_time <- Sys.time()
        cli::cli_status_clear(id = sb)
        cli::cli_alert_success(paste(
          "{.strong Step 2b:} Modules allocated",
          cli::col_blue(
            "(in ", prettyunits::pretty_dt(end_time - start_time), ")"
          )
        ))

        target_counts <- sapply(
          c(-1, seq_len(n_modules)), function(s) sum(k == s)
        )
        reg_counts <- c(NA_integer_, colSums(models))

        cat(count_table(
          list(target_counts, reg_counts),
          title = "Current counts",
          row_names = c("Targets", "Regulators")
        ))
      }

      matching_clustering <- which(
        sapply(k_history, function(k_) all(k_ == k))
      )

      k_history[[n_k]] <- k

      if (length(matching_clustering) > 0L) {
        if (verbose) {
          cat("\n")
        }
        converged <- TRUE
        if (n_k - matching_clustering == 1L) {
          if (verbose) {
            cli::cli_alert_info(
              "No change since last cycle."
            )
            cli::cli_text(
              cli::col_green(cli::symbol$square_small_filled),
              " Clustering converged"
            )
          }
        } else {
          if (verbose) {
            cli::cli_alert_info(paste(
              "Repetitive cycle of length {n_k - matching_clustering}",
              "discovered."
            ))
            cli::cli_text(
              cli::col_yellow(cli::symbol$square_small_filled),
              " Stopping"
            )
          }

          final_cycle_length <- n_k - matching_clustering
        }

        last_cycle <- TRUE
      }

      n_k <- n_k + 1L
    }

    if (verbose) {
      start_time <- Sys.time()
      sb <- cli::cli_status(
        "{cli::symbol$arrow_right} Post-processing and packaging"
      )
    }

    module <- vector("list", final_cycle_length)
    module_all <- vector("list", final_cycle_length)
    r2 <- vector("list", final_cycle_length)
    r2_idx <- vector("list", final_cycle_length)
    reg_table <- vector("list", final_cycle_length)
    r2_removed <- vector("list", final_cycle_length)
    r2_module <- vector("list", final_cycle_length)
    importance <- vector("list", final_cycle_length)
    r2_cross_module_per_target <- vector("list", final_cycle_length)
    silhouette <- vector("list", final_cycle_length)

    for (m in seq_len(final_cycle_length)) {
      k <- k_history[[length(k_history) - m + 1L]]

      if (compute_predictive_r2) {
        r2_removed[[m]] <- vector("list", n_modules)
        r2_module[[m]] <- rep(NA_real_, n_modules)
        importance[[m]] <- matrix(
          NA_real_, nrow = sum(is_regulator), ncol = n_modules
        )

        for (j in seq_len(n_modules)) {
          # Get regulators used in module j
          reg_cl <- which(models_final[[m]][, j] == TRUE)
          # Get genes in module j
          target_cl <- which(k == j)

          if (length(target_cl) > 0L && length(reg_cl) > 0L) {
            ssq <- matrix(
              NA_real_, nrow = length(target_cl), ncol = length(reg_cl) + 1
            )
            remove_reg <- c(list(integer(0)), lapply(reg_cl, function(r) r))

            for (r in seq_along(remove_reg)) {
              reg_cl_mod <- setdiff(reg_cl, remove_reg[r])

              z1_reg_scaled_cl <- z1_reg_scaled[, reg_cl_mod, drop = FALSE]
              z2_reg_scaled_cl <- z2_reg_scaled[, reg_cl_mod, drop = FALSE]

              signs_cl <- signs_final[[m]][reg_cl_mod, j]
              # Adjust the sign of the predicting regulators and estimate
              # coefficients using non-negative least squares for
              # sign-consistent estimates (following Meinshausen, 2013)
              z1_reg_scaled_cl_sign_corrected <- (
                z1_reg_scaled_cl
                %*% diag(
                  signs_cl,
                  nrow = length(signs_cl),
                  ncol = length(signs_cl)
                )
              )

              beta_hat_nnls <- coef_nnls(
                z1_reg_scaled_cl_sign_corrected,
                z1_target_centered[, target_cl, drop = FALSE] %*% diag(
                  1 / z1_target_sds[target_cl],
                  nrow = length(target_cl),
                  ncol = length(target_cl)
                ),
                eps = tol_nnls, max_iter = max_optim_iter
              )$beta * signs_cl * z1_target_sds[target_cl]

              ssq[, r] <- colSums((
                z2_target_centered[, target_cl, drop = FALSE]
                - z2_reg_scaled_cl %*% beta_hat_nnls
              )^2)
            }

            # Compute R2
            r2_removed_vec <- 1 - (
              colSums(ssq) / sum(scale(
                z2_target_centered[, target_cl, drop = FALSE],
                center = TRUE,
                scale = FALSE
              )^2)
            )

            r2_removed[[m]][[j]] <- 1 - t(
              t(ssq) / colSums(scale(
                z2_target_centered[, target_cl, drop = FALSE],
                center = TRUE,
                scale = FALSE
              )^2)
            )
            colnames(r2_removed[[m]][[j]]) <- c("None", genesymbols_reg[reg_cl])
            rownames(r2_removed[[m]][[j]]) <- genesymbols_target[target_cl]

            # Save module specific R2
            r2_module[[m]][j] <- r2_removed_vec[1]

            # Compute importances
            importance[[m]][reg_cl, j] <- (
              1 - (
                r2_removed_vec[2L:(length(reg_cl) + 1L)]
                / r2_module[[m]][j]
              )
            )
          }
        }
      }

      if (compute_silhouette) {
        # # Compute cross-module R2
        # sum_squares_test <- matrix(
        #   rep.int(colSums(z2_target_centered^2), n_modules),
        #   nrow = n_target,
        #   ncol = n_modules
        # )

        # for (j in seq_len(n_modules)) {
        #   # Get regulators used in model for module j
        #   reg_cl <- which(models_final[[m]][, j] == TRUE)
        #   if (length(reg_cl) > 0L) {
        #     # z1_reg_scaled_cl <- z1_reg_scaled[, reg_cl, drop = FALSE]
        #     # z2_reg_scaled_cl <- z2_reg_scaled[, reg_cl, drop = FALSE]

        #     # signs_cl <- signs[reg_cl, j]
        #     # # Adjust the sign of the predicting regulators and estimate
        #     # # coefficients using non-negative least squares for
        #     # # sign-consistent estimates (following Meinshausen, 2013)
        #     # z1_reg_scaled_cl_sign_corrected <- (
        #     #   z1_reg_scaled_cl
        #     #   %*% diag(
        #     #     signs_cl,
        #     #     nrow = length(signs_cl),
        #     #     ncol = length(signs_cl)
        #     #   )
        #     # )

        #     # beta_hat_nnls <- coef_nnls(
        #     #   z1_reg_scaled_cl_sign_corrected,
        #     #   z1_target_centered %*% diag(
        #     #     1 / z1_target_sds,
        #     #     nrow = n_target,
        #     #     ncol = n_target
        #     #   ),
        #     #   eps = tol_nnls, max_iter = max_optim_iter
        #     # )$beta * signs_cl * z1_target_sds

        #     sum_squares_test[, j] <- colSums((
        #       z2_target_centered - (
        #         z2_reg_scaled[, reg_cl, drop = FALSE] %*% coeffs_final[[m]][[j]]
        #       # z2_reg_scaled_cl %*% beta_hat_nnls
        #       )
        #     )^2)

        #   }
        # }

        r2_cross_module_per_target[[m]] <- r2_final[[m]]
        # (
        #   1 - sum_squares_test / colSums(z2_target_centered^2)
        # )

        silhouette[[m]] <- sapply(seq_along(k), function(i) {
          c <- k[i]
          if (c != -1) {
            b <- max(r2_cross_module_per_target[[m]][i, ][-c])
            if (b < 0) {
              b <- 0
            }
            a <- r2_cross_module_per_target[[m]][i, c]
            (a - b) / max(a, b)
          } else {
            NA
          }
        })
      }

      # # "Hack" to put regulators in tentative module
      # k_ <- k
      # k_[k == -1L] <- n_modules + 1L
      # non_empty_modules <- which(sapply(
      #   seq_len(n_modules + 1L), function(j) sum(k_ == j) > 0
      # ))
      # module_indicator <- Matrix::sparseMatrix(i = seq_along(k_), j = k_)

      # idx <- non_empty_modules[apply(
      #   diag(
      #     1 / Matrix::colSums(module_indicator)[non_empty_modules],
      #     nrow = length(non_empty_modules)
      #   )
      #   %*% t(module_indicator[, non_empty_modules])
      #   %*% cross_corr2,
      #   2,
      #   which.max
      # )]
      # idx[idx == n_modules + 1L] <- -1L

      module[[m]] <- rep.int(
        NA_integer_, n_reg + n_target
      )
      # module[[m]][is_regulator == 1L] <- idx
      module[[m]][is_regulator == 0L] <- k

      # Put rag bag genes into closest module dependent on R2 value.
      # Leave all other genes in their allocated modules.
      module_all[[m]] <- rep.int(
        NA_integer_, n_reg + n_target
      )
      module_all[[m]][is_regulator == 0L] <- k
      module_all[[m]][is_regulator == 0L][k == -1L] <- (
        best_r2_idx_final[[m]][k == -1L]
      )

      r2[[m]] <- rep.int(NA_real_, n_reg + n_target)
      r2[[m]][is_regulator == 0] <- best_r2_final[[m]]

      r2_idx[[m]] <- rep.int(
        NA_integer_, n_reg + n_target
      )
      r2_idx[[m]][is_regulator == 0] <- best_r2_idx_final[[m]]

      max_corr <- matrix(0, nrow = ncol(cross_corr1), ncol = n_modules)
      for (i in seq_len(n_modules)) {
        max_corr[, i] <- apply(
          cross_corr1[k == i, , drop = FALSE], 2, stats::median
        )
      }

      tag <- vector(mode = "character", n_modules)
      for (i in seq_len(n_modules)) {
        tag[i] <- paste0("module ", i)
      }

      reg_table[[m]] <- as.data.frame(
        models_final[[m]] * weights_final[[m]] * (abs(max_corr) > 0.025)
      )
      rownames(reg_table[[m]]) <- genesymbols_reg
      colnames(reg_table[[m]]) <- tag
    }

    if (verbose) {
      end_time <- Sys.time()
      cli::cli_status_clear(id = sb)
      cli::cli_alert_success(paste0(
        "Post-processing and packaging done ",
        cli::col_blue(
          "(in ", prettyunits::pretty_dt(end_time - start_time), ")"
        )
      ))
      cat("\n")
      cli::cli_alert_success(
        paste0(
          "Finished clustering with penalization = ",
          cli::col_blue("{penalization[l]}")
        )
      )
    }

    results[[l]] <- structure(list(
      penalization = penalization[l],
      genesymbols = genesymbols,
      is_regulator = is_regulator,
      n_modules = n_modules,
      converged = converged,
      output = lapply(seq_len(final_cycle_length), function(m) {
        structure(list(
          reg_table = reg_table[[m]],
          module = module[[m]],
          module_all = module_all[[m]],
          r2 = r2_final[[m]], # predictive R2 for each target gene and module
          best_r2 = r2[[m]], # best predictive R2 for each gene
          best_r2_idx = r2_idx[[m]], # module index of best r2
          r2_module = r2_module[[m]],
          importance = importance[[m]],
          r2_cross_module_per_target = r2_cross_module_per_target[[m]],
          silhouette = silhouette[[m]],
          models = models_final[[m]],
          signs = signs_final[[m]],
          weights = weights_final[[m]],
          coeffs = coeffs_final[[m]],
          sigmas = sigmas_final[[m]]
        ), class = "scregclust_output")
      })
    ), class = "scregclust_result")
  }

  ###############################
  # END main algorithm
  ###############################

  if (verbose) {
    cat("\n")
    cli::cli_alert_success("Finished!")
  }

  structure(
    list(
      penalization = penalization,
      results = results,
      initial_target_modules = k_start_cl,
      split_indices = split_indices
    ),
    class = "scregclust"
  )
}

format_scregclust_result <- function(res) {
  n_configs <- length(res$output)

  paste0(
    cli::col_grey("# scRegClust result"),
    "\n",
    cli::col_grey("# Clustering with penalization = "),
    cli::col_blue(res$penalization),
    "\n",
    do.call(
      function(...) paste0(..., collapse = "\n\n"),
      lapply(seq_along(res$output), function(i) {
        config <- res$output[[i]]
        target_counts <- sapply(c(-1, seq_len(res$n_modules)), function(s) {
          sum(config$module[!res$is_regulator] == s)
        })
        reg_counts <- c(NA_integer_, colSums(config$models))
        tbl_title <- if (n_configs > 1) {
          sprintf("Final counts for configuration #%d", i)
        } else {
          "Final counts"
        }
        count_table(
          list(target_counts, reg_counts),
          title = tbl_title,
          row_names = c("Targets", "Regulators")
        )
      })
    )
  )
}

#' @export
format.scregclust_result <- function(x, ...) {
  cli::ansi_strip(format_scregclust_result(x))
}

#' @export
print.scregclust_result <- function(x, ...) {
  cat(format_scregclust_result(x), "\n")
}

format_scregclust <- function(fit) {
  update_cs_str_widths <- function(str_widths) {
    if (length(str_widths) == 1) {
      str_widths[1]
    } else {
      cumsum(str_widths + c(0, rep(2, length(str_widths) - 1)))
    }
  }

  width <- cli::console_width()
  pen_strs <- as.character(fit$penalization)
  str_widths <- sapply(pen_strs, nchar)

  cs_str_widths <- update_cs_str_widths(str_widths)
  strs <- list(which(cs_str_widths < width))
  str_widths <- str_widths[cs_str_widths >= width]
  while (length(str_widths) > 0) {
    cs_str_widths <- update_cs_str_widths(str_widths)
    strs <- c(strs, list(
      max(strs[[length(strs)]]) + which(cs_str_widths < width)
    ))
    str_widths <- str_widths[cs_str_widths >= width]
  }

  pen_str <- paste(sapply(strs, function(idx) {
    paste(pen_strs[idx], collapse = ", ")
  }), collapse = ",\n")

  paste0(
    cli::col_grey("# scRegClust fit object"),
    "\n\n",
    cli::col_grey("# Modules: "),
    cli::col_blue(fit$results[[1]]$n_modules),
    "\n",
    cli::col_grey("# Penalization parameters:"),
    "\n",
    cli::col_grey("# "),
    cli::col_blue(pen_str)
  )
}

#' @export
format.scregclust <- function(x, ...) {
  cli::ansi_strip(format_scregclust(x))
}

#' @export
print.scregclust <- function(x, ...) {
  cat(format_scregclust(x), "\n")
}

format_scregclust_output <- function(output) {
  update_cs_str_widths <- function(str_widths) {
    if (length(str_widths) == 1) {
      str_widths[1]
    } else {
      cumsum(str_widths + c(0, rep(2, length(str_widths) - 1)))
    }
  }

  width <- cli::console_width()
  nms <- names(output)
  nms <- nms[!sapply(output, is.null)]
  str_widths <- sapply(nms, nchar)

  cs_str_widths <- update_cs_str_widths(str_widths)
  strs <- list(which(cs_str_widths < width))
  str_widths <- str_widths[cs_str_widths >= width]
  while (length(str_widths) > 0) {
    cs_str_widths <- update_cs_str_widths(str_widths)
    strs <- c(strs, list(
      max(strs[[length(strs)]]) + which(cs_str_widths < width)
    ))
    str_widths <- str_widths[cs_str_widths >= width]
  }

  nms_str <- paste(sapply(strs, function(idx) {
    paste(nms[idx], collapse = ", ")
  }), collapse = ",\n")

  paste0(
    cli::col_grey("# scRegClust output object"),
    "\n\n",
    cli::col_grey("# Contains: "),
    "\n",
    cli::col_grey("# "),
    cli::col_blue(nms_str)
  )
}

#' @export
format.scregclust_output <- function(x, ...) {
  cli::ansi_strip(format_scregclust_output(x))
}

#' @export
print.scregclust_output <- function(x, ...) {
  cat(format_scregclust_output(x), "\n")
}

#' Split Sample
#'
#' Splits sample in train and test set
#'
#' @param z matrix of single cell data with rows as genes and columns as cells.
#' @param stratification a vector by which the sampling will be stratified
#'                       of length `ncol(z)`
#' @param is_regulator an indicator vector, telling which rows in `z` are
#'                     candidate regulators
#' @param split_indices a vector of given split indices. can be `NULL`
#' @param split1_proportion proportion to include in first data split
#' @param total_proportion proportion of data to include overall in splitting
#' @param center TRUE if data should be row-centered. Set to FALSE otherwise.
#'
#' @return a list containing
#'   \item{z1_reg}{first data split, TF-part}
#'   \item{z2_reg}{second data split, TF-part}
#'   \item{z1_target}{first data split, non-TF part}
#'   \item{z2_target}{second data split, non-TF part}
#'   \item{split_indices}{either verbatim the vector given as input or
#'                        a vector encoding the splits as NA = not included,
#'                        1 = split 1 or 2 = split 2. Allows reproduciblity
#'                        of data splits.}
#'
#' @keywords internal
split_sample <- function(z, stratification, is_regulator, split_indices,
                         split1_proportion, total_proportion, center) {
  if (is.null(split_indices)) {
    is_included <- sort(do.call(
      c,
      lapply(
        unname(split(seq_len(ncol(z)), stratification)),
        function(s) sample(s, floor(length(s) * total_proportion))
      )
    ))

    is_split1 <- sort(do.call(
      c,
      lapply(
        unname(split(seq_along(is_included), stratification[is_included])),
        function(s) sample(s, floor(length(s) * split1_proportion))
      )
    ))

    split_indices <- rep.int(NA, ncol(z))
    split_indices[is_included] <- 2
    split_indices[is_included[is_split1]] <- 1
  } else {
    is_included <- which(!is.na(split_indices))
    is_split1 <- which(split_indices[is_included] == 1L)
  }

  z_ <- z[, is_included]
  if (center) {
    z_ <- do.call(cbind, lapply(
      unname(split(seq_along(is_included), stratification[is_included])),
      function(s) {
        # Use the means from split 1 to normalize data
        idx <- intersect(s, is_split1)
        z_[, s] - rowMeans(z_[, idx])
      }
    ))
  }

  list(
    z1_reg = t(z_[, is_split1][is_regulator == 1L, ]),
    z2_reg = t(z_[, -is_split1][is_regulator == 1L, ]),
    z1_target = t(z_[, is_split1][is_regulator == 0L, ]),
    z2_target = t(z_[, -is_split1][is_regulator == 0L, ]),
    split_indices = split_indices
  )
}

#' Package data before clustering
#'
#' @param expression_matrix The p x n gene expression matrix with gene symbols
#'                          as rownames.
#' @param mode Determines which genes are considered to be regulators.
#'
#' @return A list with
#'   \item{genesymbols}{The gene symbols extracted from the expression matrix}
#'   \item{sample_assignment}{A vector filled with `1`'s of the same length as
#'                            there are columns in the gene expression matrix.}
#'   \item{is_regulator}{Whether a gene is considered to be a regulator or not,
#'                       determined dependent on `mode`.}
#'
#' @seealso [get_regulator_list()]
#'
#' @concept main
#'
#' @export
scregclust_format <- function(expression_matrix, mode = c("TF", "kinase")) {
  regulators <- get_regulator_list(mode)
  genesymbols <- rownames(expression_matrix)

  is_regulator <- rep(0, nrow(expression_matrix))
  idx <- which(genesymbols %in% regulators)
  is_regulator[idx] <- 1

  sample_assignment <- rep(1, ncol(expression_matrix))

  list(
    genesymbols = genesymbols,
    sample_assignment = sample_assignment,
    is_regulator = is_regulator
  )
}

#' Return list of regulator genes
#'
#' @param mode Determines which genes are considered to be regulators.
#'             Currently supports TF=transcription factors and kinases.
#' @return a list of gene symbols
#' @seealso [scregclust_format()]
#'
#' @concept utilities
#'
#' @export
get_regulator_list <- function(mode = c("TF", "kinase")) {
  mode <- match.arg(mode)
  if (mode == "TF") {
    return(human_tfs_v3)
  } else if (mode == "kinase") {
    return(human_kinases)
  }
}
