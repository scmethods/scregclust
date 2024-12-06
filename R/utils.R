#' Compute OLS coefficients
#'
#' If the design matrix has full column-rank, then use the normal
#' least squares estimate. Otherwise, use the Moore-Penrose inverse
#' to compute the least squares estimate.
#'
#' @param y Target vector (n x 1)/matrix (n x m)
#' @param x Design matrix (n x p)
#'
#' @return Vector of OLS coefficients
#'
#' @keywords internal
coef_ols <- function(y, x) {
  # Pre-compute quantities
  n <- nrow(x)
  p <- ncol(x)
  xtx <- crossprod(x)
  xty <- crossprod(x, y)

  if (n < p) {
    # Compute the pseudo-inverse of xtx
    xtx_svd <- svd(xtx)
    d <- xtx_svd$d
    idx <- which(d > .Machine$double.eps * p * max(d))
    d[idx] <- 1 / d[idx]
    d[setdiff(seq_len(p), idx)] <- 0
    xtx_inv <- xtx_svd$v %*% diag(d, nrow = p, ncol = p) %*% t(xtx_svd$u)

    # Compute least squares solution using pseudo-inverse
    beta_ols <- xtx_inv %*% xty
  } else {
    # Compute least squares solution directly
    beta_ols <- solve(xtx, xty)
  }

  beta_ols
}

#' Compute ridge regression coefficients
#'
#'
#' @param y Target vector (n x 1)/matrix (n x m)
#' @param x Design matrix (n x p)
#' @param lambda Positive parameter for ridge penalty
#'
#' @return Vector of ridge regression coefficients
#'
#' @keywords internal
coef_ridge <- function(y, x, lambda) {
  # Pre-compute quantities
  p <- ncol(x)
  xtx <- crossprod(x)
  xty <- crossprod(x, y)

  # Compute ridge regression solution directly
  solve(xtx + diag(lambda, p, p), xty)
}

#' Quick'n'dirty progress bar
#'
#' Creates a progress bar and returns it as a string.
#'
#' @param step current step being worked on
#' @param n_steps total number of steps
#' @param name name of the process
#' @param finished whether the process is finished
#' @param progress_length length of the progress bar in ascii signs
#'
#' @return A string formatted as a progress bar
#'
#' @keywords internal
progstr <- function(step, n_steps, name,
                    finished = FALSE, progress_length = 20L) {
  steps_done <- floor(progress_length * (step - 1) / n_steps)

  parts <- c("|", rep.int(cli::col_blue(cli::symbol$square), steps_done))
  if (!finished) {
    parts <- c(
      parts,
      rep.int(cli::symbol$line, progress_length - steps_done)
    )
  } else {
    parts <- c(parts, rep.int(
      cli::col_blue(cli::symbol$square), progress_length - steps_done
    ))
  }
  parts <- c(parts, "| ", cli::col_grey("%d/%d ", name))

  sprintf(paste0(parts, collapse = ""), step, n_steps)
}

#' Format count table nicely
#'
#' @param counts a list of count vectors with `1 + n_cl` entries each.
#'               `NA` values are replaced with `-`
#' @param title title above the table
#' @param row_names a vector of row names, one for each count vector
#' @param col_width minimum width for columns
#'
#' @return A string formatted as a table
#'
#' @keywords internal
count_table <- function(counts,
                        title,
                        row_names,
                        col_width = 5) {
  nms <- c("Noise", as.character(seq_len(length(counts[[1]]) - 1)))
  counts_chr <- lapply(counts, as.character)
  # Replace NA with `-`
  counts_chr <- lapply(counts_chr, function(cn) {
    cn[is.na(cn)] <- "-"
    cn
  })

  stopifnot(length(row_names) == length(counts_chr))

  cws <- c(
    max(
      c(
        nchar("Noise"),
        sapply(counts_chr, function(cn) nchar(cn[1])),
        col_width
      )
    ),
    do.call(pmax, c(
      list(unname(sapply(nms[-1], nchar))),
      lapply(counts_chr, function(cn) unname(sapply(cn[-1], nchar))),
      list(col_width)
    ))
  )

  # longest row name
  width_row_nms <- max(sapply(c("Module", row_names), nchar))

  fmt_strs <- sprintf("%%%ds", cws)
  fmt_row_str <- sprintf("  %%%ds  ", width_row_nms)

  width <- cli::console_width()

  # Two spaces + longest row name + two spaces
  tbl_widths <- (
    2 + width_row_nms + 2 + cumsum(cws + c(0, rep(3, length(cws) - 1)))
  )
  tbl_rows <- list(which(tbl_widths <= width))
  cws_tmp <- cws[tbl_widths > width]
  while (length(cws_tmp) > 0) {
    tbl_widths <- (
      2 + width_row_nms + 2
      + cumsum(cws_tmp + c(0, rep(3, length(cws_tmp) - 1)))
    )

    tbl_rows <- c(tbl_rows, list(
      max(tbl_rows[[length(tbl_rows)]]) + which(tbl_widths <= width))
    )
    cws_tmp <- cws_tmp[tbl_widths > width]
  }

  # Grey-out `-` and `0`s
  counts_chr <- lapply(counts_chr, function(cn) {
    cn_out <- sprintf(fmt_strs, cn)
    cn_out[cn == "-"] <- cli::col_grey(sprintf(fmt_strs[cn == "-"], "-"))
    cn_out[cn == "0"] <- cli::col_grey(sprintf(fmt_strs[cn == "0"], "0"))
    cn_out
  })

  do.call(function(...) paste(..., sep = "\n"), c(
    list(cli::col_grey(sprintf("# %s", title))),
    lapply(
      tbl_rows,
      function(elems) {
        paste0(paste(
          paste0(
            cli::col_blue(sprintf(fmt_row_str, "Module")),
            cli::col_grey(
              paste0(
                sprintf(fmt_strs[elems], nms[elems]),
                collapse = cli::col_grey(" | ")
              )
            )
          ),
          do.call(
            function(...) paste(..., sep = "\n"),
            lapply(seq_along(counts_chr), function(i) {
              paste0(
                cli::col_blue(sprintf(fmt_row_str, row_names[i])),
                paste0(
                  sprintf(fmt_strs[elems], counts_chr[[i]][elems]),
                  collapse = cli::col_grey(" | ")
                )
              )
            })
          ),
          sep = "\n"
        ), "\n")
      }
    )
  ))
}

#' Compute indicator matrix of pairwise distances smaller than threshold
#'
#' Computes the Jaccard distance between rows of a matrix and returns a
#' sparse symmetric indicator matrix containing the entries with a distance
#' of less than a given upper bound. Note that the diagonal is always 1.
#'
#' @param x the input matrix with vectors to be compared in the rows.
#' @param upper_bnd pairs with a Jaccard distance below this upper bound are
#'                  returned as 1 while all others receive the entry 0.
#'
#' @return A list of vectors describing a sparse lower triangular pattern matrix
#'      \item{i}{Row indices}
#'      \item{j}{Column indices}
#'
#' @keywords internal
jaccard_indicator <- function(x, upper_bnd = 0.8) {
  # Treat matrix as sparse pattern matrix
  x <- methods::as(x, "ngCMatrix")

  # Dimension along which pairwise distances are computed
  n <- x@Dim[1]

  # Retrieve row and column indices of non-zero entries
  xs <- Matrix::summary(x)
  i <- xs$i
  j <- xs$j

  # Split column indices by row indices
  # -> jsplit will have exactly n entries
  # -> This is almost equivalent to the call `split(j[iord], i[iord] + 1)`
  #    except that rows with zero ones result in an empty vector
  #    whereas they would not appear in the `split` call.
  iord <- order(i)
  iord_rle <- rle(i[iord] + 1L)
  iord_rle_cs <- c(1L, cumsum(iord_rle$lengths))
  jord <- j[iord]
  jsplit <- vector(mode = "list", n)
  m <- 1L
  len_iuniq <- length(iord_rle$values)
  for (l in seq_len(n)) {
    if (m > len_iuniq) {
      break
    }

    if (iord_rle$values[m] == l) {
      jsplit[[l]] <- jord[iord_rle_cs[m]:iord_rle_cs[m + 1L]]
      m <- m + 1L
    } else {
      jsplit[[l]] <- vector("integer", 0L)
    }
  }

  # Run actual computation of Jaccard distances and save those
  # entries that have distance below the upper_bnd.
  out <- jaccard_indicator_comp(
    jsplit,
    eps = upper_bnd
  )

  # Form the indicator matrix
  methods::as(Matrix::sparseMatrix(
      c(out$i, out$j),
      c(out$j, out$i),
      dims = c(n, n)
  ) + Matrix::Diagonal(n), "ngCMatrix")
}

#' Determine initial centers for the kmeans++ algorithm
#'
#' @param x data matrix to be clustered
#' @param dm distance matrix (between rows of x; of class "dist")
#'
#' @return Row indices of initial cluster centers of x
#'
#' @keywords internal
kmeanspp_init <- function(n_cluster, x = NULL, dm = NULL) {
  if (sum(c(is.null(x), is.null(dm))) %in% c(0L, 2L)) {
    stop("Exactly one of x or dm needs to be supplied")
  }

  if (!is.null(x)) {
    dm <- dist(x)
  }

  n <- attr(dm, "Size")

  centers <- sample(n, size = 1L)
  for (i in 2L:n_cluster) {
    remaining_obs <- setdiff(seq_len(n), centers)
    log_ws_sq <- log(apply(do.call(
      cbind, lapply(
        centers,
        function(c) {
          lower_idx <- remaining_obs[remaining_obs < c]
          upper_idx <- remaining_obs[remaining_obs > c]

          c(
            dm[
              n * (lower_idx - 1)
              - lower_idx * (lower_idx - 1) / 2
              + c
              - lower_idx
            ],
            dm[
              n * (c - 1)
              - c * (c - 1) / 2
              + upper_idx
              - c
            ]
          )
          # # More straight-forward but less memory efficient method
          # dist_vals2 <- as.vector(as.matrix(dm)[remaining_obs, c])
        }
      )
    ), 1, min)^2)
    max_log_ws_sq <- max(log_ws_sq)
    ps <- (
      exp(log_ws_sq - max_log_ws_sq) / sum(exp(log_ws_sq - max_log_ws_sq))
    )
    centers <- c(centers, remaining_obs[
      sample.int(length(remaining_obs), size = 1, prob = ps)
    ])
  }

  centers
}

#' Perform the k-means++ algorithm
#'
#' Performs the k-means++ algorithm to cluster the rows of the input matrix.
#'
#' Estimation is repeated
#'
#' @param x Input matrix (n x p)
#' @param n_cluster Number of clusters
#' @param n_init_clusterings Number of repeated random initializations
#'                           to perform
#' @param n_max_iter Number of maximum iterations to perform in the k-means
#'                   algorithm
#'
#' @return An object of class [`stats::kmeans`].
#'
#' @references
#' David Arthur and Sergei Vassilvitskii. K-Means++: The advantages
#' of careful seeding. In Proceedings of the Eighteenth Annual ACM-SIAM
#' Symposium on Discrete Algorithms, SODA '07, pages 1027––1035.
#' Society for Industrial and Applied Mathematics, 2007.
#'
#' @concept helpers
#'
#' @export
kmeanspp <- function(x, n_cluster, n_init_clusterings = 10L, n_max_iter = 10L) {
  dm <- dist(x)
  initial_center_indices <- lapply(
    seq_len(n_init_clusterings),
    function(i) {
      kmeanspp_init(n_cluster, dm = dm)
    }
  )
  # Remove reference to dm
  dm <- NULL

  clusterings <- lapply(
    initial_center_indices,
    function(center_idx) {
      stats::kmeans(
        x,
        centers = x[center_idx, , drop = FALSE],
        iter.max = n_max_iter
      )
    }
  )

  min_idx <- which.min(sapply(clusterings, function(cl) cl$tot.withinss))
  clusterings[[min_idx]]$cluster
}

#' Determine module sizes
#'
#' @param module Vector of module indices
#' @param n_modules Total number of modules
#'
#' @return A named vector containing the name of the module (its index or
#'         `"Noise"`) and the number of elements in that module
#'
#' @concept helpers
#'
#' @export
find_module_sizes <- function(module, n_modules) {
  sapply(c(-1L, seq_len(n_modules)), function(i) {
    v <- sum(module == i)
    if (i == -1) {
      names(v) <- "Noise"
    } else {
      names(v) <- i
    }
    v
  })
}

#' Remove empty modules
#'
#' @details
#' Only iterates through modules with positive index, leaving the noise
#' module untouched.
#'
#' @param module Vector of module indices
#'
#' @return The updated vector of module indices with empty modules removed.
#'
#' @keywords internal
remove_empty_modules <- function(module) {
  module_ <- module
  if (max(module) > length(unique(module[module > 0]))) {
    unique_module <- unique(module[module > 0])
    for (i in seq_len(length(unique_module))) {
      module_[which(module == unique_module[i])] <- i
    }
  }

  module_
}

#' Extract target gene modules for given penalization parameters
#'
#' @param fit An object of class `scregclust`
#' @param penalization A numeric vector of penalization parameters.
#'                     The penalization parameters specified here must have
#'                     been used used during fitting of the `fit` object.
#'
#' @return A list of lists of final target modules. One list for each
#'         parameter in `penalization`. The lists contain the modules of
#'         target genes for each final configuration.
#'
#' @concept utilities
#'
#' @export
get_target_gene_modules <- function(fit, penalization = NULL) {
  if (!all(penalization %in% fit$penalization)) {
    cli::cli_abort(c(
      "Not all parameter values in {.var penalization} have been fitted.",
      "i" = paste(
        "Penalization parameters in {.class scregclust} object:",
        "{fit$penalization}"
      ),
      "i" = "Penalization parameters provided: {penalization}"
    ))
  }

  if (is.null(penalization)) {
    idx <- seq_along(fit$penalization)
  } else {
    idx <- which(fit$penalization %in% penalization)
  }

  lapply(idx, function(i) {
    lapply(
      fit$results[[i]]$output,
      function(o) {
        o$module[!fit$results[[i]]$is_regulator]
      }
    )
  })
}

#' Create a table of module overlap for two clusterings
#'
#' Compares two clusterings and creates a table of overlap between them.
#' Module labels do not have to match.
#'
#' @param k1 First clustering
#' @param k2 Second clustering
#'
#' @return A matrix showing the module overlap with the labels of `k1` in
#'         the columns and the labels of `k2` in the rows.
#'
#' @concept helpers
#'
#' @export
cluster_overlap <- function(k1, k2) {
  if (length(k1) != length(k2)) {
    cli::cli_abort(c(
      "Clusterings are not the same length.",
      "i" = "Length of {.var k1}: {length(k1)}",
      "i" = "Length of {.var k2}: {length(k2)}"
    ))
  }

  e_k1 <- sort(unique(k1))
  e_k2 <- sort(unique(k2))

  out <- do.call(cbind, lapply(e_k1, function(i1) {
    stats::setNames(vapply(e_k2, function(i2) {
      sum((k1 == i1) & (k2 == i2))
    }, 1L), e_k2)
  }))
  colnames(out) <- e_k1

  out
}

#' Extract final configurations into a data frame
#'
#' @param obj An object of class `scregclust`
#'
#' @return A [`data.frame`] containing penalization parameters and
#'         final configurations for those penalizations.
#'
#' @concept helpers
#'
#' @export
available_results <- function(obj) {
  data.frame(
    penalization = obj$penalization,
    final_configurations = sapply(obj$results, function(res) length(res$output))
  )
}

#' Fast computation of correlation
#'
#' This uses a more memory-intensive but much faster algorithm than
#' the built-in `cor` function.
#'
#' Computes the correlation between the columns of `x` and `y`.
#'
#' @param x first input matrix
#' @param y second input matrix
#'
#' @return Correlations matrix between the columns of `x` and `y`
#'
#' @concept helpers
#'
#' @export
fast_cor <- function(x, y) {
  xv <- scale(x, center = TRUE, scale = FALSE)
  yv <- scale(y, center = TRUE, scale = FALSE)
  xvss <- colSums(xv * xv)
  yvss <- colSums(yv * yv)
  result <- crossprod(xv, yv) / sqrt(outer(xvss, yvss))

  pmax(pmin(result, 1), -1)
}

#' Return the number of final configurations
#'
#' Returns the number of final configurations per penalization parameter in an
#' scRegClust object.
#'
#' @param fit An object of class `scRegClust`
#'
#' @return An integer vector containing the number of final configurations
#'         for each penalization parameter.
#'
#' @concept utilities
#'
#' @export
get_num_final_configs <- function(fit) {
  sapply(fit$results, function(r) length(r$output))
}

#' Get the average number of active regulators per module
#'
#' @param fit An object of class `scRegClust`
#'
#' @return A [`data.frame`] containing the average number of active regulators
#'         per module for each penalization parameter.
#'
#' @concept utilities
#'
#' @export
get_avg_num_regulators <- function(fit) {
  as.data.frame(do.call(rbind, lapply(fit$results, function(r) {
    c(
      penalization = r$penalization,
      colMeans(
        do.call(rbind, lapply(r$output, function(o) {
          stats::setNames(
            colSums(o$models),
            seq_len(ncol(o$models))
          )
        }))
      )
    )
  })))
}

#' Compute the Rand index
#'
#' @param k1 First clustering as vector of integers
#' @param k2 Second clustering as vector of integers
#'
#' @return The Rand index as a numeric value
#'
#' @references
#' W. M. Rand (1971). "Objective criteria for the evaluation of clustering
#' methods". Journal of the American Statistical Association 66 (336): 846–850.
#' DOI:10.2307/2284239
#'
#' @keywords internal
compute_rand_index <- function(k1, k2) {
  n <- length(k1)

  # Assertion
  stopifnot(length(k2) == n)
  stopifnot(is.numeric(k1), all(as.integer(k1) == k1))
  stopifnot(is.numeric(k2), all(as.integer(k2) == k2))

  # Requires that k1 and k2 are integer vectors (or integers in numeric format)
  m1 <- do.call(c, lapply(
    seq_len(n - 1L), function(i) abs(k1[i] - k1[(i + 1):n])
  ))
  m2 <- do.call(c, lapply(
    seq_len(n - 1L), function(i) abs(k2[i] - k2[(i + 1):n])
  ))

  # Compute Rand index
  (sum(!m1 & !m2) + sum(m1 & m2)) / choose(n, 2)
}

#' Compute Hubert's and Arabie's Adjusted Rand index
#'
#' @param k1 First clustering as vector of integers
#' @param k2 Second clustering as vector of integers
#'
#' @return The Adjusted Rand index as a numeric value
#'
#' @references
#' Lawrence Hubert and Phipps Arabie (1985). "Comparing partitions".
#' Journal of Classification. 2 (1): 193–218. DOI:10.1007/BF01908075
#'
#' @keywords internal
compute_adjusted_rand_index <- function(k1, k2) {
  n <- length(k1)

  # Assertion
  stopifnot(length(k2) == n)
  stopifnot(is.numeric(k1), all(as.integer(k1) == k1))
  stopifnot(is.numeric(k2), all(as.integer(k2) == k2))

  # Construct contingency table
  ct <- table(k1, k2)

  # Compute binomial pair sums
  sum_as <- sum(choose(rowSums(ct), 2))
  sum_bs <- sum(choose(colSums(ct), 2))
  sum_ns <- sum(choose(as.vector(ct), 2))
  denom <- choose(n, 2)

  # Compute adjusted Rand index
  (
    (sum_ns - (sum_as * sum_bs) / denom)
    / (
      0.5 * (sum_as + sum_bs)
      - (sum_as * sum_bs) / denom
    )
  )
}

#' Compute Rand indices
#'
#' Compute Rand indices for fitted scregclust object
#'
#' @param fit An object of class `scregclust`
#' @param groundtruth A known clustering of the target genes (integer vector)
#' @param adjusted If TRUE, the Adjusted Rand index is computed. Otherwise the
#'                 ordinary Rand index is computed.
#'
#' @return A [`data.frame`] containing the Rand indices. Since there can
#'         be more than one final configuration for some penalization
#'         parameters, Rand indices are averaged for each fixed penalization
#'         parameter. Returned are the mean, standard deviation and number
#'         of final configurations that were averaged.
#'
#' @references
#' W. M. Rand (1971). "Objective criteria for the evaluation of clustering
#' methods". Journal of the American Statistical Association 66 (336): 846–850.
#' DOI:10.2307/2284239
#'
#' Lawrence Hubert and Phipps Arabie (1985). "Comparing partitions".
#' Journal of Classification. 2 (1): 193–218. DOI:10.1007/BF01908075
#'
#' @concept utilities
#'
#' @export
get_rand_indices <- function(fit, groundtruth, adjusted = TRUE) {
  df <- do.call(rbind, lapply(get_target_gene_modules(fit), function(cs) {
    indices <- sapply(cs, function(cl) {
      noise_idx <- which(cl == -1)
      if (length(noise_idx) > 0) {
        cl_ <- cl[-noise_idx]
        gt_ <- groundtruth[-noise_idx]
      } else {
        cl_ <- cl
        gt_ <- groundtruth
      }

      if (length(cl_) > 0) {
        if (adjusted) {
          compute_adjusted_rand_index(gt_, cl_) * length(cl_) / length(cl)
        } else {
          compute_rand_index(gt_, cl_) * length(cl_) / length(cl)
        }
      } else {
        c(0)
      }
    })

    data.frame(mean = mean(indices), sd = sd(indices), n = length(indices))
  }))

  cbind(data.frame(penalization = fit$penalization), df)
}
