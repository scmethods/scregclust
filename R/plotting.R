#' Plotting the regulatory table from scregclust as a directed graph
#'
#' @param output Object of type `scregclust_output` from a fit of the
#'               scregclust algorithm.
#' @param arrow_size Size of arrow head
#' @param edge_scaling Scaling factor for edge width
#' @param no_links Threshold value (0-10) for number of edges to show,
#'                 higher value = more stringent threshold = less edges
#' @param col color
#'
#' @return Graph with gene modules and regulators as nodes
#'
#' @concept plotting
#'
#' @export
plot_regulator_network <- function(output,
                                   arrow_size = 0.3,
                                   edge_scaling = 30,
                                   no_links = 6,
                                   col = c(
                                     "gray80",
                                     "#FC7165",
                                     "#BD828C",
                                     "#9D8A9F",
                                     "#7D92B2",
                                     "#BDA88C",
                                     "#FCBD65",
                                     "#F2BB90",
                                     "#E7B9BA",
                                     "#BDB69C",
                                     "#92B27D",
                                     "#9B8BA5",
                                     "#9D7DB2",
                                     "#94A5BF"
                                   )) {
  reg_table <- output$reg_table
  idx <- !is.na(colSums(reg_table))
  reg_table <- reg_table[, idx]

  regulators <- c()
  for (i in seq_len(ncol(reg_table))) {
    tmp1 <- head(rownames(
      reg_table[order(reg_table[, i], decreasing = TRUE), ]
    ))
    regulators <- append(regulators, tmp1)
    tmp2 <- tail(rownames(
      reg_table[order(reg_table[, i], decreasing = TRUE), ]
    ))
    regulators <- append(regulators, tmp2)
  }

  regulators <- unique(regulators)

  f <- which(rownames(reg_table) %in% regulators)
  reg_table <- reg_table[f, ]

  reg_table$regulator <- rownames(reg_table)
  rownames(reg_table) <- NULL

  links <- reshape::melt(reg_table, id.vars = "regulator")
  colnames(links) <- c("from", "to", "weight")
  f <- which(links$weight == 0)
  links <- links[-f, ]

  m <- which(links$weight < 0)
  p <- which(links$weight > 0)

  links$mode <- array(0, dim = c(nrow(links), 1))
  links$mode[m] <- "Repress"
  links$mode[p] <- "Activate"
  links$color <- array(0, dim = c(nrow(links), 1))
  links$color[m] <- "#2B278C"
  links$color[p] <- "#BD111F"
  links$weight <- abs(links$weight)

  links <- as.data.frame(links)

  rownames(reg_table) <- reg_table$regulator
  reg_table <- reg_table[, -ncol(reg_table)]

  nodes <- array(0, dim = c((nrow(reg_table) + ncol(reg_table)), 2))
  colnames(nodes) <- c("id", "type")

  nodes[seq_len(nrow(reg_table)), 1] <- rownames(reg_table)
  nodes[seq_len(nrow(reg_table)), 2] <- "Regulator"
  nodes[(nrow(reg_table) + 1):nrow(nodes), 1] <- colnames(reg_table)
  nodes[(nrow(reg_table) + 1):nrow(nodes), 2] <- "TargetState"
  nodes <- as.data.frame(nodes)

  net <- igraph::graph_from_data_frame(
    d = links, vertices = nodes, directed = TRUE
  )

  igraph::V(net)[which(igraph::V(net)$type == "Regulator")]$shape <- 1
  igraph::V(net)[which(igraph::V(net)$type == "TargetState")]$shape <- 2

  igraph::V(net)[which(igraph::V(net)$type == "Regulator")]$type <- 1
  igraph::V(net)[which(igraph::V(net)$type == "TargetState")]$type <- (
    seq_len(ncol(reg_table))
  )

  colrs <- col
  igraph::V(net)$color <- colrs[as.numeric(igraph::V(net)$type)]

  cut.off <- quantile(links$weight, probs = seq(0, 1, 0.1))[no_links]
  net <- igraph::delete_edges(net, igraph::E(net)[links$weight < cut.off])

  isolated <- which(igraph::degree(net) == 0)
  net <- igraph::delete_vertices(net, isolated)

  igraph::E(net)$arrow.size <- arrow_size
  igraph::V(net)$shape <- c("vrectangle", "circle")[
    as.numeric(igraph::V(net)$shape)
  ]
  igraph::E(net)$width <- igraph::E(net)$weight * edge_scaling

  l <- igraph::layout_with_fr(net)

  plot(
    net,
    layout = l,
    edge.curved = 0.3,
    vertex.label.cex = .6,
    vertex.label.color = "black",
    alpha = 0.5
  )
  legend(
    x = -1.1,
    y = -0.8,
    c("Activating", "Repressing"),
    pch = 21,
    col = "#777777",
    pt.bg = c("#BD111F", "#2B278C"),
    pt.cex = 2,
    cex = .8,
    bty = "n",
    ncol = 1
  )
}

#' @concept plotting
#'
#' @export
plot.scregclust <- function(x, ...) {
  r2_module_data <- do.call(rbind, lapply(x$results, function(r) {
    do.call(rbind, lapply(r$output, function(o) {
      idx <- !is.na(o$r2_module)

      data.frame(
        penalization = r$penalization,
        module = seq_along(o$r2_module)[idx],
        value = o$r2_module[idx]
      )
    }))
  }))
  r2_module_data$penalization <- factor(
    r2_module_data$penalization, levels = x$penalization
  )
  r2_module_data$variable <- "r2-per-module"

  importance_data <- do.call(rbind, lapply(x$results, function(r) {
    do.call(rbind, lapply(seq_along(r$output), function(j) {
      o <- r$output[[j]]
      do.call(rbind, lapply(seq_len(ncol(o$models)), function(i) {
        idx <- !is.na(o$importance[, i])
        if (sum(idx) == 0) {
          return(NULL)
        }

        data.frame(
          penalization = r$penalization,
          module = i,
          value = o$importance[idx, i]
        )
      }))
    }))
  }))
  importance_data$penalization <- factor(
    importance_data$penalization, levels = x$penalization
  )
  importance_data$variable <- "importance"

  rbind(r2_module_data, importance_data) |>
    ggplot2::ggplot() +
    ggplot2::facet_wrap(
      variable ~ .,
      nrow = 2,
      scales = "free_y",
      strip.position = "left",
      labeller = ggplot2::label_bquote(
        .(
          if (variable == "importance") {
            "Regulator Importance"
          } else {
            "Predictive" ~ R^2 ~ "per module"
          }
        )
      ),
    ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(x = .data$penalization, y = .data$value),
      outlier.size = 0.5,
      lwd = 0.25,
    ) +
    ggplot2::labs(x = "Penalization", y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(
        arrow = grid::arrow(length = grid::unit(1, "mm")),
      ),
      strip.background = ggplot2::element_blank(),
      strip.placement = "outside",
      line = ggplot2::element_line(linewidth = 0.25),
      plot.margin = ggplot2::margin(t = 2, unit = "mm"),
    )
}

collect_silhouette_data <- function(list_of_fits) {
  do.call(rbind, lapply(list_of_fits, function(fit) {
    do.call(rbind, lapply(seq_along(fit$results), function(i) {
      r <- fit$results[[i]]
      do.call(rbind, lapply(seq_along(r$output), function(j) {
        o <- r$output[[j]]
        k <- o$module[!r$is_regulator]

        order_list <- lapply(seq_len(r$n_modules), function(cl) {
          if (sum(k == cl) > 0) {
            order(o$silhouette[k == cl])
          } else {
            integer(0)
          }
        })
        gene <- do.call(c, lapply(seq_len(r$n_modules), function(cl) {
          seq_along(k)[k == cl][order_list[[cl]]]
        }))

        data.frame(
          order = seq_len(sum(k != -1)),
          gene = gene,
          silhouette = o$silhouette[gene],
          module = as.factor(k[gene]),
          n_modules = r$n_modules,
          output = j,
          penalization = r$penalization
        )
      }))
    }))
  }))
}

#' Plot individual silhouette scores
#'
#' @param list_of_fits A list of `scregclust` objects each fit to the same
#'                     dataset across a variety of module counts (varying
#'                     `n_modules` when running [`scregclust`]).
#' @param penalization Either a single numeric value requesting the results
#'                     for the same penalty parameter across all fits in
#'                     `list_of_fits`, or one for each individual fit.
#' @param final_config The final configuration that should be visualized.
#'                     Either a single number to be used for all fits in
#'                     `list_of_fits`, or one for each individual fit.
#'
#' @return A ggplot2 plot showing the the silhouette scores for each
#'         supplied fit.
#'
#' @concept plotting
#'
#' @export
plot_silhouettes <- function(list_of_fits, penalization, final_config = 1L) {
  if (!(
    is.numeric(penalization)
    && (
      (
        length(penalization) == 1L
        && all(sapply(list_of_fits, function(fit) {
          penalization %in% fit$penalization
        }))
      ) || (
        length(penalization) == length(list_of_fits)
        && all(mapply(function(fit, p) {
          p %in% fit$penalization
        }, list_of_fits, penalization))
      )
    )
  )) {
    cli::cli_abort(c(
      "{.var penalization} is not supplied correctly.",
      "x" = "It needs to be one of the following two:",
      "*" = "A single penalization parameter used in all fits.",
      "*" = (
        "A list of penalization parameters, exactly one for each supplied fit."
      )
    ))
  }

  #### TODO: Checking the correctness of this is a bit of a pain
  ####       Do soon-ish!
  # if (!(
  #   is.numeric(final_config)
  #   && all(as.integer(final_config) == final_config)
  #   && (
  #     (
  #       length(final_config) == 1L
  #       && all(sapply(list_of_fits, function(fit) {
  #         final_config %in% fit$final_config
  #       }))
  #     ) || (
  #       length(final_config) == length(list_of_fits)
  #       && all(mapply(function(fit, p) {
  #         p %in% fit$final_config
  #       }, list_of_fits, final_config))
  #     )
  #   )
  # )) {
  #   cli::cli_abort(c(
  #     "{.var final_config} is not supplied correctly.",
  #     "x" = "It needs to be one of the following two:",
  #     "*" = "A single final_config parameter used in all fits.",
  #     "*" = (
  #       "A list of final_config parameters, exactly one for each supplied fit."
  #     )
  #   ))
  # }

  if (any(
    do.call(c, lapply(list_of_fits, function(fit) {
      do.call(c, lapply(fit$results, function(res) {
        sapply(res$output, function(o) {
          is.null(o$silhouette)
        })
      }))
    }))
  )) {
    cli::cli_abort(c(
      "Silhouette scores were not computed during fitting.",
      "i" = "Set `compute_silhouette = TRUE` in `scregclust`"
    ))
  }

  silhouette_data <- collect_silhouette_data(list_of_fits)
  module_counts <- sapply(
    list_of_fits, function(fit) fit$results[[1]]$n_modules
  )

  silhouette_data$n_modules_lbl <- as.factor(
    sprintf("K = %d", silhouette_data$n_modules)
  )

  if (length(penalization) == 1L) {
    silhouette_data <- silhouette_data[
      silhouette_data$penalization == penalization,
    ]
  } else {
    silhouette_data <- do.call(rbind, lapply(
      seq_along(module_counts),
      function(i) {
        df <- silhouette_data[silhouette_data$n_modules == module_counts[i], ]
        df[df$penalization == penalization[i]]
      }
    ))
  }

  if (length(final_config) == 1L) {
    silhouette_data <- silhouette_data[
      silhouette_data$output == final_config,
    ]
  } else {
    silhouette_data <- do.call(rbind, lapply(
      seq_along(module_counts),
      function(i) {
        df <- silhouette_data[silhouette_data$n_modules == module_counts[i], ]
        df[df$output == final_config[i]]
      }
    ))
  }

  module_centers <- do.call(rbind, lapply(module_counts, function(n_modules) {
    df <- silhouette_data[silhouette_data$n_modules == n_modules, ]
    contained_modules <- unique(df$module)

    data.frame(
      n_modules = n_modules,
      module = contained_modules,
      order =  sapply(contained_modules, function(cl) {
        mean(df[df$module == cl, ]$order)
      })
    )
  }))
  module_centers$n_modules_lbl <- as.factor(
    sprintf("K = %d", module_centers$n_modules)
  )

  avg_silhouette <- data.frame(
    n_modules = module_counts,
    silhouette = sapply(module_counts, function(n_modules) {
      df <- silhouette_data[silhouette_data$n_modules == n_modules, ]
      mean(df$silhouette)
    })
  )
  avg_silhouette$n_modules_lbl <- as.factor(
    sprintf("K = %d", avg_silhouette$n_modules)
  )

  silhouette_data |>
    ggplot2::ggplot() +
    ggplot2::facet_wrap(n_modules_lbl ~ .) +
    ggplot2::geom_bar(
      ggplot2::aes(x = .data$order, y = .data$silhouette, fill = .data$module),
      stat = "identity",
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = .data$order, y = -0.1, label = .data$module),
      data = module_centers,
    ) +
    ggplot2::geom_hline(
      ggplot2::aes(yintercept = .data$silhouette),
      data = avg_silhouette,
      linetype = "dashed",
      color = "red",
      linewidth = 0.25,
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_discrete(guide = "none") +
    ggplot2::labs(x = "Module", y = "Silhouette score") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
    )
}

#' Plot average silhouette scores and average predictive \eqn{R^2}
#'
#' @param list_of_fits A list of `scregclust` objects each fit to the same
#'                     dataset across a variety of module counts (varying
#'                     `n_modules` while running [`scregclust`]).
#' @param penalization Either a single numeric value requesting the results
#'                     for the same penalty parameter across all fits in
#'                     `list_of_fits`, or one for each individual fit.
#'
#' @return A ggplot2 plot showing the average silhouette score and the
#'         average predictive \eqn{R^2}
#'
#' @concept plotting
#'
#' @export
plot_module_count_helper <- function(list_of_fits, penalization) {
  if (!(
    is.list(list_of_fits)
    && all(sapply(list_of_fits, function(f) "scregclust" %in% class(f)))
  )) {
    cli::cli_abort(c(
      "{.var list_of_fits} is not supplied correctly.",
      "x" = "It needs to be a list of {.class scregclust} objects."
    ))
  }

  if (!(
    is.numeric(penalization)
    && (
      (
        length(penalization) == 1L
        && all(sapply(list_of_fits, function(fit) {
          penalization %in% fit$penalization
        }))
      ) || (
        length(penalization) == length(list_of_fits)
        && all(mapply(function(fit, p) {
          p %in% fit$penalization
        }, list_of_fits, penalization))
      )
    )
  )) {
    cli::cli_abort(c(
      "{.var penalization} is not supplied correctly.",
      "x" = "It needs to be one of the following two:",
      "*" = "A single penalization parameter used in all fits.",
      "*" = (
        "A list of penalization parameters, exactly one for each supplied fit."
      )
    ))
  }

  if (any(
    do.call(c, lapply(list_of_fits, function(fit) {
      do.call(c, lapply(fit$results, function(res) {
        sapply(res$output, function(o) {
          is.null(o$silhouette)
        })
      }))
    }))
  )) {
    cli::cli_abort(c(
      "Silhouette scores were notcomputed during fitting.",
      "i" = "Set `compute_silhouette = TRUE` in `scregclust`"
    ))
  }

  silhouette_data <- collect_silhouette_data(list_of_fits)

  avg_r2_module_data <- do.call(rbind, lapply(list_of_fits, function(fit) {
    do.call(rbind, lapply(seq_along(fit$results), function(i) {
      r <- fit$results[[i]]
      r2_module <- do.call(c, lapply(seq_along(r$output), function(j) {
        r$output[[j]]$r2_module
      })) # average across different configurations

      # If a module is empty then r2_module is NA, so use NA remove
      value <- mean(r2_module, na.rm = TRUE)
      # If all modules turn out to be empty (e.g. too high penalization) then
      # mean(...) above will evaluate to NaN. Do not return a data.frame
      # in that case.
      if (is.nan(value)) {
        return(NULL)
      }

      data.frame(
        n_modules = r$n_modules,
        penalization = r$penalization,
        value = value,
        variable = "avg-r2-module"
      )
    }))
  }))

  module_counts <- sapply(
    list_of_fits, function(fit) fit$results[[1]]$n_modules
  )

  if (length(penalization) == 1) {
    silhouette_data <- silhouette_data[
      silhouette_data$penalization == penalization,
    ]
    avg_r2_module_data <- avg_r2_module_data[
      avg_r2_module_data$penalization == penalization,
    ]
  } else {
    silhouette_data <- do.call(rbind, lapply(
      seq_along(module_counts),
      function(i) {
        df <- silhouette_data[silhouette_data$n_modules == module_counts[i], ]
        df[df$penalization == penalization[i]]
      }
    ))
    avg_r2_module_data <- do.call(rbind, lapply(
      seq_along(module_counts),
      function(i) {
        df <- avg_r2_module_data[
          avg_r2_module_data$n_modules == module_counts[i],
        ]
        df[df$penalization == penalization[i]]
      }
    ))
  }

  avg_silhouette <- sapply(seq_along(module_counts), function(i) {
    df <- silhouette_data[silhouette_data$n_modules == module_counts[i], ]
    mean(df$silhouette) # average across different configurations
  })

  rbind(
    data.frame(
      n_modules = module_counts,
      penalization = penalization,
      value = avg_silhouette,
      variable = "avg-silhouette"
    ),
    avg_r2_module_data
  ) |>
    ggplot2::ggplot() +
    ggplot2::facet_wrap(
      variable ~ .,
      nrow = 2,
      scales = "free_y",
      strip.position = "left",
      labeller = ggplot2::label_bquote(
        .(
          if (variable == "avg-silhouette") {
            "Average silhouette score"
          } else {
            "Avg. pred." ~ R^2 ~ "per module"
          }
        )
      ),
    ) +
    ggplot2::geom_line(
      ggplot2::aes(.data$n_modules, .data$value), linewidth = 0.25
    ) +
    ggplot2::geom_point(
      ggplot2::aes(.data$n_modules, .data$value), size = 0.5
    ) +
    ggplot2::labs(x = "# of modules (K)", y = NULL) +
    ggplot2::scale_x_continuous(breaks = module_counts) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(
        arrow = grid::arrow(length = grid::unit(1, "mm")),
      ),
      strip.background = ggplot2::element_blank(),
      strip.placement = "outside",
      line = ggplot2::element_line(linewidth = 0.25),
      plot.margin = ggplot2::margin(t = 2, unit = "mm"),
    )
}
