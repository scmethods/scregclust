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
  REGtable <- output$reg_table
  idx <- !is.na(colSums(REGtable))
  REGtable <- REGtable[, idx]

  regulators <- c()
  for (i in seq_len(ncol(REGtable))) {
    tmp1 <- head(rownames(REGtable[order(REGtable[, i], decreasing = TRUE), ]))
    regulators <- append(regulators, tmp1)
    tmp2 <- tail(rownames(REGtable[order(REGtable[, i], decreasing = TRUE), ]))
    regulators <- append(regulators, tmp2)
  }

  regulators <- unique(regulators)

  f <- which(rownames(REGtable) %in% regulators)
  REGtable <- REGtable[f, ]

  REGtable$regulator <- rownames(REGtable)
  rownames(REGtable) <- NULL

  links <- reshape::melt(REGtable, id.vars = "regulator")
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

  rownames(REGtable) <- REGtable$regulator
  REGtable <- REGtable[, -ncol(REGtable)]

  nodes <- array(0, dim = c((nrow(REGtable) + ncol(REGtable)), 2))
  colnames(nodes) <- c("id", "type")

  nodes[seq_len(nrow(REGtable)), 1] <- rownames(REGtable)
  nodes[seq_len(nrow(REGtable)), 2] <- "Regulator"
  nodes[(nrow(REGtable) + 1):nrow(nodes), 1] <- colnames(REGtable)
  nodes[(nrow(REGtable) + 1):nrow(nodes), 2] <- "TargetState"
  nodes <- as.data.frame(nodes)

  net <- igraph::graph_from_data_frame(
    d = links, vertices = nodes, directed = TRUE
  )

  igraph::V(net)[which(igraph::V(net)$type == "Regulator")]$shape <- 1
  igraph::V(net)[which(igraph::V(net)$type == "TargetState")]$shape <- 2

  igraph::V(net)[which(igraph::V(net)$type == "Regulator")]$type <- 1
  igraph::V(net)[which(igraph::V(net)$type == "TargetState")]$type <- (
    seq_len(ncol(REGtable))
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

#' @export
plot.scregclust <- function(x, ...) {
  r2_cluster_data <- do.call(rbind, lapply(x$results, function(r) {
    do.call(rbind, lapply(r$output, function(o) {
      idx <- !is.na(o$r2_cluster)

      data.frame(
        penalization = r$penalization,
        cluster = seq_along(o$r2_cluster)[idx],
        value = o$r2_cluster[idx]
      )
    }))
  }))
  r2_cluster_data$penalization <- factor(
    r2_cluster_data$penalization, levels = x$penalization
  )
  r2_cluster_data$variable <- "r2-per-cluster"

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
          cluster = i,
          value = o$importance[idx, i]
        )
      }))
    }))
  }))
  importance_data$penalization <- factor(
    importance_data$penalization, levels = x$penalization
  )
  importance_data$variable <- "importance"

  rbind(r2_cluster_data, importance_data) |>
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
        k <- o$cluster[!r$is_regulator]

        order_list <- lapply(seq_len(r$n_cl), function(cl) {
          if (sum(k == cl) > 0) {
            order(o$silhouette[k == cl])
          } else {
            integer(0)
          }
        })
        gene <- do.call(c, lapply(seq_len(r$n_cl), function(cl) {
          seq_along(k)[k == cl][order_list[[cl]]]
        }))

        data.frame(
          order = seq_len(sum(k != -1)),
          gene = gene,
          silhouette = o$silhouette[gene],
          cluster = as.factor(k[gene]),
          n_cl = r$n_cl,
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
#'                     dataset across a variety of cluster counts.
#' @param penalization Either a single numeric value requesting the results
#'                     for the same penalty parameter across all fits in
#'                     `list_of_fits`, or one for each individual fit.
#' @param final_config The final configuration that should be visualised.
#'                     Either a single number to be used for all fits in
#'                     `list_of_fits`, or one for each individual fit.
#'
#' @return A [`ggplot2`] plot showing the the silhouette scores for each
#'         supplied fit.
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

  silhouette_data <- collect_silhouette_data(list_of_fits)
  cluster_counts <- sapply(list_of_fits, function(fit) fit$results[[1]]$n_cl)

  silhouette_data$n_cl_lbl <- as.factor(sprintf("K = %d", silhouette_data$n_cl))

  if (length(penalization) == 1L) {
    silhouette_data <- silhouette_data[
      silhouette_data$penalization == penalization,
    ]
  } else {
    silhouette_data <- do.call(rbind, lapply(
      seq_along(cluster_counts),
      function(i) {
        df <- silhouette_data[silhouette_data$n_cl == cluster_counts[i], ]
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
      seq_along(cluster_counts),
      function(i) {
        df <- silhouette_data[silhouette_data$n_cl == cluster_counts[i], ]
        df[df$output == final_config[i]]
      }
    ))
  }

  cluster_centers <- do.call(rbind, lapply(cluster_counts, function(n_cl) {
    df <- silhouette_data[silhouette_data$n_cl == n_cl, ]
    contained_clusters <- unique(df$cluster)

    data.frame(
      n_cl = n_cl,
      cluster = contained_clusters,
      order =  sapply(contained_clusters, function(cl) {
        mean(df[df$cluster == cl, ]$order)
      })
    )
  }))

  avg_silhouette <- data.frame(
    n_cl = cluster_counts,
    silhouette = sapply(cluster_counts, function(n_cl) {
      df <- silhouette_data[silhouette_data$n_cl == n_cl, ]
      mean(df$silhouette)
    })
  )
  avg_silhouette$n_cl_lbl <- as.factor(sprintf("K = %d", avg_silhouette$n_cl))

  silhouette_data |>
    ggplot2::ggplot() +
    ggplot2::facet_wrap(n_cl_lbl ~ .) +
    ggplot2::geom_bar(
      ggplot2::aes(x = .data$order, y = .data$silhouette, fill = .data$cluster),
      stat = "identity",
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = .data$order, y = -0.1, label = .data$cluster),
      data = cluster_centers,
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
#'                     dataset across a variety of cluster counts.
#' @param penalization Either a single numeric value requesting the results
#'                     for the same penalty parameter across all fits in
#'                     `list_of_fits`, or one for each individual fit.
#'
#' @return A [ggplot2] plot showing the average silhouette score and the
#'         average predictive \eqn{R^2}
#'
#' @export
plot_cluster_count_helper <- function(list_of_fits, penalization) {
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

  silhouette_data <- collect_silhouette_data(list_of_fits)

  avg_r2_cluster_data <- do.call(rbind, lapply(list_of_fits, function(fit) {
    do.call(rbind, lapply(seq_along(fit$results), function(i) {
      r <- fit$results[[i]]
      r2_cluster <- do.call(c, lapply(seq_along(r$output), function(j) {
        r$output[[j]]$r2_cluster
      })) # average across different configurations

      # If a cluster is empty then r2_cluster is NA, so use NA remove
      value <- mean(r2_cluster, na.rm = TRUE)
      # If all clusters turn out to be empty (e.g. too high penalization) then
      # mean(...) above will evaluate to NaN. Do not return a data.frame
      # in that case.
      if (is.nan(value)) {
        return(NULL)
      }

      data.frame(
        n_cl = r$n_cl,
        penalization = r$penalization,
        value = value,
        variable = "avg-r2-cluster"
      )
    }))
  }))

  cluster_counts <- sapply(list_of_fits, function(fit) fit$results[[1]]$n_cl)

  if (length(penalization) == 1) {
    silhouette_data <- silhouette_data[
      silhouette_data$penalization == penalization,
    ]
    avg_r2_cluster_data <- avg_r2_cluster_data[
      avg_r2_cluster_data$penalization == penalization,
    ]
  } else {
    silhouette_data <- do.call(rbind, lapply(
      seq_along(cluster_counts),
      function(i) {
        df <- silhouette_data[silhouette_data$n_cl == cluster_counts[i], ]
        df[df$penalization == penalization[i]]
      }
    ))
    avg_r2_cluster_data <- do.call(rbind, lapply(
      seq_along(cluster_counts),
      function(i) {
        df <- avg_r2_cluster_data[
          avg_r2_cluster_data$n_cl == cluster_counts[i],
        ]
        df[df$penalization == penalization[i]]
      }
    ))
  }

  avg_silhouette <- sapply(seq_along(cluster_counts), function(i) {
    df <- silhouette_data[silhouette_data$n_cl == cluster_counts[i], ]
    mean(df$silhouette) # average across different configurations
  })

  rbind(
    data.frame(
      n_cl = cluster_counts,
      penalization = penalization,
      value = avg_silhouette,
      variable = "avg-silhouette"
    ),
    avg_r2_cluster_data
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
      ggplot2::aes(.data$n_cl, .data$value), linewidth = 0.25
    ) +
    ggplot2::geom_point(
      ggplot2::aes(.data$n_cl, .data$value), size = 0.5
    ) +
    ggplot2::labs(x = "# of modules (K)", y = NULL) +
    ggplot2::scale_x_continuous(breaks = cluster_counts) +
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
