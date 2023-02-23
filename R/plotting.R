#' Plotting the regulatory table from scregclust as a directed graph
#'
#' @param fit Output list from scregclust
#' @param arrow_size Size of arrow head
#' @param edge_scaling Scaling factor for edge width
#' @param no_links Threshold value (0-10) for number of edges to show,
#'                 higher value = more stringent threshold = less edges
#' @param col color
#'
#' @return Graph with gene modules and regulators as nodes
#'
#' @export
plot_regulator_network <- function(fit,
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
  REGtable <- fit[[1]]
  idx <- which(is.na(colSums(fit[[1]])))
  REGtable <- REGtable[, -idx]

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

  links <- reshape::melt(as.matrix(REGtable))
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

  igraph::V(net)$shape[which(igraph::V(net)$type == "Regulator")] <- 1
  igraph::V(net)$shape[which(igraph::V(net)$type == "TargetState")] <- 2

  igraph::V(net)$type[which(igraph::V(net)$type == "Regulator")] <- 1
  igraph::V(net)$type[which(igraph::V(net)$type == "TargetState")] <- (
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
    x = -1.5,
    y = -1.1,
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
