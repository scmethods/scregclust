#' @details
#' Computational methods for the scregclust algorithm
#' @keywords internal
#' @aliases scregclust-package
#' @author Ida Larsson, Felix Held, Sven Nelander
#' @import Rcpp
#' @import cli
#' @import ggplot2
#' @importFrom prettyunits pretty_dt
#' @importFrom Matrix Matrix sparseMatrix Diagonal t rowSums colSums summary
#' @importFrom stats cor coef predict na.omit kmeans quantile sd dist setNames
#' @importFrom utils read.table head tail
#' @importFrom graphics legend
#' @importFrom methods is as
#' @importFrom reshape melt
#' @importFrom igraph graph_from_data_frame delete_edges delete_vertices layout_with_fr V E degree
#' @importFrom grid arrow unit
#' @useDynLib scregclust, .registration = TRUE
"_PACKAGE"

.onUnload <- function(libpath) {
    library.dynam.unload("scregclust", libpath)
}