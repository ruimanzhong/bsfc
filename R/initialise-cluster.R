#' Initialize Graph and Compute Minimum Spanning Tree
#'
#' This function creates a graph from spatial data using various methods and computes its minimum spanning tree (MST).
#' It is designed to handle spatial data for clustering or network analysis, supporting multiple graph construction methods.
#'
#' @param x An sf object representing the spatial data, from which coordinates are extracted.
#' @param nclust Integer, the initial number of clusters
#' @param weights Optional numeric vector of weights for the edges in the graph. If NULL, random weights are assigned.
#'
#' @return A list with three elements: 'graph' containing the initial undirected graph
#' object, 'mst' containing the initial minimum spanning tree, and 'cluster' containing
#' the membership of `x`.
#'
#' @examples
#'
#'   librar(sf)
#'   x <- st_make_grid(cellsize = c(1, 1), offset = c(0, 0), n = c(3, 2))
#'   cluster_ini <- initial_cluster(x, nclust = 3, weights = 1:6))
#'   print(cluster_ini)
#'
#' @export
initial_cluster <- function(x, nclust = 10, weights = NULL){

  # create adjacency, initial checks and weights if required
  if (inherits(x, c("sf", "sfc"))) {
    x <- as(st_touches(st_geometry(x)), "matrix")
  } else if (!inherits(x, c("matrix", "Matrix"))) {
    stop("`x` must be of class `sf`, `sfc`, `matrix` or `Matrix`.")
  }

  if (nclust > dim(x)[1]) {
    stop("`nclust` must be smaller that number of regions.")
  }

  if (is.null(weights)) weights <- runif(length(x))

  # create weighted graph and minimum spanning tree
  graph <- graph_from_adjacency_matrix(x * weights, mode = "upper", weighted = TRUE)
  mstgraph <- mst(graph)
  V(mstgraph)$vid <- 1:vcount(mstgraph)

  # partition mst into nclust groups
  rmid <- order(E(mstgraph)$weight, decreasing = TRUE)[1:(nclust - 1)]
  partition <- components(delete_edges(mstgraph, rmid))

  return(list(graph = graph, mst = mstgraph, cluster = partition$membership))
}
