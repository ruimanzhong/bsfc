#' Initialize Graph and Compute Minimum Spanning Tree
#'
#' This function creates a graph from spatial data using various methods and computes its minimum spanning tree (MST).
#' It is designed to handle spatial data for clustering or network analysis, supporting multiple graph construction methods.
#'
#' @param map An sf object representing the spatial data, from which coordinates are extracted.
#' @param method Character string indicating the method used to construct the graph.
#'        Options include 'knn' (k-nearest neighbors), 'knn_geo' (geographical k-nearest neighbors),
#'        'rnn' (range nearest neighbors), 'rnn_geo' (geographical range nearest neighbors),
#'        'deltri' (Delaunay triangulation), and 'adjmat' (adjacency matrix based on spatial contiguity).
#' @param para Numeric value used as a parameter for the graph construction method, such as the number
#'        of neighbors (k) or a distance threshold (eps), depending on the method.
#' @param weights Optional numeric vector of weights for the edges in the graph. If NULL, random weights are assigned.
#' @param seed Integer value to set the random seed for reproducibility.
#'
#' @return A list with two elements: 'graph0' containing the initial graph object and 'mstgraph.ini'
#'         containing the initial minimum spanning tree computed from the graph.
#'
#' @examples
#' \dontrun{
#'   library(sf)
#'   library(igraph)
#'   # Assuming 'map' is an sf object containing spatial data
#'   result <- initial.mst.build(map, method = 'knn', para = 5)
#'   print(result$graph0)
#'   print(result$mstgraph.ini)
#' }
#'
#' @export
initial.mst.build <- function(map, method = 'knn', para = 10, weights = NULL, seed = 1234){
  set.seed(seed)
  graph0 = ConstructGraph0(map, method = method, para = para)
  if(is.null(weights)) { weights <- runif(length(E(graph0)), 0, 1) }
  mstgraph.ini = mst(graph0, weights = weights)
  return(list(graph0 = graph0, mstgraph.ini = mstgraph.ini))
}

#' Construct Initial Graph from Spatial Data
#'
#' This function constructs a graph from spatial coordinates using various methods,
#' including k-nearest neighbors, geographical k-nearest neighbors, and Delaunay triangulation.
#'
#' @param map An sf object from which coordinates are extracted to construct the graph.
#' @param method Character string specifying the method to use for graph construction.
#'        Available methods are 'knn', 'knn_geo', 'rnn', 'rnn_geo', 'deltri', and 'adjmat'.
#' @param para Parameter specific to the chosen method; varies by method, such as the number of neighbors or a distance threshold.
#' @param seed Integer value to set the random seed for reproducibility.
#'
#' @return An igraph object representing the constructed graph.
#'
#' @examples
#' \dontrun{
#'   library(sf)
#'   # Assuming 'map' is an sf object with spatial data
#'   graph <- ConstructGraph0(map, method = 'knn', para = 3)
#'   print(graph)
#' }
#'
#' @export
#'
ConstructGraph0=function(map,method='knn',para = 10, seed = 1234){
  set.seed(seed)
  coords <- st_coordinates(st_centroid(map)[,3])
  if(method=='knn')
  {
    k=para
    graph0=nng(coords,k=k);
    edgelist=ends(graph0,E(graph0))
    weight= sqrt(apply((coords[edgelist[,1],]-coords[edgelist[,2],])^2,1,sum));
    E(graph0)$weight=weight;}


  if(method=='knn_geo')
  {
    k=para
    graph0=nng(coords,k=k);
    edgelist=ends(graph0,E(graph0))
    weight= rdist.earth(coords[edgelist[,1],], coords[edgelist[,2],], miles=FALSE);
    E(graph0)$weight=weight;}

  if(method=='rnn')
  {
    eps=para
    adj=rdist(coords)
    adj[adj>eps]=0;
    graph0=graph_from_adjacency_matrix(adj,mode='upper',weighted=TRUE)

  }

  if(method=='rnn_geo')
  {
    eps=para
    adj=rdist.earth(coords, miles=FALSE)
    adj[adj>eps]=0;
    graph0=graph_from_adjacency_matrix(adj,mode='upper',weighted=TRUE)

  }

  if(method=='deltri'){
    threshold=para;
    graph0=dentrigraph(coords,threshold=threshold)

  }
  if(method == 'adjmat'){
    adjmat <- sf::st_touches(map) %>% as("matrix")
    graph0 <- graph.adjacency(adjmat, mode = "undirected", diag = FALSE)
  }
  return(graph0)
}
