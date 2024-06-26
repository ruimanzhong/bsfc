% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initial_build.R
\name{initial.mst.build}
\alias{initial.mst.build}
\title{Initialize Graph and Compute Minimum Spanning Tree}
\usage{
initial.mst.build(
  map,
  method = "knn",
  para = 10,
  nclust = 10,
  weights = NULL,
  seed = 1234
)
}
\arguments{
\item{map}{An sf object representing the spatial data, from which coordinates are extracted.}

\item{method}{Character string indicating the method used to construct the graph.
Options include 'knn' (k-nearest neighbors), 'knn_geo' (geographical k-nearest neighbors),
'rnn' (range nearest neighbors), 'rnn_geo' (geographical range nearest neighbors),
'deltri' (Delaunay triangulation), and 'adjmat' (adjacency matrix based on spatial contiguity).}

\item{para}{Numeric value used as a parameter for the graph construction method, such as the number
of neighbors (k) or a distance threshold (eps), depending on the method.}

\item{nclust}{Integer, the initial number of clusters}

\item{weights}{Optional numeric vector of weights for the edges in the graph. If NULL, random weights are assigned.}

\item{seed}{Integer value to set the random seed for reproducibility.}
}
\value{
A list with two elements: 'graph0' containing the initial graph object and 'mstgraph.ini'
containing the initial minimum spanning tree computed from the graph.
}
\description{
This function creates a graph from spatial data using various methods and computes its minimum spanning tree (MST).
It is designed to handle spatial data for clustering or network analysis, supporting multiple graph construction methods.
}
\examples{
\dontrun{
  library(sf)
  library(igraph)
  # Assuming 'map' is an sf object containing spatial data
  result <- initial.mst.build(map, method = 'knn', para = 5, nclust = 5)
  print(result$graph0)
  print(result$mstgraph.ini)
}

}
