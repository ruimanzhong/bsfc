#' Plot Clustered Functions
#'
#' This function plots the results of a clustering analysis, displaying each cluster with its respective time series.
#' It uses pre-computed model results to add a line representing the mean fitted values for each cluster.
#' The function saves the plots to a specified file.
#'
#' @param Y Numeric matrix containing the data to be plotted, with each column representing a different time series.
#' @param nt Integer, the number of time points.
#' @param clust_res Numeric vector containing cluster assignments for each time series in Y.
#' @param final_model List of model objects for each cluster, each containing a summary of fitted values.
#' @param a,b plot arrange number
#' @param filepath Character string, the name of the file where the plot will be saved.
#'
#' @details The function generates a jpeg file with a series of plots arranged in a grid. Each plot corresponds
#' to a cluster and shows the individual time series in the cluster as well as the mean fitted values across
#' the cluster. The function dynamically adjusts the number of rows in the plot grid based on the number
#' of clusters, aiming to arrange the plots in roughly three columns.
#'
#' The individual time series are plotted in a semi-transparent gray to allow for overlap visibility,
#' while the mean of the model's fitted values for each cluster is plotted in red for emphasis.
#'
#' @return Invisible NULL. The function's primary output is a graphical file saved to the specified location.
#' The function does not return any value within R.
#'
#' @examples
#' # Assuming Y is a matrix of time series data, nt is the number of time points,
#' # clust_res contains cluster assignments, and final_model is a list of model summaries:
#' \dontrun{
#' plotClusterFun(Y, 50, cluster_results, models, "path/to/save/cluster_plot.jpeg")
#' }
#' @export
plotClusterFun <- function(Y,E, nt,ns, clust_res, final_model,page,  filepath) {
  p = max(clust_res)
  ydf = setNames(as.data.frame(Y / E), 1:ns) |>
    mutate(time = 1:nt) |>
    tidyr::pivot_longer(1:ns, names_to = "region", names_transform = as.numeric) |>
    mutate(cluster = clust_res[region])
  cluster2 = rep(1:p, table(clust_res))
  preddf = purrr::map(1:p, function(x) final_model[[x]]$summary.fitted.values$mean / exp(final_model[[x]]$summary.random[['id']]$mean)) %>%
    purrr::map(~ matrix(., nrow = nt)) %>%
    do.call(cbind, .) |>
    as.data.frame() |>
    setNames(1:ns) |>
    mutate(time =  1:nt) |>
    tidyr::pivot_longer(1:ns, names_to = "newregion", names_transform = as.numeric) |>
    mutate(cluster = factor(cluster2[newregion]))
  pdf(filepath, height = 11, width = 8)
  for (i in 1:page) {
    ysdf <- ydf %>% filter(cluster %in% ((i-1)*9 +1):(i*9))
    predsdf <- preddf %>% filter(cluster %in% ((i-1)*9 +1):(i*9))
    p <- ggplot(ysdf) +
      geom_line(aes(x = time, y = value, group = region), color = "gray", linewidth = 0.5) +
      geom_line(data = predsdf, aes(x = time, y = value, group = newregion),
                color = "red", linewidth = 0.6, linetype = 2) +
      facet_wrap(~ cluster, ncol = 3, scales = "free_y") +
      theme_bw() +
      theme(legend.position = "none", legend.text = element_text(size = 7)) +
      labs(title = "Relative Risk by cluster per region", color = "Cluster",
           y = "RR", x = "Time")
    print(p)
    # Print plot to the PDF device
  }
  dev.off()
}

#' Plot True and Estimated Cluster Maps
#'
#' This function generates a side-by-side visualization of true and estimated cluster maps based on spatial data.
#' The function uses `ggplot2` to create the maps and `ggpubr` for arranging the plots side by side. The output
#' is saved as a jpeg file.
#'
#' @param clust_res Numeric vector of estimated cluster assignments for each area in the map.
#' @param cluster_true Numeric vector of true cluster assignments for each area in the map.
#' @param map An object of class `sf` (simple features), representing the spatial framework for the cluster data.
#' @param filepath Character string, the path and name of the file where the plot will be saved.
#' @param title Character string, title of the plot
#'
#' @details The function creates two maps using the simple features data provided in `map`. The first map
#' visualizes the true cluster assignments using different colors, and the second map displays the estimated
#' cluster assignments. These maps are helpful for comparing the accuracy of cluster predictions against
#' ground truth in spatial clustering analyses.
#'
#' The plots are arranged side by side for easy comparison, with the left map showing true clusters and the
#' right map showing estimated clusters. Both maps are enhanced with `theme_bw()` for a clean background
#' and labeled appropriately.
#'
#' @return The function saves the plot to a file and returns an invisible copy of the combined plot object.
#' The primary output is the jpeg file specified by `filepath`.
#'
#' @examples
#' # Assuming `clust_res`, `cluster_true` are available and `map` is an sf object:
#' \dontrun{
#' plotClusterMap(clust_res, cluster_true, map, "path/to/save/cluster_maps.jpeg")
#' }
#' @export
plotClusterMap <- function(clust_res, cluster_true = NULL, map, filepath = file.path(path_im, 'fun_us_p1.pdf'), title = "Estimated clusters", fill = "Estimated Cluster",palette = NULL) {
  p2 <- ggplot(map) +
    geom_sf(aes(fill = factor(clust_res))) +  # Ensure clust_res is treated as a factor
    scale_fill_manual(values = palette, na.translate = F ) +  # Apply your custom color palette
    theme_bw() + scale_fill_npg(na.translate = F) +
    labs(title = title, fill = "Estimated Cluster")
  if(!is.null(cluster_true)){
    p1 <- ggplot(map) +
      geom_sf(aes(fill = cluster_true)) +  # Ensure clust_res is treated as a factor
      scale_fill_manual(values = palette) +  # Apply your custom color palette+
      theme_bw() +
      labs(title = "True clusters", fill = "True cluster")
    p <- ggpubr::ggarrange(p1, p2, ncol = 2)
  } else {p <- p2}

  jpeg(width = 960, height = 480, filepath)
  print(p)
  dev.off()
  print(p)
  return(p)
}

#' Plot Clustering Progress Over Iterations
#'
#' This function visualizes the evolution of cluster memberships over iterations.
#' It creates a heatmap using the `fields` package to show how clusters change across
#' different iterations, saving the plot as a jpeg file.
#'
#' @param clust_res Matrix with rows representing iterations and columns representing membership indices,
#'        showing the cluster assignment of each member at each iteration.
#' @param filepath Character string specifying the path and name of the jpeg file to save the plot.
#'
#' @details The function uses `fields::image.plot` to create a heatmap where each column represents a member
#' and each row represents an iteration. The color in the heatmap represents the cluster to which
#' a member belongs in a given iteration. This visual representation helps in understanding
#' the stability and convergence of the clustering algorithm over time.
#'
#' @return Invisibly returns NULL. The function's primary output is the generation of a jpeg file specified
#' by the `filepath`. The plot is saved directly to the file, and nothing is returned to the R environment.
#'
#' @examples
#' # Assuming clust_res is a matrix where rows are iterations and columns are memberships:
#' \dontrun{
#' plotClusterIter(clust_res, "path/to/save/cluster_iterations.jpeg")
#' }
#' @export
plotClusterIter <- function(cluster_out, filepath) {
  jpeg(filepath)
  print(fields::image.plot(1:ncol(cluster_out), 1:nrow(cluster_out), t(cluster_out), main = "Cluster Procedure", axes = TRUE, xlab = "Memembership", ylab = "Iterations"))
  dev.off()
  print(fields::image.plot(1:ncol(cluster_out), 1:nrow(cluster_out), t(cluster_out), main = "Cluster Procedure", axes = TRUE, xlab = "Memembership", ylab = "Iterations"))
}
#' Plot Log Marginal Likelihood Over Iterations
#'
#' This function visualizes the log marginal likelihood over iterations using a line plot.
#' It is designed to help evaluate the convergence and performance of statistical models
#' over iterative processes.
#'
#' @param log_mlike Numeric vector containing the log marginal likelihood values at each iteration.
#'
#' @details The function converts the log marginal likelihood data into a data frame and
#' then uses `ggplot2` to plot these values. The x-axis represents the iteration number,
#' and the y-axis represents the log marginal likelihood at that iteration. This plot can
#' be especially helpful in monitoring the convergence of Bayesian models or any iterative
#' statistical modeling that provides log marginal likelihood as an output.
#'
#' @return The function returns a ggplot object representing the line plot of the log marginal
#' likelihood over iterations. This object can be further modified with additional `ggplot2`
#' layers or printed directly to view the plot.
#'
#' @examples
#' # Assuming log_mlike is a vector of log marginal likelihood values:
#' \dontrun{
#' plot <- plotmlikeIter(restult$log_mlike)
#' print(plot)
#' }
#' @export
plotmlikeIter <- function(log_mlike,filepath) {
  p <- ggplot(data.frame(id = 1:length(log_mlike), log_mlike = log_mlike)) +
    geom_line(aes(id, log_mlike)) +
    labs(x = "iteration", y = "log marginal likelihood")
  jpeg(width = 960, height = 480, filepath)
  print(p)
  dev.off()
  print(p)
}

funInterval <- function(k,final_model, fun){
  r = final_model[[k]]
  r.samples =INLA::inla.posterior.sample(100, final_model[[k]])
  f1 = inla.posterior.sample.eval(fun, r.samples)
  return(f1)
}

#' Plot a Graph on Spatial Data
#'
#' Visualizes a graph over a spatial map using specified coordinates and graph object.
#' The function uses the ggraph package to create a custom layout for displaying the graph nodes and edges.
#'
#' @param map An sf object representing the spatial data.
#' @param coords Matrix of coordinates where each row corresponds to a node in the graph.
#' @param graph0 A graph object, typically an igraph object, representing the graph to be plotted.
#' @param title Character string specifying the title of the plot.
#' @param filepath Character string, the name of the file where the plot will be saved.
#'
#' @return A ggplot object representing the graph plotted over the spatial data.
#'
#' @examples
#' \dontrun{
#'   library(sf)
#'   library(igraph)
#'   library(ggraph)
#'   # Assuming 'map' is an sf object and 'coords' are extracted from 'map'
#'   # and 'graph0' is a previously constructed graph object
#'   plot <- plotGraph(map, coords, graph0, "Sample Graph")
#'   print(plot)
#' }
#'
#' @export
plotGraph <- function(map, coords, graph0, title, filepath){
  p = ggraph(graph0, layout = 'manual', x = coords[,1], y = coords[,2]) +
    geom_edge_link() +
    geom_node_point(color = 'red', size = 3) +
    geom_sf(data = map, inherit.aes = FALSE, fill = NA) +
    ggtitle(title) +
    theme_minimal()
  if(!is.null(filepath)){
    jpeg(width = 960, height = 480, filepath)
    print(p)
    dev.off()
  }
  return(p)
}

#' Plot a Minimum Spanning Tree on Spatial Data
#'
#' Visualizes a minimum spanning tree (MST) over a spatial map using specified coordinates and graph object.
#' Nodes are colored by cluster assignments. The function uses the ggraph package to create a custom layout
#' for displaying the MST with nodes and edges.
#'
#' @param map An sf object representing the spatial data.
#' @param coords Matrix of coordinates where each row corresponds to a node in the graph.
#' @param graph0 A graph object, typically an igraph object, that includes the MST to be plotted.
#' @param cluster A vector indicating the cluster assignment for each node in the graph.
#' @param title Character string specifying the title of the plot.
#' @param filepath Character string, the name of the file where the plot will be saved.
#'
#' @return A ggplot object representing the MST plotted over the spatial data with nodes colored by cluster.
#'
#' @examples
#' \dontrun{
#'   library(sf)
#'   library(igraph)
#'   library(ggraph)
#'   # Assuming 'map' is an sf object, 'coords' are coordinates, 'graph0' is a graph object,
#'   # and 'cluster' is a vector of cluster assignments
#'   plot <- plotMST(map, coords, graph0, cluster, "Sample MST")
#'   print(plot)
#' }
#'
#' @export
plotMST <- function(map, coords, graph0, cluster, title, filepath = NULL){
  V(graph0)$cluster <- as.factor(cluster)
  p = ggraph(graph0, layout = 'manual', x = coords[,1], y = coords[,2]) +
    geom_edge_link() +
    geom_node_point(aes(color = cluster), size = 3) +
    geom_sf(data = map, inherit.aes = FALSE, fill = NA) +
    ggtitle(title)+
    theme_minimal()+
    theme(legend.position="none")
  if(!is.null(filepath)){
    jpeg(width = 960, height = 480, filepath)
    print(p)
    dev.off()
  }
  return(p)
}

