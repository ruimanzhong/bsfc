#' Plot Clustered Functions
#'
#' This function plots the results of a clustering analysis, displaying each cluster with its respective time series.
#' It uses pre-computed model results to add a line representing the mean fitted values for each cluster.
#' The function saves the plots to a specified file.
#'
#' @param Y Numeric matrix containing the data to be plotted, with each column representing a different time series.
#' @param E Numeric matrix of expected cases
#' @param nt Integer, the number of time points.
#' @param clust_res Numeric vector containing cluster assignments for each time series in Y.
#' @param final_model List of model objects for each cluster, each containing a summary of fitted values.
#' @param height plot height
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
plotClusterFun <- function(Y, N = NULL, link, nt,ns, clust_res,cluster_k, final_model,height = 8,ncol = 3,  linetype = 1, filepath) {
  # Define link functions
  link_function <- switch(link,
                          "exp" = exp,
                          "sigmoid" = function(x) 1 / (1 + exp(-x)),
                          "identity" = function(x) x,
                          stop("Invalid link function provided. Use 'exp', 'sigmoid', or 'identity'."))

  # Convert data to long format
  ydf <- setNames(as.data.frame(Y / E), 1:ns) |>
    mutate(time = nt) |>
    tidyr::pivot_longer(1:ns, names_to = "region", names_transform = as.numeric) |>
    mutate(cluster = clust_res[region])

  # Filter the data to plot only the selected cluster
  ydf <- ydf |> filter(cluster == cluster_k)

  # Apply the link function to the fitted values for the selected cluster
  preddf <- final_model[[cluster_k]]$summary.fitted.values$mean |>
    link_function() / exp(final_model[[cluster_k]]$summary.random[['id']]$mean) |>
    matrix(nrow = length(nt)) |>
    as.data.frame() |>
    setNames(1:ns) |>
    mutate(time = nt) |>
    tidyr::pivot_longer(1:ns, names_to = "newregion", names_transform = as.numeric) |>
    mutate(cluster = factor(cluster_k))

  # Define the custom labeller function for clusters
  custom_labeller <- function(cluster) {
    return(paste("Cluster", cluster))
  }

  # Plot the results
  p <- ggplot() +
    geom_line(data = ydf, aes(x = time, y = value, group = region, color = as.factor(region)), linetype = linetype) +
    geom_line(data = preddf, aes(x = time, y = value, group = newregion, color = as.factor(newregion)), linetype = linetype, alpha = 0.6) +
    labs(title = paste("Cluster", cluster_k), x = "Time", y = "Response") +
    theme_minimal() +
    facet_wrap(~ cluster, ncol = ncol, labeller = as_labeller(custom_labeller)) +
    theme(legend.position = "none")


  return(p)
}
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
plotClusterFun2 <- function(Y,E, nt,ns, clust_res, final_model,page,limit = max(clust_res),height = 8,ncol= 3,filepath) {
  p = max(clust_res)
  ydf = setNames(as.data.frame( Y / (E+0.001)), 1:ns) |>
    mutate(time = nt) |>
    tidyr::pivot_longer(1:ns, names_to = "region", names_transform = as.numeric) |>
    mutate(cluster = clust_res[region])
  ydf[ydf$value > 45,"value"] <- 0
  cluster2 = rep(1:p, table(clust_res))
  mean = purrr::map(1:p, function(x) final_model[[x]]$summary.linear.predictor$mean - final_model[[x]]$summary.random[['id']]$mean) %>%
    purrr::map(~ matrix(., nrow = length(nt))) %>%
    do.call(cbind, .) |>
    as.data.frame() |>
    setNames(1:ns) |>
    mutate(time =  nt) |>
    tidyr::pivot_longer(1:ns, names_to = "newregion", names_transform = as.numeric) |>
    mutate(cluster = cluster2[newregion])
  mean[is.infinite(mean$value),"value"] <- 0
  mean[mean$value > 100,"value"] <- 0
  plot_list = NULL
  # Define the custom labeller function
  custom_labeller <- function(cluster) {
    return(paste("Cluster", cluster))
  }
  for (i in 1:page) {
    ysdf <- ydf %>% filter(cluster %in% ((i-1)*15 +1):(i*15) & cluster <= limit)
    predsdf <- mean %>% filter(cluster %in% ((i-1)*15 +1):(i*15) & cluster <= limit)
    predsdf[predsdf$value > 3000, "value"] <- 0

    p <- ggplot(ysdf) +
      geom_line(aes(x = time, y = value, group = region), color = "gray", linewidth = 0.2, alpha = 0.8) +
      geom_line(data = predsdf, aes(x = time, y = exp(value), group = newregion),
                color = "red", linewidth = 0.6, linetype = 1) +
      # geom_ribbon(data = predsdf,aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
      facet_wrap(~ cluster, ncol = ncol, scales = "free_y", labeller = labeller(cluster = as_labeller(custom_labeller))) +
      theme_bw() +
      theme(legend.position = "none", legend.text = element_text(size = 7),
            axis.title = element_blank()) +
      labs( color = "Cluster",
           y = "Relative Risk")
    plot_list[[i]] = p
    print(p + theme_Publication())
  }
  dev.off()
  return(plot_list)
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
plotClusterMap <- function(clust_res, cluster_true = NULL, map, height = NULL, filepath = NULL, title = "Estimated Clusters", fill = "Estimated",palette = NULL) {
  # Plot for estimated clusters
  letters <- LETTERS
  p2 <- ggplot(map) +
    geom_sf(aes(fill = factor(clust_res)))  +  # Apply your custom color palette
    theme_bw() +
    theme( # Larger title font size
      axis.title = element_blank(),         # No axis titles
      axis.text = element_blank(),          # No axis text
      axis.ticks = element_blank(),         # No axis ticks
      panel.grid = element_blank(),
      panel.border = element_blank(),       # Remove panel border
      plot.margin = unit(c(0, 0, 0, 0), "cm") # No panel grid
    ) +
    labs(title = title, fill = fill)

  # Plot for true clusters (if provided)
  if (!is.null(cluster_true)) {
    p1 <- ggplot(map) +
      geom_sf(aes(fill = factor(cluster_true)))  +  # Apply your custom color palette
      theme_bw() +
      theme( # Larger title font size
        axis.title = element_blank(),         # No axis titles
        axis.text = element_blank(),          # No axis text
        axis.ticks = element_blank(),         # No axis ticks
        panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm")# No panel grid
      ) +
      labs(title = "  (A) True Clusters", fill = "True")

    # Arrange the two plots side by side
    p <- patchwork::wrap_plots(p1,p2)
  } else {
    p <- p2
  }
  if(!is.null(filepath)){
    ggsave(filename = filepath, plot = p, device = "pdf", width = 8.27, height = height)
  }

  return(p)
}
#' Plot True and Estimated Cluster Maps
#'
#' This function generates a side-by-side visualization of true and estimated cluster maps based on spatial data.
#' The function uses `ggplot2` to create the maps and `ggpubr` for arranging the plots side by side. The output
#' is saved as a jpeg file.
#'
#' @param cluster_res Numeric vector of estimated cluster assignments for each area in the map.
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
#' plotClusterMap2(clust_res,3, map, "path/to/save/cluster_maps.jpeg")
#' }
#' @export
#'
plotClusterMap2 <- function(cluster_res, map,height = 8, filepath = file.path(path_im, 'fun_us_p1.pdf'), title = "Estimated clusters", fill = "Cluster", palette = NULL) {
  p = max(cluster_res)
  plot_list = list()
  num_labels <- length(unique(cluster_res))

  for (i in 1:ceiling(p/9)) {
    start_label <- (i-1) * 9 + 1
    end_label <- min(i *9, num_labels)

    current_labels <- seq(start_label,end_label, by = 1)
    clust_res = cluster_res
    clust_res[!clust_res %in% current_labels] <- NA
    map$clust_res = clust_res

    plot_list[[i]] <- ggplot(map) +
      geom_sf(aes(fill = as.factor(clust_res))) +  # Ensure clust_res is treated as a factor
      scale_fill_npg(na.translate = F) +  # Apply your custom color palette
      theme_bw() +
      theme(
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title = element_blank(),         # No axis titles
        axis.text = element_blank(),          # No axis text
        axis.ticks = element_blank(),         # No axis ticks
        panel.grid = element_blank()          # No panel grid
      ) + labs(fill = "Cluster")
  }
  # Combine the plots into a single layout
  combined_plot <- patchwork::wrap_plots(plot_list, ncol = 2)

  # Save the combined plot as a PDF using ggsave
  return(combined_plot)
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
#'
plotClusterIter <- function(cluster_out, height, filepath) {
  # Finding the highest label to set the level
  level = max(cluster_out)
  melted_data <- reshape2::melt(cluster_out)
  colnames(melted_data) <- c("Iteration", "Membership", "Cluster")

  # Ensuring that labels are treated as factors and checking the levels
  melted_data$Cluster <- factor(melted_data$Cluster, levels = 1:level)

  # Creating a color palette
  mycolor <- Polychrome::createPalette(level, c("#00ff00", "#ff0040", "#0000ff"))
  names(mycolor) <- as.character(1:level)  # Ensure color names match factor levels

  # Assigning the manual color scale
  colscl = scale_fill_manual(name = "Cluster", values = mycolor)

  # Creating the plot
  p <- ggplot(melted_data, aes(x = Membership, y = Iteration, fill = Cluster)) +
    geom_tile() +
    colscl +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 20),
      axis.title = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20),
          axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15)) +
    labs(title = "Cluster Procedure", x = "Regions", y = "Iterations")
  # Saving the plot
  ggsave(filename = filepath, plot = p, device = "jpeg", width = 8.27, height = height, units = "in", dpi = 300)

  return(p)
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
plotmlikeIter <- function(log_mlike,height, filepath) {
  p <- ggplot(data.frame(id = 1:length(log_mlike), log_mlike = log_mlike)) +
    geom_line(aes(id, log_mlike)) +
    theme_minimal() +     theme(
      axis.title = element_text(size = 18),         # No axis titles
      axis.text = element_text(size = 18),          # No axis text      # No panel grid
    ) +
    labs(x = "iteration", y = "log marginal likelihood")
  ggsave(filename = filepath, plot = p, device = "jpeg", width = 8.27, height = height, units = "in", dpi = 300)
  return(p)
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
plotGraph <- function(map, coords, graph0, title, filepath = NULL){
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
plotMST <- function(map, coords, graph0, cluster, title, height = NULL, filepath = NULL){
  V(graph0)$cluster <- as.factor(cluster)
  p = ggraph(graph0, layout = 'manual', x = coords[,1], y = coords[,2]) +
    geom_edge_link() +
    geom_node_point(aes(color = cluster), size = 3) +
    geom_sf(data = map, inherit.aes = FALSE, fill = NA) +
    ggtitle(title)+  theme_minimal()+
    theme(legend.position="none") + labs(x = NULL, y = NULL)
  if(!is.null(filepath)){
    ggsave(filename = filepath, plot = p, device = "pdf", width = 8.27, height = height)
  }
  return(p)
}

