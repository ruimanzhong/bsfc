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
#' @param filename Character string, the name of the file where the plot will be saved.
#'
#' @details The function generates a PNG file with a series of plots arranged in a grid. Each plot corresponds
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
#' plotClusterFun(Y, 50, cluster_results, models, "path/to/save/cluster_plot.png")
#' }
#' @export
plotClusterFun <- function(Y, nt, clust_res, final_model, a, b, filename) {
  time <- 1:nt
  p <- max(clust_res)
  col_n <- gray.colors(80, alpha = 0.5)
  ymin <- min(Y)
  ymax <- max(Y)
  par(mfrow = c(a, b))
  time <- seq(0, 1, length.out = nt)
  for (k in 1:p) {
    IND <- which(clust_res == k)
    plot(time, Y[, IND[1]], "l",
      col = col_n[40], xlab = "t", ylab = "f(t)",
      main = sprintf("Cluster %d", k), ylim = c(ymin, ymax)
    )
    m <- t(matrix(final_model[[k]]$summary.fitted.values$mean, nrow = nt, byrow = FALSE))
    for (i in 2:length(IND)) {
      lines(time, Y[, IND[i]], col = col_n[40])
    }
    lines(time, colMeans(m), col = "red")
  }
  png(width = 960, height = 480, filename)
  par(mfrow = c(3, p / 3))
  time <- seq(0, 1, length.out = nt)
  for (k in 1:p) {
    IND <- which(clust_res == k)
    plot(time, Y[, IND[1]], "l",
      col = col_n[40], xlab = "t", ylab = "f(t)",
      main = sprintf("Cluster %d", k), ylim = c(ymin, ymax)
    )
    m <- t(matrix(final_model[[k]]$summary.fitted.values$mean, nrow = nt, byrow = FALSE))
    for (i in 2:length(IND)) {
      lines(time, Y[, IND[i]], col = col_n[40])
    }
    lines(time, colMeans(m), col = "red")
  }
  dev.off()
}

#' Plot True and Estimated Cluster Maps
#'
#' This function generates a side-by-side visualization of true and estimated cluster maps based on spatial data.
#' The function uses `ggplot2` to create the maps and `ggpubr` for arranging the plots side by side. The output
#' is saved as a PNG file.
#'
#' @param clust_res Numeric vector of estimated cluster assignments for each area in the map.
#' @param cluster_true Numeric vector of true cluster assignments for each area in the map.
#' @param map An object of class `sf` (simple features), representing the spatial framework for the cluster data.
#' @param filename Character string, the path and name of the file where the plot will be saved.
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
#' The primary output is the PNG file specified by `filename`.
#'
#' @examples
#' # Assuming `clust_res`, `cluster_true` are available and `map` is an sf object:
#' \dontrun{
#' plotClusterMap(clust_res, cluster_true, map, "path/to/save/cluster_maps.png")
#' }
#' @export
plotClusterMap <- function(clust_res, cluster_true, map, filename) {
  p1 <- ggplot(map) +
    geom_sf(aes(fill = cluster_true)) +
    theme_bw() +
    labs(title = "True clusters", fill = "True cluster")
  p2 <- ggplot(map) +
    geom_sf(aes(fill = clust_res)) +
    theme_bw() +
    labs(title = "Estimated clusters", fill = "Estimated Cluster")
  p <- ggpubr::ggarrange(p1, p2, ncol = 2)
  png(width = 960, height = 480, filename)
  print(p)
  dev.off()
  print(p)
}

#' Plot Clustering Progress Over Iterations
#'
#' This function visualizes the evolution of cluster memberships over iterations.
#' It creates a heatmap using the `fields` package to show how clusters change across
#' different iterations, saving the plot as a PNG file.
#'
#' @param clust_res Matrix with rows representing iterations and columns representing membership indices,
#'        showing the cluster assignment of each member at each iteration.
#' @param filename Character string specifying the path and name of the PNG file to save the plot.
#'
#' @details The function uses `fields::image.plot` to create a heatmap where each column represents a member
#' and each row represents an iteration. The color in the heatmap represents the cluster to which
#' a member belongs in a given iteration. This visual representation helps in understanding
#' the stability and convergence of the clustering algorithm over time.
#'
#' @return Invisibly returns NULL. The function's primary output is the generation of a PNG file specified
#' by the `filename`. The plot is saved directly to the file, and nothing is returned to the R environment.
#'
#' @examples
#' # Assuming clust_res is a matrix where rows are iterations and columns are memberships:
#' \dontrun{
#' plotClusterIter(clust_res, "path/to/save/cluster_iterations.png")
#' }
#' @export
plotClusterIter <- function(clust_res, filename) {
  png(filename)
  print(fields::image.plot(t(clust_res), main = "Cluster Procedure", axes = TRUE, xlab = "Memembership", ylab = "Iterations"))
  dev.off()
  print(fields::image.plot(t(clust_res), main = "Cluster Procedure", axes = TRUE, xlab = "Memembership", ylab = "Iterations"))
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
plotmlikeIter <- function(mod) {
  data.frame(id = 1:length(mod$log_mlike), log_mlike = mod$log_mlike) |>
    ggplot() +
    geom_line(aes(id, log_mlike)) +
    labs(x = "iteration", y = "log marginal likelihood")
}
