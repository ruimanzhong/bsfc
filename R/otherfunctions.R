#' Print method for sfclust objects
#'
#' This function prints details of an sfclust object, including the formula, counts, hyperparameter q, and the final log marginal likelihood.
#'
#' @param x An object of class 'sfclust'.
#' @param ... Additional arguments (currently unused).
#' @export
print.sfclust <- function(x, ...) {
  
  # Ensure that x is an sfclust object
  if (!inherits(x, "sfclust")) {
    stop("The object must be of class 'sfclust'.")
  }
  
  # Print the formula
  cat("Formula:\n")
  print(x$formula)
  
  # Print the counts (births, deaths, changes, hypers)
  cat("\nCounts:\n")
  print(x$counts)
  
  # Print the hyperparameter q
  cat("\nHyperparameter q:\n")
  print(x$q)
  
  # Print the final log marginal likelihood
  cat("\nFinal Log Marginal Likelihood:\n")
  print(x$log_mlike[length(x$log_mlike)])
  
  invisible(x)  # Return the object invisibly
}


#' Summary method for sfclust objects
#'
#' This function summarizes the cluster assignments from the a-th row of the cluster matrix.
#'
#' @param x An object of class 'sfclust'.
#' @param a An integer specifying the row of the cluster matrix to summarize (default is the last row).
#' @param ... Additional arguments (currently unused).
#' @export
summary.sfclust <- function(x, a = nrow(x$cluster), ...) {
  
  # Ensure that x is an sfclust object
  if (!inherits(x, "sfclust")) {
    stop("The object must be of class 'sfclust'.")
  }
  
  # Ensure that a is valid (within the range of available rows in the cluster matrix)
  if (a < 1 || a > nrow(x$cluster)) {
    stop("The value of a must be between 1 and the number of rows in the cluster matrix.")
  }
  
  # Extract the a-th row of the cluster matrix
  cluster_assignments <- x$cluster[a, ]
  
  # Generate a table summarizing the cluster assignments
  cluster_summary <- table(cluster_assignments)
  
  # Print the summary table
  cat("Summary of Cluster Assignments (Row", a, "):\n")
  print(cluster_summary)
  
  invisible(cluster_summary)  # Return the summary invisibly
}


#' Plot function for sfclust objects with a conditional legend and log marginal likelihood
#'
#' This function plots the map with estimated clusters from an sfclust object and adds a legend if the number of clusters is less than 10. 
#' It also plots the log marginal likelihood.
#'
#' @param map An sf object provided by the user.
#' @param output An sfclust object containing the clustering results.
#' @param k The row of the cluster matrix to use for plotting (default is the last row).
#' @param title A title for the plot (default is "Estimated Clusters").
#' @export
plot.sfclust <- function(output,map,  k = nrow(output$cluster), title = "Estimated Clusters") {
  
  # Ensure that output is an sfclust object
  if (!inherits(output, "sfclust")) {
    stop("The output must be an object of class 'sfclust'.")
  }
  
  # Ensure that map is an sf object
  if (!inherits(map, "sf")) {
    stop("The map must be an sf object.")
  }
  
  # Ensure that k is valid (within the range of available rows in the cluster matrix)
  if (k < 1 || k > nrow(output$cluster)) {
    stop("The value of k must be between 1 and the number of rows in the cluster matrix.")
  }
  # Extract the k-th row of the cluster matrix
  cluster_estimated <- output$cluster[k, ]
  
  # Get the unique clusters
  unique_clusters <- sort(unique(cluster_estimated))
  
  ## First Plot: The map with the estimated clusters
  # Adjust the margins to prevent the "figure margins too large" error
  par(mar = c(4, 4, 2, 2))
  
  # Plot the map with the estimated clusters
  plot(map$geometry, col = factor(cluster_estimated), main = title, 
       xlab = "Longitude", ylab = "Latitude", border = "black")
  
  # Add a legend if the number of clusters is less than 10
  if (length(unique_clusters) < 10) {
    legend("topright", legend = paste("Cluster", unique_clusters),
           fill = unique(factor(cluster_estimated)), border = "black",
           title = "Legend", cex = 0.8, bty = "n")
  }
  
  ## Second Plot: The log marginal likelihood
  par(mar = c(4, 4, 2, 2))  # Adjust margins again
  
  # Plot the log marginal likelihood
  plot(output$log_mlike, type = "l", main = "Log Marginal Likelihood",
       xlab = "Iteration", ylab = "Log Marginal Likelihood", col = "blue", lwd = 2)
  
  invisible(NULL)  # Return nothing, just plot
}

#' Plot function for sfclustm objects
#'
#' This function plots the fitted values (mean) from the a-th model in an sfclustm object.
#'
#' @param x An sfclustm object containing the models for each cluster.
#' @param a The index of the model to plot (default is 1).
#' @param title A title for the plot (default is "Fitted Values").
#' @export
plot.sfclustm <- function(x, a = 1, title = "Fitted Values") {
  
  # Ensure that x is an sfclustm object
  if (!inherits(x, "sfclustm")) {
    stop("The object must be of class 'sfclustm'.")
  }
  
  # Ensure that a is valid (within the range of available models)
  if (a < 1 || a > length(x)) {
    stop("The value of a must be between 1 and the number of models in the sfclustm object.")
  }
  
  # Extract the fitted values (mean) from the a-th model
  fitted_values <- x[[a]]$summary.fitted.values$mean
  
  # Ensure that fitted_values exists
  if (is.null(fitted_values)) {
    stop("Fitted values not found in the selected model.")
  }
  
  # Plot the fitted values
  plot(fitted_values, type = "l", col = "blue", lwd = 2, 
       xlab = "Index", ylab = "Fitted Mean", 
       main = paste(title, "- Model", a))
  
  invisible(NULL)  # Return nothing, just plot
}


#' Fit function for sfclust objects with cluster sorting
#'
#' This function applies the log_mlik_each function to fit models for each cluster in sfclust.
#' The clusters are sorted in descending order based on the number of members in each cluster.
#'
#' @param output An sfclust object containing the clustering results.
#' @param k The row of the cluster matrix to use (default is the last row).
#' @param Y The response variable.
#' @param X Covariates or predictor variables.
#' @param N Number of observations or other input (optional).
#' @param formula A model formula.
#' @param family A family object specifying the distribution (default is "normal").
#' @param correction Logical, indicating whether to apply corrections (default is FALSE).
#' @param detailed Logical, indicating whether detailed results should be returned (default is TRUE).
#' @param ... Additional arguments passed to log_mlik_each.
#' @return A list of models, each of class 'sfclustm', fit for each sorted cluster.
#' @export
fit.sfclust <- function(output, k = nrow(output$cluster), Y, X = NULL, N = NULL, formula, family = "normal", correction = FALSE, detailed = TRUE, ...) {
  
  # Ensure that output is an sfclust object
  if (!inherits(output, "sfclust")) {
    stop("The output must be an object of class 'sfclust'.")
  }
  
  # Ensure that k is valid (within the range of available rows in the cluster matrix)
  if (k < 1 || k > nrow(output$cluster)) {
    stop("The value of k must be between 1 and the number of rows in the cluster matrix.")
  }
  
  # Extract the cluster assignments from the k-th row
  cluster <- output$cluster[k, ]
  
  # Determine the number of clusters (p) based on the max of the cluster vector
  p <- max(cluster)
  
  # Count the number of members in each cluster
  cluster_sizes <- table(cluster)
  
  # Sort the clusters by size (descending order)
  sorted_clusters <- as.numeric(names(sort(cluster_sizes, decreasing = TRUE)))
  
  # Apply log_mlik_each to each cluster (based on the sorted clusters)
  models <- lapply(sorted_clusters, function(c) {
    log_mlik_each(Y, cluster == c, X, N, formula, family, correction = correction, detailed = detailed, ...)
  })
  
  # Assign the class 'sfclustm' to the result
  class(models) <- "sfclustm"
  
  return(models)
}

#' Print method for sfclustm objects
#'
#' This function prints the summary of each model within an sfclustm object.
#'
#' @param x An object of class 'sfclustm'.
#' @param ... Additional arguments (currently unused).
#' @export
print.sfclustm <- function(x, ...) {
  cat("sfclustm object containing models for each cluster:\n")
  for (i in seq_along(x)) {
    cat("\nModel for cluster", i, ":\n")
    print(summary(x[[i]]))
  }
}

#' Summary method for sfclustm objects
#'
#' This function provides a summary of each model within an sfclustm object.
#'
#' @param x An object of class 'sfclustm'.
#' @param ... Additional arguments (currently unused).
#' @return A list of summaries for each model.
#' @export
summary.sfclustm <- function(x, ...) {
  # Ensure that x is an sfclustm object
  if (!inherits(x, "sfclustm")) {
    stop("The object must be of class 'sfclustm'.")
  }
  
  # Initialize a list to store the summaries
  summaries <- lapply(seq_along(x), function(i) {
    cat("\nSummary of Model for Cluster", i, ":\n")
    summary(x[[i]])  # Call the summary function for each model in the sfclustm object
  })
  
  invisible(summaries)  # Return the summaries invisibly
}
