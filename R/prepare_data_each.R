#' Prepare data for a specific cluster without dplyr
#'
#' This function extracts data corresponding to a specific cluster from a given \code{sftime} object and reshapes it into a long format, using only base R.
#'
#' @param k An integer representing the cluster index (e.g., 1, 2, etc.). This identifies the regions that belong to the cluster.
#' @param data An \code{sftime} object containing spatial and temporal data, including the response and covariates.
#' @param membership A numeric vector where each element corresponds to the region membership for each spatial unit. The length of this vector should be equal to the number of regions.
#' @param formula A model formula specifying the response variable and the covariates (e.g., \code{case ~ temperature + precipitation + f(idt, model = "iid")}).
#'
#' @return A long-format data frame for the regions in cluster \code{k}, with additional columns for:
#' \item{id}{Unique identifier for each observation.}
#' \item{idt}{Time index for each observation.}
#' \item{ids}{Region index for each observation.}
#'
#' @details The function identifies the regions belonging to cluster \code{k} based on the \code{membership} vector. It then filters the \code{sftime} object for those regions and reshapes the data into a long format. The response variable and covariates are dynamically extracted from the \code{formula}. 
#'
#' @examples
#' # Example usage
#' # Assuming 'data' is an sftime object, 'membership' is a vector of region memberships, 
#' # and 'formula' specifies the model with response and covariates:
#' prepared_data <- prepare_data_each(k = 1, data = data, membership = membership, formula = formula)
#'
#' @export
prepare_data_each <- function(k, data, membership, formula) {
  # Identify the regions that belong to cluster 'k'
  ind <- which(membership == k)
  nk <- length(ind)  # Number of regions in this cluster
  
  # Extract the number of time points (ntime)
  ntime <- nrow(data) / length(membership)  # Total rows divided by the number of regions
  
  # Extract variable names from the formula
  response_var <- all.vars(formula[[2]])  # Extract the response variable
  predictor_vars <- all.vars(formula[[3]])  # Extract the predictor variables
  
  # Filter data for the relevant regions (ind)
  rows_to_keep <- rep(1:length(membership), each = ntime) %in% ind
  data_filtered <- data[rows_to_keep, ]
  
  # Create id, idt, ids columns
  id <- seq_len(nk * ntime)
  idt <- rep(seq_len(ntime), times = nk)
  ids <- rep(seq_len(nk), each = ntime)
  
  # Add the new columns to the filtered data
  data_filtered$id <- id
  data_filtered$idt <- idt
  data_filtered$ids <- ids
  
  return(data_filtered)
}
