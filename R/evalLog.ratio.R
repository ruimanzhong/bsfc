#' Evaluate Log-Likelihood Ratio for a Proposed Move
#'
#' This function calculates the log-likelihood ratio for a proposed "split" or "merge"
#' move within a clustering or partitioning algorithm using the Integrated Nested
#' Laplace Approximation (INLA) approach. It updates local likelihoods based on
#' whether a split or merge operation has been proposed and supports multiple family
#' distributions.
#'
#' @param move A character string specifying the type of move; either "split" or "merge".
#' @param log_like_vec A numeric vector of existing log likelihoods for all clusters.
#' @param control.move A list containing control parameters and information about the move,
#'        such as `clust_old`, `k`, `cluster_newid`, and `cluster_rm`.
#' @param Y Numeric matrix of data observations.
#' @param membership Integer vector indicating the cluster membership of each observation.
#' @param X Numeric matrix representing the model's design matrix or covariates.
#' @param N The number of trials or cases, applicable in certain statistical families like binomial.
#' @param family The family of distributions to use, e.g., "normal", "poisson", etc.
#' @param formula An object of class \code{\link[stats]{formula}} specifying the model to be used in INLA.
#' @param correction Logical, specifying whether correction calculations should be executed.
#' @param detailed Logical, indicating whether detailed output from INLA is required.
#'
#' @return A list containing the updated log likelihood vector and the log likelihood ratio (`llratio`).
#' @export
evalLogMLike.ratio <- function(move, log_like_vec, control.move, Y, membership, X = NULL, N = NULL, family = "normal", formula = Yk ~ 1 + Xk, correction = FALSE, detailed = FALSE, ...){

  # update local likelihoods for split move
  if(move == 'split'){
    log_like_vec_new <- log_like_vec;
    M1 <- evalLogMLike_each(control.move$clust_old, Y, membership, X, N, family, formula, correction, detailed, ...)
    M2 <- evalLogMLike_each(control.move$k, Y, membership, X, N, family, formula, correction, detailed, ...)
    log_like_vec_new[control.move$clust_old] <- M1
    log_like_vec_new[control.move$k] <- M2
    llratio <- M1 + M2 - log_like_vec[control.move$clust_old]
  }

  # update local likelihoods for merge move
  if(move == 'merge'){
    M1 <- evalLogMLike_each(control.move$cluster_newid, Y, membership, X, N, family, formula, correction, detailed, ...)
    log_like_vec_new <- log_like_vec[-control.move$cluster_rm]
    log_like_vec_new[control.move$cluster_newid] <- M1
    llratio <-  M1 - sum(log_like_vec[c(control.move$cluster_rm, control.move$cluster_comb)])
  }

  return(list(ratio = llratio, log_mlike_vec = log_like_vec_new))
}
