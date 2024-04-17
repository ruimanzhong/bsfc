#' Compute Likelihood Ratio for a Proposed Move
#'
#' This function calculates the likelihood ratio for a proposed "split" or "merge" move
#' within a clustering or partitioning algorithm, using the Integrated Nested Laplace
#' Approximation (INLA). It is designed to work with spatial or model-based clustering
#' where likelihood needs to be recomputed for different configurations.
#'
#' @param move A character string specifying the type of move; either "split" or "merge".
#' @param log_like_vec A numeric vector of existing log likelihoods for all clusters.
#' @param control.move A list containing control parameters and information about the move,
#'        such as `clust_old`, `k`, `cluster_newid`, and `cluster_rm`.
#' @param Y Numeric matrix of data observations.
#' @param X Numeric matrix representing the model's design matrix or covariates.
#' @param inla.extra Additional parameters or data needed for running the INLA model.
#' @param membership Integer vector indicating the cluster membership of each observation.
#' @param formula An object of class \code{\link[stats]{formula}} specifying the model to be used in INLA.
#' @param detailed Logical, indicating whether detailed output from INLA is required.
#' @param correction Logical, specifying whether correction calculations should be executed.
#'
#' @return A list containing updated log likelihood vector and the log likelihood ratio (`llratio`).
#' @examples
#' # Assuming appropriate data and control setup
#' @export
evalLogLike.ratio=function(move,log_like_vec, control.move,Y, X, inla.extra, membership, formula, detailed, correction){


  ### update local likelihoods for split move
  if(move=='split'){
    log_like_vec_new=log_like_vec;
    M1 <- evalLogLike_each_INLA (control.move$clust_old, Y, X, inla.extra, membership, formula, detailed, correction)  ## new clust 1
    M2 <- evalLogLike_each_INLA (control.move$k,  Y, X, inla.extra, membership, formula, detailed, correction)  ## new clust 2
    log_like_vec_new[control.move$clust_old]=M1
    log_like_vec_new[control.move$k]=M2
    llratio=M1+M2-log_like_vec[control.move$clust_old]
  }

  ### update local likelihoods for merge move
  if(move=='merge'){
    M1 <- evalLogLike_each_INLA (control.move$cluster_newid, Y, X, inla.extra, membership, formula, detailed, correction)
    log_like_vec_new=log_like_vec[-control.move$cluster_rm]
    log_like_vec_new[control.move$cluster_newid]=M1
    llratio= M1-sum(log_like_vec[c(control.move$cluster_rm,control.move$cluster_comb)])
  }

  return(list(ratio=llratio,log_like_vec=log_like_vec_new))
}
