#' Detect Spatial Functional Clusters Based on Bayesian Spanning Tree
#'
#' This function implements a Reversible Jump Markov Chain Monte Carlo (RJMCMC) algorithm
#' for detecting spatial functional clusters based on a Bayesian analysis of spanning trees.
#' It handles various types of moves including birth, death, change, and hyperparameter updates
#' to explore the space of possible cluster configurations.
#'
#' @param data A list containing data matrices 'Y' and 'X', and 'N' for the number of trials or cases.
#' @param formula A formula object representing the model to be fitted.
#' @param graph0 Initial spanning tree used for the Bayesian model.
#' @param init_val List of initial values for parameters 'trees', 'beta', and 'cluster'.
#' @param hyperpar List containing hyperparameters used in the model.
#' @param MCMC Integer, number of MCMC iterations to perform.
#' @param BURNIN Integer, number of burn-in iterations to discard.
#' @param THIN Integer, thinning interval for recording the results.
#' @param path_save Character, the path where results should be saved.
#' @param seed Integer, seed value for random number generation, defaults to 1234.
#'
#' @return NULL The function primarily outputs results to a specified path and does not return anything.
#'
#' @examples
#' # Example setup (note: actual data and parameters need to be defined)
#' \dontrun{
#' data <- list(Y = matrix(rnorm(100), ncol=10), X = matrix(rnorm(100), ncol=10))
#' formula <- Y ~ X1 + X2
#' inla.extra <- list(correction = TRUE)
#' graph0 <- matrix(sample(0:1, 100, replace=TRUE), ncol=10)
#' init_val <- list(trees = graph0, beta = runif(10), cluster = sample(1:5, 10, replace=TRUE))
#' hyperpar <- list(c = 0.5)
#' path_save <- "path/to/save/results/"
#'
#' BayesClust(data,family, formula, graph0, init_val, hyperpar, MCMC, BURNIN, THIN, path_save, seed = 1234)
#'}
#' @export
BayesClust <- function(data, formula = Yk ~ 1 + Xk, family = "normal", graph0, init_val, hyperpar, correction = F, detailed = F, MCMC, BURNIN, THIN, path_save, seed = 1234, ...) {
  set.seed(seed)


  # initial values

  mstgraph = init_val[['trees']]
  beta = init_val[['beta']]
  cluster = init_val[['cluster']]
  k = max(cluster) # number of clusters
  Y = data$Y
  X = data$X
  N = data$N
  #For recording the acceptance ratio
  hyper_cnt = 0; birth_cnt=0; death_cnt=0; change_cnt=0;
  nt = ncol(Y); ns = nrow(Y);
  ################# RJMCMC ####################

  # hyper-parameter
  rhy = 0.05
  c = hyperpar$c



  ### initialize log likelihood vector
  p = max(cluster)
  llik_res = lapply(1:p, evalLogMLike_each_INLA, Y, cluster, X, N, family, formula, correction, FALSE)
  log_like_vec = unlist(llik_res)
  log_like = sum(log_like_vec)



  # whether an edge in graph0 is within a cluster or bewteen two clusters
  # n*p matrix
  edge_status = getEdgeStatus(cluster, mstgraph)

  ## MCMC results
  beta_out = list();
  mstgraph_out=list(); cluster_out = array(0, dim = c((MCMC-BURNIN)/THIN, ns))
  MST_out = list()
  loglike_out = numeric((MCMC-BURNIN)/THIN)

  #time_start = Sys.time()
  ## MCMC iteration
  for(iter in 1:MCMC) {

    if(k == 1) {rb = 0.95; rd = 0; rc = 0
    } else if(k == ns) {rb = 0; rd = 0.85; rc = 0.1
    } else {rb = 0.425; rd = 0.425; rc = 0.1}

    move = sample(4, 1, prob = c(rb, rd, rc, rhy))

    if(move == 1) { ## Birth move
      # split an existing cluster
      split_res = splitCluster(mstgraph, k, cluster)
      membership_new = split_res$cluster
      split_res$k=k+1

      # compute log-proposal ratio
      if(k == ns-1) {
        rd_new = 0.85
      } else {rd_new = 0.425}
      # if(k_m == n-1) {
      #   rd_new = 0.6
      # } else {rd_new = 0.3}
      log_P = log(rd_new) - log(rb)
      # compute log-likelihood ratio

      ind_k=which(membership_new==k+1); ns1=length(ind_k)
      Y_clust=Y[ind_k,]

      # compute log-prior ratio
      log_A = log(1-c)
      ## calculate loglikelihood ratio by only comparing local likelihood of the two clusters that changed.

      log_L_new=evalLogLike.ratio('split',log_like_vec, split_res, Y, membership_new, X, N, family, formula, correction, detailed = F, ...)
      log_L=log_L_new$ratio

      #acceptance probability
      acc_prob = min(0, log_A + log_P + log_L)
      acc_prob = exp(acc_prob)
      if(runif(1) < acc_prob){
        # accept
        cluster = membership_new
        k = k + 1

        log_like_vec=log_L_new$log_like_vec
        log_like = sum(log_like_vec);

        edge_status = getEdgeStatus(cluster, mstgraph)
        birth_cnt=birth_cnt+1

      }

    }

    if(move == 2) { ## Death move
      # merge two existing clusters (c1, c2) -> c2
      merge_res = mergeCluster(mstgraph, edge_status, cluster)
      membership_new = merge_res$cluster

      ##index of cluster that is removed
      cid_rm = merge_res$cluster_rm


      # # compute log-proposal ratio
      if(k == 2) {rb_new = 0.85
      }else {rb_new = 0.425}

      log_P = log(rb_new) - log(rd)


      ##For hetero variances
      ##index of the new merged cluster
      cid_comb = merge_res$cluster_comb
      cid_newid = merge_res$cluster_newid

      ind_k=which(membership_new == cid_newid); ns1=length(ind_k)
      Y_clust=Y[ind_k,]


      # compute log-prior ratio

      log_A = -log(1-c)

      # compute log-likelihood ratio

      log_L_new=evalLogLike.ratio('merge',log_like_vec, merge_res, Y, membership_new, X, N, family, formula, correction, detailed = F, ...)
      log_L=log_L_new$ratio



      #acceptance probability
      acc_prob = min(0, log_A + log_P + log_L)
      acc_prob = exp(acc_prob)
      if(runif(1) < acc_prob){
        # accept
        cluster = membership_new
        k = k - 1

        log_like_vec=log_L_new$log_like_vec
        log_like = sum(log_like_vec);
        # ind_gam=ind_gam_new;

        edge_status = getEdgeStatus(cluster, mstgraph)

        death_cnt=death_cnt+1

      }

    }


    if(move == 3) { ## change move
      # first perform death move: (c1, c2) -> c2
      merge_res = mergeCluster(mstgraph, edge_status, cluster)
      cid_rm = merge_res$cluster_rm
      membership_new = merge_res$cluster

      # ind_gam_new=ind_gam[-cid_rm,]

      ##For hetero variances

      cid_comb = merge_res$cluster_comb
      cid_newid = merge_res$cluster_newid

      ind_k=which(membership_new == cid_newid); ns1=length(ind_k)

      Y_clust=Y[ind_k,]

      log_prior_merge = 0

      k = k-1

      log_L_new_merge=evalLogLike.ratio('merge',log_like_vec, merge_res, Y, membership_new, X, N, family, formula, correction, detailed = F, ...)

      # then perform birth move
      split_res = splitCluster(mstgraph, k, merge_res$cluster);
      split_res$k=k+1;membership_new=split_res$cluster

      k = k+1

      # For hetero variances
      ind_k=which(membership_new==k); ns1=length(ind_k)
      Y_clust=Y[ind_k,]


      log_prior_split=0


      # compute log-likelihood ratio

      log_L_new=evalLogLike.ratio('split',log_L_new_merge$log_like_vec, split_res, Y, X, inla.extra, membership_new, formula, detailed  = F, correction)
      log_L = log_L_new$ratio + log_L_new_merge$ratio + log_prior_merge + log_prior_split




      # acceptance probability
      acc_prob = min(0, log_L)
      acc_prob = exp(acc_prob)
      if(runif(1) < acc_prob){
        # accept
        cluster = membership_new

        log_like_vec=log_L_new$log_like_vec
        log_like = sum(log_like_vec);
        # ind_gam=ind_gam_new

        edge_status = getEdgeStatus(cluster, mstgraph)

        change_cnt=change_cnt+1

      }

    }

    if(move == 4) { ## Hyper move
      hyper_cnt = hyper_cnt + 1

      #################################################################################

      # update MST

      edge_status_G = getEdgeStatus(cluster, graph0)
      mstgraph = proposeMST(graph0, edge_status_G)

      V(mstgraph)$vid=1:ns

      edge_status = getEdgeStatus(cluster, mstgraph)

      #######################################################################################

      ##Update log-likelihood value
      p = max(cluster)
      llik_res = lapply(1:p, evalLogMLike_each_INLA, Y, cluster, X, N, family, formula, correction, FALSE)
      log_like_vec = unlist(llik_res)
      log_like = sum(log_like_vec)
    }

    ##############################################################################

    ###store estimates
    if(iter %% 10 == 0) {
      cat('move', move, 'k', k, 'ncluster', length(unique(cluster)), '\n')
    }
    ## save result
    if(iter > BURNIN & (iter - BURNIN) %% THIN == 0) {

      MST_out[[(iter-BURNIN)/THIN]] = mstgraph
      cluster_out[(iter-BURNIN)/THIN, ] = cluster
      loglike_out[(iter-BURNIN)/THIN] = log_like
      if(iter %% 10 == 0) {
        cat('Iteration', iter, 'birth cnt', birth_cnt, 'death cnt', death_cnt, 'change cnt', change_cnt, 'hyper cnt', hyper_cnt, '\n')
        # beta_new = postSample_par(k, Y, X,  a0, b0, cluster, cl)
        save( MST_out, cluster_out, loglike_out, hyperpar, birth_cnt, death_cnt, change_cnt, file=path_save)
        print(cluster_out[(iter-BURNIN)/THIN, ])
      }
    }

  }

  ##Update log-likelihood value
  p=max(cluster)
  llik_res=lapply(1:p, evalLogLike_each_INLA, Y, X, inla.extra, cluster, formula, detailed = FALSE, correction)
  log_like_vec=unlist(llik_res)
  log_like=sum(log_like_vec)

}





