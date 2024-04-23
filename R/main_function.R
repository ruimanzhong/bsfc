#' Detect Spatial Functional Clusters Based on Bayesian Spanning Tree
#'
#' This function implements a Reversible Jump Markov Chain Monte Carlo (RJMCMC) algorithm
#' for detecting spatial functional clusters based on a Bayesian analysis of spanning trees.
#' It handles various types of moves including birth, death, change, and hyperparameter updates
#' to explore the space of possible cluster configurations.
#'
#' @param data A list containing data matrices 'Y' and 'X', and 'N' for the number of trials or cases.
#' @param formula A formula object representing the model to be fitted.
#' @param graph Initial spanning tree used for the Bayesian model.
#' @param init_val List of initial values for parameters 'trees', 'beta', and 'cluster'.
#' @param hyperpar List containing hyperparameters used in the model.
#' @param MCMC Integer, number of MCMC iterations to perform.
#' @param burnin Integer, number of burn-in iterations to discard.
#' @param thin Integer, thinning interval for recording the results.
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
#' graph <- matrix(sample(0:1, 100, replace=TRUE), ncol=10)
#' init_val <- list(trees = graph, beta = runif(10), cluster = sample(1:5, 10, replace=TRUE))
#' hyperpar <- list(c = 0.5)
#' path_save <- "path/to/save/results/"
#'
#' BayesClust(data,family, formula, graph, init_val, hyperpar, MCMC, burnin, THIN, path_save, seed = 1234)
#'}
#' @export
spfc <- function(Y, graph, X = NULL, N = NULL, formula = Yk ~ 1 + Xk, family = "normal", initial, hyperpar,
                 correction = FALSE, niter = 100, burnin = 0, thin = 1, path_save, ...) {

  ## Initial setup

  # dimensions
  ns <- nrow(Y)

  # initial values
  mstgraph <- initial[['trees']]
  cluster <- initial[['cluster']]
  k = max(cluster)

  # hyperparameters
  rhy = 0.05
  c = hyperpar$c

  # movement counts
  hyper_cnt = 0; birth_cnt = 0; death_cnt = 0; change_cnt = 0

  ## Initialize

  # initialize log likelihood vector
  log_mlike_vec = log_mlik_all(Y, cluster, X, N, formula, family, correction)
  log_mlike = sum(log_mlike_vec)

  # whether an edge in graph is within a cluster or bewteen two clusters
  # n*p matrix
  edge_status = getEdgeStatus(cluster, mstgraph)

  ## Prepare output

  cluster_out = array(0, dim = c((niter-burnin)/thin, ns))
  MST_out = list()
  log_mlike_out = numeric((niter-burnin)/thin)

  ## MCMC sampling

  for(iter in 1:niter) {

    if(k == 1) {rb = 0.95; rd = 0; rc = 0
    } else if(k == ns) {rb = 0; rd = 0.85; rc = 0.1
    } else {rb = 0.425; rd = 0.425; rc = 0.1}

    move = sample(4, 1, prob = c(rb, rd, rc, rhy))

    if(move == 1) { ## Birth move

      ## Propose a split movement
      split_res = splitCluster(mstgraph, k, cluster)
      split_res$k = k+1
      membership_new = split_res$cluster

      ## Compute log-ratios for probability of acceptance

      # movement ratio
      if(k == ns-1) {
        rd_new = 0.85
      } else {rd_new = 0.425}
      # if(k_m == n-1) {
      #   rd_new = 0.6
      # } else {rd_new = 0.3}
      log_P = log(rd_new) - log(rb)

      # prior ratio
      log_A = log(1-c)

      # marginal likelihood ratio
      log_L_new = log_mlik_ratio('split', log_mlike_vec, split_res, Y, X, N, formula, family, correction, FALSE, ...)
      log_L = log_L_new$ratio

      ## Accept with probability
      acc_prob = min(0, log_A + log_P + log_L)
      acc_prob = exp(acc_prob)
      if(runif(1) < acc_prob){
        cluster = membership_new
        k = k + 1

        log_mlike_vec = log_L_new$log_mlike_vec
        log_mlike = sum(log_mlike_vec)

        edge_status = getEdgeStatus(cluster, mstgraph)
        birth_cnt = birth_cnt + 1
      }
    }

    if(move == 2) { ## Death move

      ## Propose a merge movement (c1, c2) -> c2
      merge_res = mergeCluster(mstgraph, edge_status, cluster)
      membership_new = merge_res$cluster
      cid_rm = merge_res$cluster_rm

      ## Compute log-ratios for probability of acceptance

      # movement ratio
      if(k == 2) {rb_new = 0.85
      }else {rb_new = 0.425}
      log_P = log(rb_new) - log(rd)

      # prior ratio
      log_A = -log(1-c)

      # marginal likelihood ratio
      log_L_new = log_mlik_ratio('merge', log_mlike_vec, merge_res, Y, X, N, formula, family, correction, FALSE, ...)
      log_L = log_L_new$ratio

      ## Accept with probability
      acc_prob = min(0, log_A + log_P + log_L)
      acc_prob = exp(acc_prob)
      if(runif(1) < acc_prob){
        cluster = membership_new
        k = k - 1

        log_mlike_vec = log_L_new$log_mlike_vec
        log_mlike = sum(log_mlike_vec);

        edge_status = getEdgeStatus(cluster, mstgraph)
        death_cnt=death_cnt+1
      }
    }


    if(move == 3) { ## change move

      ## First: Propose a merge movement (c1, c2) -> c2
      merge_res = mergeCluster(mstgraph, edge_status, cluster)
      membership_new = merge_res$cluster
      cid_rm = merge_res$cluster_rm
      k = k-1

      log_L_new_merge = log_mlik_ratio('merge', log_mlike_vec, merge_res, Y, X, N, formula, family, correction, FALSE, ...)

      # Second: Propose a split movement
      split_res = splitCluster(mstgraph, k, merge_res$cluster);
      split_res$k = k+1
      membership_new = split_res$cluster
      k = k+1

      log_L_new = log_mlik_ratio('split', log_L_new_merge$log_mlike_vec, split_res, Y, X, N,
        formula, family, correction, FALSE, ...)
      log_L = log_L_new$ratio + log_L_new_merge$ratio

      ## Accept with probability
      acc_prob = min(0, log_L)
      acc_prob = exp(acc_prob)
      if(runif(1) < acc_prob){
        cluster = membership_new

        log_mlike_vec = log_L_new$log_mlike_vec
        log_mlike = sum(log_mlike_vec)

        edge_status = getEdgeStatus(cluster, mstgraph)
        change_cnt=change_cnt+1
      }
    }

    if(move == 4) { ## Hyper move

      hyper_cnt = hyper_cnt + 1

      ## Update MST

      edge_status_G = getEdgeStatus(cluster, graph)
      mstgraph = proposeMST(graph, edge_status_G)
      V(mstgraph)$vid = 1:ns
      edge_status = getEdgeStatus(cluster, mstgraph)
    }

    ## Store estimates

    if(iter %% 10 == 0) {
      cat('move', move, 'k', k, 'ncluster', length(unique(cluster)), '\n')
    }

    if(iter > burnin & (iter - burnin) %% thin == 0) {

      MST_out[[(iter-burnin)/thin]] = mstgraph
      cluster_out[(iter-burnin)/thin, ] = cluster
      log_mlike_out[(iter-burnin)/thin] = log_mlike
      if(iter %% 10 == 0) {
        cat('Iteration', iter, 'birth cnt', birth_cnt, 'death cnt', death_cnt, 'change cnt', change_cnt, 'hyper cnt', hyper_cnt, '\n')
        # beta_new = postSample_par(k, Y, X,  a0, b0, cluster, cl)
        save(MST_out, cluster_out, log_mlike_out, hyperpar, birth_cnt, death_cnt, change_cnt, file = path_save)
        print(cluster_out[(iter-burnin)/thin, ])
      }
    }

  }
}





