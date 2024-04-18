#' Evaluate Log-Likelihood for Each Group in INLA
#'
#' This function computes the log-likelihood for each group based on the membership
#' indicator using the Integrated Nested Laplace Approximation (INLA) approach.
#' It allows flexible definition of the model formula and can be customized
#' to include detailed output or just the marginal likelihood.
#'
#' @param k integer, specifying the group index to compute.
#' @param Y A ns by nt numeric matrix, containing the response matrix.
#' @param X A numeric matrix, the inverse discrete wavelet transform matrix used in the model.
#' @param inla.extra A list of within cluster model, the model type, and related parameters
#' @param membership integer vector, indicating the membership of each observation in `Y`.
#' @param formula an object of class \code{\link[stats]{formula}}, defining the model to be fit.
#' @param detailed logical, indicating whether to return a detailed output including
#'        full INLA object (`TRUE`) or just the computed marginal likelihood (`FALSE`).
#' @param correction logical, if a correction should be applied when computing marginal loglikelihood, for more details, please see inla.doc("rw1")
#' @return If `detailed` is TRUE, returns the full INLA object; otherwise, returns
#'         a numeric value of the marginal likelihood.
#' @examples
#' \dontrun{# Define data matrices and parameters
#' Y <- matrix(rpois(100, lambda = 10), ncol=10)
#' population <- matrix(rpois(100, lambda = 5), ncol=10)
#' inla.extra <- list(family = "poisson", N = population)
#' membership <- sample(1:2, 10, replace=TRUE)
#' X <- seq(0,1, length.out = 10)
#' a0 <- 0.1; b0 <- 0.1
#'
#' # Define a formula, the name of the iid random effec should be idx
#' formula <- Yk ~ 1 + COV + f(idx, model = 'iid',
#'                  hyper = list(prec = list(prior = "loggamma", param = c(a0, b0))))
#' formula <- Yk ~ 1 + COV + f(idx, model = 'iid',hyper = list(prec = list(prior = "loggamma", param = c(a0, b0))))
#' formula= Yk~ 1  +f(COV, model = "rw1", scale.model = TRUE,hyper = list(theta = list(prior="pc.prec", param=c(1,0.01))))+
#' f(idx, model = 'iid',hyper = list(prec = list(prior = "loggamma", param = c(a0, b0))))
#' formula= Yk~ 1  + f(ids,model="iid",group=idt, control.group=list(model="ar1"),hyper = list(prec = list(prior = "loggamma", param = c(a0, b0))))
#' # Call the function
#' result <- evalLogLike_each_INLA(1, Y, X, inla.extra, membership, custom_formula, TRUE)
#' }
#' @export
evalLogLike_each_INLA <- function(k, Y, X, inla.extra, membership, formula, detailed = FALSE, correction = T) {
  N <- inla.extra$N

  nt <- dim(Y)[2]
  data_prep <- prepare_data(Y, N, membership, k, nt)
  #
  if(is.vector(X)){COV = rep(times = data_prep$n_k, X)}
  if(is.matrix(X)){COV = do.call(rbind, replicate(data_prep$n_k, X, simplify = FALSE))}
  ns = data_prep$n_k
  idx = 1:(ns*nt)
  idt = rep(1:nt,ns)
  ids = rep(1:ns, each = nt)
  inla.data <- as.data.frame(cbind(Yk = data_prep$vec_Yk, COV = COV, N = data_prep$vec_N, idx = idx, idt = idt, ids = ids))
  family <- inla.extra$family
  if(family == "poisson"){
    model <- INLA::inla(formula, family = "poisson", E = inla.data$N,
                        data = inla.data, control.predictor = list(compute = TRUE),
                        control.compute = list(config=TRUE))
  } else if(family == "binomial"){
    model <- INLA::inla(formula, family = "binomial", Ntrials = inla.data$N,
                        data = inla.data, control.predictor = list(compute = TRUE),
                        control.compute = list(config=TRUE))
  }else if(family == "nbinomial"){
    model = INLA::inla(formula, family = "binomial",
                 control.family = list(variant = 1), Ntrials = inla.data$N,
                 data = inla.data, control.predictor = list(compute = TRUE),
                 control.compute = list(config=TRUE))
  }else if(family == "normal" | is.null(family)){
    model = INLA::inla(formula, family = "normal",
                 data = inla.data, control.predictor = list(compute = TRUE),
                 control.compute = list(config=TRUE))
  }


  if (detailed) {
    return(model)
  } else {
    if (correction) {
      # correction = T
      theta <- exp(as.numeric(model[["misc"]][["configs"]][["config"]][[1]][["theta"]][[1]]))
      Q <- model[["misc"]][["configs"]][["config"]][[1]][["Qprior"]]
      dim <- dim(Q)[1] - 1
      det <- sum(diag(as.matrix(SparseM::chol(Q[1:dim,1:dim])))^2)
      mlik <- as.numeric(model[["mlik"]][[1]]) + log(det/theta^dim) * 0.5
      return(mlik)
    } else {
      return(model[["mlik"]][[1]])}
  }
}


###############################

## Auxiliary function to create data.frame for cluster `k`

prepare_data <- function(k, Y, X, N, membership) {
  ind <- which(membership == k)
  nk <- length(ind)
  nt <- nrow(Y)

  # response
  Yk <- as.vector(Y[,ind])

  # size
  if (is.vector(N)) {
    Nk <- rep(N[ind], each = nt)
  } else if (is.matrix(N)) {
    Nk <- as.vector(N[,ind])
  }

  # predictors
  if(is.vector(X)) {
    COV <- rep(X, times = n_k)
  } else if (is.matrix(X)) {
    COV <- kronecker(rep(1, nk), X)
  }

  list(
    inla_data = data.frame(Yk, Nk, id = 1:(nk*nt), idt = rep(1:nt, nk), ids = rep(1:nk, each = nt)),
    COV = COV
  )
}

