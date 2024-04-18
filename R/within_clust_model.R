#' Evaluate Log-Marginal-Likelihood for Each Group in INLA
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
evalLogMLike_each_INLA <- function(k, Y, membership, X = NULL, N = NULL, family = "normal",
                                  formula = Yk ~ 1 + Xk, correction = FALSE, detailed = FALSE, ...) {

  inla_data <- prepare_data_each(k, Y, membership, X, N)

  if(family == "poisson"){
    model <- INLA::inla(formula, family = "poisson", E = inla_data$N,
                        data = inla_data, control.predictor = list(compute = TRUE),
                        control.compute = list(config=TRUE), ...)
  } else if (family == "binomial"){
    model <- INLA::inla(formula, family = "binomial", Ntrials = inla_data$N,
                        data = inla_data, control.predictor = list(compute = TRUE),
                        control.compute = list(config=TRUE), ...)
  } else if (family == "nbinomial"){
    model <- INLA::inla(formula, family = "binomial",
                 control.family = list(variant = 1), Ntrials = inla_data$N,
                 data = inla_data, control.predictor = list(compute = TRUE),
                 control.compute = list(config=TRUE), ...)
  } else if (family == "normal" | is.null(family)){
    model <- INLA::inla(formula, family = "normal",
                 data = inla_data, control.predictor = list(compute = TRUE),
                 control.compute = list(config=TRUE), ...)
  }

  if (detailed) {
    model
  } else {
    if (correction) {
      mlik_corrected(model)
    } else {
      model[["mlik"]][[1]]
    }
  }
}


###############################

## Auxiliary function to create data.frame for cluster `k`

prepare_data_each <- function(k, Y, membership, X = NULL, N = NULL) {
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
  } else {
    Nk <- NULL
  }

  # predictors
  if(is.vector(X)) {
    Xk <- rep(X, times = n_k)
  } else if (is.matrix(X)) {
    Xk <- kronecker(rep(1, nk), X)
  } else {
    Xk <- NULL
  }

  list(
    Yk = Yk, Nk = Nk, id = 1:(nk*nt), idt = rep(1:nt, nk), ids = rep(1:nk, each = nt),
    Xk = Xk
  )
}

## Auxiliary function to correct the marginal likelihood of INLA model

mlik_corrected <- function(inla_model) {
  theta <- exp(as.numeric(inla_model[["misc"]][["configs"]][["config"]][[1]][["theta"]][[1]]))
  Q <- inla_model[["misc"]][["configs"]][["config"]][[1]][["Qprior"]]
  dim <- dim(Q)[1] - 1
  det <- sum(diag(as.matrix(SparseM::chol(Q[1:dim,1:dim])))^2)
  as.numeric(inla_model[["mlik"]][[1]]) + log(det/theta^dim) * 0.5
}

