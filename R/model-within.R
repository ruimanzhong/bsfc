#' Evaluate Log Marginal Likelihood for Each Cluster Using INLA
#'
#' This function computes the log marginal likelihood for a specific cluster `k` within a dataset. It supports multiple family distributions and allows for the inclusion of exposure or trial sizes if applicable. The function is flexible, accepting custom formulas and additional arguments for the INLA model fitting.
#'
#' @param k Integer, the cluster index for which to compute the log marginal likelihood.
#' @param Y Numeric vector or matrix, representing the response variable(s) in the model.
#' @param membership Integer vector, indicating the cluster membership of each observation.
#' @param X Numeric matrix or NULL, optional covariates for the model; defaults to NULL if not provided.
#' @param N Numeric vector or NULL, representing the number of trials or exposures for each observation (necessary for binomial and Poisson models); defaults to NULL.
#' @param family Character string specifying the family of the model, with support for "normal", "poisson", "binomial", and "nbinomial" (negative binomial variant 1); defaults to "normal".
#' @param formula Object of class \code{\link[stats]{formula}}, specifying the model to be fitted; defaults to `Yk ~ 1 + Xk`.
#' @param correction Logical, indicating whether a correction should be applied to the computed log marginal likelihood; defaults to FALSE.
#' @param detailed Logical, specifying whether to return the full model object or just the log marginal likelihood; defaults to FALSE.
#' @param ... Additional arguments passed to the INLA function.
#'
#' @return Depending on the `detailed` parameter:
#'   - If `detailed` is TRUE, returns the full `INLA` model object.
#'   - If `detailed` is FALSE, returns the log marginal likelihood (corrected if `correction` is TRUE).
#'
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
#' formula <- Yk ~ 1 + Xk + f(idx, model = 'iid',
#'                  hyper = list(prec = list(prior = "loggamma", param = c(a0, b0))))
#' formula <- Yk ~ 1 + Xk + f(idx, model = 'iid',hyper = list(prec = list(prior = "loggamma", param = c(a0, b0))))
#' formula= Yk~ 1  +f(Xk, model = "rw1", scale.model = TRUE,hyper = list(theta = list(prior="pc.prec", param=c(1,0.01))))+
#' f(idx, model = 'iid',hyper = list(prec = list(prior = "loggamma", param = c(a0, b0))))
#' formula= Yk ~ 1  + f(ids,model="iid",group=idt, control.group=list(model="ar1"),hyper = list(prec = list(prior = "loggamma", param = c(a0, b0))))
#' # Call the function
#' result <- log_mlik_each(1, Y, X, inla.extra, membership, custom_formula, TRUE)
#' }
#' @export
log_mlik_each <- function(k, Y, membership, X = NULL, N = NULL, formula = Yk ~ 1 + Xk,
                          family = "normal", correction = FALSE, detailed = FALSE, ...) {
  inla_data <- prepare_data_each(k, Y, membership, X, N)
  if(family == "poisson"){
    model <- INLA::inla(formula, family, E = Nk,
                        data = inla_data, control.predictor = list(compute = TRUE),
                        control.compute = list(config=TRUE), ...)
  } else if (family == "binomial"){
    model <- INLA::inla(formula, family, Ntrials = Nk,
                        data = inla_data, control.predictor = list(compute = TRUE),
                        control.compute = list(config=TRUE), ...)
  } else if (family == "nbinomial"){
    model <- INLA::inla(formula, family,
                 control.family = list(variant = 1), Ntrials = Nk,
                 data = inla_data, control.predictor = list(compute = TRUE),
                 control.compute = list(config=TRUE), ...)
  } else if (family == "normal" | is.null(family)){
    model <- INLA::inla(formula, family,
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

log_mlik_all <- function(Y, membership, X = NULL, N = NULL, formula = Yk ~ 1 + Xk,
                         family = "normal", correction = FALSE, ...) {

  k = max(membership)
  sapply(1:k, log_mlik_each, Y, membership, X, N, formula, family, correction, FALSE, ...)
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
    Xk <- rep(X, times = nk)
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
  theta <- exp(as.numeric(inla_model[["misc"]][["configs"]][["config"]][[1]][["theta"]][[2]]))
  Q <- inla_model[["misc"]][["configs"]][["config"]][[1]][["Qprior"]]
  dim <- dim(Q)[1] - 1
  det <- sum(diag(as.matrix(SparseM::chol(Q[1:dim,1:dim])))^2)
  as.numeric(model[["mlik"]][[1]]) + log(det) * 0.5 - dim * log(theta)
}

