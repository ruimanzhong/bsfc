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
#' \dontrun{
#' # Define data matrices and parameters
#' Y <- matrix(rpois(100, lambda = 10), ncol = 10)
#' population <- matrix(rpois(100, lambda = 5), ncol = 10)
#' inla.extra <- list(family = "poisson", N = population)
#' membership <- sample(1:2, 10, replace = TRUE)
#' X <- seq(0, 1, length.out = 10)
#' a0 <- 0.1
#' b0 <- 0.1
#'
#' # Define a formula, the name of the iid random effec should be idx
#' formula <- Yk ~ 1 + Xk + f(idx,
#'   model = "iid",
#'   hyper = list(prec = list(prior = "loggamma", param = c(a0, b0)))
#' )
#' formula <- Yk ~ 1 + Xk + f(idx, model = "iid", hyper = list(prec = list(prior = "loggamma", param = c(a0, b0))))
#' formula <- Yk~ 1 + f(Xk, model = "rw1", scale.model = TRUE, hyper = list(theta = list(prior = "pc.prec", param = c(1, 0.01)))) +
#'   f(idx, model = "iid", hyper = list(prec = list(prior = "loggamma", param = c(a0, b0))))
#' formula <- Yk ~ 1 + f(ids, model = "iid", group = idt, control.group = list(model = "ar1"), hyper = list(prec = list(prior = "loggamma", param = c(a0, b0))))
#' # Call the function
#' result <- log_mlik_each(1, Y, X, inla.extra, membership, custom_formula, TRUE)
#' }
#' @export
log_mlik_each <- function(k, Y, membership, X = NULL, N = NULL, formula = Yk ~ 1 + Xk,
                          family = "normal", correction = FALSE, detailed = FALSE, ...) {
  inla_data <- prepare_data_each(k, Y, membership, X, N)
  if (family == "poisson") {
    model <- INLA::inla(formula, family,
      E = Nk,
      data = inla_data, control.predictor = list(compute = TRUE),
      control.compute = list(config = TRUE), ...
    )
  } else if (family == "binomial") {
    model <- INLA::inla(formula, family,
      Ntrials = Nk,
      data = inla_data, control.predictor = list(compute = TRUE),
      control.compute = list(config = TRUE), ...
    )
  } else if (family == "nbinomial") {
    model <- INLA::inla(formula, family,
      control.family = list(variant = 1), Ntrials = Nk,
      data = inla_data, control.predictor = list(compute = TRUE),
      control.compute = list(config = TRUE), ...
    )
  } else if (family == "normal" | is.null(family)) {
    model <- INLA::inla(formula, family,
      data = inla_data, control.predictor = list(compute = TRUE),
      control.compute = list(config = TRUE), ...
    )
  }

  if (detailed) {
    model
  } else {
    if (correction) {
      log_mlik_corrected(model, formula)
    } else {
      model[["mlik"]][[1]]
    }
  }
}

#' Calculate Log Marginal Likelihood for All Clusters
#'
#' This function computes the log marginal likelihood for each cluster specified in the membership vector.
#' It applies a specified model to each cluster to calculate the likelihood.
#'
#' @param Y Numeric matrix or data frame of observations, where rows typically represent observational units
#'        and columns represent different variables.
#' @param membership Integer vector indicating the cluster membership for each observation.
#' @param X Optional numeric matrix or data frame of covariates used in the model (if applicable).
#' @param N Optional numeric vector indicating the number of trials or cases, relevant for
#'        distributions like binomial (default is NULL).
#' @param formula An object of class \code{\link[stats]{formula}} specifying the model to be used in INLA.
#'        Default is `Yk ~ 1 + Xk`, which can be adjusted based on the model requirements.
#' @param family Character string specifying the family of distributions to use for the model.
#'        Defaults to "normal". Other possible values include "binomial", "poisson", etc.
#' @param correction Logical indicating whether a correction for dispersion or other factors
#'        should be applied. Defaults to FALSE.
#' @param ... Additional arguments passed to the underlying `log_mlik_each` function.
#'
#' @details This function iterates over the unique clusters defined in the `membership` vector and
#'          calculates the log marginal likelihood for each cluster using the `log_mlik_each` function.
#'          The function is flexible and allows for different families of distributions and model specifications
#'          by adjusting the formula, family, and other parameters.
#'
#' @return A numeric vector containing the log marginal likelihood for each cluster.
#'
#' @examples
#' \dontrun{
#' # Assume Y is a matrix of observations, membership defines cluster memberships
#' # and X is a matrix of covariates:
#' log_mlikelihoods <- log_mlik_all(Y, membership, X = X, family = "poisson")
#' }
#'
#' @export
log_mlik_all <- function(Y, membership, X = NULL, N = NULL, formula = Yk ~ 1 + Xk,
                         family = "normal", correction = FALSE, ...) {
  k <- max(membership)
  sapply(1:k, log_mlik_each, Y, membership, X, N, formula, family, correction, FALSE, ...)
}

##############################################################

## Auxiliary function to create data.frame for cluster `k`

prepare_data_each <- function(k, Y, membership, X = NULL, N = NULL) {
  ind <- which(membership == k)
  nk <- length(ind)
  nt <- nrow(Y)

  # response
  Yk <- as.vector(Y[, ind])

  # size
  if (is.vector(N)) {
    Nk <- rep(N[ind], each = nt)
  } else if (is.matrix(N)) {
    Nk <- as.vector(N[, ind])
  } else {
    Nk <- NULL
  }

  # predictors
  if (is.vector(X)) {
    Xk <- rep(X, times = nk)
  } else if (is.matrix(X)) {
    Xk <- kronecker(rep(1, nk), X)
  } else {
    Xk <- NULL
  }

  list(
    Yk = Yk, Nk = Nk, id = 1:(nk * nt), idt = rep(1:nt, nk), ids = rep(1:nk, each = nt),
    Xk = Xk
  )
}

##############################################################

## Auxiliary function to correct the marginal likelihood of INLA model
get_structure_matrix = function (model, formula) {

  model = model[["misc"]][["configs"]]

  # effects dimension information
  x_info = model[["contents"]]
  ef_start = setNames(x_info$start[-1] - x_info$length[1], x_info$tag[-1])
  ef_end = ef_start + x_info$length[-1] - 1

  # select effect that requires correction
  fs = as.list(attr(terms(formula), "variables"))[c(-1,-2)]
  pattern <- "f\\((\\w+)(?=[^,]*?, (?!model = \"iid\"))"
   results <- stringr::str_extract_all(sapply(fs, deparse), pattern)
  fs_vars <- sapply(unlist(results), function(x) gsub("f\\(|\\)", "", x))

  # provide structure matrix for selected effects
  out = list()
  for (x in fs_vars) {
    i = ef_start[x]; j = ef_end[x]
    out[[x]] = model[["config"]][[1]][["Qprior"]][i:j, i:j] /
      exp(model[["config"]][[1]][["theta"]][paste0("Log precision for ", x)])
  }

  return(out)
}

log_mlik_corrected <- function(model, formula) {
  Slist <- get_structure_matrix(model, formula)
  Slogdet <- sapply(Slist, function(x) 2 * sum(log(Matrix::diag(SparseM::chol(x)))))
  model[["mlik"]][[1]] + 0.5 * sum(Slogdet)
}
