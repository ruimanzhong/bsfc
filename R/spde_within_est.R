evalLogLike_each_INLA <- function(k, Y, X, population, membership, a0 = NULL, b0 = NULL, detailed = F) {
  nt <- as.numeric(dim(Y)[2])
  ind <- which(membership == k)
  Yk <- matrix(Y[ind, ], ncol = nt)
  ns <- dim(Yk)[1]
  popk <- matrix(population[ind, ], ncol = nt)
  vec_Yk <- as.vector(t(Yk))
  vec_population <- as.vector(t(popk))
  knots <- seq(round(min(X, na.rm = T) - 0.5), round(max(X, na.rm = T) + 0.5), length = 12)
  mesh1d <- INLA::inla.mesh.1d(loc = knots, degree = 2, boundary = "free")
  spde1d <- INLA::inla.spde2.pcmatern(
    mesh = mesh1d, alpha = 2, ### mesh and smoothness parameter
    prior.range = c(a0, NA), ### Prange == NA means range is fixed at 20
    prior.sigma = c(b0, 0.1), ### P(sigma>1)=0.1

    constr = TRUE
  )
  COV <- rep(times = nrow(Yk), X)
  A.1d <- INLA::inla.spde.make.A(mesh1d, COV)
  covar.index <- INLA::inla.spde.make.index("Cov", n.spde = spde1d$n.spde)
  idx <- 1:length(vec_Yk)
  # stack
  formula <- y ~ 0 + intercept + f(Cov, model = spde1d) + f(idx, model = "iid", hyper = list(prec = list(prior = "loggamma", param = c(0.1, 0.1))))
  data <- as.data.frame(cbind(vec_Yk = vec_Yk, E = vec_population))
  stack <- INLA::inla.stack(
    data = list(y = data$vec_Yk, e = data$E),
    A = list(1, A.1d, 1),
    effects = list(intercept = rep(1, nt * ns), covar.index, idx = 1:length(vec_Yk)),
    tag = "est"
  )

  m.spde <- INLA::inla(formula,
    data = INLA::inla.stack.data(stack),
    family = "poisson", E = INLA::inla.stack.data(stack)$e,
    control.predictor = list(A = INLA::inla.stack.A(stack), compute = TRUE)
  )
  summary(m.spde)
  if (detailed) {
    return(m.spde)
  } else {
    return(m.spde[["mlik"]][[1]])
  }
}
