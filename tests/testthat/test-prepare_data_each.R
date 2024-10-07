library(testthat)
test_that("prepare_data_each works as expected", {

  # Example inputs
  Y <- matrix(1:12, nrow = 4, ncol = 3)  # Response matrix (4 time points, 3 regions)
  membership <- c(1, 2, 1)  # Membership vector
  X <- matrix(1:8, nrow = 4, ncol = 2)  # Predictor matrix (4 time points, 2 predictors)
  N <- c(100, 200, 300)  # Population size vector

  # Case 1: Test with matrix X and vector N
  result <- prepare_data_each(k = 1, Y = Y, membership = membership, X = X, N = N)

  expect_equal(result$Y, as.vector(Y[, membership == 1]))  # Check Yk
  expect_equal(result$N, rep(N[membership == 1], each = 4))  # Check Nk (size)
  expect_equal(result$id, 1:(sum(membership == 1) * 4))  # Check id
  expect_equal(result$idt, rep(1:4, sum(membership == 1)))  # Check idt
  expect_equal(result$ids, rep(1:sum(membership == 1), each = 4))  # Check ids
  expect_equal(result$X, kronecker(rep(1, sum(membership == 1)), X))  # Check Xk

  # Case 2: Test with NULL X and vector N
  result <- prepare_data_each(k = 2, Y = Y, membership = membership, X = NULL, N = N)

  expect_equal(result$Y, as.vector(Y[, membership == 2]))  # Check Yk
  expect_equal(result$N, rep(N[membership == 2], each = 4))  # Check Nk (size)
  expect_equal(result$id, 1:(sum(membership == 2) * 4))  # Check id
  expect_equal(result$idt, rep(1:4, sum(membership == 2)))  # Check idt
  expect_equal(result$ids, rep(1:sum(membership == 2), each = 4))  # Check ids
  expect_null(result$X)  # X should be NULL

  # Case 3: Test with NULL X and NULL N
  result <- prepare_data_each(k = 1, Y = Y, membership = membership, X = NULL, N = NULL)

  expect_equal(result$Y, as.vector(Y[, membership == 1]))  # Check Yk
  expect_null(result$N)  # N should be NULL
  expect_equal(result$id, 1:(sum(membership == 1) * 4))  # Check id
  expect_equal(result$idt, rep(1:4, sum(membership == 1)))  # Check idt
  expect_equal(result$ids, rep(1:sum(membership == 1), each = 4))  # Check ids
  expect_null(result$X)  # X should be NULL

  # Case 4: Test with vector X and matrix N
  N_matrix <- matrix(1:12, nrow = 4, ncol = 3)  # Population matrix
  result <- prepare_data_each(k = 1, Y = Y, membership = membership, X = c(5, 10), N = N_matrix)

  expect_equal(result$Y, as.vector(Y[, membership == 1]))  # Check Yk
  expect_equal(result$N, as.vector(N_matrix[, membership == 1]))  # Check Nk (matrix)
  expect_equal(result$id, 1:(sum(membership == 1) * 4))  # Check id
  expect_equal(result$idt, rep(1:4, sum(membership == 1)))  # Check idt
  expect_equal(result$ids, rep(1:sum(membership == 1), each = 4))  # Check ids
  expect_equal(result$X, rep(c(5, 10), times = sum(membership == 1)))  # Check Xk

})
