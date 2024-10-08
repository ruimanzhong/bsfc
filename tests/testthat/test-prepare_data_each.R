
library(testthat)
library(sf)
library(sftime)
library(dplyr)
library(tidyr)

# Create an example test for the function `prepare_data_each`

test_that("prepare_data_each works as expected", {

  # Spatial coordinates (3 regions)
  coords <- st_sfc(st_point(c(0, 0)), st_point(c(1, 1)), st_point(c(2, 2)))

  # Time variable (2 time points)
  time_var <- as.POSIXct(c('2022-01-01', '2022-01-02'))

  # Case data (for each region and time point)
  case <- matrix(runif(6, 50, 100), nrow = 2, ncol = 3)  # 2 time points x 3 regions
  temperature <- matrix(runif(6, 10, 20), nrow = 2, ncol = 3)
  precipitation <- matrix(runif(6, 0, 50), nrow = 2, ncol = 3)
  population <- matrix(runif(6, 100, 200), nrow = 2, ncol = 3)

  # Create sftime object
  data <- st_sftime(
    data.frame(
      case = as.vector(case),
      temperature = as.vector(temperature),
      precipitation = as.vector(precipitation),
      population = as.vector(population),
      time = rep(time_var, each = 3)
    ),
    geometry = coords
  )

  # Example membership vector (3 regions, 2 clusters)
  membership <- c(1, 1, 2)

  # Define the formula (as used in the question)
  formula <- case ~ temperature + precipitation + f(idt, model = "iid")

  # Call the prepare_data_each function
  result <- prepare_data_each(k = 1, data = data, membership = membership, formula = formula)

  # Expected values
  expected_nrows <- 2 * sum(membership == 1)  # 2 time points * number of regions in cluster 1

  # Tests
  expect_equal(nrow(result), expected_nrows)  # Ensure correct number of rows
  expect_equal(result$case, as.vector(case[, membership == 1]))  # Check that case values are correct
  expect_equal(result$temperature, as.vector(temperature[, membership == 1]))  # Check temperature
  expect_equal(result$precipitation, as.vector(precipitation[, membership == 1]))  # Check precipitation
  expect_equal(result$id, 1:expected_nrows)  # Check the id column
  expect_equal(result$idt, rep(1:2, sum(membership == 1)))  # Ensure time index is correct
  expect_equal(result$ids, rep(1:sum(membership == 1), each = 2))  # Check region ids
})
