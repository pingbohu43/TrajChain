test_that("built-in kernels are densities supported on [-1, 1]", {
  for (k in c("epanechnikov", "biweight", "triweight", "triangular", "uniform")) {
    K <- kernel_function(k)
    expect_equal(integrate(K, -1, 1)$value, 1, tolerance = 1e-8, label = k)
    expect_equal(K(c(-3, -1.0001, 1.0001, 3)), rep(0, 4), label = k)
    expect_true(all(K(seq(-0.99, 0.99, by = 0.01)) > 0), label = k)
    expect_equal(K(0.3), K(-0.3), label = k)
  }
})

test_that("kernels preserve matrix dimensions", {
  M <- matrix(seq(-2, 2, length.out = 12), 3, 4)
  for (k in c("epanechnikov", "biweight", "triweight", "triangular", "uniform")) {
    expect_equal(dim(kernel_function(k)(M)), dim(M))
  }
})

test_that("the Epanechnikov kernel matches its formula", {
  u <- c(-1, -0.5, 0, 0.25, 1, 1.5)
  expect_equal(kernel_function()(u), ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0))
})

test_that("user-supplied kernels are validated", {
  quartic <- function(u) {
    out <- (15 / 16) * (1 - u^2)^2
    out[abs(u) > 1] <- 0
    out
  }
  expect_identical(kernel_function(quartic), quartic)
  expect_error(kernel_function(dnorm), "vanish outside")
  expect_error(kernel_function(function(u) -abs(u)), "non-negative")
  expect_error(kernel_function(function(u) sum(u)), "one value per input")
  expect_error(kernel_function("gaussian"))
})
