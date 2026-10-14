test_that("nonconvergent full steps fall back to the same ridge optimum", {
  set.seed(108)
  X <- matrix(rnorm(120 * 4), 120, 4)
  y <- rep(0:1, each = 60)
  X[, 1] <- X[, 1] + y
  predictor <- getFromNamespace(".pairwise_ridge_predict", "scNPM")
  reference <- predictor(X, y, X)
  # Fault injection represents a failed full-step solve; the fallback must
  # recover the unchanged objective instead of returning an unfinished fit.
  calls <- 0L
  controlled_solve <- function(a, b, ...) {
    calls <<- calls + 1L
    if (calls <= 100L) rep(.1, ncol(a)) else base::solve(a, b, ...)
  }
  environment(predictor) <- list2env(list(solve = controlled_solve),
    parent = environment(predictor))
  rng <- .Random.seed
  recovered <- predictor(X, y, X)
  expect_gt(calls, 100L)
  expect_identical(.Random.seed, rng)
  expect_equal(recovered, reference, tolerance = 1e-7)
  expect_true(all(is.finite(recovered) & recovered > 0 & recovered < 1))

  A <- cbind(1, scale(X))
  w <- ifelse(y == 1, .5 / mean(y == 1), .5 / mean(y == 0))
  penalty <- c(1e-8, rep(.001 * sum(w), ncol(X)))
  objective <- function(b) {
    eta <- drop(A %*% b)
    sum(w * (pmax(eta, 0) + log1p(exp(-abs(eta))) - y * eta)) +
      sum(penalty * b^2) / 2
  }
  gradient <- function(b) drop(crossprod(A, w * (plogis(A %*% b) - y))) + penalty * b
  fit <- stats::optim(numeric(ncol(A)), objective, gradient, method = "BFGS",
    control = list(maxit = 1000, reltol = 1e-14))
  expect_equal(fit$convergence, 0L)
  expect_equal(recovered, plogis(drop(A %*% fit$par)), tolerance = 1e-6)
})
