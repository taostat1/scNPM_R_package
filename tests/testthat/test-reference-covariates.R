# Independent checks for the common, centered reference feature space.
scnpm_reference_internal <- function(name) getFromNamespace(name, "scNPM")

scnpm_reference_names <- function(K) {
  c(paste0("qf_ref_T", seq_len(K)), paste0("qf_ref_C", seq_len(K)))
}

scnpm_reference_fixture <- function(K = 4L) {
  theta <- rbind(c(0, 0, 0), c(1, 2, 3), c(-1, 1, 2), c(2, -1, 0))
  mu_treat <- cbind(c(0, 0, 0), c(1, 0, -1), c(-1, 2, 0), c(2, -1, 1))
  mu_ctrl <- cbind(c(0, 1, 0), c(-1, 0, 1), c(1, 1, 1), c(0, -2, 2))
  invcov_treat <- array(0, c(3L, 3L, 4L))
  invcov_ctrl <- array(0, c(3L, 3L, 4L))
  invcov_treat[, , 1] <- diag(c(1, 2, 3))
  invcov_treat[, , 2] <- diag(c(2, 1, 1))
  invcov_treat[, , 3] <- rbind(c(2, 1, 0), c(1, 3, 0), c(0, 0, 1))
  invcov_treat[, , 4] <- rbind(c(1, 0, 0), c(0, 2, 0.5), c(0, 0.5, 2))
  invcov_ctrl[, , 1] <- diag(c(3, 1, 2))
  invcov_ctrl[, , 2] <- rbind(c(2, 0, 1), c(0, 1, 0), c(1, 0, 2))
  invcov_ctrl[, , 3] <- diag(c(1, 4, 2))
  invcov_ctrl[, , 4] <- rbind(c(3, -1, 0), c(-1, 2, 0), c(0, 0, 1))
  list(theta_t = theta,
       mu_treat = mu_treat[, seq_len(K), drop = FALSE],
       mu_ctrl = mu_ctrl[, seq_len(K), drop = FALSE],
       invcov_treat = invcov_treat[, , seq_len(K), drop = FALSE],
       invcov_ctrl = invcov_ctrl[, , seq_len(K), drop = FALSE])
}

scnpm_reference_match_args <- function(K = 4L) {
  reference <- scnpm_reference_fixture(K)
  n <- 2L * K
  rows <- rep(seq_len(nrow(reference$theta_t)), length.out = n)
  treatment <- reference$theta_t[rows, , drop = FALSE]
  treatment <- treatment + rep(c(0, 0.25), each = K)
  control <- treatment + 0.5
  control[1, ] <- treatment[1, ]  # The same cell state across conditions.
  theta <- matrix(0, 2L * n, 3L)
  it <- seq.int(1L, 2L * n, by = 2L)
  ic <- it + 1L
  theta[it, ] <- treatment
  theta[ic, ] <- control
  reference$theta_t <- theta
  reference$cell_t <- cbind(treatment = rep(c(1, 0), n),
                           quality = seq_len(2L * n) + rep(c(10, 0), n),
                           constant = 5)
  reference$group_treat <- rep(seq_len(K) - 1L, each = 2L)
  reference$group_ctrl <- rev(reference$group_treat)
  reference$treat_idx <- it - 1L
  reference$ctrl_idx <- ic - 1L
  reference$G <- 3L
  reference$K <- K
  reference$network_covariates <- "all_reference"
  reference$propensity_ridge <- 1e-6
  reference
}

scnpm_reference_match <- function(args) {
  do.call(scnpm_reference_internal("network_propensity_matching_by_difficulty"), args)
}

test_that("K=4 scores match hand-calculated centered quadratic forms in reference order", {
  actual <- do.call(scnpm_reference_internal("calculate_reference_quadratic_forms"),
                    scnpm_reference_fixture())
  # Columns are T1..T4, C1..C4. Nonzero off-diagonal entries test cross terms.
  expected <- rbind(c(0, 3, 10, 7, 1, 2, 7, 12),
                    c(36, 20, 17, 33, 22, 28, 12, 28),
                    c(15, 18, 7, 21, 11, 3, 6, 27),
                    c(6, 4, 27, 2, 16, 15, 19, 14))
  expect_equal(dim(actual), c(4L, 8L))
  expect_equal(unname(actual), expected, tolerance = 1e-12)

  two <- do.call(scnpm_reference_internal("calculate_reference_quadratic_forms"),
                scnpm_reference_fixture(2L))
  expect_equal(dim(two), c(4L, 4L))
  expect_equal(unname(two), expected[, c(1, 2, 5, 6)], tolerance = 1e-12)
})

test_that("identical cell states share reference features regardless of condition or labels", {
  args <- scnpm_reference_match_args()
  first <- scnpm_reference_match(args)
  expect_equal(as.numeric(first$network_covariates[1, ]),
               as.numeric(first$network_covariates[2, ]), tolerance = 1e-12)

  # Change every assignment without changing theta, references, or condition.
  args$group_treat <- (args$group_treat + 1L) %% args$K
  args$group_ctrl <- (args$group_ctrl + 2L) %% args$K
  relabelled <- scnpm_reference_match(args)
  expect_equal(first$network_covariates, relabelled$network_covariates,
               tolerance = 1e-12)
  expect_equal(first$propensity_scores, relabelled$propensity_scores,
               tolerance = 1e-12)
})

test_that("reference centering is invariant to a joint translation of cells and means", {
  reference <- scnpm_reference_fixture()
  calculate <- scnpm_reference_internal("calculate_reference_quadratic_forms")
  original <- do.call(calculate, reference)
  translation <- c(10, -7, 2.5)
  reference$theta_t <- sweep(reference$theta_t, 2L, translation, "+")
  reference$mu_treat <- sweep(reference$mu_treat, 1L, translation, "+")
  reference$mu_ctrl <- sweep(reference$mu_ctrl, 1L, translation, "+")
  expect_equal(do.call(calculate, reference), original, tolerance = 1e-12)
})

test_that("propensity scaling uses one pooled population for every QC and reference feature", {
  args <- scnpm_reference_match_args()
  result <- scnpm_reference_match(args)
  raw <- cbind(args$cell_t[, -1L, drop = FALSE], result$network_covariates)
  pooled_center <- colMeans(raw)
  pooled_scale <- apply(raw, 2L, stats::sd)
  pooled_scale[pooled_scale < 1e-12] <- 1
  expect_equal(as.numeric(result$propensity_scaling$center),
               as.numeric(pooled_center), tolerance = 1e-12)
  expect_equal(as.numeric(result$propensity_scaling$scale),
               as.numeric(pooled_scale), tolerance = 1e-12)

  standardized <- sweep(sweep(raw, 2L, pooled_center, "-"),
                        2L, pooled_scale, "/")
  expect_equal(as.numeric(colMeans(standardized)), rep(0, ncol(raw)),
               tolerance = 1e-12)
  # A condition-specific centering would erase this deliberate QC mean gap.
  expect_gt(abs(mean(standardized[args$treat_idx + 1L, 1L])), 0.1)
  expect_gt(abs(mean(standardized[args$ctrl_idx + 1L, 1L])), 0.1)
  expect_true(all(is.finite(result$propensity_scores)))
  expected_scores <- scnpm_reference_internal("calculate_propensity_score")(
    raw, args$cell_t[, 1L])
  expect_equal(result$propensity_scores, expected_scores, tolerance = 1e-12)
})

test_that("duplicated references remain finite and retain their shared feature values", {
  args <- scnpm_reference_match_args()
  args$mu_treat[,] <- args$mu_ctrl[,] <- 0
  args$invcov_treat <- args$invcov_ctrl <- array(diag(args$G), c(args$G, args$G, args$K))
  result <- scnpm_reference_match(args)
  expect_equal(dim(result$network_covariates), c(nrow(args$theta_t), 2L * args$K))
  for (j in seq_len(2L * args$K)) {
    expect_equal(as.numeric(result$network_covariates[, j]),
                 rowSums(args$theta_t^2), tolerance = 1e-12)
  }
  expect_true(all(is.finite(result$propensity_scores)))
  expect_true(all(result$propensity_scores > 0 & result$propensity_scores < 1))
})

test_that("explicit legacy mode reproduces the original own-cluster uncentered formula", {
  args <- scnpm_reference_match_args()
  args$network_covariates <- "legacy"
  result <- scnpm_reference_match(args)
  expected <- numeric(nrow(args$theta_t))
  for (side in c("treat", "ctrl")) {
    index <- args[[paste0(side, "_idx")]] + 1L
    cluster <- args[[paste0("group_", side)]] + 1L
    precision <- args[[paste0("invcov_", side)]]
    for (i in seq_along(index)) {
      x <- args$theta_t[index[i], ]
      expected[index[i]] <- drop(t(x) %*% precision[, , cluster[i]] %*% x)
    }
  }
  expect_equal(result$network_covariate_mode, "legacy")
  expect_equal(dim(result$network_covariates), c(length(expected), 1L))
  expect_equal(as.numeric(result$network_covariates), expected, tolerance = 1e-12)
  expect_equal(as.numeric(result$quadratic_form_global), expected, tolerance = 1e-12)
  old_probabilities <- scnpm_reference_internal("calculate_propensity_score")(
    cbind(args$cell_t[, -1L, drop = FALSE], expected), args$cell_t[, 1L])
  expect_equal(result$propensity_scores, old_probabilities, tolerance = 1e-12)
  # Old low-level callers omit means and mode; their result must stay legacy.
  args$mu_treat <- args$mu_ctrl <- args$network_covariates <- args$propensity_ridge <- NULL
  compatible <- scnpm_reference_match(args)
  expect_equal(compatible$network_covariates, result$network_covariates)
  expect_equal(compatible$propensity_scores, result$propensity_scores)
})

scnpm_reference_public_args <- function(K) {
  n <- 24L
  cell <- seq_len(n)
  expression <- rbind(A = 2 + cell %% 7,
                      B = 3 + (cell * 3) %% 11,
                      C = 1 + (cell * 5) %% 13)
  colnames(expression) <- paste0("cell", cell)
  list(inference_mode = "coupled", goal_gene_mt = expression,
       cell_t = cbind(treatment = rep(c(1, 0), n / 2L),
                      quality = sin(cell / 3)),
       K = K, ds_seed = 193L, num_iter = 2L, num_save = 1L,
       epsilon_theta = 1e-8, num_step_theta = 1L,
       epsilon_lam = 0.001, num_step_lam = 1L)
}

test_that("public default fits expose all 2K reference features and featurewise balance", {
  set_threads <- scnpm_reference_internal("set_scNPM_threads")
  previous <- set_threads(1L)
  on.exit(set_threads(previous), add = TRUE)
  for (K in c(2L, 4L)) {
    args <- scnpm_reference_public_args(K)
    result <- NULL
    invisible(capture.output(result <- suppressMessages(do.call(run_scNPM, args))))
    feature_names <- scnpm_reference_names(K)
    expect_equal(result$network_covariate_mode, "all_reference")
    expect_equal(dim(result$network_covariates), c(ncol(args$goal_gene_mt), 2L * K))
    expect_equal(colnames(result$network_covariates), feature_names)
    expect_true(all(is.finite(result$network_covariates)))
    expect_true(all(is.finite(result$propensity_scores)))
    expect_equal(colnames(result$smd_covariate_results), c("quality", feature_names))
    expect_equal(length(result$propensity_scaling$center), 1L + 2L * K)
    expect_equal(length(result$propensity_scaling$scale), 1L + 2L * K)
    expect_true(all(is.finite(result$propensity_scaling$center)))
    expect_true(all(is.finite(result$propensity_scaling$scale)))
    expect_equal(result$provenance$settings$network_covariates, "all_reference")
    expect_equal(result$provenance$settings$n_network_covariates, 2L * K)
    expect_equal(result$provenance$settings$network_feature_names, feature_names)

    balance <- compute_qf_smd(result)
    expect_equal(attr(balance, "mode"), "exact-all-reference")
    expect_equal(nrow(balance), 2L * K * K)
    expect_setequal(unique(balance$network_feature), feature_names)
    for (feature in feature_names) {
      rows <- balance[balance$network_feature == feature, , drop = FALSE]
      expect_equal(sort(rows$treat_cluster), seq_len(K))
      expect_equal(rows$qf_smd,
                   unname(round(result$smd_covariate_results[rows$treat_cluster, feature], 4)),
                   tolerance = 1e-4)
    }
  }
})

test_that("featurewise SMD uses the final common feature matrix rather than legacy diagnostics", {
  feature_names <- scnpm_reference_names(2L)
  features <- cbind(c(1, 2, 2, 1, 3, 3),
                    c(2, 1, 3, 2, 2, 4),
                    c(1, 2, 4, 3, 2, 6),
                    c(0, 1, 5, 0, 3, 3))
  colnames(features) <- feature_names
  result <- list(K = 2L, network_covariate_mode = "all_reference",
                 network_covariates = features,
                 matching_result_df = data.frame(treatment_cluster = 1:2,
                                                 matched_control_cluster = 1:2),
                 final_group_label = c(0, 1, 1, 2, 3, 3),
                 qf_global = rep(999, 6))
  actual <- compute_qf_smd(result)
  expect_equal(attr(actual, "mode"), "exact-all-reference")
  expect_equal(nrow(actual), 8L)
  expected <- rbind(c(0, 0, Inf, 0), c(Inf, 1 / sqrt(2), 1 / sqrt(5), 0))
  for (j in seq_along(feature_names)) {
    rows <- actual[actual$network_feature == feature_names[j], , drop = FALSE]
    expect_equal(rows$qf_smd[match(1:2, rows$treat_cluster)], round(expected[, j], 4))
  }
})
