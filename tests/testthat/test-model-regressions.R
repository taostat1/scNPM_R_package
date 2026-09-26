# Small deterministic regression cases for the probabilistic model and matching.
# These exercise scientific invariants and independently derived reference values.
scnpm_internal <- function(name) getFromNamespace(name, "scNPM")

with_scnpm_single_thread <- function(code) {
  set_threads <- scnpm_internal("set_scNPM_threads")
  previous <- set_threads(1L)
  on.exit(set_threads(previous), add = TRUE)
  force(code)
}

scnpm_tiny_chain_args <- function(edge_prior = 0.2, missing = TRUE,
                                  num_iter = 24L, num_save = 12L) {
  G <- 3L
  N <- 12L
  theta <- matrix(seq(-0.6, 0.9, length.out = N * G), N, G)
  is_zero <- matrix(FALSE, N, G)
  if (missing) {
    is_zero[1, 1] <- TRUE
    is_zero[9, 2] <- TRUE
    theta[is_zero] <- 0
  }
  prior <- matrix(edge_prior, G, G)
  diag(prior) <- 0
  list(
    num_iter = num_iter, num_save = num_save,
    theta_t = theta, ind_zero = is_zero,
    mu_treat = matrix(colMeans(theta[1:6, , drop = FALSE]), G, 1),
    mu_ctrl = matrix(colMeans(theta[7:12, , drop = FALSE]), G, 1),
    invcov_treat = array(diag(G), c(G, G, 1)),
    cov_treat = array(diag(G), c(G, G, 1)),
    invcov_ctrl = array(diag(G), c(G, G, 1)),
    cov_ctrl = array(diag(G), c(G, G, 1)),
    edge_treat = array(0, c(G, G, 1)),
    edge_ctrl = array(0, c(G, G, 1)),
    group_treat = rep(1, 6), group_ctrl = rep(1, 6),
    lambda0_t = rep(0.3, G), lambda1_t = rep(-1, G),
    pi_treat = 1, pi_ctrl = 1, gam = 1,
    cell_t = cbind(treatment = rep(c(1, 0), each = 6),
                   quality = rep(seq(-1, 1, length.out = 6), 2)),
    G = G, N_int = N, K = 1L,
    ssp_v0 = 0.1, ssp_v1 = 1, ssp_l = 1, ssp_xi = 0.2,
    epsilon_theta = 0.02, num_step_theta = 3L,
    epsilon_lam = 0.01, num_step_lam = 3L,
    edge_prob_threshold = 0.6, edge_prob_adj = 1,
    edge_prior = prior
  )
}

scnpm_run_tiny_chain <- function(args, seed = 1749L) {
  # Rcpp may use mutable cube views; each run must start from a fresh state.
  args <- unserialize(serialize(args, NULL))
  set.seed(seed)
  result <- NULL
  invisible(capture.output(
    result <- do.call(scnpm_internal("MCMC_with_difficulty_based_matching"), args)
  ))
  result
}

scnpm_match_fixture <- function(treat_quality, ctrl_quality, labels) {
  n <- length(labels)
  K <- max(labels) + 1L
  G <- 2L
  scnpm_internal("network_propensity_matching_by_difficulty")(
    theta_t = matrix(0, 2L * n, G),
    cell_t = cbind(treatment = rep(c(1, 0), each = n),
                   quality = c(treat_quality, ctrl_quality)),
    group_treat = labels, group_ctrl = labels,
    treat_idx = seq_len(n) - 1L,
    ctrl_idx = n + seq_len(n) - 1L,
    invcov_treat = array(diag(G), c(G, G, K)),
    invcov_ctrl = array(diag(G), c(G, G, K)),
    G = G, K = K
  )
}

expect_scnpm_bijection <- function(result, K) {
  matching <- result$matching_result_df
  expect_equal(result$matched_pairs, K)
  expect_false(anyNA(matching$matched_control_cluster))
  expect_equal(sort(as.integer(matching$matched_control_cluster)), seq_len(K))
  expect_equal(sort(as.integer(matching$treatment_cluster)), seq_len(K))
}

test_that("edge probabilities agree with a Bernoulli spike-and-slab model", {
  probability <- scnpm_internal("calculate_edge_probability")
  cases <- expand.grid(omega = c(0, 0.03, 0.2, -0.2, 3),
                       xi = c(0.01, 0.2, 0.8))
  for (i in seq_len(nrow(cases))) {
    omega <- cases$omega[i]
    xi <- cases$xi[i]
    log_slab <- log(xi) + dnorm(omega, sd = 1, log = TRUE)
    log_spike <- log1p(-xi) + dnorm(omega, sd = 0.05, log = TRUE)
    expected <- plogis(log_slab - log_spike)
    expect_equal(probability(omega, 0.05, 1, xi), expected,
                 tolerance = 1e-12)
  }
  # With indistinguishable components, observing omega conveys no edge evidence.
  expect_equal(probability(0.7, 0.4, 0.4, 0.17), 0.17, tolerance = 1e-12)
  expect_equal(probability(-0.7, 0.4, 0.4, 0.17), 0.17, tolerance = 1e-12)
})

test_that("short chains reproduce and preserve latent and precision constraints", {
  with_scnpm_single_thread({
    args <- scnpm_tiny_chain_args()
    first <- scnpm_run_tiny_chain(args)
    second <- scnpm_run_tiny_chain(args)
    expect_identical(first, second)
    observed <- !args$ind_zero
    expect_equal(first$theta_est[observed], args$theta_t[observed],
                 tolerance = 1e-12)
    expect_true(all(is.finite(first$theta_est)))
    expect_true(all(first$lam1_iter < 0))
    for (side in c("treat", "ctrl")) {
      precision_samples <- first[[paste0("invcov_", side, "_iter")]]
      for (sample in precision_samples) {
        for (k in seq_len(dim(sample)[3])) {
          precision <- sample[, , k]
          expect_equal(precision, t(precision), tolerance = 1e-12)
          expect_gt(min(eigen(precision, symmetric = TRUE,
                              only.values = TRUE)$values), 0)
        }
      }
      edges <- first[[paste0("edge_", side, "_est")]][, , 1]
      expect_equal(edges, t(edges))
      expect_true(all(edges >= 0 & edges <= 1))
    }
  })
})

test_that("edge priors act during sampling even when initial graphs are identical", {
  with_scnpm_single_thread({
    low <- scnpm_tiny_chain_args(edge_prior = 1e-9, missing = FALSE,
                                num_iter = 80L, num_save = 40L)
    high <- low
    high$edge_prior[upper.tri(high$edge_prior) | lower.tri(high$edge_prior)] <- 1 - 1e-9
    # Initial Z, all continuous states, random seed, and matching adjustment agree.
    expect_identical(low$edge_treat, high$edge_treat)
    expect_identical(low$edge_ctrl, high$edge_ctrl)
    low_result <- scnpm_run_tiny_chain(low, seed = 8521L)
    high_result <- scnpm_run_tiny_chain(high, seed = 8521L)
    for (side in c("treat", "ctrl")) {
      low_edges <- low_result[[paste0("edge_", side, "_est")]][, , 1]
      high_edges <- high_result[[paste0("edge_", side, "_est")]][, , 1]
      expect_lt(mean(low_edges[upper.tri(low_edges)]), 0.2)
      expect_gt(mean(high_edges[upper.tri(high_edges)]), 0.8)
    }
  })
})

test_that("matching keeps zero-priority clusters when finite candidates are available", {
  with_scnpm_single_thread({
    # Cluster 1 has the overall treatment mean score, hence zero priority.
    # Every control cluster contains both score values, giving finite SMDs.
    # Version 1.0.1 omitted cluster 1 on this fixture.
    treat <- c(-2, 0.5, -2, -2, 0.5, 0.5)
    result <- scnpm_match_fixture(treat, rep(c(-4, -1.5), 3), rep(0:2, each = 2))
    expect_scnpm_bijection(result, 3L)
    scores <- as.numeric(result$propensity_scores)[1:6]
    priority <- abs(tapply(scores, rep(0:2, each = 2), mean) - mean(scores))
    expect_lt(priority[1], 1e-12)
    expect_gt(max(priority[-1]), 0.1)
    expect_false(anyNA(result$smd_results))
    expect_true(all(result$smd_results >= 0))
  })
})

test_that("tied scores and singleton clusters still give deterministic bijections", {
  with_scnpm_single_thread({
    tied <- scnpm_match_fixture(rep(0, 6), rep(0, 6), rep(0:2, each = 2))
    repeated <- scnpm_match_fixture(rep(0, 6), rep(0, 6), rep(0:2, each = 2))
    expect_scnpm_bijection(tied, 3L)
    expect_equal(as.numeric(tied$smd_results), rep(0, 3))
    expect_identical(tied$matching_result_df, repeated$matching_result_df)
    singleton <- scnpm_match_fixture(rep(0, 3), rep(0, 3), 0:2)
    expect_scnpm_bijection(singleton, 3L)
    expect_equal(as.numeric(singleton$smd_results), rep(0, 3))
  })
})


scnpm_public_fixture <- function() {
  expression <- rbind(A = seq_len(16),
                      B = rep(c(3, 7, 2, 8), 4),
                      C = rep(c(5, 2, 9, 4), 4),
                      D = rep(c(1, 6, 8, 3), 4))
  colnames(expression) <- paste0("cell", seq_len(16))
  expression[1, c(11, 12)] <- 0
  list(
    inference_mode = "coupled",
    goal_gene_mt = expression,
    cell_t = cbind(treatment = rep(c(1, 0), 8), quality = seq_len(16)),
    K = 1L, ds_seed = 711L, num_iter = 1L, num_save = 1L,
    epsilon_theta = 1e-10, num_step_theta = 1L,
    epsilon_lam = 0.01, num_step_lam = 1L,
    ppi_prior_low = 0.05, ppi_prior_high = 0.5
  )
}

scnpm_run_public_fixture <- function(args) {
  result <- NULL
  invisible(capture.output(
    result <- suppressMessages(do.call(run_scNPM, args))
  ))
  result
}

test_that("interleaved conditions initialize missing expression within the right group", {
  with_scnpm_single_thread({
    args <- scnpm_public_fixture()
    result <- scnpm_run_public_fixture(args)
    expression <- t(args$goal_gene_mt)
    observed <- expression > 0
    expect_equal(result$theta_est[observed], log(expression[observed]),
                 tolerance = 1e-12)
    for (cell in c(11L, 12L)) {
      condition <- args$cell_t[cell, 1]
      same_group <- args$cell_t[, 1] == condition
      nonzero <- expression[same_group & expression[, 1] > 0, 1]
      expected <- log(as.numeric(quantile(nonzero, 0.05)))
      # A near-zero HMC step exposes the actual initialization without adding
      # another production hook or assuming that a chain converges in one sweep.
      expect_equal(result$theta_est[cell, 1], expected, tolerance = 1e-8)
    }
    expect_equal(result$provenance$package_version,
                 as.character(utils::packageVersion("scNPM")))
    expect_equal(result$provenance$algorithm, "scNPM-1.0.5-coupled")
    expect_equal(result$provenance$settings$ds_seed, args$ds_seed)
    expect_equal(result$provenance$settings$num_iter, args$num_iter)
    expect_equal(result$provenance$settings$num_save, args$num_save)
    expect_equal(result$provenance$settings$ppi_prior_low, 0.05)
    expect_equal(result$provenance$settings$ppi_prior_high, 0.5)
  })
})

test_that("PPI input forms give the same persistent prior after edge cleanup", {
  with_scnpm_single_thread({
    ppi <- data.frame(node1 = c("A", "B", "C", "unknown"),
                      node2 = c("B", "A", "C", "D"),
                      combined_score = c(0.8, 0.9, 0.7, 0.8))
    text <- paste(c("#node1\tnode2\tcombined_score", "", "A\tB\t0.8",
                    "B\tA\t0.9", "C\tC\t0.7", "unknown\tD\t0.8", ""),
                  collapse = "\n")
    path <- tempfile(fileext = ".tsv")
    on.exit(unlink(path), add = TRUE)
    writeLines(text, path)
    results <- lapply(list(ppi, text, path), function(input) {
      args <- scnpm_public_fixture()
      args$protein_input <- input
      suppressWarnings(scnpm_run_public_fixture(args))
    })
    expected <- matrix(0, 4, 4, dimnames = list(LETTERS[1:4], LETTERS[1:4]))
    expected["A", "B"] <- expected["B", "A"] <- 1
    for (result in results) {
      expect_equal(unname(result$protein_adjacency), unname(expected))
      expect_equal(result$edge_prior["A", "B"], 0.5)
      expect_equal(result$edge_prior["A", "C"], 0.05)
      expect_equal(result$edge_prior["D", "C"], 0.05)
      expect_equal(unname(diag(result$edge_prior)), rep(0, 4))
    }
    expect_equal(results[[1]]$edge_prior, results[[2]]$edge_prior)
    expect_equal(results[[1]]$edge_prior, results[[3]]$edge_prior)
  })
})


test_that("all-zero cells and empty mixture components still receive updates", {
  with_scnpm_single_thread({
    all_zero <- scnpm_tiny_chain_args()
    all_zero$ind_zero[2, ] <- TRUE
    all_zero$theta_t[2, ] <- 0
    result <- scnpm_run_tiny_chain(all_zero)
    expect_true(all(is.finite(result$theta_est[2, ])))
    expect_gt(max(abs(result$theta_est[2, ])), 1e-6)

    empty <- scnpm_tiny_chain_args(num_iter = 1L, num_save = 1L)
    empty$K <- 2L
    empty$mu_treat <- cbind(empty$mu_treat, 0)
    empty$mu_ctrl <- cbind(empty$mu_ctrl, 0)
    empty$invcov_treat <- empty$cov_treat <- array(diag(3), c(3, 3, 2))
    empty$invcov_ctrl <- empty$cov_ctrl <- array(diag(3), c(3, 3, 2))
    empty$edge_treat <- empty$edge_ctrl <- array(0, c(3, 3, 2))
    empty$pi_treat <- empty$pi_ctrl <- c(0.999, 0.001)
    empty$gam <- c(1, 1)
    # All initial assignments are to component 1; component 2 has no data.
    result <- scnpm_run_tiny_chain(empty)
    for (side in c("treat", "ctrl")) {
      expect_gt(max(abs(result[[paste0("mu_", side, "_est")]][, 2])), 1e-8)
      precision <- result[[paste0("invcov_", side, "_est")]][, , 2]
      expect_gt(max(abs(precision - diag(3))), 1e-8)
      expect_gt(min(eigen(precision, symmetric = TRUE,
                          only.values = TRUE)$values), 0)
      weights <- result[[paste0("pi_", side, "_iter")]][, 1]
      expect_true(all(weights > 0))
      expect_equal(sum(weights), 1, tolerance = 1e-12)
    }
  })
})

test_that("quadratic-form SMD distinguishes singleton, empty and separated pairs", {
  exact <- list(
    K = 4L,
    matching_result_df = data.frame(treatment_cluster = 1:4,
                                    matched_control_cluster = c(1, 2, 3, NA)),
    qf_global = c(1, 2, 2, 1, 3, 3),
    final_group_label = c(0, 1, 1, 4, 5, 5)
  )
  result <- compute_qf_smd(exact)
  expect_equal(attr(result, "mode"), "exact")
  expect_equal(result$qf_smd, c(0, Inf, NA, NA))

  approximate <- list(
    K = 1L,
    matching_result_df = data.frame(treatment_cluster = 1,
                                    matched_control_cluster = 1),
    theta_est = matrix(c(1, 3, 2, 4), ncol = 1),
    invcov_treat_est = array(1, c(1, 1, 1)),
    invcov_ctrl_est = array(1, c(1, 1, 1)),
    group_est = c(1, 1, 2, 2)
  )
  result <- compute_qf_smd(approximate)
  expect_equal(attr(result, "mode"), "posterior-approx")
  expect_equal(result$qf_smd, round(5 / sqrt((32 + 72) / 2), 4))
  expect_error(compute_qf_smd(list()), "matching_result_df")
  approximate$group_est <- c(1, 2)
  expect_error(compute_qf_smd(approximate), "same number of rows")
})

test_that("invalid expression, sampling and prior inputs fail before inference", {
  args <- scnpm_public_fixture()
  bad <- args
  bad$goal_gene_mt[1, 1] <- -1
  expect_error(do.call(run_scNPM, bad), "nonnegative")
  bad$goal_gene_mt[1, 1] <- NA_real_
  expect_error(do.call(run_scNPM, bad), "finite")
  bad <- args
  bad$cell_t[, 1] <- 1
  expect_error(do.call(run_scNPM, bad), "Both treatment and control")
  bad <- args
  bad$num_save <- bad$num_iter + 1L
  expect_error(do.call(run_scNPM, bad), "num_save")
  bad <- args
  bad$ppi_prior_low <- 0
  expect_error(do.call(run_scNPM, bad), "ppi_prior_low")
  bad <- args
  bad$ppi_prior_high <- bad$ppi_prior_low / 2
  expect_error(do.call(run_scNPM, bad), "ppi_prior_high")
  probability <- scnpm_internal("calculate_edge_probability")
  expect_error(probability(0, 0.1, 1, 0), "xi")
  expect_error(probability(0, 0.1, 1, 1), "xi")
})


test_that("public fits reproduce from the same seed without modifying expression", {
  with_scnpm_single_thread({
    args <- scnpm_public_fixture()
    args$num_iter <- 8L
    args$num_save <- 4L
    original_expression <- unserialize(serialize(args$goal_gene_mt, NULL))
    first <- scnpm_run_public_fixture(args)
    second <- scnpm_run_public_fixture(args)
    expect_identical(first$Result, second$Result)
    expect_identical(first$matching_result_df, second$matching_result_df)
    expect_identical(first$edge_diff_networks, second$edge_diff_networks)
    expect_identical(args$goal_gene_mt, original_expression)
  })
})

test_that("the low-level sampler does not mutate caller-owned initial states", {
  with_scnpm_single_thread({
    args <- scnpm_tiny_chain_args(num_iter = 3L, num_save = 2L)
    diag(args$edge_prior) <- 0.7  # Ignored diagonals must not be cleared in place.
    original <- unserialize(serialize(args, NULL))
    set.seed(519L)
    invisible(capture.output(
      do.call(scnpm_internal("MCMC_with_difficulty_based_matching"), args)
    ))
    expect_identical(args, original)
  })
})

