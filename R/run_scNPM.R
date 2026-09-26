#' Run Single-Cell Network Propensity Matching (scNPM)
#'
#' Estimate condition-specific gene networks using staged stochastic inference.
#' The default first adapts networks and matching, aligns component labels to
#' fixed initial mean anchors, and freezes a complete mapping after sustained
#' stability. It then disables counterpart feedback, completes post-switch
#' burn-in and retains final baseline-model samples. Failure to stabilize by
#' the adaptation limit stops the run without producing final discoveries.
#'
#' Network features are centered distances to all 2*K current references,
#' with pooled scaling; the legacy uncentered feature remains available.
#' Default differential selection bounds the estimated posterior directional
#' error over every possible control counterpart. This is model-based Bayesian
#' error selection, not a finite-sample frequentist FDR guarantee. The coupled
#' mode reproduces version 1.0.3 and uses descriptive frequency contrasts.
#' The J_aug trace uses an augmented-data likelihood, a state-wise network
#' prior kernel and total SMD. It is a diagnostic, not an optimized objective
#' or the observed-data marginal likelihood; stochastic decreases are allowed.
#'
#' @param goal_gene_mt Finite nonnegative numeric matrix or data frame, with
#'   genes in rows and cells in columns. At least two genes are required.
#'   Row names must be unique gene identifiers; missing names are generated.
#' @param cell_t Finite numeric matrix or data frame, with cells in the same
#'   order as columns of `goal_gene_mt`. The first column is 0 (control) or
#'   1 (treatment); both groups must be present. Remaining columns are numeric
#'   covariates, standardized across both groups. Encode categorical covariates
#'   numerically before calling this function.
#' @param K Positive integer number of clusters in each condition, no larger
#'   than the number of distinct expression profiles in either condition.
#' @param ds_seed Integer seed between 0 and 2147483646. Reproducibility requires
#'   the same inputs, software environment and OpenMP thread count.
#' @param num_iter Positive integer total number of sampler iterations.
#' @param num_save Positive integer number of final iterations retained,
#'   no greater than `num_iter`.
#' @param protein_input Optional PPI data frame, local TSV path, or TSV text.
#'   Required columns are `node1`, `node2` and `combined_score`. A leading `#`
#'   on the header is accepted. Positive scores indicate an interaction;
#'   duplicate pairs are combined, self-interactions and genes outside the
#'   expression matrix are ignored. `NULL` assigns the lower prior to all edges.
#' @param ssp_v0 Positive spike standard deviation, smaller than `ssp_v1`.
#' @param ssp_v1 Positive slab standard deviation.
#' @param ssp_l Positive diagonal shrinkage parameter; the exponential-prior
#'   rate for precision diagonals is `ssp_l / 2`.
#' @param edge_prob_threshold Confidence threshold in (0.5, 1) for deciding
#'   whether an edge receives matching feedback; this is separate from the
#'   differential-network cutoff.
#' @param edge_prob_adj Feedback multiplier in (0, 1]. A value of 1 leaves
#'   baseline edge probabilities unchanged by matching feedback.
#' @param epsilon_theta Positive fixed HMC step size for latent expression.
#' @param num_step_theta Positive integer number of HMC steps for latent expression.
#' @param epsilon_lam Positive fixed HMC step size for dropout parameters.
#' @param num_step_lam Positive integer number of HMC steps for dropout parameters.
#' @param ppi_prior_low Baseline inclusion probability for edges without a PPI
#'   interaction (default 0.05).
#' @param ppi_prior_high Baseline inclusion probability for PPI-supported edges
#'   (default 0.5). Requires `0 < ppi_prior_low <= ppi_prior_high <= 1 - 1e-6`.
#'
#' @param network_covariates Network features in the propensity model.
#'   `"all_reference"` (default) evaluates every cell against all current
#'   treatment networks followed by all control networks, using
#'   \eqn{(\theta_i-\mu_r)^T\Omega_r(\theta_i-\mu_r)}.
#'   `"legacy"` uses the single uncentered quadratic form from version 1.0.2.
#'   Other cell covariates are retained; matching_method selects the assignment rule.
#' @param propensity_ridge Positive diagonal ridge stabilization for the
#'   pooled logistic IRLS fit (default 1e-6, also applied to the intercept).
#'   This is a numerical setting, not an empirically tuned penalty.
#'
#' @param inference_mode `staged` (default) or `coupled` (continuous feedback; use greedy matching for version 1.0.3 behavior).
#' @param adapt_min Earliest freezing iteration. NULL uses min(100, floor(adapt_max/2)), at least 2.
#' @param adapt_max Latest adaptation iteration. NULL uses num_iter-num_save-final_burnin.
#' @param freeze_patience Consecutive identical complete maps needed to freeze,
#'   at least 2. NULL uses min(20, floor(adapt_max/2)), at least 2.
#' @param final_burnin Minimum iterations after disabling feedback and before
#'   retention. NULL uses floor((num_iter-num_save)/3), at least 1.
#' @param diagnostic_every Positive interval for final-stage J_aug diagnostics;
#'   adaptation diagnostics and matching histories are recorded every iteration.
#' @param differential_method `auto` chooses `bayes_fdr` for staged fits and
#'   descriptive `frequency` selection for coupled fits. Explicit `frequency`
#'   is available for staged comparisons; `bayes_fdr` requires staged sampling.
#' @param fdr_level Posterior directional error target in (0,1), default 0.05.
#' @param fdr_scope `frozen` (default) uses the frozen correspondence.
#'   `all_matchings` uses a conservative bound over all control counterparts.
#' @param fdr_min_draws Minimum retained draws for FDR discoveries, at least 20,
#'   default 200. Fewer draws return probabilities and no selected discoveries.
#'
#' @param matching_method `lexicographic` (default) minimizes maximum pair SMD,
#'   then total SMD, with deterministic ties and no SMD cap. `greedy` reproduces
#'   the previous assignment rule. A frozen map is not reoptimized.
#' @param fdr_rule `delta_threshold` (default) scans all distinct absolute
#'   contrast cutoffs and retains the largest complete threshold set whose mean
#'   posterior directional error is at most fdr_level. `posterior_error` retains
#'   the version 1.0.4 error-sorted prefix rule.
#'
#' @return A named list containing all key analysis results (can be saved with `save()`):
#' \item{runtime}{MCMC runtime in minutes (difftime object)}
#' \item{K}{Number of cell clusters used in analysis}
#' \item{G}{Number of genes (rows in goal_gene_mt)}
#' \item{N}{Number of cells (columns in goal_gene_mt)}
#' \item{lam0_est}{Retained-draw means of lambda0 parameter (length = G)}
#' \item{lam1_est}{Retained-draw means of lambda1 parameter (length = G)}
#' \item{mu_treat_est}{Retained-draw summaries for treatment group (G x K matrix)}
#' \item{mu_ctrl_est}{Retained-draw summaries for control group (G x K matrix)}
#' \item{invcov_treat_est}{Retained-draw precision summaries for treatment group (G x G x K array)}
#' \item{invcov_ctrl_est}{Retained-draw precision summaries for control group (G x G x K array)}
#' \item{edge_treat_est}{Retained-draw edge frequencies for treatment group (G x G x K array)}
#' \item{edge_ctrl_est}{Retained-draw edge frequencies for control group (G x G x K array)}
#' \item{group_treat_est}{Cluster labels for treatment cells (length = number of treatment cells)}
#' \item{group_ctrl_est}{Cluster labels for control cells (length = number of control cells)}
#' \item{pi_treat_est}{Retained-draw mixing-weight means for treatment group (length = K)}
#' \item{pi_ctrl_est}{Retained-draw mixing-weight means for control group (length = K)}
#' \item{theta_treat_est}{Retained-draw latent-expression summaries for treatment group}
#' \item{theta_ctrl_est}{Retained-draw latent-expression summaries for control group}
#' \item{theta_est}{Retained-draw latent-expression summaries (cells x genes)}
#' \item{matching_result_df}{Cluster matching table (treatment_cluster vs matched_control_cluster)}
#' \item{matching_method, smd_cost_matrix, frozen_smd_cost_matrix, matching_objective}{Assignment rule, all K by K final-state pair costs, and selected maximum/total SMD. Frozen maps are retained even if final costs change.}
#' \item{smd_results}{Standardized Mean Difference (SMD) values for each treatment cluster (computed on the propensity score)}
#' \item{smd_covariate_results}{Per-covariate SMD matrix for matched cluster pairs (rows = treatment clusters, columns = the QC covariates in cell_t[, -1] plus the selected network features; NA for unmatched pairs). A value below 0.1 is a conventional balance diagnostic}
#' \item{match_map}{Named vector mapping treatment clusters to matched control clusters}
#' \item{qf_global}{Legacy uncentered per-cell quadratic form retained for compatibility. It enters the propensity model only in legacy mode.}
#' \item{network_covariates}{Raw network-feature matrix from the final matching call. Columns are qf_ref_T1 through qf_ref_TK then qf_ref_C1 through qf_ref_CK in all-reference mode, or qf_network in legacy mode.}
#' \item{network_covariate_mode}{Selected network feature mode.}
#' \item{propensity_scaling}{Pooled centers and scales for the full predictor matrix, with other cell covariates before network features.}
#' \item{propensity_scores}{Per-cell propensity scores from the final matching call}
#' \item{final_group_label}{Final MCMC cluster label per cell (0-based: treatment cells 0..K-1, control cells K..2K-1; used by compute_qf_smd())}
#' \item{edge_diff_networks}{Per-cluster networks selected by differential_method; solid/dashed indicates treatment/control enrichment. The matrix uses Rao-Blackwell contrasts for bayes_fdr and empirical frequency contrasts for frequency.}
#' \item{inference_mode, differential_method}{Resolved inference and selection modes.}
#' \item{sampling_schedule}{Freeze time, actual post-switch burn-in and retained interval.}
#' \item{stage_history, matching_history}{Per-iteration phase, feedback status and aligned mapping history.}
#' \item{criterion_trace}{Recorded augmented likelihood, network prior kernel, SMD loss and J_aug.}
#' \item{joint_edge_probabilities}{Final all-pair 00/10/01/11 Rao-Blackwell probabilities, empirical 10/01 frequencies and batch means; slice (k-1)*K+c.}
#' \item{differential_fdr}{Full directional error table, selected rows, estimated bound, scope and status; NULL for frequency mode.}
#' \item{Result}{Retained sampler output, including parameter draws and selected trace summaries}
#'
#' \item{protein_adjacency}{Symmetric binary PPI adjacency matrix; diagonal entries are zero.}
#' \item{edge_prior}{Per-edge baseline inclusion probabilities used at every edge update. Diagonal entries are unused.}
#' \item{provenance}{Package version, input fingerprint, native-library fingerprint, execution environment, times and analysis settings.}
#'
#' @examples
#' \donttest{
#' set.seed(7)
#' x <- matrix(rpois(6 * 40, lambda = 4), nrow = 6,
#'             dimnames = list(paste0("G", 1:6), NULL))
#' cells <- cbind(condition = rep(0:1, each = 20))
#' ppi <- data.frame(node1 = "G1", node2 = "G2", combined_score = 0.9)
#' # Short demonstration; increase the run length for an analysis.
#' fit <- run_scNPM(x, cells, K = 2, num_iter = 20, num_save = 10,
#'                  inference_mode = "coupled",
#'                  protein_input = ppi)
#' fit$matching_result_df
#' compute_qf_smd(fit)
#' }
#' @export
#' @importFrom stats kmeans quantile rnorm setNames
run_scNPM <- function(
    goal_gene_mt,
    cell_t,
    K,
    ds_seed = 20250102,
    num_iter = 1000,
    num_save = 500,
    protein_input = NULL,
    ssp_v0 = 0.02,
    ssp_v1 = 1,
    ssp_l = 1,
    edge_prob_threshold = 0.6,
    edge_prob_adj = 0.75,
    epsilon_theta = 0.2,
    num_step_theta = 20,
    epsilon_lam = 0.01,
    num_step_lam = 10,
    ppi_prior_low = 0.05,
    ppi_prior_high = 0.5,
    network_covariates = c("all_reference", "legacy"),
    propensity_ridge = 1e-6,
    inference_mode = c("staged", "coupled"),
    adapt_min = NULL, adapt_max = NULL, freeze_patience = NULL,
    final_burnin = NULL, diagnostic_every = 10L,
    differential_method = c("auto", "bayes_fdr", "frequency"),
    fdr_level = 0.05, fdr_scope = c("frozen", "all_matchings"),
    fdr_min_draws = 200L,
    matching_method = c("lexicographic", "greedy"),
    fdr_rule = c("delta_threshold", "posterior_error")
) {
  network_covariates <- match.arg(network_covariates)
  inference_mode <- match.arg(inference_mode)
  differential_method <- match.arg(differential_method)
  fdr_scope <- match.arg(fdr_scope)
  matching_method <- match.arg(matching_method)
  fdr_rule <- match.arg(fdr_rule)
  if (differential_method=="auto") differential_method <- if(inference_mode=="staged") "bayes_fdr" else "frequency"
  if (inference_mode=="coupled" && differential_method=="bayes_fdr")
    stop("bayes_fdr requires staged final baseline sampling; coupled frequencies are not posterior probabilities.")
  if (!is.matrix(goal_gene_mt) && !is.data.frame(goal_gene_mt)) {
    stop("'goal_gene_mt' must be a matrix or data frame (rows=genes, columns=cells).")
  }
  X <- as.matrix(goal_gene_mt)
  G <- nrow(X)
  N <- ncol(X)
  if (G < 2L || N < 2L) stop("'goal_gene_mt' requires at least two genes and two cells.")
  if (!is.numeric(X) || any(!is.finite(X)) || any(X < 0)) {
    stop("'goal_gene_mt' must contain finite nonnegative numeric expression values.")
  }
  gene_names <- rownames(X)
  if (is.null(gene_names)) {
    warning("'goal_gene_mt' has no row names; generating Gene_1, Gene_2, ...")
    gene_names <- paste0("Gene_", seq_len(G))
    rownames(X) <- gene_names
  }
  if (anyNA(gene_names) || any(!nzchar(trimws(gene_names))) || anyDuplicated(gene_names)) {
    stop("Gene row names must be nonempty and unique.")
  }
  if (!is.matrix(cell_t) && !is.data.frame(cell_t)) {
    stop("'cell_t' must be a numeric matrix or data frame.")
  }
  cell_t <- as.matrix(cell_t)
  if (nrow(cell_t) != N || ncol(cell_t) < 1L) {
    stop("'cell_t' must have one row per cell and at least one column.")
  }
  if (!is.numeric(cell_t) || any(!is.finite(cell_t))) {
    stop("'cell_t' must contain finite numeric values; encode categorical covariates first.")
  }
  treatment_col <- cell_t[, 1L]
  if (!all(treatment_col %in% c(0, 1))) {
    stop("First column of 'cell_t' must contain only 0 (control) and 1 (treatment).")
  }
  treat_cell_idx <- which(treatment_col == 1)
  ctrl_cell_idx <- which(treatment_col == 0)
  if (!length(treat_cell_idx) || !length(ctrl_cell_idx)) {
    stop("Both treatment and control cells are required.")
  }
  check_scalar <- function(x, name, lower, upper = Inf, integer = FALSE,
                           lower_open = FALSE, upper_open = FALSE) {
    ok <- is.numeric(x) && length(x) == 1L && is.finite(x)
    if (ok) ok <- (if (lower_open) x > lower else x >= lower) &&
      (if (upper_open) x < upper else x <= upper) &&
      (!integer || x == floor(x))
    if (!ok) stop("Invalid '", name, "'; see help('run_scNPM') for its allowed range.")
  }
  check_scalar(K, "K", 1, min(length(treat_cell_idx), length(ctrl_cell_idx)), TRUE)
  check_scalar(ds_seed, "ds_seed", 0, .Machine$integer.max - 1, TRUE)
  check_scalar(num_iter, "num_iter", 1, .Machine$integer.max, TRUE)
  check_scalar(num_save, "num_save", 1, num_iter, TRUE)
  check_scalar(diagnostic_every,"diagnostic_every",1,.Machine$integer.max,TRUE)
  check_scalar(fdr_level,"fdr_level",0,1,lower_open=TRUE,upper_open=TRUE)
  check_scalar(fdr_min_draws,"fdr_min_draws",20,.Machine$integer.max,TRUE)
  if (inference_mode=="staged") {
    available <- num_iter-num_save
    if (available<3) stop("Staged inference requires at least three pre-retention iterations; increase num_iter or use inference_mode='coupled'.")
    if (is.null(final_burnin)) final_burnin <- max(1L,floor(available/3))
    check_scalar(final_burnin,"final_burnin",1,available-2L,TRUE)
    if (is.null(adapt_max)) adapt_max <- available-final_burnin
    check_scalar(adapt_max,"adapt_max",2,available-final_burnin,TRUE)
    if (is.null(adapt_min)) adapt_min <- max(2L,min(100L,floor(adapt_max/2)))
    if (is.null(freeze_patience)) freeze_patience <- max(2L,min(20L,floor(adapt_max/2)))
    check_scalar(adapt_min,"adapt_min",1,adapt_max,TRUE)
    check_scalar(freeze_patience,"freeze_patience",2,adapt_max,TRUE)
  } else {
    adapt_min <- adapt_max <- final_burnin <- 1L; freeze_patience <- 2L
  }
  check_scalar(ssp_v0, "ssp_v0", 0, lower_open = TRUE)
  check_scalar(ssp_v1, "ssp_v1", ssp_v0, lower_open = TRUE)
  check_scalar(ssp_l, "ssp_l", 0, lower_open = TRUE)
  check_scalar(edge_prob_threshold, "edge_prob_threshold", 0.5, 1,
               lower_open = TRUE, upper_open = TRUE)
  check_scalar(edge_prob_adj, "edge_prob_adj", 0, 1, lower_open = TRUE)
  check_scalar(epsilon_theta, "epsilon_theta", 0, lower_open = TRUE)
  check_scalar(epsilon_lam, "epsilon_lam", 0, lower_open = TRUE)
  check_scalar(num_step_theta, "num_step_theta", 1, .Machine$integer.max, TRUE)
  check_scalar(num_step_lam, "num_step_lam", 1, .Machine$integer.max, TRUE)
  check_scalar(ppi_prior_low, "ppi_prior_low", 0, 1 - 1e-6, lower_open = TRUE)
  check_scalar(ppi_prior_high, "ppi_prior_high", ppi_prior_low, 1 - 1e-6)
  check_scalar(propensity_ridge, "propensity_ridge", 0, lower_open = TRUE)
  if (is.null(colnames(cell_t))) {
    colnames(cell_t) <- c("condition", if (ncol(cell_t) > 1L)
      paste0("covariate_", seq_len(ncol(cell_t) - 1L)))
  }
  raw_cell_t <- cell_t
  if (ncol(cell_t) > 1L) {
    for (i in seq.int(2L, ncol(cell_t))) {
      cov_sd <- stats::sd(cell_t[, i])
      cell_t[, i] <- if (cov_sd > 0) (cell_t[, i] - mean(cell_t[, i])) / cov_sd else 0
    }
  }

  protein_mat <- .scnpm_ppi_adjacency(protein_input, gene_names)
  edge_prior <- ppi_prior_low + (ppi_prior_high - ppi_prior_low) * protein_mat
  diag(edge_prior) <- 0
  input_fingerprint <- .scnpm_fingerprint(list(expression = X, cells = raw_cell_t,
                                               protein_adjacency = protein_mat))

  # --------------------------
  # 4. Independent K-means clustering
  # --------------------------
  set.seed(ds_seed)
  X_log <- log1p(X)  # Log-transform expression matrix for better clustering

  # Clustering for treatment group (cells x genes matrix)
  X_treat_log <- t(X_log[, treat_cell_idx, drop = FALSE])
  group_treat <- .scnpm_initial_clusters(X_treat_log, K, "treatment")

  # Clustering for control group (different seed to avoid correlation)
  set.seed(ds_seed + 1)
  X_ctrl_log <- t(X_log[, ctrl_cell_idx, drop = FALSE])
  group_ctrl <- .scnpm_initial_clusters(X_ctrl_log, K, "control")

  # --------------------------
  # 5. Theta_t initialization
  # --------------------------
  ind_zero <- (X == 0)  # Zero value indicator matrix
  theta_t <- X          # Initialize theta_t with original expression values

  # Stratified zero imputation (by treatment/control and cluster)
  for (g in 1:G) {
    ind_zero_g <- ind_zero[g, ]
    X_g <- X[g, ]

    # Impute treatment group cells
    for (k in 1:K) {
      treat_k_idx <- treat_cell_idx[which(group_treat == k)]
      if (length(treat_k_idx) == 0) next

      X_gk <- X_g[treat_k_idx]
      ind_zero_gk <- ind_zero_g[treat_k_idx]

      # Impute zeros with 5th percentile of non-zero values (if available)
      if (sum(!ind_zero_gk) > 0) {
        theta_gk <- X_gk
        theta_gk[ind_zero_gk] <- quantile(X_gk[!ind_zero_gk], na.rm = TRUE, probs = 0.05)
        theta_t[g, treat_k_idx] <- theta_gk
      }
    }

    # Impute control group cells
    for (k in 1:K) {
      ctrl_k_idx <- ctrl_cell_idx[which(group_ctrl == k)]
      if (length(ctrl_k_idx) == 0) next

      X_gk <- X_g[ctrl_k_idx]
      ind_zero_gk <- ind_zero_g[ctrl_k_idx]

      # Impute zeros with 5th percentile of non-zero values (if available)
      if (sum(!ind_zero_gk) > 0) {
        theta_gk <- X_gk
        theta_gk[ind_zero_gk] <- quantile(X_gk[!ind_zero_gk], na.rm = TRUE, probs = 0.05)
        theta_t[g, ctrl_k_idx] <- theta_gk
      }
    }
  }

  # Log transformation and NA/Inf handling
  theta_t <- log(theta_t)
  theta_t[is.na(theta_t) | theta_t == -Inf] <- 0  # Replace NA/-Inf with 0
  message("Theta_t initialization completed (stratified zero imputation)")

  # --------------------------
  # 6. Initialize mean/covariance/precision matrices
  # --------------------------
  # Treatment group matrices
  mu_treat <- matrix(0, nrow = G, ncol = K)
  cov_treat <- array(diag(1, G), dim = c(G, G, K))
  invcov_treat <- array(diag(1, G), dim = c(G, G, K))

  for (k in 1:K) {
    treat_k_idx <- treat_cell_idx[which(group_treat == k)]
    if (length(treat_k_idx) == 0) next

    theta_k <- theta_t[, treat_k_idx, drop = FALSE]
    mu_treat[, k] <- rowMeans(theta_k)  # Cluster-specific mean

    # Gene-specific variance (diagonal covariance matrix)
    gene_var <- if (ncol(theta_k) > 1L) {
      rowSums((theta_k - mu_treat[, k])^2) / (ncol(theta_k) - 1L)
    } else rep(1, G)
    gene_var[!is.finite(gene_var) | gene_var <= 0] <- 1
    diag(cov_treat[, , k]) <- gene_var
    diag(invcov_treat[, , k]) <- 1.0 / gene_var  # Precision = inverse variance
  }

  # Control group matrices
  mu_ctrl <- matrix(0, nrow = G, ncol = K)
  cov_ctrl <- array(diag(1, G), dim = c(G, G, K))
  invcov_ctrl <- array(diag(1, G), dim = c(G, G, K))

  for (k in 1:K) {
    ctrl_k_idx <- ctrl_cell_idx[which(group_ctrl == k)]
    if (length(ctrl_k_idx) == 0) next

    theta_k <- theta_t[, ctrl_k_idx, drop = FALSE]
    mu_ctrl[, k] <- rowMeans(theta_k)  # Cluster-specific mean

    # Gene-specific variance (diagonal covariance matrix)
    gene_var <- if (ncol(theta_k) > 1L) {
      rowSums((theta_k - mu_ctrl[, k])^2) / (ncol(theta_k) - 1L)
    } else rep(1, G)
    gene_var[!is.finite(gene_var) | gene_var <= 0] <- 1
    diag(cov_ctrl[, , k]) <- gene_var
    diag(invcov_ctrl[, , k]) <- 1.0 / gene_var  # Precision = inverse variance
  }

  # --------------------------
  # 7. Initialize edge matrices
  # --------------------------
  edge_treat <- array(0, dim = c(G, G, K))
  edge_ctrl <- array(0, dim = c(G, G, K))

  # Assign protein matrix to all clusters
  for (k in 1:K) {
    edge_treat[, , k] <- protein_mat
    edge_ctrl[, , k] <- protein_mat
  }
  message("Edge matrices initialized (based on protein interaction matrix)")

  # --------------------------
  # 8. Final parameter preparation
  # --------------------------
  theta_t <- t(theta_t)       # Transpose to cells x genes (C++ input format)
  ind_zero <- t(ind_zero)     # Match theta_t dimensions
  gam <- rep(1, K)            # Dirichlet prior weights for cell type probabilities

  # Cell type probability initialization
  pi_treat <- update_pi_R(group_treat - 1L, gam, K)
  pi_ctrl <- update_pi_R(group_ctrl - 1L, gam, K)

  # Spike-slab prior parameters for lambda
  lambda0_t <- rnorm(G, mean = 3, sd = 0.2)
  lambda1_t <- rnorm(G, mean = -2, sd = 0.2)
  lambda1_t <- -abs(lambda1_t)
  ssp_xi <- ppi_prior_low  # Scalar fallback for the internal interface

  # --------------------------
  # 9. Run C++ MCMC function
  # --------------------------
  # The parallel sampler seeds its workers from R's RNG. Thread count and
  # software environment must also be fixed for reproducibility.
  set.seed(ds_seed)
  start_time <- Sys.time()
  Result <- MCMC_with_difficulty_based_matching(
    num_iter = num_iter,
    num_save = num_save,
    theta_t = theta_t,
    ind_zero = ind_zero,
    mu_treat = mu_treat,
    mu_ctrl = mu_ctrl,
    invcov_treat = invcov_treat,
    cov_treat = cov_treat,
    invcov_ctrl = invcov_ctrl,
    cov_ctrl = cov_ctrl,
    edge_treat = edge_treat,
    edge_ctrl = edge_ctrl,
    group_treat = group_treat,
    group_ctrl = group_ctrl,
    lambda0_t = lambda0_t,
    lambda1_t = lambda1_t,
    pi_treat = pi_treat,
    pi_ctrl = pi_ctrl,
    gam = gam,
    cell_t = cell_t,
    G = G,
    N_int = N,
    K = K,
    ssp_v0 = ssp_v0,
    ssp_v1 = ssp_v1,
    ssp_l = ssp_l,
    ssp_xi = ssp_xi,
    epsilon_theta = epsilon_theta,
    num_step_theta = num_step_theta,
    eta_mu = 0,
    tau_sq_mu = 1,
    lam0_0 = 3,
    lam1_0 = -2,
    sigma2_lam0 = 0.2,
    sigma2_lam1 = 0.2,
    edge_prob_threshold = edge_prob_threshold,
    edge_prob_adj = edge_prob_adj,
    epsilon_lam = epsilon_lam,
    num_step_lam = num_step_lam,
    edge_prior = edge_prior,
    network_covariates = network_covariates,
    propensity_ridge = propensity_ridge,
    inference_mode = inference_mode, adapt_min = adapt_min, adapt_max = adapt_max,
    freeze_patience = freeze_patience, final_burnin = final_burnin,
    diagnostic_every = diagnostic_every, matching_method = matching_method
  )

  # Calculate MCMC runtime
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "mins")
  message(paste0("MCMC Runtime: ", round(runtime, 2), " minutes"))

  # --------------------------
  # 10. Extract retained-draw summaries
  # --------------------------
  # Retained-draw summaries
  lam0_est <- rowMeans(Result$lam0_iter)
  lam1_est <- rowMeans(Result$lam1_iter)
  mu_treat_est <- Result$mu_treat_est
  mu_ctrl_est <- Result$mu_ctrl_est
  invcov_treat_est <- Result$invcov_treat_est
  invcov_ctrl_est <- Result$invcov_ctrl_est
  edge_treat_est <- Result$edge_treat_est
  edge_ctrl_est <- Result$edge_ctrl_est
  group_treat_est <- Result$group_treat_est
  group_ctrl_est <- Result$group_ctrl_est
  group_est <- Result$group_est
  pi_treat_est <- rowMeans(Result$pi_treat_iter)
  pi_ctrl_est <- rowMeans(Result$pi_ctrl_iter)
  theta_treat_est <- Result$theta_treat_est
  theta_ctrl_est <- Result$theta_ctrl_est
  theta_est <- Result$theta_est

  # Cluster matching results
  matching_result_df <- Result$matching_results$matching_result_df
  smd_results <- Result$matching_results$smd_results

  # Final-state balance diagnostics use exactly the selected network features.
  network_names <- if (network_covariates == "all_reference") {
    c(paste0("qf_ref_T", seq_len(K)), paste0("qf_ref_C", seq_len(K)))
  } else "qf_network"
  network_values <- as.matrix(Result$matching_results$network_covariates)
  colnames(network_values) <- network_names
  predictor_names <- c(colnames(cell_t)[-1], network_names)
  propensity_scaling <- Result$matching_results$propensity_scaling
  if (!is.null(propensity_scaling)) {
    propensity_scaling$center <- setNames(as.numeric(propensity_scaling$center), predictor_names)
    propensity_scaling$scale <- setNames(as.numeric(propensity_scaling$scale), predictor_names)
  }
  smd_covariate_results <- Result$matching_results$smd_covariate_results
  if (!is.null(smd_covariate_results) && is.matrix(smd_covariate_results) &&
      ncol(smd_covariate_results) > 0) {
    cov_names <- predictor_names
    colnames(smd_covariate_results) <- cov_names[seq_len(ncol(smd_covariate_results))]
    rownames(smd_covariate_results) <- paste0("treat_cluster_",
                                              seq_len(nrow(smd_covariate_results)))
  }

  # Create cluster matching map (treatment → control)
  match_map <- setNames(
    matching_result_df$matched_control_cluster,
    matching_result_df$treatment_cluster
  )
  match_map <- match_map[!is.na(match_map)]  # Remove unmatched clusters

  stage_history <- as.data.frame(Result$stage_history)
  names(stage_history) <- c("iteration","stage_code","feedback_active","frozen_after_iteration","stable_count")
  stage_history$stage <- c("coupled","adaptation","final_burnin","retained")[stage_history$stage_code+1L]
  matching_history <- Result$matching_history
  colnames(matching_history) <- paste0("T",seq_len(K))
  criterion_trace <- as.data.frame(Result$criterion_trace)
  names(criterion_trace) <- c("iteration","log_augmented_likelihood","network_prior_kernel","matching_loss","J_aug","stage_code")
  criterion_trace <- criterion_trace[!is.na(criterion_trace$iteration),,drop=FALSE]
  schedule <- list(frozen_at=if(inference_mode=="staged") Result$frozen_at else NA_integer_,
    adapt_min=adapt_min,adapt_max=adapt_max,freeze_patience=freeze_patience,
    minimum_final_burnin=final_burnin,
    actual_final_burnin=if(inference_mode=="staged") num_iter-num_save-Result$frozen_at else NA_integer_,
    retained_start=num_iter-num_save+1L,retained_draws=num_save,
    feedback_disabled_for_retained=inference_mode=="staged",
    labeling=if(inference_mode=="staged") "minimum squared-mean assignment to fixed initial anchors" else "legacy")
  fdr <- NULL
  if (differential_method=="bayes_fdr") {
    fdr <- select_differential_edges(list(K=K,G=G,gene_names=gene_names,
      matching_result_df=matching_result_df,inference_mode=inference_mode,
      sampling_schedule=schedule,joint_edge_probabilities=Result$joint_edge_probabilities),
      fdr_level=fdr_level,scope=fdr_scope,min_draws=fdr_min_draws,rule=fdr_rule)
    if (fdr$status=="insufficient_draws")
      warning("Fewer than fdr_min_draws retained iterations: probabilities returned, but no FDR discoveries selected.")
  }

  # --------------------------
  # 11. Per-cluster differential network analysis
  # --------------------------
  edge_diff_networks <- list()

  # Iterate over each matched treatment-control cluster pair
  for (treat_clust in names(match_map)) {
    treat_clust <- as.numeric(treat_clust)
    ctrl_clust <- as.numeric(match_map[as.character(treat_clust)])

    # Extract cluster-specific edge matrices
    edge_treat_clust <- edge_treat_est[, , treat_clust]
    edge_ctrl_clust <- edge_ctrl_est[, , ctrl_clust]

    # Calculate differential network (treatment - control)
    edge_diff <- edge_treat_clust - edge_ctrl_clust
    selected_directions <- matrix(0,G,G)
    if (!is.null(fdr)) {
      h <- (treat_clust-1L)*K+ctrl_clust
      edge_diff <- Result$joint_edge_probabilities$p10[,,h]-Result$joint_edge_probabilities$p01[,,h]
      hits <- fdr$selected[fdr$selected$treatment_cluster==treat_clust,,drop=FALSE]
      if (nrow(hits)) selected_directions[cbind(hits$gene_i,hits$gene_j)] <-
        ifelse(hits$direction=="treatment",1,-1)
    }
    dimnames(edge_diff) <- list(gene_names, gene_names)

    # Initialize solid/dashed edge matrices and edge lists
    solid_mat <- matrix(0, nrow = G, ncol = G, dimnames = list(gene_names, gene_names))  # Solid edges (difference > 0.5)
    dashed_mat <- matrix(0, nrow = G, ncol = G, dimnames = list(gene_names, gene_names)) # Dashed edges (difference < -0.5)
    solid_edges <- data.frame(from = character(), to = character(), value = numeric(), stringsAsFactors = FALSE)
    dashed_edges <- data.frame(from = character(), to = character(), value = numeric(), stringsAsFactors = FALSE)

    # Identify solid/dashed edges (upper triangle only for undirected network)
    for (i in seq_len(G - 1L)) {
      for (j in seq.int(i + 1L, G)) {
        val <- edge_diff[i, j]

        # Solid edge (positive difference > 0.5)
        if (if(is.null(fdr)) val > 0.5 else selected_directions[i,j]==1) {
          solid_mat[i, j] <- 1
          solid_mat[j, i] <- 1  # Symmetric for undirected network
          solid_edges <- rbind(solid_edges, data.frame(
            from = gene_names[i],
            to = gene_names[j],
            value = val,
            stringsAsFactors = FALSE
          ))
        }
        # Dashed edge (negative difference < -0.5)
        else if (if(is.null(fdr)) val < -0.5 else selected_directions[i,j]==-1) {
          dashed_mat[i, j] <- 1
          dashed_mat[j, i] <- 1  # Symmetric for undirected network
          dashed_edges <- rbind(dashed_edges, data.frame(
            from = gene_names[i],
            to = gene_names[j],
            value = val,
            stringsAsFactors = FALSE
          ))
        }
      }
    }

    # Carry the selection with the contrast so plotting cannot reintroduce
    # rejected edges through the legacy amplitude thresholds.
    if (!is.null(fdr)) attr(edge_diff,"scNPM_selection") <- list(
      solid=solid_mat,dashed=dashed_mat,method="bayes_fdr",
      scope=fdr$scope,target=fdr$target,rule=fdr$rule,threshold=fdr$threshold)

    # Store per-cluster differential network results
    edge_diff_networks[[paste0("treat_", treat_clust, "_ctrl_", ctrl_clust)]] <- list(
      pair_info = data.frame(
        treatment_cluster = treat_clust,
        matched_control_cluster = ctrl_clust,
        smd_value = smd_results[treat_clust]
      ),
      diff_matrix = edge_diff,
      solid_network = list(
        mat = solid_mat,
        edges = solid_edges,
        count = nrow(solid_edges)
      ),
      dashed_network = list(
        mat = dashed_mat,
        edges = dashed_edges,
        count = nrow(dashed_edges)
      )
    )
  }

  # --------------------------
  # 12. Compile final results list
  # --------------------------
  result_list <- list(
    # Basic metadata
    inference_mode = inference_mode,
    differential_method = differential_method,
    matching_method = matching_method,
    smd_cost_matrix = Result$matching_results$smd_cost_matrix,
    frozen_smd_cost_matrix = Result$frozen_smd_cost_matrix,
    matching_objective = list(max_smd=Result$matching_results$matching_max_smd,
      total_smd=Result$matching_results$matching_total_smd,
      frozen=Result$matching_results$matching_is_frozen),
    sampling_schedule = schedule,
    stage_history = stage_history,
    matching_history = matching_history,
    criterion_trace = criterion_trace,
    joint_edge_probabilities = Result$joint_edge_probabilities,
    differential_fdr = fdr,
    runtime = runtime,
    K = K,
    G = G,
    N = N,
    # Retained-draw parameter summaries
    lam0_est = lam0_est,
    lam1_est = lam1_est,
    gene_names = gene_names,
    mu_treat_est = mu_treat_est,
    mu_ctrl_est = mu_ctrl_est,
    invcov_treat_est = invcov_treat_est,
    invcov_ctrl_est = invcov_ctrl_est,
    edge_treat_est = edge_treat_est,
    edge_ctrl_est = edge_ctrl_est,
    group_treat_est = group_treat_est,
    group_ctrl_est = group_ctrl_est,
    group_est = group_est,
    pi_treat_est = pi_treat_est,
    pi_ctrl_est = pi_ctrl_est,
    theta_treat_est = theta_treat_est,
    theta_ctrl_est = theta_ctrl_est,
    theta_est = theta_est,
    # Cluster matching results
    matching_result_df = matching_result_df,
    smd_results = smd_results,
    smd_covariate_results = smd_covariate_results,
    match_map = match_map,
    # Propensity-matching internals (lightweight copies of the final matching
    # call, kept even if the heavy `Result` element is stripped before saving;
    # used by compute_qf_smd() for exact post-hoc qf SMD recomputation)
    qf_global = Result$matching_results$quadratic_form_global,
    network_covariates = network_values,
    network_covariate_mode = network_covariates,
    propensity_scaling = propensity_scaling,
    propensity_scores = Result$matching_results$propensity_scores,
    final_group_label = Result$group_iter[, ncol(Result$group_iter)],
    # Per-cluster differential networks
    edge_diff_networks = edge_diff_networks,
    protein_adjacency = protein_mat,
    edge_prior = edge_prior,
    provenance = list(
      package_version = as.character(utils::packageVersion("scNPM")),
      algorithm = paste0("scNPM-1.0.5-",inference_mode),
      R_version = R.version.string,
      platform = R.version$platform,
      threads = sampler_thread_info(),
      started_at = format(start_time, tz = "UTC", usetz = TRUE),
      finished_at = format(end_time, tz = "UTC", usetz = TRUE),
      input_md5 = input_fingerprint,
      native_library_md5 = .scnpm_native_fingerprint(),
      settings = list(
        ds_seed = ds_seed, num_iter = num_iter, num_save = num_save, K = K,
        ssp_v0 = ssp_v0, ssp_v1 = ssp_v1, ssp_l = ssp_l,
        ppi_prior_low = ppi_prior_low, ppi_prior_high = ppi_prior_high,
        edge_prob_threshold = edge_prob_threshold, edge_prob_adj = edge_prob_adj,
        epsilon_theta = epsilon_theta, num_step_theta = num_step_theta,
        epsilon_lam = epsilon_lam, num_step_lam = num_step_lam,
        eta_mu = 0, tau_sq_mu = 1, lam0_0 = 3, lam1_0 = -2,
        sigma2_lam0 = 0.2, sigma2_lam1 = 0.2, gam = rep(1, K),
        differential_cutoff = if(differential_method=="frequency") 0.5 else fdr$threshold, covariate_scaling = "pooled",
        network_covariates = network_covariates,
        n_network_covariates = length(network_names),
        network_feature_names = network_names,
        propensity_ridge = propensity_ridge,
        reference_update = "current-matching-state",
        quadratic_form = if (network_covariates == "all_reference")
          "centered-shared-2K-references" else "uncentered",
        continuous_feedback = inference_mode=="coupled",
        inference_mode = inference_mode, differential_method = differential_method,
        adapt_min=adapt_min,adapt_max=adapt_max,freeze_patience=freeze_patience,
        final_burnin=final_burnin,diagnostic_every=diagnostic_every,
        fdr_level=fdr_level,fdr_scope=fdr_scope,fdr_min_draws=fdr_min_draws,
        fdr_rule=fdr_rule,matching_method=matching_method,smd_cap=NULL
      )
    ),
    # Retained sampler output, including selected trace summaries
    Result = Result
  )

  # Return complete results list
  return(result_list)
}


.scnpm_initial_clusters <- function(x, K, condition) {
  if (nrow(unique(x)) < K) {
    stop("'K' exceeds the number of distinct expression profiles in ", condition, ".")
  }
  if (K == 1L) return(rep.int(1L, nrow(x)))
  if (K == nrow(x)) return(seq_len(K))
  kmeans(x, centers = K, nstart = 5)$cluster
}

.scnpm_ppi_adjacency <- function(protein_input, gene_names) {
  G <- length(gene_names)
  adj <- matrix(0, G, G, dimnames = list(gene_names, gene_names))
  if (is.null(protein_input)) return(adj)
  required <- c("node1", "node2", "combined_score")
  if (is.data.frame(protein_input)) {
    dat <- protein_input
    names(dat) <- sub("^#", "", trimws(names(dat)))
  } else if (is.character(protein_input) && length(protein_input) && !anyNA(protein_input)) {
    if (length(protein_input) == 1L && file.exists(protein_input)) {
      lines <- readLines(protein_input, warn = FALSE, encoding = "UTF-8")
    } else {
      lines <- unlist(strsplit(protein_input, "\n", fixed = TRUE), use.names = FALSE)
    }
    lines <- sub("\r$", "", lines)
    lines <- lines[nzchar(trimws(lines))]
    candidates <- sub("^\\s*#\\s*", "", lines)
    header <- which(vapply(strsplit(candidates, "\t", fixed = TRUE),
                            function(z) all(required %in% trimws(z)), logical(1)))
    if (!length(header)) stop("PPI TSV requires columns node1, node2 and combined_score.")
    first <- header[1L]
    tail_lines <- if (first < length(lines)) lines[seq.int(first + 1L, length(lines))] else character()
    tail_lines <- tail_lines[!grepl("^\\s*#", tail_lines)]
    dat <- utils::read.delim(text = paste(c(candidates[first], tail_lines), collapse = "\n"),
                             check.names = FALSE, stringsAsFactors = FALSE,
                             comment.char = "", quote = "")
    names(dat) <- trimws(names(dat))
  } else {
    stop("'protein_input' must be a data frame, local TSV path or TSV text.")
  }
  if (!all(required %in% names(dat))) {
    stop("PPI input requires columns node1, node2 and combined_score.")
  }
  node1 <- trimws(as.character(dat$node1))
  node2 <- trimws(as.character(dat$node2))
  score <- suppressWarnings(as.numeric(as.character(dat$combined_score)))
  if (anyNA(node1) || anyNA(node2) || any(!nzchar(node1)) || any(!nzchar(node2)) ||
      any(!is.finite(score)) || any(score < 0)) {
    stop("PPI nodes must be nonempty and scores must be finite nonnegative numbers.")
  }
  i <- match(node1, gene_names)
  j <- match(node2, gene_names)
  keep <- !is.na(i) & !is.na(j) & i != j & score > 0
  if (any(keep)) {
    adj[cbind(i[keep], j[keep])] <- 1
    adj[cbind(j[keep], i[keep])] <- 1
  }
  adj
}

.scnpm_fingerprint <- function(object) {
  path <- tempfile("scnpm-input-", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(object, path, version = 2, compress = FALSE)
  unname(tools::md5sum(path))
}

.scnpm_native_fingerprint <- function() {
  dll <- getLoadedDLLs()[["scNPM"]]
  if (is.null(dll)) return(NA_character_)
  unname(tools::md5sum(dll[["path"]]))
}
