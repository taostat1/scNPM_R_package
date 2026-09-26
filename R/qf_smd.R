#' Compute the SMD of the Network Quadratic Form for Matched Cluster Pairs
#'
#' In all-reference mode (version 1.0.3), this helper returns one row per
#' matched pair and reference network, using the exact final network features.
#' The additional column `network_feature` identifies each reference.
#' The legacy vector interface below remains supported for older results.
#'
#' The legacy network quadratic form \eqn{qf_i = \theta_i' \Omega_{k(i)} \theta_i}
#' enters the propensity score model as an additional covariate. This helper
#' reports the matched-pair SMD of the quadratic form using the same formula
#' as `smd_covariate_results`:
#' \deqn{SMD = |mean_t - mean_c| / sqrt((var_t + var_c) / 2)}
#'
#' Two modes are supported (returned in the attribute `mode`):
#' \itemize{
#'   \item **exact**: `result_list` contains `qf_global` and
#'   `final_group_label` (scNPM >= 1.0.1). The helper uses the quadratic
#'   forms and labels from the internal final matching call.
#'   \item **posterior-approx**: for results saved by older versions (with the
#'   heavy `Result` element stripped), the quadratic form is recomputed from
#'   stored sampling summaries (`theta_est`, `invcov_treat_est`,
#'   `invcov_ctrl_est`, `group_est`). The legacy mode name is retained for
#'   compatibility; these summaries need not approximate the final matching
#'   state closely.
#' }
#'
#' @param result_list Output of [run_scNPM()] (as saved with [save()]).
#'
#' @return A data frame with columns `treat_cluster`,
#'   `matched_control_cluster` and `qf_smd`, rounded to four decimal places.
#'   All-reference results additionally include `network_feature` and contain
#'   K * (2 * K) rows; their mode attribute is `exact-all-reference`.
#'   Unmatched or empty pairs have NA. Singleton samples have variance zero.
#'   At the numerical zero-variance boundary, equal means give zero and
#'   unequal means give infinity, using the matching routine's tolerances.
#' @export
#'
#' @examples
#' # Small synthetic example (exact mode), runs in well under 5 seconds
#' set.seed(1)
#' G <- 5; K <- 2; nT <- 30; nC <- 30
#' invT <- array(diag(G), c(G, G, K)); invC <- array(diag(G), c(G, G, K))
#' thT <- matrix(rnorm(nT * G), nT, G)
#' thC <- matrix(rnorm(nC * G, mean = 0.5), nC, G)
#' qfT <- rowSums((thT %*% invT[, , 1]) * thT)
#' qfC <- rowSums((thC %*% invC[, , 2]) * thC)
#' rl <- list(
#'   K = K,
#'   matching_result_df = data.frame(treatment_cluster = 1:K,
#'                                   matched_control_cluster = c(2, 1)),
#'   qf_global = c(qfT, qfC),
#'   final_group_label = c(rep(0, nT), rep(K + 1, nC))
#' )
#' compute_qf_smd(rl)
compute_qf_smd <- function(result_list) {
  if (is.null(result_list$matching_result_df)) {
    stop("'result_list' must be the output of run_scNPM() (matching_result_df missing).")
  }
  if (identical(result_list$network_covariate_mode, "all_reference")) {
    features <- result_list$network_covariates
    if (!is.matrix(features) || is.null(result_list$final_group_label) ||
        nrow(features) != length(result_list$final_group_label) ||
        ncol(features) != 2L * result_list$K || any(!is.finite(features))) {
      stop("All-reference diagnostics require the final network_covariates matrix and final_group_label.")
    }
    feature_names <- colnames(features)
    if (is.null(feature_names)) {
      feature_names <- c(paste0("qf_ref_T", seq_len(result_list$K)),
                         paste0("qf_ref_C", seq_len(result_list$K)))
    }
    out <- lapply(seq_len(ncol(features)), function(j) {
      legacy <- list(K = result_list$K,
                     matching_result_df = result_list$matching_result_df,
                     qf_global = features[, j],
                     final_group_label = result_list$final_group_label)
      x <- compute_qf_smd(legacy)
      x$network_feature <- feature_names[j]
      x[, c("treat_cluster", "matched_control_cluster", "network_feature", "qf_smd")]
    })
    out <- do.call(rbind, out)
    rownames(out) <- NULL
    attr(out, "mode") <- "exact-all-reference"
    return(out)
  }
  K <- result_list$K
  mdf <- result_list$matching_result_df
  n_pairs <- nrow(mdf)
  smd_out <- rep(NA_real_, n_pairs)

  has_exact <- !is.null(result_list$qf_global) &&
               !is.null(result_list$final_group_label)

  if (has_exact) {
    # ---- exact mode: reuse the qf and the final cluster labels of the last
    #      matching call (labels: treat 1..K stored as 0..K-1, ctrl stored as
    #      K..2K-1 in final_group_label) ----
    qf <- as.numeric(result_list$qf_global)
    gl <- as.numeric(result_list$final_group_label)
    is_treat <- gl < K
    grp_t <- gl[is_treat] + 1        # 1-based treatment cluster
    grp_c <- gl[!is_treat] - K + 1   # 1-based control cluster
    qf_t <- qf[is_treat]
    qf_c <- qf[!is_treat]
    mode <- "exact"
  } else {
    # ---- legacy posterior-approx mode: recompute qf from stored sampling
    #      summaries (mean theta and precision matrices, modal labels) ----
    theta <- as.matrix(result_list$theta_est)   # N x G (cells x genes)
    invT  <- result_list$invcov_treat_est        # G x G x K
    invC  <- result_list$invcov_ctrl_est
    gN    <- as.numeric(result_list$group_est)   # N-length, 1-based:
    # treat cells 1..K, control cells K+1..2K (indexes rows of theta_est)
    if (length(gN) != nrow(theta)) {
      stop("posterior-approx mode requires 'group_est' (N-length) with the ",
           "same number of rows as 'theta_est'.")
    }
    qf_all <- numeric(nrow(theta))
    for (k in 1:K) {
      idx_t <- which(gN == k)
      if (length(idx_t) > 0) {
        Tk <- theta[idx_t, , drop = FALSE]
        qf_all[idx_t] <- rowSums((Tk %*% invT[, , k]) * Tk)
      }
      idx_c <- which(gN == K + k)
      if (length(idx_c) > 0) {
        Ck <- theta[idx_c, , drop = FALSE]
        qf_all[idx_c] <- rowSums((Ck %*% invC[, , k]) * Ck)
      }
    }
    if (any(!is.finite(qf_all))) stop("Non-finite network quadratic form in stored summaries.")
    is_treat <- gN >= 1 & gN <= K
    qf_t  <- qf_all[is_treat]
    qf_c  <- qf_all[!is_treat]
    grp_t <- gN[is_treat]
    grp_c <- gN[!is_treat] - K
    mode <- "posterior-approx"
  }

  for (r in seq_len(n_pairs)) {
    kt <- mdf$treatment_cluster[r]
    kc <- mdf$matched_control_cluster[r]
    if (is.na(kc)) next  # unmatched pair
    vt <- qf_t[grp_t == kt]
    vc <- qf_c[grp_c == kc]
    if (length(vt) == 0L || length(vc) == 0L) next
    m1 <- mean(vt); m2 <- mean(vc)
    v1 <- if (length(vt) < 2L) 0 else stats::var(vt)
    v2 <- if (length(vc) < 2L) 0 else stats::var(vc)
    vs <- max(0, (v1 + v2) / 2)
    difference <- abs(m1 - m2)
    smd_out[r] <- if (!is.finite(difference) || !is.finite(vs)) {
      Inf
    } else if (vs < 1e-12) {
      if (difference < 1e-10) 0 else Inf
    } else {
      difference / sqrt(vs)
    }
  }

  out <- data.frame(
    treat_cluster = mdf$treatment_cluster,
    matched_control_cluster = mdf$matched_control_cluster,
    qf_smd = round(smd_out, 4)
  )
  attr(out, "mode") <- mode
  out
}
