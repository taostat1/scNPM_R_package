#' Select Group Networks by Posterior Expected False Discovery Proportion
#'
#' Reconstructs marginal Rao-Blackwell slab probabilities from final-stage joint
#' probabilities: p10 + p11 for treatment and p01 + p11 for control. Each
#' condition and subtype is a separate family of G*(G-1)/2 unordered pairs.
#' Probabilities are sorted decreasingly, complete tie blocks are evaluated,
#' and the largest prefix with mean error 1-p at most fdr_level is retained.
#' Selection uses p >= cutoff; an empty set has cutoff Inf and BFDR zero.
#'
#' Selection is performed on the full gene panel before any subnetwork is
#' extracted. Returned masks should be shared by plotting and sensitivity
#' calculations. BFDR applies to each original network, not to a subsequent
#' union over subtypes. Edge states are continuous spike/slab allocations.
#' This estimates posterior expected FDP under the baseline model; it is not
#' a uniform frequentist FDR guarantee. Existing staged 1.0.5 fits with saved
#' joint probabilities can be processed without rerunning the sampler.
#'
#' @param result Output of a staged [run_scNPM()] fit with joint probabilities.
#' @param fdr_level BFDR target in (0,1), default 0.05, separately for each network.
#' @param min_draws Minimum retained iterations, at least 20, default 200.
#'   Below this value, probabilities and curves are returned with empty masks.
#' @return A list with `networks`, `summary`, `table`, `selected`, `probabilities`,
#'   `masks`, `target`, `status` and draw counts. Each named network (for example
#'   `treatment_1`) contains `probability_matrix`, `mask`, `cutoff`,
#'   `estimated_bfdr`, `threshold_curve`, `table` and `selected`.
#'   `probabilities` and `masks` each contain treatment/control G x G x K arrays.
#'   `probability_matrix` carries a saved selection attribute for plotting.
#' @export
#' @examples
#' # Synthetic retained posterior summaries for one matched pair.
#' # These values illustrate selection; no MCMC is run in this example.
#' gn <- paste0("G", 1:3)
#' p <- matrix(c(0,.99,.97,.99,0,.60,.97,.60,0), 3,
#'             dimnames = list(gn, gn))
#' z <- array(0, c(3,3,1))
#' p10 <- array(p, c(3,3,1)); p00 <- array(1-p, c(3,3,1))
#' diag(p00[, , 1]) <- 0
#' example_fit <- list(
#'   K = 1L, G = 3L, gene_names = gn, inference_mode = "staged",
#'   matching_result_df = data.frame(treatment_cluster=1L, matched_control_cluster=1L),
#'   sampling_schedule = list(feedback_disabled_for_retained=TRUE, retained_draws=500L),
#'   joint_edge_probabilities = list(p00=p00, p10=p10, p01=z, p11=z)
#' )
#' group <- select_group_edges(example_fit)
#' group$summary
select_group_edges <- function(result, fdr_level = 0.05, min_draws = 200L) {
  scalar <- function(x, lo, hi, integer = FALSE) is.numeric(x) && length(x) == 1L &&
    is.finite(x) && x >= lo && x <= hi && (!integer || x == floor(x))
  if (!scalar(fdr_level, 0, 1) || fdr_level == 0 || fdr_level == 1 ||
      !scalar(min_draws, 20, .Machine$integer.max, TRUE))
    stop("Invalid FDR level or min_draws.")
  if (!identical(result$inference_mode, "staged") ||
      !isTRUE(result$sampling_schedule$feedback_disabled_for_retained))
    stop("Group BFDR requires staged final sampling with feedback disabled.")
  G <- result$G; K <- result$K; gn <- result$gene_names
  saved <- result$sampling_schedule$retained_draws
  if (!scalar(G, 2, 1e5, TRUE) || !scalar(K, 1, 1e4, TRUE) ||
      !scalar(saved, 1, .Machine$integer.max, TRUE)) stop("Invalid G, K or retained draw count.")
  if (length(gn) != G || anyNA(gn) || anyDuplicated(gn) || any(!nzchar(gn)))
    stop("Invalid gene names.")
  m <- result$matching_result_df
  if (!is.data.frame(m) || nrow(m) != K ||
      !identical(as.integer(m$treatment_cluster), seq_len(K)) ||
      length(m$matched_control_cluster) != K || anyNA(m$matched_control_cluster) ||
      !setequal(m$matched_control_cluster, seq_len(K)))
    stop("A complete frozen bijection is required.")
  map <- as.integer(m$matched_control_cluster)
  ij <- which(upper.tri(matrix(FALSE, G, G)), arr.ind = TRUE)
  joint <- result$joint_edge_probabilities
  total <- array(0, c(G, G, K*K))
  for (nm in c("p00", "p10", "p01", "p11")) {
    a <- joint[[nm]]
    if (!is.numeric(a) || !is.array(a) ||
        !identical(as.integer(dim(a)), as.integer(c(G, G, K*K))) ||
        any(!is.finite(a)) || any(a < -1e-12 | a > 1+1e-12))
      stop("Invalid joint edge probabilities: ", nm)
    for (h in seq_len(K*K)) if (max(abs(a[, , h] - t(a[, , h]))) > 1e-12)
      stop("Joint edge probabilities must be symmetric.")
    total <- total + a
  }
  for (h in seq_len(K*K)) if (any(abs(total[, , h][ij] - 1) > 1e-10))
    stop("Joint probabilities must sum to one off diagonal.")
  probs <- list(treatment = array(0, c(G, G, K), dimnames = list(gn, gn, NULL)),
                control = array(0, c(G, G, K), dimnames = list(gn, gn, NULL)))
  masks <- lapply(probs, function(a) array(FALSE, dim(a), dimnames = dimnames(a)))
  networks <- list(); rows <- list(); status <- if (saved < min_draws) "insufficient_draws" else "estimated"
  counterpart_maxdiff <- 0
  for (condition in c("treatment", "control")) for (k in seq_len(K)) {
    tk <- if (condition == "treatment") k else match(k, map)
    ck <- if (condition == "control") k else map[k]
    h <- (tk-1L)*K + ck
    marginal <- function(h) {
      if (condition == "treatment") joint$p10[, , h] + joint$p11[, , h]
      else joint$p01[, , h] + joint$p11[, , h]
    }
    raw <- marginal(h)[ij]
    for (other in seq_len(K)) {
      oh <- if (condition == "treatment") (k-1L)*K + other else (other-1L)*K + k
      counterpart_maxdiff <- max(counterpart_maxdiff, max(abs(raw - marginal(oh)[ij])))
    }
    if (counterpart_maxdiff > 1e-10)
      stop("Inconsistent marginal probabilities across counterpart clusters.")
    p <- pmin(1, pmax(0, raw))
    ord <- order(-p, seq_along(p)); scores <- p[ord]
    ends <- c(which(diff(scores) != 0), length(scores))
    avg <- cumsum(1-scores)/seq_along(scores)
    feasible <- which(avg[ends] <= fdr_level)
    pick <- if (length(feasible) && saved >= min_draws) utils::tail(feasible, 1L) else 0L
    cutoff <- if (pick) scores[ends[pick]] else Inf
    keep <- p >= cutoff
    curve <- data.frame(cutoff = scores[ends], discoveries = ends,
                        estimated_bfdr = avg[ends], feasible = avg[ends] <= fdr_level,
                        chosen = seq_along(ends) == pick)
    pm <- matrix(0, G, G, dimnames = list(gn, gn)); pm[ij] <- p; pm <- pm + t(pm)
    mask <- matrix(FALSE, G, G, dimnames = list(gn, gn)); mask[ij] <- keep; mask <- mask | t(mask)
    probs[[condition]][, , k] <- pm; masks[[condition]][, , k] <- mask
    tab <- data.frame(condition = condition, cluster = k,
                      gene_i = ij[, 1], gene_j = ij[, 2], from = gn[ij[, 1]], to = gn[ij[, 2]],
                      probability = p, posterior_error = 1-p, selected = keep,
                      stringsAsFactors = FALSE)
    bfdr <- if (any(keep)) mean(1-p[keep]) else 0
    attr(pm, "scNPM_group_selection") <- list(mask = mask, method = "bayes_fdr",
       target = fdr_level, cutoff = cutoff, estimated_bfdr = bfdr, status = status,
       condition = condition, cluster = k, retained_draws = saved, min_draws = min_draws)
    name <- paste0(condition, "_", k)
    networks[[name]] <- list(condition = condition, cluster = k,
      probability_matrix = pm, mask = mask, cutoff = cutoff, threshold_operator = ">=",
      estimated_bfdr = bfdr, target = fdr_level, status = status,
      table = tab, selected = tab[keep, , drop = FALSE], threshold_curve = curve)
    rows[[name]] <- data.frame(condition = condition, cluster = k, family_size = nrow(ij),
      selected_edges = sum(keep), cutoff = cutoff, estimated_bfdr = bfdr,
      expected_false_allocations = sum(1-p[keep]), target = fdr_level,
      retained_draws = saved, status = status)
  }
  tab <- do.call(rbind, lapply(networks, function(x) x$table)); rownames(tab) <- NULL
  summary <- do.call(rbind, rows); rownames(summary) <- NULL
  list(networks = networks, summary = summary, table = tab,
       selected = tab[tab$selected, , drop = FALSE], probabilities = probs, masks = masks,
       target = fdr_level, status = status, retained_draws = saved, min_draws = min_draws,
       counterpart_max_difference = counterpart_maxdiff,
       family = "Each condition and subtype separately, over all unordered gene pairs",
       probability_source = "Rao-Blackwell marginals: treatment p10+p11; control p01+p11",
       interpretation = "Monte Carlo estimate of posterior expected FDP under the baseline model",
       frequentist_guarantee = FALSE)
}
