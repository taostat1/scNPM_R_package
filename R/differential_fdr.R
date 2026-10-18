#' Select Differential Edges by Posterior Directional Error
#'
#' Uses the final final-model Rao-Blackwell probabilities from a staged
#' run. The selection unit is one treatment cluster and unordered gene pair;
#' all K clusters form a single family. The most supported direction under
#' the frozen map is reported. For `all_matchings`, its error is bounded by
#' the largest posterior error over all K possible control counterparts.
#' This is a conservative bound for every fixed bijection, not a posterior
#' distribution over matchings and not a frequentist FDR guarantee.
#'
#' The selected set has mean estimated directional error at most
#' `fdr_level`. Monte Carlo approximation, label identification and model
#' assumptions still apply. Edge states refer to continuous spike/slab
#' allocations, not exact zero precision entries. Batch MCSE is a diagnostic,
#' not a confidence bound or convergence certificate.
#'
#' @param result Output of a staged [run_scNPM()] fit.
#' @param fdr_level Target posterior expected directional FDP in (0, 1).
#' @param scope `frozen` (default) uses the selected matching. `all_matchings`
#'   bounds errors over every control counterpart for sensitivity analysis.
#' @param min_draws Minimum retained iterations before producing discoveries,
#'   at least 20, default 200. This guard does not establish convergence.
#' @param rule `delta_threshold` scans strict absolute contrast thresholds without
#'   splitting tied scores; `posterior_error` uses the version 1.0.4 error prefix.
#' @return A list with threshold (strict |Delta P| > threshold), threshold_curve,
#'   rule, and the full edge table, selected rows, estimated error
#'   bound, target, scope, sample count and status. Each table row has an
#'   `error_bound` for that edge and a `frozen_error` for its reported pairing.
#'   `estimated_error_bound` is the mean bound across selected rows.
#' @export
select_differential_edges <- function(result, fdr_level = 0.05,
    scope = c("frozen", "all_matchings"), min_draws = 200L,
    rule = c("delta_threshold", "posterior_error")) {
  scope <- match.arg(scope)
  rule <- match.arg(rule)
  scalar <- function(x, lo, hi, integer=FALSE) is.numeric(x) && length(x)==1L &&
    is.finite(x) && x>=lo && x<=hi && (!integer || x==floor(x))
  if (!scalar(fdr_level,0,1) || fdr_level==0 || fdr_level==1 ||
      !scalar(min_draws,20,.Machine$integer.max,TRUE)) stop("Invalid FDR level or min_draws.")
  if (!identical(result$inference_mode,"staged") || is.null(result$sampling_schedule) ||
      !isTRUE(result$sampling_schedule$feedback_disabled_for_retained))
    stop("Posterior error selection requires staged final sampling with feedback disabled.")
  K <- result$K; G <- result$G; joint <- result$joint_edge_probabilities
  if (!scalar(K,1,1e4,TRUE) || !scalar(G,2,1e5,TRUE)) stop("Invalid G or K.")
  wanted <- c(G,G,K*K)
  for (nm in c("p00","p10","p01","p11")) {
    a <- joint[[nm]]
    if (!is.array(a) || !identical(as.integer(dim(a)),as.integer(wanted)) ||
        any(!is.finite(a)) || any(a < -1e-10 | a > 1+1e-10))
      stop("Invalid joint edge probabilities: ",nm)
  }
  # Native loop order is i first, j second; retain it for batch diagnostics.
  ij <- do.call(rbind,lapply(seq_len(G-1L),function(i) cbind(i,seq.int(i+1L,G))))
  ind <- ij[,1]+(ij[,2]-1L)*G
  total <- joint$p00+joint$p10+joint$p01+joint$p11
  if (any(abs(matrix(total,G*G,K*K)[ind,,drop=FALSE]-1)>1e-7))
    stop("Joint probabilities must sum to one off diagonal.")
  m <- result$matching_result_df
  if (nrow(m)!=K || !identical(as.integer(m$treatment_cluster),seq_len(K)) ||
      anyNA(m$matched_control_cluster) ||
      !setequal(m$matched_control_cluster,seq_len(K))) stop("A complete frozen bijection is required.")
  names_g <- result$gene_names
  if (length(names_g)!=G || anyNA(names_g) || anyDuplicated(names_g)) stop("Invalid gene names.")
  saved <- result$sampling_schedule$retained_draws
  if (!scalar(saved,1,.Machine$integer.max,TRUE)) stop("Invalid retained draw count.")
  prob10 <- matrix(joint$p10,G*G,K*K)[ind,,drop=FALSE]
  prob01 <- matrix(joint$p01,G*G,K*K)[ind,,drop=FALSE]
  blocks <- lapply(seq_len(K),function(k) {
    c <- m$matched_control_cluster[k]
    h <- (k-1L)*K+c
    plus <- prob10[,h]; minus <- prob01[,h]
    direction <- ifelse(plus>minus,"treatment",ifelse(minus>plus,"control","undetermined"))
    frozen_error <- 1-pmax(plus,minus)
    bound <- frozen_error
    if (scope=="all_matchings") {
      hs <- (k-1L)*K+seq_len(K)
      least_plus <- apply(prob10[,hs,drop=FALSE],1,min)
      least_minus <- apply(prob01[,hs,drop=FALSE],1,min)
      bound <- 1-ifelse(direction=="treatment",least_plus,least_minus)
    }
    mcse <- rep(NA_real_,nrow(ij))
    batch <- joint$batch_p10
    if (is.array(batch) && length(dim(batch))==3L && dim(batch)[3]>=2L) {
      for (e in seq_len(nrow(ij))) {
        if (direction[e]=="undetermined") next
        source <- if (direction[e]=="treatment") joint$batch_p10 else joint$batch_p01
        col <- h
        if (scope=="all_matchings") {
          vals <- if (direction[e]=="treatment") prob10[e,(k-1L)*K+seq_len(K)] else prob01[e,(k-1L)*K+seq_len(K)]
          col <- (k-1L)*K+which.min(vals)
        }
        mcse[e] <- stats::sd(source[e,col,])/sqrt(dim(source)[3])
      }
    }
    data.frame(treatment_cluster=k,matched_control_cluster=c,
      gene_i=ij[,1],gene_j=ij[,2],from=names_g[ij[,1]],to=names_g[ij[,2]],
      p_treatment_only=plus,p_control_only=minus,p_same=1-plus-minus,
      direction=direction,delta_p=plus-minus,abs_delta_p=abs(plus-minus),
      frozen_error=frozen_error,error_bound=pmin(1,pmax(0,bound)),
      batch_mcse=mcse,stringsAsFactors=FALSE)
  })
  tab <- do.call(rbind,blocks); rownames(tab) <- NULL
  tab$selected <- FALSE; tab$prefix_mean_error <- NA_real_
  eligible <- which(tab$direction!="undetermined")
  threshold <- if(rule=="delta_threshold") Inf else NA_real_
  curve <- data.frame(threshold=numeric(),min_included_abs_delta=numeric(),
    discoveries=integer(),mean_error=numeric(),mean_frozen_error=numeric(),
    feasible=logical(),chosen=logical())
  if(length(eligible) && rule=="delta_threshold") {
    ord <- eligible[order(-tab$abs_delta_p[eligible],tab$treatment_cluster[eligible],
                         tab$gene_i[eligible],tab$gene_j[eligible])]
    scores <- tab$abs_delta_p[ord]
    ends <- c(which(diff(scores)!=0),length(ord))
    means <- cumsum(tab$error_bound[ord])/seq_along(ord)
    frozen_means <- cumsum(tab$frozen_error[ord])/seq_along(ord)
    # The next lower score is a canonical strict cutoff; 0 includes every
    # nonzero contrast. No rounding tolerance may split an exact score tie.
    cuts <- c(scores[ends[-length(ends)]+1L],0)
    curve <- data.frame(threshold=cuts,min_included_abs_delta=scores[ends],
      discoveries=ends,mean_error=means[ends],mean_frozen_error=frozen_means[ends],
      feasible=means[ends]<=fdr_level,chosen=FALSE)
    # Error means need not be monotone in contrast order: scan every cutoff.
    feasible <- which(curve$feasible)
    if(length(feasible) && saved>=min_draws) {
      pick <- feasible[length(feasible)]
      threshold <- curve$threshold[pick]; curve$chosen[pick] <- TRUE
      tab$selected <- tab$abs_delta_p>threshold
    }
    # Each row receives its whole tie-block's prefix mean, not a split prefix.
    tab$prefix_mean_error[ord] <- rep(means[ends],diff(c(0L,ends)))
  } else if(length(eligible)) {
    ord <- eligible[order(tab$error_bound[eligible],tab$treatment_cluster[eligible],
                         tab$gene_i[eligible],tab$gene_j[eligible])]
    avg <- cumsum(tab$error_bound[ord])/seq_along(ord)
    tab$prefix_mean_error[ord] <- avg
    r <- if(any(avg<=fdr_level)) max(which(avg<=fdr_level)) else 0L
    if(r && saved>=min_draws) tab$selected[ord[seq_len(r)]] <- TRUE
  }
  selected <- tab[tab$selected,,drop=FALSE]
  list(table=tab,selected=selected,target=fdr_level,scope=scope,
    rule=rule,threshold=threshold,threshold_operator=if(rule=="delta_threshold") ">" else NA_character_,
    threshold_curve=curve,
    estimated_error_bound=if(nrow(selected)) mean(selected$error_bound) else 0,
    estimated_frozen_error=if(nrow(selected)) mean(selected$frozen_error) else 0,
    status=if(saved<min_draws) "insufficient_draws" else "estimated",
    retained_draws=saved,min_draws=min_draws,
    family="all treatment clusters x unordered gene pairs; one direction per item",
    interpretation="Monte Carlo estimate of posterior expected directional FDP under the selected final model",
    frequentist_guarantee=FALSE)
}
