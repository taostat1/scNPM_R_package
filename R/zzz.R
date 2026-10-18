#' Single-Cell Network Propensity Matching
#'
#' Estimate and compare heterogeneous gene conditional-dependence networks
#' across two conditions. A latent Gaussian mixture and probit observation
#' mechanism represent cellular heterogeneity and zero observations. Protein
#' interaction information supplies edge-specific prior weights, and network
#' propensity matching guides stochastic updates of uncertain edges.
#'
#' @details
#' Matching uses centered quadratic forms against all current treatment and
#' control reference networks (2 * K features). Pairwise propensity regression
#' uses the observed condition indicator, balances the two conditions' total
#' weights and standardizes predictors within each training fold. Global scores
#' remain available with pooled scaling. Default staged inference freezes a
#' stable complete matching or selects a recent consensus map at the adaptation
#' limit. It then disables feedback and discards a post-switch burn-in before retaining final-model
#' draws. Staged fits use a four-state paired prior by default, with optional
#' independent priors for compatibility. Paired precision updates and joint
#' relabeling preserve frozen partners. No gene modules are introduced.
#' Group networks use per-subtype BFDR selection from marginal slab
#' probabilities. Directional edge selection controls an estimated posterior expected
#' FDP, with a conservative all-counterparts option; it is not a uniform
#' frequentist FDR guarantee. Coupled mode with global scores and greedy
#' assignment retains the earlier continuous-feedback workflow.
#'
#' The augmented criterion trace records likelihood, network-prior and SMD
#' contributions. It is diagnostic and need not increase. HMC uses fixed step
#' sizes, leapfrog counts and Metropolis acceptance correction.
#'
#' @seealso \code{\link{run_scNPM}}, \code{\link{compute_qf_smd}}
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @importClassesFrom Matrix dgTMatrix
#' @useDynLib scNPM, .registration = TRUE
## usethis namespace: end
NULL

# ggplot2 uses non-standard evaluation for aesthetics; declare the variables
# used inside aes() so R CMD check does not flag "no visible binding".
# ".data" is the rlang pronoun used for dynamic aesthetic mapping.
utils::globalVariables(c(
  ".data",
  "xcoord", "ycoord",
  "xcoord.from", "ycoord.from", "xcoord.to", "ycoord.to",
  "size_norm", "color_norm", "gene"
))
