#include "scNPM_types.h"
#include <random>
#include <iomanip>
#include <algorithm>
#include <chrono>
#include "global_matching.h"
#include "freeze_consensus.h"

// [[Rcpp::export]]
int set_scNPM_threads(int threads) {
  if (threads < 1) Rcpp::stop("threads must be a positive integer.");
  const int previous = omp_get_max_threads();
#ifdef _OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(threads);
#else
  if (threads != 1) Rcpp::stop("This build does not support OpenMP.");
#endif
  return previous;
}

// [[Rcpp::export]]
Rcpp::List sampler_thread_info() {
#ifdef _OPENMP
  const bool enabled = true;
#else
  const bool enabled = false;
#endif
  return Rcpp::List::create(Rcpp::_["openmp"] = enabled,
                            Rcpp::_["max_threads"] = omp_get_max_threads());
}

inline double thread_safe_rnorm(std::mt19937& rng) {
  std::normal_distribution<double> norm_dist(0.0, 1.0);
  return norm_dist(rng);
}

inline double thread_safe_randu(std::mt19937& rng) {
  std::uniform_real_distribution<double> unif_dist(0.0, 1.0);
  return unif_dist(rng);
}

vec mvnrnd(const vec& mu, const mat& sigma, int n) {
  int d = mu.n_elem;
  mat Z = randn(d, n);
  mat L = chol(sigma, "lower");
  mat X = repmat(mu, 1, n) + L * Z;
  return X.col(0);
}

double U_theta(const vec &theta_miss, const vec &mu, const mat &prec, const vec &lambda0_miss, const vec &lambda1_miss) {
  vec linear_term = lambda0_miss + lambda1_miss % theta_miss;
  double U = as_scalar((theta_miss - mu).t() * prec * (theta_miss - mu)) / 2.0 -
    sum(rcpp_log_normcdf(linear_term));
  return U;
}

// [[Rcpp::export]]
double calculate_edge_probability(double omega, double v0, double v1, double xi) {
  if (!std::isfinite(omega) || !std::isfinite(v0) || !std::isfinite(v1) ||
      v0 <= 0.0 || v1 <= 0.0 || !std::isfinite(xi) || xi <= 0.0 || xi >= 1.0) {
    Rcpp::stop("Edge probability requires finite omega, positive scales, and 0 < xi < 1.");
  }
  // Slab/spike posterior odds include both prior masses, xi and 1 - xi.
  // Equal scales leave the probability equal to its prior, even in the tails.
  if (v0 == v1) return xi;
  const long double w = omega, s0 = v0, s1 = v1;
  const long double log_odds = std::log(xi) - std::log1p(-xi) +
    std::log(v0) - std::log(v1) +
    0.5L * w * w * (1.0L / (s0 * s0) - 1.0L / (s1 * s1));
  if (log_odds >= 0.0L) return static_cast<double>(1.0L / (1.0L + std::exp(-log_odds)));
  const long double e = std::exp(log_odds);
  return static_cast<double>(e / (1.0L + e));
}

inline double sample_variance_or_zero(const vec& x) {
  return x.n_elem < 2 ? 0.0 : var(x);
}

inline double standardized_mean_difference(double m1, double m2, double v1, double v2) {
  const double difference = std::abs(m1 - m2);
  const double pooled_variance = std::max(0.0, (v1 + v2) / 2.0);
  if (!std::isfinite(difference) || !std::isfinite(pooled_variance)) return datum::inf;
  if (pooled_variance < 1e-12) return difference < 1e-10 ? 0.0 : datum::inf;
  return difference / std::sqrt(pooled_variance);
}

void U_theta_grad(vec &U_grad, const vec &theta_miss, const vec &mu, const mat &prec, const vec &lambda0_miss, const vec &lambda1_miss) {
  vec linear_term = lambda0_miss + lambda1_miss % theta_miss;
  U_grad = prec * (theta_miss - mu) - 
    lambda1_miss % rcpp_normal_hazard(linear_term, true);
}

void HMC_theta(vec &theta_new, double &log_r, const vec &theta_miss, const vec &mu, const mat &prec, const vec &lambda0_miss,
               const vec &lambda1_miss, const vec &p_rnorm, const double &epsilon, const int &num_step) {
  vec U_grad;
  
  theta_new = theta_miss;
  vec p_new = p_rnorm;
  
  U_theta_grad(U_grad, theta_new, mu, prec, lambda0_miss, lambda1_miss);
  p_new -= epsilon *  U_grad / 2.0;
  
  for (int i = 0; i < (num_step-1); i++) {
    theta_new += epsilon * p_new;
    U_theta_grad(U_grad, theta_new, mu, prec, lambda0_miss, lambda1_miss);
    p_new -= epsilon * U_grad;
  }
  theta_new += epsilon * p_new;
  
  U_theta_grad(U_grad, theta_new, mu, prec, lambda0_miss, lambda1_miss);
  p_new -= epsilon * U_grad / 2.0;
  
  p_new = -p_new;
  log_r = U_theta(theta_miss, mu, prec, lambda0_miss, lambda1_miss) - U_theta(theta_new, mu, prec, lambda0_miss, lambda1_miss) +
    sum(p_rnorm % p_rnorm) / 2.0 - sum(p_new % p_new) / 2.0;
}

void update_theta(mat &theta_t, const mat &ind_zero, const mat &mu_t, const cube &invcov_t, const cube &cov_t, const vec &lambda0_t, 
                  const vec &lambda1_t, const vec &group_t, const int &N, const unsigned long &seed, const double &epsilon = 0.01, 
                  const int &num_step = 50) {
  int n_threads = std::min(omp_get_max_threads(), 24); 
  
#pragma omp parallel shared(theta_t) num_threads(n_threads)
{
  int rank = omp_get_thread_num();
  std::mt19937 rng(seed + rank * 10000);
  
#pragma omp for schedule(static)
  for (int i = 0; i < N; i++) {
    vec ind_zero_i = ind_zero.row(i).t();
    uvec ind_0 = find(ind_zero_i == true);
    
    if (ind_0.n_elem > 0) {
      vec theta_i = theta_t.row(i).t();
      uvec ind_obs = find(ind_zero_i == false);

      vec theta_i_obs = theta_i(ind_obs);
      vec theta_i_0 = theta_i(ind_0);
      
      int group_idx = static_cast<int>(group_t(i));
      if (group_idx < 0 || group_idx >= (int)invcov_t.n_slices) continue;
      
      vec mu_i = mu_t.col(static_cast<uword>(group_idx)); 
      mat cov_i = cov_t.slice(static_cast<uword>(group_idx));
      mat invcov_i = invcov_t.slice(static_cast<uword>(group_idx));
      
      if (ind_0.max() >= mu_i.n_elem ||
          (!ind_obs.is_empty() && ind_obs.max() >= mu_i.n_elem)) continue;
      const mat prec_cond = invcov_i.submat(ind_0, ind_0);
      vec mu_cond = mu_i(ind_0);
      // An all-zero cell still has a proper latent conditional and must move.
      if (!ind_obs.is_empty()) {
        mu_cond -= solve(prec_cond, invcov_i.submat(ind_0, ind_obs) *
                         (theta_i_obs - mu_i(ind_obs)), solve_opts::likely_sympd);
      }

      vec p_rnorm(ind_0.n_elem);
      for (uword t = 0; t < ind_0.n_elem; t++) {
        p_rnorm(t) = thread_safe_rnorm(rng);
      }
      double tmp_unif = thread_safe_randu(rng);
      
      vec lambda0_miss = lambda0_t(ind_0);
      vec lambda1_miss = lambda1_t(ind_0);
      
      vec theta_star_0;
      double log_r;
      HMC_theta(theta_star_0, log_r, theta_i_0, mu_cond, prec_cond, lambda0_miss, lambda1_miss, p_rnorm, epsilon, num_step);
      
      if (tmp_unif < exp(log_r)) {
        theta_i(ind_0) = theta_star_0;
        theta_t.row(i) = theta_i.t();
      }
    }
  }
}
}

void update_mu(mat &mu_t, const mat &theta_t, const cube &invcov_t, const vec &group_t, const int &G, const int &K, const double &eta_mu = 0, const double &tau_sq_mu = 1) {
  mat I_tau_sq = eye(G,G);
  I_tau_sq /= tau_sq_mu;
  
  for (int k = 0; k < K; k++) {
    uvec ind_k = find(group_t == k);  
    vec tmp1(G, fill::zeros);
    if (!ind_k.is_empty()) tmp1 = sum(theta_t.rows(ind_k), 0).t();
    mat invcov_k = invcov_t.slice(static_cast<uword>(k));
    
    mat COV = inv(ind_k.n_elem * invcov_k + I_tau_sq);
    vec MU = COV * (invcov_k * tmp1 + eta_mu/tau_sq_mu);
    
    vec mu_k = mvnrnd(MU, COV, 1);
    mu_t.col(k) = mu_k;
  }
}

void update_invcov(cube &invcov_t, cube &cov_t, const cube &edge_t, const mat &theta_t, const mat &mu_t, const vec &group_t, 
                   const double &ssp_v0, const double &ssp_v1, const double &ssp_l, const int &G, const int &K, const unsigned long &seed) {
  int n_threads = std::min(omp_get_max_threads(), K);
  
#pragma omp parallel shared(invcov_t, cov_t) num_threads(n_threads)
{
  int rank = omp_get_thread_num();
  std::mt19937 rng(seed + rank * 10000);
  
#pragma omp for schedule(static)
  for (int k = 0; k < K; k++) {
    uvec ind_k = find(group_t == k);  
    mat S(G, G, fill::zeros);
    if (!ind_k.is_empty()) {
      const mat theta_k = theta_t.rows(ind_k).t();
      const mat theta_mu = theta_k.each_col() - mu_t.col(k);
      S = theta_mu * theta_mu.t();
    }
    
    mat edge_k = edge_t.slice(static_cast<uword>(k));
    mat V = edge_k * ssp_v1 * ssp_v1;
    uvec ind_n_v1 = find(edge_k == 0);
    V(ind_n_v1).fill(ssp_v0 * ssp_v0); 
    V.diag().fill(0);
    
    mat invcov_k = invcov_t.slice(static_cast<uword>(k));
    mat cov_k = cov_t.slice(static_cast<uword>(k));
    
    vec G_vec = regspace(0, G-1);
    
    for (int g = 0; g < G; g++) {
      uvec ind_ng = find(G_vec != g);
      vec v12_g = V.col(g);
      mat v12_inv = diagmat(1/v12_g(ind_ng));
      vec cov_k_g = cov_k.col(g);
      
      mat w11_inv = cov_k.submat(ind_ng, ind_ng) - cov_k_g(ind_ng) * cov_k_g(ind_ng).t() / cov_k_g(g);
      mat C = (S(g,g) + ssp_l) * w11_inv + v12_inv;
      
      mat C_chol_inv = inv(chol(C));
      vec s12_g = S.col(g);
      vec mu_w12 = -C_chol_inv * C_chol_inv.t() * s12_g(ind_ng);
      
      vec rnorm_vec(G-1);
      for (int t = 0; t < (G-1); t++) {
        rnorm_vec(t) = thread_safe_rnorm(rng);
      }
      vec w12 = mu_w12 + C_chol_inv * rnorm_vec;
      
      // Thread-safe gamma draw (R::rgamma uses R's global RNG, which is neither
      // thread-safe nor reproducible inside OpenMP regions)
      std::gamma_distribution<double> gam_dist(
        ind_k.n_elem / 2.0 + 1.0,
        2.0 / (S(g, g) + ssp_l));
      double w_v = gam_dist(rng);
      
      vec w11_inv_w12 = w11_inv * w12;
      double w22 = w_v + as_scalar(w12.t() * w11_inv_w12); 
      
      vec w_update(G);
      w_update(ind_ng) = w12;
      w_update(g) = w22;
      
      invcov_k.col(g) = w_update;
      invcov_k.row(g) = w_update.t();
      
      cov_k.submat(ind_ng, ind_ng) = w11_inv + w11_inv_w12 * w11_inv_w12.t() / w_v;
      vec cov_k_ng = - w11_inv_w12 / w_v;
      
      vec cov_update(G);
      cov_update(ind_ng) = cov_k_ng;
      cov_update(g) = 1 / w_v;
      
      cov_k.col(g) = cov_update;
      cov_k.row(g) = cov_update.t();
    }
    
    cov_t.slice(static_cast<uword>(k)) = cov_k;
    invcov_t.slice(static_cast<uword>(k)) = invcov_k;
  }
}
}

// ---------------------------------------------------------------------------
// Logistic propensity score via IRLS with configurable ridge stabilisation.
// Every QC and network predictor uses its own pooled (all-cell) mean and sample
// standard deviation. The ridge adds to all coefficients, including intercept,
// preserving the v1.0.1/1.0.2 numerical convention at the default 1e-6.
// ---------------------------------------------------------------------------
vec fit_logistic_propensity(const mat& Xraw, const vec& y,
                            const double ridge = 1e-6,
                            vec* center_out = nullptr, vec* scale_out = nullptr) {
  const uword N = Xraw.n_rows;
  const uword p = Xraw.n_cols;
  vec prob(N, fill::value(0.5));

  if (N == 0 || y.n_elem != N || !Xraw.is_finite() || !y.is_finite() ||
      any((y != 0.0) % (y != 1.0))) {
    Rcpp::stop("Propensity inputs must be finite and condition labels must be binary.");
  }
  if (!std::isfinite(ridge) || ridge <= 0.0)
    Rcpp::stop("propensity_ridge must be finite and positive.");

  mat X(N, p + 1, fill::ones);
  vec centers(p), scales(p);
  for (uword j = 0; j < p; j++) {
    const double m = mean(Xraw.col(j));
    double s = N > 1 ? stddev(Xraw.col(j)) : 0.0;
    if (!std::isfinite(m) || !std::isfinite(s))
      Rcpp::stop("Non-finite propensity scaling; check covariate scales.");
    if (s < 1e-12) s = 1.0;
    centers(j) = m;
    scales(j) = s;
    X.col(j + 1) = (Xraw.col(j) - m) / s;
  }
  if (!X.is_finite()) Rcpp::stop("Non-finite standardized propensity covariates.");
  if (center_out) *center_out = centers;
  if (scale_out) *scale_out = scales;
  if (N < 2 || p == 0 || all(y == 0) || all(y == 1)) {
    prob.fill(mean(y));
    return prob;
  }

  vec beta(p + 1, fill::zeros);
  const int max_iter = 25;
  const double tol = 1e-6;

  for (int it = 0; it < max_iter; it++) {
    vec eta = clamp(X * beta, -30.0, 30.0);
    vec pi_hat = 1.0 / (1.0 + exp(-eta));
    vec w = clamp(pi_hat % (1.0 - pi_hat), 1e-10, 0.25);
    vec z = eta + (y - pi_hat) / w;

    mat Xw = X;
    Xw.each_col() %= w;
    mat A = X.t() * Xw;
    A.diag() += ridge;
    vec b = X.t() * (w % z);

    vec beta_new;
    if (!solve(beta_new, A, b) || !beta_new.is_finite()) {
      beta_new = pinv(A) * b;
    }
    if (!beta_new.is_finite()) break;

    double delta = max(abs(beta_new - beta));
    beta = beta_new;
    if (delta < tol) break;
  }

  prob = clamp(1.0 / (1.0 + exp(-clamp(X * beta, -30.0, 30.0))), 1e-8, 1.0 - 1e-8);
  if (!prob.is_finite()) Rcpp::stop("Non-finite fitted propensity scores.");
  return prob;
}

// Exported wrapper retained for compatibility with v1.0.0.
// [[Rcpp::export]]
vec calculate_propensity_score(const mat& X, const vec& treatment_status) {
  return fit_logistic_propensity(X, treatment_status);
}

// Both conditions use every reference in the same order and from the same
// state: T1..TK followed by C1..CK. No cell's cluster assignment enters here.
// Each BLAS multiplication processes every cell for one centered reference.
// [[Rcpp::export]]
mat calculate_reference_quadratic_forms(
  const mat& theta_t, const mat& mu_treat, const mat& mu_ctrl,
  const cube& invcov_treat, const cube& invcov_ctrl) {
  const uword G = theta_t.n_cols;
  const uword K = mu_treat.n_cols;
  if (theta_t.n_rows == 0 || G == 0 || K == 0 ||
      mu_treat.n_rows != G || mu_ctrl.n_rows != G || mu_ctrl.n_cols != K ||
      invcov_treat.n_rows != G || invcov_treat.n_cols != G || invcov_treat.n_slices != K ||
      invcov_ctrl.n_rows != G || invcov_ctrl.n_cols != G || invcov_ctrl.n_slices != K ||
      !theta_t.is_finite() || !mu_treat.is_finite() || !mu_ctrl.is_finite() ||
      !invcov_treat.is_finite() || !invcov_ctrl.is_finite()) {
    Rcpp::stop("Reference quadratic forms require finite N x G expression, G x K means, and G x G x K precision arrays.");
  }
  mat features(theta_t.n_rows, 2 * K);
  for (uword k = 0; k < K; ++k) {
    mat centered = theta_t.each_row() - mu_treat.col(k).t();
    features.col(k) = sum((centered * invcov_treat.slice(k)) % centered, 1);
    centered = theta_t.each_row() - mu_ctrl.col(k).t();
    features.col(K + k) = sum((centered * invcov_ctrl.slice(k)) % centered, 1);
  }
  if (!features.is_finite())
    Rcpp::stop("Non-finite reference quadratic form; check expression, mean, and precision scales.");
  return features;
}

vec calculate_quadratic_form(
  const mat& theta, const cube& invcov, const vec& group, 
  const uvec& idx, const int& K, const int& G
) {
  vec qf(idx.n_elem, fill::zeros);
  int n_threads = std::min(omp_get_max_threads(), 16);
  
  #pragma omp parallel for schedule(dynamic, 32) num_threads(n_threads)
  for (uword i = 0; i < idx.n_elem; i++) {
    const uword k = static_cast<uword>(std::clamp(static_cast<int>(group(i)), 0, K-1));
    if (k >= static_cast<uword>(K) || k >= invcov.n_slices) continue;
    
    const uword cell_global_idx = idx(i);
    if (cell_global_idx >= theta.n_rows) continue;
    
    const rowvec theta_row = theta.row(cell_global_idx);
    qf(i) = dot(theta_row, theta_row * invcov.slice(k));
    
  }
  // Report numerical failure on the main thread instead of changing a cell's
  // network covariate to zero or throwing from an OpenMP worker.
  if (!qf.is_finite()) Rcpp::stop("Non-finite network quadratic form; check expression and precision scales.");
  return qf;
}

// The sampler needs only the map; retain full diagnostics for the exported
// helper and final state while avoiding repeated R allocations in the loop.
Rcpp::List network_propensity_matching_impl(
  const mat& theta_t, const mat& cell_t, 
  const vec& group_treat, const vec& group_ctrl,
  const uvec& treat_idx, const uvec& ctrl_idx,
  const cube& invcov_treat, const cube& invcov_ctrl,
  const int& G, const int& K,
  const mat* mu_treat, const mat* mu_ctrl,
  const std::string& network_covariates, const double propensity_ridge,
  const bool diagnostics, const uvec* fixed_map = nullptr,
  const std::string& matching_method = "greedy",
  SEXP matching_cost_function = R_NilValue) {

  if(matching_method!="lexicographic" && matching_method!="greedy")
    Rcpp::stop("matching_method must be lexicographic or greedy.");

  if (network_covariates != "legacy" && network_covariates != "all_reference")
    Rcpp::stop("network_covariates must be 'legacy' or 'all_reference'.");
  if (!std::isfinite(propensity_ridge) || propensity_ridge <= 0.0)
    Rcpp::stop("propensity_ridge must be finite and positive.");
  if (network_covariates == "all_reference" && (!mu_treat || !mu_ctrl))
    Rcpp::stop("all_reference requires both mu_treat and mu_ctrl.");
  const uword N = cell_t.n_rows;
  if (G < 1 || K < 1 || N == 0 || cell_t.n_cols == 0 || theta_t.n_rows != N ||
      theta_t.n_cols != static_cast<uword>(G) || !theta_t.is_finite() || !cell_t.is_finite() ||
      group_treat.n_elem != treat_idx.n_elem || group_ctrl.n_elem != ctrl_idx.n_elem ||
      treat_idx.is_empty() || ctrl_idx.is_empty() || treat_idx.max() >= N || ctrl_idx.max() >= N ||
      invcov_treat.n_rows != static_cast<uword>(G) || invcov_treat.n_cols != static_cast<uword>(G) ||
      invcov_ctrl.n_rows != static_cast<uword>(G) || invcov_ctrl.n_cols != static_cast<uword>(G) ||
      invcov_treat.n_slices != static_cast<uword>(K) || invcov_ctrl.n_slices != static_cast<uword>(K) ||
      !invcov_treat.is_finite() || !invcov_ctrl.is_finite() ||
      !group_treat.is_finite() || !group_ctrl.is_finite() ||
      any(group_treat < 0) || any(group_treat >= K) || any(group_treat != floor(group_treat)) ||
      any(group_ctrl < 0) || any(group_ctrl >= K) || any(group_ctrl != floor(group_ctrl))) {
    Rcpp::stop("Invalid dimensions, values, or zero-based cluster labels in NPM.");
  }
  vec treatment_status = cell_t.col(0);
  uvec treatment_indices = treat_idx;
  uvec control_indices = ctrl_idx;

  vec smd_results(K);
  smd_results.fill(datum::nan);
  // Preserve the historical assigned-cluster, uncentered vector as a returned
  // diagnostic. It is not a propensity predictor in all_reference mode.
  vec qf_global;
  if (diagnostics || network_covariates == "legacy") {
    const vec qf_treat = calculate_quadratic_form(theta_t, invcov_treat, group_treat, treat_idx, K, G);
    const vec qf_ctrl = calculate_quadratic_form(theta_t, invcov_ctrl, group_ctrl, ctrl_idx, K, G);
    qf_global.zeros(N);
    qf_global.elem(treat_idx) = qf_treat;
    qf_global.elem(ctrl_idx) = qf_ctrl;
  }
  mat network_features;
  if (network_covariates == "all_reference") {
    network_features = calculate_reference_quadratic_forms(
      theta_t, *mu_treat, *mu_ctrl, invcov_treat, invcov_ctrl);
  } else {
    network_features = qf_global;
  }
  mat Xcov;
  if (cell_t.n_cols > 1) {
    Xcov = join_horiz(cell_t.cols(1, cell_t.n_cols - 1), network_features);
  } else {
    Xcov = network_features;
  }
  vec centers, scales;
  vec base_propensity = fit_logistic_propensity(Xcov, treatment_status, propensity_ridge,
    diagnostics ? &centers : nullptr, diagnostics ? &scales : nullptr);

  if (base_propensity.n_elem != N) {
    Rcpp::stop("Propensity score length mismatch with cell count!");
  }
  
  vec prop_treat = base_propensity.elem(treatment_indices);
  vec prop_ctrl = base_propensity.elem(control_indices);

  vec cluster_mean_treat(K, fill::zeros), cluster_mean_ctrl(K, fill::zeros);
  vec difficulty_treat(K, fill::zeros);
  vec cluster_var_treat(K, fill::zeros), cluster_var_ctrl(K, fill::zeros);
  uvec count_treat(K, fill::zeros), count_ctrl(K, fill::zeros);
  const double overall_treatment_mean = mean(prop_treat);
  for (uword k = 0; k < static_cast<uword>(K); ++k) {
    const uvec it = find(group_treat == k), ic = find(group_ctrl == k);
    count_treat(k) = it.n_elem;
    count_ctrl(k) = ic.n_elem;
    if (!it.is_empty()) {
      const vec scores = prop_treat.elem(it);
      cluster_mean_treat(k) = mean(scores);
      cluster_var_treat(k) = sample_variance_or_zero(scores);
      difficulty_treat(k) = std::abs(cluster_mean_treat(k) - overall_treatment_mean);
    }
    if (!ic.is_empty()) {
      const vec scores = prop_ctrl.elem(ic);
      cluster_mean_ctrl(k) = mean(scores);
      cluster_var_ctrl(k) = sample_variance_or_zero(scores);
    }
  }

  mat smd_cost(K,K,fill::value(datum::nan));
  for(int k=0;k<K;++k) for(int c=0;c<K;++c)
    if(count_treat(k) && count_ctrl(c))
      smd_cost(k,c)=standardized_mean_difference(cluster_mean_treat(k),cluster_mean_ctrl(c),
        cluster_var_treat(k),cluster_var_ctrl(c));

  mat pair_scores;
  if (!Rf_isNull(matching_cost_function)) {
    Rcpp::Function callback(matching_cost_function);
    Rcpp::NumericMatrix callback_X = Rcpp::wrap(Xcov);
    callback_X.attr("matching_K")=K;
    Rcpp::List custom = callback(callback_X, treatment_status, group_treat,
                                group_ctrl, treat_idx, ctrl_idx);
    smd_cost = Rcpp::as<mat>(custom["cost"]);
    pair_scores = Rcpp::as<mat>(custom["scores"]);
    if (smd_cost.n_rows != (uword)K || smd_cost.n_cols != (uword)K ||
        pair_scores.n_rows != N || pair_scores.n_cols != (uword)K)
      Rcpp::stop("Invalid pairwise propensity dimensions.");
  }
  std::vector<uword> sorted_clusters;
  for (uword k = 0; k < static_cast<uword>(K); ++k)
    if (count_treat(k) > 0) sorted_clusters.push_back(k);
  std::stable_sort(sorted_clusters.begin(), sorted_clusters.end(),
    [&](uword a, uword b) { return difficulty_treat(a) > difficulty_treat(b); });
  // Stable sorting preserves the ascending index tie rule, including zero priorities.
  uvec ctrl_matched_flag(K, fill::zeros);
  uvec treat2ctrl_map(K);
  treat2ctrl_map.fill(K);
  if (fixed_map) {
    if (fixed_map->n_elem != static_cast<uword>(K)) Rcpp::stop("Invalid fixed map length.");
    treat2ctrl_map = *fixed_map;
    for (uword k = 0; k < static_cast<uword>(K); ++k) {
      const uword c = treat2ctrl_map(k);
      if (c >= static_cast<uword>(K) || ctrl_matched_flag(c))
        Rcpp::stop("A frozen map must be a complete bijection.");
      ctrl_matched_flag(c) = 1;
      if (count_treat(k) && count_ctrl(c))
        smd_results(k) = standardized_mean_difference(cluster_mean_treat(k), cluster_mean_ctrl(c),
          cluster_var_treat(k), cluster_var_ctrl(c));
    }
  } else if(matching_method=="lexicographic") {
    // Empty clusters cannot define an SMD. Keep a complete provisional map,
    // expose missing costs, and let the nonempty-cluster freeze guard reject it.
    mat assignment_cost=smd_cost;
    assignment_cost.replace(datum::nan,datum::inf);
    treat2ctrl_map=lexicographic_smd_assignment(assignment_cost);
    ctrl_matched_flag.ones();
    for(int k=0;k<K;++k) smd_results(k)=smd_cost(k,treat2ctrl_map(k));
  } else {
  for (const uword k : sorted_clusters) {
    double min_smd = datum::inf;
    uword best_c = K;
    for (uword c = 0; c < static_cast<uword>(K); ++c) {
      if (ctrl_matched_flag(c) || count_ctrl(c) == 0) continue;
      const double smd = standardized_mean_difference(cluster_mean_treat(k), cluster_mean_ctrl(c),
        cluster_var_treat(k), cluster_var_ctrl(c));
      if (smd < min_smd) {
        min_smd = smd;
        best_c = c;
      }
    }
    if (best_c < static_cast<uword>(K)) {
      treat2ctrl_map(k) = best_c;
      ctrl_matched_flag(best_c) = 1;
      smd_results(k) = min_smd;
    }
  }
  }

  if (!Rf_isNull(matching_cost_function)) {
    for (uword k=0; k<(uword)K; ++k) {
      uword c=treat2ctrl_map(k);
      if(c>=(uword)K) continue;
      smd_results(k)=smd_cost(k,c);
      for(uword i=0;i<treat_idx.n_elem;++i)
        if(group_treat(i)==k) base_propensity(treat_idx(i))=pair_scores(treat_idx(i),c);
      for(uword i=0;i<ctrl_idx.n_elem;++i)
        if(group_ctrl(i)==c) base_propensity(ctrl_idx(i))=pair_scores(ctrl_idx(i),k);
    }
  }
  if (!diagnostics)
    return Rcpp::List::create(Rcpp::_["treat2ctrl_cluster_map"] = treat2ctrl_map);

  int unmatched_treat_clusters = sum(treat2ctrl_map == K);
  int matched_pairs = K - unmatched_treat_clusters;

  vec treat_cluster(K);
  vec ctrl_cluster(K);
  
  for (uword k = 0; k < (uword)K; k++) {
    treat_cluster(k) = k + 1;
    if (treat2ctrl_map(k) < (uword)K) {
      ctrl_cluster(k) = treat2ctrl_map(k) + 1;  
    } else {
      ctrl_cluster(k) = NA_REAL;  
    }
  }
  
  Rcpp::DataFrame matching_result_df = Rcpp::DataFrame::create(
    Rcpp::Named("treatment_cluster") = treat_cluster,
    Rcpp::Named("matched_control_cluster") = ctrl_cluster
  );

  // SMD uses the same raw QC and network columns as the propensity fit.
  // Empty/unmatched clusters remain NaN; identical constants have SMD 0 and
  // unequal constants have infinite SMD, as in the legacy matcher.
  const uword n_cov = cell_t.n_cols - 1;
  const uword n_cov_tot = Xcov.n_cols;
  mat smd_cov(K, n_cov_tot);
  smd_cov.fill(datum::nan);
  {
    mat csumT(K, n_cov_tot, fill::zeros), csqT(K, n_cov_tot, fill::zeros);
    mat csumC(K, n_cov_tot, fill::zeros), csqC(K, n_cov_tot, fill::zeros);
    vec nT(K, fill::zeros), nC(K, fill::zeros);
    for (uword t = 0; t < group_treat.n_elem; t++) {
      const uword k = static_cast<uword>(group_treat(t));
      const uword row = treatment_indices(t);
      for (uword j = 0; j < n_cov_tot; j++) {
        const double v = Xcov(row, j);
        csumT(k, j) += v; csqT(k, j) += v * v;
      }
      nT(k)++;
    }
    for (uword t = 0; t < group_ctrl.n_elem; t++) {
      const uword k = static_cast<uword>(group_ctrl(t));
      const uword row = control_indices(t);
      for (uword j = 0; j < n_cov_tot; j++) {
        const double v = Xcov(row, j);
        csumC(k, j) += v; csqC(k, j) += v * v;
      }
      nC(k)++;
    }
    for (uword k = 0; k < static_cast<uword>(K); k++) {
      if (treat2ctrl_map(k) >= static_cast<uword>(K)) continue;
      const uword c = treat2ctrl_map(k);
      if (nT(k) < 1 || nC(c) < 1) continue;
      for (uword j = 0; j < n_cov_tot; j++) {
        const double m1 = csumT(k, j) / nT(k);
        const double m2 = csumC(c, j) / nC(c);
        double v1 = nT(k) > 1 ? (csqT(k, j) - nT(k) * m1 * m1) / (nT(k) - 1) : 0.0;
        double v2 = nC(c) > 1 ? (csqC(c, j) - nC(c) * m2 * m2) / (nC(c) - 1) : 0.0;
        if (v1 < 0) v1 = 0;
        if (v2 < 0) v2 = 0;
        smd_cov(k, j) = standardized_mean_difference(m1, m2, v1, v2);
      }
    }
  }
  Rcpp::CharacterVector network_names(network_features.n_cols), covariate_names(n_cov_tot);
  if (network_covariates == "all_reference") {
    for (int k = 0; k < K; ++k) {
      network_names[k] = "qf_ref_T" + std::to_string(k + 1);
      network_names[K + k] = "qf_ref_C" + std::to_string(k + 1);
    }
  } else {
    network_names[0] = "qf_network";
  }
  for (uword j = 0; j < n_cov; ++j)
    covariate_names[j] = "covariate_" + std::to_string(j + 1);
  for (uword j = 0; j < network_features.n_cols; ++j)
    covariate_names[n_cov + j] = network_names[j];
  Rcpp::NumericVector center_result(centers.begin(), centers.end());
  Rcpp::NumericVector scale_result(scales.begin(), scales.end());
  center_result.attr("names") = covariate_names;
  scale_result.attr("names") = covariate_names;
  Rcpp::NumericMatrix network_result = Rcpp::wrap(network_features);
  Rcpp::colnames(network_result) = network_names;

  return Rcpp::List::create(
    Rcpp::_["smd_results"] = smd_results,
    Rcpp::_["smd_cost_matrix"] = smd_cost,
    Rcpp::_["matching_method"] = matching_method,
    Rcpp::_["matching_is_frozen"] = fixed_map != nullptr,
    Rcpp::_["matching_max_smd"] = smd_results.has_nan() ? datum::nan : smd_results.max(),
    Rcpp::_["matching_total_smd"] = sum(smd_results),
    Rcpp::_["smd_covariate_results"] = smd_cov,
    Rcpp::_["propensity_scores"] = base_propensity,
    Rcpp::_["quadratic_form_global"] = qf_global,
    Rcpp::_["network_covariates"] = network_result,
    Rcpp::_["network_covariate_names"] = network_names,
    Rcpp::_["network_covariate_mode"] = network_covariates,
    Rcpp::_["propensity_scaling"] = Rcpp::List::create(
      Rcpp::_["center"] = center_result, Rcpp::_["scale"] = scale_result),
    Rcpp::_["treat2ctrl_cluster_map"] = treat2ctrl_map,
    Rcpp::_["matched_pairs"] = matched_pairs,
    Rcpp::_["ctrl_matched_flag"] = ctrl_matched_flag,
    Rcpp::_["matching_result_df"] = matching_result_df
  );
}

// The legacy default preserves calls to this internal exported helper from
// earlier versions. The public run_scNPM/MCMC default is all_reference.
// [[Rcpp::export]]
Rcpp::List network_propensity_matching_by_difficulty(
  const mat& theta_t, const mat& cell_t,
  const vec& group_treat, const vec& group_ctrl,
  const uvec& treat_idx, const uvec& ctrl_idx,
  const cube& invcov_treat, const cube& invcov_ctrl,
  const int& G, const int& K,
  Rcpp::Nullable<Rcpp::NumericMatrix> mu_treat = R_NilValue,
  Rcpp::Nullable<Rcpp::NumericMatrix> mu_ctrl = R_NilValue,
  std::string network_covariates = "legacy", double propensity_ridge = 1e-6,
  std::string matching_method = "greedy") {
  mat means_treat, means_ctrl;
  if (mu_treat.isNotNull()) means_treat = Rcpp::as<mat>(mu_treat.get());
  if (mu_ctrl.isNotNull()) means_ctrl = Rcpp::as<mat>(mu_ctrl.get());
  if ((mu_treat.isNotNull() && (means_treat.n_rows != static_cast<uword>(G) ||
       means_treat.n_cols != static_cast<uword>(K) || !means_treat.is_finite())) ||
      (mu_ctrl.isNotNull() && (means_ctrl.n_rows != static_cast<uword>(G) ||
       means_ctrl.n_cols != static_cast<uword>(K) || !means_ctrl.is_finite())))
    Rcpp::stop("mu_treat and mu_ctrl must be finite G x K matrices when supplied.");
  return network_propensity_matching_impl(theta_t, cell_t, group_treat, group_ctrl,
    treat_idx, ctrl_idx, invcov_treat, invcov_ctrl, G, K,
    mu_treat.isNotNull() ? &means_treat : nullptr,
    mu_ctrl.isNotNull() ? &means_ctrl : nullptr,
    network_covariates, propensity_ridge, true, nullptr, matching_method);
}

void update_edge_with_matching_by_difficulty(
  cube &edge_treat, cube &edge_ctrl,
  const cube &invcov_treat, const cube &invcov_ctrl,
  double ssp_v0, double ssp_v1, double ssp_xi, double edge_prob_threshold, double edge_prob_adj,
  const int &G, const int &K, const uword &N,
  const unsigned long &seed, 
  const mat &theta_t, const mat &cell_t, 
  const vec &group_treat, const vec &group_ctrl,
  const uvec &treat_idx, const uvec &ctrl_idx, const mat& edge_prior,
  const mat& mu_treat, const mat& mu_ctrl,
  const std::string& network_covariates, double propensity_ridge, const uvec* fixed_map = nullptr,
  const std::string& matching_method = "greedy",
  SEXP matching_cost_function = R_NilValue) {
  
  double edge_prob_threshold_up = std::max(edge_prob_threshold, 1.0 - edge_prob_threshold);
  double edge_prob_threshold_down = std::min(edge_prob_threshold, 1.0 - edge_prob_threshold);
  
  uvec treat2ctrl_map;
  if (fixed_map) treat2ctrl_map = *fixed_map;
  else {
    Rcpp::List matching_result = network_propensity_matching_impl(
      theta_t, cell_t, group_treat, group_ctrl, treat_idx, ctrl_idx,
      invcov_treat, invcov_ctrl, G, K,
      &mu_treat, &mu_ctrl, network_covariates, propensity_ridge, false, nullptr, matching_method, matching_cost_function);
    treat2ctrl_map = Rcpp::as<uvec>(matching_result["treat2ctrl_cluster_map"]);
  }
  uvec ctrl2treat_map(K);
  ctrl2treat_map.fill(K); 
  
  for (uword k = 0; k < (uword)K; k++) {
    if (treat2ctrl_map(k) < (uword)K) {
      ctrl2treat_map(treat2ctrl_map(k)) = k;
    }
  }

  const uword G_u = static_cast<uword>(G);
  const uword K_u = static_cast<uword>(K);
  // Each unordered pair is owned by one g iteration, so parallel workers
  // never update the same edge or read a counterpart owned by another worker.
  // Within each pair preserve the existing interleaved k, treatment/control scan.
  int n_threads = std::min(omp_get_max_threads(), 8); 
  #pragma omp parallel for schedule(dynamic) num_threads(n_threads)
  for (uword g = 0; g < G_u; g++) {
    std::mt19937 rng(seed + g * 10000);
    
    for (uword k = 0; k < K_u; k++) {
      
      if (k < edge_treat.n_slices && k < invcov_treat.n_slices) {
        const mat& invcov_k_treat = invcov_treat.slice(k);
        int matched_ctrl_cluster = treat2ctrl_map(k);
        bool is_matched = (matched_ctrl_cluster >= 0 && matched_ctrl_cluster < K);
        
        for (uword i = g+1; i < G_u; i++) {
          double omega_treat = 0.0;
          if (g < invcov_k_treat.n_rows && i < invcov_k_treat.n_cols) {
            omega_treat = invcov_k_treat(g, i);
          }
          
          const double xi = edge_prior(g, i);
          const double xi_present = std::max(xi * edge_prob_adj, std::numeric_limits<double>::min());
          const double xi_absent = std::min(xi / edge_prob_adj, std::max(xi, 1.0 - 1e-6));
          const double original_prob = calculate_edge_probability(omega_treat, ssp_v0, ssp_v1, xi);
          const double tmp = thread_safe_randu(rng);
          
          if (is_matched) {
            const uword c = static_cast<uword>(matched_ctrl_cluster);
            if (original_prob > edge_prob_threshold_up || original_prob < edge_prob_threshold_down) {
              edge_treat(g, i, k) = (tmp < original_prob) ? 1 : 0;
            } else if (c < edge_ctrl.n_slices) {
              double ctrl_edge_val = 0.0;
              if (g < edge_ctrl.n_rows && i < edge_ctrl.n_cols && c < edge_ctrl.n_slices) {
                ctrl_edge_val = edge_ctrl(g, i, c);
              }
              const double adjusted_xi = (ctrl_edge_val == 1) ? xi_present : xi_absent;
              const double adjusted_prob = calculate_edge_probability(omega_treat, ssp_v0, ssp_v1, adjusted_xi);
              edge_treat(g, i, k) = (tmp < adjusted_prob) ? 1 : 0;
            } else {
              edge_treat(g, i, k) = (tmp < original_prob) ? 1 : 0;
            }
          } else {
            edge_treat(g, i, k) = (tmp < original_prob) ? 1 : 0;
          }
          edge_treat(i, g, k) = edge_treat(g, i, k); // 对称矩阵
        }
      }

      if (k < edge_ctrl.n_slices && k < invcov_ctrl.n_slices) {
        const mat& invcov_k_ctrl = invcov_ctrl.slice(k);
        int matched_treat_cluster = ctrl2treat_map(k);
        bool is_matched_ctrl = (matched_treat_cluster >= 0 && matched_treat_cluster < K);
        
        for (uword i = g+1; i < G_u; i++) {
          double omega_ctrl = 0.0;
          if (g < invcov_k_ctrl.n_rows && i < invcov_k_ctrl.n_cols) {
            omega_ctrl = invcov_k_ctrl(g, i);
          }
          const double xi = edge_prior(g, i);
          const double xi_present = std::max(xi * edge_prob_adj, std::numeric_limits<double>::min());
          const double xi_absent = std::min(xi / edge_prob_adj, std::max(xi, 1.0 - 1e-6));
          const double original_prob = calculate_edge_probability(omega_ctrl, ssp_v0, ssp_v1, xi);
          const double tmp = thread_safe_randu(rng);

          if (is_matched_ctrl) {
            const uword t = static_cast<uword>(matched_treat_cluster);
            if (original_prob > edge_prob_threshold_up || original_prob < edge_prob_threshold_down) {
              edge_ctrl(g, i, k) = (tmp < original_prob) ? 1 : 0;
            } else if (t < edge_treat.n_slices) {
              double treat_edge_val = 0.0;
              if (g < edge_treat.n_rows && i < edge_treat.n_cols && t < edge_treat.n_slices) {
                treat_edge_val = edge_treat(g, i, t);
              }
              const double adjusted_xi = (treat_edge_val == 1) ? xi_present : xi_absent;
              const double adjusted_prob = calculate_edge_probability(omega_ctrl, ssp_v0, ssp_v1, adjusted_xi);
              edge_ctrl(g, i, k) = (tmp < adjusted_prob) ? 1 : 0;
            } else {
              edge_ctrl(g, i, k) = (tmp < original_prob) ? 1 : 0;
            }
          } else {
            edge_ctrl(g, i, k) = (tmp < original_prob) ? 1 : 0;
          }
          edge_ctrl(i, g, k) = edge_ctrl(g, i, k); // 对称矩阵
        }
      }
    }
  }
}

double U_lam(const vec &lambda, const vec &theta_miss, const vec &theta_obs, 
             const double &lam0_0, const double &lam1_0, const double &sigma2_lam0, const double &sigma2_lam1) {
  vec linear_obs = lambda(0) + lambda(1) * theta_obs;
  vec linear_miss = lambda(0) + lambda(1) * theta_miss;
  double U = - sum(rcpp_log_normsf(linear_obs)) - sum(rcpp_log_normcdf(linear_miss)) +
    (lambda(0) - lam0_0) * (lambda(0) - lam0_0) / (2.0 * sigma2_lam0) + 
    (lambda(1) - lam1_0) * (lambda(1) - lam1_0) / (2.0 * sigma2_lam1);
  return U;
}

void U_lam_grad(vec &U_grad, const vec &lambda, const vec &theta_miss, const vec &theta_obs, 
                const double &lam0_0, const double &lam1_0, const double &sigma2_lam0, const double &sigma2_lam1) {
  vec linear_obs = lambda(0) + lambda(1) * theta_obs;
  vec linear_miss = lambda(0) + lambda(1) * theta_miss;
  
  vec part1 = rcpp_normal_hazard(linear_obs, false);
  vec part2 = - rcpp_normal_hazard(linear_miss, true);
  
  U_grad(0) = sum(part1) + sum(part2) + (lambda(0) - lam0_0) / sigma2_lam0;
  U_grad(1) = sum(part1 % theta_obs) + sum(part2 % theta_miss) + (lambda(1) - lam1_0) / sigma2_lam1;
}

void HMC_lam(vec &lambda_new, double &log_r, const vec &lambda, const vec &theta_miss, const vec &theta_obs, const vec &p_rnorm,
             const double &lam0_0, const double &lam1_0, const double &sigma2_lam0, const double &sigma2_lam1,
             const double &epsilon, const int &num_step) {
  vec U_grad(2);
  
  lambda_new = lambda;
  vec p_new = p_rnorm;
  
  U_lam_grad(U_grad, lambda_new, theta_miss, theta_obs, lam0_0, lam1_0, sigma2_lam0, sigma2_lam1);
  p_new -= epsilon *  U_grad / 2.0;
  
  for (int i = 0; i < (num_step-1); i++) {
    lambda_new += epsilon * p_new;
    U_lam_grad(U_grad, lambda_new, theta_miss, theta_obs, lam0_0, lam1_0, sigma2_lam0, sigma2_lam1);
    p_new -= epsilon * U_grad;
  }
  lambda_new += epsilon * p_new;
  
  U_lam_grad(U_grad, lambda_new, theta_miss, theta_obs, lam0_0, lam1_0, sigma2_lam0, sigma2_lam1);
  p_new -= epsilon * U_grad / 2.0;
  
  p_new = -p_new;
  log_r = U_lam(lambda, theta_miss, theta_obs, lam0_0, lam1_0, sigma2_lam0, sigma2_lam1) - 
    U_lam(lambda_new, theta_miss, theta_obs, lam0_0, lam1_0, sigma2_lam0, sigma2_lam1) +
    sum(p_rnorm % p_rnorm) / 2.0 - sum(p_new % p_new) / 2.0;
}

void update_lambda(vec &lambda0_t, vec &lambda1_t, const mat &theta_t, const mat &ind_zero, const int &G, const unsigned long &seed,
                   const double &lam0_0=2, const double &lam1_0=-2, const double &sigma2_lam0=0.25, const double &sigma2_lam1=0.25,
                   const double &epsilon = 0.01, const int &num_step = 50) {
  int n_threads = std::min(omp_get_max_threads(), 24);
  
#pragma omp parallel shared(lambda0_t, lambda1_t) num_threads(n_threads)
{
  int rank = omp_get_thread_num();
  std::mt19937 rng(seed + rank * 10000);
  
#pragma omp for schedule(static)
  for (int g = 0; g < G; g++) {
    vec lambda_g(2);
    lambda_g(0) = lambda0_t(g);
    lambda_g(1) = lambda1_t(g);
    
    vec theta_g = theta_t.col(g);
    uvec ind_0= find(ind_zero.col(g) == true);
    uvec ind_obs = find(ind_zero.col(g) == false);
    
    vec theta_g_obs = theta_g(ind_obs);
    vec theta_g_miss = theta_g(ind_0);
    
    vec p_rnorm(2);
    p_rnorm(0) = thread_safe_rnorm(rng);
    p_rnorm(1) = thread_safe_rnorm(rng);
    
    double tmp_unif = thread_safe_randu(rng);
    
    vec lambda_g_new;
    double log_r;
    HMC_lam(lambda_g_new, log_r, lambda_g, theta_g_miss, theta_g_obs, p_rnorm, lam0_0, lam1_0, sigma2_lam0, sigma2_lam1, epsilon, num_step);
    
    if (tmp_unif < exp(log_r) && lambda_g_new(1) < 0) {
      lambda0_t(g) = lambda_g_new(0);
      lambda1_t(g) = lambda_g_new(1);
    }
  }
}
}

vec rDirichlet(const vec &alpha_vec) {
  if (alpha_vec.is_empty() || !alpha_vec.is_finite() || any(alpha_vec <= 0.0))
    Rcpp::stop("Dirichlet concentrations must be finite and strictly positive.");
  // This routine runs on the R main thread. Direct R gamma draws avoid
  // Armadillo RNG backend differences and preserve set.seed reproducibility.
  std::vector<long double> log_draws(alpha_vec.n_elem);
  long double largest = -std::numeric_limits<long double>::infinity();
  for (uword i = 0; i < alpha_vec.n_elem; ++i) {
    const double a = alpha_vec(i);
    const double gamma = R::rgamma(a < 1.0 ? a + 1.0 : a, 1.0);
    if (!std::isfinite(gamma) || gamma <= 0.0) Rcpp::stop("Non-finite Dirichlet gamma draw.");
    long double value = std::log(gamma);
    if (a < 1.0) {
      double u;
      do { u = R::runif(0.0, 1.0); } while (u <= 0.0);
      value += static_cast<long double>(std::log(u)) / a;
    }
    log_draws[i] = value;
    largest = std::max(largest, value);
  }
  vec result(alpha_vec.n_elem);
  for (uword i = 0; i < alpha_vec.n_elem; ++i)
    result(i) = std::max(static_cast<double>(std::exp(log_draws[i] - largest)),
                         std::numeric_limits<double>::min());
  return result / sum(result);
}

void update_pi(vec &pi_t, const vec &group_t, const vec &gam, const int &K) {
  vec tmp0(K, fill::zeros);
  for (int k = 0; k < K; k++) {
    uvec ind_j = find(group_t == k);
    tmp0(k) = ind_j.n_elem;
  }
  pi_t = rDirichlet(tmp0+gam);
}

void update_group(vec &group_t, const cube &invcov_t, const mat &theta_t, const mat &mu_t, const vec &pi_t, const int &G, const int &N, const int &K, const unsigned long &seed) {
  uvec pi_n0 = find(pi_t > 0);
  vec invcov_logdet(K);
  
  for (int k = 0; k < K; k++) {
    const mat chol_precision = chol(invcov_t.slice(static_cast<uword>(k)));
    invcov_logdet(k) = sum(log(chol_precision.diag()));
  }
  
  int n_threads = std::min(omp_get_max_threads(), 24);
  
#pragma omp parallel shared(group_t, pi_n0, invcov_logdet) num_threads(n_threads)
{
  int rank = omp_get_thread_num();
  std::mt19937 rng(seed + rank * 10000);
  
#pragma omp for schedule(static) 
  for (int i = 0; i < N; i++) {
    vec theta_i = theta_t.row(i).t();
    vec tmp(K);
    tmp.fill(- datum::inf);
    
    for (uword k_pi = 0; k_pi < pi_n0.n_elem; k_pi++) {
      int k = pi_n0(k_pi);
      vec tmp_k_pi = invcov_logdet(k) - 
        (theta_i - mu_t.col(k)).t() * invcov_t.slice(static_cast<uword>(k)) * (theta_i - mu_t.col(k))/2.0;
      tmp(k) = tmp_k_pi(0) + std::log(pi_t(k));
    }
    
    vec prob = exp(tmp - max(tmp));
    prob /= sum(prob);
    
    vec cum_prob = cumsum(prob);
    double u = thread_safe_randu(rng);
    int new_group = K - 1;
    for (int g = 0; g < K; g++) {
      if (u < cum_prob(g)) {
        new_group = g;
        break;
      }
    }
    group_t(i) = new_group;
  }
}
}

#include "staged_inference.h"

// [[Rcpp::export]]
Rcpp::List MCMC_with_difficulty_based_matching(
  const int num_iter, const int num_save, 
  mat theta_t, mat ind_zero,
  mat mu_treat, mat mu_ctrl,
  cube invcov_treat, cube cov_treat,
  cube invcov_ctrl, cube cov_ctrl,
  cube edge_treat, cube edge_ctrl,
  vec group_treat, vec group_ctrl,
  vec lambda0_t, vec lambda1_t,
  vec pi_treat, vec pi_ctrl,
  vec gam, mat cell_t,
  const int G, const int N_int, const int K,
  double ssp_v0, double ssp_v1, double ssp_l, double ssp_xi,
  double epsilon_theta = 0.2, int num_step_theta = 20,
  double eta_mu = 0.0, double tau_sq_mu = 1.0,
  double lam0_0 = 2.0, double lam1_0 = -2.0,
  double sigma2_lam0 = 0.25, double sigma2_lam1 = 0.25,
  double edge_prob_threshold = 0.6, double edge_prob_adj = 0.75,
  double epsilon_lam = 0.01, int num_step_lam = 10,
  Rcpp::Nullable<Rcpp::NumericMatrix> edge_prior = R_NilValue,
  std::string network_covariates = "all_reference", double propensity_ridge = 1e-6,
  std::string inference_mode = "coupled", int adapt_min = 1, int adapt_max = 1,
  int freeze_patience = 2, int final_burnin = 1, int diagnostic_every = 10,
  std::string matching_method = "greedy", bool verbose = false,
  std::string freeze_policy = "consensus", int freeze_window = 100,
  SEXP matching_cost_function = R_NilValue, double progress_every = 60.0,
  SEXP progress_callback = R_NilValue) {

  if (!std::isfinite(progress_every) || progress_every <= 0)
    Rcpp::stop("progress_every must be a positive number of seconds.");

  if(matching_method!="lexicographic" && matching_method!="greedy")
    Rcpp::stop("matching_method must be lexicographic or greedy.");

  if (network_covariates != "legacy" && network_covariates != "all_reference")
    Rcpp::stop("network_covariates must be 'legacy' or 'all_reference'.");
  if (!std::isfinite(propensity_ridge) || propensity_ridge <= 0.0)
    Rcpp::stop("propensity_ridge must be finite and positive.");

  // RcppArmadillo cube inputs can view the caller's R array storage even
  // when passed by value. Every mutable cube below owns an explicit copy;
  // assigning a temporary back into the original view would not detach it.
  cube invcov_treat_owned(invcov_treat.memptr(), invcov_treat.n_rows, invcov_treat.n_cols, invcov_treat.n_slices, true);
  cube cov_treat_owned(cov_treat.memptr(), cov_treat.n_rows, cov_treat.n_cols, cov_treat.n_slices, true);
  cube invcov_ctrl_owned(invcov_ctrl.memptr(), invcov_ctrl.n_rows, invcov_ctrl.n_cols, invcov_ctrl.n_slices, true);
  cube cov_ctrl_owned(cov_ctrl.memptr(), cov_ctrl.n_rows, cov_ctrl.n_cols, cov_ctrl.n_slices, true);
  cube edge_treat_owned(edge_treat.memptr(), edge_treat.n_rows, edge_treat.n_cols, edge_treat.n_slices, true);
  cube edge_ctrl_owned(edge_ctrl.memptr(), edge_ctrl.n_rows, edge_ctrl.n_cols, edge_ctrl.n_slices, true);

  if (G < 2 || K < 1 || N_int < 2 || num_iter < 1 || num_save < 1 || num_save > num_iter ||
      !std::isfinite(ssp_v0) || !std::isfinite(ssp_v1) || ssp_v0 <= 0 || ssp_v1 <= ssp_v0 ||
      !std::isfinite(ssp_l) || ssp_l <= 0 || !std::isfinite(ssp_xi) || ssp_xi <= 0 || ssp_xi >= 1 ||
      !std::isfinite(edge_prob_adj) || edge_prob_adj <= 0 || edge_prob_adj > 1 ||
      !std::isfinite(edge_prob_threshold) || edge_prob_threshold <= 0 || edge_prob_threshold >= 1 ||
      !std::isfinite(epsilon_theta) || epsilon_theta <= 0 || num_step_theta < 1 ||
      !std::isfinite(epsilon_lam) || epsilon_lam <= 0 || num_step_lam < 1 ||
      !std::isfinite(tau_sq_mu) || tau_sq_mu <= 0 ||
      !std::isfinite(sigma2_lam0) || sigma2_lam0 <= 0 ||
      !std::isfinite(sigma2_lam1) || sigma2_lam1 <= 0 ||
      gam.n_elem != static_cast<uword>(K) || !gam.is_finite() || any(gam <= 0)) {
    Rcpp::stop("Invalid MCMC dimensions, iterations, probabilities, or positive hyperparameters.");
  }
  mat edge_prior_mat(G, G, fill::value(ssp_xi));
  if (edge_prior.isNotNull()) {
    Rcpp::NumericMatrix supplied_prior(edge_prior);
    // Force owned storage before clearing the unused diagonal.
    edge_prior_mat = mat(supplied_prior.begin(), supplied_prior.nrow(), supplied_prior.ncol(), true);
    if (edge_prior_mat.n_rows != static_cast<uword>(G) || edge_prior_mat.n_cols != static_cast<uword>(G))
      Rcpp::stop("edge_prior must be a G x G matrix.");
    for (uword j = 0; j < static_cast<uword>(G); ++j) {
      for (uword t = j + 1; t < static_cast<uword>(G); ++t) {
        const double a = edge_prior_mat(j, t), b = edge_prior_mat(t, j);
        if (!std::isfinite(a) || !std::isfinite(b) || a <= 0 || a >= 1 || b <= 0 || b >= 1 ||
            std::abs(a - b) > 1e-12)
          Rcpp::stop("Off-diagonal edge_prior weights must be symmetric, finite, and strictly between 0 and 1.");
      }
    }
  }
  edge_prior_mat.diag().zeros();

  // NOTE: RcppArmadillo routes Armadillo's RNG (randi/randn) through R's RNG.
  // Reproducibility is therefore controlled at the R level: run_scNPM() calls
  // set.seed(ds_seed) immediately before invoking this function.
  
  if (theta_t.n_rows != (uword)N_int || theta_t.n_cols != (uword)G) {
    Rcpp::stop("theta_t dimension mismatch! Must be N x G");
  }
  if (ind_zero.n_rows != (uword)N_int || ind_zero.n_cols != (uword)G) {
    Rcpp::stop("ind_zero dimension mismatch! Must be N x G");
  }
  if (mu_treat.n_rows != (uword)G || mu_treat.n_cols != (uword)K) {
    Rcpp::stop("mu_treat dimension mismatch! Must be G x K");
  }
  if (invcov_treat_owned.n_rows != (uword)G || invcov_treat_owned.n_cols != (uword)G || invcov_treat_owned.n_slices != (uword)K) {
    Rcpp::stop("invcov_treat_owned dimension mismatch! Must be G x G x K");
  }
  if (cell_t.n_rows != (uword)N_int) {
    Rcpp::stop("cell_t row count mismatch! Must be N rows");
  }
  
  if (cell_t.n_cols < 1 || !cell_t.is_finite() || !theta_t.is_finite() || !ind_zero.is_finite() ||
      any((cell_t.col(0) != 0.0) % (cell_t.col(0) != 1.0)) ||
      any(vectorise((ind_zero != 0.0) % (ind_zero != 1.0))) ||
      mu_ctrl.n_rows != static_cast<uword>(G) || mu_ctrl.n_cols != static_cast<uword>(K) ||
      !mu_treat.is_finite() || !mu_ctrl.is_finite() ||
      lambda0_t.n_elem != static_cast<uword>(G) || lambda1_t.n_elem != static_cast<uword>(G) ||
      !lambda0_t.is_finite() || !lambda1_t.is_finite() || any(lambda1_t >= 0.0) ||
      pi_treat.n_elem != static_cast<uword>(K) || pi_ctrl.n_elem != static_cast<uword>(K) ||
      !pi_treat.is_finite() || !pi_ctrl.is_finite() || any(pi_treat <= 0.0) || any(pi_ctrl <= 0.0) ||
      !group_treat.is_finite() || !group_ctrl.is_finite() ||
      any(group_treat < 1) || any(group_treat > K) || any(group_treat != floor(group_treat)) ||
      any(group_ctrl < 1) || any(group_ctrl > K) || any(group_ctrl != floor(group_ctrl)))
    Rcpp::stop("Invalid initial states, binary condition indicators, or one-based MCMC cluster labels.");
  const cube* state_cubes[] = {&invcov_treat_owned, &invcov_ctrl_owned, &cov_treat_owned, &cov_ctrl_owned, &edge_treat_owned, &edge_ctrl_owned};
  for (const cube* state : state_cubes) {
    if (state->n_rows != static_cast<uword>(G) || state->n_cols != static_cast<uword>(G) ||
        state->n_slices != static_cast<uword>(K) || !state->is_finite())
      Rcpp::stop("Initial precision, covariance, and edge arrays must be finite G x G x K arrays.");
  }
  for (int k = 0; k < K; ++k) {
    if (!invcov_treat_owned.slice(k).is_sympd() || !invcov_ctrl_owned.slice(k).is_sympd() ||
        !cov_treat_owned.slice(k).is_sympd() || !cov_ctrl_owned.slice(k).is_sympd())
      Rcpp::stop("Initial precision and covariance matrices must be symmetric positive definite.");
    if (!approx_equal(edge_treat_owned.slice(k), edge_treat_owned.slice(k).t(), "absdiff", 1e-12) ||
        !approx_equal(edge_ctrl_owned.slice(k), edge_ctrl_owned.slice(k).t(), "absdiff", 1e-12) ||
        any(vectorise((edge_treat_owned.slice(k) != 0.0) % (edge_treat_owned.slice(k) != 1.0))) ||
        any(vectorise((edge_ctrl_owned.slice(k) != 0.0) % (edge_ctrl_owned.slice(k) != 1.0))))
      Rcpp::stop("Initial edge matrices must be symmetric and binary.");
  }
  pi_treat /= sum(pi_treat);
  pi_ctrl /= sum(pi_ctrl);
  group_treat = group_treat - 1;
  group_ctrl = group_ctrl - 1;
  
  if ((freeze_policy != "consensus" && freeze_policy != "error") || freeze_window < 1)
    Rcpp::stop("Invalid freeze policy or consensus window.");
  const bool staged = inference_mode == "staged";
  if (!staged && inference_mode != "coupled") Rcpp::stop("Unknown inference_mode.");
  int save_start = num_iter - num_save;
  if (diagnostic_every < 1 || (staged && (adapt_min < 1 || adapt_max < adapt_min ||
      freeze_patience < 2 || final_burnin < 1 ||
      adapt_max + final_burnin > save_start || any(abs(gam - gam(0)) > 1e-12))))
    Rcpp::stop("Invalid staged budget or non-exchangeable component prior.");
  const mat anchor_treat = mu_treat, anchor_ctrl = mu_ctrl;
  bool frozen = false;
  std::string freeze_reason = "not_frozen";
  bool matching_stable = false;
  double freeze_consensus_fraction = datum::nan;
  int freeze_window_used = 0;
  std::deque<FreezeRecord> freeze_records;
  int frozen_at = -1, stable_count = 0;
  uvec frozen_map(K), previous_map(K); previous_map.fill(K); frozen_map.fill(K);
  mat map_history(num_iter, K, fill::value(datum::nan));
  mat stage_history(num_iter, 5, fill::zeros);
  mat criterion_trace(num_iter, 6, fill::value(datum::nan));
  mat frozen_smd_cost(K,K,fill::value(datum::nan));
  const uword pair_count = static_cast<uword>(K) * K;
  cube joint00(G, G, staged ? pair_count : 0, fill::zeros);
  cube joint10 = joint00, joint01 = joint00, joint11 = joint00;
  cube empirical10 = joint00, empirical01 = joint00;
  const int batch_count = std::min(10, num_save);
  const int edge_count = G * (G - 1) / 2;
  cube batch10(staged ? edge_count : 0, pair_count, staged ? batch_count : 0, fill::zeros);
  cube batch01 = batch10;
  vec batch_sizes(staged ? batch_count : 0, fill::zeros);

  const uword N = static_cast<uword>(N_int);
  const uword G_u = static_cast<uword>(G);
  const uword K_u = static_cast<uword>(K);

  uvec treat_idx = find(cell_t.col(0) == 1);
  uvec ctrl_idx = find(cell_t.col(0) == 0);
  const uword n_treat = treat_idx.n_elem;
  const uword n_ctrl = ctrl_idx.n_elem;
  if (n_treat == 0 || n_ctrl == 0 || group_treat.n_elem != n_treat || group_ctrl.n_elem != n_ctrl)
    Rcpp::stop("Both conditions require cells and matching-length group label vectors.");

  mat theta_treat_t = theta_t.rows(treat_idx);
  mat theta_ctrl_t  = theta_t.rows(ctrl_idx);
  mat ind_zero_treat = ind_zero.rows(treat_idx);
  mat ind_zero_ctrl  = ind_zero.rows(ctrl_idx);

  // All R callbacks occur on the main thread, outside OpenMP regions.
  Rcpp::Function progress_message = Rcpp::Environment::base_env()["message"];
  if (verbose) progress_message("MCMC started...");
  const auto progress_start = std::chrono::steady_clock::now();
  auto last_progress = progress_start;
  auto report_progress = [&](int completed, const std::string& stage, bool force) {
    if (!verbose && Rf_isNull(progress_callback)) return;
    const auto now = std::chrono::steady_clock::now();
    if (!force && std::chrono::duration<double>(now-last_progress).count() < progress_every) return;
    const double elapsed = std::chrono::duration<double>(now-progress_start).count();
    const double eta = completed > 0 ? elapsed * (num_iter-completed)/completed : NA_REAL;
    Rcpp::List record = Rcpp::List::create(
      Rcpp::_["iteration"]=completed, Rcpp::_["total"]=num_iter,
      Rcpp::_["stage"]=stage, Rcpp::_["elapsed_seconds"]=elapsed,
      Rcpp::_["eta_seconds"]=eta, Rcpp::_["stable_count"]=stable_count,
      Rcpp::_["stable_required"]=freeze_patience,
      Rcpp::_["freeze_reason"]=freeze_reason,
      Rcpp::_["frozen_at"]=frozen ? frozen_at : NA_INTEGER,
      Rcpp::_["retained_draws"]=std::max(0,completed-save_start),
      Rcpp::_["retained_target"]=num_save);
    if (!Rf_isNull(progress_callback)) {
      Rcpp::Function callback(progress_callback); callback(record);
    } else if (verbose) {
      progress_message("MCMC iteration " + std::to_string(completed) + "/" +
        std::to_string(num_iter) + " | " + stage + " | elapsed " +
        std::to_string(elapsed) + " s | ETA " + std::to_string(eta) + " s");
    }
    last_progress=now;
  };
  report_progress(0,staged ? "adaptation" : "coupled",true);

  cube invcov_treat_save(G_u, G_u, K_u, fill::zeros);
  cube invcov_ctrl_save(G_u, G_u, K_u, fill::zeros);
  cube edge_treat_save(G_u, G_u, K_u, fill::zeros);
  cube edge_ctrl_save(G_u, G_u, K_u, fill::zeros);
  mat mu_treat_save(G_u, K_u, fill::zeros);
  mat mu_ctrl_save(G_u, K_u, fill::zeros);

  mat theta_treat_save(n_treat, G, fill::zeros);
  mat theta_ctrl_save(n_ctrl, G, fill::zeros);
  mat theta_save(N, G, fill::zeros);
  mat group_treat_save(K_u, n_treat, fill::zeros);
  mat group_ctrl_save(K_u, n_ctrl, fill::zeros);  
  
  mat group_iter(N, num_save, fill::zeros);
  cube mu_treat_iter(G_u, K_u, num_save, fill::zeros);
  cube mu_ctrl_iter(G_u, K_u, num_save, fill::zeros);
  arma::field<cube> invcov_treat_iter(num_save);
  arma::field<cube> invcov_ctrl_iter(num_save);
  
  for (int s = 0; s < num_save; s++) {
    invcov_treat_iter(s) = cube(G_u, G_u, K_u, fill::zeros);
    invcov_ctrl_iter(s) = cube(G_u, G_u, K_u, fill::zeros);
  }
  
  mat edge_treat_iter(G_u, num_save, fill::zeros);
  mat edge_ctrl_iter(G_u, num_save, fill::zeros);
  mat theta_treat_iter(n_treat, num_save, fill::zeros);
  mat theta_ctrl_iter(n_ctrl, num_save, fill::zeros);
  mat theta_iter(N, num_save, fill::zeros);
  mat lam0_iter(G_u, num_save, fill::zeros);
  mat lam1_iter(G_u, num_save, fill::zeros);
  mat pi_treat_iter(K_u, num_save, fill::zeros);
  mat pi_ctrl_iter(K_u, num_save, fill::zeros);

  vec group_t(N, fill::zeros);
  group_t(treat_idx) = group_treat;
  group_t(ctrl_idx) = group_ctrl + static_cast<double>(K);
  
  unsigned long iter_seed = static_cast<unsigned long>(randi(distr_param(1,1000000)));
  
  for (int t_iter = 0; t_iter < num_iter; t_iter++) {
    Rcpp::checkUserInterrupt();
    const bool feedback_this_iter = !staged || !frozen;
    const int stage_code = !staged ? 0 : (!frozen ? 1 : (t_iter < save_start ? 2 : 3));

    iter_seed = static_cast<unsigned long>(randi(distr_param(1, 1000000)));
    update_theta(theta_treat_t, ind_zero_treat, mu_treat, invcov_treat_owned, cov_treat_owned, 
                 lambda0_t, lambda1_t, group_treat, n_treat, iter_seed, epsilon_theta, num_step_theta);
    
    iter_seed = static_cast<unsigned long>(randi(distr_param(1, 1000000)));
    update_theta(theta_ctrl_t, ind_zero_ctrl, mu_ctrl, invcov_ctrl_owned, cov_ctrl_owned, 
                 lambda0_t, lambda1_t, group_ctrl, n_ctrl, iter_seed, epsilon_theta, num_step_theta);
    
    theta_t.rows(treat_idx) = theta_treat_t;
    theta_t.rows(ctrl_idx) = theta_ctrl_t;
    
    iter_seed = static_cast<unsigned long>(randi(distr_param(1,1000000)));
    update_lambda(lambda0_t, lambda1_t, theta_t, ind_zero, G, iter_seed, lam0_0, lam1_0, sigma2_lam0, sigma2_lam1, epsilon_lam, num_step_lam);
    
    update_mu(mu_treat, theta_treat_t, invcov_treat_owned, group_treat, G, K, eta_mu, tau_sq_mu);
    update_mu(mu_ctrl, theta_ctrl_t, invcov_ctrl_owned, group_ctrl, G, K, eta_mu, tau_sq_mu);
    
    iter_seed = static_cast<unsigned long>(randi(distr_param(1, 1000000)));
    update_invcov(invcov_treat_owned, cov_treat_owned, edge_treat_owned, theta_treat_t, mu_treat, group_treat, 
                  ssp_v0, ssp_v1, ssp_l, G, K, iter_seed);
    
    iter_seed = static_cast<unsigned long>(randi(distr_param(1, 1000000)));
    update_invcov(invcov_ctrl_owned, cov_ctrl_owned, edge_ctrl_owned, theta_ctrl_t, mu_ctrl, group_ctrl, 
                  ssp_v0, ssp_v1, ssp_l, G, K, iter_seed);
    
    iter_seed = static_cast<unsigned long>(randi(distr_param(1,1000000))); 
    update_edge_with_matching_by_difficulty(
      edge_treat_owned, edge_ctrl_owned, invcov_treat_owned, invcov_ctrl_owned, ssp_v0, ssp_v1, ssp_xi, edge_prob_threshold, feedback_this_iter ? edge_prob_adj : 1.0,
      G, K, N, iter_seed,
      theta_t, cell_t, group_treat, group_ctrl, treat_idx, ctrl_idx, edge_prior_mat,
      mu_treat, mu_ctrl, network_covariates, propensity_ridge, frozen ? &frozen_map : nullptr, matching_method, matching_cost_function
    );
    
    iter_seed = static_cast<unsigned long>(randi(distr_param(1, 1000000)));
    update_group(group_treat, invcov_treat_owned, theta_treat_t, mu_treat, pi_treat, G, n_treat, K, iter_seed);
    
    iter_seed = static_cast<unsigned long>(randi(distr_param(1, 1000000)));
    update_group(group_ctrl, invcov_ctrl_owned, theta_ctrl_t, mu_ctrl, pi_ctrl, G, n_ctrl, K, iter_seed);
    
    group_t(treat_idx) = group_treat;
    group_t(ctrl_idx) = group_ctrl + static_cast<double>(K);
    
    update_pi(pi_treat, group_treat, gam, K);
    update_pi(pi_ctrl, group_ctrl, gam, K);
    
    if (staged) {
      align_components(anchor_treat, mu_treat, invcov_treat_owned, cov_treat_owned,
                       edge_treat_owned, pi_treat, group_treat);
      align_components(anchor_ctrl, mu_ctrl, invcov_ctrl_owned, cov_ctrl_owned,
                       edge_ctrl_owned, pi_ctrl, group_ctrl);
      group_t(treat_idx) = group_treat;
      group_t(ctrl_idx) = group_ctrl + static_cast<double>(K);
    }
    const bool diagnose = staged && (!frozen || (t_iter + 1) % diagnostic_every == 0 ||
                                     t_iter == num_iter - 1 || t_iter == save_start);
    if (diagnose) {
      Rcpp::List current_match = network_propensity_matching_impl(
        theta_t, cell_t, group_treat, group_ctrl, treat_idx, ctrl_idx,
        invcov_treat_owned, invcov_ctrl_owned, G, K, &mu_treat, &mu_ctrl,
        network_covariates, propensity_ridge, true, frozen ? &frozen_map : nullptr, matching_method, matching_cost_function);
      uvec current_map = Rcpp::as<uvec>(current_match["treat2ctrl_cluster_map"]);
      if (!frozen) {
        bool complete = all(current_map < static_cast<uword>(K));
        for (int k = 0; k < K; ++k)
          complete = complete && any(group_treat == k) && any(group_ctrl == k);
        mat current_cost = Rcpp::as<mat>(current_match["smd_cost_matrix"]);
        bool valid_cost = current_cost.is_finite();
        if(complete && valid_cost) {
          freeze_records.push_back(FreezeRecord{current_map,current_cost});
          if((int)freeze_records.size()>freeze_window) freeze_records.pop_front();
        }
        complete = complete && valid_cost;
        stable_count = complete ? (all(current_map == previous_map) ? stable_count + 1 : 1) : 0;
        previous_map = current_map;
        if (t_iter + 1 >= adapt_min && stable_count >= freeze_patience) {
          frozen = true; frozen_at = t_iter + 1; frozen_map = current_map;
          frozen_smd_cost=current_cost;
          freeze_reason="stable"; matching_stable=true;
          freeze_window_used=(int)freeze_records.size();
          int votes=0; for(const auto& r:freeze_records) if(all(r.map==frozen_map)) ++votes;
          freeze_consensus_fraction=double(votes)/freeze_window_used;
        }
        if (!frozen && t_iter + 1 >= adapt_max) {
          if(freeze_policy=="error")
            Rcpp::stop("Matching did not stabilize by adapt_max (strict freeze_policy=error).");
          if(freeze_records.empty())
            Rcpp::stop("No valid complete finite-cost matching was observed; inspect empty clusters or numerical inputs.");
          FreezeChoice choice=select_freeze_consensus(freeze_records);
          frozen=true; frozen_at=t_iter+1; frozen_map=choice.map;
          frozen_smd_cost=choice.average_cost;
          freeze_reason="budget_consensus"; matching_stable=false;
          freeze_consensus_fraction=choice.fraction;
          freeze_window_used=(int)freeze_records.size();
          current_map=frozen_map;
          if(verbose) Rcpp::message(Rcpp::wrap("Adaptation budget reached: froze consensus map; stability not established. Final burn-in and sampling continue."));
        }
      }
      map_history.row(t_iter) = conv_to<rowvec>::from(current_map + 1);
      vec smd = Rcpp::as<vec>(current_match["smd_results"]);
      if(freeze_reason=="budget_consensus" && frozen_at==t_iter+1) {
        mat cost_now=Rcpp::as<mat>(current_match["smd_cost_matrix"]);
        for(int k=0;k<K;++k) smd(k)=cost_now(k,current_map(k));
      }
      double match_loss = smd.is_finite() ? sum(smd) : datum::nan;
      const double ll = augmented_likelihood(theta_t, ind_zero, mu_treat, mu_ctrl,
        invcov_treat_owned, invcov_ctrl_owned, group_treat, group_ctrl,
        treat_idx, ctrl_idx, pi_treat, pi_ctrl, lambda0_t, lambda1_t);
      const double lp = network_prior_kernel(invcov_treat_owned, invcov_ctrl_owned,
        edge_treat_owned, edge_ctrl_owned, edge_prior_mat, current_map,
        ssp_v0, ssp_v1, ssp_l, edge_prob_threshold,
        feedback_this_iter ? edge_prob_adj : 1.0);
      criterion_trace.row(t_iter) = rowvec({double(t_iter+1), ll, lp, match_loss,
                                           ll+lp-match_loss, double(stage_code)});
    }
    if (staged && frozen && !diagnose)
      map_history.row(t_iter) = conv_to<rowvec>::from(frozen_map + 1);
    stage_history.row(t_iter) = rowvec({double(t_iter+1), double(stage_code),
      feedback_this_iter && edge_prob_adj < 1.0 ? 1.0 : 0.0, double(frozen), double(stable_count)});
    if (t_iter >= save_start) {
      int save_i = t_iter - save_start;
      if (staged) {
        if (!frozen || save_start - frozen_at < final_burnin)
          Rcpp::stop("Final samples must follow freezing and post-feedback burn-in.");
        accumulate_joint_probabilities(invcov_treat_owned, invcov_ctrl_owned,
          edge_treat_owned, edge_ctrl_owned, edge_prior_mat, ssp_v0, ssp_v1,
          joint00, joint10, joint01, joint11, empirical10, empirical01,
          batch10, batch01, batch_sizes, save_i * batch_count / num_save);
      }

      
      for (uword i = 0; i < n_treat; i++) {
        int k_t = static_cast<int>(group_treat(i));
        if (k_t >= 0 && k_t < K) {
          group_treat_save(k_t, i) += 1;
        }
      }
      
      for (uword i = 0; i < n_ctrl; i++) {
        int k_c = static_cast<int>(group_ctrl(i));
        if (k_c >= 0 && k_c < K) {
          group_ctrl_save(k_c, i) += 1;
        }
      }
      
      mu_treat_save   += mu_treat;
      mu_ctrl_save    += mu_ctrl;
      invcov_treat_save += invcov_treat_owned;
      invcov_ctrl_save += invcov_ctrl_owned;
      edge_treat_save  += edge_treat_owned;
      edge_ctrl_save   += edge_ctrl_owned;
      
      theta_treat_save += theta_treat_t;
      theta_ctrl_save  += theta_ctrl_t;
      theta_save       += theta_t;

      group_iter.col(save_i)      = group_t;
      mu_treat_iter.slice(save_i) = mu_treat;
      mu_ctrl_iter.slice(save_i)  = mu_ctrl;
      invcov_treat_iter(save_i)   = invcov_treat_owned;
      invcov_ctrl_iter(save_i)    = invcov_ctrl_owned;
      
      if (G_u > 0) {
        edge_treat_iter.col(save_i) = edge_treat_owned.slice(0).col(0);
        edge_ctrl_iter.col(save_i)  = edge_ctrl_owned.slice(0).col(0);
      }
      
      lam0_iter.col(save_i)       = lambda0_t;
      lam1_iter.col(save_i)       = lambda1_t;
      pi_treat_iter.col(save_i)   = pi_treat;
      pi_ctrl_iter.col(save_i)    = pi_ctrl;
      
      if (n_treat > 0) theta_treat_iter.col(save_i) = theta_treat_t.col(0);
      if (n_ctrl > 0) theta_ctrl_iter.col(save_i)  = theta_ctrl_t.col(0);
      if (N > 0) theta_iter.col(save_i)            = theta_t.col(0);
    }
    const bool just_frozen = frozen && frozen_at==t_iter+1;
    const std::string progress_stage = t_iter+1==num_iter ? "completed" :
      (stage_code==0 ? "coupled" : (stage_code==1 ? "adaptation" :
        (stage_code==2 ? "final_burnin" : "retained")));
    report_progress(t_iter+1,progress_stage,
      just_frozen || t_iter==save_start || t_iter+1==num_iter);
  }
  
  if (verbose) progress_message("MCMC completed!");
  
  mu_treat_save   /= num_save;
  mu_ctrl_save    /= num_save;
  invcov_treat_save /= num_save;
  invcov_ctrl_save /= num_save;
  edge_treat_save  /= num_save;
  edge_ctrl_save   /= num_save;
  theta_treat_save /= num_save;
  theta_ctrl_save  /= num_save;
  theta_save       /= num_save;
  
  vec group_treat_est(n_treat);
  vec group_ctrl_est(n_ctrl);
  
  for (uword i = 0; i < n_treat; i++) {
    group_treat_est(i) = group_treat_save.col(i).index_max() + 1;
  }
  
  for (uword i = 0; i < n_ctrl; i++) {
    group_ctrl_est(i) = group_ctrl_save.col(i).index_max() + 1;
  }
  
  vec group_est(N, fill::zeros);
  group_est(treat_idx) = group_treat_est;
  group_est(ctrl_idx)  = group_ctrl_est + static_cast<double>(K);
  
  Rcpp::List final_matching = network_propensity_matching_impl(
    theta_t, cell_t, group_treat, group_ctrl, treat_idx, ctrl_idx,
    invcov_treat_owned, invcov_ctrl_owned, G, K,
    &mu_treat, &mu_ctrl, network_covariates, propensity_ridge, true, frozen ? &frozen_map : nullptr, matching_method, matching_cost_function);
  if (staged) {
    joint00 /= num_save; joint10 /= num_save; joint01 /= num_save; joint11 /= num_save;
    empirical10 /= num_save; empirical01 /= num_save;
    for (int b = 0; b < batch_count; ++b) { batch10.slice(b) /= batch_sizes(b); batch01.slice(b) /= batch_sizes(b); }
  }
  return Rcpp::List::create(
    Rcpp::Named("group_est")        = group_est, 
    Rcpp::Named("group_treat_est")  = group_treat_est, 
    Rcpp::Named("group_ctrl_est")   = group_ctrl_est,
    Rcpp::Named("mu_treat_est")     = mu_treat_save, 
    Rcpp::Named("mu_ctrl_est")      = mu_ctrl_save,
    Rcpp::Named("invcov_treat_est") = invcov_treat_save, 
    Rcpp::Named("invcov_ctrl_est")  = invcov_ctrl_save,
    Rcpp::Named("edge_treat_est")   = edge_treat_save, 
    Rcpp::Named("edge_ctrl_est")    = edge_ctrl_save,
    Rcpp::Named("theta_treat_est")  = theta_treat_save, 
    Rcpp::Named("theta_ctrl_est")   = theta_ctrl_save, 
    Rcpp::Named("theta_est")        = theta_save,
    Rcpp::Named("group_iter")       = group_iter,
    Rcpp::Named("mu_treat_iter")    = mu_treat_iter, 
    Rcpp::Named("mu_ctrl_iter")     = mu_ctrl_iter,
    Rcpp::Named("invcov_treat_iter")= invcov_treat_iter, 
    Rcpp::Named("invcov_ctrl_iter") = invcov_ctrl_iter,
    Rcpp::Named("edge_treat_iter")  = edge_treat_iter, 
    Rcpp::Named("edge_ctrl_iter")  = edge_ctrl_iter,
    Rcpp::Named("theta_treat_iter") = theta_treat_iter, 
    Rcpp::Named("theta_ctrl_iter")  = theta_ctrl_iter, 
    Rcpp::Named("theta_iter")       = theta_iter,
    Rcpp::Named("lam0_iter")        = lam0_iter, 
    Rcpp::Named("lam1_iter")        = lam1_iter, 
    Rcpp::Named("pi_treat_iter")    = pi_treat_iter, 
    Rcpp::Named("pi_ctrl_iter")     = pi_ctrl_iter,
    Rcpp::Named("group_treat_save") = group_treat_save, 
    Rcpp::Named("group_ctrl_save")  = group_ctrl_save,
    Rcpp::Named("matching_results") = final_matching,
    Rcpp::Named("stage_history") = stage_history,
    Rcpp::Named("matching_history") = map_history,
    Rcpp::Named("criterion_trace") = criterion_trace,
    Rcpp::Named("frozen_at") = frozen_at,
    Rcpp::Named("freeze_reason") = freeze_reason,
    Rcpp::Named("matching_stable") = matching_stable,
    Rcpp::Named("freeze_consensus_fraction") = freeze_consensus_fraction,
    Rcpp::Named("freeze_window_used") = freeze_window_used,
    Rcpp::Named("frozen_smd_cost_matrix") = frozen_smd_cost,
    Rcpp::Named("label_anchors") = Rcpp::List::create(Rcpp::_["treatment"] = anchor_treat, Rcpp::_["control"] = anchor_ctrl),
    Rcpp::Named("joint_edge_probabilities") = Rcpp::List::create(
      Rcpp::_["p00"] = joint00, Rcpp::_["p10"] = joint10,
      Rcpp::_["p01"] = joint01, Rcpp::_["p11"] = joint11,
      Rcpp::_["empirical_p10"] = empirical10, Rcpp::_["empirical_p01"] = empirical01,
      Rcpp::_["batch_p10"] = batch10, Rcpp::_["batch_p01"] = batch01,
      Rcpp::_["batch_sizes"] = batch_sizes)
  );
}

// [[Rcpp::export]]
mat update_mu_R(mat theta_t, cube invcov_t, vec group_t, int G, int K, double eta_mu = 0, double tau_sq_mu = 1) {
  mat mu_updated(G, K);
  mat I_tau_sq = eye(G,G);
  I_tau_sq /= tau_sq_mu;
  
  for (int k = 0; k < K; k++) {
    uvec ind_k = find(group_t == k);  
    vec tmp1(G, fill::zeros);
    if (!ind_k.is_empty()) tmp1 = sum(theta_t.rows(ind_k), 0).t();
    mat invcov_k = invcov_t.slice(static_cast<uword>(k));
    
    mat COV = inv(ind_k.n_elem * invcov_k + I_tau_sq);
    vec MU = COV * (invcov_k * tmp1 + eta_mu/tau_sq_mu);
    
    vec mu_k = mvnrnd(MU, COV, 1);
    mu_updated.col(k) = mu_k;
  }
  
  return mu_updated;
}

// [[Rcpp::export]]
vec update_pi_R(vec group_t, vec gam, int K) {
  if (K < 1 || gam.n_elem != static_cast<uword>(K) || !gam.is_finite() || any(gam <= 0.0) ||
      !group_t.is_finite() || any(group_t < 0) || any(group_t >= K) || any(group_t != floor(group_t)))
    Rcpp::stop("update_pi_R requires positive concentrations and zero-based labels in 0, ..., K-1.");
  vec pi_new(K);
  vec tmp0(K, fill::zeros);
  
  for (int k = 0; k < K; k++) {
    uvec ind_j = find(group_t == k);
    tmp0(k) = ind_j.n_elem;
  }
  
  pi_new = rDirichlet(tmp0 + gam);
  return pi_new;
}

// [[Rcpp::export]]
arma::uvec staged_assignment(const arma::mat& cost) {
  return minimum_assignment(cost) + 1;
}

// [[Rcpp::export]]
arma::uvec lexicographic_assignment(const arma::mat& cost) {
  return lexicographic_smd_assignment(cost)+1;
}
