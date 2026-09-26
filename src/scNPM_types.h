// scNPM_types.h (最终版：无 rTRNG 依赖)
#ifndef SCNPM_TYPES_H
#define SCNPM_TYPES_H

// Rcpp 属性标签
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

// 宏定义
#define ARMA_NO_DEBUG
#define RCPP_USE_GLOBAL_STRINGSTREAM

// 标准库
#ifdef _OPENMP
#include <omp.h>
#else
inline int omp_get_max_threads() { return 1; }
inline int omp_get_thread_num() { return 0; }
#endif
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>
#include <vector>

using namespace arma;

// 支持 Armadillo 向量的 normpdf/normcdf
inline vec rcpp_normpdf(const vec& x) {
  vec res(x.n_elem);
  for (uword i = 0; i < x.n_elem; i++) {
    res(i) = R::dnorm(x(i), 0.0, 1.0, 0);
  }
  return res;
}

inline vec rcpp_normcdf(const vec& x) {
  vec res(x.n_elem);
  for (uword i = 0; i < x.n_elem; i++) {
    res(i) = R::pnorm(x(i), 0.0, 1.0, 1, 0);
  }
  return res;
}

// 保留单个 double 版本
inline double rcpp_normpdf(double x) {
  return R::dnorm(x, 0.0, 1.0, 0);
}

inline double rcpp_normcdf(double x) {
  return R::pnorm(x, 0.0, 1.0, 1, 0);
}

// Pure Rmath density/tail evaluations use no R random-number generator.
inline vec rcpp_log_normcdf(const vec& x) {
  vec out(x.n_elem);
  for (uword i = 0; i < x.n_elem; ++i) out(i) = R::pnorm(x(i), 0.0, 1.0, 1, 1);
  return out;
}

inline vec rcpp_log_normsf(const vec& x) {
  vec out(x.n_elem);
  for (uword i = 0; i < x.n_elem; ++i) out(i) = R::pnorm(x(i), 0.0, 1.0, 0, 1);
  return out;
}

inline vec rcpp_normal_hazard(const vec& x, bool lower_tail) {
  vec out(x.n_elem);
  for (uword i = 0; i < x.n_elem; ++i)
    out(i) = std::exp(R::dnorm(x(i), 0.0, 1.0, 1) -
                     R::pnorm(x(i), 0.0, 1.0, lower_tail ? 1 : 0, 1));
  return out;
}

#endif // SCNPM_TYPES_H
