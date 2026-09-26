#ifndef SCNPM_GLOBAL_MATCHING_H
#define SCNPM_GLOBAL_MATCHING_H
#include <functional>
// Minimum-cost assignment; ascending traversal fixes all exact ties.
inline uvec minimum_assignment(const mat& cost) {
  if (!cost.is_square() || cost.n_rows == 0 || !cost.is_finite())
    Rcpp::stop("Assignment costs must be a nonempty finite square matrix.");
  const int n = cost.n_rows;
  std::vector<double> u(n+1), v(n+1);
  std::vector<int> p(n+1), way(n+1);
  for (int i=1; i<=n; ++i) {
    p[0]=i; int j0=0;
    std::vector<double> minv(n+1, datum::inf);
    std::vector<bool> used(n+1, false);
    do {
      used[j0]=true; int i0=p[j0], j1=0; double delta=datum::inf;
      for (int j=1; j<=n; ++j) if (!used[j]) {
        double cur=cost(i0-1,j-1)-u[i0]-v[j];
        if (cur<minv[j]) {minv[j]=cur; way[j]=j0;}
        if (minv[j]<delta) {delta=minv[j]; j1=j;}
      }
      for (int j=0; j<=n; ++j) {
        if (used[j]) {u[p[j]]+=delta; v[j]-=delta;}
        else minv[j]-=delta;
      }
      j0=j1;
    } while (p[j0]!=0);
    do {int j1=way[j0]; p[j0]=p[j1]; j0=j1;} while (j0);
  }
  uvec out(n);
  for (int j=1; j<=n; ++j) out[p[j]-1]=j-1;
  return out;
}


// Perfect matching feasibility on edges at or below a candidate bottleneck.
inline bool bottleneck_feasible(const mat& cost, double limit) {
  const uword n=cost.n_rows;
  std::vector<int> owner(n,-1);
  std::function<bool(uword,std::vector<bool>&)> augment;
  augment = [&](uword r, std::vector<bool>& seen) {
    for(uword c=0;c<n;++c) if(!seen[c] && cost(r,c)<=limit) {
      seen[c]=true;
      if(owner[c]<0 || augment(static_cast<uword>(owner[c]),seen)) {
        owner[c]=static_cast<int>(r); return true;
      }
    }
    return false;
  };
  for(uword r=0;r<n;++r) {
    std::vector<bool> seen(n,false);
    if(!augment(r,seen)) return false;
  }
  return true;
}

// Minimize (maximum SMD, total SMD), then resolve ties deterministically.
// For K<=8 enumerate permutations in ascending lexicographic order (K=4: 24).
// For larger K find the bottleneck by feasibility, then a minimum-sum assignment
// restricted to that bottleneck. Ascending Hungarian traversal fixes exact ties.
inline uvec lexicographic_smd_assignment(const mat& cost) {
  if(!cost.is_square() || cost.n_rows==0 || cost.has_nan() || any(vectorise(cost)<0))
    Rcpp::stop("SMD costs must be a nonempty nonnegative square matrix without NA/NaN.");
  const uword n=cost.n_rows;
  uvec best=regspace<uvec>(0,n-1);
  if(n<=8) {
    std::vector<uword> p(n); std::iota(p.begin(),p.end(),0);
    double best_b=datum::inf; long double best_s=std::numeric_limits<long double>::infinity();
    do {
      double b=0; long double total=0;
      for(uword k=0;k<n;++k) { b=std::max(b,cost(k,p[k])); total+=cost(k,p[k]); }
      if(b<best_b || (b==best_b && total<best_s)) {
        best_b=b; best_s=total;
        for(uword k=0;k<n;++k) best(k)=p[k];
      }
    } while(std::next_permutation(p.begin(),p.end()));
    return best;
  }
  std::vector<double> levels;
  for(uword i=0;i<cost.n_elem;++i) if(std::isfinite(cost(i))) levels.push_back(cost(i));
  std::sort(levels.begin(),levels.end());
  levels.erase(std::unique(levels.begin(),levels.end()),levels.end());
  // If every complete permutation contains infinity, B=S=infinity for all:
  // the deterministic tertiary rule is the identity permutation.
  if(levels.empty() || !bottleneck_feasible(cost,levels.back())) return best;
  size_t lo=0,hi=levels.size()-1;
  while(lo<hi) {const size_t mid=lo+(hi-lo)/2;
    if(bottleneck_feasible(cost,levels[mid])) hi=mid; else lo=mid+1;}
  const double bottleneck=levels[lo];
  mat restricted(n,n,fill::value(static_cast<double>(n)+1.0));
  for(uword k=0;k<n;++k) for(uword c=0;c<n;++c)
    if(cost(k,c)<=bottleneck) restricted(k,c)=bottleneck>0?cost(k,c)/bottleneck:0;
  best=minimum_assignment(restricted);
  for(uword k=0;k<n;++k) if(cost(k,best(k))>bottleneck)
    Rcpp::stop("Internal bottleneck assignment failure.");
  return best;
}
#endif
