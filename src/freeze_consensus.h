#pragma once
#include <deque>
#include <vector>
struct FreezeRecord { arma::uvec map; arma::mat cost; };
struct FreezeChoice { arma::uvec map; arma::mat average_cost; double fraction; };
inline FreezeChoice select_freeze_consensus(const std::deque<FreezeRecord>& records) {
  if(records.empty()) Rcpp::stop("No consensus records");
  const arma::uword K=records.front().map.n_elem;
  arma::mat avg(K,K,arma::fill::zeros);
  for(const auto& r:records)avg+=r.cost;
  avg/=records.size();
  int best_count=-1; double best_max=arma::datum::inf,best_sum=arma::datum::inf;
  arma::uvec best;
  for(const auto& candidate:records) {
    int count=0;for(const auto& r:records)if(arma::all(r.map==candidate.map))++count;
    double mx=0,total=0;for(arma::uword k=0;k<K;++k){double v=avg(k,candidate.map(k));mx=std::max(mx,v);total+=v;}
    bool lex=best.n_elem==K && std::lexicographical_compare(candidate.map.begin(),candidate.map.end(),best.begin(),best.end());
    if(count>best_count || (count==best_count && (mx<best_max || (mx==best_max && (total<best_sum || (total==best_sum && lex)))))) {
      best_count=count;best_max=mx;best_sum=total;best=candidate.map;
    }
  }
  return FreezeChoice{best,avg,double(best_count)/records.size()};
}
