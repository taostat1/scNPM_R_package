# Pairwise propensity matching: condition labels are observed responses, never imputed.
# Does not touch R's RNG, so cross-fitting cannot reseed the MCMC sampler.
.pairwise_ridge_predict <- function(X, y, newX, lambda=.001) {
 center<-colMeans(X);ss<-apply(X,2,stats::sd)
 ss[!is.finite(ss)|ss<1e-10]<-1
 A<-cbind(1,sweep(sweep(X,2,center,"-"),2,ss,"/"))
 B<-cbind(1,sweep(sweep(newX,2,center,"-"),2,ss,"/"))
 weights<-ifelse(y==1,.5/mean(y==1),.5/mean(y==0))
 penalty<-c(1e-8,rep(lambda*sum(weights),ncol(X)))
 beta<-numeric(ncol(A));converged<-FALSE
 for(it in 1:100) {
  p<-stats::plogis(pmax(-30,pmin(30,drop(A%*%beta))))
  w<-pmax(p*(1-p),1e-10)*weights
  H<-crossprod(A,A*w)+diag(penalty)
  step<-drop(solve(H,crossprod(A,weights*(y-p))-penalty*beta))
  beta<-beta+step
  if(max(abs(step))<1e-7){converged<-TRUE;break}
 }
 if(!converged || any(!is.finite(beta))) {
  # Preserve the original fast path. Restart an oscillating Newton fit on the
  # same penalized likelihood, using Armijo backtracking instead of full steps.
  # Labels, data, ridge strength and the sampler random-number state are unchanged.
  objective<-function(b) {
   eta<-drop(A%*%b)
   sum(weights*(pmax(eta,0)+log1p(exp(-abs(eta)))-y*eta))+sum(penalty*b*b)/2
  }
  beta<-numeric(ncol(A));value<-objective(beta)
  for(it in 1:100) {
   p<-stats::plogis(drop(A%*%beta))
   gradient<-drop(crossprod(A,weights*(p-y)))+penalty*beta
   if(max(abs(gradient))/sum(weights)<1e-10)break
   w<-pmax(p*(1-p),1e-10)*weights
   H<-crossprod(A,A*w)+diag(penalty)
   direction<-drop(solve(H,gradient));decrement<-sum(gradient*direction)
   if(!all(is.finite(direction)) || !is.finite(decrement) || decrement<=0)break
   alpha<-1;accepted<-FALSE
   for(backtrack in 0:50) {
    candidate<-beta-alpha*direction;candidate_value<-objective(candidate)
    slack<-10*.Machine$double.eps*max(1,abs(value))
    if(is.finite(candidate_value) && candidate_value<=value-1e-4*alpha*decrement+slack) {
     accepted<-TRUE;break
    }
    alpha<-alpha/2
   }
   if(!accepted)break
   beta<-candidate;value<-candidate_value
  }
  # Verify stationarity, including a solution reached on the final step.
  gradient<-drop(crossprod(A,weights*(stats::plogis(drop(A%*%beta))-y)))+penalty*beta
  converged<-all(is.finite(beta)) && all(is.finite(gradient)) &&
    max(abs(gradient))/sum(weights)<1e-10
 }
 if(!converged || any(!is.finite(beta)))stop("Pairwise logistic solver did not converge; inspect numerical inputs")
 pmin(1-1e-8,pmax(1e-8,stats::plogis(pmax(-30,pmin(30,drop(B%*%beta))))))
}
.pairwise_propensity_cost <- function(X,condition,group_treat,group_ctrl,treat_idx,ctrl_idx) {
 X<-as.matrix(X);condition<-as.numeric(condition)
 gt<-as.integer(group_treat)+1L;gc<-as.integer(group_ctrl)+1L
 ti<-as.integer(treat_idx)+1L;ci<-as.integer(ctrl_idx)+1L
 # The native caller adds K as an attribute below; do not infer it from missing labels.
 K<-attr(X,"matching_K")
 if(is.null(K))stop("Missing matching dimension")
 cost<-matrix(NaN,K,K);scores<-matrix(NaN,nrow(X),K)
 for(k in seq_len(K))for(l in seq_len(K)) {
  a<-ti[gt==k];b<-ci[gc==l]
  if(!length(a)||!length(b))next
  ii<-c(a,b);y<-condition[ii];v<-X[ii,,drop=FALSE]
  nf<-min(5L,length(a),length(b))
  if(nf<2L)pred<-.pairwise_ridge_predict(v,y,v)else {
   # Cell-index fold allocation is invariant to arbitrary cluster numbering.
   folds<-((ii-1L)%%nf)+1L;pred<-numeric(length(ii))
   for(f in sort(unique(folds))) {
    test<-folds==f;train<-!test
    if(length(unique(y[train]))<2L) {
     pred[test]<-.pairwise_ridge_predict(v,y,v[test,,drop=FALSE])
    } else pred[test]<-.pairwise_ridge_predict(v[train,,drop=FALSE],y[train],v[test,,drop=FALSE])
   }
  }
  va<-if(length(a)>1)stats::var(pred[y==1])else 0
  vb<-if(length(b)>1)stats::var(pred[y==0])else 0
  cost[k,l]<-abs(mean(pred[y==1])-mean(pred[y==0]))/sqrt(max((va+vb)/2,1e-12))
  scores[a,l]<-pred[y==1];scores[b,k]<-pred[y==0]
 }
 list(cost=cost,scores=scores)
}
