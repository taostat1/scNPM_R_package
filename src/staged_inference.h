#ifndef SCNPM_STAGED_INFERENCE_H
#define SCNPM_STAGED_INFERENCE_H

inline void align_components(const mat& anchor, mat& mu, cube& precision,
                             cube& covariance, cube& edges, vec& pi, vec& group) {
  const uword K=mu.n_cols;
  mat cost(K,K);
  for (uword a=0; a<K; ++a) for (uword k=0; k<K; ++k)
    cost(a,k)=accu(square(anchor.col(a)-mu.col(k)));
  const uvec order=minimum_assignment(cost);
  uvec inverse(K);
  for (uword a=0; a<K; ++a) inverse(order(a))=a;
  mu=mu.cols(order); pi=pi.elem(order);
  const cube old_p=precision, old_c=covariance, old_e=edges;
  for (uword a=0; a<K; ++a) {
    precision.slice(a)=old_p.slice(order(a)); covariance.slice(a)=old_c.slice(order(a));
    edges.slice(a)=old_e.slice(order(a));
  }
  for (uword i=0; i<group.n_elem; ++i) group(i)=inverse(static_cast<uword>(group(i)));
}

// Complete-data augmented likelihood, including the positive-observation
// Jacobian. This is not the intractable observed-data marginal likelihood.
inline double augmented_likelihood(const mat& theta, const mat& zero,
  const mat& mt, const mat& mc, const cube& ot, const cube& oc,
  const vec& gt, const vec& gc, const uvec& ti, const uvec& ci,
  const vec& pit, const vec& pic, const vec& l0, const vec& l1) {
  const double log2pi=std::log(2.0*arma::datum::pi);
  double out=0;
  for (int d=0; d<2; ++d) {
    const mat& mu=d==0?mt:mc; const cube& om=d==0?ot:oc;
    const vec& grp=d==0?gt:gc; const uvec& idx=d==0?ti:ci;
    const vec& pi=d==0?pit:pic;
    for (uword k=0; k<mu.n_cols; ++k) {
      const uvec local=find(grp==k);
      if (local.is_empty()) continue;
      const mat chol_o=chol(om.slice(k));
      const double logdet=2*sum(log(chol_o.diag()));
      mat centered=theta.rows(idx.elem(local)); centered.each_row()-=mu.col(k).t();
      out+=local.n_elem*(std::log(pi(k))+0.5*logdet-0.5*theta.n_cols*log2pi)
        -0.5*accu((centered*om.slice(k))%centered);
    }
  }
  for (uword g=0; g<theta.n_cols; ++g) {
    const vec eta=l0(g)+l1(g)*theta.col(g);
    const vec log_zero=rcpp_log_normcdf(eta), log_obs=rcpp_log_normsf(eta);
    for (uword i=0; i<theta.n_rows; ++i)
      out+=zero(i,g)>0.5?log_zero(i):log_obs(i)-theta(i,g);
  }
  return out;
}

// A state-wise diagnostic kernel. During adaptation the adjusted factors need
// not be compatible joint conditionals; J_aug is never used to accept/reject.
inline double network_prior_kernel(const cube& ot, const cube& oc,
  const cube& zt, const cube& zc, const mat& base, const uvec& map,
  double v0, double v1, double rate2, double threshold, double adjustment) {
  const uword K=ot.n_slices, G=ot.n_rows;
  uvec inverse(K); inverse.fill(K);
  for (uword k=0;k<K;++k) if(map(k)<K) inverse(map(k))=k;
  double out=0;
  for(int d=0;d<2;++d) {
    const cube& om=d==0?ot:oc; const cube& z=d==0?zt:zc; const cube& other=d==0?zc:zt;
    const uvec& partners=d==0?map:inverse;
    for(uword k=0;k<K;++k) {
      for(uword i=0;i<G;++i) {
        out+=std::log(rate2/2)-rate2*om(i,i,k)/2;
        for(uword j=i+1;j<G;++j) {
          const double xi=base(i,j), w=om(i,j,k); double a=xi;
          const double p=calculate_edge_probability(w,v0,v1,xi);
          if(partners(k)<K && p>=1-threshold && p<=threshold) {
            a=other(i,j,partners(k))>0.5 ? std::max(xi*adjustment,std::numeric_limits<double>::min())
              : std::min(xi/adjustment,std::max(xi,1-1e-6));
          }
          const bool included=z(i,j,k)>0.5;
          const double sd=included?v1:v0;
          out+=-0.5*std::log(2*arma::datum::pi)-std::log(sd)-0.5*std::pow(w/sd,2)
            +(included?std::log(a):std::log1p(-a));
        }
      }
    }
  }
  return out;
}

// Rao-Blackwellize within each draw, never multiply marginal edge means.
// Slice k*K+c refers to treatment k / control c, including non-frozen pairs.
inline void accumulate_joint_probabilities(const cube& ot,const cube& oc,
  const cube& zt,const cube& zc,const mat& base,double v0,double v1,
  cube& p00,cube& p10,cube& p01,cube& p11,cube& e10,cube& e01,
  cube& b10,cube& b01,vec& sizes,int batch) {
  const uword K=ot.n_slices,G=ot.n_rows;
  uword edge=0;
  for(uword i=0;i<G;++i) for(uword j=i+1;j<G;++j,++edge) {
    vec pt(K),pc(K);
    for(uword k=0;k<K;++k) {
      pt(k)=calculate_edge_probability(ot(i,j,k),v0,v1,base(i,j));
      pc(k)=calculate_edge_probability(oc(i,j,k),v0,v1,base(i,j));
    }
    for(uword k=0;k<K;++k) for(uword c=0;c<K;++c) {
      uword h=k*K+c;
      const double a=(1-pt(k))*(1-pc(c)), b=pt(k)*(1-pc(c));
      const double f=(1-pt(k))*pc(c), g=pt(k)*pc(c);
      p00(i,j,h)+=a;p10(i,j,h)+=b;p01(i,j,h)+=f;p11(i,j,h)+=g;
      p00(j,i,h)+=a;p10(j,i,h)+=b;p01(j,i,h)+=f;p11(j,i,h)+=g;
      const double eb=zt(i,j,k)>0.5 && zc(i,j,c)<0.5?1:0;
      const double ef=zt(i,j,k)<0.5 && zc(i,j,c)>0.5?1:0;
      e10(i,j,h)+=eb;e10(j,i,h)+=eb;e01(i,j,h)+=ef;e01(j,i,h)+=ef;
      b10(edge,h,batch)+=b;b01(edge,h,batch)+=f;
    }
  }
  sizes(batch)+=1;
}
#endif
