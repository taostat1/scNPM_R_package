#ifndef SCNPM_PAIRED_PRIOR_H
#define SCNPM_PAIRED_PRIOR_H

// Four support-state factors (00,10,01,11), in one globally normalized prior
// restricted to positive-definite precision matrices. These are not marginal
// edge probabilities after that restriction.
inline vec paired_weights(double pi, double base, double discordance, int reference) {
  const double d = discordance * base * (1-base);
  vec q = {1-base-d, d, d, base-d};
  const double odds = (pi/(1-pi))/(base/(1-base));
  if (reference == 0 || reference == 2) { q(2)*=odds; q(3)*=odds; }
  if (reference == 1 || reference == 2) { q(1)*=odds; q(3)*=odds; }
  return q/sum(q);
}
inline mat paired_precision(int h, double v0, double v1, double correlation) {
  const double a=(h==1||h==3)?v1:v0, b=(h==2||h==3)?v1:v0;
  const double r=h==3?correlation:0, den=1-r*r;
  mat J(2,2); J(0,0)=1/(a*a*den); J(1,1)=1/(b*b*den);
  J(0,1)=J(1,0)=-r/(a*b*den); return J;
}
inline double paired_logdet(const mat& A) {
  const double d=A(0,0)*A(1,1)-A(0,1)*A(1,0);
  if (!(d>0) || !std::isfinite(d)) throw std::runtime_error("Invalid paired precision determinant.");
  return std::log(d);
}
inline mat paired_inverse(const mat& A) {
  const double d=std::exp(paired_logdet(A));
  mat V(2,2); V(0,0)=A(1,1)/d; V(1,1)=A(0,0)/d;
  V(0,1)=V(1,0)=-A(0,1)/d; return V;
}
inline vec paired_softmax(vec lp) {
  lp-=lp.max(); vec p=exp(lp); p/=sum(p);
  if (!p.is_finite()) throw std::runtime_error("Non-finite four-state probabilities.");
  return p;
}
inline vec paired_posterior(double wt,double wc,double pi,double base,double discordance,
                           double v0,double v1,double correlation,int reference) {
  vec w={wt,wc}, lp=log(paired_weights(pi,base,discordance,reference));
  for(int h=0;h<4;++h) {
    mat J=paired_precision(h,v0,v1,correlation);
    lp(h)+=.5*paired_logdet(J)-.5*dot(w,J*w);
  }
  return paired_softmax(lp);
}
inline int paired_draw(const vec& p,std::mt19937& rng) {
  const double u=thread_safe_randu(rng); double s=0;
  for(int h=0;h<3;++h) { s+=p(h); if(u<s) return h; } return 3;
}

// Integrate one T/C coefficient pair when updating its state, then immediately
// restore both coefficients. Other column entries stay conditioned upon.
inline vec paired_column(mat& H,vec& zt,vec& zc,const vec& s,const vec& prior,
                         double base,double discordance,double v0,double v1,
                         double correlation,int reference,std::mt19937& rng,vec beta) {
  const uword p=zt.n_elem; vec residual=-s-H*beta;
  mat Js[4]; double ld[4];
  for(int h=0;h<4;++h) { Js[h]=paired_precision(h,v0,v1,correlation);ld[h]=paired_logdet(Js[h]); }
  for(uword j=0;j<p;++j) {
    uvec ix={j,p+j}; const int old=(zt(j)>.5?1:0)+(zc(j)>.5?2:0);
    mat Hjj=H.submat(ix,ix),L=Hjj-Js[old]; vec r=residual(ix)+Hjj*beta(ix);
    vec lp=log(paired_weights(prior(j),base,discordance,reference)); mat Vs[4];vec ms[4];
    for(int h=0;h<4;++h) {
      mat P=L+Js[h];Vs[h]=paired_inverse(P);ms[h]=Vs[h]*r;
      lp(h)+=.5*ld[h]-.5*paired_logdet(P)+.5*dot(r,ms[h]);
    }
    const int h=paired_draw(paired_softmax(lp),rng);
    vec z={thread_safe_rnorm(rng),thread_safe_rnorm(rng)};
    vec draw=ms[h]+chol(Vs[h]).t()*z,previous=beta(ix);
    mat delta=Js[h]-Js[old]; H.submat(ix,ix)+=delta;residual(ix)-=delta*previous;
    residual-=H.cols(ix)*(draw-previous);beta(ix)=draw;
    zt(j)=h==1||h==3;zc(j)=h==2||h==3;
  }
  return beta;
}
inline void paired_restore(mat& W,mat& Sigma,mat& Z,const mat& Wi,const vec& beta,
                           const vec& z,const uvec& other,int g,int n,double Sgg,
                           double l,std::mt19937& rng) {
  std::gamma_distribution<double> gamma(n/2.0+1.0,2.0/(Sgg+l));const double schur=gamma(rng);
  if (!(schur>0)) throw std::runtime_error("Nonpositive paired Schur complement.");
  vec t=Wi*beta,wcol(W.n_rows);wcol(other)=beta;wcol(g)=schur+dot(beta,t);
  W.col(g)=wcol;W.row(g)=wcol.t();
  Sigma.submat(other,other)=Wi+t*t.t()/schur;
  vec cc(W.n_rows);cc(other)=-t/schur;cc(g)=1/schur;Sigma.col(g)=cc;Sigma.row(g)=cc.t();
  for(uword j=0;j<other.n_elem;++j)Z(other(j),g)=Z(g,other(j))=z(j);
}
inline void paired_graph(cube& ot,cube& ct,cube& et,cube& oc,cube& cc,cube& ec,
 const mat& xt,const mat& mt,const vec& gt,const mat& xc,const mat& mc,const vec& gc,
 const mat& prior,const uvec& map,double base,double discordance,double v0,double v1,
 double correlation,int reference,double l,int G,int K,unsigned long seed) {
  std::vector<std::string> errors(K);
  #pragma omp parallel for schedule(static) num_threads(std::min(omp_get_max_threads(),K))
  for(int k=0;k<K;++k) {
    try {
      std::mt19937 rng(seed+10000UL*k);const int c=map(k);
      uvec it=find(gt==k),ic=find(gc==c);
      mat Xt=xt.rows(it).t(),Xc=xc.rows(ic).t();Xt.each_col()-=mt.col(k);Xc.each_col()-=mc.col(c);
      mat St=Xt*Xt.t(),Sc=Xc*Xc.t(),Wt=ot.slice(k),Wc=oc.slice(c),Sigt=ct.slice(k),Sigc=cc.slice(c),Zt=et.slice(k),Zc=ec.slice(c);
      uvec all=regspace<uvec>(0,G-1);const int p=G-1;
      for(int g=0;g<G;++g) {
        uvec other=all(find(all!=(uword)g));vec cgT=Sigt.col(g),cgC=Sigc.col(g);
        mat Wit=Sigt.submat(other,other)-cgT(other)*cgT(other).t()/cgT(g);
        mat Wic=Sigc.submat(other,other)-cgC(other)*cgC(other).t()/cgC(g);
        vec zt=Zt.col(g),zc=Zc.col(g);zt=zt(other);zc=zc(other);
        mat H(2*p,2*p,fill::zeros);
        H.submat(0,0,p-1,p-1)=(St(g,g)+l)*Wit;H.submat(p,p,2*p-1,2*p-1)=(Sc(g,g)+l)*Wic;
        for(int j=0;j<p;++j){uvec ix={(uword)j,(uword)(p+j)};H.submat(ix,ix)+=paired_precision((zt(j)>.5?1:0)+(zc(j)>.5?2:0),v0,v1,correlation);}
        vec st=St.col(g),sc=Sc.col(g),pp=prior.col(g),bt=Wt.col(g),bc=Wc.col(g);
        vec beta=paired_column(H,zt,zc,join_cols(st(other),sc(other)),pp(other),base,discordance,v0,v1,correlation,reference,rng,join_cols(bt(other),bc(other)));
        paired_restore(Wt,Sigt,Zt,Wit,beta.head(p),zt,other,g,it.n_elem,St(g,g),l,rng);
        paired_restore(Wc,Sigc,Zc,Wic,beta.tail(p),zc,other,g,ic.n_elem,Sc(g,g),l,rng);
      }
      ot.slice(k)=Wt;ct.slice(k)=Sigt;et.slice(k)=Zt;oc.slice(c)=Wc;cc.slice(c)=Sigc;ec.slice(c)=Zc;
    } catch(const std::exception& ex) { errors[k]=ex.what(); }
  }
  for(int k=0;k<K;++k)if(!errors[k].empty())Rcpp::stop("Paired network update: "+errors[k]);
}
inline void paired_support(cube& zt,cube& zc,const cube& ot,const cube& oc,const mat& prior,const uvec& map,
 double base,double discordance,double v0,double v1,double correlation,int reference,unsigned long seed) {
  std::mt19937 rng(seed);
  for(uword k=0;k<ot.n_slices;++k)for(uword i=0;i<ot.n_rows;++i)for(uword j=i+1;j<ot.n_rows;++j){
    uword c=map(k);int h=paired_draw(paired_posterior(ot(i,j,k),oc(i,j,c),prior(i,j),base,discordance,v0,v1,correlation,reference),rng);
    zt(i,j,k)=zt(j,i,k)=h==1||h==3;zc(i,j,c)=zc(j,i,c)=h==2||h==3;
  }
}
inline double paired_prior_kernel(const cube& ot,const cube& oc,const cube& zt,const cube& zc,
 const mat& prior,const uvec& map,double base,double discordance,double v0,double v1,double correlation,int reference,double l) {
  double value=0;
  for(uword k=0;k<ot.n_slices;++k)for(uword i=0;i<ot.n_rows;++i){
    value+=2*std::log(l/2)-l*(ot(i,i,k)+oc(i,i,map(k)))/2;
    for(uword j=i+1;j<ot.n_rows;++j){
      int h=(zt(i,j,k)>.5?1:0)+(zc(i,j,map(k))>.5?2:0);
      mat J=paired_precision(h,v0,v1,correlation);vec w={ot(i,j,k),oc(i,j,map(k))};
      value+=std::log(paired_weights(prior(i,j),base,discordance,reference)(h))-.5*dot(w,J*w)+.5*paired_logdet(J)-std::log(2*arma::datum::pi);
    }
  }return value;
}
inline void paired_reorder(const uvec& order,mat& mu,cube& w,cube& s,cube& z,vec& pi,vec& group){
  uvec inverse(order.n_elem);for(uword k=0;k<order.n_elem;++k)inverse(order(k))=k;
  mu=mu.cols(order);pi=pi.elem(order);cube ow=w,os=s,oz=z;
  for(uword k=0;k<order.n_elem;++k){w.slice(k)=ow.slice(order(k));s.slice(k)=os.slice(order(k));z.slice(k)=oz.slice(order(k));}
  for(uword i=0;i<group.n_elem;++i)group(i)=inverse((uword)group(i));
}
// Relabel whole matched pairs together; never independently change partners.
inline void paired_align(const mat& at,const mat& ac,mat& mt,mat& mc,cube& wt,cube& wc,cube& st,cube& sc,
 cube& zt,cube& zc,vec& pt,vec& pc,vec& gt,vec& gc,const uvec& map){
  uword K=map.n_elem;mat cost(K,K);for(uword k=0;k<K;++k)for(uword r=0;r<K;++r)
    cost(k,r)=accu(square(at.col(k)-mt.col(r)))+accu(square(ac.col(map(k))-mc.col(map(r))));
  uvec order=minimum_assignment(cost),co(K);for(uword k=0;k<K;++k)co(map(k))=map(order(k));
  paired_reorder(order,mt,wt,st,zt,pt,gt);paired_reorder(co,mc,wc,sc,zc,pc,gc);
}
inline void paired_accumulate(const cube& ot,const cube& oc,const cube& zt,const cube& zc,const mat& prior,
 const uvec& map,double base,double discordance,double v0,double v1,double correlation,int reference,
 cube& p00,cube& p10,cube& p01,cube& p11,cube& e10,cube& e01,cube& b10,cube& b01,vec& sizes,int batch){
  uword K=ot.n_slices,G=ot.n_rows,edge=0;
  for(uword i=0;i<G;++i)for(uword j=i+1;j<G;++j,++edge){
    mat pp(4,K);vec pt(K),pc(K);
    for(uword k=0;k<K;++k){pp.col(k)=paired_posterior(ot(i,j,k),oc(i,j,map(k)),prior(i,j),base,discordance,v0,v1,correlation,reference);pt(k)=pp(1,k)+pp(3,k);pc(map(k))=pp(2,k)+pp(3,k);}
    for(uword k=0;k<K;++k)for(uword c=0;c<K;++c){uword h=k*K+c;
      vec p={(1-pt(k))*(1-pc(c)),pt(k)*(1-pc(c)),(1-pt(k))*pc(c),pt(k)*pc(c)};
      if(c==map(k))p=pp.col(k);
      p00(i,j,h)+=p(0);p00(j,i,h)+=p(0);p10(i,j,h)+=p(1);p10(j,i,h)+=p(1);p01(i,j,h)+=p(2);p01(j,i,h)+=p(2);p11(i,j,h)+=p(3);p11(j,i,h)+=p(3);
      double a=zt(i,j,k)>.5&&zc(i,j,c)<.5,b=zt(i,j,k)<.5&&zc(i,j,c)>.5;
      e10(i,j,h)+=a;e10(j,i,h)+=a;e01(i,j,h)+=b;e01(j,i,h)+=b;
      b10(edge,h,batch)+=p(1);b01(edge,h,batch)+=p(2);
    }
  }sizes(batch)+=1;
}
#endif
