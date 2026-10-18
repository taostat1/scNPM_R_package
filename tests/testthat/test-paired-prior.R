# Independent R formulae, not calls back into the native probability helper.
paired_reference_weights <- function(pi, base=.05, d=.75, reference="control") {
  a<-d*base*(1-base);q<-c(1-base-a,a,a,base-a)
  odds<-(pi/(1-pi))/(base/(1-base))
  if(reference %in% c("control","both"))q[c(3,4)]<-q[c(3,4)]*odds
  if(reference %in% c("treatment","both"))q[c(2,4)]<-q[c(2,4)]*odds
  q/sum(q)
}
paired_reference_cov <- function(h,v0=.02,v1=1,r=.995) {
  ss<-switch(as.character(h),"0"=c(v0,v0),"1"=c(v1,v0),"2"=c(v0,v1),"3"=c(v1,v1))
  z<-diag(ss^2);if(h==3)z[1,2]<-z[2,1]<-r*v1^2;z
}
paired_reference_post <- function(w,pi,reference="control",d=.75,r=.995) {
  lp<-log(paired_reference_weights(pi,d=d,reference=reference))+
    vapply(0:3,function(h){v<-paired_reference_cov(h,r=r)
      -.5*as.numeric(determinant(v,logarithm=TRUE)$modulus)-.5*drop(t(w)%*%solve(v,w))},numeric(1))
  p<-exp(lp-max(lp));p/sum(p)
}
paired_fixture <- function(K=1L) {
  set.seed(1100);n<-20*K;g<-4L
  x<-exp(matrix(rnorm(g*n,sd=.1),g,n)+matrix(rep(rep(seq_len(K)*2,each=20),each=g),g,n))
  x<-cbind(x,x);rownames(x)<-paste0("g",seq_len(g));x[1,c(1,n+1)]<-0
  list(goal_gene_mt=x,cell_t=cbind(condition=rep(c(1,0),each=n)),K=K,
       protein_input=data.frame(node1="g1",node2="g2",combined_score=1),
       num_iter=75L,num_save=40L,adapt_min=5L,adapt_max=15L,freeze_patience=3L,
       final_burnin=15L,fdr_min_draws=20L,diagnostic_every=1L,
       num_step_theta=2L,num_step_lam=2L,epsilon_theta=.01,epsilon_lam=.01)
}

test_that("paired probabilities use the full four-state conditional, for each PPI source", {
  for(ref in c("control","treatment","both")) {
    a<-paired_fixture();a$ppi_reference<-ref;f<-do.call(run_scNPM,a)
    expect_identical(f$network_prior,"paired")
    expect_false(f$paired_prior$modules)
    expect_identical(f$paired_prior$observation_sampler,"fixed-step HMC")
    expect_identical(f$provenance$settings$ppi_reference,ref)
    for(i in 1:3)for(j in (i+1):4) {
      p<-vapply(seq_along(f$Result$invcov_treat_iter),function(s)
        paired_reference_post(c(f$Result$invcov_treat_iter[[s]][i,j,1],
          f$Result$invcov_ctrl_iter[[s]][i,j,1]),f$edge_prior[i,j],ref),numeric(4))
      actual<-vapply(c("p00","p10","p01","p11"),function(n)f$joint_edge_probabilities[[n]][i,j,1],numeric(1))
      expect_equal(unname(actual),rowMeans(p),tolerance=2e-10)
    }
    for(w in c(f$Result$invcov_treat_iter,f$Result$invcov_ctrl_iter)) {
      expect_equal(w[,,1],t(w[,,1]),tolerance=1e-12)
      expect_gt(min(eigen(w[,,1],symmetric=TRUE,only.values=TRUE)$values),0)
    }
    expect_true(all(f$stage_history$feedback_active[f$stage_history$iteration>f$sampling_schedule$frozen_at]==0))
    expect_true(all(is.finite(f$criterion_trace$network_prior_kernel)))
    expect_lte(f$differential_fdr$estimated_error_bound,.05+1e-12)
  }
})

test_that("paired updates keep all-pair marginals and frozen correspondence consistent", {
  a<-paired_fixture(2L);f<-do.call(run_scNPM,a);z<-f$joint_edge_probabilities
  for(i in 1:3)for(j in (i+1):4) {
    expect_equal(unname((z$p00+z$p10+z$p01+z$p11)[i,j,]),rep(1,4),tolerance=1e-12)
    pt<-matrix((z$p10+z$p11)[i,j,],2,byrow=TRUE)
    pc<-matrix((z$p01+z$p11)[i,j,],2,byrow=TRUE)
    expect_equal(pt[,1],pt[,2],tolerance=1e-12)
    expect_equal(pc[1,],pc[2,],tolerance=1e-12)
  }
  post<-f$matching_history[(f$sampling_schedule$frozen_at+1):a$num_iter,,drop=FALSE]
  expect_true(all(sweep(post,2,as.numeric(post[1,]),"-")==0))
  expect_match(f$sampling_schedule$labeling,"joint matched-pair")
  expect_identical(f$group_fdr$status,"estimated")
})

test_that("paired integrated column sampler agrees with independent 16-state enumeration", {
  z<-getFromNamespace("paired_column_diagnostic","scNPM")(25000L)
  hs<-as.matrix(expand.grid(0:3,0:3));lp<-numeric(16);mus<-matrix(0,16,4);vs<-vector("list",16)
  for(k in seq_len(16)) {
    P<-z$likelihood_precision;score<-0
    for(j in 1:2){ix<-c(j,j+2);v<-paired_reference_cov(hs[k,j],.2,.8,.8)
      J<-solve(v);P[ix,ix]<-P[ix,ix]+J
      score<-score+log(paired_reference_weights(z$ppi_weights[j])[hs[k,j]+1])-.5*as.numeric(determinant(v,TRUE)$modulus)}
    V<-solve(P);mu<-drop(-V%*%z$linear_term);mus[k,]<-mu;vs[[k]]<-V
    lp[k]<-score-.5*as.numeric(determinant(P,TRUE)$modulus)+.5*sum(mu*(P%*%mu))
  }
  p<-exp(lp-max(lp));p<-p/sum(p);expected<-colSums(mus*p)
  empirical<-tabulate(1+z$states[,1]+4*z$states[,2],nbins=16)/nrow(z$states)
  expect_lt(max(abs(empirical-p)),.025)
  # Absolute Monte Carlo tolerance, not a relative tolerance around a near-zero mean.
  expect_lt(max(abs(colMeans(z$coefficients)-expected)),.025)
  second<-Reduce(`+`,lapply(seq_len(16),function(k)p[k]*(vs[[k]]+tcrossprod(mus[k,]))))
  expect_lt(max(abs(crossprod(z$coefficients)/nrow(z$coefficients)-second)),.025)
})

test_that("old independent prior remains selectable and invalid paired parameters fail", {
  a<-paired_fixture();a$network_prior<-"independent";f<-do.call(run_scNPM,a)
  expect_identical(f$network_prior,"independent");expect_null(f$paired_prior)
  a$network_prior<-"paired";a$shared_edge_correlation<-1
  expect_error(do.call(run_scNPM,a),"shared_edge_correlation")
  a$shared_edge_correlation<-.995;a$paired_discordance<-0
  expect_error(do.call(run_scNPM,a),"paired_discordance")
  a$paired_discordance<-.75;a$inference_mode<-"coupled"
  expect_error(do.call(run_scNPM,a),"requires staged")
})
