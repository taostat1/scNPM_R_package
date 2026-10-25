# Final manuscript settings: cells per condition, 10 replicates, eight methods.
perms3 <- rbind(c(1,2,3),c(1,3,2),c(2,1,3),c(2,3,1),c(3,1,2),c(3,2,1))
METHODS <- c('scNPM','scNPM_noPPI','scNPM_independent_fixedNPM',
             'BLGGM','glasso','mixtureGGM','AJGM','JGL_fused')
CFG <- list(num_iter=5000L,num_save=1000L,adapt_min=100L,adapt_max=3000L,
  freeze_patience=5L,final_burnin=1000L,bfdr=.05,threads=8L,
  lambda_grid=exp(seq(log(.6),log(.02),length.out=10)),lambda2_grid=c(.01,.05,.1),ebic_gamma=.5,
  mixture_lambda=c(.01,.03,.1),mixture_starts=3L,mixture_maxiter=150L,mixture_tol=1e-5,
  ajgm_python=Sys.getenv('SCNPM_AJGM_PYTHON','python3'),
  wsl_distro=Sys.getenv('SCNPM_AJGM_WSL',''),wsl_user=Sys.getenv('SCNPM_AJGM_USER',''),
  ajgm_maxiter=100L,ajgm_edge_cut=.05)
DATA_CFG <- list(seed=68261018L,n_per_condition=3000L,proportions=c(.4,.3,.3),n_replicates=10L,
  zero_target=.40,change_fraction=.20,change_mode='replace',ppi_coverage=.60,ppi_false_fraction=.10,
  min_eigenvalue=.10,max_proposals=20000L,control_label_permutation=c(2L,3L,1L),gamma_shape=5,gamma_rate=2)
replicate_seeds <- function(r,n) {
  network <- 68261018L+1000L*r
  base <- if(n==3000L)network else 73261020L+1000L*r
  list(network=network,means=network+100L,dropout=network+200L,
       expression=base+200L,ppi=network+300L,fit=base+500L)
}
