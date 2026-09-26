staged_single_thread <- function(code) {
  set <- getFromNamespace("set_scNPM_threads","scNPM")
  old <- set(1L); on.exit(set(old),add=TRUE); force(code)
}
# Final-stage tests use an exactly balanced expression fixture: K=1 has a
# complete map independent of numerical state, so phase boundaries are known.
staged_fixture <- function(K=1L) {
  set.seed(28)
  G <- max(4L,K); n <- 12L*K
  mu <- if(K==1L) matrix(0,G,1) else 3*diag(G)[,seq_len(K),drop=FALSE]
  labels <- rep(seq_len(K),each=12L)
  theta <- t(mu[,labels])+matrix(rnorm(n*G,sd=.12),n,G)
  one <- t(exp(theta)); x <- cbind(one,one)
  rownames(x) <- paste0("G",seq_len(G))
  list(goal_gene_mt=x,cell_t=cbind(condition=rep(c(1,0),each=n)),K=K,
       num_iter=70L,num_save=40L,adapt_min=6L,adapt_max=20L,
       freeze_patience=4L,final_burnin=10L,diagnostic_every=1L,
       fdr_min_draws=20L,num_step_theta=2L,num_step_lam=2L,
       epsilon_theta=.01,epsilon_lam=.01,ds_seed=918L)
}
staged_run <- function(args) {
  out <- NULL
  invisible(capture.output(out <- suppressMessages(do.call(run_scNPM,args))))
  out
}

fdr_fixture <- function(K=1L,G=4L) {
  a <- array(0,c(G,G,K*K))
  list(K=K,G=G,gene_names=paste0("G",seq_len(G)),inference_mode="staged",
    matching_result_df=data.frame(treatment_cluster=seq_len(K),matched_control_cluster=seq_len(K)),
    sampling_schedule=list(feedback_disabled_for_retained=TRUE,retained_draws=500L),
    joint_edge_probabilities=list(p00=a,p10=a,p01=a,p11=a))
}


# Independent enumerating oracle for small complete assignment problems.
all_permutations <- function(x) {
  if(length(x)==1L) return(matrix(x,1L))
  do.call(rbind,lapply(x,function(a)cbind(a,all_permutations(x[x!=a]))))
}
assignment_oracle <- function(D) {
  ps <- all_permutations(seq_len(nrow(D)))
  obj <- t(apply(ps,1,function(p) {v<-D[cbind(seq_len(nrow(D)),p)];c(max(v),sum(v))}))
  as.numeric(ps[order(obj[,1],obj[,2])[1],])
}

test_that("lexicographic matching prioritizes bottleneck, then total, then fixed ties", {
  fn <- getFromNamespace("lexicographic_assignment","scNPM")
  # Sum minimization picks (1,2): sum 10, max 9; minimax picks (2,1): sum 12, max 6.
  D <- matrix(c(1,6,6,9),2,byrow=TRUE)
  expect_equal(as.numeric(fn(D)),c(2,1))
  D <- matrix(c(4,4,1,3),2,byrow=TRUE)
  expect_equal(as.numeric(fn(D)),c(2,1)) # both max 4; sums 5 versus 7
  expect_equal(as.numeric(fn(matrix(1,4,4))),1:4)
  set.seed(59)
  for(n in 2:6) for(rep in 1:8) {
    D <- matrix(sample(0:15,n*n,TRUE),n)
    expect_equal(as.numeric(fn(D)),assignment_oracle(D))
  }
  expect_equal(as.numeric(fn(matrix(c(Inf,1,1,Inf),2))),c(2,1))
  expect_equal(as.numeric(fn(matrix(Inf,4,4))),1:4)
  expect_error(fn(matrix(NA_real_,2,2)),"without NA")
  expect_error(fn(matrix(-1,2,2)),"nonnegative")
})

test_that("large-K solver enforces the bottleneck before minimum sum", {
  fn <- getFromNamespace("lexicographic_assignment","scNPM")
  set.seed(73)
  for(rep in 1:10) {
    block <- matrix(sample(0:20,16,TRUE),4)
    D <- matrix(100,9,9);D[1:4,1:4] <- block
    for(k in 5:9) D[k,k] <- 0
    got <- as.numeric(fn(D));want <- c(assignment_oracle(block),5:9)
    vg <- D[cbind(1:9,got)]; vw <- D[cbind(1:9,want)]
    expect_equal(c(max(vg),sum(vg)),c(max(vw),sum(vw)))
    expect_equal(sort(got),1:9)
    expect_equal(as.numeric(fn(D)),got)
  }
  expect_equal(as.numeric(fn(matrix(0,9,9))),1:9)
  expect_equal(as.numeric(fn(matrix(Inf,9,9))),1:9)
})

make_threshold_fit <- function(delta,error) {
  # G=4 creates six unordered edges, with remaining directions undetermined.
  f <- fdr_fixture(); delta <- c(delta,rep(0,6-length(delta)))
  error <- c(error,rep(1,6-length(error)))
  e<-0L
  for(i in 1:3) for(j in (i+1):4) {
    e<-e+1L
    plus<-if(delta[e]>0) 1-error[e] else 0
    minus<-plus-delta[e]
    vals<-c(p00=1-plus-minus,p10=plus,p01=minus,p11=0)
    stopifnot(all(vals>=-1e-12))
    for(nm in names(vals)) f$joint_edge_probabilities[[nm]][i,j,1] <-
      f$joint_edge_probabilities[[nm]][j,i,1] <- vals[nm]
  }
  f
}

test_that("delta threshold matches hand calculations and never splits equal contrasts", {
  f<-make_threshold_fit(c(.98,.96,.92,.88,.69),c(.0199,.0397,.0793,.1189,.307))
  z<-select_differential_edges(f)
  expect_equal(z$scope,"frozen")
  expect_equal(z$rule,"delta_threshold")
  expect_equal(z$threshold,.88,tolerance=1e-12)
  expect_equal(nrow(z$selected),3)
  expect_equal(z$estimated_error_bound,.0463,tolerance=1e-12)
  expect_identical(z$table$selected,z$table$abs_delta_p>z$threshold)
  expect_equal(sum(z$threshold_curve$chosen),1)
  # Exact score tie at 0.5, different errors: an error-prefix can select one,
  # but no single strict contrast threshold can split the two candidates.
  f<-make_threshold_fit(c(.5,.5),c(.25,.5))
  a<-select_differential_edges(f,fdr_level=.3)
  b<-select_differential_edges(f,fdr_level=.3,rule="posterior_error")
  expect_equal(nrow(a$selected),0);expect_identical(a$threshold,Inf)
  expect_equal(nrow(b$selected),1)
  expect_equal(nrow(a$threshold_curve),1)
  expect_equal(a$table$prefix_mean_error[1:2],rep(.375,2))
})

test_that("all cutoffs are scanned even when the first prefix fails", {
  f<-make_threshold_fit(c(.875,.75,.625),c(.125,.125,.1875))
  # Construct a nonmonotone curve: higher contrast can have higher error.
  f<-make_threshold_fit(c(.9,.89,.88),c(.1,.055,.06))
  z<-select_differential_edges(f,fdr_level=.08)
  expect_false(z$threshold_curve$feasible[1])
  expect_true(z$threshold_curve$feasible[2])
  expect_equal(nrow(z$selected),3)
  expect_equal(z$threshold,0)
  f$sampling_schedule$retained_draws<-30L
  z<-select_differential_edges(f,fdr_level=.08)
  expect_equal(z$status,"insufficient_draws")
  expect_equal(nrow(z$selected),0)
  expect_identical(z$threshold,Inf)
  expect_false(any(z$threshold_curve$chosen))
})

test_that("default public matching freezes a complete globally optimal assignment", {
  staged_single_thread({
    a<-staged_fixture();f<-staged_run(a)
    expect_equal(f$matching_method,"lexicographic")
    expect_true(f$matching_objective$frozen)
    expect_equal(dim(f$smd_cost_matrix),c(1L,1L))
    expect_equal(f$matching_objective$max_smd,max(f$smd_results))
    expect_equal(f$matching_objective$total_smd,sum(f$smd_results))
    expect_equal(f$matching_result_df$matched_control_cluster,
      as.numeric(getFromNamespace("lexicographic_assignment","scNPM")(f$frozen_smd_cost_matrix)))
    expect_identical(f$provenance$settings$smd_cap,NULL)
    expect_equal(f$provenance$settings$fdr_rule,"delta_threshold")
  })
})
