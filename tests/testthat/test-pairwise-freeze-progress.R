test_that("pairwise costs preserve labels, observed condition, and the sampler RNG", {
  cost <- getFromNamespace(".pairwise_propensity_cost","scNPM")
  set.seed(61)
  n <- 90L; k <- rep(1:3,each=30)
  a <- cbind(k*5,matrix(rnorm(n*5,sd=.1),n,5))
  X <- rbind(a,a); attr(X,"matching_K") <- 3L
  y <- rep(c(1,0),each=n)
  before <- .Random.seed
  z <- cost(X,y,k-1L,k-1L,0:(n-1L),n:(2*n-1L))
  expect_identical(.Random.seed,before)
  expect_equal(diag(z$cost),rep(0,3),tolerance=1e-8)
  expect_true(all(z$cost[row(z$cost)!=col(z$cost)]>1))
  p <- c(3L,1L,2L);q <- c(2L,3L,1L)
  zz <- cost(X,y,p[k]-1L,q[k]-1L,0:(n-1L),n:(2*n-1L))
  expect_equal(zz$cost[p,q],z$cost,tolerance=1e-10)
  expect_true(all(is.finite(z$scores)))
})

new_freeze_fixture <- function() {
  set.seed(31)
  x <- exp(matrix(rnorm(4*32,sd=.2),4,32))
  rownames(x) <- paste0("g",1:4)
  list(goal_gene_mt=x,cell_t=cbind(condition=rep(c(1,0),each=16)),K=1L,
    ds_seed=51,num_iter=40L,num_save=20L,adapt_min=3L,adapt_max=8L,
    freeze_patience=20L,final_burnin=10L,fdr_min_draws=20L,
    num_step_theta=2L,num_step_lam=2L,epsilon_theta=.01,epsilon_lam=.01)
}

test_that("consensus fallback completes sampling and remains distinct from stability", {
  threads <- getFromNamespace("set_scNPM_threads","scNPM")
  old <- threads(1L);on.exit(threads(old))
  a <- new_freeze_fixture()
  f <- do.call(run_scNPM,a)
  expect_equal(f$sampling_schedule$freeze_reason,"budget_consensus")
  expect_false(f$sampling_schedule$matching_stable)
  expect_equal(f$sampling_schedule$frozen_at,8)
  expect_gte(f$sampling_schedule$actual_final_burnin,10)
  expect_equal(f$sampling_schedule$retained_draws,20)
  expect_true(all(f$stage_history$feedback_active[f$stage_history$iteration>8]==0))
  expect_true(all(f$matching_history[8:40,]==1))
  expect_equal(f$sampling_schedule$freeze_consensus_fraction,1)
  expect_error(do.call(run_scNPM,c(a,list(freeze_policy="error"))),"strict freeze_policy")
  a$freeze_patience <- 2L
  g <- do.call(run_scNPM,a)
  expect_true(g$sampling_schedule$matching_stable)
  expect_equal(g$sampling_schedule$freeze_reason,"stable")
  expect_null(g$propensity_scaling)
  expect_length(g$global_reference_scaling$center,2L)
})

test_that("timestamped progress records match completed draws without changing them", {
  threads <- getFromNamespace("set_scNPM_threads","scNPM")
  old <- threads(1L);on.exit(threads(old))
  a <- new_freeze_fixture();a$freeze_patience <- 2L
  quiet <- do.call(run_scNPM,a);seed <- .Random.seed
  path <- tempfile(fileext=".csv");on.exit(unlink(path),add=TRUE)
  logged <- do.call(run_scNPM,c(a,list(progress_file=path,progress_every=.00001)))
  expect_identical(logged$Result,quiet$Result)
  expect_identical(.Random.seed,seed)
  z <- utils::read.csv(path)
  expect_equal(z$iteration[1],0)
  expect_equal(tail(z$iteration,1),40)
  expect_equal(tail(z$retained_draws,1),20)
  expect_equal(tail(z$stage,1),"completed")
  expect_true(all(diff(z$iteration)>0))
  expect_true(all(diff(z$elapsed_seconds)>=0))
  expect_true(any(z$stage=="final_burnin"))
  expect_true(any(z$stage=="retained"))
  mid <- z$iteration>0 & z$iteration<a$num_iter
  expect_identical(z$stage[mid],logged$stage_history$stage[z$iteration[mid]])
  expect_true(all(z$stable_required==a$freeze_patience))
  expect_error(do.call(run_scNPM,c(a,list(progress_every=0))),"progress_every")
  expect_error(do.call(run_scNPM,c(a,list(progress_file=tempdir()))),"progress_file")
})
