test_that('R and native progress are quiet or suppressible without changing draws', {
  threads <- getFromNamespace('set_scNPM_threads','scNPM')
  old <- threads(1L); on.exit(threads(old))
  set.seed(92)
  x <- matrix(rpois(4*24,4),4,dimnames=list(paste0('G',1:4),NULL))
  args <- list(goal_gene_mt=x,cell_t=cbind(condition=rep(0:1,each=12)),
    K=1,num_iter=35,num_save=20,adapt_min=3,adapt_max=5,freeze_patience=2,
    final_burnin=5,fdr_min_draws=20,num_step_theta=2,num_step_lam=2,
    epsilon_theta=.01,epsilon_lam=.01,ds_seed=357)
  capture_run <- function(verbose,suppress=FALSE) {
    msgs <- character(); fit <- NULL
    job <- function() do.call(run_scNPM,c(args,list(verbose=verbose)))
    stdout <- capture.output(fit <- withCallingHandlers(
      if(suppress) suppressMessages(job()) else job(),
      message=function(m) {msgs <<- c(msgs,conditionMessage(m)); invokeRestart('muffleMessage')}
    ))
    list(fit=fit,messages=msgs,stdout=stdout,seed=.Random.seed)
  }
  quiet <- capture_run(FALSE)
  loud <- capture_run(TRUE)
  suppressed <- capture_run(TRUE,TRUE)
  expect_length(quiet$stdout,0); expect_length(quiet$messages,0)
  expect_length(loud$stdout,0)
  expect_true(any(grepl('MCMC iteration',loud$messages)))
  expect_true(any(grepl('MCMC Runtime',loud$messages)))
  expect_length(suppressed$stdout,0); expect_length(suppressed$messages,0)
  expect_identical(quiet$fit$Result,loud$fit$Result)
  expect_identical(quiet$fit$Result,suppressed$fit$Result)
  expect_identical(quiet$seed,loud$seed)
  expect_identical(quiet$seed,suppressed$seed)
  expect_error(do.call(run_scNPM,c(args,list(verbose=NA))),'verbose')
  expect_error(do.call(run_scNPM,c(args,list(verbose=1))),'verbose')
})
