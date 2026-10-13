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

test_that("staged sampling freezes, turns feedback off, and retains only final draws", {
  staged_single_thread({
    f <- staged_run(staged_fixture())
    s <- f$sampling_schedule; h <- f$stage_history
    expect_identical(f$inference_mode,"staged")
    expect_identical(f$differential_method,"bayes_fdr")
    expect_equal(s$frozen_at,6)
    expect_gte(s$actual_final_burnin,s$minimum_final_burnin)
    expect_true(all(h$feedback_active[h$iteration>s$frozen_at]==0))
    expect_equal(sum(h$stage=="retained"),40)
    expect_equal(which(h$stage=="retained")[1],s$retained_start)
    expect_equal(unique(as.numeric(f$matching_history[s$frozen_at:70,])),1)
    expect_equal(f$matching_result_df$matched_control_cluster,1)
    expect_equal(dim(f$joint_edge_probabilities$p10),c(4,4,1))
    q <- f$criterion_trace
    expect_equal(q$J_aug,q$log_augmented_likelihood+q$network_prior_kernel-q$matching_loss)
    expect_true(any(diff(q$J_aug[q$stage_code==3])<0))
    expect_true(all(is.finite(q$J_aug)))
    expect_equal(f$differential_fdr$scope,"frozen")
    expect_equal(f$differential_fdr$rule,"delta_threshold")
    expect_equal(f$matching_method,"lexicographic")
    expect_lte(f$differential_fdr$estimated_error_bound,.05)
    expect_false(f$differential_fdr$frequentist_guarantee)
    net <- f$edge_diff_networks[[1]]
    expect_equal(attr(net$diff_matrix,"scNPM_selection")$solid,net$solid_network$mat)
    expect_equal(attr(net$diff_matrix,"scNPM_selection")$dashed,net$dashed_network$mat)
  })
})

test_that("joint probabilities equal within-draw conditional products", {
  staged_single_thread({
    f <- staged_run(staged_fixture())
    prob <- getFromNamespace("calculate_edge_probability","scNPM")
    draws <- f$Result$invcov_treat_iter
    pt <- vapply(seq_along(draws),function(s) prob(draws[[s]][1,2,1],.02,1,.05),numeric(1))
    pc <- vapply(seq_along(draws),function(s) prob(f$Result$invcov_ctrl_iter[[s]][1,2,1],.02,1,.05),numeric(1))
    expect_equal(f$joint_edge_probabilities$p10[1,2,1],mean(pt*(1-pc)),tolerance=1e-12)
    expect_equal(f$joint_edge_probabilities$p01[1,2,1],mean((1-pt)*pc),tolerance=1e-12)
    expect_gt(abs(mean(pt*(1-pc))-mean(pt)*(1-mean(pc))),1e-12)
    z <- f$joint_edge_probabilities
    expect_equal((z$p00+z$p10+z$p01+z$p11)[,,1][upper.tri(diag(4))],rep(1,6))
    expect_equal(sum(z$batch_sizes),40)
    expect_equal(weighted.mean(z$batch_p10[1,1,],z$batch_sizes),z$p10[1,2,1])
  })
})

test_that("K=4 staged inference retains eight references and all sixteen pairs", {
  staged_single_thread({
    a <- staged_fixture(4L)
    # Keep the original global-score regression distinct from new-score tests.
    a$matching_score <- "global"
    set.seed(28)
    a$goal_gene_mt <- do.call(cbind,lapply(1:4,function(k) {
      exp(matrix(3*diag(4)[,k],4,12)+matrix(rnorm(48,sd=.12),4))
    }))
    a$goal_gene_mt <- cbind(a$goal_gene_mt,a$goal_gene_mt)
    rownames(a$goal_gene_mt) <- paste0("G",1:4)
    a$num_iter <- 140L; a$adapt_min <- 20L; a$adapt_max <- 60L
    a$final_burnin <- 40L; a$freeze_patience <- 5L
    f <- staged_run(a)
    expect_equal(ncol(f$network_covariates),8)
    expect_equal(dim(f$joint_edge_probabilities$p10),c(4,4,16))
    expect_equal(nrow(f$differential_fdr$table),4*choose(4,2))
    expect_equal(f$matching_result_df$matched_control_cluster,1:4)
    expect_equal(f$matching_history[a$num_iter,],f$matching_history[f$sampling_schedule$frozen_at,],ignore_attr=TRUE)
  })
})

fdr_fixture <- function(K=1L,G=4L) {
  a <- array(0,c(G,G,K*K))
  list(K=K,G=G,gene_names=paste0("G",seq_len(G)),inference_mode="staged",
    matching_result_df=data.frame(treatment_cluster=seq_len(K),matched_control_cluster=seq_len(K)),
    sampling_schedule=list(feedback_disabled_for_retained=TRUE,retained_draws=500L),
    joint_edge_probabilities=list(p00=a,p10=a,p01=a,p11=a))
}

test_that("posterior directional prefix rule matches a hand calculation", {
  f <- fdr_fixture(); probs <- c(.99,.98,.97,.96,.80,0)
  e <- 0L
  for(i in 1:3) for(j in (i+1):4) {
    e <- e+1L
    f$joint_edge_probabilities$p10[i,j,1] <- f$joint_edge_probabilities$p10[j,i,1] <- probs[e]
    f$joint_edge_probabilities$p00[i,j,1] <- f$joint_edge_probabilities$p00[j,i,1] <- 1-probs[e]
  }
  x <- select_differential_edges(f)
  expect_equal(nrow(x$selected),4)
  expect_equal(x$estimated_error_bound,.025)
  expect_equal(x$selected$error_bound,c(.01,.02,.03,.04))
  expect_equal(x$table$prefix_mean_error[5],.06)
  f$sampling_schedule$retained_draws <- 50L
  y <- select_differential_edges(f)
  expect_equal(y$status,"insufficient_draws")
  expect_equal(nrow(y$selected),0)
  f$joint_edge_probabilities$p10[1,2,1] <- 2
  expect_error(select_differential_edges(f),"Invalid joint")
})

test_that("all-matchings bound protects against an alternative correspondence", {
  f <- fdr_fixture(2L,2L)
  # Posterior mass .98 on (T1,T2,C1,C2)=(1,0,0,1), .02 on its complement.
  ps <- rbind(c(0,.98,.02,0),c(.02,0,0,.98),c(.98,0,0,.02),c(0,.02,.98,0))
  for(h in 1:4) for(j in seq_along(c("p00","p10","p01","p11"))) {
    nm <- c("p00","p10","p01","p11")[j]
    f$joint_edge_probabilities[[nm]][1,2,h] <- f$joint_edge_probabilities[[nm]][2,1,h] <- ps[h,j]
  }
  a <- select_differential_edges(f,scope="frozen")
  b <- select_differential_edges(f,scope="all_matchings")
  expect_equal(nrow(a$selected),2)
  expect_equal(a$estimated_error_bound,.02)
  expect_equal(nrow(b$selected),0)
  expect_equal(b$table$error_bound,c(1,1))
  expect_true(all(b$table$error_bound>=b$table$frozen_error))
  f$inference_mode <- "coupled"
  expect_error(select_differential_edges(f),"staged final sampling")
})

test_that("label assignment is globally minimal and deterministic", {
  assign <- getFromNamespace("staged_assignment","scNPM")
  expect_equal(as.numeric(assign(matrix(c(1,2,.1,3.1),2,byrow=TRUE))),c(2,1))
  expect_equal(as.numeric(assign(matrix(0,4,4))),1:4)
  expect_error(assign(matrix(NA_real_,2,2)),"finite square")
})

test_that("invalid stage budgets and opt-in strict stability cannot produce discoveries", {
  args <- staged_fixture()
  args$final_burnin <- 29
  expect_error(staged_run(args),"final_burnin")
  args <- staged_fixture();args$freeze_patience <- 21
  args$freeze_policy <- "error"
  expect_error(staged_run(args),"Matching did not stabilize")
  args <- staged_fixture();args$inference_mode <- "coupled";args$differential_method <- "bayes_fdr"
  expect_error(staged_run(args),"requires staged")
  staged_single_thread({
    a <- staged_fixture(4L); a$matching_method <- "greedy"
    a$matching_score <- "global"; a$freeze_policy <- "error"
    expect_error(staged_run(a),"Matching did not stabilize")
  })
})


test_that("differential plots preserve Bayesian decisions instead of amplitude thresholds", {
  gn <- paste0("G",1:3)
  dm <- matrix(0,3,3,dimnames=list(gn,gn))
  dm[1,2] <- dm[2,1] <- .99 # strong contrast rejected by the error rule
  dm[1,3] <- dm[3,1] <- .10 # weak contrast selected in a large low-error family
  dm[2,3] <- dm[3,2] <- -.90
  solid <- dashed <- 0*dm
  solid[1,3] <- solid[3,1] <- 1
  dashed[2,3] <- dashed[3,2] <- 1
  attr(dm,"scNPM_selection") <- list(solid=solid,dashed=dashed)
  p <- plot_diff_network(dm,gn)
  expect_equal(nrow(p$layers[[1]]$data),1)
  expect_equal(nrow(p$layers[[2]]$data),1)
  expect_equal(sort(c(p$layers[[1]]$data$from,p$layers[[1]]$data$to)),c("G1","G3"))
  attr(dm,"scNPM_selection")$dashed[] <- 0
  p <- plot_diff_network(dm,gn)
  expect_equal(nrow(p$layers[[2]]$data),0)
  expect_no_error(ggplot2::ggplot_build(p))
  attr(dm,"scNPM_selection")$solid <- matrix(1,2,2)
  expect_error(plot_diff_network(dm,gn),"Invalid scNPM_selection")
})
