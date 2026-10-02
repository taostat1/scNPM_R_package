group_fixture <- function(pt = c(.99,.98,.97,.96,.8,0), pc = rep(0,6), K = 1L) {
  G <- 4L; a <- array(0,c(G,G,K*K)); joint <- list(p00=a,p10=a,p01=a,p11=a)
  ij <- which(upper.tri(matrix(0,G,G)),arr.ind=TRUE)
  pt <- matrix(pt,6,K); pc <- matrix(pc,6,K)
  for(k in seq_len(K)) for(c in seq_len(K)) {
    h <- (k-1L)*K+c
    vals <- list(p00=(1-pt[,k])*(1-pc[,c]),p10=pt[,k]*(1-pc[,c]),
                 p01=(1-pt[,k])*pc[,c],p11=pt[,k]*pc[,c])
    for(nm in names(vals)) {
      m <- matrix(0,G,G); m[ij] <- vals[[nm]]; m <- m+t(m); joint[[nm]][,,h] <- m
    }
  }
  list(K=K,G=G,gene_names=paste0('G',1:G),inference_mode='staged',
       matching_result_df=data.frame(treatment_cluster=seq_len(K),matched_control_cluster=rev(seq_len(K))),
       sampling_schedule=list(feedback_disabled_for_retained=TRUE,retained_draws=500L),
       joint_edge_probabilities=joint)
}

test_that('group BFDR matches hand calculation and inclusive threshold', {
  f <- group_fixture(); g <- select_group_edges(f)
  x <- g$networks$treatment_1
  expect_equal(nrow(x$selected),4)
  expect_equal(x$cutoff,.96)
  expect_equal(x$estimated_bfdr,.025)
  expect_identical(x$table$selected,x$table$probability>=x$cutoff)
  expect_equal(sum(x$mask),8)
  expect_true(all(diag(x$mask)==0)); expect_equal(x$mask,t(x$mask))
  expect_equal(g$masks$treatment[,,1],x$mask)
  expect_equal(g$summary$selected_edges,c(4,0))
  expect_identical(g$networks$control_1$cutoff,Inf)
  expect_identical(g$networks$control_1$estimated_bfdr,0)
  expect_false(g$frequentist_guarantee)
})

test_that('ties stay together, full and empty sets work, and draw guard is preserved', {
  g <- select_group_edges(group_fixture(c(.99,.9,.9,0,0,0)))
  x <- g$networks$treatment_1
  expect_equal(nrow(x$selected),1) # splitting the .9 pair would incorrectly admit two
  expect_equal(x$threshold_curve$discoveries,c(1,3,6))
  expect_equal(sum(x$threshold_curve$chosen),1)
  a <- select_group_edges(group_fixture(rep(.99,6)))
  expect_equal(a$summary$selected_edges,c(6,0))
  f <- group_fixture(); f$sampling_schedule$retained_draws <- 100L
  z <- select_group_edges(f)
  expect_equal(z$status,'insufficient_draws'); expect_equal(nrow(z$selected),0)
  expect_true(any(z$networks$treatment_1$threshold_curve$feasible))
  expect_false(any(z$networks$treatment_1$threshold_curve$chosen))
  expect_true(all(is.infinite(z$summary$cutoff)))
})

test_that('independent threshold enumeration agrees and each network has its own family', {
  set.seed(881)
  for (i in 1:25) {
    p <- sample(seq(0,1,by=.01),6,replace=TRUE); alpha <- .073
    x <- select_group_edges(group_fixture(p),fdr_level=alpha)$networks$treatment_1
    cuts <- sort(unique(p),decreasing=TRUE)
    counts <- vapply(cuts,function(cut) {
      s <- p>=cut; if(mean(1-p[s])<=alpha) sum(s) else 0L
    },integer(1))
    expect_equal(nrow(x$selected),max(counts))
    expect_lte(x$estimated_bfdr,alpha)
  }
  f <- group_fixture(c(rep(.99,6),rep(.8,6)),c(rep(.4,6),rep(.98,6)),K=2L)
  g <- select_group_edges(f)
  expect_equal(g$summary$selected_edges,c(6,0,0,6))
  expect_equal(nrow(g$table),24)
  expect_equal(g$networks$treatment_1$table$probability,rep(.99,6))
  expect_equal(g$networks$control_2$table$probability,rep(.98,6))
  # A different bijection reconstructs the same marginals.
  f$matching_result_df$matched_control_cluster <- 1:2
  z <- select_group_edges(f)
  expect_equal(z$masks,g$masks)
  expect_equal(z$probabilities,g$probabilities,tolerance=1e-12)
})

test_that('group selection rejects invalid posterior inputs', {
  f <- group_fixture()
  expect_error(select_group_edges(f,fdr_level=1),'Invalid')
  expect_error(select_group_edges(f,min_draws=19),'Invalid')
  f$inference_mode <- 'coupled'
  expect_error(select_group_edges(f),'staged final sampling')
  f <- group_fixture(); f$joint_edge_probabilities$p00[1,2,1] <- NA
  expect_error(select_group_edges(f),'Invalid joint')
  f <- group_fixture(); f$joint_edge_probabilities$p10[1,2,1] <- .8
  expect_error(select_group_edges(f),'symmetric')
  f <- group_fixture(); f$joint_edge_probabilities$p00[1,2,1] <- .5
  f$joint_edge_probabilities$p00[2,1,1] <- .5
  expect_error(select_group_edges(f),'sum to one')
  f <- group_fixture(K=2L)
  f$joint_edge_probabilities$p10[1,2,2] <- f$joint_edge_probabilities$p10[2,1,2] <- .8
  f$joint_edge_probabilities$p00[1,2,2] <- f$joint_edge_probabilities$p00[2,1,2] <- .2
  expect_error(select_group_edges(f),'Inconsistent marginal')
})

test_that('group plotting uses saved masks and never silently applies .5', {
  g <- select_group_edges(group_fixture(c(.99,.98,.97,.96,.8,0)))
  net <- g$networks$treatment_1
  p <- plot_group_network(net)
  expect_s3_class(p,'ggplot')
  expect_equal(nrow(p$layers[[1]]$data),4)
  expect_s3_class(plot_group_network(net$probability_matrix),'ggplot')
  pm <- net$probability_matrix; attr(pm,'scNPM_group_selection') <- NULL
  expect_error(plot_group_network(pm),'saved group selection')
  old <- plot_group_network(pm,selection='frequency')
  expect_equal(nrow(old$layers[[1]]$data),5)
  expect_error(plot_group_network(net,gene_names=rev(rownames(pm))),'Gene order')
  expect_null(suppressWarnings(plot_group_network(g$networks$control_1)))
  bad <- net$probability_matrix
  attr(bad,'scNPM_group_selection')$mask[1,1] <- TRUE
  expect_error(plot_group_network(bad),'Invalid saved')
})

test_that('staged wrapper returns group decisions from the actual conditional draws', {
  # Tiny baseline run checks the integration, not convergence of scientific data.
  set <- getFromNamespace('set_scNPM_threads','scNPM'); old <- set(1L); on.exit(set(old))
  set.seed(2026); x <- exp(matrix(rnorm(4*24,sd=.2),4))
  rownames(x) <- paste0('G',1:4)
  f <- suppressMessages(run_scNPM(x,cbind(condition=rep(0:1,each=12)),K=1,
        num_iter=70,num_save=40,adapt_min=6,adapt_max=20,freeze_patience=4,
        final_burnin=10,fdr_min_draws=20,num_step_theta=2,num_step_lam=2,
        ds_seed=19,group_fdr_level=.075))
  expect_equal(f$group_fdr$target,.075)
  expect_equal(f$group_fdr$summary$target,rep(.075,2))
  prob <- getFromNamespace('calculate_edge_probability','scNPM')
  pt <- vapply(f$Result$invcov_treat_iter,function(a) prob(a[1,2,1],.02,1,.05),numeric(1))
  pc <- vapply(f$Result$invcov_ctrl_iter,function(a) prob(a[1,2,1],.02,1,.05),numeric(1))
  expect_equal(f$group_fdr$probabilities$treatment[1,2,1],mean(pt),tolerance=1e-12)
  expect_equal(f$group_fdr$probabilities$control[1,2,1],mean(pc),tolerance=1e-12)
  expect_equal(f$group_fdr,select_group_edges(f,fdr_level=.075,min_draws=20))
  expect_equal(f$provenance$settings$group_selection,'bayes_fdr')
  expect_equal(f$provenance$settings$group_fdr_level,.075)
})
