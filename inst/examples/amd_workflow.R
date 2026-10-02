# Complete retinal example for scNPM 1.0.6.
# This script runs MCMC; set the working directory for its output folder.

library(scNPM)
data("top50_goal_gene_mt", package = "scNPM")
expression <- as.matrix(goal_gene_mt)
stopifnot(identical(as.character(cell$Barcode), colnames(expression)))
stopifnot(all(cell$PR_or_MR %in% c("No", "Yes")))

qc_names <- c(
  "scater_qc.all.total_features_by_counts",
  "scater_qc.all.total_counts",
  "scater_qc.endogenous.pct_counts",
  "scater_qc.feature_control_MT.pct_counts",
  "scater_qc.all.pct_counts_in_top_50_features"
)
cell_t <- cbind(
  treatment = as.integer(cell$PR_or_MR == "Yes"),
  as.matrix(cell[, qc_names, drop = FALSE])
)
colnames(cell_t) <- c("treatment", "detected_genes", "total_counts",
                      "endogenous_pct", "mitochondrial_pct", "top50_pct")
stopifnot(nrow(cell_t) == ncol(expression), all(is.finite(cell_t)))
ppi <- system.file("extdata", "top50_string_interactions_short.tsv",
                   package = "scNPM")
stopifnot(nzchar(ppi))

fit <- run_scNPM(
  goal_gene_mt = expression, cell_t = cell_t, protein_input = ppi, K = 4,
  num_iter = 5000, num_save = 2000,
  inference_mode = "staged", adapt_min = 500, adapt_max = 2000,
  freeze_patience = 20, final_burnin = 1000,
  network_covariates = "all_reference", matching_method = "lexicographic",
  differential_method = "bayes_fdr", fdr_level = 0.05,
  fdr_scope = "frozen", fdr_rule = "delta_threshold",
  group_fdr_level = 0.05, diagnostic_every = 10, ds_seed = 123
)

out_dir <- "scNPM_example_output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_csv <- function(x, name, row_names = FALSE) {
  write.csv(x, file.path(out_dir, name), row.names = row_names)
}
saveRDS(fit, file.path(out_dir, "fit.rds"))
write_csv(fit$matching_result_df, "matching.csv")
write_csv(compute_qf_smd(fit), "network_feature_smd.csv")
write_csv(fit$smd_covariate_results, "qc_covariate_smd.csv", row_names = TRUE)
write_csv(fit$group_fdr$summary, "group_bfdr_summary.csv")
write_csv(fit$group_fdr$table, "group_all.csv")
write_csv(fit$group_fdr$selected, "group_selected.csv")
set.seed(123)
for (network_name in names(fit$group_fdr$networks)) {
  net <- fit$group_fdr$networks[[network_name]]
  write_csv(net$mask, paste0("group_", network_name, "_mask.csv"), row_names = TRUE)
  if (!any(net$mask)) next
  p <- plot_group_network(net)
  if (!is.null(p)) ggplot2::ggsave(
    file.path(out_dir, paste0("group_", network_name, ".pdf")),
    plot = p, width = 7, height = 5
  )
}
write_csv(fit$differential_fdr$table, "differential_all.csv")
write_csv(fit$differential_fdr$selected, "differential_selected.csv")
write_csv(fit$differential_fdr$threshold_curve, "bfdr_threshold_curve.csv")
write_csv(fit$stage_history, "stage_history.csv")
write_csv(data.frame(iteration = seq_len(nrow(fit$matching_history)),
                     fit$matching_history), "matching_history.csv")
write_csv(fit$criterion_trace, "criterion_trace.csv")

maps <- fit$matching_history
mapping_changes <- if (nrow(maps) < 2L) 0L else
  sum(rowSums(maps[-1, , drop = FALSE] !=
                maps[-nrow(maps), , drop = FALSE]) > 0)
write_csv(data.frame(
  runtime_minutes = as.numeric(fit$runtime, units = "mins"),
  frozen_at = fit$sampling_schedule$frozen_at,
  adaptation_iterations = sum(fit$stage_history$stage == "adaptation"),
  mapping_changes = mapping_changes,
  retained_draws = fit$sampling_schedule$retained_draws,
  cutoff = fit$differential_fdr$threshold,
  selected_edges = nrow(fit$differential_fdr$selected),
  estimated_BFDR = fit$differential_fdr$estimated_error_bound,
  selection_status = fit$differential_fdr$status
), "run_summary.csv")

trace <- fit$criterion_trace
if (sum(is.finite(trace$J_aug)) >= 2L) {
  pdf(file.path(out_dir, "criterion_trace.pdf"), width = 7, height = 4)
  plot(trace$iteration, trace$J_aug, type = "l",
       xlab = "Iteration", ylab = "Diagnostic criterion")
  abline(v = fit$sampling_schedule$frozen_at, lty = 2, col = "grey50")
  dev.off()
}
set.seed(123)  # Reproducible plotting layouts
for (pair_name in names(fit$edge_diff_networks)) {
  net <- fit$edge_diff_networks[[pair_name]]
  if (net$solid_network$count + net$dashed_network$count == 0L) next
  p <- plot_diff_network(net$diff_matrix, fit$gene_names,
                         output_mode = "combined")
  if (!is.null(p)) ggplot2::ggsave(
    file.path(out_dir, paste0("differential_", pair_name, ".pdf")),
    plot = p, width = 7, height = 5
  )
}
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
