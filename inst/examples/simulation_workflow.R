# Small, self-contained scNPM simulation: generate data, fit, and export results.
# Run in a fresh R session. Set OMP_NUM_THREADS before starting R to choose threads.
library(scNPM)
stopifnot(packageVersion("scNPM") >= "1.1.1")
set.seed(123)
G <- 12L                         # genes
K <- 2L                          # clusters per condition
n_per_condition <- 160L          # 320 newly generated cells in total
genes <- sprintf("Gene%02d", seq_len(G))

# Symmetric, strictly diagonally dominant precision matrices are positive definite.
make_precision <- function(edges) {
  omega <- matrix(0, G, G, dimnames = list(genes, genes))
  omega[edges] <- -0.35
  omega <- omega + t(omega)
  diag(omega) <- rowSums(abs(omega)) + 0.5
  stopifnot(min(eigen(omega, symmetric = TRUE, only.values = TRUE)$values) > 0)
  omega
}
edges1 <- cbind(1:5, 2:6)
edges2 <- cbind(7:11, 8:12)
omega_control <- list(make_precision(edges1), make_precision(edges2))
omega_treatment <- list(
  make_precision(rbind(edges1[-1, ], c(1, 6))),
  make_precision(rbind(edges2[-1, ], c(7, 12)))
)
mu <- matrix(1.5, G, K)
mu[1:6, 1] <- 3.5
mu[7:12, 2] <- 3.5

# Latent Gaussian expression, followed by expression-dependent probit dropout.
# Observed values are nonnegative continuous expression, not UMI counts.
generate_condition <- function(precision, prefix) {
  labels <- sample(rep(seq_len(K), length.out = n_per_condition))
  theta <- matrix(NA_real_, G, n_per_condition)
  for (k in seq_len(K)) {
    ix <- which(labels == k)
    sigma <- solve(precision[[k]])
    theta[, ix] <- mu[, k] + t(
      matrix(rnorm(length(ix) * G), ncol = G) %*% chol(sigma)
    )
  }
  x <- exp(theta)
  x[matrix(runif(length(x)), G) < pnorm(1 - theta)] <- 0
  dimnames(x) <- list(genes, paste0(prefix, seq_len(n_per_condition)))
  list(x = x, labels = labels)
}
control <- generate_condition(omega_control, "C")
treatment <- generate_condition(omega_treatment, "T")
expression <- cbind(control$x, treatment$x)
cell_t <- cbind(treatment = rep(0:1, each = n_per_condition),
                measured_covariate = rnorm(2 * n_per_condition))
rownames(cell_t) <- colnames(expression)

# A simulated, imperfect PPI reference; it is not a downloaded biological network.
support <- Reduce(`|`, lapply(c(omega_control, omega_treatment), function(x) x != 0))
ppi_mask <- matrix(runif(G * G), G) < ifelse(support, 0.8, 0.05)
ij <- which(upper.tri(ppi_mask) & ppi_mask, arr.ind = TRUE)
ppi <- data.frame(node1 = genes[ij[, 1]], node2 = genes[ij[, 2]], combined_score = 1)
out_dir <- getOption("scNPM.simulation_output", "scNPM_simulation_output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(expression = expression, cell_t = cell_t, ppi = ppi,
             truth = list(control_precision = omega_control,
                          treatment_precision = omega_treatment,
                          control_labels = control$labels,
                          treatment_labels = treatment$labels,
                          correspondence = data.frame(T = 1:2, C = 1:2))),
        file.path(out_dir, "simulated_inputs_and_truth.rds"))

fit <- run_scNPM(
  goal_gene_mt = expression, cell_t = cell_t, protein_input = ppi, K = K,
  ds_seed = 123, num_iter = 1200, num_save = 400,
  inference_mode = "staged", network_prior = "paired", ppi_reference = "union",
  adapt_min = 100, adapt_max = 400, freeze_patience = 5,
  final_burnin = 400, freeze_policy = "consensus",
  matching_score = "pairwise", network_covariates = "all_reference",
  matching_method = "lexicographic", differential_method = "bayes_fdr",
  fdr_level = 0.05, group_fdr_level = 0.05,
  fdr_scope = "frozen", fdr_rule = "delta_threshold",
  verbose = TRUE, diagnostic_every = 10, progress_every = 30,
  progress_file = file.path(out_dir, "iteration_progress.csv")
)
saveRDS(fit, file.path(out_dir, "fit.rds"))
for (name in c("matching_result_df", "criterion_trace", "stage_history")) {
  write.csv(fit[[name]], file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
}
write.csv(fit$group_fdr$selected, file.path(out_dir, "group_selected.csv"), row.names = FALSE)
write.csv(fit$differential_fdr$selected, file.path(out_dir, "differential_selected.csv"), row.names = FALSE)
write.csv(fit$differential_fdr$table, file.path(out_dir, "differential_all.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
print(fit$matching_result_df)
print(fit$sampling_schedule)
print(fit$group_fdr$summary)
print(head(fit$differential_fdr$selected))
