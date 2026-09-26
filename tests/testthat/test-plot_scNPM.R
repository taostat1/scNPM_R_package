# Synthetic gene names and edge matrices used to test the plotting logic
# without shipping the 40+ MB posterior result files in data/
make_test_gene_names <- function(n = 10) paste0("Gene", seq_len(n))

make_test_edge_mat <- function(gene_names, p_high = 0.3) {
  G <- length(gene_names)
  m <- matrix(runif(G * G), nrow = G, ncol = G,
              dimnames = list(gene_names, gene_names))
  m[lower.tri(m)] <- t(m)[lower.tri(m)]  # symmetric
  diag(m) <- 0
  # Force some edges above the 0.5 binarization threshold
  high_idx <- which(m > (1 - p_high))
  m[high_idx] <- runif(length(high_idx), 0.6, 1)
  m
}

make_test_diff_mat <- function(gene_names, n_gain = 5, n_loss = 5) {
  G <- length(gene_names)
  m <- matrix(0, nrow = G, ncol = G,
              dimnames = list(gene_names, gene_names))
  ut <- which(upper.tri(m), arr.ind = TRUE)
  set.seed(42)
  gain <- ut[sample(seq_len(nrow(ut)), n_gain), , drop = FALSE]
  m[gain] <- runif(n_gain, 0.6, 0.9)
  loss <- ut[sample(seq_len(nrow(ut)), n_loss), , drop = FALSE]
  m[loss] <- runif(n_loss, -0.9, -0.6)
  m <- m + t(m)  # symmetric
  m
}

# Test 1: plot_group_network functionality
test_that("plot_group_network works with synthetic edge matrices", {
  set.seed(123)
  gene_names <- make_test_gene_names(10)
  treat_edge_mat <- make_test_edge_mat(gene_names)
  ctrl_edge_mat <- make_test_edge_mat(gene_names)

  p_treat <- plot_group_network(
    edge_mat = treat_edge_mat,
    gene_names = gene_names,
    group_type = "treatment"
  )
  expect_true("ggplot" %in% class(p_treat))

  p_ctrl <- plot_group_network(
    edge_mat = ctrl_edge_mat,
    gene_names = gene_names,
    group_type = "control"
  )
  expect_true("ggplot" %in% class(p_ctrl))
})

# Test 2: plot_diff_network functionality
test_that("plot_diff_network works with synthetic differential matrices", {
  gene_names <- make_test_gene_names(10)
  diff_mat <- make_test_diff_mat(gene_names)

  p_combined <- plot_diff_network(
    diff_mat = diff_mat,
    gene_names = gene_names,
    output_mode = "combined",
    solid_threshold = 0.5,
    dashed_threshold = -0.5
  )
  expect_true("ggplot" %in% class(p_combined) || is.null(p_combined))

  p_separate <- plot_diff_network(
    diff_mat = diff_mat,
    gene_names = gene_names,
    output_mode = "separate",
    solid_threshold = 0.5,
    dashed_threshold = -0.5
  )
  expect_type(p_separate, "list")
  expect_true(all(c("solid", "dashed") %in% names(p_separate)))
  if (!is.null(p_separate$solid)) expect_true("ggplot" %in% class(p_separate$solid))
  if (!is.null(p_separate$dashed)) expect_true("ggplot" %in% class(p_separate$dashed))
})

# Test 3: empty matrices return NULL with a warning (no edges)
test_that("plot functions handle empty networks gracefully", {
  gene_names <- make_test_gene_names(5)
  empty_mat <- matrix(0, nrow = 5, ncol = 5,
                      dimnames = list(gene_names, gene_names))

  expect_null(suppressWarnings(plot_group_network(
    edge_mat = empty_mat, gene_names = gene_names, group_type = "treatment"
  )))
  expect_null(suppressWarnings(plot_diff_network(
    diff_mat = empty_mat, gene_names = gene_names, output_mode = "combined"
  )))
})
