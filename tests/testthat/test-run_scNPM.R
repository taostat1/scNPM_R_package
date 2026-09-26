test_that("run_scNPM The main function runs normally (tested with real data in short iterations).", {
  data(top50_goal_gene_mt)
  tsv_path <- system.file(
    "extdata",
    "top50_string_interactions_short.tsv",
    package = "scNPM"
  )
  tsv_lines <- readLines(tsv_path)
  tsv_lines_processed <- gsub("^\\s*#", "", tsv_lines)
  tsv_lines_processed <- tsv_lines_processed[tsv_lines_processed != ""]
  top50_string_interactions_short <- read.table(
    text = paste(tsv_lines_processed, collapse = "\n"),
    sep = "\t", header = TRUE, stringsAsFactors = FALSE,
    check.names = FALSE
  )

  test_gene_mt <- as.matrix(goal_gene_mt)
  N_cells <- ncol(test_gene_mt)
  cell_t <- matrix(
    c(
      ifelse(cell$PR_or_MR == "Yes", 1, 0),
      scale(cell$scater_qc.all.total_counts),
      scale(cell$scater_qc.all.total_features_by_counts)
    ),
    ncol = 3,
    dimnames = list(NULL, c("treatment", "total_counts", "total_features"))
  )

  # Run the main function
  result <- run_scNPM(inference_mode = "coupled", 
    goal_gene_mt = test_gene_mt,
    cell_t = cell_t,
    K = 4,
    protein_input = top50_string_interactions_short,
    ds_seed = 123,
    num_iter = 10,
    num_save = 5
  )

  expect_type(result, "list")
  expect_true(all(c("edge_diff_networks", "mu_treat_est", "edge_treat_est") %in% names(result)))
})
