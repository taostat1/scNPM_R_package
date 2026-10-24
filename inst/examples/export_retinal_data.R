# Export the packaged public-data subset. No download or MCMC is needed.
# Source and reuse information: inst/extdata/DATA_SOURCES.md.
library(scNPM)
data("top50_goal_gene_mt", package = "scNPM")
qc <- c("scater_qc.all.total_features_by_counts", "scater_qc.all.total_counts",
        "scater_qc.endogenous.pct_counts", "scater_qc.feature_control_MT.pct_counts",
        "scater_qc.all.pct_counts_in_top_50_features")
counts <- as.matrix(goal_gene_mt)
cells <- data.frame(cell_id = as.character(cell$Barcode),
  region = ifelse(cell$PR_or_MR == "Yes", "macular", "peripheral"),
  treatment = as.integer(cell$PR_or_MR == "Yes"), donor = cell$individual,
  sample = cell$tissue, cell[, qc], check.names = FALSE)
names(cells)[6:10] <- c("detected_genes", "total_counts", "endogenous_pct",
                      "mitochondrial_pct", "top50_pct")
stopifnot(identical(cells$cell_id, colnames(counts)))
out_dir <- "retina_GSE137537_export"
dir.create(out_dir, showWarnings = FALSE)
con <- gzfile(file.path(out_dir, "counts.csv.gz"), "wt")
write.csv(counts, con)
close(con)
write.csv(cells, file.path(out_dir, "cells.csv"), row.names = FALSE)
write.csv(data.frame(gene = rownames(counts)), file.path(out_dir, "genes.csv"), row.names = FALSE)
