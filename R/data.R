#' Top-50 Gene Expression Matrix (AMD Example Data)
#'
#' Sparse UMI count matrix (genes x cells) for the 50 target genes used in the
#' AMD (age-related macular degeneration) example analysis. This is a stratified
#' random subset (200 PR control and 200 MR treatment cells) of the full
#' 20091-cell dataset, drawn with a fixed seed. Loaded via
#' \code{data(top50_goal_gene_mt)} together with the companion objects
#' \code{cell}, \code{goal_gene}, \code{mt}, \code{cols_to_PR} and
#' \code{cols_to_MR}.
#'
#' @format A \link[Matrix:dgTMatrix-class]{dgTMatrix} with 50 rows (genes) and
#'   400 columns (cells); entries are raw UMI counts.
#' @source Menon et al. (2019), Nature Communications 10:4902,
#'   \doi{10.1038/s41467-019-12780-8}; GEO \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137537}{GSE137537}.
#'   Normal human retinas, macular versus peripheral regions.
#'   Subset preparation and third-party reuse terms: installed file
#'   \code{extdata/DATA_SOURCES.md}.
#' @docType data
#' @name goal_gene_mt
#' @usage data(top50_goal_gene_mt)
NULL

#' Cell-Level Metadata (AMD Example Data)
#'
#' A data frame with per-cell metadata (quality-control metrics computed with
#' \pkg{scater}, sample provenance and the \code{PR_or_MR} region indicator)
#' for the 400 cells of the AMD example subset (200 PR, 200 MR). Rows are
#' aligned with the columns of \code{goal_gene_mt}. The first example
#' covariate column in \code{\link{run_scNPM}} is typically built from
#' \code{PR_or_MR} (1 = MR, treatment; 0 = PR, control).
#'
#' @format A data frame with 400 rows and 52 columns: cell barcodes
#'   (\code{Barcode}), cell type labels (\code{Labels}), the region indicator
#'   \code{PR_or_MR}, donor identifier (\code{individual}), tissue, spatial
#'   coordinates (\code{x}, \code{y}), and \pkg{scater} quality-control metrics
#'   (columns prefixed \code{scater_qc.}: library size, detected features,
#'   percentage of counts in top features, mitochondrial percentage, etc.).
#' @source Menon et al. (2019), Nature Communications 10:4902,
#'   \doi{10.1038/s41467-019-12780-8}; GEO \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137537}{GSE137537}.
#'   Normal human retinas, macular versus peripheral regions.
#'   Subset preparation and third-party reuse terms: installed file
#'   \code{extdata/DATA_SOURCES.md}.
#' @docType data
#' @name cell
#' @usage data(top50_goal_gene_mt)
NULL

#' Full Gene Expression Matrix (AMD Example Data)
#'
#' Sparse UMI count matrix (all detected genes x cells) from which the top-50
#' target gene matrix \code{goal_gene_mt} was extracted. Provided so users can
#' derive alternative gene subsets. This object covers the same 400 cells as
#' \code{goal_gene_mt}; genes detected in fewer than 3 of these cells were
#' removed.
#'
#' @format A \link[Matrix:dgTMatrix-class]{dgTMatrix} with 12530 rows (genes)
#'   and 400 columns (cells); entries are raw UMI counts.
#' @source Menon et al. (2019), Nature Communications 10:4902,
#'   \doi{10.1038/s41467-019-12780-8}; GEO \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137537}{GSE137537}.
#'   Normal human retinas, macular versus peripheral regions.
#'   Subset preparation and third-party reuse terms: installed file
#'   \code{extdata/DATA_SOURCES.md}.
#' @docType data
#' @name mt
#' @usage data(top50_goal_gene_mt)
NULL

#' Target Gene Names (AMD Example Data)
#'
#' Character vector with the symbols of the 50 target genes (rows of
#' \code{goal_gene_mt}).
#'
#' @format A character vector of length 50 (e.g., \code{"FTH1"}, \code{"SAG"},
#'   \code{"TF"}, \code{"CLU"}, \code{"PDE6G"}, \code{"GLUL"}, ...).
#' @source Selected in the scNPM reference publication.
#' @docType data
#' @name goal_gene
#' @usage data(top50_goal_gene_mt)
NULL

#' Column Indices of PR Cells (AMD Example Data)
#'
#' Integer vector of column indices (in \code{goal_gene_mt}/\code{mt}) of the
#' cells belonging to the PR group, used as the control arm in the example
#' analysis.
#'
#' @format An integer vector of length 200.
#' @source Derived from the \code{PR_or_MR} metadata of the example dataset.
#' @docType data
#' @name cols_to_PR
#' @usage data(top50_goal_gene_mt)
NULL

#' Column Indices of MR Cells (AMD Example Data)
#'
#' Integer vector of column indices (in \code{goal_gene_mt}/\code{mt}) of the
#' cells belonging to the MR group, used as the treatment arm in the example
#' analysis.
#'
#' @format An integer vector of length 200.
#' @source Derived from the \code{PR_or_MR} metadata of the example dataset.
#' @docType data
#' @name cols_to_MR
#' @usage data(top50_goal_gene_mt)
NULL
