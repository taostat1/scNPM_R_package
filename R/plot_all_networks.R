#' Core Network Plotting Function
#'
#' Generates publication-quality gene regulatory network plots from adjacency matrices.
#'
#' @param adj_matrix Square adjacency matrix (rows/columns = gene names).
#' @param mode Network type: "undirected" (default) or "directed".
#' @param layout_fun Igraph layout function (default = igraph::layout_nicely; currently informational).
#' @param node_id_col Column name for node IDs (default = "gene").
#' @param size_var Node attribute for size mapping (default = "degree").
#' @param color_var Node attribute for color mapping (default = "degree").
#' @param size_range Numeric vector (length 2) for node size range (default = c(1,5)).
#' @param color_palette Color gradient for nodes (default = c("steelblue", "white", "orangered")).
#' @param label_size Font size for node labels (default = 2).
#' @param edge_color Edge color (default = "#ADD8E6").
#' @param edge_alpha Edge transparency (default = 0.6).
#' @param edge_size Edge thickness (default = 0.3).
#' @param curvature Edge curvature (0 = straight, default = 0).
#' @param edge_linetype Edge line type ("solid" (default) or "dashed").
#' @param legend_position Legend position ("right" (default), "bottom", "none").
#'
#' @return ggplot object (network plot) or NULL if no edges exist.
#' @importFrom igraph graph_from_adjacency_matrix degree subgraph vertex_attr as_edgelist vcount ecount V layout_nicely
#' @importFrom ggplot2 ggplot geom_curve geom_point scale_color_gradientn scale_size_continuous
#' @importFrom ggplot2 theme_minimal theme element_blank element_rect unit guide_colorbar
#' @importFrom ggplot2 guide_legend element_text
#' @importFrom ggrepel geom_text_repel
#' @importFrom scales rescale number_format
#' @keywords internal  # Mark as internal (not exported to users)
plot_network_core <- function(
    adj_matrix,
    mode = "undirected",
    layout_fun = igraph::layout_nicely,
    node_id_col = "gene",
    size_var = "degree",
    color_var = "degree",
    size_range = c(1, 5),
    color_palette = c("steelblue", "white", "orangered"),
    label_size = 2,
    edge_color = "#ADD8E6",
    edge_alpha = 0.6,
    edge_size = 0.3,
    curvature = 0,
    edge_linetype = "solid",
    legend_position = "right"
) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required! Install with: install.packages('igraph')")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required! Install with: install.packages('ggplot2')")
  }
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required! Install with: install.packages('ggrepel')")
  }

  if (!is.matrix(adj_matrix) || nrow(adj_matrix) != ncol(adj_matrix)) {
    stop("adj_matrix must be a square matrix!")
  }
  if (is.null(rownames(adj_matrix)) || is.null(colnames(adj_matrix))) {
    stop("adj_matrix must have row/column names!")
  }

  g <- igraph::graph_from_adjacency_matrix(
    adjmatrix = adj_matrix,
    mode = mode,
    weighted = TRUE,
    add.rownames = node_id_col
  )
  igraph::V(g)$degree <- igraph::degree(g)
  g <- igraph::subgraph(g, v = igraph::V(g)[igraph::V(g)$degree > 0])

  if (igraph::vcount(g) == 0 || igraph::ecount(g) == 0) {
    warning("No edges found in adjacency matrix - returning NULL")
    return(NULL)
  }

  layout <- igraph::layout_nicely(g)
  igraph::V(g)$xcoord <- layout[, 1]
  igraph::V(g)$ycoord <- layout[, 2]

  nodes_df <- as.data.frame(igraph::vertex_attr(g), stringsAsFactors = FALSE)
  nodes_df$size_norm <- scales::rescale(nodes_df[[size_var]], size_range)
  nodes_df$color_norm <- scales::rescale(nodes_df[[color_var]], c(0, 1))

  edges_mat <- igraph::as_edgelist(g, names = TRUE)
  edges_df <- data.frame(
    from = edges_mat[, 1],
    to = edges_mat[, 2],
    stringsAsFactors = FALSE
  )
  edges_df <- merge(edges_df, nodes_df[, c(node_id_col, "xcoord", "ycoord")],
                    by.x = "from", by.y = node_id_col, all.x = TRUE)
  edges_df <- merge(edges_df, nodes_df[, c(node_id_col, "xcoord", "ycoord")],
                    by.x = "to", by.y = node_id_col, all.x = TRUE,
                    suffixes = c(".from", ".to"))

  p <- ggplot2::ggplot() +
    ggplot2::geom_curve(
      data = edges_df,
      ggplot2::aes(x = xcoord.from, y = ycoord.from, xend = xcoord.to, yend = ycoord.to),
      color = edge_color,
      alpha = edge_alpha,
      linewidth  = edge_size,
      curvature = curvature,
      linetype = edge_linetype
    ) +
    ggplot2::geom_point(
      data = nodes_df,
      ggplot2::aes(x = xcoord, y = ycoord, size = size_norm, color = color_norm),
      alpha = 0.9
    ) +
    ggrepel::geom_text_repel(
      data = nodes_df,
      ggplot2::aes(x = xcoord, y = ycoord, label = .data[[node_id_col]]),
      size = label_size,
      color = "black",
      alpha = 0.95,
      box.padding = 0.1,
      force = 0.01
    ) +
    ggplot2::scale_color_gradientn(
      colors = color_palette,
      name = color_var,
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = scales::number_format(accuracy = 0.1),
      guide = ggplot2::guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        title.theme = ggplot2::element_text(face = "bold", size = 10),
        barwidth = 1.5,
        barheight = 15
      )
    ) +
    ggplot2::scale_size_continuous(
      name = size_var,
      range = size_range,
      guide = ggplot2::guide_legend(
        title.position = "top",
        title.hjust = 0.5,
        title.theme = ggplot2::element_text(face = "bold", size = 10),
        keyheight = ggplot2::unit(0.8, "cm"),
        keywidth = ggplot2::unit(0.8, "cm")
      )
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      legend.position = legend_position,
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.key = ggplot2::element_rect(fill = "white"),
      legend.spacing.y = ggplot2::unit(0.5, "cm")
    )

  return(p)
}

#' Plot Group-Specific Gene Regulatory Network (No File Saving)
#'
#' Plots treatment/control group gene regulatory networks from 3D edge estimate arrays
#' (output by run_scNPM) without saving files.
#'
#' @param edge_mat 2D edge matrix (rows/columns = gene names) extracted from 3D estimate array.
#' @param gene_names Vector of gene names (matches row/column names of edge_mat).
#' @param group_type Group type: "treatment" (default) or "control".
#' @param edge_color Edge color for network edges (default = "#ADD8E6", the original scNPM light blue; set e.g. "#E74C3C" for treatment-red or "#3498DB" for control-blue).
#' @param ... Additional parameters passed to plot_network_core().
#'
#' @return ggplot object (network plot) or NULL if no edges exist.
#' @importFrom igraph graph_from_adjacency_matrix
#' @export
#'
#' @examples
#' # Self-contained example with a small random edge matrix (runs in < 5 sec)
#' set.seed(1)
#' gn <- paste0("Gene_", 1:10)
#' em <- matrix(runif(100), 10, 10, dimnames = list(gn, gn))
#' em <- (em + t(em)) / 2
#'
#' # Treatment group network (default light-blue edges)
#' p_treat <- plot_group_network(
#'   edge_mat = em,
#'   gene_names = gn,
#'   group_type = "treatment",
#'   label_size = 2.5,
#'   size_range = c(1.5, 5),
#'   edge_alpha = 0.6,
#'   legend_position = "bottom"
#' )
#'
#' # Custom edge color (e.g., treatment-red)
#' p_treat_red <- plot_group_network(
#'   edge_mat = em, gene_names = gn, group_type = "treatment",
#'   edge_color = "#E74C3C"
#' )
plot_group_network <- function(
    edge_mat,
    gene_names,
    group_type = c("treatment", "control"),
    edge_color = "#ADD8E6",
    ...
) {
  # Validate inputs
  group_type <- match.arg(group_type)
  rownames(edge_mat) <- gene_names
  colnames(edge_mat) <- gene_names

  # Binarize edge matrix (core threshold: weight >= 0.5)
  edge_mat_bin <- edge_mat
  edge_mat_bin[edge_mat_bin >= 0.5] <- 1
  edge_mat_bin[edge_mat_bin < 0.5] <- 0

  # Generate core network plot
  p <- plot_network_core(
    adj_matrix = edge_mat_bin,
    edge_color = edge_color,
    ...
  )

  # Return early if no valid edges
  if (is.null(p)) {
    warning(paste("No edges found for", group_type, "group - skipping plot"))
    return(NULL)
  }

  return(p)
}

#' Plot Differential Gene Regulatory Network (Combined/Separate Mode)
#'
#' Plots differential networks (treatment - control) with flexible output:
#' - Combined mode: solid/dashed edges in a single plot (default)
#' - Separate mode: solid/dashed edges as two independent plots
#'
#' @param diff_mat Differential matrix (treatment - control, rows/columns = gene names).
#'   Matrices returned by Bayesian-mode run_scNPM carry a scNPM_selection
#'   attribute; its selected edge masks are used directly, without re-thresholding.
#' @param gene_names Vector of gene names (matches diff_mat row/column names).
#' @param solid_threshold Minimum positive difference for unannotated matrices
#'   (default = 0.5); ignored for matrices carrying Bayesian selection masks.
#' @param dashed_threshold Maximum negative difference for unannotated matrices
#'   (default = -0.5); ignored for matrices carrying Bayesian selection masks.
#' @param solid_color Color for solid edges (treatment-enhanced, default = "#ADD8E6", the original scNPM light blue; set e.g. "#2ECC71" for green).
#' @param dashed_color Color for dashed edges (treatment-reduced, default = "#ADD8E6", the original scNPM light blue; set e.g. "#F39C12" for orange).
#' @param output_mode Output type: "combined" (single plot, default) or "separate" (two plots).
#' @param ... Additional parameters passed to plot_network_core() (separate mode) or internal plotting (combined mode).
#'
#' @return
#' - If output_mode = "combined": ggplot object (single plot with solid/dashed edges) or NULL (no edges)
#' - If output_mode = "separate": list of ggplot objects (solid/dashed) or NULL for empty edges
#' @importFrom igraph graph_from_adjacency_matrix degree subgraph vertex_attr as_edgelist vcount ecount V layout_nicely
#' @importFrom ggplot2 ggplot geom_curve geom_point scale_color_gradientn scale_size_continuous
#' @importFrom ggplot2 theme_minimal theme element_blank element_rect unit guide_colorbar
#' @importFrom ggplot2 guide_legend element_text aes
#' @importFrom ggrepel geom_text_repel
#' @importFrom scales rescale number_format
#' @export
#'
#' @examples
#' # Self-contained example with a small random differential matrix (< 5 sec)
#' set.seed(1)
#' gn <- paste0("Gene_", 1:10)
#' dm <- matrix(rnorm(100), 10, 10, dimnames = list(gn, gn))
#' dm <- (dm + t(dm)) / 2
#' diag(dm) <- 0
#'
#' # Mode 1 - Combined (solid+dashed in one plot, default colors: light blue)
#' p_combined <- plot_diff_network(
#'   diff_mat = dm,
#'   gene_names = gn,
#'   output_mode = "combined"
#' )
#'
#' # Mode 2 - Separate (solid/dashed as two plots; custom colors optional)
#' p_separate <- plot_diff_network(
#'   diff_mat = dm,
#'   gene_names = gn,
#'   output_mode = "separate",
#'   solid_color = "#2ECC71",
#'   dashed_color = "#F39C12",
#'   label_size = 3
#' )
plot_diff_network <- function(
    diff_mat,
    gene_names,
    solid_threshold = 0.5,
    dashed_threshold = -0.5,
    solid_color = "#ADD8E6",
    dashed_color = "#ADD8E6",
    output_mode = c("combined", "separate"),
    ...
) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required! Install with: install.packages('igraph')")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required! Install with: install.packages('ggplot2')")
  }
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Package 'ggrepel' is required! Install with: install.packages('ggrepel')")
  }

  output_mode <- match.arg(output_mode)
  rownames(diff_mat) <- gene_names
  colnames(diff_mat) <- gene_names

  solid_mat <- matrix(0, nrow = nrow(diff_mat), ncol = ncol(diff_mat),
                      dimnames = list(gene_names, gene_names))
  dashed_mat <- matrix(0, nrow = nrow(diff_mat), ncol = ncol(diff_mat),
                       dimnames = list(gene_names, gene_names))
  selected <- attr(diff_mat,"scNPM_selection")
  if (!is.null(selected)) {
    valid_mask <- function(x) is.matrix(x) && identical(dim(x),dim(diff_mat)) &&
      all(is.finite(x)) && all(x %in% c(0,1)) && all(diag(x)==0) &&
      isTRUE(all.equal(unname(x),unname(t(x))))
    if (!is.list(selected) || !valid_mask(selected$solid) ||
        !valid_mask(selected$dashed) || any(selected$solid*selected$dashed != 0))
      stop("Invalid scNPM_selection edge masks; retain the original matrix dimensions.")
    solid_mat[] <- selected$solid
    dashed_mat[] <- selected$dashed
  } else {
    solid_mat[diff_mat > solid_threshold] <- 1
    dashed_mat[diff_mat < dashed_threshold] <- 1
  }

  if (output_mode == "combined") {
    combined_mat <- solid_mat + dashed_mat
    if (sum(combined_mat) == 0) {
      warning(if (!is.null(selected)) "No selected differential edges - returning NULL" else
        paste("No edges found (solid >", solid_threshold, "or dashed <", dashed_threshold, ") - returning NULL"))
      return(NULL)
    }

    g <- igraph::graph_from_adjacency_matrix(combined_mat, mode = "undirected", weighted = TRUE, add.rownames = "gene")
    igraph::V(g)$degree <- igraph::degree(g)
    g <- igraph::subgraph(g, v = igraph::V(g)[igraph::V(g)$degree > 0])

    if (igraph::vcount(g) == 0 || igraph::ecount(g) == 0) {
      warning("No valid edges for differential network - returning NULL")
      return(NULL)
    }

    layout <- igraph::layout_nicely(g)
    igraph::V(g)$xcoord <- layout[, 1]
    igraph::V(g)$ycoord <- layout[, 2]
    nodes_df <- as.data.frame(igraph::vertex_attr(g), stringsAsFactors = FALSE)
    nodes_df$size_norm <- scales::rescale(nodes_df$degree, c(1, 5))
    nodes_df$color_norm <- scales::rescale(nodes_df$degree, c(0, 1))

    edges_solid_df <- if (sum(solid_mat) > 0) {
      g_solid <- igraph::graph_from_adjacency_matrix(solid_mat, mode = "undirected", add.rownames = "gene")
      edges_solid_mat <- igraph::as_edgelist(g_solid, names = TRUE)
      df <- data.frame(from = edges_solid_mat[,1], to = edges_solid_mat[,2], stringsAsFactors = FALSE)
      merge(merge(df, nodes_df[,c("gene","xcoord","ycoord")], by.x="from", by.y="gene"),
            nodes_df[,c("gene","xcoord","ycoord")], by.x="to", by.y="gene", suffixes = c(".from", ".to"))
    } else {
      data.frame(xcoord.from=numeric(), ycoord.from=numeric(),
                 xcoord.to=numeric(), ycoord.to=numeric())
    }

    edges_dashed_df <- if (sum(dashed_mat) > 0) {
      g_dashed <- igraph::graph_from_adjacency_matrix(dashed_mat, mode = "undirected", add.rownames = "gene")
      edges_dashed_mat <- igraph::as_edgelist(g_dashed, names = TRUE)
      df <- data.frame(from = edges_dashed_mat[,1], to = edges_dashed_mat[,2], stringsAsFactors = FALSE)
      merge(merge(df, nodes_df[,c("gene","xcoord","ycoord")], by.x="from", by.y="gene"),
            nodes_df[,c("gene","xcoord","ycoord")], by.x="to", by.y="gene", suffixes = c(".from", ".to"))
    } else {
      data.frame(xcoord.from=numeric(), ycoord.from=numeric(),
                 xcoord.to=numeric(), ycoord.to=numeric())
    }

    label_size <- if ("label_size" %in% names(list(...))) list(...)$label_size else 2
    p <- ggplot2::ggplot() +

      ggplot2::geom_curve(
        data = edges_solid_df,
        ggplot2::aes(x = xcoord.from, y = ycoord.from, xend = xcoord.to, yend = ycoord.to),
        color = solid_color, alpha = 0.7, linewidth  = 0.3, curvature = 0, linetype = "solid"
      ) +

      ggplot2::geom_curve(
        data = edges_dashed_df,
        ggplot2::aes(x = xcoord.from, y = ycoord.from, xend = xcoord.to, yend = ycoord.to),
        color = dashed_color, alpha = 0.7, linewidth  = 0.3, curvature = 0, linetype = "dashed"
      ) +

      ggplot2::geom_point(data = nodes_df, ggplot2::aes(x = xcoord, y = ycoord, size = size_norm, color = color_norm), alpha = 0.9) +
      ggrepel::geom_text_repel(data = nodes_df, ggplot2::aes(x = xcoord, y = ycoord, label = gene),
                               size = label_size, color = "black", alpha = 0.95, box.padding = 0.1, force = 0.01) +

      ggplot2::scale_color_gradientn(colors = c("steelblue", "white", "orangered"), name = "Degree",
                                     limits = c(0,1), breaks = seq(0,1,0.2), labels = scales::number_format(accuracy = 0.1),
                                     guide = ggplot2::guide_colorbar(title.position = "top", title.hjust = 0.5,
                                                                     title.theme = ggplot2::element_text(face = "bold", size = 10),
                                                                     barwidth = 1.5, barheight = 15)) +
      ggplot2::scale_size_continuous(name = "Degree", range = c(1,5),
                                     guide = ggplot2::guide_legend(title.position = "top", title.hjust = 0.5,
                                                                   title.theme = ggplot2::element_text(face = "bold", size = 10),
                                                                   keyheight = ggplot2::unit(0.8, "cm"), keywidth = ggplot2::unit(0.8, "cm"))) +
      ggplot2::theme_minimal() +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
                     legend.position = "right", plot.background = ggplot2::element_rect(fill = "white", color = NA),
                     legend.key = ggplot2::element_rect(fill = "white"), legend.spacing.y = ggplot2::unit(0.5, "cm"))

    return(p)
  }

  else if (output_mode == "separate") {
    plots <- list()
    if (sum(solid_mat) > 0) {
      plots$solid <- plot_network_core(
        adj_matrix = solid_mat,
        edge_color = solid_color,
        edge_linetype = "solid",
        ...
      )
    } else {
      warning(if (!is.null(selected)) "No selected solid edges - solid plot is NULL" else
        paste("No solid edges (difference >", solid_threshold, ") - solid plot is NULL"))
      plots$solid <- NULL
    }
    if (sum(dashed_mat) > 0) {
      plots$dashed <- plot_network_core(
        adj_matrix = dashed_mat,
        edge_color = dashed_color,
        edge_linetype = "dashed",
        ...
      )
    } else {
      warning(if (!is.null(selected)) "No selected dashed edges - dashed plot is NULL" else
        paste("No dashed edges (difference <", dashed_threshold, ") - dashed plot is NULL"))
      plots$dashed <- NULL
    }
    return(plots)
  }
}
