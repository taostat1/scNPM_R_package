# Extracted from the production simulation functions; see PROVENANCE.md.

align_labels <- function(pred, truth) {
    acc <- apply(perms3, 1, function(p) sum(p[pred] == truth))
    as.integer(perms3[which.max(acc), ])
}

ari <- function(x, y) {
    tb <- table(x, y)
    choose2 <- function(v) sum(v * (v - 1)/2)
    a <- choose2(tb)
    b <- choose2(rowSums(tb))
    c <- choose2(colSums(tb))
    expected <- b * c/choose2(length(x))
    denom <- (b + c)/2 - expected
    if (denom == 0)
        return(if (all(outer(x, x, "==") == outer(y, y, "=="))) 1 else 0)
    (a - expected)/denom
}

evaluate_fit <- function(f, truth) {
    maps <- lapply(c("T", "C"), function(d) align_labels(f$labels[[d]], truth$cells[[d]]$labels))
    correct <- vapply(1:3, function(k) maps[[2]][f$mapping[k]] == truth$networks$config$control_label_permutation[maps[[1]][k]],
        TRUE)
    rows <- list()
    for (i in 1:2) for (k in 1:3) {
        d <- c("T", "C")[i]
        truek <- maps[[i]][k]
        A <- adjacency(truth$networks$precision[[d]][[truek]])
        rows[[length(rows) + 1L]] <- data.frame(condition = d, inferred_cluster = k, true_cluster = truek, as.list(edge_metrics(f$masks[[d]][[k]],
            f$scores[[d]][[k]], A)))
    }
    graphs <- do.call(rbind, rows)
    summary <- c(matching_accuracy = mean(correct), exact_matching = as.numeric(all(correct)), ARI_T = ari(f$labels$T,
        truth$cells$T$labels), ARI_C = ari(f$labels$C, truth$cells$C$labels), frozen_at = f$frozen_at, colMeans(graphs[,
        c("precision", "recall", "F1", "FDP", "FPR", "AP")]))
    list(summary = summary, graphs = graphs, alignment = maps)
}

average_precision <- function(y, s) {
    stopifnot(length(y) == length(s), all(is.finite(s)))
    if (!sum(y))
        return(NA_real_)
    o <- order(s, decreasing = TRUE)
    ends <- which(!duplicated(s[o], fromLast = TRUE))
    rr <- cumsum(y[o])[ends]/sum(y)
    pp <- cumsum(y[o])[ends]/ends
    sum(diff(c(0, rr)) * pp)
}

metrics <- function(mask, score, truth) {
    mask <- as.logical(mask)
    truth <- as.logical(truth)
    tp <- sum(mask & truth)
    fp <- sum(mask & !truth)
    fn <- sum(!mask & truth)
    tn <- sum(!mask & !truth)
    c(precision = if (tp + fp) tp/(tp + fp) else NA_real_, recall = if (tp + fn) tp/(tp + fn) else NA_real_, F1 = if (2 *
        tp + fp + fn) 2 * tp/(2 * tp + fp + fn) else NA_real_, FDP = fp/max(1, tp + fp), FPR = fp/max(1, fp + tn),
        AP = average_precision(truth, score), TP = tp, FP = fp, FN = fn, TN = tn, selected = tp + fp, true_edges = tp +
            fn)
}

edge_metrics <- function(mask, score, truth) {
    up <- upper.tri(truth)
    metrics(mask[up], score[up], truth[up])
}

evaluate_all <- function(f, truth) {
    e <- evaluate_fit(f, truth)
    d <- f$differential
    d$truth_gain <- d$truth_loss <- FALSE
    for (k in 1:3) {
        rows <- which(d$treatment_cluster == k)
        t <- e$alignment[[1]][k]
        c <- truth$networks$config$control_label_permutation[t]
        a <- adjacency(truth$networks$precision$T[[t]])
        b <- adjacency(truth$networks$precision$C[[c]])
        ix <- cbind(d$gene_i[rows], d$gene_j[rows])
        d$truth_gain[rows] <- (a & !b)[ix]
        d$truth_loss[rows] <- (!a & b)[ix]
    }
    gain <- d$selected & d$direction == "treatment"
    loss <- d$selected & d$direction == "control"
    out <- list()
    add <- function(scope, m, s, y) {
        out[[length(out) + 1L]] <<- data.frame(scope = scope, as.list(metrics(m, s, y)))
    }
    add("differential_directional", c(gain, loss), c(d$p_treatment_only, d$p_control_only), c(d$truth_gain, d$truth_loss))
    add("differential_undirected", d$selected, d$abs_delta_p, d$truth_gain | d$truth_loss)
    add("gain", gain, d$p_treatment_only, d$truth_gain)
    add("loss", loss, d$p_control_only, d$truth_loss)
    xor_gain <- xor_loss <- numeric(nrow(d))
    for (k in 1:3) {
        rows <- which(d$treatment_cluster == k)
        ix <- cbind(d$gene_i[rows], d$gene_j[rows])
        c <- f$mapping[k]
        a <- f$masks$T[[k]][ix]
        b <- f$masks$C[[c]][ix]
        xor_gain[rows] <- a & !b
        xor_loss[rows] <- !a & b
    }
    add("differential_support_xor", c(xor_gain, xor_loss), c(d$p_treatment_only, d$p_control_only), c(d$truth_gain,
        d$truth_loss))
    for (k in 1:3) {
        ix <- d$treatment_cluster == k
        add(paste0("differential_T", e$alignment[[1]][k]), c(gain[ix], loss[ix]), c(d$p_treatment_only[ix], d$p_control_only[ix]),
            c(d$truth_gain[ix], d$truth_loss[ix]))
    }
    mm <- ss <- yy <- list()
    i <- 0L
    for (dd in 1:2) for (k in 1:3) {
        i <- i + 1L
        label <- c("T", "C")[dd]
        a <- adjacency(truth$networks$precision[[label]][[e$alignment[[dd]][k]]])
        up <- upper.tri(a)
        mm[[i]] <- f$masks[[label]][[k]][up]
        ss[[i]] <- f$scores[[label]][[k]][up]
        yy[[i]] <- a[up]
    }
    add("overall_micro", unlist(mm), unlist(ss), unlist(yy))
    metric_names <- names(metrics(unlist(mm), unlist(ss), unlist(yy)))
    out[[length(out) + 1L]] <- data.frame(scope = "overall_macro", as.list(colMeans(e$graphs[, metric_names, drop = FALSE])))
    e$metrics <- do.call(rbind, out)
    e$edge_audit <- d
    stopifnot(sum(d$truth_gain) == 33, sum(d$truth_loss) == 33)
    e
}
