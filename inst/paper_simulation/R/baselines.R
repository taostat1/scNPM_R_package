# Extracted from the production simulation functions; see PROVENANCE.md.

edge_table <- function(A, weight = NULL) {
    ix <- which(upper.tri(A) & A, arr.ind = TRUE)
    data.frame(gene1 = rownames(A)[ix[, 1]], gene2 = colnames(A)[ix[, 2]], weight = if (is.null(weight))
        rep(1, nrow(ix))
    else weight[ix])
}

as_graphs <- function(a) lapply(seq_len(dim(a)[3]), function(k) a[, , k])

partial_cor <- function(P) {
    a <- -P/sqrt(outer(diag(P), diag(P)))
    diag(a) <- 0
    a
}

bfdr_mask <- function(P, alpha = 0.05) {
    up <- upper.tri(P)
    p <- pmin(1, pmax(0, P[up]))
    o <- order(p, decreasing = TRUE)
    ends <- which(!duplicated(p[o], fromLast = TRUE))
    ok <- ends[(cumsum(1 - p[o])/seq_along(p))[ends] <= alpha]
    a <- matrix(FALSE, nrow(P), ncol(P))
    if (length(ok))
        a[up] <- p >= p[o[max(ok)]]
    a | t(a)
}

match_graphs <- function(T, C) {
    D <- outer(1:3, 1:3, Vectorize(function(i, j) mean((T[[i]] - C[[j]])^2)))
    v <- t(apply(perms3, 1, function(p) {
        d <- D[cbind(1:3, p)]
        c(max(d), sum(d))
    }))
    as.integer(perms3[order(v[, 1], v[, 2])[1], ])
}

standardize <- function(y) {
    z <- scale(y)
    bad <- !is.finite(attr(z, "scaled:scale")) | attr(z, "scaled:scale") == 0
    z[, bad] <- 0
    if (any(!is.finite(z)))
        stop("Nonfinite standardized data")
    unname(z)
}

cluster_input <- function(obs, seed) {
    set.seed(seed)
    lapply(c("T", "C"), function(d) {
        y <- standardize(t(log1p(obs[[paste0("X_", d)]])))
        pc <- prcomp(y, center = FALSE, scale. = FALSE, rank. = min(10, ncol(y)))
        lab <- kmeans(pc$x, centers = 3, nstart = 20, iter.max = 100)$cluster
        list(labels = lab, Y = lapply(1:3, function(k) {
            v <- y[lab == k, , drop = FALSE]
            if (nrow(v) < 5) stop("Too few cells in estimated cluster")
            standardize(v)
        }))
    })
}

ebic_graph <- function(P, S, n, gamma) {
    det <- determinant(P, logarithm = TRUE)
    if (det$sign <= 0)
        return(Inf)
    e <- sum(abs(P[upper.tri(P)]) > 1e-05)
    p <- nrow(P)
    n * (sum(S * P) - as.numeric(det$modulus)) + e * log(n) + 4 * gamma * e * log(p)
}

fit_glasso <- function(Y, cfg) {
    S <- crossprod(Y)/nrow(Y)
    diag(S) <- pmax(diag(S), 1e-06)
    fits <- lapply(cfg$lambda_grid, function(lam) {
        f <- glasso::glasso(S, rho = lam, penalize.diagonal = FALSE)
        if (!is.null(f$errflag) && f$errflag != 0)
            stop("glasso solver error")
        f$wi
    })
    val <- vapply(fits, ebic_graph, 0, S = S, n = nrow(Y), gamma = cfg$ebic_gamma)
    fits[[which.min(val)]]
}

fit_jgl_pair <- function(Y, cfg, penalty) {
    best <- Inf
    out <- NULL
    S <- lapply(Y, function(y) crossprod(y)/nrow(y))
    for (l1 in cfg$lambda_grid) for (l2 in cfg$lambda2_grid) {
        f <- JGL::JGL(Y, penalty = penalty, lambda1 = l1, lambda2 = l2, weights = "sample.size", penalize.diagonal = FALSE,
            maxiter = 1000, tol = 1e-05, return.whole.theta = TRUE)
        value <- sum(vapply(1:2, function(i) ebic_graph(f$theta[[i]], S[[i]], nrow(Y[[i]]), cfg$ebic_gamma), 0))
        if (value < best) {
            best <- value
            out <- f$theta
        }
    }
    if (is.null(out))
        stop("No valid JGL fit")
    out
}

fit_mixture_ggm <- function(Y, cfg, seed) {
    n <- nrow(Y)
    p <- ncol(Y)
    K <- 3L
    best <- NULL
    best_bic <- Inf
    records <- list()
    estep <- function(mu, P, pi) {
        lp <- vapply(1:K, function(k) {
            D <- sweep(Y, 2, mu[k, ], "-")
            log(pi[k]) + 0.5 * (as.numeric(determinant(P[[k]], logarithm = TRUE)$modulus) - p * log(2 * pi_const) -
                rowSums((D %*% P[[k]]) * D))
        }, numeric(n))
        mx <- apply(lp, 1, max)
        a <- exp(lp - mx)
        den <- rowSums(a)
        list(w = a/den, ll = sum(mx + log(den)))
    }
    pi_const <- base::pi
    for (lam in cfg$mixture_lambda) for (start in seq_len(cfg$mixture_starts)) {
        set.seed(seed + start)
        lab <- kmeans(Y, 3, nstart = 1, iter.max = 100)$cluster
        w <- sapply(1:K, function(k) as.numeric(lab == k))
        old <- -Inf
        converged <- FALSE
        invalid <- FALSE
        for (it in seq_len(cfg$mixture_maxiter)) {
            nk <- colSums(w)
            if (any(nk < 5)) {
                invalid <- TRUE
                break
            }
            mu <- sweep(t(w) %*% Y, 1, nk, "/")
            mix <- nk/n
            P <- lapply(1:K, function(k) {
                D <- sweep(Y, 2, mu[k, ], "-")
                S <- crossprod(D * sqrt(w[, k]))/nk[k]
                diag(S) <- diag(S) + 1e-06
                f <- glasso::glasso(S, rho = lam * n/nk[k], penalize.diagonal = FALSE)
                if (f$errflag != 0)
                  stop("Mixture GGM glasso solver failed")
                f$wi
            })
            E <- estep(mu, P, mix)
            penalty <- n * lam/2 * sum(vapply(P, function(x) sum(abs(x)) - sum(abs(diag(x))), 0))
            objective <- E$ll - penalty
            w <- E$w
            if (is.finite(old) && abs(objective - old) < cfg$mixture_tol * (1 + abs(old))) {
                converged <- TRUE
                break
            }
            old <- objective
        }
        if (invalid)
            next
        edges <- sum(vapply(P, function(x) sum(abs(x[upper.tri(x)]) > 1e-05), 0L))
        df <- K * (2 * p) + (K - 1) + edges
        bic <- -2 * E$ll + df * log(n) + 4 * cfg$ebic_gamma * edges * log(p)
        records[[length(records) + 1L]] <- data.frame(lambda = lam, start = start, iterations = it, converged = converged,
            EBIC = bic)
        if (converged && bic < best_bic) {
            best_bic <- bic
            best <- list(P = P, labels = max.col(w, ties.method = "first"), lambda = lam, iterations = it)
        }
    }
    if (is.null(best))
        stop("No converged mixture GGM candidate")
    best$tuning <- do.call(rbind, records)
    best
}

probability_diff <- function(f) {
    G <- nrow(f$scores$T[[1]])
    joint <- lapply(1:4, function(x) array(0, c(G, G, 9)))
    names(joint) <- c("p00", "p10", "p01", "p11")
    for (k in 1:3) for (c in 1:3) {
        h <- (k - 1) * 3 + c
        pt <- f$scores$T[[k]]
        pc <- f$scores$C[[c]]
        joint$p00[, , h] <- (1 - pt) * (1 - pc)
        joint$p10[, , h] <- pt * (1 - pc)
        joint$p01[, , h] <- (1 - pt) * pc
        joint$p11[, , h] <- pt * pc
    }
    dummy <- list(K = 3L, G = G, gene_names = rownames(f$scores$T[[1]]), inference_mode = "staged", joint_edge_probabilities = joint,
        sampling_schedule = list(retained_draws = CFG$num_save, feedback_disabled_for_retained = TRUE), matching_result_df = data.frame(treatment_cluster = 1:3,
            matched_control_cluster = f$mapping))
    if (is.null(dummy$gene_names))
        dummy$gene_names <- sprintf("G%03d", 1:G)
    scNPM::select_differential_edges(dummy, fdr_level = CFG$bfdr, min_draws = 20L, rule = "delta_threshold")$table
}
