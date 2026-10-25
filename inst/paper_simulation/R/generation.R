# Extracted from the production simulation functions; see PROVENANCE.md.

take <- function(x, n) {
    if (!n)
        return(x[integer(0)])
    if (n > length(x))
        stop("Too few eligible pairs.")
    x[sample.int(length(x), n)]
}

min_eigen <- function(x) min(eigen(x, symmetric = TRUE, only.values = TRUE)$values)

adjacency <- function(x) {
    a <- x != 0
    diag(a) <- FALSE
    a
}

wu_module <- function(type) {
    B <- diag(2, 10)
    for (j in 2:10) for (t in seq_len(j - 1L)) {
        w <- switch(type, d = if (j - t <= 5) 0.7 else 0, c = if (j - t == 1 || (j == 10 && t == 1)) 0.9 else 0,
            s = if (t == 1) 0.6 else 0, nd = if (j - t <= 5 && t <= 5) {
                if (t <= 2) -0.6 else 0.6
            } else 0, stop("Unknown module"))
        B[j, t] <- B[t, j] <- w
    }
    stopifnot(min_eigen(B) > 0)
    B
}

stress_overlap_templates <- function() {
    ring <- matrix(0, 100L, 100L)
    diag(ring) <- 2
    out <- vector("list", 3L)
    for (b in seq_len(10L)) {
        ids <- (10L * (b - 1L) + 1L):(10L * b)
        block <- wu_module("c")
        block[row(block) != col(block) & block != 0] <- 0.6
        ring[ids, ids] <- block
    }
    out[] <- list(ring)
    for (b in seq_len(10L)) {
        ids <- (10L * (b - 1L) + 1L):(10L * b)
        local <- ring[ids, ids]
        choices <- take(which(upper.tri(local) & local == 0), 3L)
        for (k in seq_len(3L)) {
            jt <- arrayInd(choices[k], dim(local))
            j <- ids[jt[1L]]
            t <- ids[jt[2L]]
            out[[k]][j, t] <- out[[k]][t, j] <- 0.25
        }
    }
    up <- upper.tri(ring)
    aa <- lapply(out, adjacency)
    stopifnot(all(vapply(out, min_eigen, 0) > 0.1), all(vapply(aa, function(a) sum(a & up), 0L) == 110L), sum(Reduce(`&`,
        aa) & up) == 100L)
    list(precision = out, modules = rep(list(rep("ring_plus_chord", 10L)), 3L))
}

stress_perturb_all_blocks <- function(B, cfg) {
    up <- upper.tri(B)
    blocks <- rep(seq_len(10L), each = 10L)
    eligible <- up & outer(blocks, blocks, "==")
    on <- which(eligible & B != 0)
    off <- which(eligible & B == 0)
    nchange <- round(cfg$change_fraction * sum(up & B != 0))
    nd <- floor(nchange/2)
    na <- nchange - nd
    actions <- c(rep(c("add", "delete"), nd), rep("add", na - nd))
    C <- B
    removed <- added <- integer(0L)
    proposals <- 0L
    for (action in actions) {
        accepted <- FALSE
        while (!accepted) {
            proposals <- proposals + 1L
            if (proposals > cfg$max_proposals)
                stop("Positive-definite perturbation budget exhausted.")
            pool <- if (action == "delete")
                setdiff(on, removed)
            else setdiff(off, added)
            ix <- take(pool, 1L)
            jt <- arrayInd(ix, dim(B))
            j <- jt[1L]
            t <- jt[2L]
            value <- 0
            if (action == "add") {
                ids <- which(blocks == blocks[j])
                block <- B[ids, ids]
                value <- take(block[upper.tri(block) & block != 0], 1L)
            }
            trial <- C
            trial[j, t] <- trial[t, j] <- value
            accepted <- min_eigen(trial) > cfg$min_eigenvalue
            if (accepted) {
                C <- trial
                if (action == "delete")
                  removed <- c(removed, ix)
                else added <- c(added, ix)
            }
        }
    }
    list(precision = C, removed = removed, added = added, proposals = proposals)
}

make_networks <- function(seed) {
    cfg <- DATA_CFG
    set.seed(seed)
    tpl <- stress_overlap_templates()
    TT <- CC <- vector("list", 3)
    audit <- list()
    for (k in 1:3) {
        C <- tpl$precision[[k]]
        c <- cfg$control_label_permutation[k]
        z <- stress_perturb_all_blocks(C, cfg)
        T <- z$precision
        TT[[k]] <- T
        CC[[c]] <- C
        up <- upper.tri(C)
        a <- adjacency(T)
        b <- adjacency(C)
        audit[[k]] <- data.frame(T = k, C = c, edges_T = sum(a & up), edges_C = sum(b & up), gains = sum(a & !b &
            up), losses = sum(!a & b & up), change_fraction = sum(xor(a, b) & up)/sum(b & up), min_eigen_T = min_eigen(T),
            min_eigen_C = min_eigen(C))
    }
    genes <- sprintf("G%03d", 1:100)
    nm <- function(x) {
        dimnames(x) <- list(genes, genes)
        x
    }
    TT <- lapply(TT, nm)
    CC <- lapply(CC, nm)
    names(TT) <- paste0("T", 1:3)
    names(CC) <- paste0("C", 1:3)
    A <- lapply(c(TT, CC), adjacency)
    eligible <- Reduce(`|`, lapply(CC, adjacency))
    up <- upper.tri(eligible)
    set.seed(seed + 300L)
    ids <- take(which(up & eligible), round(cfg$ppi_coverage * sum(up & eligible)))
    noise <- take(which(up & !eligible), round(length(ids) * cfg$ppi_false_fraction/(1 - cfg$ppi_false_fraction)))
    ppi <- matrix(FALSE, 100, 100, dimnames = list(genes, genes))
    ppi[c(ids, noise)] <- TRUE
    ppi <- ppi | t(ppi)
    audit <- do.call(rbind, audit)
    stopifnot(all(audit$edges_T == 110), all(audit$edges_C == 110), all(audit$gains == 11), all(audit$losses ==
        11), all(audit$change_fraction == 0.2), all(eligible[ids]), !any(eligible[noise]))
    list(precision = list(T = TT, C = CC), adjacency = A, ppi = ppi, intersection = Reduce(`&`, A), union = Reduce(`|`,
        A), ppi_eligible = eligible, ppi_true_ids = ids, ppi_false_ids = noise, matching = data.frame(T = 1:3,
        C = cfg$control_label_permutation), modules = tpl$modules, audit = audit, config = cfg, seed = seed)
}

generate_original_means <- function(cfg, seed) {
    set.seed(seed)
    base <- rgamma(100, shape = cfg$gamma_shape, rate = cfg$gamma_rate)
    mu <- cbind(base, base, base)
    mu[1:10, 2] <- rnorm(10, 4, 0.2)
    mu[11:20, 2] <- rnorm(10, 5, 0.2)
    mu[16:30, 3] <- rnorm(15, 6, 0.2)
    dimnames(mu) <- list(sprintf("G%03d", 1:100), paste0("type", 1:3))
    mu
}

generate_expression <- function(networks, means, seed) {
    cfg <- networks$config
    set.seed(seed)
    lambda0 <- rnorm(100, 1, 0.1)
    lambda1 <- rnorm(100, -1, 0.1)
    while (any(lambda1 >= 0)) {
        ii <- which(lambda1 >= 0)
        lambda1[ii] <- rnorm(length(ii), -1, 0.1)
    }
    expected_zero <- function(shift) {
        vals <- numeric(2)
        for (di in 1:2) {
            d <- c("T", "C")[di]
            vals[di] <- sum(vapply(1:3, function(k) {
                j <- if (d == "T") k else cfg$control_label_permutation[k]
                variance <- diag(chol2inv(chol(networks$precision[[d]][[j]])))
                cfg$proportions[k] * mean(pnorm((lambda0 + shift + lambda1 * means[, k])/sqrt(1 + lambda1^2 * variance)))
            }, 0))
        }
        mean(vals)
    }
    shift <- uniroot(function(a) expected_zero(a) - cfg$zero_target, c(-100, 100), tol = 1e-10)$root
    lambda0 <- lambda0 + shift
    observed <- truth <- list()
    summaries <- list()
    for (d in c("T", "C")) {
        n <- cfg$n_per_condition
        canonical <- sample.int(3, n, replace = TRUE, prob = cfg$proportions)
        local <- if (d == "T")
            canonical
        else cfg$control_label_permutation[canonical]
        theta <- matrix(NA_real_, 100, n)
        for (k in 1:3) {
            ii <- which(canonical == k)
            if (!length(ii))
                next
            j <- if (d == "T")
                k
            else cfg$control_label_permutation[k]
            Sigma <- chol2inv(chol(networks$precision[[d]][[j]]))
            theta[, ii] <- means[, k] + t(matrix(rnorm(length(ii) * 100), length(ii), 100) %*% chol(Sigma))
        }
        prob_zero <- pnorm(lambda0 + lambda1 * theta)
        X <- exp(theta)
        X[matrix(runif(100 * n), 100, n) < prob_zero] <- 0
        dimnames(X) <- dimnames(theta) <- list(sprintf("G%03d", 1:100), sprintf("%s_cell%04d", d, 1:n))
        stopifnot(all(is.finite(X)), all(X >= 0))
        observed[[d]] <- X
        truth[[d]] <- list(labels = local, canonical = canonical, theta = theta)
        summaries[[d]] <- data.frame(condition = d, n = n, n1 = sum(local == 1), n2 = sum(local == 2), n3 = sum(local ==
            3), zero_median = median(colMeans(X == 0)), zero_mean = mean(X == 0))
    }
    list(observed = list(X_T = observed$T, X_C = observed$C, PPI = networks$ppi), truth = list(cells = truth, means = means,
        lambda0 = lambda0, lambda1 = lambda1, networks = networks, zero_target = cfg$zero_target, dropout_intercept_shift = shift),
        summary = do.call(rbind, summaries))
}

fresh_expression <- function(old_truth, seed) {
    set.seed(seed)
    n <- 1500L
    G <- 100L
    net <- old_truth$networks
    net$config$n_per_condition <- n
    means <- old_truth$means
    l0 <- old_truth$lambda0
    l1 <- old_truth$lambda1
    obs <- cells <- audit <- list()
    for (d in c("T", "C")) {
        canonical <- sample.int(3L, n, replace = TRUE, prob = net$config$proportions)
        labels <- if (d == "T")
            canonical
        else net$config$control_label_permutation[canonical]
        theta <- matrix(NA_real_, G, n)
        for (k in 1:3) {
            ix <- which(canonical == k)
            j <- if (d == "T")
                k
            else net$config$control_label_permutation[k]
            sigma <- chol2inv(chol(net$precision[[d]][[j]]))
            theta[, ix] <- means[, k] + t(matrix(rnorm(length(ix) * G), length(ix), G) %*% chol(sigma))
        }
        X <- exp(theta)
        X[matrix(runif(G * n), G, n) < pnorm(l0 + l1 * theta)] <- 0
        dimnames(X) <- dimnames(theta) <- list(sprintf("G%03d", 1:G), sprintf("%s_new1500_%04d", d, 1:n))
        obs[[d]] <- X
        cells[[d]] <- list(labels = labels, canonical = canonical, theta = theta)
        audit[[d]] <- data.frame(condition = d, n = n, n1 = sum(labels == 1), n2 = sum(labels == 2), n3 = sum(labels ==
            3), zero_mean = mean(X == 0))
    }
    list(observed = list(X_T = obs$T, X_C = obs$C, PPI = net$ppi), truth = list(cells = cells, means = means, lambda0 = l0,
        lambda1 = l1, networks = net, zero_target = 0.4, dropout_intercept_shift = old_truth$dropout_intercept_shift),
        summary = do.call(rbind, audit))
}
