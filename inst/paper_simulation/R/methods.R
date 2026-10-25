# Production adapter with the final union-PPI and prior settings resolved.
fit_adapter <- function(obs, method, seed, folder, cfg = CFG) {
    set.seed(seed)
    G <- nrow(obs$X_T)
    if (grepl("^scNPM", method)) {
        stage(folder, paste(method, "joint sampling"))
        protein <- edge_table(obs$PPI)
        protein <- data.frame(node1 = protein$gene1, node2 = protein$gene2, combined_score = protein$weight)
        f <- scNPM::run_scNPM(cbind(obs$X_T, obs$X_C), cbind(condition = c(rep(1, ncol(obs$X_T)), rep(0, ncol(obs$X_C)))),
            K = 3L, ds_seed = seed, num_iter = cfg$num_iter, num_save = cfg$num_save, protein_input = if (method ==
                "scNPM_noPPI")
                NULL
            else protein, ppi_prior_low = 0.05, ppi_prior_high = if (method == "scNPM_noPPI")
                0.05
            else 0.5, edge_prob_adj = 1, adapt_min = cfg$adapt_min, adapt_max = cfg$adapt_max, freeze_patience = cfg$freeze_patience,
            final_burnin = cfg$final_burnin, inference_mode = "staged", matching_score = "pairwise", matching_method = "lexicographic",
            freeze_policy = "consensus", fdr_level = 0.05, group_fdr_level = 0.05, fdr_min_draws = min(200L, cfg$num_save),
            network_prior = if (method == "scNPM_independent_fixedNPM")
                "independent"
            else "paired", paired_discordance = 0.75, shared_edge_correlation = 0.995, ppi_reference = "union",
            verbose = TRUE, progress_every = 60, progress_file = file.path(folder, "iterations.csv"))
        stopifnot(f$network_prior == if (method == "scNPM_independent_fixedNPM")
            "independent"
        else "paired", f$provenance$settings$ppi_reference == "union", f$provenance$settings$edge_prob_adj == 1,
            all(f$stage_history$feedback_active == 0))
        ans <- list(labels = list(T = as.integer(f$group_treat_est), C = as.integer(f$group_ctrl_est)), scores = list(T = as_graphs(f$group_fdr$probabilities$treatment),
            C = as_graphs(f$group_fdr$probabilities$control)), masks = list(T = as_graphs(f$group_fdr$masks$treatment),
            C = as_graphs(f$group_fdr$masks$control)), mapping = as.integer(f$matching_result_df$matched_control_cluster),
            frozen_at = f$sampling_schedule$frozen_at, differential = f$differential_fdr$table, selection = "directional BFDR 0.05",
            schedule = f$sampling_schedule, provenance = f$provenance, diagnostics = list(freeze_reason = f$sampling_schedule$freeze_reason,
                final_max_smd = max(f$smd_cost_matrix[cbind(1:3, f$matching_result_df$matched_control_cluster)]),
                batch_mcse = f$differential_fdr$table$batch_mcse))
        {
            ans$adaptation_audit <- list(history = f$matching_history[seq_len(f$sampling_schedule$frozen_at), ,
                drop = FALSE], frozen_cost = f$frozen_smd_cost_matrix)
            ans$matching_history <- f$matching_history
            ans$stage_history <- f$stage_history
            rm(f)
        }
        gc(FALSE)
    }
    else if (method == "BLGGM") {
        fs <- list()
        for (d in c("T", "C")) {
            stage(folder, paste("BLGGM condition", d, "iterations", cfg$num_iter))
            set.seed(seed + match(d, c("T", "C")))
            fs[[d]] <- BLGGM::BLGGM(obs[[paste0("X_", d)]], n_celltype = 3, num_iterations = cfg$num_iter, num_burnin = cfg$num_iter -
                cfg$num_save, collect_post_sample = FALSE, print_label = TRUE, print_per_iteration = 100, num_threads = cfg$threads)
        }
        scores <- lapply(fs, function(f) as_graphs(f$cell_type_edge_post_prob))
        pc <- lapply(fs, function(f) lapply(as_graphs(f$cell_type_precision_matr), partial_cor))
        ans <- list(labels = lapply(fs, function(f) as.integer(f$cell_labels)), scores = scores, masks = lapply(scores,
            function(pp) lapply(pp, bfdr_mask, alpha = cfg$bfdr)), mapping = match_graphs(pc$T, pc$C), frozen_at = NA_integer_)
        ans$differential <- probability_diff(ans)
        ans$selection <- "independent posterior product; directional BFDR 0.05"
    }
    else if (method == "AJGM") {
        stage(folder, "AJGM optimization")
        ans <- fit_ajgm(obs, cfg, seed)
    }
    else if (method == "mixtureGGM") {
        stage(folder, "Mixture GGM EM and EBIC selection")
        ff <- lapply(c("T", "C"), function(d) fit_mixture_ggm(standardize(t(log1p(obs[[paste0("X_", d)]]))), cfg,
            seed))
        pc <- lapply(ff, function(z) lapply(z$P, partial_cor))
        ans <- list(labels = lapply(ff, function(z) as.integer(z$labels)), scores = lapply(pc, function(z) lapply(z,
            abs)), masks = lapply(ff, function(z) lapply(z$P, function(a) {
            m <- abs(a) > 1e-05
            diag(m) <- FALSE
            m
        })), mapping = match_graphs(pc[[1]], pc[[2]]), frozen_at = NA_integer_, tuning = lapply(ff, function(z) z$tuning))
    }
    else {
        stage(folder, paste(method, "PCA/kmeans clustering and EBIC selection"))
        cl <- cluster_input(obs, seed)
        P <- lapply(cl, function(a) lapply(a$Y, fit_glasso, cfg = cfg))
        mp <- match_graphs(lapply(P[[1]], partial_cor), lapply(P[[2]], partial_cor))
        if (method == "JGL_fused")
            for (k in 1:3) {
                z <- fit_jgl_pair(list(cl[[1]]$Y[[k]], cl[[2]]$Y[[mp[k]]]), cfg, "fused")
                P[[1]][[k]] <- z[[1]]
                P[[2]][[mp[k]]] <- z[[2]]
            }
        ans <- list(labels = lapply(cl, function(z) z$labels), scores = lapply(P, function(a) lapply(a, function(p) abs(partial_cor(p)))),
            masks = lapply(P, function(a) lapply(a, function(p) {
                m <- abs(p) > 1e-05
                diag(m) <- FALSE
                m
            })), mapping = mp, frozen_at = NA_integer_)
    }
    names(ans$labels) <- names(ans$scores) <- names(ans$masks) <- c("T", "C")
    stopifnot(identical(sort(ans$mapping), 1:3))
    if (is.null(ans$differential)) {
        ij <- which(upper.tri(obs$PPI), arr.ind = TRUE)
        ans$differential <- do.call(rbind, lapply(1:3, function(k) {
            c <- ans$mapping[k]
            pt <- ans$scores$T[[k]][ij]
            pc <- ans$scores$C[[c]][ij]
            mt <- ans$masks$T[[k]][ij]
            mc <- ans$masks$C[[c]][ij]
            delta <- pt - pc
            data.frame(treatment_cluster = k, matched_control_cluster = c, gene_i = ij[, 1], gene_j = ij[, 2],
                selected = xor(mt, mc), direction = ifelse(mt & !mc, "treatment", ifelse(mc & !mt, "control", ifelse(delta >
                  0, "treatment", "control"))), abs_delta_p = abs(delta), p_treatment_only = pmax(delta, 0), p_control_only = pmax(-delta,
                  0))
        }))
        ans$selection <- "support difference; continuous partial-correlation contrast for ranking"
    }
    ans
}
