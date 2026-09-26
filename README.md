# scNPM: Single-Cell Network Propensity Matching

Version **1.0.5** estimates heterogeneous gene conditional-dependence networks
and compares matched clusters across two conditions.

## Installation

```r
install.packages("scNPM_1.0.5.tar.gz", repos = NULL, type = "source")
# On Windows, install the supplied binary with type = "win.binary".
library(scNPM)
```

Source installation requires a compatible C++17 toolchain. See DESCRIPTION for
dependencies and NEWS.md for earlier releases.

## The default staged algorithm

1. Update the latent expression, mixture and network states with matching feedback.
   Align component labels to fixed initial mean anchors. The NPM matching is a
   complete assignment minimizing maximum cluster propensity-score SMD, then total SMD.
   Exact ties use a fixed traversal; no SMD cap is applied.
2. Freeze a complete bijection after `freeze_patience` consecutive identical
   assignments, no earlier than `adapt_min`. Every cluster must be nonempty at
   these checks. Stop with an error if this has not occurred by `adapt_max`.
3. From the next iteration, switch the edge adjustment multiplier to 1. This
   restores the edge-specific baseline PPI weights. Finish at least
   `final_burnin` iterations under this model before retaining the last
   `num_save` draws. The frozen matching remains the reported correspondence.
4. Select directional differential edges from final-stage joint edge-state
   probabilities, across all K clusters as one family.

The final sampling target is the baseline model without counterpart feedback.
Matching affects adaptation and the reported comparisons; it does not keep
altering the final network prior. Re-run analyses when moving from 1.0.3.

The observation model is unchanged: positive values fix latent expression at
log X; zero values have sampled latent expression. HMC uses fixed step sizes,
fixed leapfrog counts and Metropolis acceptance correction.

## Network propensity covariates

For every cell in either condition, the default uses centered distances to all
2K current reference networks, ordered T1,...,TK,C1,...,CK:

$$q_{i,r}^{d;s}=(\theta_i^d-\mu_r^s)^\top\Omega_r^s(\theta_i^d-\mu_r^s),\quad s\in\{T,C\}.$$

K = 4 therefore uses eight network covariates. These and the optional QC
covariates enter the logistic propensity model. Each column is standardized
using its pooled mean and sample standard deviation over both conditions.
Constant columns become zero. `propensity_ridge = 1e-6` stabilizes IRLS.
`network_covariates = "legacy"` requests the original assigned-network,
uncentered quadratic form. Inspect both propensity-score and feature SMDs.

## Inputs and a staged analysis

`goal_gene_mt` is a finite, nonnegative gene-by-cell matrix on the original
expression scale, with unique gene names. Rows of `cell_t` follow the cell
columns: its first column codes control 0 and treatment 1; remaining numeric
columns are covariates. `protein_input` accepts a data frame, TSV text or a
local TSV file with node1, node2 and combined_score, or NULL.

```r
# Supply expression, cell_t and optional ppi for your data.
fit <- run_scNPM(
  goal_gene_mt = expression, cell_t = cell_t, protein_input = ppi, K = 4,
  num_iter = 5000, num_save = 2000,
  inference_mode = "staged", adapt_min = 500, adapt_max = 2000,
  freeze_patience = 20, final_burnin = 1000,
  differential_method = "bayes_fdr", fdr_level = 0.05,
  fdr_scope = "frozen", fdr_rule = "delta_threshold",
  matching_method = "lexicographic", diagnostic_every = 10, ds_seed = 123
)
fit$sampling_schedule
fit$matching_result_df
fit$differential_fdr$selected
fit$differential_fdr$estimated_error_bound
compute_qf_smd(fit)
plot(fit$criterion_trace$iteration, fit$criterion_trace$J_aug, type = "l")
saveRDS(fit, "scNPM_1.0.5_fit.rds")
```

Budgets satisfy `adapt_min <= adapt_max` and
`adapt_max + final_burnin <= num_iter - num_save`. Earlier freezing gives a
longer actual final burn-in; it does not move the retained window. Defaults
resolve from this available pre-retention budget; see `?run_scNPM`. Stability
of the matching and completion of a burn-in budget are separate from MCMC
convergence. Examine multiple chains and the retained diagnostic traces.

## Differential edges and error interpretation

For a treatment cluster k and an unordered gene pair, a treatment-direction
claim means slab allocation in k and spike allocation in its control
counterpart. A control-direction claim means the reverse. These are continuous
spike-and-slab states, not exact zero precision entries or gene absence.

The package averages within-draw conditional products for all K x K pairs.
It does not multiply marginal edge-frequency estimates. The more supported
direction under the frozen pairing supplies the reported direction.

- `fdr_scope = "frozen"` (default): use errors for the frozen counterpart.
- `fdr_scope = "all_matchings"`: sensitivity analysis using the largest error
  over all K control counterparts for the chosen direction.

The default `fdr_rule = "delta_threshold"` computes `delta_p = p10 - p01`,
then scans every distinct strict cutoff `abs(delta_p) > tau`. Equal contrasts
are never split. For each whole threshold set, average its directional errors
`1 - max(p10,p01)` (or the all-counterpart bound). Choose the set with the most
edges among those whose average is at most `fdr_level`. These averages need not
be monotone in contrast order; the scan does not stop at the first failure.

`fit$differential_fdr$threshold` records tau; `threshold_curve` records every
candidate's size and mean error. Infinity selects no edges. The cutoff is the
next excluded score, or zero if all nonzero scores are included; the operator
is strictly `>`. The plot uses selection masks and does not apply +/-0.5 again.
`fdr_rule = "posterior_error"` retains the error-sorted prefix from 1.0.4.
There is one direction per unordered edge and treatment cluster, with all K
clusters sharing one budget. Below `fdr_min_draws` (default 200), return the
exploratory curve but no threshold discoveries. This minimum is not a mixing
guarantee. The batch MCSE is a diagnostic, not a confidence bound.

This is an estimated **posterior expected directional FDP**, conditional on the
baseline model and label convention. A target of 0.05 is not a uniform
frequentist 5% guarantee across fixed network configurations, donors, model
misspecification or incomplete mixing. The supplied simulation report evaluates
this distinction explicitly. Donor dependence is not modeled by adding a donor
column to the propensity model.

```r
conditional <- select_differential_edges(fit, scope = "frozen")
conservative <- select_differential_edges(fit, scope = "all_matchings")
```

## Criterion diagnostics

`criterion_trace` records `log_augmented_likelihood`, `network_prior_kernel`,
`matching_loss` and their combination `J_aug = likelihood + prior - loss`.
The likelihood includes sampled missing latent expression and component labels;
it is not the integrated observed-data likelihood. The prior column is a
state-wise network factor, not the complete posterior density. SMDs are
recomputed against the frozen map after the switch. Undefined SMDs produce NA
loss/criterion values. Changes in stages also change the prior factor.

This trace describes interactions between modules. It does not define a joint
optimization objective. Random draws can decrease any of these values; no
updates are rejected merely because J_aug decreases. Global matching optimizes
(maximum SMD, total SMD) for the current cost matrix, not a common objective
across sampling iterations. Freezing preserves that map as costs change.

## Compatibility and outputs

Use `inference_mode = "coupled", differential_method = "frequency",
matching_method = "greedy"` for the
1.0.3 continuous-feedback algorithm. Its +/-0.5 selections are descriptive.
Posterior error selection rejects coupled-mode results. Old and new draws
must not be pooled. K=4 still has four clusters per condition, not eight.

Existing network, plotting, SMD and provenance outputs are retained. New fields:

| Field | Contents |
| --- | --- |
| `sampling_schedule` | Freeze iteration, minimum/actual final burn-in and retained interval. |
| `stage_history` | Per-iteration phase and feedback status. |
| `matching_history` | Per-iteration map, fixed after freezing. |
| `criterion_trace` | Augmented criterion components at recorded iterations. |
| `joint_edge_probabilities` | p00, p10, p01, p11 for every cluster pair, empirical checks and batch means. |
| `differential_fdr` | Edge table, threshold, full threshold curve, selections, target, scope and status. |
| `smd_cost_matrix`, `frozen_smd_cost_matrix` | All pair costs at the final state and at the freezing iteration. |
| `matching_objective` | Final-state maximum/total SMD under the reported map; frozen flag. |

`edge_diff_networks` and `plot_diff_network()` use the selected differential
edges in Bayesian mode. Keep the `scNPM_selection` attribute on the returned
`diff_matrix`: the plotting function reads these masks directly. Plain matrices
without that attribute retain the descriptive +/-0.5 plotting behavior.
`match_map` and `matching_result_df` use the frozen map,
while SMDs refer to the final sampled state under that map. The main edge
summaries are retained allocation frequencies; the new joint probabilities
use Rao-Blackwellization. Raw draws remain in `Result`.

The baseline PPI inclusion weight is
`ppi_prior_low + (ppi_prior_high - ppi_prior_low) * protein_adjacency`, with
defaults 0.05 and 0.5. It is used at every edge update, including the final
phase. The SPD matrix constraint is part of the joint network prior.

For formulae, assumptions and the error-bound proof, read
[algorithm-1.0.5.md](inst/algorithm-1.0.5.md). This file is also installed in the
package directory. Development repository:
[taostat1/scNPM_R_package](https://github.com/taostat1/scNPM_R_package).
License: GPL >= 3.
