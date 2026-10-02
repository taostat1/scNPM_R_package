# scNPM 1.0.6

* Address CRAN review: keep the 39-character title, make all help examples
  self-contained and unwrapped, and add quiet-by-default verbosity across R
  and C++ with optional progress emitted as suppressible message conditions.

* Select group networks by BFDR separately for each condition and subtype,
  using Rao-Blackwell marginal probabilities (p10+p11 or p01+p11). Evaluate
  complete probability tie blocks and retain the largest set with mean 1-p
  at most the target (default 0.05); use p >= the smallest selected probability.
* Add select_group_edges(), automatically return group_fdr in staged fits,
  and expose summary/edge tables, thresholds, full curves and symmetric masks.
  Existing staged 1.0.5 fits can be post-processed without new sampling.
* Group plots now require saved BFDR selections by default. The earlier 0.5
  display is available only with selection = "frequency". Allocation-frequency
  arrays retain their original meaning for compatibility.
* Add the full retinal example, including group and differential exports.
  Consolidate algorithm documentation into algorithm.md and exclude historical
  version files and unrelated figure-generation scripts from distribution.
* Matching, freezing, baseline sampling and differential selection are unchanged.

# scNPM 1.0.5

* Default to complete lexicographic matching: minimize the largest SMD first,
  then total SMD, with deterministic ties and no SMD cap. K<=8 uses exact
  permutation enumeration; larger K uses bottleneck feasibility and restricted
  minimum-sum assignment. Preserve frozen maps and expose freeze/final costs.
* Select differential networks by strict absolute-contrast thresholds, calibrated
  by the mean posterior directional error across all K cluster-edge families.
  Scan every distinct cutoff, preserve entire tied-score blocks, and choose the
  largest feasible network. Return the threshold and full calibration curve.
* Default FDR scope to the frozen correspondence. Retain all_matchings for
  sensitivity analysis, posterior_error for 1.0.4 selection, and greedy for
  matching compatibility. Keep the 200-draw guard and plotting decision masks.
* Add independent assignment-oracle, infinity/tie, nonmonotone cutoff and
  stage-integration tests. Runtime comparisons use full AMD data with matched
  seeds, threads and iteration budgets; short runs do not establish convergence.

# scNPM 1.0.4

* Default staged inference: align labels to fixed initial mean anchors, freeze only after complete-map stability, disable feedback, and discard post-switch warmup before retention. Unstable adaptation fails explicitly.
* Compute joint 00/10/01/11 probabilities for all K*K reference pairs by averaging conditional products within final baseline draws; retain empirical and batch diagnostics.
* Add posterior directional error selection across all K gene-pair families, with a conservative all-matchings default and conditional frozen-map option. This is an estimated Bayesian error bound, not a frequentist guarantee.
* Return J_aug diagnostics and stage/mapping histories. The augmented likelihood differs from the observed marginal likelihood and is not forced to increase.
* Keep inference_mode="coupled" and differential_method="frequency" for 1.0.3 compatibility. Retain eight shared reference features when K=4.

# scNPM 1.0.3

## Shared network reference features

* Change the default propensity network features to centered quadratic forms
  against all current treatment and control networks, ordered T1..TK,C1..CK.
  Every cell uses the same 2K reference models; K = 4 yields eight features.
  The references are refreshed from the current state at each matching step.
* Standardize each predictor using the pooled cell mean and sample standard
  deviation, retaining the existing other cell covariates and SMD matching.
* Add `network_covariates = "legacy"` for the previous single uncentered
  feature, and expose `propensity_ridge` (default 1e-6) for logistic IRLS
  stabilization. The feature mode and setting are recorded in provenance.
* Return final `network_covariates`, `network_covariate_mode` and
  `propensity_scaling`. Report each reference feature in the covariate SMD
  matrix. `compute_qf_smd()` returns a long feature-specific table for new
  all-reference fits, while retaining the old format for legacy results.
* Retain `qf_global` as a compatibility diagnostic. It is not the propensity
  network predictor in all-reference mode.

## Implementation and validation

* Compute centered reference features with batched matrix operations.
* Skip output-only balance diagnostics during internal sampling iterations;
  compute the complete diagnostics for the final matching call.
* Add independently calculated quadratic-form, pooled-scaling, translation,
  reference-order, duplicate-reference and public-interface regression tests.
* Keep the observation and sampling mechanisms from 1.0.2. Changing feature
  mode can change matching and subsequent edge feedback; compare modes using
  separately run chains. Short runtime benchmarks are not convergence tests.

# scNPM 1.0.2

## Sampling and network prior corrections

* Correct the spike-and-slab allocation probability by including the missing
  slab-weight term `log(xi)` in its log odds. The probability now equals
  `xi * slab / (xi * slab + (1 - xi) * spike)`.
* Add `ppi_prior_low` and `ppi_prior_high` to `run_scNPM()`, with defaults
  0.05 and 0.5, constrained by `0 < low <= high <= 1 - 1e-6`. These define the baseline inclusion weights for absent and
  present PPI edges. The resulting edge-specific weights are used at every
  network update, rather than using PPI adjacency only for initialization.
  The defaults are working analysis settings, not empirically optimized values.
* Keep the uncentered network covariate `q = theta' Omega theta` and common
  column standardization across the pooled cells in the propensity fit.
  Matching feedback remains active throughout retained sampling. HMC uses the
  supplied fixed step sizes and numbers of leapfrog steps.

* Include treatment clusters with zero matching-priority scores in the
  greedy matching pass, rather than leaving them unmatched solely because
  their priority is zero.
* Update zero-associated latent expression for all-zero cells, and draw
  empty-cluster parameters from their corresponding prior updates instead
  of leaving stale parameter values in place.
* Use seeded thread-local random generation in parallel updates and R's
  seeded random stream for Dirichlet draws. Reproducibility requires a fixed
  software environment and thread count.
* Evaluate the probit observation terms on the log scale and stabilize the
  associated likelihood calculations in the tails.

## Input handling and reproducibility

* Correct the use of global versus within-condition cell indices during
  initialization and the indexing used to initialize mixture proportions.
* Check that expression values, covariates, and numeric settings are finite
  and admissible, and validate small or duplicated clustering inputs.
* PPI inputs now support local files, TSV text, and data frames; the earlier
  implementation accepted only data frames despite documenting all formats.
* Add `protein_adjacency`, `edge_prior`, and `provenance` to the result.
  Provenance records analysis settings, the package version, and fingerprints
  of serialized inputs and the native library when available.

* Align `compute_qf_smd()` with the C++ boundary convention: singleton samples
  have variance zero, equal constant samples have SMD zero, and unequal
  constant samples have infinite SMD. Empty or unmatched pairs remain missing.
  Summary-based approximations are no longer described as necessarily close
  to the final matching state.

## Interpretation and migration

* Describe `edge_treat_est` and `edge_ctrl_est` as sampled edge inclusion
  frequencies. The matching-guided transition is not established to preserve
  a compatible joint posterior.
* Clarify that differential cutoffs of +0.5 and -0.5 are descriptive selection
  thresholds, not FDR control. The reported matching map comes from the final
  state and need not stabilize as iterations increase.
* Results from 1.0.1 and 1.0.2 should not be pooled or treated as the same run.
  The allocation correction and persistent PPI prior change the sampling
  transition; rerun an analysis to obtain results for 1.0.2. Save the package
  version, inputs, seed, and settings with the output.
* Validation for this revision uses small synthetic examples and package
  checks. The archived AMD, IPF, and DN analyses have not been rerun for 1.0.2.

# scNPM 1.0.1

## New features

* New exported function `compute_qf_smd()`: standardized mean difference of
  the network quadratic-form covariate for each matched cluster pair.
  Two modes: exact (when the new `qf_global` / `final_group_label` outputs
  of `run_scNPM()` are available) and posterior-approximation (for result
  files saved with earlier versions).
* `run_scNPM()` now reports the network quadratic form as an additional
  covariate in `smd_covariate_results` (last column, `qf_network`), so
  balance of the matching-guiding network statistic can be assessed
  together with the QC covariates. The C++ change is output-only: the
  mathematical sampling updates and matching logic were unchanged by these
  additional output fields.
* Network plots (`plot_group_network`, `plot_diff_network`) now default to
  the original light-blue edges (`"#ADD8E6"`, as in the reference
  implementation) instead of red/blue or green/orange. Both colours remain
  fully user-controllable via the `edge_color`, `solid_color` and
  `dashed_color` arguments.
* `run_scNPM()` additionally returns `qf_global` (per-cell network
  quadratic form), `propensity_scores` and `final_group_label` from the
  final matching call.

## Bug fixes

* `network_propensity_matching_by_difficulty()` now computes a real logistic
  propensity score by IRLS (Newton-Raphson with ridge stabilisation) from the
  QC covariates and the network quadratic form. In 1.0.0 the propensity score
  was an i.i.d. `runif()` placeholder rather than a fitted propensity score.
* New output `smd_covariate_results`: per-covariate standardized mean
  difference for each matched cluster pair (rows = treatment clusters,
  columns = the QC covariates). This is the recommended balance diagnostic
  (values < 0.1 indicate good covariate balance). Also exposed:
  `propensity_scores`.
* `update_theta()` skips cells with zero counts in every gene instead of
  crashing with `Mat::max(): object has no elements` (the conditional HMC
  update previously called `.max()` on an empty observed-index vector).
* Spike-slab prior adjustment in `update_edge_with_matching_by_difficulty()`
  clamps the adjusted prior to `min(xi / rho, 1 - 1e-6)`, keeping
  `log(1 - xi)` finite in `calculate_edge_probability()` and preventing silent
  NaN edge probabilities when `xi / rho > 1`.

## CRAN compliance

* Package title shortened to below 65 characters.
* Package name quoted as 'scNPM' in the DESCRIPTION description field.
* All `\dontrun{}` examples replaced by self-contained `\donttest{}`
  examples that run in a few seconds.
* Progress messages in `run_scNPM()` switched from `cat()` to `message()`
  so output goes to `stderr()` and can be suppressed.

# scNPM 1.0.0

* Initial CRAN submission.
* Core framework: two-stage analysis combining K-means cell clustering with an
  MCMC-based network propensity matching model (`run_scNPM`), differential
  network construction, and protein-protein interactome integration.
* Visualization functions for per-group and differential regulatory networks
  (`plot_group_network`, `plot_diff_network`).
* Example dataset: stratified 400-cell subset (200 PR control, 200 MR
  treatment) of a public human retina single-cell RNA-seq study.
* The `ds_seed` argument controls random initialization and sampling in
  `run_scNPM()`.
