# scNPM algorithm: staged inference and BFDR-selected networks

## Sampling stages

The state includes latent expression, cluster memberships, mixture weights,
component means, precisions, spike/slab allocations and observation parameters.
Nonzero expression fixes its latent value to log X; missing latent values at
zero observations are sampled using fixed-step HMC with Metropolis correction.
After every sweep, component labels are aligned to fixed initial mean anchors.
Both conditions use the same 2K network references in the propensity regression,
with pooled scaling of these distances and the measured cell covariates.

Matching and edge feedback operate during adaptation. All clusters must be
nonempty and the complete mapping must repeat for freeze_patience sweeps,
no earlier than adapt_min. Failure to stabilize by adapt_max stops the run.
Starting with the next sweep after freezing, the feedback multiplier is one
and the PPI-specific baseline weights are restored. With L=num_iter and
S=num_save, retained sweeps are L-S+1 through L. A freeze at F leaves L-S-F
baseline warmup sweeps, required to be at least final_burnin. The default
minimum retained sample count for BFDR selection is 200.

The diagnostic criterion is augmented likelihood plus a network-prior kernel
minus the matching SMD sum. It is not a jointly optimized objective, and
stochastic updates need not increase it.

## Group-network BFDR

For each condition d and subtype k, use all H=G*(G-1)/2 unordered gene pairs.
Let p[e] be the retained mean of the baseline conditional slab probability.
It is recovered from the saved joint quantities as p10+p11 for treatment or
p01+p11 for control. The marginal is independent of which counterpart is used;
select_group_edges checks this numerical consistency. It uses the frozen-pair
slice, as in the manuscript's group-network analysis, without resampling.

For E(c)={e:p[e]>=c}, estimate BFDR as mean(1-p[e]) over E(c). Sort probabilities
in decreasing order and evaluate complete tie-block endpoints. Keep the largest
prefix with estimated BFDR <= group_fdr_level (default 0.05). Its smallest
selected probability is c; c=Inf denotes an empty set. This is an inclusive
threshold, unlike the strict contrast threshold for differential edges.
For fewer than min_draws retained draws, keep the curve but return an empty set.

For any reported set, conditional on the data, linearity of expectation gives
E[FDP|X] = sum_{e in E} Pr(z[e]=0|X) / max(1,|E|). Independence across edges is
not needed. With estimated probabilities, if their maximum absolute error is
epsilon, the selected posterior expected FDP is at most alpha+epsilon.
The bound concerns spike/slab allocations under the baseline posterior.

Every condition/subtype and gene-panel fit is a separate family. Select before
restriction to shared genes; reuse the same saved mask for figures, exports
and sensitivity analysis. A union across clusters is a summary of those
selections and does not inherit a new per-union BFDR bound. Empirical allocation
frequencies remain in edge_treat_est/edge_ctrl_est; group_fdr contains the
Rao-Blackwell probabilities and BFDR masks. Default plotting requires these
saved selections; selection="frequency" explicitly enables descriptive plots.

## Matching

At each unfrozen matching call, construct D[k,c] from the cluster means and
sample variances of the current cell propensity scores. Minimize the ordered
pair (max_k D[k,phi(k)], sum_k D[k,phi(k)]) over complete bijections.
For K<=8, examine permutations in ascending order and update only for a strict
improvement in one of the two objectives. K=4 examines 24 possibilities.
For K>8, binary-search the sorted finite cost levels for a perfect matching,
then use minimum-sum assignment restricted to that minimum bottleneck. Normalize
allowed costs to [0,1] and use K+1 for forbidden entries, so a feasible allowed
assignment always beats one containing a forbidden edge. Ascending solver
traversal deterministically handles equal costs. Computation uses floating-point
costs; no SMD upper bound or user-chosen near-tie tolerance is imposed.

If no all-finite permutation exists, every complete permutation has infinite
maximum and total cost. Return the identity by the fixed tie rule. Missing SMDs
from empty clusters are kept missing in diagnostics and treated as infinite
for provisional assignment. Empty clusters prevent freezing. Finite but large
SMDs are not discarded. Stable nonempty complete maps freeze after the configured patience interval; feedback is disabled in subsequent iterations and final warmup precedes
retention. Label alignment still uses minimum squared distance to fixed anchors,
which is a separate assignment problem from matching T to C clusters.

The optimizer's guarantee refers to each current cost matrix. A frozen map is
not reoptimized at the final state. Its final maximum SMD need not remain
optimal, nor must the maximum SMD decrease from one iteration to the next.

## Differential-network threshold selection

Let e=(k,j,t), j<t. Using the final baseline-model draws, average within-draw
conditional products to estimate p10[e] and p01[e] for the frozen counterpart.
Set Delta[e]=p10[e]-p01[e], direction[e]=sign(Delta[e]), and
ell[e]=1-max(p10[e],p01[e]). The direction asserts unequal spike/slab allocation
states in the indicated groups, not gene absence or directed regulation.

For S(tau)={e:|Delta[e]|>tau}, define

    BFDR(tau) = sum_{e in S(tau)} ell[e] / |S(tau)|.

All K*G*(G-1)/2 candidates are a single family, with at most one direction per
candidate. For any data-determined nonempty set S, conditional on the data and
the baseline posterior,

    E[FDP(S)|X] = (1/|S|) * sum_{e in S} Pr(reported direction is false|X).

This is linearity of conditional expectation; edge independence is unnecessary.
It applies to data-selected thresholds because the selection is fixed once X
is conditioned on. It concerns the selected pairing and model allocation
states; it does not provide posterior uncertainty for biological correspondence.
Monte Carlo estimates, rather than exact probabilities, are used in code.

Evaluate every distinct absolute-contrast cutoff with exact score ties kept
together. Among sets with estimated mean error <= alpha, take the largest.
The average error can be nonmonotone in this contrast ordering. Do not stop
at the first failure. The canonical strict cutoff is the next excluded score,
or zero if all nonzero directions are included. Infinity means no discoveries.
Below the retained-draw minimum, no set is selected even if a candidate curve
passes alpha. Recompute the rule if defining a different inferential family.

The previous error-sorted prefix remains available as posterior_error.
all_matchings replaces ell by the worst error over all control counterparts;
it is a conservative sensitivity rule, not a posterior over mappings.
Neither rule is claimed to give uniform frequentist FDR control. Empirical
calibration requires repeated simulations with an explicitly specified truth.

Reference for the posterior expected error principle: Mueller, Parmigiani and
Rice, FDR and Bayesian Multiple Comparisons Rules,
https://web.ma.utexas.edu/users/pmueller/pap/MPR07.pdf .
