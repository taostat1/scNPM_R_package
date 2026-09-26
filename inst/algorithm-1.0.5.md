# scNPM 1.0.5: complete global matching and calibrated contrast thresholds

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
SMDs are not discarded. Stable nonempty complete maps freeze as in version
1.0.4; feedback is disabled in subsequent iterations and final warmup precedes
retention. Label alignment still uses minimum squared distance to fixed anchors,
which is a separate assignment problem from matching T to C clusters.

The optimizer's guarantee refers to each current cost matrix. A frozen map is
not reoptimized at the final state. Its final maximum SMD need not remain
optimal, nor must the maximum SMD decrease from one iteration to the next.

## Threshold selection

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
The earlier algorithm-1.0.4.md records that version's defaults and sampling
mechanism; the matching and threshold defaults above supersede it in 1.0.5.
