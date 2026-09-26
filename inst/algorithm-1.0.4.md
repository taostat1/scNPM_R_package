# scNPM 1.0.4: staged inference and directional error selection

## 1. State, stages and label convention

The state comprises latent expression theta, component memberships C, mixture
weights pi, means mu, precisions Omega, edge allocations z and observation
parameters Lambda. Each condition has K components. Initial component means
are fixed labeling anchors. At every complete sweep, a minimum-cost bijection
aligns current means to these anchors, using squared Euclidean distance. The
same permutation acts on all component-specific parameters and memberships.
The implementation requires exchangeable mixture priors. Anchors name
components; they do not establish biological subtype identities.

During adaptation the NPM assignment phi is recomputed after each sweep. All
K clusters in both conditions must be nonempty. A complete assignment identical
for `freeze_patience` consecutive sweeps, after `adapt_min`, is frozen. On the
next sweep, the counterpart adjustment multiplier becomes one. A run stops if
stability has not been reached by `adapt_max`; it does not silently freeze the
last assignment. A stable correspondence need not be a correct correspondence.

Let L = num_iter, S = num_save and F be the 1-based freeze iteration. Retained
iterations are L-S+1,...,L; the post-switch discarded interval is F+1,...,L-S.
Thus actual final burn-in is L-S-F, checked to be at least final_burnin.
Neither adaptation nor post-switch burn-in contributes to reported samples.
Final matching diagnostics use the frozen bijection, not another greedy match.

The final chain has baseline PPI weights xi_jt and no counterpart feedback.
Its intended target is the baseline augmented posterior modulo the deterministic
label convention. With valid component kernels, a proper target and adequate
ergodic sampling, averages estimate expectations under this target. Finite
burn-in and minimum draw counts do not prove convergence. In particular, the
invariant target no longer depends on the frozen map through the edge prior.
This is a substantive inferential change from the coupled 1.0.3 algorithm.

## 2. Joint edge-state probabilities

For a retained draw s and precision entry omega, the conditional probability is

p(omega,xi) = xi * Normal(omega;0,v1^2) /
             [xi * Normal(omega;0,v1^2) + (1-xi) * Normal(omega;0,v0^2)].

The implementation evaluates the log odds, including both prior weights.
Given precisions, final-stage edge allocations in the two conditions have
independent Bernoulli conditionals under the baseline joint prior. The SPD
constraint is fixed when conditioning on Omega. For treatment cluster k and
control cluster c, set pt_s = p(Omega_T[j,t,k],xi_jt) and similarly pc_s.
Then estimate

p10(k,c) = mean_s pt_s*(1-pc_s),
p01(k,c) = mean_s (1-pt_s)*pc_s,
p00(k,c) = mean_s (1-pt_s)*(1-pc_s),
p11(k,c) = mean_s pt_s*pc_s.

The product is taken WITHIN the same draw; shared observation parameters and
mixture uncertainty can induce marginal dependence. Multiplying separate
marginal averages is not generally equivalent. Raw paired allocations and
contiguous batch averages are returned for numerical diagnostics.

Array slice (k-1)*K+c uses R's 1-based indices. Diagonal entries are zero
placeholders, not probabilities; off-diagonal entries are symmetric and sum
to one across the four states. K=4 produces sixteen pair slices.

## 3. Directional selection and proof

One item h=(k,j,t), j<t, can yield at most one directional claim. Choose d_h as
the more probable of states (1,0) and (0,1) under frozen counterpart phi(k).
Equal probabilities give no directional candidate. Let A_h(c) be the event that
this chosen direction is correct against control c. For exact posterior
probabilities define

ell_h(c) = P(A_h(c)^c | X),
ell_h,frozen = ell_h(phi(k)),
ell_h,robust = max_c ell_h(c).

For any fixed bijection psi, ell_h(psi(k)) <= ell_h,robust. Let S(X) be a
selected set of items and R=|S|. Conditional on the observed data, S and all
reported directions are fixed. Linearity of expectation gives

E[FDP_directional(S,psi) | X]
 = sum_{h in S} P(A_h(psi(k))^c | X) / max(1,R)
 <= sum_{h in S} ell_h,robust / max(1,R).

For frozen scope the same calculation is an equality with ell_h,frozen and
psi=phi. Independence across candidate edges is not required. This argument
bounds every fixed correspondence using the same baseline posterior; it is
not an integration over a random, jointly modeled correspondence. No posterior
distribution for phi is fitted in this version.

Sort the chosen errors ell_(1)<=...<=ell_(M), and take the largest r for which
(sum_{m=1}^r ell_(m))/r <= alpha. If no such r exists, take the empty set and FDP
zero. This yields an exact posterior expected directional FDP bound alpha when
the probabilities are exact. Sorting minimizes the sum of errors at each
fixed cardinality, so this prefix has the largest feasible cardinality for the
chosen candidate directions. All K*G*(G-1)/2 items share this one family.

In code, probabilities are Monte Carlo estimates; the computed prefix mean is
an estimate, not an exact bound on the unknown expectations. At least 200
retained draws are required by default (configurable, lower limit 20). This
protects against accidental tiny runs, not poor mixing. Batch MCSE estimates
sampling variability without simultaneous confidence coverage and is not used
as an upper confidence bound in the selection rule. The max over counterparts
can also be noisy; examine independent chains and longer runs.

The theoretical posterior statement concerns continuous spike/slab allocations.
It is not uniform frequentist FDR control for every fixed precision matrix or
for exact nonzero precision entries. In particular, under a graph-equality
null, one group's uncertain edge allocation can differ from the other group's
allocation. Cell-level independence is also a modeling assumption: biological
replicate dependence requires an appropriate observation model. The separate
simulation report includes exact-zero graph recovery as a stress test rather
than treating it as a proof of the posterior statement.

## 4. Meaning of the criterion trace

The proposed observed-data criterion requires integrating missing latent
expression and summing memberships. The native trace instead records

J_aug = log p(X,theta_missing,C | pi,mu,Omega,Lambda)
        + network_prior_kernel(z,Omega;current_weights)
        - sum_k SMD(k,phi(k)).

The likelihood includes Gaussian densities, mixing probabilities, the probit
zero/positive observation probabilities and the positive-observation Jacobian
-log X. The network factor includes spike/slab normal densities, current
allocation weights and exponential diagonal terms. It omits the other priors
and common normalization constants. During adaptation, adjusted edge factors
are state-wise diagnostics; they are not asserted to be compatible joint full
conditionals. After feedback is disabled, the factor uses baseline xi.

Missing/undefined SMD produces NA loss and J_aug, preserving the likelihood and
prior columns. Diagnostic evaluation does not consume random numbers, change
updates, or decide Metropolis acceptance. Acceptance correction uses the
corresponding HMC Hamiltonian. The criterion can decrease even in stationarity;
changes at the stage boundary must not be read as objective improvement.
Greedy NPM makes sequential assignments and need not minimize total SMD.

The defensible wording is: "The components are updated conditional on the
current state; the criterion records their changing contributions." No
monotonic likelihood, prior or common-objective ascent is claimed.

## 5. Compatibility, validation and references

`inference_mode='coupled'` retains 1.0.3 update order, RNG consumption and
continuous feedback; it reports descriptive frequency differences. Posterior
error selection refuses those coupled outputs. Version 1.0.4 keeps 1.0.3's
2K reference covariates and pooled column scaling in both modes.

Tests cover stage boundaries, failure to stabilize, fixed correspondence,
label assignment, per-draw joint probabilities, all sixteen K=4 comparisons,
a hand-computed prefix, and an alternative-matching counterexample. The release
folder includes the full regression suite, package check and repeated simulations.

The posterior expected error sorting principle is standard Bayesian multiple
testing, rather than a newly claimed theorem of scNPM. References:

- Mueller, Parmigiani and Rice (2007), FDR and Bayesian multiple comparisons:
  https://web.ma.utexas.edu/users/pmueller/pap/MPR07.pdf
- Peterson, Stingo and Vannucci (2015), Bayesian inference of multiple Gaussian
  graphical models, JASA: https://pmc.ncbi.nlm.nih.gov/articles/PMC4465207/

The joint-probability estimator, fixed-map versus worst-counterpart comparison
and stage bookkeeping above document their implementation in this package.
