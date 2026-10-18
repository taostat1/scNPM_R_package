## Update: scNPM 1.1.0

Addresses the reported 1.0.6 macOS/Clang test-staged-fdr.R failure. The old
test expected a random multi-cluster chain not to stabilize, an assumption
that depended on the numerical trajectory. The replacement retains the
multi-cluster/global/greedy path but uses freeze_patience > adapt_max, which
deterministically exercises the strict-stability error.
Test fixtures use one native thread to respect constrained check environments.

This release also adds the documented four-state matched-network prior and
regression tests for its conditional probabilities, sampling and labels.
The previous network prior remains available with network_prior="independent".
