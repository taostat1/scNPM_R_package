# Code and result provenance

The statistical functions were extracted from the production simulation
scripts used in the manuscript. File handling and process launch code were
consolidated into `run.R` and `R/portable.R` so paths and output locations can
be specified by the reader.

Original project-relative scripts:

- `Simulation/scNPM_110_Cunion_N3000_N1500_20261020/RUN_COMPARISON.R`:
  shared generation, baseline fitting and evaluation functions.
- `Simulation/scNPM_111_N3000_publicPPI_20261020/RUN_EXPERIMENT.R`:
  final 3000-cell scNPM settings and six-network union PPI.
- `Simulation/scNPM_111_N1500_publicPPI_20261020/RUN_EXPERIMENT.R`:
  final 1500-cell settings and fresh expression generation.
- `Simulation/scNPM_110_N3000_Tunion_20261018/RUN_COMPARISON.R` and
  `Simulation/scNPM_110_N1500_Tunion_20261020`: baseline fits reused in the
  final comparison; these five methods do not consume the PPI reference.
- `outputs/simulation_review_20261022/audit_results.R`: final result audit
  and tables. Reference CSVs retain its values, with local source paths removed.

## Validation

- [Data equivalence](validation/data_equivalence.csv): all 20 regenerated
  datasets exactly matched the archived inputs, including both expression
  matrices, PPI, labels, latent expression and generating parameters.
- [Function audit](validation/function_audit.csv): extracted functions were
  compared as parsed R expressions; baseline/evaluation functions were also
  checked across the original experiment scripts.
- [Metric equivalence](validation/metric_equivalence.csv): evaluation of all
  160 archived fits reproduced their original network and differential metrics
  to a tolerance of 1e-12.

The release validation regenerated data and reevaluated existing fits. Full
5000-iteration fits are produced by the supplied `run` command.

The command-line data preparation was exercised for both cell sizes. The
reference-summary command reproduced all 32 summary rows and 112 paired
comparisons to a tolerance of 1e-12. Short integration runs checked the three
scNPM adapters, fixed-matching replay, and the AJGM WSL bridge with stopping
diagnostics. These short runs are separate from the manuscript estimates.

The AJGM adapter preserves the production preprocessing, seed, optimization
arguments and output conversion. The convergence recorder observes the
upstream stopping state without changing the returned estimate. Upstream
Python source fingerprints are supplied in
`python/upstream_AJGM_sha256.json`; line endings are normalized before hashing.
BLGGM and AJGM are installed from their authors' repositories; their source
and license information remain with those repositories.
