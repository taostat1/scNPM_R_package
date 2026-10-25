# Manuscript simulation study

This directory contains the data generator, method adapters, evaluation and
summary code used for the manuscript's 3000- and 1500-cell experiments.
Cell counts are **per condition**: each experiment contains 6000 or 3000 cells
in total. Each size has ten replicates and eight methods (160 fits).

## Design and seeds

Each condition has 100 genes and three populations with proportions 0.4, 0.3
and 0.3. The generator constructs positive-definite precision matrices with
shared and population-specific edges. Each matched pair has 11 gained and
11 lost edges, giving 22 differing edges relative to 110 baseline edges
(20%). Latent Gaussian expression
is transformed and subjected to expression-dependent dropout, targeting 40%
zeros. The PPI reference samples 60% of the union of all six true networks;
10% of its selected edges are drawn from the union's complement.

The 1500-cell observations are generated anew using the same network, mean
and dropout parameters as their corresponding 3000-cell replicate.
[seeds.csv](seeds.csv) lists every seed, including expression, PPI and fitting
seeds. [R/config.R](R/config.R) defines the design and all tuning settings.

## Dependencies

The analysis uses R 4.5.2, scNPM 1.1.1, BLGGM, glasso, JGL and jsonlite.
[R_package_versions.csv](R_package_versions.csv) records the local validation
environment. Install this repository's R package first, using the installation
instructions in the main README. Source builds require a C++17 toolchain.

From the repository root, install the baseline dependencies with:

```sh
Rscript inst/paper_simulation/install_dependencies.R
```

BLGGM is installed from the authors' repository at commit
`2acabd4c583b49801460865ef7ccb8b944edc7fa`. The installer changes its build
standard to C++14 for compatibility with the Rcpp headers; its statistical
source is unchanged.

AJGM uses the authors' Python implementation and proxTV. On Linux or macOS,
create a Python environment, install the dependencies below, and install
[AJGM](https://github.com/yyytim/AJGM) following its repository instructions.
The authors' source is obtained directly from their repository.

```sh
python3 -m venv ajgm-env
. ajgm-env/bin/activate
pip install numpy==2.5.3 scipy==1.18.1 scikit-learn==1.9.1 prox-tv==3.3.0
git clone https://github.com/yyytim/AJGM.git
pip install ./AJGM
python inst/paper_simulation/python/verify_upstream.py AJGM
export SCNPM_AJGM_PYTHON="$PWD/ajgm-env/bin/python"
```

proxTV source builds require BLAS/LAPACK and compiler development libraries.
On Windows, run AJGM in WSL and set these environment variables before Rscript:

```powershell
$env:SCNPM_AJGM_WSL = 'Ubuntu'
$env:SCNPM_AJGM_PYTHON = '/absolute/path/to/ajgm-env/bin/python'
# Optional: choose a WSL user instead of the distribution default.
# $env:SCNPM_AJGM_USER = 'your-user'
```

## Run the study

Run these commands from the downloaded repository root. The output directory
is created automatically. Fits run serially with eight threads by default;
`--threads` controls the requested thread count.

```sh
# Generate all twenty datasets.
Rscript inst/paper_simulation/run.R --phase prepare --n both --output paper_simulation_output

# Fit all methods, evaluate each fit, and summarize both cell sizes.
Rscript inst/paper_simulation/run.R --phase run --n both --threads 8 --output paper_simulation_output

# Rebuild summaries from completed fits.
Rscript inst/paper_simulation/run.R --phase summarize --output paper_simulation_output
```

Select a size, replicate or method when running individual jobs:

```sh
Rscript inst/paper_simulation/run.R --phase run --n 3000 --replicates 1 --methods scNPM,scNPM_independent_fixedNPM --threads 8 --output paper_simulation_output
```

Completed fits are reused on restart. The independent-prior method follows
the corresponding scNPM fit; the runner verifies identical adaptation,
frozen matching and freeze iteration. Failed jobs record their error and stop
the runner. Rerunning the command retries the incomplete job.

## Methods and output conversion

| Method identifier | Fitting and selection |
| --- | --- |
| `scNPM` | Staged NPM, frozen four-state paired prior, union PPI weights; BFDR 0.05 |
| `scNPM_noPPI` | Same paired model with uniform base edge prior |
| `scNPM_independent_fixedNPM` | Same adaptation and frozen matching, independent network prior |
| `BLGGM` | Separate condition fits; posterior-product differential probabilities and BFDR 0.05 |
| `glasso` | PCA/k-means populations, graphical lasso, EBIC tuning |
| `mixtureGGM` | Penalized mixture EM, three starts, EBIC tuning |
| `AJGM` | Authors' algorithm, 100-iteration cap, absolute partial-correlation cutoff 0.05 |
| `JGL_fused` | Matched populations, fused graphical lasso, EBIC tuning |

scNPM and BLGGM use 5000 iterations and 1000 saved draws. AJGM uses the
100-iteration settings of the final tables. Its stopping diagnostic records
the attained objective change and whether the 0.1 tolerance was reached;
the iteration cap and convergence status are separate fields. The reference
table includes the 100-iteration estimates. Tuning grids, preprocessing,
precision-to-score conversion and edge selection are implemented in
[R/baselines.R](R/baselines.R) and [R/methods.R](R/methods.R).

## Metrics and summaries

[R/evaluation.R](R/evaluation.R) calculates matching accuracy, precision,
recall, F1, AP, false-discovery proportion and edge counts. Differential
evaluation requires both the correct gene pair and direction. Overall
precision/recall/F1 pool the six networks; overall mAP averages their six
AP values. Differential mAP averages AP over the three matched pairs.
AP uses grouped score ties and a stepwise precision-recall integral.

Each fit writes `metrics.csv`, `overall_by_network.csv`,
`differential_edges.csv`, `fit_summary.rds`, settings and session information.
scNPM also writes `iterations.csv`; `progress.log` records job completion.
The summary contains replicate means and standard errors, the number of
finite values contributing to each metric, and paired differences against
scNPM. Metrics with no defined precision remain missing.

The small [reference_results](reference_results) tables contain the actual
per-replicate manuscript results. Recompute their summary without fitting:

```sh
Rscript inst/paper_simulation/run.R --phase reference --output paper_reference_summary
```

This command summarizes the supplied results; `--phase run` generates and
fits new results. [PROVENANCE.md](PROVENANCE.md) describes the original scripts
and the data/function/metric equivalence checks.

## Method sources

- [scNPM](https://github.com/taostat1/scNPM_R_package)
- [BLGGM implementation](https://github.com/WgitU/BLGGM) and
  [Wu and Luo, Annals of Applied Statistics](https://doi.org/10.1214/21-AOAS1582)
- [AJGM implementation and publication information](https://github.com/yyytim/AJGM)
- [glasso package](https://cran.r-project.org/package=glasso)
- [JGL package](https://cran.r-project.org/package=JGL)
- [proxTV](https://github.com/albarji/proxTV)
