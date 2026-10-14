## Update: scNPM 1.0.8

This patch adds a damped Newton fallback for nonconvergent pairwise ridge
logistic fits. The statistical objective and existing convergent path are
unchanged.

## Test environment

* Windows, R 4.5.2, Rtools45
* R CMD check --no-manual --no-build-vignettes: Status OK
* Tests: 576 passed; 0 failures, warnings or skips
* Source build and Windows binary installation succeeded
* Separate Windows binary installation and solver smoke test passed

The PDF manual was not built in this local check.
