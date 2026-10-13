## Update: scNPM 1.0.7

This update introduces pairwise propensity matching, consensus freezing at the
adaptation limit, and optional elapsed-time progress and CSV logging.
The earlier global matching and strict freezing are available through arguments.
Examples remain self-contained; console progress is quiet by default and uses
suppressible message conditions.

## Test environment

Windows 11 x64, R 4.5.2 UCRT, Rtools45 (GCC 14.3.0).

## Check results

R CMD check --no-manual --no-build-vignettes scNPM_1.0.7.tar.gz:
0 errors, 0 warnings, 0 notes (Status: OK).
All examples and testthat regression tests passed.

The Windows binary was also installed into a separate library and exercised
with a small sampling run. Package version, consensus freezing, retained draws,
and progress CSV output were verified.
