## Resubmission

This is a revised submission of scNPM 1.0.6 addressing the review of version 1.0.0 by Leonore Hochhauser.

* The Title is "Single-Cell Network Propensity Matching" (39 characters).
* Title and Description have been checked for software/package/API names; no unquoted software names remain in these fields.
* All help examples are self-contained and unwrapped. The three executable help examples took 0.08, 0.01 and 0.11 seconds on the test machine. There are no dontrun or donttest wrappers.
* run_scNPM() is quiet by default (verbose = FALSE). With verbose = TRUE, R and native sampler progress use message conditions and are suppressible by suppressMessages(). Warnings and errors are unaffected.
* Tests verify identical sampled results and random-number state with informational output enabled, disabled or suppressed.

## Test environment

Windows 11 x64, R 4.5.2, GCC 14.3.0; OpenMP/BLAS threads set to 1 for testing.

All 538 test assertions passed (44 test blocks).

## R CMD check results

R CMD check --as-cran --run-donttest scNPM_1.0.6.tar.gz:
0 ERRORs, 0 WARNINGs, 3 NOTEs.

* Future file timestamps: "unable to verify current time"; the remote time check was unavailable.
* Top-level files: checking README.md and NEWS.md was skipped because pandoc is not installed locally.
* HTML manual: checking mathematical rendering was skipped because the V8 package is not installed locally. The PDF manual check passed.

The check also reported timeouts downloading the CRAN package index. Package installation, compiled code, R code, documentation, examples and tests passed. These results are from the available local Windows R 4.5.2 installation; no R-devel or other-platform check is claimed.
