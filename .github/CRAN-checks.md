# CRAN preparation

The workflows follow the structure of the local `clmplus` package and use
[r-lib/actions](https://github.com/r-lib/actions/tree/v2/check-r-package).

`r-checkrelease.yml` runs on pull requests, pushes to main/master, and manual
dispatch. It provides:

- A documentation audit pinned to the version in `RoxygenNote`, unit tests,
  URL checks, and explicit COX/XGB smoke checks that fail on fitting/prediction
  errors. Torch tests are skipped on CI and CRAN using `skip_on_ci()` and
  `skip_on_cran()`. Run `testthat::test_local()` locally to include torch tests,
  with the native torch runtime installed. For local `R CMD check`, set
  `NOT_CRAN=true` to include them. Run `Rscript scripts/smoke_backends.R` locally
  to smoke-test all three backends; CI passes `--skip-torch`.
- R release checks on Windows and macOS, and R devel/oldrel-1 on Linux.
  These fail on warnings and retain check logs as artifacts.
- A separate Linux release check with TinyTeX, HTML Tidy, vignettes, and the
  PDF manual. This fails on notes as well as warnings and errors.

`pkgdown.yaml` builds the website on pull requests and main/master pushes.
Only a successful main/master build can deploy to `gh-pages`; the build job
has read-only repository permissions. Configure GitHub Pages to serve the
root of `gh-pages` if it does not already do so. Manual runs on feature
branches build an artifact without publishing it.

`rhub.yaml` retains manual R-hub checks; the obsolete Python installation
steps have been removed.

## Fixes awaiting validation

The issues from the previous check have been addressed without changing the
model-fitting, simulation, or prediction formulas or their public defaults:

- Added missing namespace imports and removed the unused `tibble` dependency.
- Registered data.table/dplyr column names for static code analysis.
- Removed an obsolete `Om.df` argument passed to `hazard_data_frame()`.
- Made `time_unit` and `ref_claim` explicit required inputs of two unused
  internal inflation helpers; their formulas are unchanged.
- Removed the package-load side effect that disabled R's system-clock check.
- Marked three internal dotted-name helpers as non-exported helpers for roxygen.

No tests, model runs, vignette execution, or package checks were run after
these fixes, as requested. Earlier results below describe the pre-fix code;
they are not validation of these changes. Run the tests and a full package
check before submitting to CRAN.

## Previous local validation on 2026-09-14

Windows 11, R 4.5.2:

- The pre-fix `devtools::check(document = FALSE, manual = TRUE,
  error_on = "never")` completed with **0 errors, 1 warning, 2 notes**.
  Existing findings include an obsolete `Om.df` argument in an internal
  helper, missing namespace imports, unregistered data-mask column names,
  and unused Imports declarations.
- PDF/HTML manuals, examples, tests, and vignette rebuilding passed with the
  original implementation.
- Unit tests: 72 passed, 0 failed, 0 warnings, 2 skipped. The older small-data
  XGB/NN tests catch fitting errors and skip; the separate backend smoke
  checks use simulated data and fail on errors instead.
- Explicit COX, XGB, and native torch smoke checks passed.
- Re-running roxygen2 7.3.3 reproduced `NAMESPACE` and `man/` exactly.
- All workflow YAML files parsed; the pkgdown website built successfully.
- All package documentation URLs passed `urlchecker::url_check()` after
  updating the redirected neural-network reference.

The local check uses devtools defaults, including disabled incoming checks
and optional Suggests enforcement. It does not replace the stricter GitHub
release job or CRAN's own checks. The Windows environment required `LC_ALL`,
`LC_CTYPE`, and `LANG` to be set to `C`; its inherited `C.UTF-8` locale was
invalid. A local Quarto command-wrapper warning was emitted after the check;
it did not affect the package's R Markdown vignette or manual results.

## Before submission

1. Push the changes and require successful package-audit, platform, and full
   release jobs. Download and review the check artifacts, including skipped
   tests and any notes.
2. Verify the website deployment and documentation links.
3. Confirm maintainer details and update the submission comments with the
   actual GitHub/CRAN preflight results. Follow the current
   [CRAN policies](https://cran.r-project.org/web/packages/policies.html).

No CRAN submission or remote workflow execution was performed during the
local preparation. Historical research sources remain in `articles/historical/`
and are excluded from the source package. The current executable vignette is
`vignettes/getting-started.Rmd`.
