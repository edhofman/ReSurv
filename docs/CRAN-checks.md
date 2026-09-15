# CRAN preparation

The workflows follow the structure of the local `clmplus` package and
use
[r-lib/actions](https://github.com/r-lib/actions/tree/v2/check-r-package).

`r-checkrelease.yml` runs on pull requests, pushes to main/master, and
manual dispatch. It provides:

- A documentation audit pinned to the version in `RoxygenNote`, unit
  tests, URL checks, and explicit COX/XGB/NN smoke checks. The smoke
  checks install and verify the native torch runtime and fail on
  fitting/prediction errors.
- R release checks on Windows and macOS, and R devel/oldrel-1 on Linux.
  These fail on warnings and retain check logs as artifacts.
- A separate Linux release check with TinyTeX, HTML Tidy, vignettes, and
  the PDF manual. This fails on notes as well as warnings and errors.

`pkgdown.yaml` builds the website on pull requests and main/master
pushes. Only a successful main/master build can deploy to `gh-pages`;
the build job has read-only repository permissions. Configure GitHub
Pages to serve the root of `gh-pages` if it does not already do so.
Manual runs on feature branches build an artifact without publishing it.

`rhub.yaml` retains manual R-hub checks; the obsolete Python
installation steps have been removed.

## Local validation on 2026-09-14

Windows 11, R 4.5.2:

- Implementation changes are outside this task’s scope. Existing R
  function bodies and `NAMESPACE` are unchanged; only roxygen
  documentation was edited in existing R files.
- The initial documentation-only check completed with **0 errors, 1
  warning, 3 notes**. Existing findings include an obsolete `Om.df`
  argument in an internal helper, missing namespace imports,
  unregistered data-mask column names, unused Imports declarations, and
  inability to verify current time.
- A subsequent full check passed the PDF/HTML manuals, examples, tests,
  and vignette rebuilding. Its code fixes were reverted to retain the
  requested documentation/workflow-only scope; its cleaner result is not
  a validation of the final unchanged implementation.
- Unit tests: 72 passed, 0 failed, 0 warnings, 2 skipped. The older
  small-data XGB/NN tests catch fitting errors and skip; the separate
  backend smoke checks use simulated data and fail on errors instead.
- Explicit COX, XGB, and native torch smoke checks passed.
- Re-running roxygen2 7.3.3 reproduced `NAMESPACE` and `man/` exactly.
- All workflow YAML files parsed; the pkgdown website built
  successfully.
- All package documentation URLs passed
  [`urlchecker::url_check()`](https://rdrr.io/pkg/urlchecker/man/url_check.html)
  after updating the redirected neural-network reference.

The strict CI jobs will report the existing implementation findings and
are not expected to be green until those are addressed in a separate
code change. The package currently sets `_R_CHECK_SYSTEM_CLOCK_` when
loaded; this existing behavior is also unchanged and should be reviewed
before submission.

The local check uses devtools defaults, including disabled incoming
checks and optional Suggests enforcement. It does not replace the
stricter GitHub release job or CRAN’s own checks. The Windows
environment required `LC_ALL`, `LC_CTYPE`, and `LANG` to be set to `C`;
its inherited `C.UTF-8` locale was invalid. A local Quarto
command-wrapper warning was emitted after the check; it did not affect
the package’s R Markdown vignette or manual results.

## Before submission

1.  Push the changes and require successful package-audit, platform, and
    full release jobs. Download and review the check artifacts,
    including skipped tests and any notes.
2.  Verify the website deployment and documentation links.
3.  Confirm maintainer details and update the submission comments with
    the actual GitHub/CRAN preflight results. Follow the current [CRAN
    policies](https://cran.r-project.org/web/packages/policies.html).

No CRAN submission or remote workflow execution was performed during the
local preparation. Historical research sources remain in
`articles/historical/` and are excluded from the source package. The
current executable vignette is `vignettes/getting-started.Rmd`.
