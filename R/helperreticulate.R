#' Install Python Environment for ReSurv
#'
#' Install a Python environment for optional Python-based features (e.g., SHAP for NN models).
#' The core NN backend now uses native R torch and does not require Python.
#'
#' @param ... Additional arguments for `virtualenv_create`.
#' @param envname `character`. Name of the environment created. Default `pyresurv`.
#' @param new_env `logical`. If `TRUE`, any existing Python virtual environment and/or `conda` environment specified by `envname` is deleted first.
#'
#' @return No return value.
#'
#' @export
install_pyresurv <- function(...,
                             envname = "pyresurv",
                             new_env = identical(envname, "pyresurv")) {

  if(!requireNamespace("reticulate", quietly = TRUE))
    stop("Package 'reticulate' is required for install_pyresurv(). Install it with install.packages('reticulate').")

  if(new_env && reticulate::virtualenv_exists(envname)){

    reticulate::virtualenv_remove(envname)


    }

  packages_list <- c("pip",
                     "scipy",
                     "numpy",
                     "torch",
                     "torchtuples",
                     "shap")


  reticulate::virtualenv_create(envname = "pyresurv",
                    packages = packages_list,
                    force = TRUE,
                    ...)



}

