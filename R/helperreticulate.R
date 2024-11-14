#' Install Python Environment for ReSurv
#'
#' Install a Python environment that allows the user to apply the Neural Network (NN) models.
#'
#' @param ... Additional arguments for `virtualenv_create`.
#' @param envname `character`. Name of the environment created. Default `pyresurv`.
#' @param new_env `logical`. If `TRUE`, any existing Python virtual environment and/or `conda` environment specified by `envname` is deleted first.
#'
#' @return No return value.
#'
#' @import reticulate
#' @export
install_pyresurv <- function(...,
                             envname = "pyresurv",
                             new_env = identical(envname, "pyresurv")) {

  if(new_env && virtualenv_exists(envname)){

    virtualenv_remove(envname)


    }

  packages_list <- c("pip",
                     "scipy",
                     "numpy",
                     "torch",
                     "torchtuples",
                     "shap")


  virtualenv_create(envname = "pyresurv",
                    packages = packages_list,
                    force = TRUE,
                    ...)



}


.onLoad <- function(...) {
  use_virtualenv("pyresurv", required = FALSE)
}


# global reference to scipy (will be initialized in .onLoad)
torch <- NULL
torchtuple <- NULL
shap <- NULL

.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to scipy
  torch <<- reticulate::import("torch", delay_load = TRUE)
  torchtuples <<- reticulate::import('torchtuples', delay_load = TRUE)
  shap <<- reticulate::import('shap', delay_load = TRUE)
}
