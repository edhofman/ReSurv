#'@import reticulate

if(!("pyresurv"%in%virtualenv_list())){
virtualenv_create(envname = "pyresurv",packages = c("scipy",
                                                 "numpy",
                                                 "torch",
                                                 "torchtuples",
                                                 "shap"))}
use_virtualenv("pyresurv")

# global reference to scipy (will be initialized in .onLoad)
torch <- NULL
torchtuple <- NULL

.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to scipy
  torch <<- reticulate::import("torch", delay_load = TRUE)
  torchtuples <<- reticulate::import('torchtuples', delay_load = TRUE)
  shap <<- reticulate::import('shap', delay_load = TRUE)
}


