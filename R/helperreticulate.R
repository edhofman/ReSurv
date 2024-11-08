#'@import reticulate

if(!("pyresurv"%in%virtualenv_list())){


  packages_list <- c("pip",
                     "scipy",
                     "numpy",
                     "torch",
                     "torchtuples",
                     "shap")

  if(grepl("ubuntu",tolower(Sys.info()[["sysname"]])) | grepl("debian",tolower(Sys.info()[["sysname"]])) ){
    packages_list <- c("venv",packages_list)
  }

virtualenv_create(envname = "pyresurv",
                  packages = packages_list,
                  force = TRUE)}
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


