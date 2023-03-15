# global reference to scipy (will be initialized in .onLoad)
torch <- NULL
torchtuple <- NULL

.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to scipy
  scipy <<- reticulate::import("torch", delay_load = TRUE)
  hdbscan <<- reticulate::import('torchtuples', delay_load = TRUE)
}
