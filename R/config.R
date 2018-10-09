# Tell roxygen to document c++ functions
#' @useDynLib rsurvival
#' @importFrom Rcpp sourceCpp
NULL

# Unload the DLL when package is unloaded
.onUnload <- function (libpath) {
  library.dynam.unload("rsurvival", libpath)
}