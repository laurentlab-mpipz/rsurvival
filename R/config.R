# Tell roxygen to document c++ functions
#' @useDynLib rsurvival
#' @importFrom Rcpp sourceCpp
NULL

# Unload the DLL when package is unloaded
.onUnload <- function (libpath) {
  library.dynam.unload("rsurvival", libpath)
}

"%dopar%" <- foreach::"%dopar%"
"%do%" <- foreach::"%do%"

# Print when pocakge loads
.onLoad <- function(...) {
  packageStartupMessage("Welcome to rsurvival, enjoy your stay!")
  packageStartupMessage("Can't believe it works...")
  packageStartupMessage("Have you ever noticed Samuel is awesome?")
}