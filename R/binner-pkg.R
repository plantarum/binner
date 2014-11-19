##' Import and process ABI data files for AFLP analyses
##' 
##' \tabular{ll}{
##' Package: \tab binner\cr
##' Type: \tab Package\cr
##' Version: \tab 0.1\cr
##' Date: \tab 2013-10-09\cr
##' License: \tab GPL (>= 3)\cr
##' LazyLoad: \tab yes\cr
##' }
##'
##' Provides functions for reading fsa files 
##' 
##' @name binner
##' @aliases binner
##' @docType package
##' @title Import and process ABI data files for AFLP analyses
##' @author Tyler Smith \email{tyler@@plantarum.ca}
##' @keywords package aflp abi fsa
##' @import abind seqinr caTools
NA

binnerVersion <- 0.1

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("binner version ", binnerVersion)
  packageStartupMessage("See ?readFsa to get started")
}
