.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Loading MGDrivE: Mosquito Gene Drive Explorer")
}

.onUnload <- function (libpath) {
  library.dynam.unload("MGDrivE", libpath)
}

#' @importFrom utils globalVariables

# CRAN Note avoidance
if(getRversion() >= "2.15.1"){
  utils::globalVariables(
    c("private", # network object access
      "self"     # patch/network function call
    )
  )
}
