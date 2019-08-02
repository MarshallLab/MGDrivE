.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Loading MGDrivE: Mosquito Gene Drive Explorer")
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
