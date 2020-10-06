.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Loading MGDrivE2: Mosquito Gene Drive Explorer Version 2")
}

#' @importFrom utils globalVariables

# CRAN Note avoidance
#  all of these variables come from set_populations(), which is used in
#  equilibrium_SEI_SIS() and equilibrium_SEI_SEIR()
#  it makes use of cute scoping tricks to reuse code and supply the correct
#  variables necessary
if(getRversion() >= "2.15.1"){
  utils::globalVariables(
    c("fSEI",
      "mosyHList",
      "hPopAq",
      "hPopF",
      "hPopM",
      "mosyList",
      "mPopAq",
      "mPopF",
      "mPopM"
    )
  )
}



