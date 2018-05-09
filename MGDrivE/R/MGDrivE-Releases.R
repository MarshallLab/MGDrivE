########################################################################
#       __  _____________       _       ______
#      /  |/  / ____/ __ \_____(_)   __/ ____/
#     / /|_/ / / __/ / / / ___/ / | / / __/
#    / /  / / /_/ / /_/ / /  / /| |/ / /___
#   /_/  /_/\____/_____/_/  /_/ |___/_____/
#
#   MGDrivE Releases
#   Marshall Lab
#   November 2017
#
########################################################################

#' Make List of Modified Mosquito Releases
#'
#' Sets up a release schedule for a single patch, returns a list to be used in
#' \code{\link{oneDay_maleReleases_Patch}} or \code{\link{oneDay_femaleReleases_Patch}}.
#'
#' @param genotypes possible genotypes
#' @param releaseStart day releases start
#' @param releaseEnd day releases end
#' @param releaseInterval interval between releases
#' @param releaseVector named character vector of release composition
#' @param sex character in 'M','F'
#'
#' @examples
#' # to setup for 3 patches but only release in the first with a defined release schedule:
#'
#' patchReleases = replicate(n = 3,expr = {
#'   list(maleReleases = NULL,femaleReleases = NULL)
#' },simplify = FALSE)
#'
#' patchReleases[[1]]$femaleReleases = Release_basicRepeatedReleases(genotypes = Cube_Homing1RA$genotypesID,releaseStart = 5,releaseEnd = 30,releaseInterval = 5,releaseVector = c("HH"=100,"Hh"=0,"HR"=0,"hh"=0,"hR"=0,"RR"=0),sex = "F")
#' patchReleases[[1]]$maleReleases = Release_basicRepeatedReleases(genotypes = Cube_Homing1RA$genotypesID,releaseStart = 50,releaseEnd = 60,releaseInterval = 1,releaseVector = c("HH"=100,"Hh"=0,"HR"=0,"hh"=0,"hR"=0,"RR"=0),sex = "M")
#'
#' @export
Release_basicRepeatedReleases <- function(genotypes, releaseStart, releaseEnd, releaseInterval, releaseVector, sex="M"){

  # check timing of releases
  if(releaseInterval > (releaseEnd - releaseStart)){
    stop("interval between releases cannot be greater than time between start and end of releases")
  }

  # name and check releaseVector. Initialize release times. Initialize return list
  releaseTimes = seq(from=releaseStart,to = releaseEnd,by = floor(releaseInterval))
  releaseList = vector(mode="list",length=length(releaseTimes))

  # check for male or female. Fill appropriate list.
  if(sex=="M"){
    for(tx in 1:length(releaseTimes)){
      releaseList[[tx]]$nuM = releaseVector
      releaseList[[tx]]$tRelease = releaseTimes[[tx]]
    }
  } else if(sex=="F"){
    for(tx in 1:length(releaseTimes)){
      releaseList[[tx]]$nuF = releaseVector
      releaseList[[tx]]$tRelease = releaseTimes[[tx]]
    }
  } else if(sex=="L"){
    for(tx in 1:length(releaseTimes)){
      releaseList[[tx]]$larvae = releaseVector
      releaseList[[tx]]$tRelease = releaseTimes[[tx]]
    }
  } else {
    stop(paste0("expected character in 'M','F','L' in argument 'sex', got: ",sex))
  }
  return(releaseList)
}

#' Make List of Modified Mosquito Releases
#'
#' Sets up a release schedule for a single patch, calls \code{\link{Release_basicRepeatedReleases}} internally.
#'
#' @param driveCube gene-drive cube
#' @param releasesParameters A list containing the releasesStart, releasesNumber
#' releasesInterval, and releaseProportion named values.
#' @param sex character in 'M','F'
#'
#' @export
generateReleaseVector=function(driveCube=driveCube,releasesParameters=releasesParameters,sex="M"){
  genotypes=driveCube$genotypesID
  releaseVector=rep(0,length=length(genotypes))
  names(releaseVector)=genotypes
  releaseVector[driveCube$releaseType]=releasesParameters$releaseProportion
  Release_basicRepeatedReleases(
    genotypes=genotypes,
    releaseStart=releasesParameters$releasesStart,
    releaseEnd=releasesParameters$releasesStart+releasesParameters$releasesInterval*(releasesParameters$releasesNumber-1),
    releaseInterval=releasesParameters$releasesInterval,
    releaseVector=releaseVector,
    sex=sex
  )
}
