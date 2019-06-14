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
#' @param sex character in 'M','F','E'
#'
#' @examples
#' # to setup for 3 patches but only release in the first with a defined release schedule:
#'
#' patchReleases = replicate(n = 3,expr = {
#'   list(maleReleases = NULL,femaleReleases = NULL)
#' },simplify = FALSE)
#'
#' patchReleases[[1]]$femaleReleases = basicRepeatedReleases(genotypes = cubeHoming1RA$genotypesID,
#'                                                           releaseStart = 5,
#'                                                           releaseEnd = 30,
#'                                                           releaseInterval = 5,
#'                                                           releaseVector = c("HH"=100,
#'                                                                             "Hh"=0,
#'                                                                             "HR"=0,
#'                                                                             "hh"=0,
#'                                                                             "hR"=0,
#'                                                                             "RR"=0),
#'                                                           sex = "F")
#' patchReleases[[1]]$maleReleases = basicRepeatedReleases(genotypes = cubeHoming1RA$genotypesID,
#'                                                         releaseStart = 50,
#'                                                         releaseEnd = 60,
#'                                                         releaseInterval = 1,
#'                                                         releaseVector = c("HH"=100,
#'                                                                           "Hh"=0,
#'                                                                           "HR"=0,
#'                                                                           "hh"=0,
#'                                                                           "hR"=0,
#'                                                                           "RR"=0),
#'                                                         sex = "M")
#'
#' @export
basicRepeatedReleases <- function(genotypes, releaseStart, releaseEnd,
                                  releaseInterval, releaseVector, sex="M"){

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
  } else if(sex=="E"){
    for(tx in 1:length(releaseTimes)){
      releaseList[[tx]]$nuE = releaseVector
      releaseList[[tx]]$tRelease = releaseTimes[[tx]]
      releaseList[[tx]]$releaseFlag = FALSE # trip to true to release into larval window
    }
  } else {
    stop(paste0("expected character in 'M','F','L' in argument 'sex', got: ",sex))
  }
  return(releaseList)
}

#' Make List of Modified Mosquito Releases
#'
#' Sets up a release schedule for a single patch, calls \code{\link{basicRepeatedReleases}} internally.
#'
#' @param driveCube gene-drive cube
#' @param releasesParameters A list containing the releasesStart, releasesNumber
#' releasesInterval, and releaseProportion named values.
#' @param sex character in 'M','F'
#'
#' @examples
#' # setup a drive cube, using Mendelian as the example
#' cube <- cubeMendelian()
#'
#' # setup release parameter list
#' #  releasesStart is the time of first release
#' #  releasesNumber is the number of releases
#' #  releasesInterval is the number of days between releases
#' #  releaseProportion is the number of mosquitoes released
#' relParams <- list(releasesStart = 25, releasesNumber = 1,
#'                   releasesInterval = 0, releaseProportion = 10)
#'
#' # generate male releases
#' mRelVec <- generateReleaseVector(driveCube = cube,
#'                                  releasesParameters = relParams,
#'                                  sex = "M")
#'
#' # generate female releases
#' fRelVec <- generateReleaseVector(driveCube = cube,
#'                                  releasesParameters = relParams,
#'                                  sex = "F")
#'
#' @export
generateReleaseVector <- function(driveCube=driveCube,
                                  releasesParameters=releasesParameters,
                                  sex="M"){

  # edge cases
  if(releasesParameters$releasesNumber==0L){
    return(NULL)
  }

  # generate timing parameters
  if(releasesParameters$releasesNumber==1L){
    start = releasesParameters$releasesStart
    end = releasesParameters$releasesStart
    interval = 0
  } else {
    start = releasesParameters$releasesStart
    end = releasesParameters$releasesStart + releasesParameters$releasesInterval * (releasesParameters$releasesNumber-1)
    interval = releasesParameters$releasesInterval
  }

  # generate other parameters
  genotypes=driveCube$genotypesID
  releaseVector=rep(0,length=length(genotypes))
  names(releaseVector)=genotypes
  releaseVector[driveCube$releaseType]=releasesParameters$releaseProportion

  # make releases vector
  basicRepeatedReleases(
    genotypes=genotypes,
    releaseStart=start,
    releaseEnd=end,
    releaseInterval=interval,
    releaseVector=releaseVector,
    sex=sex
  )
}

#' Make List of Batch Migration Parameters
#'
#' Sets up a list containing the probability of a batch migration, the fractional amount of males/females
#' that migrate, and the weighted probabilities for where to migrate.
#'
#' @param batchProbs Probability of a batch migration, either 1 number or vector of length equal to the number of patches
#' @param sexProbs Population fraction of males and females that migration. Either vector c(M,F) or matrix of 2 columns
#' @param numPatches Number of patches in the simulation
#'
#' @examples
#' # to setup for 3 patches
#' batchMigration = basicBatchMigration(batchProbs = 1e-5, sexProbs = c(0.1, 0.01), numPatches = 3)
#'
#' @export
basicBatchMigration <- function(batchProbs = 1e-5, sexProbs = c(0.01, 0.01),
                                numPatches = 1){

  # check length of probs
  if(!all(batchProbs<1)){
    stop("Probability of batch migration must be less than 1")
  }
  if(length(batchProbs) != numPatches){
    batchProbs = rep(x = batchProbs, numPatches)
  }

  # check length of sexes, make sure less than 1
  if(!all(sexProbs<1)){
    stop("Sex specific movement fraction must be less than 1")
  }
  if(is.null(dim(sexProbs))){
    sexProbsMat = matrix(data = sexProbs, nrow = numPatches, ncol = 2,
                      byrow = TRUE, dimnames = list(NULL, c("M","F")))
  }
  sexProbsCube = array(data = c(1-sexProbsMat, sexProbsMat), dim = c(numPatches, 2, 2))

  # setup movement matrix
  moveMat <- matrix(data = 1L, nrow = numPatches, ncol = numPatches)
  diag(moveMat) <- 0L
  moveMat <- moveMat/rowSums(x = moveMat)

  # return basic batch migration
  return(list("batchProbs" = batchProbs,
              "sexProbs" = sexProbsCube,
              "moveProbs" = moveMat)
         )
}
