###############################################################################
#       __  _____________       _       ______
#      /  |/  / ____/ __ \_____(_)   __/ ____/
#     / /|_/ / / __/ / / / ___/ / | / / __/
#    / /  / / /_/ / /_/ / /  / /| |/ / /___
#   /_/  /_/\____/_____/_/  /_/ |___/_____/
#
#   MGDrivE Releases
#   Marshall Lab
#   jared_bennett@berkeley.edu
#   December 2019
#
###############################################################################

#' Make List of Modified Mosquito Releases
#'
#' Sets up a release schedule for a single patch, returns a list to be used in
#' \code{\link{oneDay_releases_Patch}} or \code{\link{oneDay_eggReleases_Patch}}.
#' This function is no longer intended to be used alone, please use the standard
#' interface, \code{\link{generateReleaseVector}}.
#'
#' @param releaseStart Day releases start
#' @param releaseEnd Day releases end
#' @param releaseInterval Interval between releases
#' @param releaseMatrix Numeric matrix specifying the genotype and release amount
#'
#' @examples
#' \dontrun{
#' # Setup for 3 patches but only release in the first with a defined release
#' #  schedule, for the cube cubeHomingDrive:
#'
#' patchReleases = replicate(n = 3, expr = {
#'   list(maleReleases = NULL, femaleReleases = NULL, eggReleases = NULL, matedFemaleReleases = NULL)
#' },simplify = FALSE)
#'
#' patchReleases[[1]]$femaleReleases = MGDrivE::basicRepeatedReleases(releaseStart = 5,
#'                                                           releaseEnd = 30,
#'                                                           releaseInterval = 5,
#'                                                           releaseMatrix = matrix(c(5,100),1,2))
#'
#' patchReleases[[1]]$maleReleases = MGDrivE::basicRepeatedReleases(releaseStart = 50,
#'                                                         releaseEnd = 60,
#'                                                         releaseInterval = 1,
#'                                                         releaseMatrix = matrix(c(5,100),1,2))
#' }
#'
basicRepeatedReleases <- function(releaseStart, releaseEnd,
                                  releaseInterval, releaseMatrix){

  # check timing of releases
  if(releaseInterval > (releaseEnd - releaseStart)){
    stop("interval between releases cannot be greater than time between start and end of releases")
  }

  # create release times and setup list with releases
  releaseList <- lapply(X = seq(from=releaseStart,to = releaseEnd,by = floor(releaseInterval)),
                        FUN = function(x){
                          list("nRelease"=releaseMatrix,"tRelease"=x)
                        })

  # return
  return(releaseList)
}

#' Make List of Modified Mosquito Releases
#'
#' Sets up a release schedule for a single patch, calls
#' \code{\link{basicRepeatedReleases}} internally.
#'
#' @param driveCube Gene-drive cube
#' @param releasesParameters A list containing the releasesStart, releasesNumber
#' releasesInterval, and releaseProportion named values.
#' @param nameGenotypes Optional list to specify different genotypes for egg/male/female
#' releases. This is required for mated female releases. This parameter overrides
#' the default release type.
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
#'                                  releasesParameters = relParams)
#'
#' # generate mated female releases
#' fRelVec <- generateReleaseVector(driveCube = cube,
#'                                  releasesParameters = relParams,
#'                                  nameGenotypes = list(c("AA","AA", 10),
#'                                                       c("AA","aa", 10)))
#'
#' @export
generateReleaseVector <- function(driveCube, releasesParameters, nameGenotypes = NULL){

  ##########
  # edge case
  ##########
  if(releasesParameters$releasesNumber==0L){
    return(NULL)
  }

  ##########
  # generate timing parameters
  ##########
  start = releasesParameters$releasesStart
  end = releasesParameters$releasesStart + releasesParameters$releasesInterval * (releasesParameters$releasesNumber-1)
  interval = releasesParameters$releasesInterval
  if(start == end) interval = 0

  ##########
  # generate release matrix
  ##########
  if(is.null(nameGenotypes)){
    # default behaviour
    # This generates a release matrix for the pre-specified release type in the cube
    releaseMat <- matrix(data = c(which(driveCube$releaseType == driveCube$genotypesID),
                                        releasesParameters$releaseProportion),
                         byrow = TRUE, nrow = 1, ncol = 2)

  } else if(all(lengths(nameGenotypes) == 2) ) {
    # user-specified released for eggs, males, unmated females
    #  Takes the input list of genotypes and amount, builds the releases from that
    #  Allows the user to release several types at once, not necessarily the default release type

    # check for correctly specified genotypes
    if(!all(lapply(X = nameGenotypes, FUN = "[[", 1) %in% driveCube$genotypesID)){
      stop(paste0("Genotypes specified for egg/male/female release is not one of the
                  genotypes in the specified inheritance pattern."))
    }

    # setup and fill releases
    releaseMat <- matrix(data = 0, nrow = length(nameGenotypes), ncol = 2)
    for(i in 1:length(nameGenotypes)){
      releaseMat[i,1] <- which(nameGenotypes[[i]][1] == driveCube$genotypesID)
      releaseMat[i,2] <- as.numeric(nameGenotypes[[i]][2])
    } # end loop over what to release

  } else if(all(lengths(nameGenotypes) == 3) ) {
    # mated female releases
    # takes input list of female genotype, mate genotype, and number to release
    #  builds the releases from that.

    # check for correctly specified genotypes
    if(!all(lapply(X = nameGenotypes, FUN = "[[", 1) %in% driveCube$genotypesID)){
      stop(paste0("Genotypes specified for female in mated-female release is not
                  one of the genotypes in the specified inheritance pattern."))
    }
    if(!all(lapply(X = nameGenotypes, FUN = "[[", 2) %in% driveCube$genotypesID)){
      stop(paste0("Genotypes specified for male mate in mated-female release is not
                  one of the genotypes in the specified inheritance pattern."))
    }

    # setup and fill releases
    releaseMat <- matrix(data = 0, nrow = length(nameGenotypes), ncol = 3)
    for(i in 1:length(nameGenotypes)){
      # female genotype
      releaseMat[i,1] <- which(nameGenotypes[[i]][1] == driveCube$genotypesID)
      # male genotype
      releaseMat[i,2] <- which(nameGenotypes[[i]][2] == driveCube$genotypesID)
      # number to release
      releaseMat[i,3] <- as.numeric(nameGenotypes[[i]][3])
    } # end loop over releaseList

  } else {
    stop("nameGenotypes is incorrectly formatted.\n
         Check the examples or set as NULL for default behaviour.")
  }

  ##########
  # make/return releases list
  ##########
  basicRepeatedReleases(releaseStart=start,
                        releaseEnd=end,
                        releaseInterval=interval,
                        releaseMatrix=releaseMat)
}

#' Make List of Batch Migration Parameters
#'
#' Sets up a list containing the probability of a batch migration, the fractional
#' amount of males/females that migrate, and the weighted probabilities for where
#' to migrate. The default weights for migration are equal for all patches.
#' These can be changed after running the function. This is only used in
#' \code{\link{oneDay_migrationOut_stochastic_Patch}}.
#'
#' @param batchProbs Probability of a batch migration, either 1 number or a vector
#' of length equal to the number of patches
#' @param sexProbs Population fraction of males and females that migrate. Either
#' a vector c(M,F) or matrix of 2 columns
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
  if(length(batchProbs) == 1){
    batchProbs = rep(x = batchProbs, numPatches)
  } else if(length(batchProbs) != numPatches){
    stop("Probability of batch migration should be length 1 (all patches specified the same)
         or length numPatches (all patches specified individually).")
  }

  # check length of sexes, make sure less than 1
  if(!all(sexProbs<1)){
    stop("Sex specific movement fraction must be less than 1")
  }
  if(is.null(dim(sexProbs))){
    sexProbs = matrix(data = sexProbs, nrow = numPatches, ncol = 2,
                      byrow = TRUE, dimnames = list(NULL, c("M","F")))
  } else if(dim(sexProbs)[1] != numPatches){
    stop("Probability of male/female migration should have 1 row (all patches specified the same)
     or numPatches of rows (all patches specified individually).")
  }

  # setup movement matrix
  #  This is uniform everywhere. IDK a better way to let users specify.
  moveMat <- matrix(data = 1L, nrow = numPatches, ncol = numPatches)
  diag(moveMat) <- 0L
  moveMat <- moveMat/rowSums(x = moveMat)

  # return basic batch migration
  return(list("batchProbs" = batchProbs,
              "sexProbs" = sexProbs,
              "moveProbs" = moveMat)
         )
}
