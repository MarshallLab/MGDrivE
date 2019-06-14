## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE,
  hold = TRUE,
  fig.width = 7,
  fig.height = 10.5,
  eval = TRUE
)

# set seed for reproduibility
set.seed(seed = 42)

## ---- eval=FALSE---------------------------------------------------------
#  ####################
#  # Load libraries
#  ####################
#  library(MGDrivE)
#  
#  ####################
#  # Output Folder
#  ####################
#  outFolder <- "mgdrive"
#  dir.create(path = outFolder)
#  
#  ####################
#  # Simulation Parameters
#  ####################
#  # days to run the simulation
#  tMax <- 365
#  
#  # number of Monte Carlo iterations
#  nRep <- 5
#  
#  # each Monte Carlo iteration gets its own folder
#  folderNames <- file.path(outFolder,
#                          formatC(x = 1:nRep, width = 3, format = "d", flag = "0"))
#  
#  # entomological parameters
#  bioParameters <- list(betaK=20,tEgg=5,tLarva=6,tPupa=4,popGrowth=1.175,muAd=0.09)
#  
#  # a 5-node network with 5% per day migration rate
#  sitesNumber <- 5
#  adultPopEquilibrium <- 500
#  patchPops <- rep(adultPopEquilibrium,sitesNumber)
#  
#  # auxiliary function
#  triDiag <- function(upper, lower){
#  
#    # return matrix
#    retMat <- matrix(data = 0, nrow = length(upper) + 1, ncol = length(upper) + 1)
#  
#    # set index values for upper/lower triangles
#    indx <- 1:length(upper)
#  
#    # set forward/backward migration using matrix access
#    retMat[cbind(indx+1,indx)] <- lower
#    retMat[cbind(indx,indx+1)] <- upper
#  
#    # set stay probs
#    diag(x = retMat) <- 1-rowSums(x = retMat)
#  
#    return(retMat)
#  }
#  
#  # fill movement matrix
#  #  Remember, rows need to sum to 1.
#  moveMat <- triDiag(upper = rep.int(x = 0.05, times = sitesNumber-1),
#                     lower = rep.int(x = 0.05, times = sitesNumber-1))
#  
#  # batch migration is disabled by setting the probability to 0
#  batchMigration <- basicBatchMigration(batchProbs=0,
#                                       sexProbs=c(.5,.5),
#                                       numPatches=sitesNumber)
#  
#  ####################
#  # Basic Inheritance pattern
#  ####################
#  # Mendelian cube with standard (default) paraameters
#  cube <- cubeMendelian()
#  
#  
#  ####################
#  # Setup releases and batch migration
#  ####################
#  # set up the empty release vector
#  #  MGDrivE pulls things out by name
#  patchReleases <- replicate(n=sitesNumber,
#                            expr={list(maleReleases=NULL,femaleReleases=NULL)},
#                            simplify=FALSE)
#  
#  # choose release parameters
#  #  Releases start at time 25, occur every day, for 5 day.
#  #  There are 50 mosquitoes relaeses every time.
#  releasesParameters <- list(releasesStart=25,
#                            releasesNumber=10,
#                            releasesInterval=7,
#                            releaseProportion=100)
#  
#  # generate male release vector
#  maleReleasesVector <- generateReleaseVector(driveCube=cube,
#                                             releasesParameters=releasesParameters,
#                                             sex="M")
#  
#  # generate female release vector
#  femaleReleasesVector <- generateReleaseVector(driveCube=cube,
#                                               releasesParameters=releasesParameters,
#                                               sex="F")
#  
#  # put releases into the proper place in the release list
#  patchReleases[[1]]$maleReleases <- maleReleasesVector
#  patchReleases[[1]]$femaleReleases <- femaleReleasesVector
#  
#  
#  ####################
#  # Combine parameters and run!
#  ####################
#  # setup parameters for the network. This builds a list of parameters required for
#  #  every population in the network.
#  netPar <- parameterizeMGDrivE(runID=1, simTime=tMax, nPatch=sitesNumber,
#                               beta=bioParameters$betaK, muAd=bioParameters$muAd,
#                               popGrowth=bioParameters$popGrowth, tEgg=bioParameters$tEgg,
#                               tLarva=bioParameters$tLarva, tPupa=bioParameters$tPupa,
#                               AdPopEQ=patchPops)
#  
#  # set MGDrivE to run stochastic
#  setupMGDrivE(stochasticityON = TRUE)
#  
#  # build network prior to run
#  MGDrivESim <- Network$new(params=netPar,
#                           driveCube=cube,
#                           patchReleases=patchReleases,
#                           migrationMale=moveMat,
#                           migrationFemale=moveMat,
#                           migrationBatch=batchMigration,
#                           directory=folderNames)
#  # run simulation
#  MGDrivESim$multRun(verbose = FALSE)
#  
#  
#  ####################
#  # Post Analysis
#  ####################
#  # First, split output by patch
#  # Second, aggregate females by their mate choice
#  for(i in 1:nRep){
#   splitOutput(readDir=folderNames[i],remFile = TRUE, verbose = FALSE)
#   aggregateFemales(readDir=folderNames[i],genotypes=cube$genotypesID,
#                    remFile=TRUE,verbose = FALSE)
#  }
#  
#  # plot output of first run to see effect
#  plotMGDrivESingle(readDir=folderNames[1],totalPop = TRUE,lwd=3.5,alpha=1)
#  
#  # plot all 5 repetitions together
#  plotMGDrivEMult(readDir=outFolder,lwd=0.35,alpha=0.75)

## ------------------------------------------------------------------------
####################
# Load libraries
####################
library(MGDrivE)


## ------------------------------------------------------------------------
####################
# Output Folder
####################
outFolder <- "mgdrive"
dir.create(path = outFolder)


## ------------------------------------------------------------------------

####################
# Simulation Parameters
####################
# days to run the simulation
tMax <- 365

# number of Monte Carlo iterations
nRep <- 5

# each Monte Carlo iteration gets its own folder
folderNames <- file.path(outFolder,
                        formatC(x = 1:nRep, width = 3, format = "d", flag = "0"))

folderNames

## ------------------------------------------------------------------------
# entomological parameters
bioParameters <- list(betaK=20,tEgg=5,tLarva=6,tPupa=4,popGrowth=1.175,muAd=0.09)

## ------------------------------------------------------------------------
# a 5-node network with 5% per day migration rate
sitesNumber <- 5
adultPopEquilibrium <- 500
patchPops <- rep(adultPopEquilibrium,sitesNumber)

# landscape
# auxiliary function
triDiag <- function(upper, lower){
  
  # return matrix
  retMat <- matrix(data = 0, nrow = length(upper) + 1, ncol = length(upper) + 1)
  
  # set index values for upper/lower triangles
  indx <- 1:length(upper)
  
  # set forward/backward migration using matrix access
  retMat[cbind(indx+1,indx)] <- lower
  retMat[cbind(indx,indx+1)] <- upper
  
  # set stay probs
  diag(x = retMat) <- 1-rowSums(x = retMat)
  
  return(retMat)
}

# fill movement matrix
#  Remember, rows need to sum to 1.
moveMat <- triDiag(upper = rep.int(x = 0.05, times = sitesNumber-1),
                   lower = rep.int(x = 0.05, times = sitesNumber-1))

moveMat

## ------------------------------------------------------------------------
# batch migration is disabled by setting the probability to 0
batchMigration <- basicBatchMigration(batchProbs=0, numPatches=sitesNumber)

## ------------------------------------------------------------------------
####################
# Basic Inheritance pattern
####################
# Mendelian cube with standard (default) paraameters
cube <- cubeMendelian()

## ------------------------------------------------------------------------
####################
# Setup releases and batch migration
####################
# set up the empty release vector
#  MGDrivE pulls things out by name
patchReleases <- replicate(n=sitesNumber,
                          expr={list(maleReleases=NULL,femaleReleases=NULL)},
                          simplify=FALSE)

# choose release parameters
#  Releases start at time 25, occur every day, for 5 day.
#  There are 50 mosquitoes relaeses every time.
releasesParameters <- list(releasesStart=25,
                          releasesNumber=10,
                          releasesInterval=7,
                          releaseProportion=100)

## ------------------------------------------------------------------------
# generate male release vector
maleReleasesVector <- generateReleaseVector(driveCube=cube, 
                                           releasesParameters=releasesParameters,
                                           sex="M")
maleReleasesVector[[1]]

## ------------------------------------------------------------------------
# generate female release vector
femaleReleasesVector <- generateReleaseVector(driveCube=cube,
                                             releasesParameters=releasesParameters,
                                             sex="F")

# put releases into the proper place in the release list
patchReleases[[1]]$maleReleases <- maleReleasesVector
patchReleases[[1]]$femaleReleases <- femaleReleasesVector


## ------------------------------------------------------------------------
####################
# Combine parameters and run!
####################
# setup parameters for the network. This builds a list of parameters required for
#  every population in the network.
netPar <- parameterizeMGDrivE(runID=1, simTime=tMax, nPatch=sitesNumber,
                             beta=bioParameters$betaK, muAd=bioParameters$muAd,
                             popGrowth=bioParameters$popGrowth, tEgg=bioParameters$tEgg,
                             tLarva=bioParameters$tLarva, tPupa=bioParameters$tPupa,
                             AdPopEQ=patchPops)

## ------------------------------------------------------------------------
# set MGDrivE to run stochastic
setupMGDrivE(stochasticityON = TRUE)

## ------------------------------------------------------------------------
# build network prior to run
MGDrivESim <- Network$new(params=netPar,
                         driveCube=cube,
                         patchReleases=patchReleases,
                         migrationMale=moveMat,
                         migrationFemale=moveMat,
                         migrationBatch=batchMigration,
                         directory=folderNames)

# list folders to show that they have been created
list.files(path = outFolder)

## ------------------------------------------------------------------------
# run simulation
MGDrivESim$multRun(verbose = FALSE)

## ------------------------------------------------------------------------
# first and last repetitions
list.files(path = outFolder)[1]
list.files(path = list.files(path = outFolder, full.names = TRUE)[1], recursive = TRUE)
list.files(path = outFolder)[5]
list.files(path = list.files(path = outFolder, full.names = TRUE)[5], recursive = TRUE)

## ------------------------------------------------------------------------
# read in male and female files
mMat <- as.matrix(read.csv(file = list.files(path = outFolder, full.names = TRUE,
                                             recursive = TRUE)[2],
                           header = TRUE, sep = ","))
fMat <- as.matrix(read.csv(file = list.files(path = outFolder, full.names = TRUE,
                                             recursive = TRUE)[1],
                           header = TRUE, sep = ","))

# look at male file header
colnames(mMat)

## ------------------------------------------------------------------------
head(x = mMat, n = 2*sitesNumber)

## ------------------------------------------------------------------------
# look at female file header
colnames(fMat)

## ------------------------------------------------------------------------
head(x = fMat, n = 2*sitesNumber)

## ------------------------------------------------------------------------
# First, split output by patch
for(i in 1:nRep){
 splitOutput(readDir=folderNames[i],remFile = TRUE, verbose = FALSE)
}

## ------------------------------------------------------------------------
# first and last repetitions
list.files(path = outFolder)[1]
list.files(path = list.files(path = outFolder, full.names = TRUE)[1], recursive = TRUE)
list.files(path = outFolder)[5]
list.files(path = list.files(path = outFolder, full.names = TRUE)[5], recursive = TRUE)

# read in examples of new male/female files
twoFiles <- list.files(path = outFolder, full.names = TRUE,
                       recursive = TRUE)[c(1,6)]

# read in male and female files
mMat <- as.matrix(read.csv(file = twoFiles[2], header = TRUE, sep = ","))
fMat <- as.matrix(read.csv(file = twoFiles[1], header = TRUE, sep = ","))


## ------------------------------------------------------------------------
head(x = mMat, n = 5)

## ------------------------------------------------------------------------
head(x = fMat, n = 5)

## ------------------------------------------------------------------------
# Second, aggregate females by their mate choice
for(i in 1:nRep){
 aggregateFemales(readDir=folderNames[i],genotypes=cube$genotypesID,
                  remFile=TRUE,verbose = FALSE)
}

## ------------------------------------------------------------------------
# changed female files in first repetition
list.files(path = outFolder)[1]
list.files(path = list.files(path = outFolder, full.names = TRUE)[1], recursive = TRUE)

# read in examples of new female file
fMat2 <- as.matrix(read.csv(file = list.files(path = outFolder,
                                              recursive = TRUE,
                                              full.names = TRUE)[1],
                            header = TRUE, sep = ","))

## ------------------------------------------------------------------------
head(x = fMat2, n = 5)

## ------------------------------------------------------------------------
# show non-aggregated female file split by patch
#  This is for patch number 1
fMat[(releasesParameters$releasesStart-2):(releasesParameters$releasesStart+2), ]
cat("\n")

# show aggregated female file
#  This is for patch number 1
fMat2[(releasesParameters$releasesStart-2):(releasesParameters$releasesStart+2), ]

## ------------------------------------------------------------------------
# plot the first repetition
plotMGDrivESingle(readDir=folderNames[1],totalPop = TRUE,lwd=3.5,alpha=1)

## ------------------------------------------------------------------------
# plot all 5 repetitions together
plotMGDrivEMult(readDir=outFolder,lwd=0.35,alpha=0.75)

## ---- echo=FALSE---------------------------------------------------------
####################
# Cleanup before next run
####################
unlink(x = outFolder, recursive = TRUE)
rm(list=ls())

