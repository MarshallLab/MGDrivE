## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  echo = TRUE,
  hold = TRUE,
  fig.width = 7,
  fig.height = 4,
  eval = TRUE
)

# set seed for reproducibility
set.seed(seed = 42)

## -----------------------------------------------------------------------------
# setup movement matrix for 1 node
moveMat <- matrix(data = 1, nrow = 1, ncol = 1)
moveMat

## -----------------------------------------------------------------------------
# setup movement matrix for 2 nodes

####################
# 2 nodes, no migration
####################
moveMat <- matrix(data = c(1,0,0,1), nrow = 2, ncol = 2, byrow = TRUE)
moveMat

####################
# 2 nodes, with migration
####################
# 5% migration per day from population 1
# 10% migraton per day from population 2
#  Notice that the rows sum to 1
moveMat <- matrix(data = c(0.95, 0.05, 0.10, 0.90),
                  nrow = 2, ncol = 2, byrow = TRUE)
moveMat

## -----------------------------------------------------------------------------
# setup random movement matrix for 5 nodes

####################
# 5 nodes
####################
nNodes <- 5

# fill with random data
moveMat <- matrix(data = runif(n = nNodes*nNodes), nrow = nNodes, ncol = nNodes)

# normalize
moveMat <- moveMat/rowSums(x = moveMat)
moveMat

## -----------------------------------------------------------------------------
# setup line with 10 nodes

####################
# 10 nodes in a line
####################
nNodes <- 10

# define function for use
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
moveMat <- triDiag(upper = rep.int(x = 0.05, times = nNodes-1),
                   lower = rep.int(x = 0.05, times = nNodes-1))

moveMat

## -----------------------------------------------------------------------------
# realistic landscape

# matrix of coordinates as latitude/longitude pairs
lat_longs <- matrix(data = c(37.873507, -122.268181,
                             37.873578, -122.254430,
                             37.869806, -122.267639),
                    nrow = 3, ncol = 2, byrow = TRUE,
                    dimnames = list(NULL, c('Lat','Lon')))

# calculate distance matrix between points
# dmat <- MGDrivE::calcHaversine(latLongs = lat_longs)
# dmat <- MGDrivE::calcVinSph(latLongs = lat_longs)
distMat <- MGDrivE::calcVinEll(latLongs = lat_longs)

# calculate a zero-inflated movement kernal over the distances
# p0 is the probability, per day, that a mosquito doesn't move.
#  This is the value used in Code sample 1 from the paper, and in the examples in our
#  github repository.
# rate is the average migration rate per day, implying 1/rate is the average
#  migration distance. The average distance was estimated as ~55.5 meters per day,
#  which is the value used in Code sample 1 and in the examples on github.
p0 <- 0.991
rate <- 1/55.5

moveMat <- MGDrivE::calcHurdleExpKernel(distMat = distMat, rate = rate, p0 = p0)
moveMat

## -----------------------------------------------------------------------------
####################
# Load libraries
####################
library(MGDrivE)

####################
# Output Folder
####################
outFolder <- "mgdrive"

####################
# Simulation Parameters
####################
# days to run the simulation
tMax <- 365

# entomological parameters
bioParameters <- list(betaK=20, tEgg=5, tLarva=6, tPupa=4, popGrowth=1.175, muAd=0.09)

# a 1-node network where mosquitoes do not leave
moveMat <- matrix(data = 1, nrow = 1, ncol = 1)

# parameters of the population equilibrium
adultPopEquilibrium <- 500
sitesNumber <- nrow(moveMat)

####################
# Basic Inheritance pattern
####################
# Mendelian cube with standard (default) parameters
cube <- cubeMendelian()

####################
# Setup releases and batch migration
####################
# set up the empty release vector
#  MGDrivE pulls things out by name
patchReleases <- replicate(n=sitesNumber,
                           expr={list(maleReleases=NULL,femaleReleases=NULL,
                                      eggReleases=NULL,matedFemaleReleases=NULL)},
                           simplify=FALSE)

# choose release parameters
#  Releases start at time 25, occur every day, for 1 day.
#  There are 10 mosquitoes released every time.
releasesParameters <- list(releasesStart=25,
                           releasesNumber=1,
                           releasesInterval=0,
                           releaseProportion=10)

# generate release vector
releasesVector <- generateReleaseVector(driveCube=cube,
                                        releasesParameters=releasesParameters)

# put releases into the proper place in the release list
patchReleases[[1]]$maleReleases <- releasesVector
patchReleases[[1]]$femaleReleases <- releasesVector


# batch migration is disabled by setting the probability to 0
# This is required because of the stochastic simulations, but doesn't make sense
#  in a deterministic simulation.
batchMigration <- basicBatchMigration(batchProbs=0,
                                      sexProbs=c(.5,.5),
                                      numPatches=sitesNumber)

####################
# Combine parameters and run!
####################
# set MGDrivE to run deterministic
setupMGDrivE(stochasticityON = FALSE, verbose = FALSE)

# setup parameters for the network. This builds a list of parameters required for
#  every population in the network. In ths case, we havee a network of 1 population.
netPar <- parameterizeMGDrivE(runID=1, simTime=tMax, sampTime = 1, nPatch=sitesNumber,
                              beta=bioParameters$betaK, muAd=bioParameters$muAd,
                              popGrowth=bioParameters$popGrowth, tEgg=bioParameters$tEgg,
                              tLarva=bioParameters$tLarva, tPupa=bioParameters$tPupa,
                              AdPopEQ=adultPopEquilibrium, inheritanceCube = cube)

# build network prior to run
MGDrivESim <- Network$new(params=netPar,
                          driveCube=cube,
                          patchReleases=patchReleases,
                          migrationMale=moveMat,
                          migrationFemale=moveMat,
                          migrationBatch=batchMigration,
                          directory=outFolder,
                          verbose=FALSE)

# run simulation
MGDrivESim$oneRun(verbose = FALSE)

####################
# Post Analysis
####################
# split output by patch
#  Required for plotting later
splitOutput(readDir = outFolder, remFile = TRUE, verbose = FALSE)

# aggregate females by their mate choice
#  This reduces the female file to have the same columns as the male file
aggregateFemales(readDir = outFolder, genotypes = cube$genotypesID,
                 remFile = TRUE, verbose = FALSE)

# plot output to see effect
plotMGDrivESingle(readDir = outFolder, totalPop = TRUE, lwd = 3.5, alpha = 1)

## ---- echo=FALSE--------------------------------------------------------------
####################
# Theory Comparison
####################

# read in simulation files
totPop <- read.csv(file = file.path(outFolder, "M_Run001_Patch001.csv"),
                   header = TRUE, sep = ",")[ ,-1] +
          read.csv(file = file.path(outFolder, "F_Aggregate_Run001_Patch001.csv"),
                   header = TRUE, sep = ",")[ ,-1]

## ---- echo=FALSE--------------------------------------------------------------
####################
# Cleanup before next run
####################
unlink(x = outFolder, recursive = TRUE)
rm(list=ls())

## -----------------------------------------------------------------------------
####################
# Load libraries
####################
library(MGDrivE)

####################
# Output Folder
####################
outFolder <- "mgdrive"

####################
# Simulation Parameters
####################
# days to run the simulation
tMax <- 365*2

# entomological parameters
bioParameters <- list(betaK=20, tEgg=5, tLarva=6, tPupa=4, popGrowth=1.175, muAd=0.09)

# a 1-node network where mosquitoes do not leave
moveMat <- as.matrix(1)

# parameters of the population equilibrium
adultPopEquilibrium <- 500
sitesNumber <- nrow(moveMat)
patchPops <- rep(adultPopEquilibrium,sitesNumber)

####################
# Basic Inheritance pattern
####################
# Mendelian cube with standard (default) parameters

# This time, lets put fitness cost on the homozygotes, giving the heterozygote
#  an advantage

# These genotypes correspond to ones in the cube. Look at a base cube first,
#  then set this.
# homozygotes are 60% as fit as heterozygote over their entire lifetime
#  Since omega is the adult daily death rate, we use the built-in function to
#  calculate our desired lifetime cost as applied daily
dayOmega <- calcOmega(mu = bioParameters$muAd, lifespanReduction = 0.60)
omegaNew <- c("AA"=dayOmega, "aa"=dayOmega)

# setup cube
cube <- cubeMendelian(omega = omegaNew)

####################
# Setup releases and batch migration
####################
# set up the empty release vector
#  MGDrivE pulls things out by name
patchReleases <- replicate(n=sitesNumber,
                           expr={list(maleReleases=NULL,femaleReleases=NULL,
                                      eggReleases=NULL,matedFemaleReleases=NULL)},
                           simplify=FALSE)

# choose release parameters
#  Releases start at time 100, occur every day, for 5 days.
#  There are 50 mosquitoes released every time.
releasesParameters <- list(releasesStart=100,
                           releasesNumber=5,
                           releasesInterval=0,
                           releaseProportion=50)

# generate male release vector
maleReleasesVector <- generateReleaseVector(driveCube=cube,
                                            releasesParameters=releasesParameters)

# generate female release vector
femaleReleasesVector <- generateReleaseVector(driveCube=cube,
                                              releasesParameters=releasesParameters)

# put releases into the proper place in the release list
patchReleases[[1]]$maleReleases <- maleReleasesVector
patchReleases[[1]]$femaleReleases <- femaleReleasesVector


# batch migration is disabled by setting the probability to 0
# This is required because of the stochastic simulations, but doesn't make sense
#  in a deterministic simulation.
batchMigration <- basicBatchMigration(batchProbs=0,
                                      sexProbs=c(.5,.5),
                                      numPatches=sitesNumber)

####################
# Combine parameters and run!
####################
# set MGDrivE to run deterministic
setupMGDrivE(stochasticityON = FALSE, verbose = FALSE)

# setup parameters for the network. This builds a list of parameters required for
#  every population in the network. In ths case, we have a network of 1 population.
netPar <- parameterizeMGDrivE(runID=1, simTime=tMax, nPatch=sitesNumber,
                              beta=bioParameters$betaK, muAd=bioParameters$muAd,
                              popGrowth=bioParameters$popGrowth, tEgg=bioParameters$tEgg,
                              tLarva=bioParameters$tLarva, tPupa=bioParameters$tPupa,
                              AdPopEQ=patchPops, inheritanceCube = cube)

# build network prior to run
MGDrivESim <- Network$new(params=netPar,
                          driveCube=cube,
                          patchReleases=patchReleases,
                          migrationMale=moveMat,
                          migrationFemale=moveMat,
                          migrationBatch=batchMigration,
                          directory=outFolder,
                          verbose=FALSE)
# run simulation
MGDrivESim$oneRun(verbose = FALSE)

####################
# Post Analysis
####################
# split output by patch
#  Required for plotting later
splitOutput(readDir = outFolder, remFile = TRUE, verbose = FALSE)

# aggregate females by their mate choice
#  This reduces the female file to have the same columns as the male file
aggregateFemales(readDir = outFolder, genotypes = cube$genotypesID, remFile = TRUE,
                 verbose = FALSE)

# plot output to see effect
plotMGDrivESingle(readDir = outFolder, totalPop = TRUE, lwd = 3.5, alpha = 1)

## ---- echo=FALSE--------------------------------------------------------------
####################
# Theory Comparison
####################

# read in simulation files
totPop <- read.csv(file = file.path(outFolder, "M_Run001_Patch001.csv"),
                   header = TRUE, sep = ",")[ ,-1] +
  read.csv(file = file.path(outFolder, "F_Aggregate_Run001_Patch001.csv"),
           header = TRUE, sep = ",")[ ,-1]

## ---- echo=FALSE--------------------------------------------------------------
####################
# Cleanup before next run
####################
unlink(x = outFolder, recursive = TRUE)
rm(list=ls())

## -----------------------------------------------------------------------------
####################
# Load libraries
####################
library(MGDrivE)

####################
# Output Folder
####################
# directory
# This is slightly obtuse for vignette building reasons
#  Really, all you need is a base directory, then the repetitions in folders inside that.
#  Here, our base directory is "mgdrive", and the repetition folders are "001","002", etc.
#  So, the final structure is "mgdrive/001","mgdrive/002", etc.
outFolder <- "mgdrive"
dir.create(path = outFolder)

####################
# Simulation Parameters
####################
# days to run the simulation
tMax <- 365*2

# number of Monte Carlo iterations
nRep <- 50

# each Monte Carlo iteration gets its own folder
folderNames <- file.path(outFolder,
                         formatC(x = 1:nRep, width = 3, format = "d", flag = "0"))

# entomological parameters
bioParameters <- list(betaK=20, tEgg=5, tLarva=6, tPupa=4, popGrowth=1.175, muAd=0.09)

# a 1-node network where mosquitoes do not leave
moveMat <- as.matrix(1)

# parameters of the population equilibrium
adultPopEquilibrium <- 500
sitesNumber <- nrow(moveMat)

####################
# Basic Inheritance pattern
####################
# Mendelian cube with standard (default) parameters

# This time, lets put fitness cost on the homozygotes, giving the heterozygote
#  an advantage

# These genotypes correspond to ones in the cube. Look at a base cube first,
#  then set this.
# homozygotes are 60% as fit as heterozygote over their entire lifetime
#  Since omega is the adult daily death rate, we use the built-in function to
#  calculate our desired lifetime cost as applied daily
dayOmega <- calcOmega(mu = bioParameters$muAd, lifespanReduction = 0.60)
omegaNew <- c("AA"=dayOmega, "aa"=dayOmega)

# setup cube
cube <- cubeMendelian(omega = omegaNew)

####################
# Setup releases and batch migration
####################
# set up the empty release vector
#  MGDrivE pulls things out by name
patchReleases <- replicate(n=sitesNumber,
                           expr={list(maleReleases=NULL,femaleReleases=NULL,
                                      eggReleases=NULL,matedFemaleReleases=NULL)},
                           simplify=FALSE)

# choose release parameters
#  Releases start at time 100, occur every day, for 5 days.
#  There are 50 mosquitoes released every time.
releasesParameters <- list(releasesStart=100,
                           releasesNumber=5,
                           releasesInterval=0,
                           releaseProportion=50)

# generate male release vector
maleReleasesVector <- generateReleaseVector(driveCube=cube,
                                            releasesParameters=releasesParameters)

# generate female release vector
femaleReleasesVector <- generateReleaseVector(driveCube=cube,
                                              releasesParameters=releasesParameters)

# put releases into the proper place in the release list
patchReleases[[1]]$maleReleases <- maleReleasesVector
patchReleases[[1]]$femaleReleases <- femaleReleasesVector


# batch migration is disabled by setting the probability to 0
# This is required because of the stochastic simulations, but doesn't make sense
#  in a deterministic simulation.
batchMigration <- basicBatchMigration(batchProbs=0,
                                      sexProbs=c(.5,.5),
                                      numPatches=sitesNumber)

####################
# Combine parameters and run!
####################
# set MGDrivE to run stochastic
setupMGDrivE(stochasticityON = TRUE, verbose = FALSE)

# setup parameters for the network. This builds a list of parameters required for
#  every population in the network. In ths case, we have a network of 1 population.
netPar <- parameterizeMGDrivE(runID=1, simTime=tMax, sampTime = 5, nPatch=sitesNumber,
                              beta=bioParameters$betaK, muAd=bioParameters$muAd,
                              popGrowth=bioParameters$popGrowth, tEgg=bioParameters$tEgg,
                              tLarva=bioParameters$tLarva, tPupa=bioParameters$tPupa,
                              AdPopEQ=adultPopEquilibrium, inheritanceCube = cube)

# build network prior to run
MGDrivESim <- Network$new(params=netPar,
                          driveCube=cube,
                          patchReleases=patchReleases,
                          migrationMale=moveMat,
                          migrationFemale=moveMat,
                          migrationBatch=batchMigration,
                          directory=folderNames,
                          verbose = FALSE)
# run simulation
MGDrivESim$multRun(verbose = FALSE)

####################
# Post Analysis
####################
# First, split output by patch
# Second, aggregate females by their mate choice
for(i in 1:nRep){
  splitOutput(readDir = folderNames[i], remFile = TRUE, verbose = FALSE)
  aggregateFemales(readDir = folderNames[i], genotypes = cube$genotypesID,
                   remFile = TRUE, verbose = FALSE)
}

# plot output of first run to see effect
plotMGDrivESingle(readDir = folderNames[1], totalPop = TRUE, lwd = 3.5, alpha = 1)

# plot all 50 repetitions together
plotMGDrivEMult(readDir = outFolder, lwd = 0.35, alpha = 0.75)

## ---- echo=FALSE--------------------------------------------------------------
####################
# Cleanup before next run
####################
unlink(x = outFolder, recursive = TRUE)
rm(list=ls())

## -----------------------------------------------------------------------------
####################
# Load libraries
####################
library(MGDrivE)

####################
# Output Folder
####################
outFolder <- "mgdrive"

####################
# Simulation Parameters
####################
# days to run the simulation
tMax <- 365

# entomological parameters
bioParameters <- list(betaK=20, tEgg=5, tLarva=6, tPupa=4, popGrowth=1.175, muAd=0.09)

# a 2-node network where mosquitoes do not leave
moveMat <- matrix(data = c(1,0,0,1), nrow = 2, ncol = 2)
moveMat

# parameters of the population equilibrium
adultPopEquilibrium <- 500
sitesNumber <- nrow(moveMat)

####################
# Basic Inheritance pattern
####################
# Mendelian cube with standard (default) parameters
cube <- cubeMendelian()

####################
# Setup releases and batch migration
####################
# set up the empty release vector
#  MGDrivE pulls things out by name
patchReleases <- replicate(n=sitesNumber,
                           expr={list(maleReleases=NULL,femaleReleases=NULL,
                                      eggReleases=NULL,matedFemaleReleases=NULL)},
                           simplify=FALSE)

# choose release parameters
#  Releases start at time 25, occur every day, for 5 days.
#  There are 50 mosquitoes released every time.
releasesParameters <- list(releasesStart=25,
                           releasesNumber=5,
                           releasesInterval=0,
                           releaseProportion=50)

# generate release vector
releasesVector <- generateReleaseVector(driveCube=cube,
                                        releasesParameters=releasesParameters)

# put releases into the proper place in the release list
patchReleases[[1]]$maleReleases <- releasesVector
patchReleases[[1]]$femaleReleases <- releasesVector

# batch migration is disabled by setting the probability to 0
# This is required because of the stochastic simulations, but doesn't make sense
#  in a deterministic simulation.
batchMigration <- basicBatchMigration(batchProbs=0,
                                      sexProbs=c(.5,.5),
                                      numPatches=sitesNumber)

####################
# Combine parameters and run!
####################
# set MGDrivE to run deterministic
setupMGDrivE(stochasticityON = FALSE, verbose = FALSE)

# setup parameters for the network. This builds a list of parameters required for
#  every population in the network. In ths case, we have a network of 2 populations.
netPar <- parameterizeMGDrivE(runID=1, simTime=tMax, nPatch=sitesNumber,
                              beta=bioParameters$betaK, muAd=bioParameters$muAd,
                              popGrowth=bioParameters$popGrowth, tEgg=bioParameters$tEgg,
                              tLarva=bioParameters$tLarva, tPupa=bioParameters$tPupa,
                              AdPopEQ=adultPopEquilibrium, inheritanceCube = cube)

# build network prior to run
MGDrivESim <- Network$new(params=netPar,
                          driveCube=cube,
                          patchReleases=patchReleases,
                          migrationMale=moveMat,
                          migrationFemale=moveMat,
                          migrationBatch=batchMigration,
                          directory=outFolder,
                          verbose=FALSE)
# run simulation
MGDrivESim$oneRun(verbose = FALSE)

####################
# Post Analysis
####################
# split output by patch
#  Required for plotting later
splitOutput(readDir = outFolder, verbose = FALSE, remFile = TRUE)

# aggregate females by their mate choice
#  This reduces the female file to have the same columns as the male file
aggregateFemales(readDir = outFolder, genotypes = cube$genotypesID,
                 remFile = TRUE, verbose = FALSE)

# plot output to see effect
plotMGDrivESingle(readDir = outFolder, lwd = 3.5, alpha = 1)

## ---- echo=FALSE--------------------------------------------------------------
####################
# Cleanup before next run
####################
unlink(x = outFolder, recursive = TRUE)
rm(list=ls())

## -----------------------------------------------------------------------------
####################
# Load libraries
####################
library(MGDrivE)

####################
# Output Folder
####################
outFolder <- "mgdrive"

####################
# Simulation Parameters
####################
# days to run the simulation
tMax <- 365

# entomological parameters
bioParameters <- list(betaK=20, tEgg=5, tLarva=6, tPupa=4, popGrowth=1.175, muAd=0.09)

# a 2-node network with 1% per day migration rate
moveMat <- matrix(data = c(0.99,0.01,0.01,0.99), nrow = 2, ncol = 2)
moveMat

# parameters of the population equilibrium
adultPopEquilibrium <- 500
sitesNumber <- nrow(moveMat)
patchPops <- rep(adultPopEquilibrium,sitesNumber)

####################
# Basic Inheritance pattern
####################
# Mendelian cube with standard (default) parameters
cube <- cubeMendelian()

####################
# Setup releases and batch migration
####################
# set up the empty release vector
#  MGDrivE pulls things out by name
patchReleases <- replicate(n=sitesNumber,
                           expr={list(maleReleases=NULL,femaleReleases=NULL,
                                      eggReleases=NULL,matedFemaleReleases=NULL)},
                           simplify=FALSE)

# choose release parameters
#  Releases start at time 25, occur every day, for 5 days.
#  There are 50 mosquitoes released every time.
releasesParameters <- list(releasesStart=25,
                           releasesNumber=5,
                           releasesInterval=0,
                           releaseProportion=50)

# generate male release vector
maleReleasesVector <- generateReleaseVector(driveCube=cube,
                                            releasesParameters=releasesParameters)

# generate female release vector
femaleReleasesVector <- generateReleaseVector(driveCube=cube,
                                              releasesParameters=releasesParameters)

# put releases into the proper place in the release list
patchReleases[[1]]$maleReleases <- maleReleasesVector
patchReleases[[1]]$femaleReleases <- femaleReleasesVector

# batch migration is disabled by setting the probability to 0
# This is required because of the stochastic simulations, but doesn't make sense
#  in a deterministic simulation.
batchMigration <- basicBatchMigration(batchProbs=0,
                                     sexProbs=c(.5,.5),
                                     numPatches=sitesNumber)

####################
# Combine parameters and run!
####################
# set MGDrivE to run deterministic
setupMGDrivE(stochasticityON = FALSE, verbose = FALSE)

# setup parameters for the network. This builds a list of parameters required for
#  every population in the network. In ths case, we have a network of 2 populations.
netPar <- parameterizeMGDrivE(runID=1, simTime=tMax, nPatch=sitesNumber,
                             beta=bioParameters$betaK, muAd=bioParameters$muAd,
                             popGrowth=bioParameters$popGrowth, tEgg=bioParameters$tEgg,
                             tLarva=bioParameters$tLarva, tPupa=bioParameters$tPupa,
                             AdPopEQ=patchPops, inheritanceCube = cube)

# build network prior to run
MGDrivESim <- Network$new(params=netPar,
                         driveCube=cube,
                         patchReleases=patchReleases,
                         migrationMale=moveMat,
                         migrationFemale=moveMat,
                         migrationBatch=batchMigration,
                         directory=outFolder,
                         verbose = FALSE)
# run simulation
MGDrivESim$oneRun(verbose = FALSE)

####################
# Post Analysis
####################
# split output by patch
#  Required for plotting later
splitOutput(readDir = outFolder, verbose = FALSE, remFile = TRUE)

# aggregate females by their mate choice
#  This reduces the female file to have the same columns as the male file
aggregateFemales(readDir = outFolder, genotypes = cube$genotypesID,
                remFile = TRUE, verbose = FALSE)

# plot output to see effect
plotMGDrivESingle(readDir = outFolder, totalPop = TRUE, lwd = 3.5, alpha = 1)

## ---- echo=FALSE--------------------------------------------------------------
####################
# Cleanup before next run
####################
unlink(x = outFolder, recursive = TRUE)
rm(list=ls())

## -----------------------------------------------------------------------------
####################
# Load libraries
####################
library(MGDrivE)

####################
# Output Folder
####################
# This is slightly obtuse for vignette building reasons
#  Really, all you need is a base directory, then the repetitions in folders inside that.
#  Here, our base directory is "mgdrive", and the repetition folders are "001","002", etc.
#  So, the final structure is "mgdrive/001","mgdrive/002", etc.
outFolder <- "mgdrive"
dir.create(path = outFolder)

####################
# Simulation Parameters
####################
# days to run the simulation
tMax <- 365

# number of Monte Carlo iterations
nRep <- 25

# each Monte Carlo iteration gets its own folder
folderNames <- file.path(outFolder,
                        formatC(x = 1:nRep, width = 3, format = "d", flag = "0"))

# entomological parameters
bioParameters <- list(betaK=20, tEgg=5, tLarva=6, tPupa=4, popGrowth=1.175, muAd=0.09)

# a 2-node network with 1% per day migration rate
moveMat <- matrix(data = c(0.99,0.01,0.01,0.99), nrow = 2, ncol = 2)
moveMat

# parameters of the population equilibrium
adultPopEquilibrium <- 500
sitesNumber <- nrow(moveMat)

####################
# Basic Inheritance pattern
####################
# Mendelian cube with standard (default) parameters
cube <- cubeMendelian()

####################
# Setup releases and batch migration
####################
# set up the empty release vector
#  MGDrivE pulls things out by name
patchReleases <- replicate(n=sitesNumber,
                           expr={list(maleReleases=NULL,femaleReleases=NULL,
                                      eggReleases=NULL,matedFemaleReleases=NULL)},
                           simplify=FALSE)

# choose release parameters
#  Releases start at time 25, occur every day, for 5 days.
#  There are 50 mosquitoes released every time.
releasesParameters <- list(releasesStart=25,
                          releasesNumber=5,
                          releasesInterval=0,
                          releaseProportion=50)

# generate release vector
releasesVector <- generateReleaseVector(driveCube=cube,
                                        releasesParameters=releasesParameters)

# put releases into the proper place in the release list
patchReleases[[1]]$maleReleases <- releasesVector
patchReleases[[1]]$femaleReleases <- releasesVector

# batch migration is disabled by setting the probability to 0
# This is required because of the stochastic simulations, but doesn't make sense
#  in a deterministic simulation.
batchMigration <- basicBatchMigration(batchProbs=0,
                                     sexProbs=c(.5,.5),
                                     numPatches=sitesNumber)

####################
# Combine parameters and run!
####################
# set MGDrivE to run stochastic
setupMGDrivE(stochasticityON = TRUE, verbose = FALSE)

# setup parameters for the network. This builds a list of parameters required for
#  every population in the network. In ths case, we have a network of 2 populations.
netPar <- parameterizeMGDrivE(runID=1, simTime=tMax, nPatch=sitesNumber,
                             beta=bioParameters$betaK, muAd=bioParameters$muAd,
                             popGrowth=bioParameters$popGrowth, tEgg=bioParameters$tEgg,
                             tLarva=bioParameters$tLarva, tPupa=bioParameters$tPupa,
                             AdPopEQ=adultPopEquilibrium, inheritanceCube = cube)

# build network prior to run
MGDrivESim <- Network$new(params=netPar,
                         driveCube=cube,
                         patchReleases=patchReleases,
                         migrationMale=moveMat,
                         migrationFemale=moveMat,
                         migrationBatch=batchMigration,
                         directory=folderNames,
                         verbose = FALSE)
# run simulation
MGDrivESim$multRun(verbose = FALSE)

####################
# Post Analysis
####################
# First, split output by patch
# Second, aggregate females by their mate choice
for(i in 1:nRep){
 splitOutput(readDir = folderNames[i], remFile = TRUE, verbose = FALSE)
 aggregateFemales(readDir = folderNames[i], genotypes = cube$genotypesID,
                  remFile = TRUE, verbose = FALSE)
}

# plot output of first run to see effect
plotMGDrivESingle(readDir = folderNames[1], totalPop = TRUE, lwd = 3.5, alpha = 1)

# plot all 25 repetitions together
plotMGDrivEMult(readDir = outFolder, lwd = 0.35, alpha = 0.75)

## ---- echo=FALSE--------------------------------------------------------------
####################
# Cleanup before next run
####################
unlink(x = outFolder, recursive = TRUE)
rm(list=ls())

## -----------------------------------------------------------------------------
####################
# Load libraries
####################
library(MGDrivE)

####################
# Output Folder
####################
# This is slightly obtuse for vignette building reasons
#  Really, all you need is a base directory, then the repetitions in folders inside that.
#  Here, our base directory is "mgdrive", and the repetition folders are "001","002", etc.
#  So, the final structure is "mgdrive/001","mgdrive/002", etc.
outFolder <- "mgdrive"
dir.create(path = outFolder)

####################
# Simulation Parameters
####################
# days to run the simulation
tMax <- 365*2

# number of Monte Carlo iterations
nRep <- 25

# each Monte Carlo iteration gets its own folder
folderNames <- file.path(outFolder,
                         formatC(x = 1:nRep, width = 3, format = "d", flag = "0"))

# entomological parameters
bioParameters <- list(betaK=20, tEgg=5, tLarva=6, tPupa=4, popGrowth=1.175, muAd=0.09)

# a 2-node network with 1% per day migration rate
moveMat <- matrix(data = c(0.99,0.01,0.01,0.99), nrow = 2, ncol = 2)

# parameters of the population equilibrium
adultPopEquilibrium <- 500
sitesNumber <- nrow(moveMat)

####################
# Basic Inheritance pattern
####################
# Mendelian cube

# This time, lets put fitness cost on the homozygotes, giving the heterozygote
#  an advantage

# These genotypes correspond to ones in the cube. Look at a base cube first,
#  then set this.
# Homozygotes are 60% as fit as heterozygote over their entire lifetime
#  Since omega is a daily death rate, we use the built-in function to calculate
#  our desired lifetime cost as applied daily
dayOmega <- calcOmega(mu = bioParameters$muAd, lifespanReduction = 0.60)
omegaNew <- c("AA"=dayOmega, "aa"=dayOmega)

# setup cube
cube <- cubeMendelian(omega = omegaNew)

####################
# Setup releases and batch migration
####################
# set up the empty release vector
#  MGDrivE pulls things out by name
patchReleases <- replicate(n=sitesNumber,
                           expr={list(maleReleases=NULL,femaleReleases=NULL,
                                      eggReleases=NULL,matedFemaleReleases=NULL)},
                           simplify=FALSE)

# choose release parameters
#  Releases start at time 100, occur every day, for 5 days.
#  There are 50 mosquitoes released every time.
releasesParameters <- list(releasesStart=100,
                           releasesNumber=5,
                           releasesInterval=0,
                           releaseProportion=50)

# generate male release vector
maleReleasesVector <- generateReleaseVector(driveCube=cube,
                                            releasesParameters=releasesParameters)

# generate female release vector
femaleReleasesVector <- generateReleaseVector(driveCube=cube,
                                              releasesParameters=releasesParameters)

# put releases into the proper place in the release list
patchReleases[[1]]$maleReleases <- maleReleasesVector
patchReleases[[1]]$femaleReleases <- femaleReleasesVector

# batch migration is disabled by setting the probability to 0
# This is required because of the stochastic simulations, but doesn't make sense
#  in a deterministic simulation.
batchMigration <- basicBatchMigration(batchProbs=0,
                                     sexProbs=c(.5,.5),
                                     numPatches=sitesNumber)

####################
# Combine parameters and run!
####################
# set MGDrivE to run stochastic
setupMGDrivE(stochasticityON = TRUE, verbose = FALSE)

# setup parameters for the network. This builds a list of parameters required for
#  every population in the network. In ths case, we have a network of 2 populations.
netPar <- parameterizeMGDrivE(runID=1, simTime=tMax, nPatch=sitesNumber,
                             beta=bioParameters$betaK, muAd=bioParameters$muAd,
                             popGrowth=bioParameters$popGrowth, tEgg=bioParameters$tEgg,
                             tLarva=bioParameters$tLarva, tPupa=bioParameters$tPupa,
                             AdPopEQ=adultPopEquilibrium, inheritanceCube = cube)

# build network prior to run
MGDrivESim <- Network$new(params=netPar,
                         driveCube=cube,
                         patchReleases=patchReleases,
                         migrationMale=moveMat,
                         migrationFemale=moveMat,
                         migrationBatch=batchMigration,
                         directory=folderNames,
                         verbose = FALSE)
# run simulation
MGDrivESim$multRun(verbose = FALSE)

####################
# Post Analysis
####################
# First, split output by patch
# Second, aggregate females by their mate choice
for(i in 1:nRep){
 splitOutput(readDir = folderNames[i], remFile = TRUE, verbose = FALSE)
 aggregateFemales(readDir = folderNames[i], genotypes = cube$genotypesID,
                  remFile = TRUE, verbose = FALSE)
}

# plot output of first run to see effect
#  per the structure above, we are reading "mgdrive/001" for the single plot
plotMGDrivESingle(readDir = folderNames[1], totalPop = TRUE, lwd = 3.5, alpha = 1)

# plot all 50 repetitions together
#  Here, we feed the function "mgdrive/", and it finds all repetition folders
#   inside that.
plotMGDrivEMult(readDir = outFolder, lwd = 0.35, alpha = 0.75)

## ---- echo=FALSE--------------------------------------------------------------
####################
# Cleanup before next run
####################
unlink(x = outFolder, recursive = TRUE)
rm(list=ls())

## -----------------------------------------------------------------------------
####################
# Load libraries
####################
library(MGDrivE)

####################
# Output Folder
####################
outFolder <- "mgdrive"

####################
# Simulation Parameters
####################
# days to run the simulation, 2 years
tMax <- 365*2

# entomological parameters
bioParameters <- list(betaK=20, tEgg=5, tLarva=6, tPupa=4, popGrowth=1.175, muAd=0.09)

# a 1-node network where mosquitoes do not leave
moveMat <- as.matrix(1)

# parameters of the population equilibrium
adultPopEquilibrium <- 500
sitesNumber <- nrow(moveMat)

####################
# Basic Inheritance pattern
####################
# Reciprocal translocation cube with standard (default) parameters
cube <- cubeReciprocalTranslocations()

####################
# Setup releases and batch migration
####################
# set up the empty release vector
#  MGDrivE pulls things out by name
patchReleases <- replicate(n=sitesNumber,
                           expr={list(maleReleases=NULL,femaleReleases=NULL,
                                      eggReleases=NULL,matedFemaleReleases=NULL)},
                           simplify=FALSE)

# choose release parameters
#  Releases start at time 25, occur once a week, for 5 weeks.
#  There are 100 mosquitoes released every time.
releasesParameters <- list(releasesStart=25,
                          releasesNumber=5,
                          releasesInterval=7,
                          releaseProportion=100)

# generate male release vector
maleReleasesVector <- generateReleaseVector(driveCube=cube,
                                           releasesParameters=releasesParameters)

# generate female release vector
femaleReleasesVector <- generateReleaseVector(driveCube=cube,
                                             releasesParameters=releasesParameters)

# put releases into the proper place in the release list
patchReleases[[1]]$maleReleases <- maleReleasesVector
patchReleases[[1]]$femaleReleases <- femaleReleasesVector

# batch migration is disabled by setting the probability to 0
# This is required because of the stochastic simulations, but doesn't make sense
#  in a deterministic simulation.
batchMigration <- basicBatchMigration(batchProbs=0,
                                     sexProbs=c(.5,.5),
                                     numPatches=sitesNumber)

####################
# Combine parameters and run!
####################
# set MGDrivE to run deterministic
setupMGDrivE(stochasticityON = FALSE, verbose = FALSE)

# setup parameters for the network. This builds a list of parameters required for
#  every population in the network. In ths case, we have a network of 1 population.
netPar <- parameterizeMGDrivE(runID=1, simTime=tMax, nPatch=sitesNumber,
                             beta=bioParameters$betaK, muAd=bioParameters$muAd,
                             popGrowth=bioParameters$popGrowth, tEgg=bioParameters$tEgg,
                             tLarva=bioParameters$tLarva, tPupa=bioParameters$tPupa,
                             AdPopEQ=adultPopEquilibrium, inheritanceCube = cube)

# build network prior to run
MGDrivESim <- Network$new(params=netPar,
                         driveCube=cube,
                         patchReleases=patchReleases,
                         migrationMale=moveMat,
                         migrationFemale=moveMat,
                         migrationBatch=batchMigration,
                         directory=outFolder,
                         verbose = FALSE)
# run simulation
MGDrivESim$oneRun(verbose = FALSE)

####################
# Post Analysis
####################
# split output by patch
#  Required for plotting later
splitOutput(readDir = outFolder, remFile = TRUE, verbose = FALSE)

# aggregate females by their mate choice
#  This reduces the female file to have the same columns as the male file
aggregateFemales(readDir = outFolder, genotypes = cube$genotypesID,
                remFile = TRUE, verbose = FALSE)

# plot output to see effect
plotMGDrivESingle(readDir = outFolder, totalPop = TRUE, lwd = 3.5, alpha = 1)

## ---- echo=FALSE--------------------------------------------------------------
####################
# Cleanup before next run
####################
unlink(x = outFolder, recursive = TRUE)
rm(list=ls())

## -----------------------------------------------------------------------------
####################
# Load libraries
####################
library(MGDrivE)

####################
# Output Folder
####################
outFolder <- "mgdrive"

####################
# Simulation Parameters
####################
# days to run the simulation, 2 years
tMax <- 365*2

# entomological parameters
bioParameters <- list(betaK=20, tEgg=5, tLarva=6, tPupa=4, popGrowth=1.175, muAd=0.09)

# a 1-node network where mosquitoes do not leave
moveMat <- as.matrix(1)

# parameters of the population equilibrium
adultPopEquilibrium <- 500
sitesNumber <- nrow(moveMat)
patchPops <- rep(adultPopEquilibrium,sitesNumber)

####################
# Basic Inheritance pattern
####################
# Reciprocal translocation cube with standard (default) parameters
cube <- cubeReciprocalTranslocations()

####################
# Setup releases and batch migration
####################
# set up the empty release vector
#  MGDrivE pulls things out by name
patchReleases <- replicate(n=sitesNumber,
                           expr={list(maleReleases=NULL,femaleReleases=NULL,
                                      eggReleases=NULL,matedFemaleReleases=NULL)},
                           simplify=FALSE)

# choose release parameters
#  Releases start at time 25, occur once a week, for 6 weeks.
#  There are 100 mosquitoes released every time.
releasesParameters <- list(releasesStart=25,
                          releasesNumber=6,
                          releasesInterval=7,
                          releaseProportion=100)

# generate release vector
releasesVector <- generateReleaseVector(driveCube=cube,
                                        releasesParameters=releasesParameters)

# put releases into the proper place in the release list
patchReleases[[1]]$maleReleases <- releasesVector
patchReleases[[1]]$femaleReleases <- releasesVector

# batch migration is disabled by setting the probability to 0
# This is required because of the stochastic simulations, but doesn't make sense
#  in a deterministic simulation.
batchMigration <- basicBatchMigration(batchProbs=0,
                                     sexProbs=c(.5,.5),
                                     numPatches=sitesNumber)

####################
# Combine parameters and run!
####################
# set MGDrivE to run deterministic
setupMGDrivE(stochasticityON = FALSE, verbose = FALSE)

# setup parameters for the network. This builds a list of parameters required for
#  every population in the network. In ths case, we have a network of 1 population.
netPar <- parameterizeMGDrivE(runID=1, simTime=tMax, nPatch=sitesNumber,
                             beta=bioParameters$betaK, muAd=bioParameters$muAd,
                             popGrowth=bioParameters$popGrowth, tEgg=bioParameters$tEgg,
                             tLarva=bioParameters$tLarva, tPupa=bioParameters$tPupa,
                             AdPopEQ=patchPops, inheritanceCube = cube)

# build network prior to run
MGDrivESim <- Network$new(params=netPar,
                         driveCube=cube,
                         patchReleases=patchReleases,
                         migrationMale=moveMat,
                         migrationFemale=moveMat,
                         migrationBatch=batchMigration,
                         directory=outFolder,
                         verbose = FALSE)
# run simulation
MGDrivESim$oneRun(verbose = FALSE)

####################
# Post Analysis
####################
# split output by patch
#  Required for plotting later
splitOutput(readDir = outFolder, remFile = TRUE, verbose = FALSE)

# aggregate females by their mate choice
#  This reduces the female file to have the same columns as the male file
aggregateFemales(readDir = outFolder, genotypes = cube$genotypesID,
                 remFile = TRUE, verbose = FALSE)

# plot output to see effect
plotMGDrivESingle(readDir = outFolder, totalPop = TRUE, lwd = 3.5, alpha = 1)

## ---- echo=FALSE--------------------------------------------------------------
####################
# Cleanup before next run
####################
unlink(x = outFolder, recursive = TRUE)
rm(list=ls())

## -----------------------------------------------------------------------------
####################
# Load libraries
####################
library(MGDrivE)

####################
# Output Folder
####################
# This is slightly obtuse for vignette building reasons
#  Really, all you need is a base directory, then the repetitions in folders inside that.
#  Here, our base directory is "mgdrive", and the repetition folders are "001","002", etc.
#  So, the final structure is "mgdrive/001","mgdrive/002", etc.
outFolder <- "mgdrive"
dir.create(path = outFolder)

####################
# Simulation Parameters
####################
# days to run the simulation, 3 years
tMax <- 365*3

# number of Monte Carlo iterations
nRep <- 20

# each Monte Carlo iteration gets its own folder
folderNames <- file.path(outFolder,
                        formatC(x = 1:nRep, width = 3, format = "d", flag = "0"))

# entomological parameters
bioParameters <- list(betaK=20, tEgg=5, tLarva=6, tPupa=4, popGrowth=1.175, muAd=0.09)

# a 1-node network where mosquitoes do not leave
moveMat <- as.matrix(1)

# parameters of the population equilibrium
adultPopEquilibrium <- 500
sitesNumber <- nrow(moveMat)

####################
# Basic Inheritance pattern
####################
# Reciprocal translocation cube with standard (default) parameters
cube <- cubeReciprocalTranslocations()

####################
# Setup releases and batch migration
####################
# set up the empty release vector
#  MGDrivE pulls things out by name
patchReleases <- replicate(n=sitesNumber,
                           expr={list(maleReleases=NULL,femaleReleases=NULL,
                                      eggReleases=NULL,matedFemaleReleases=NULL)},
                           simplify=FALSE)

# choose release parameters
#  Releases start at time 25, occur once a week, for 6 weeks.
#  There are 100 mosquitoes released every time.
releasesParameters <- list(releasesStart=25,
                          releasesNumber=6,
                          releasesInterval=7,
                          releaseProportion=100)

# generate male release vector
maleReleasesVector <- generateReleaseVector(driveCube=cube,
                                           releasesParameters=releasesParameters)

# generate female release vector
femaleReleasesVector <- generateReleaseVector(driveCube=cube,
                                             releasesParameters=releasesParameters)

# put releases into the proper place in the release list
patchReleases[[1]]$maleReleases <- maleReleasesVector
patchReleases[[1]]$femaleReleases <- femaleReleasesVector

# batch migration is disabled by setting the probability to 0
# This is required because of the stochastic simulations, but doesn't make sense
#  in a deterministic simulation.
batchMigration <- basicBatchMigration(batchProbs=0,
                                     sexProbs=c(.5,.5),
                                     numPatches=sitesNumber)

####################
# Combine parameters and run!
####################
# set MGDrivE to run stochastic
setupMGDrivE(stochasticityON = TRUE, verbose = FALSE)

# setup parameters for the network. This builds a list of parameters required for
#  every population in the network. In ths case, we have a network of 1 population.
netPar <- parameterizeMGDrivE(runID=1, simTime=tMax, nPatch=sitesNumber,
                             beta=bioParameters$betaK, muAd=bioParameters$muAd,
                             popGrowth=bioParameters$popGrowth, tEgg=bioParameters$tEgg,
                             tLarva=bioParameters$tLarva, tPupa=bioParameters$tPupa,
                             AdPopEQ=adultPopEquilibrium, inheritanceCube = cube)

# build network prior to run
MGDrivESim <- Network$new(params=netPar,
                         driveCube=cube,
                         patchReleases=patchReleases,
                         migrationMale=moveMat,
                         migrationFemale=moveMat,
                         migrationBatch=batchMigration,
                         directory=folderNames,
                         verbose = TRUE)
# run simulation
MGDrivESim$multRun(verbose = FALSE)

####################
# Post Analysis
####################
# First, split output by patch
# Second, aggregate females by their mate choice
for(i in 1:nRep){
 splitOutput(readDir = folderNames[i], remFile = TRUE, verbose = FALSE)
 aggregateFemales(readDir = folderNames[i], genotypes = cube$genotypesID,
                  remFile = TRUE, verbose = FALSE)
}

# plot output of first run to see effect
plotMGDrivESingle(readDir = folderNames[1], lwd = 3.5, alpha = 1)

# plot all 50 repetitions together
plotMGDrivEMult(readDir = outFolder, lwd = 0.35, alpha = 0.75)

## ---- echo=FALSE--------------------------------------------------------------
####################
# Theory Comparison
####################

# list male files
mFiles <- list.files(path = outFolder, recursive = TRUE, pattern = "^M.*.csv$",
                    full.names = TRUE)

# read in simulation files
successCount <- 0
for(f in mFiles){
 # read in files, one by one
 hold <- matrix(data = scan(file = f, what = double(), sep = ",", skip = 1, quiet = TRUE),
                nrow = tMax, ncol = 1 + cube$genotypesN, byrow = TRUE)

 # check if the simulation was successful
 successCount <- successCount + (hold[tMax,10]!=0)
}

# percentage of success
sCPerc <- successCount / nRep * 100

## ---- echo=FALSE--------------------------------------------------------------
####################
# Cleanup before next run
####################
unlink(x = outFolder, recursive = TRUE)
rm(list=ls())

