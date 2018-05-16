###############################################################################################################
# ______  __________________       _____       __________
# ___   |/  /_  ____/__  __ \_________(_)__   ____  ____/
# __  /|_/ /_  / __ __  / / /_  ___/_  /__ | / /_  __/
# _  /  / / / /_/ / _  /_/ /_  /   _  / __ |/ /_  /___
# /_/  /_/  \____/  /_____/ /_/    /_/  _____/ /_____/
# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# Population Replacements in the: "Tale of Two Cities" Spatial Scenario
###############################################################################################################
######################### LOAD AND SETUP PACKAGES #############################################################
rm(list=ls());gc()
library(stringr)
library(MGDrivE)
###############################################################################################################
######################### SETUP PATHS #########################################################################
MGDrivE.Setup(stochasticityON=FALSE)
setwd("/Users/sanchez.hmsc/Documents/GitHub/MGDrivE_Releases/Examples/")
outputDirectory="/Users/sanchez.hmsc/Documents/GitHub/MGDrivE_Releases/Examples/OUTPUT/ReplacementDeterministic/"
###############################################################################################################
######################### SETUP SIMULATION ####################################################################
simulationTime=3000 # Number of "days" run in the simulation
repetitions=1       # Number of repetitions to run on each scenario (for stochastic version)
###############################################################################################################
######################### SETUP MOSQUITO BIOLOGY ##############################################################
bioParameters=list(betaK=2*10,tEgg=5,tLarva=6,tPupa=4,popGrowth=1.175,muAd=.09)
###############################################################################################################
######################### SETUP LANDSCAPE #####################################################################
distancesMatrix=as.matrix(read.csv("./GeoLandscapes/ATaleOfTwoCities_Distances.csv",sep=",",header=FALSE))
lifespanStayProbability=.90
pulseHeight=lifespanStayProbability^(bioParameters$muAd)
movementKernel=calc_HurdleExpKernel(distancesMatrix,MGDrivE::kernels$exp_rat,pulseHeight)
sitesNumber=nrow(movementKernel)
patchPops=rep(50,sitesNumber)
###############################################################################################################
######################### SETUP GENE-DRIVE AND RELEASES #######################################################
sH=.25
sR=.25
sB=.50
eM=0.9
eF=0.5
driveCube=Cube_HomingDrive(
  eM=eM,eF=eF,
  rM=1/3*(1-eM), bM=2/3*(1-eM),
  rF=1/3*(1-eF), bF=1/3*(1-eF),
  s=c(
    "WW"=1,"WH"=1-sH,"WR"=1-sR,"WB"=1-sB,
    "HH"=1-2*sH,"HR"=1-sH-sR,"HB"=1-sH-sB,
    "RR"=1-2*sR,"RB"=1-sR-sB,
    "BB"=1-2*sB
  ),
  eta=c(
    "WW"=1,"WH"=1-sH,"WR"=1-sR,"WB"=1-sB,
    "HH"=1-2*sH,"HR"=1-sH-sR,"HB"=1-sH-sB,
    "RR"=1-2*sR,"RB"=1-sR-sB,
    "BB"=1-2*sB
  )
)
patchReleases=replicate(n=sitesNumber,expr={list(maleReleases=NULL,femaleReleases=NULL)},simplify=FALSE)
releasesParameters=list(
  releasesStart=100,releasesNumber=2,releasesInterval=1000*(bioParameters$tEgg+bioParameters$tLarva+bioParameters$tPupa),
  releaseProportion=2*round(mean(patchPops))
)
maleReleasesVector=generateReleaseVector(driveCube=driveCube,releasesParameters=releasesParameters,sex="M")
for(i in 6:6){patchReleases[[i]]$maleReleases=maleReleasesVector}
###############################################################################################################
################################ PREPARE THE FOLDERS ##########################################################
folderNames=list()
for(i in 1:repetitions){
  folderName=paste0(outputDirectory,str_pad(i,4,"left","0"))
  dir.create(folderName)
  folderNames=c(folderNames,folderName)
}
###############################################################################################################
################################ RUN THE MODEL ################################################################
for(i in 1:repetitions){
  outputFolder=folderNames[[i]]
  netPar=Network.Parameters(
    runID=i,simTime=simulationTime,nPatch=sitesNumber,
    beta=bioParameters$betaK,muAd=bioParameters$muAd,popGrowth=bioParameters$popGrowth,
    tEgg=bioParameters$tEgg,tLarva=bioParameters$tLarva,tPupa=bioParameters$tPupa,
    AdPopEQ=patchPops
  )
  network=Network$new(
    networkParameters=netPar,
    driveCube=driveCube,
    patchReleases=patchReleases,
    migrationMale=movementKernel,
    migrationFemale=movementKernel,
    directory=outputFolder
  )
  network$oneRun()
  network$reset()
}
###############################################################################################################
############################### POST-ANALYSIS #################################################################
for(i in 1:repetitions){
  splitOutput(directory=folderNames[[i]])
  aggregateFemales(folderNames[[i]],driveCube$genotypesID,remove=FALSE)
}
