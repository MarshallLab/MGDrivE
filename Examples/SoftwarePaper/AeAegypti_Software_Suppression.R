###############################################################################################################
# ______  __________________       _____       __________
# ___   |/  /_  ____/__  __ \_________(_)__   ____  ____/
# __  /|_/ /_  / __ __  / / /_  ___/_  /__ | / /_  __/
# _  /  / / / /_/ / _  /_/ /_  /   _  / __ |/ /_  /___
# /_/  /_/  \____/  /_____/ /_/    /_/  _____/ /_____/
# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
# Population Suppression in the: "Tale of Two Cities" Spatial Scenario
###############################################################################################################
######################### LOAD AND SETUP PACKAGES #############################################################
rm(list=ls());gc()
library(foreach)
library(doMC)
library(MGDrivE)
registerDoMC(4)
###############################################################################################################
######################### SETUP PATHS #########################################################################
setupMGDrivE(stochasticityON=TRUE)
setwd("/Users/sanchez.hmsc/Downloads/MGDrivE-master/Examples/SoftwarePaper/")
outputDirectory="/Users/sanchez.hmsc/Desktop/SoftwarePaper/SuppressionStochastic/"
postAnalysisDirectory="/Users/sanchez.hmsc/Desktop/SoftwarePaper/SuppressionStochastic_Analyzed/"
###############################################################################################################
######################### SETUP SIMULATION ####################################################################
simulationTime=5000 # Number of "days" run in the simulation
repetitions=20    # Number of repetitions to run on each scenario (for stochastic version)
###############################################################################################################
######################### SETUP MOSQUITO BIOLOGY ##############################################################
#bioParameters=list(betaK=8,tEgg=6,tLarva=11,tPupa=4,popGrowth=1.175,muAd=.09)
bioParameters=list(betaK=2*10,tEgg=5,tLarva=6,tPupa=4,popGrowth=1.175,muAd=.09)
###############################################################################################################
######################### SETUP LANDSCAPE #####################################################################
distancesMatrix=as.matrix(read.csv("./GeoLandscapes/ATaleOfTwoCitiesScaled_Distances.csv",sep=",",header=FALSE))
#movementKernel=as.matrix(read.table("./GeoLandscapes/ATaleOfTwoCities_Kernel.csv",sep=",",header=TRUE))
#distancesMatrix=matrix(0)
lifespanStayProbability=.90
pulseHeight=lifespanStayProbability^(bioParameters$muAd)
movementKernel=calcHurdleExpKernel(distancesMatrix,MGDrivE::kernels$exp_rat,pulseHeight)
write.table(movementKernel,file="./GeoLandscapes/ATaleOfTwoCitiesScaled_Kernel.csv",row.names=FALSE,col.names=FALSE,sep=",")
sitesNumber=nrow(movementKernel)
patchPops=rep(50,sitesNumber)
###############################################################################################################
######################### SETUP GENE-DRIVE AND RELEASES #######################################################
sHet=.9
eM=0.999
eF=0.999
driveCube=cubeHomingDrive(cM = 1, cF = 1, chM = eM, crM = 1/3, chF = eF, crF = 1/3,
  s=c(
    "WW"=1, "WH"=1-sHet, "WR"=1, "WB"=1-sHet,
    "HH"=0, "HR"=1-sHet, "HB"=0,
    "RR"=1, "RB"=1-sHet,
    "BB"=0
  )
)
patchReleases=replicate(n=sitesNumber,expr={list(maleReleases=NULL,femaleReleases=NULL)},simplify=FALSE)
releasesParameters=list(
  releasesStart=100,releasesNumber=5,releasesInterval=2*(bioParameters$tEgg+bioParameters$tLarva+bioParameters$tPupa),
  releaseProportion=2*round(mean(patchPops))
)
maleReleasesVector=generateReleaseVector(driveCube=driveCube,releasesParameters=releasesParameters,sex="M")
for(i in 6:6){patchReleases[[i]]$maleReleases=maleReleasesVector}
batchMigration=basicBatchMigration(batchProbs=0,sexProbs=c(.5,.5),numPatches=sitesNumber)
###############################################################################################################
################################ PREPARE THE FOLDERS ##########################################################
folderNames=character(length = repetitions)
for(i in 1:repetitions){
  folderName=paste0(outputDirectory, formatC(x = i, width = 3, format = "d", flag = "0"))
  dir.create(folderName)
  folderNames[i]=folderName
}
###############################################################################################################
################################ RUN THE MODEL ################################################################
foreach(i=1:repetitions) %dopar% {
  netPar=parameterizeMGDrivE(
    runID=i,simTime=simulationTime,nPatch=sitesNumber,
    beta=bioParameters$betaK,muAd=bioParameters$muAd,popGrowth=bioParameters$popGrowth,
    tEgg=bioParameters$tEgg,tLarva=bioParameters$tLarva,tPupa=bioParameters$tPupa,
    AdPopEQ=patchPops
  )
  network=Network$new(
    params=netPar,
    driveCube=driveCube,
    patchReleases=patchReleases,
    migrationMale=movementKernel,
    migrationFemale=movementKernel,
    directory=folderNames[i],
    migrationBatch=batchMigration,
    verbose=TRUE
  )
  network$oneRun()
  network$reset()
}
###############################################################################################################
############################### POST-ANALYSIS #################################################################
for(i in 1:repetitions){
   splitOutput(readDir=folderNames[i])
   aggregateFemales(readDir=folderNames[i], genotypes=driveCube$genotypesID)
}
calcQuantiles(
  readDir=outputDirectory,
  writeDirectory=postAnalysisDirectory,
  mean=TRUE
)
###############################################################################################################
############################### PLOTS #################################################################
plotMGDrivESingle(readDir = folderName[1], whichPatches=1:10)
plotMGDrivEMult(readDir = outputDirectory)
