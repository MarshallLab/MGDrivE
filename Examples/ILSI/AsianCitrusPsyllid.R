###############################################################################################################
# ______  __________________       _____       __________
# ___   |/  /_  ____/__  __ \_________(_)__   ____  ____/
# __  /|_/ /_  / __ __  / / /_  ___/_  /__ | / /_  __/
# _  /  / / / /_/ / _  /_/ /_  /   _  / __ |/ /_  /___
# /_/  /_/  \____/  /_____/ /_/    /_/  _____/ /_____/
# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
###############################################################################################################
######################### LOAD AND SETUP PACKAGES #############################################################
rm(list=ls());gc()
library(MGDrivE)
###############################################################################################################
######################### SETUP PATHS #########################################################################
# 1. Change the working directory to the folder where the repository is stored
# 2. Setup an output folder for the results of the simulation to be stored
setupMGDrivE(stochasticityON=TRUE)
outputDirectory = "~/Desktop/OUTPUT/"
eraseDirectory(outputDirectory)
###############################################################################################################
######################### SETUP SIMULATION ####################################################################
simulationTime=2000 # Number of "days" run in the simulation
repetitions=2       # Number of repetitions to run on each scenario (for stochastic version)
###############################################################################################################
######################### SETUP LANDSCAPE #####################################################################
# 1: Close
# 2: Mid
# 3: Far
SCENARIO=1
if(SCENARIO==1){movementKernel=matrix(data=c(0.99999,0.0001,0.0001,0.9999),nrow=2,ncol=2,byrow=TRUE)}
if(SCENARIO==2){movementKernel=matrix(data=c(.999995,0.000005,0.000005,.999995),nrow=2,ncol=2,byrow=TRUE)}
if(SCENARIO==3){movementKernel=matrix(data=c(.99999999,0.00000001,0.00000001,.99999999),nrow=2,ncol=2,byrow=TRUE)}
adultPopEquilibrium=500
sitesNumber=nrow(movementKernel)
patchPops=rep(adultPopEquilibrium,sitesNumber)
###############################################################################################################
######################### SETUP MOSQUITO BIOLOGY ##############################################################
#Ref: Liu & Tsai 2000, Effects of temperature on biology and life table parameters of the Asian citrus psyllid...
# 1: Citrus Psyllid
# 2: Mosquito
# 3: Medfly
SPECIES=1
if(SPECIES==1){bioParameters=list(betaK=2,tEgg=10,tLarva=25,tPupa=14,popGrowth=1.036,muAd=0.02)}
if(SPECIES==2){bioParameters=list(betaK=16,tEgg=4,tLarva=3,tPupa=6,popGrowth=1.096,muAd=.09)}
if(SPECIES==3){bioParameters=list(betaK=20,tEgg=2,tLarva=6,tPupa=10,popGrowth=1.031,muAd=.1)}
###############################################################################################################
######################### SETUP GENE-DRIVE AND RELEASES #######################################################
scaler=1/1
sCost=list(
  "H"=1-scaler*.25,
  "R"=1-scaler*.25,
  "B"=1-scaler*.50,
  "W"=1-scaler*0.0
)
driveCube=cubeHomingDrive(chM = 0.9, crM = 1/3, chF = 0.5, crF = 1/3,
                          s=c(
                            'HH'=(sCost$H+sCost$H),
                            'WH'=(sCost$H+sCost$W),
                            'HR'=(sCost$H+sCost$R),
                            'HB'=(sCost$H+sCost$B),
                            'WW'=(sCost$W+sCost$W),
                            'WR'=(sCost$W+sCost$R),
                            'WB'=(sCost$W+sCost$B),
                            'RR'=(sCost$R+sCost$R),
                            'RB'=(sCost$R+sCost$B),
                            'BB'=(sCost$B+sCost$B)
                          )
)
#driveCube=cubeCRISPRTwoResistantAllele(s=c(1,1,1,1,1,1,1,1,1,1),e=.9,p1=(1/10 * 1/3),p2=(1/10 * 2/3))
patchReleases=replicate(n=sitesNumber,expr={list(maleReleases=NULL,femaleReleases=NULL)},simplify=FALSE)
releasesParameters=list(
  releasesStart=100,releasesNumber=1,releasesInterval=0,
  releaseProportion=round(mean(patchPops))
)
maleReleasesVector=generateReleaseVector(driveCube=driveCube,releasesParameters=releasesParameters,sex="M")
for(i in 1:1){patchReleases[[i]]$maleReleases=maleReleasesVector}
batchMigration=basicBatchMigration(batchProbs=0,sexProbs=c(.5,.5),numPatches=sitesNumber)
###############################################################################################################
################################ PREPARE THE FOLDERS ##########################################################
folderNames <- paste0(outputDirectory,formatC(x = 1:repetitions, width = 3, format = "d", flag = "0"))
###############################################################################################################
################################ RUN THE MODEL ################################################################
netPar=parameterizeMGDrivE(runID=1,simTime=simulationTime,nPatch=sitesNumber,
                           beta=bioParameters$betaK,muAd=bioParameters$muAd,
                           popGrowth=bioParameters$popGrowth,
                           tEgg=bioParameters$tEgg,tLarva=bioParameters$tLarva,
                           tPupa=bioParameters$tPupa, AdPopEQ=patchPops
)
MGDrivESim=Network$new(params=netPar,
                       driveCube=driveCube,
                       patchReleases=patchReleases,
                       migrationMale=movementKernel,
                       migrationFemale=movementKernel,
                       migrationBatch=batchMigration,
                       directory=folderNames
)
MGDrivESim$multRun()
###############################################################################################################
############################### POST-ANALYSIS #################################################################
for(i in 1:repetitions){
  splitOutput(readDir=folderNames[i])
  aggregateFemales(readDir=folderNames[i],genotypes=driveCube$genotypesID,remFile=TRUE)
}
###############################################################################################################
############################### PLOTS #################################################################
#png(filename= paste0(outputDirectory ,"/name.png"), width=5000, height=3000, pointsize=50)
plotMGDrivESingle(readDir=folderNames[1],lwd=.35,alpha=.75)
plotMGDrivEMult(readDir=outputDirectory,lwd=.35,alpha=.75)
#dev.off()
