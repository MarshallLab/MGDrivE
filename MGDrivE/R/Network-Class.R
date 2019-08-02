########################################################################
#       _   __     __                      __
#      / | / /__  / /__      ______  _____/ /__
#     /  |/ / _ \/ __/ | /| / / __ \/ ___/ //_/
#    / /|  /  __/ /_ | |/ |/ / /_/ / /  / ,<
#   /_/ |_/\___/\__/ |__/|__/\____/_/  /_/|_|
#
#   Network Class Definition
#   Marshall Lab
#   November 2017
#
########################################################################


########################################################################
# Class Definition
########################################################################

#' Network Class Definition
#'
#' A \code{Network} class object stores all the information for a simulation on a defined landscape.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#'
#' @section **Constructor**:
#'  * params: see \code{\link{parameterizeMGDrivE}}
#'  * driveCube: an inheritance cube
#'  * patchReleases: see \code{\link{basicRepeatedReleases}} for examples on how to set up release schedules
#'  * migrationMale: a stochastic matrix whose dimensions conform to the number of patches
#'  * migrationFemale: a stochastic matrix whose dimensions conform to the number of patches
#'  * directory: character string of output directory
#'
#' @section **Methods**:
#'  * get_moveVar: see \code{\link{get_moveVar_Network}}
#'  * get_timeAq: see \code{\link{get_timeAq_Network}}
#'  * get_thetaAq: see \code{\link{get_thetaAq_Network}}
#'  * get_windowSize: see \code{\link{get_windowSize_Network}}
#'  * get_beta: see \code{\link{get_beta_Network}}
#'  * get_muAd: see \code{\link{get_muAd_Network}}
#'  * get_dayPopGrowth: see \code{\link{get_dayPopGrowth_Network}}
#'  * get_AdPopEQ: see \code{\link{get_AdPopEQ_Network}}
#'  * get_g: see \code{\link{get_g_Network}}
#'  * get_genPopGrowth: see \code{\link{get_genPopGrowth_Network}}
#'  * get_muAq: see \code{\link{get_muAq_Network}}
#'  * get_alpha: see \code{\link{get_alpha_Network}}
#'  * get_Leq_Network: see \code{\link{get_Leq_Network}}
#'  * get_drivecubegenotype: see \code{\link{get_drivecubegenotype_Network}}
#'  * get_drivecubeindex: see \code{\link{get_drivecubeindex_Network}}
#'  * get_tau: see \code{\link{get_tau_Network}}
#'  * get_genotypesID: see \code{\link{get_genotypesID_Network}}
#'  * get_genotypesN: see \code{\link{get_genotypesN_Network}}
#'  * get_wildType: see \code{\link{get_wildType_Network}}
#'  * get_eta: see \code{\link{get_eta_Network}}
#'  * get_phi: see \code{\link{get_phi_Network}}
#'  * getOmega: see \code{\link{getOmega_Network}}
#'  * get_xiF: see \code{\link{get_xiF_Network}}
#'  * get_xiM: see \code{\link{get_xiM_Network}}
#'  * get_s: see \code{\link{get_s_Network}}
#'  * get_releaseType: see \code{\link{get_releaseType_Network}}
#'  * get_patch: see \code{\link{get_patch_Network}}
#'  * get_patches: see \code{\link{get_patches_Network}}
#'  * get_nPatch: see \code{\link{get_nPatch_Network}}
#'  * get_directory: see \code{\link{get_directory_Network}}
#'  * get_simTime: see \code{\link{get_simTime_Network}}
#'  * get_conADM: see \code{\link{get_conM_Network}}
#'  * get_conAF1: see \code{\link{get_conF_Network}}
#'  * close_allConnections: see \code{\link{close_allConnections_Network}}
#'  * get_tNow: see \code{\link{get_tNow_Network}}
#'  * get_migrationMale: see \code{\link{get_migrationMale_Network}}
#'  * get_migrationMaleRow: see \code{\link{get_migrationMaleRow_Network}}
#'  * set_migrationMale: see \code{\link{set_migrationMale_Network}}
#'  * get_migrationFemale: see \code{\link{get_migrationFemale_Network}}
#'  * get_migrationFemaleRow: see \code{\link{get_migrationFemaleRow_Network}}
#'  * set_migrationFemale: see \code{\link{set_migrationFemale_Network}}
#'  * get_patchReleases: see \code{\link{get_patchReleases_Network}}
#'  * oneDay_Migration: see \code{\link{oneDay_Migration_Network}}
#'  * reset: see \code{\link{reset_Network}}
#'  * oneRun: see \code{\link{oneRun_Network}}
#'  * multRun: see \code{\link{multRun_Network}}
#'  * oneDay: see \code{\link{oneDay_Network}}
#'
#' @section **Fields**:
#'  * params: see \code{\link{parameterizeMGDrivE}}
#'  * patches: a list of \code{\link{Patch}} objects
#'  * nPatch: number of patches
#'  * simTime: maximum time of simulation
#'  * driveCube: an inheritance cube
#'  * tNow: current time of simulation (time starts at 2 because time 1 is the initial equilibrium state)
#'  * runID: an identifier for the current simulation run, useful for Monte Carlo simulation
#'  * directory: a character string of where to store output
#'  * conADM: a \code{\link[base]{connection}} to write male population dynamics out to
#'  * conAF1: a \code{\link[base]{connection}} to write female population dynamics out to
#'  * migrationMale: a stochastic matrix whose dimensions conform to the number of patches
#'  * migrationFemale: a stochastic matrix whose dimensions conform to the number of patches
#'  * migrationBatch: list of items for batch migration in stochastic sim.
#'  * patchReleases: a list of release schedules for each patch
#'  * verbose: Chatty? Default is TRUE
#'
#'  @examples
#'  \dontrun{
#'  # There are no simple examples for this, so looking at the vignettes would be
#'  #  most useful.
#'
#'  # Complete manual with examples, but none explored in depth.
#'  vignette("MGDrivE-Examples", package = "MGDrivE")
#'
#'  # One example, explored in great detail. This is probably more helpful.
#'  vignette("MGDrivE-Run", package = "MGDrivE")
#'
#'  }
#'
#' @export
Network <- R6::R6Class(classname = "Network",
            portable = TRUE,
            cloneable = FALSE,
            lock_class = FALSE,
            lock_objects = FALSE,
            class = FALSE,

            # public memebers
            public = list(

                #################################################
                # Constructor
                #################################################

                initialize = function(params, driveCube, patchReleases,
                                      migrationMale, migrationFemale, migrationBatch = NULL,
                                      directory, verbose = TRUE){

                  if(length(patchReleases) != params$nPatch){
                    stop("length of patchReleases must equal number of patches in params!")
                  }

                  private$parameters = params
                  private$nPatch = params$nPatch
                  private$patches = vector(mode="list",length=private$nPatch)
                  private$simTime = params$simTime
                  private$driveCube = driveCube
                  private$directory = directory
                  private$runID = params$runID

                  private$migrationMale = migrationMale
                  private$migrationFemale = migrationFemale
                  private$migrationBatch = migrationBatch

                  private$patchReleases = patchReleases

                  # initialize patches
                  for(i in 1:private$nPatch){

                    # initialize patch i
                    if(verbose){cat("initializing patch: ", i, " of ",
                                    private$nPatch, "\n")}

                    # initial aquatic populations
                    EGGt0 = createNamedPopVector(driveCube$genotypesID)
                    EGGt0[driveCube$wildType] = params$Leq[i]

                    # initial adult male population
                    ADMt0 = createNamedPopVector(driveCube$genotypesID)
                    ADMt0[driveCube$wildType] = params$AdPopEQ[i]/2

                    # initial adult female population
                    AF1t0 = createNamedPopMatrix(driveCube$genotypesID)
                    AF1t0[driveCube$wildType,driveCube$wildType] = params$AdPopEQ[i]/2

                    # initialize patch
                    private$patches[[i]] = Patch$new(patchID = i,
                                                     genotypesID = driveCube$genotypesID,
                                                     simTime = private$simTime,
                                                     windowSize = params$windowSize,
                                                     EGGt0 = EGGt0,
                                                     LARt0 = EGGt0,
                                                     PUPt0 = EGGt0,
                                                     ADMt0 = ADMt0,
                                                     AF1t0 = AF1t0,
                                                     maleReleases = patchReleases[[i]]$maleReleases,
                                                     femaleReleases = patchReleases[[i]]$femaleReleases,
                                                     eggReleases = patchReleases[[i]]$eggReleases,
                                                     numPatches = private$nPatch
                                                   )

                    # set pointers
                    private$patches[[i]]$set_NetworkPointer(self)
                  }

                  # Output
                  if(!all(dir.exists(directory))){
                    for(f in directory){suppressWarnings(dir.create(f))}
                  } else {
                    # if running in serial remove files, else do nothing
                    if(!params$parallel){
                      dirFiles = list.files(path = directory)
                      if(length(dirFiles)>0){
                        if(verbose){cat("warning: ", length(dirFiles),
                            " files found in the output directory; please move files to avoid being overwritten\n",
                            sep="")}
                      }

                    }

                  }

                } # end constructor

              ),

            # private members
            private = list(

                parameters = NULL,
                patches = NULL, # list of patches
                nPatch = NULL, # number of patches
                simTime = NULL, # max time of sim
                driveCube = NULL, # number of genotypes in simulation
                tNow = 2L, # time starts at 2 because time 1 is the initial condition
                runID = numeric(1),

                # output
                directory = NULL, # directory to store all patch output
                conADM = NULL,
                conAF1 = NULL,

                # inter-patch migration
                migrationMale = NULL,
                migrationFemale = NULL,
                migrationBatch = NULL,

                # release schedule
                patchReleases = NULL

              )
)


########################################################################
# Getters & Setters: Parameters
########################################################################

#' Get moveVar
#'
#' Return numeric variance in Dirchlet-Multinomial movement
#'
get_moveVar_Network <- function(){return(private$parameters$moveVar)}

Network$set(which = "public",name = "get_moveVar",
  value = get_moveVar_Network,overwrite = TRUE
)

#' Get timeAq
#'
#' Return duration of aquatic stages, see \code{\link{initStagesDurations}}
#'
#' @param stage character in 'E', 'L', 'P'; if \code{NULL} return total duration
#'
get_timeAq_Network <- function(stage = NULL){
  if(is.null(stage)){
    return(sum(private$parameters$timeAq))
  } else {
    return(private$parameters$timeAq[[stage]])
  }
}

Network$set(which = "public",name = "get_timeAq",
  value = get_timeAq_Network,overwrite = TRUE
)

#' Get thetaAq
#'
#' Return aquatic stage survival probability, see \code{\link{calcAquaticStagesSurvivalProbability}} and \code{\link{calcAquaticStageSurvivalProbability}}
#'
#' @param stage character in 'E', 'L', 'P'
#'
get_thetaAq_Network <- function(stage){
  return(private$parameters$thetaAq[[stage]])
}

Network$set(which = "public",name = "get_thetaAq",
  value = get_thetaAq_Network,overwrite = TRUE
)

#' Get windowSize
#'
#' Return memory window size, see \code{\link{calcMemoryWindow}}
#'
get_windowSize_Network <- function(){return(private$parameters$windowSize)}

Network$set(which = "public",name = "get_windowSize",
  value = get_windowSize_Network,overwrite = TRUE
)

#' Get beta
#'
#' Return size of wild-type egg batch
#'
get_beta_Network <- function(){return(private$parameters$beta)}

Network$set(which = "public",name = "get_beta",
  value = get_beta_Network,overwrite = TRUE
)

#' Get muAd
#'
#' Return adult mortality
#'
get_muAd_Network <- function(){return(private$parameters$muAd)}

Network$set(which = "public",name = "get_muAd",
  value = get_muAd_Network,overwrite = TRUE
)

#' Get dayPopGrowth
#'
#' Return daily population growth rate (rm)
#'
get_dayPopGrowth_Network <- function(){return(private$parameters$dayPopGrowth)}

Network$set(which = "public",name = "get_dayPopGrowth",
  value = get_dayPopGrowth_Network,overwrite = TRUE
)

#' Get AdPopEQ
#'
#' Return equilibrium adult population
#'
#' @param ix index of patch
#'
get_AdPopEQ_Network <- function(ix){return(private$parameters$AdPopEQ[ix])}

Network$set(which = "public",name = "get_AdPopEQ",
  value = get_AdPopEQ_Network,overwrite = TRUE
)

#' Get g
#'
#' Return average generation time, see \code{\link{calcAverageGenerationTime}}
#'
get_g_Network <- function(){return(private$parameters$g)}

Network$set(which = "public",name = "get_g",
  value = get_g_Network,overwrite = TRUE
)

#' Get genPopGrowth
#'
#' Return population growth rate, see \code{\link{calcPopulationGrowthRate}}
#'
get_genPopGrowth_Network <- function(){return(private$parameters$genPopGrowth)}

Network$set(which = "public",name = "get_genPopGrowth",
  value = get_genPopGrowth_Network,overwrite = TRUE
)

#' Get muAq
#'
#' Return larval mortality, see \code{\link{calcLarvalStageMortalityRate}}
#'
get_muAq_Network <- function(){return(private$parameters$muAq)}

Network$set(which = "public",name = "get_muAq",
  value = get_muAq_Network,overwrite = TRUE
)

#' Get alpha
#'
#' Return density dependent mortality, see \code{\link{calcDensityDependentDeathRate}}
#'
#' @param ix index of patch
#'
get_alpha_Network <- function(ix){return(private$parameters$alpha[ix])}

Network$set(which = "public",name = "get_alpha",
  value = get_alpha_Network,overwrite = TRUE
)

#' Get Leq
#'
#' Return equilibrium larval population, see \code{\link{calcLarvalPopEquilibrium}}
#'
#' @param ix index of patch
#'
get_Leq_Network <- function(ix){return(private$parameters$Leq[ix])}

Network$set(which = "public",name = "get_Leq_Network",
  value = get_Leq_Network,overwrite = TRUE
)

########################################################################
# Getters & Setters: Drive Cube
########################################################################

#' Get Element(s) of Drive Cube by Genotype
#'
#' Return elements or slices of drive cube. If all \code{NULL} return entire cube.
#'
#' @param fG female genotype
#' @param mG male genotype
#' @param oG offspring genotype
#'
get_drivecubegenotype_Network <- function(fG=NULL,mG=NULL,oG=NULL){
  if(is.null(fG)){fG = private$driveCube$genotypesID}
  if(is.null(mG)){mG = private$driveCube$genotypesID}
  if(is.null(oG)){oG = private$driveCube$genotypesID}
  return(private$driveCube$ih[fG,mG,oG])
}

Network$set(which = "public",name = "get_drivecubegenotype",
  value = get_drivecubegenotype_Network,overwrite = TRUE
)

#' Get Element(s) of Drive Cube by Index
#'
#' Return elements or slices of drive cube. If all \code{NULL} return entire cube.
#'
#' @param fG female genotype index
#' @param mG male genotype index
#' @param oG offspring genotype index
#'
get_drivecubeindex_Network <- function(fG=NULL,mG=NULL,oG=NULL){
  if(is.null(fG)){fG = 1:private$driveCube$genotypesN}
  if(is.null(mG)){mG = 1:private$driveCube$genotypesN}
  if(is.null(oG)){oG = 1:private$driveCube$genotypesN}
  return(private$driveCube$ih[fG,mG,oG])
}

Network$set(which = "public",name = "get_drivecubeindex",
  value = get_drivecubeindex_Network,overwrite = TRUE
)

#' Get Female Viability Mask (tau)
#'
#' @param fG Number for which female genotype to get
#' @param mG Number for which male genotype to get
#' @param oG Number for which offspring genotype to get
#'
#' Return matrix
#'
get_tau_Network <- function(fG=NULL,mG=NULL,oG=NULL){
  if(is.null(fG)){fG = 1:private$driveCube$genotypesN}
  if(is.null(mG)){mG = 1:private$driveCube$genotypesN}
  if(is.null(oG)){oG = 1:private$driveCube$genotypesN}
  return(private$driveCube$tau[fG,mG,oG])
}

Network$set(which = "public",name = "get_tau",
  value = get_tau_Network,overwrite = TRUE
)

#' Get genotypesID
#'
#' Return character vector of possible genotypes
#'
get_genotypesID_Network <- function(){return(private$driveCube$genotypesID)}

Network$set(which = "public",name = "get_genotypesID",
  value = get_genotypesID_Network,overwrite = TRUE
)

#' Get genotypesN
#'
#' Return number of possible genotypes
#'
get_genotypesN_Network <- function(){return(private$driveCube$genotypesN)}

Network$set(which = "public",name = "get_genotypesN",
  value = get_genotypesN_Network,overwrite = TRUE
)

#' Get wildType
#'
#' Return wild-type genotype
#'
get_wildType_Network <- function(){return(private$driveCube$wildType)}

Network$set(which = "public",name = "get_wildType",
  value = get_wildType_Network,overwrite = TRUE
)

#' Get eta
#'
#' Return genotype-specific mating fitness
#'
get_eta_Network <- function(){return(private$driveCube$eta)}

Network$set(which = "public",name = "get_eta",
  value = get_eta_Network,overwrite = TRUE
)

#' Get phi
#'
#' Return genotype-specific sex ratio at emergence
#'
get_phi_Network <- function(){return(private$driveCube$phi)}

Network$set(which = "public",name = "get_phi",
  value = get_phi_Network,overwrite = TRUE
)

#' Get omega
#'
#' Return genotype-specific multiplicative modifier of adult mortality
#'
getOmega_Network <- function(){return(private$driveCube$omega)}

Network$set(which = "public",name = "getOmega",
  value = getOmega_Network,overwrite = TRUE
)

#' Get xiF
#'
#' Return genotype-specific female pupatory success
#'
get_xiF_Network <- function(){return(private$driveCube$xiF)}

Network$set(which = "public",name = "get_xiF",
  value = get_xiF_Network,overwrite = TRUE
)

#' Get xiM
#'
#' Return genotype-specific male pupatory success
#'
get_xiM_Network <- function(){return(private$driveCube$xiM)}

Network$set(which = "public",name = "get_xiM",
  value = get_xiM_Network,overwrite = TRUE
)

#' Get s
#'
#' Return genotype-specific fractional reduction(increase) in fertility
#'
get_s_Network <- function(){return(private$driveCube$s)}

Network$set(which = "public",name = "get_s",
  value = get_s_Network,overwrite = TRUE
)

#' Get releaseType
#'
#' Return genotype of release
#'
get_releaseType_Network <- function(){return(private$driveCube$releaseType)}

Network$set(which = "public",name = "get_releaseType",
  value = get_releaseType_Network,overwrite = TRUE
)


########################################################################
# Getters & Setters: Other
########################################################################

#' Get Patch
#'
#' Return a \code{\link{Patch}} object
#'
#' @param ix integer id of patch to return
#'
get_patch_Network <- function(ix){return(private$patches[[ix]])}

Network$set(which = "public",name = "get_patch",
  value = get_patch_Network,overwrite = TRUE
)

#' Get all Patches
#'
#' Return a list of \code{\link{Patch}} objects
#'
get_patches_Network <- function(){return(private$patches)}

Network$set(which = "public",name = "get_patches",
  value = get_patches_Network,overwrite = TRUE
)

#' Get nPatch
#'
#' Return number of patches
#'
get_nPatch_Network <- function(){return(private$nPatch)}

Network$set(which = "public",name = "get_nPatch",
  value = get_nPatch_Network,overwrite = TRUE
)

#' Get directory
#'
#' Return character string of directory being written to
#'
get_directory_Network <- function(){return(private$directory)}

Network$set(which = "public",name = "get_directory",
  value = get_directory_Network,overwrite = TRUE
)

#' Get simTime
#'
#' Return maximum time to run simulation
#'
get_simTime_Network <- function(){return(private$simTime)}

Network$set(which = "public",name = "get_simTime",
  value = get_simTime_Network,overwrite = TRUE
)

#' Get conADM
#'
#' Return \code{\link[base]{connection}} where adult male dynamics are written to
#'
get_conM_Network <- function(){return(private$conADM)}

Network$set(which = "public",name = "get_conADM",
  value = get_conM_Network,overwrite = TRUE
)

#' Get conAF1
#'
#' Return \code{\link[base]{connection}} where adult female dynamics are written to
#'
get_conF_Network <- function(){return(private$conAF1)}

Network$set(which = "public",name = "get_conAF1",
  value = get_conF_Network,overwrite = TRUE
)

#' Close all Output Connections
#'
#' Close \code{private$conADM} and \code{private$conAF1}
#'
close_allConnections_Network <- function(){
  close(private$conADM)
  close(private$conAF1)
}

Network$set(which = "public",name = "close_allConnections",
  value = close_allConnections_Network,overwrite = TRUE
)

#' Get tNow
#'
#' Return current simulation time
#'
get_tNow_Network <- function(){return(private$tNow)}

Network$set(which = "public",name = "get_tNow",
  value = get_tNow_Network,overwrite = TRUE
)

# migration

#' Get Male Migration Matrix
#'
#' Return a matrix object
#'
get_migrationMale_Network <- function(){return(private$migrationMale)}

Network$set(which = "public",name = "get_migrationMale",
  value = get_migrationMale_Network,overwrite = TRUE
)

#' Get Row of Male Migration Matrix
#'
#' Return a matrix object (does not drop dimensions)
#'
#' @param ix index of row
#'
get_migrationMaleRow_Network <- function(ix){return(private$migrationMale[ix,,drop=FALSE])}

Network$set(which = "public",name = "get_migrationMaleRow",
  value = get_migrationMaleRow_Network,overwrite = TRUE
)

#' Set Male Migration Matrix
#'
#' Sets \code{link{private$migrationMale}} field
#'
#' @param migrationMale matrix object (rows must sum to one)
#'
set_migrationMale_Network <- function(migrationMale){private$migrationMale = migrationMale}

Network$set(which = "public",name = "set_migrationMale",
  value = set_migrationMale_Network,overwrite = TRUE
)

#' Get Female Migration Matrix
#'
#' Return a matrix object
#'
get_migrationFemale_Network <- function(){return(private$migrationFemale)}

Network$set(which = "public",name = "get_migrationFemale",
  value = get_migrationFemale_Network,overwrite = TRUE
)

#' Get Row of Female Migration Matrix
#'
#' Return a matrix object (does not drop dimensions)
#'
#' @param ix index of row
#'
get_migrationFemaleRow_Network <- function(ix){return(private$migrationFemale[ix,,drop=FALSE])}

Network$set(which = "public",name = "get_migrationFemaleRow",
  value = get_migrationFemaleRow_Network,overwrite = TRUE
)

#' Set Female Migration Matrix
#'
#' Sets \code{link{private$migrationFemale}} field
#'
#' @param migrationFemale matrix object (rows must sum to one)
#'
set_migrationFemale_Network <- function(migrationFemale){private$migrationFemale = migrationFemale}

Network$set(which = "public",name = "set_migrationFemale",
  value = set_migrationFemale_Network,overwrite = TRUE
)

#' Get Patch Release Schedule
#'
#' Return the release schedule for a patch for male or female
#'
#' @param ix index of patch
#' @param sex character in 'M', 'F'
#'
get_patchReleases_Network <- function(ix, sex = "M"){
  switch(sex,
    M = {return(private$patchReleases[[ix]]$maleReleases)},
    F = {return(private$patchReleases[[ix]]$femaleReleases)}
  )
}

Network$set(which = "public",name = "get_patchReleases",
  value = get_patchReleases_Network,overwrite = TRUE
)

#' Get Batch Migration Probability
#'
#' Return the probability of undergoing batch migration each day
#'
#' @param ix index of patch
#'
get_batchProbs_Network <- function(ix){return(private$migrationBatch[[1L]][ix])}

Network$set(which = "public",name = "get_batchProbs",
            value = get_batchProbs_Network,overwrite = TRUE
)

#' Get Batch Migration Sex
#'
#' Return the batch migration size for each sex
#'
#' @param ix index of patch
#' @param sex number, M=1 and F=2
#'
get_batchSex_Network <- function(ix, sex){return(private$migrationBatch[[2]][ix,sex, ])}

Network$set(which = "public",name = "get_batchSex",
            value = get_batchSex_Network,overwrite = TRUE
)

#' Get Batch Location Distribution
#'
#' Return the distribution of location probabilities
#'
#' @param ix index of patch
#'
get_batchLocRow_Network <- function(ix){return(private$migrationBatch[[3]][ix, ,drop=FALSE])}

Network$set(which = "public",name = "get_batchLocations",
            value = get_batchLocRow_Network,overwrite = TRUE
)

