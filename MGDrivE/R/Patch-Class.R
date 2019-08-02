########################################################################
#        ____        __       __
#       / __ \____ _/ /______/ /_
#      / /_/ / __ `/ __/ ___/ __ \
#     / ____/ /_/ / /_/ /__/ / / /
#    /_/    \__,_/\__/\___/_/ /_/
#
#   Patch Class Definition
#   Marshall Lab
#   November 2017
#
########################################################################


#' Patch Class Definition
#'
#' A Patch is a single well-mixed population that is the smallest unit of simulation for MGDrivE.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#'
#' @importFrom stats rbinom rmultinom rpois
#'
#' @section **Constructor**:
#'  * patchID: integer ID of this patch
#'  * genotypesID: character vector of genotypes
#'  * simTime: maximum time of simulation
#'  * windowSize: necessary memory window size for model
#'  * EGGt0: initial egg population, \eqn{L_{eq}}
#'  * LARt0: initial larval population
#'  * PUPt0: initial pupae population
#'  * ADMt0: initial adult male population, \eqn{Ad_{eq}}
#'  * AF1t0: initial adult female population, \eqn{Ad_{eq}}
#'  * maleReleases: integer ID of this patch
#'  * femaleReleases: female release schedule for this patch, see \code{\link{basicRepeatedReleases}}
#'  * larvaeReleases: male release schedule for this patch, see \code{\link{basicRepeatedReleases}}
#'
#' @section **Methods**:
#'  * get_patchID: see \code{\link{get_patchID_Patch}}
#'  * get_AF1new: see \code{\link{get_AF1new_Patch}}
#'  * set_AF1new: see \code{\link{set_AF1new_Patch}}
#'  * get_ADMnew: see \code{\link{get_ADMnew_Patch}}
#'  * set_ADMnew: see \code{\link{set_ADMnew_Patch}}
#'  * accumulate_ADMnew: see \code{\link{accumulate_ADMnew_Patch}}
#'  * get_EGG: see \code{\link{get_EGG_Patch}}
#'  * get_LAR: see \code{\link{get_LAR_Patch}}
#'  * get_PUP: see \code{\link{get_PUP_Patch}}
#'  * get_ADM: see \code{\link{get_M_Patch}}
#'  * get_AF1: see \code{\link{get_F_Patch}}
#'  * get_EGGdly: see \code{\link{get_EGGdly_Patch}}
#'  * get_LARdly: see \code{\link{get_LARdly_Patch}}
#'  * get_PUPdly: see \code{\link{get_PUPdly_Patch}}
#'  * get_ADMdly: see \code{\link{get_ADMdly_Patch}}
#'  * get_AF1dly: see \code{\link{get_AF1dly_Patch}}
#'  * get_maleMigration: see \code{\link{get_maleMigration_Patch}}
#'  * get_femaleMigration: see \code{\link{get_femaleMigration_Patch}}
#'  * set_NetworkPointer: see \code{\link{set_NetworkPointer_Patch}}
#'  * get_NetworkPointer: see \code{\link{get_NetworkPointer_Patch}}
#'  * reset: see \code{\link{reset_Patch}}
#'  * oneDay_initOutput: see \code{\link{oneDay_initOutput_Patch}}
#'  * oneDay_writeOutput: see \code{\link{oneDay_writeOutput_Patch}}
#'  * oneDay_migrationIn: see \code{\link{oneDay_migrationIn_Patch}}
#'  * oneDay_maleReleases: see \code{\link{oneDay_maleReleases_Patch}}
#'  * oneDay_femaleReleases: see \code{\link{oneDay_femaleReleases_Patch}}
#'  * oneDay_PopDynamics: see \code{\link{oneDay_PopDynamics_Patch}}
#'  * oneDay_updatePopulation: see \code{\link{oneDay_updatePopulation_Patch}}
#'  * oneDay_calcLarvalDensityDependentFactor: see \code{\link{oneDay_calcLarvalDensityDependentFactor_Patch}}
#'  * oneDay_calcCumulativeLarvalDensityDependentFactor: see \code{\link{oneDay_calcCumulativeLarvalDensityDependentFactor_Patch}}
#'  * oneDay_calcCumulativePupaDensityDependentFactor: see \code{\link{oneDay_calcCumulativePupaDensityDependentFactor_Patch}}
#'  * oneDay_migrationOut: see \code{\link{oneDay_migrationOut_stochastic_Patch}} or \code{\link{oneDay_migrationOut_deterministic_Patch}}
#'  * oneDay_ovipositG1: see \code{\link{oneDay_ovipositG1_stochastic_Patch}} or \code{\link{oneDay_ovipositG1_deterministic_Patch}}
#'  * oneDay_larSurvival: see \code{\link{oneDay_larSurvival_stochastic_Patch}} or \code{\link{oneDay_larSurvival_deterministic_Patch}}
#'  * oneDay_hatchingFract: see \code{\link{oneDay_hatchingFract_stochastic_Patch}} or \code{\link{oneDay_hatchingFract_deterministic_Patch}}
#'  * oneDay_larHatching: see \code{\link{oneDay_larHatching_stochastic_Patch}} or \code{\link{oneDay_larHatching_deterministic_Patch}}
#'  * oneDay_eggsFract2: see \code{\link{oneDay_eggsFract2_stochastic_Patch}} or \code{\link{oneDay_eggsFract2_deterministic_Patch}}
#'  * oneDay_larPupating: see \code{\link{oneDay_larPupating_stochastic_Patch}} or \code{\link{oneDay_larPupating_deterministic_Patch}}
#'  * oneDay_numMaleFemale: see \code{\link{oneDay_numMaleFemale_stochastic_Patch}} or \code{\link{oneDay_numMaleFemale_deterministic_Patch}}
#'  * oneDay_admSurvival: see \code{\link{oneDay_admSurvival_stochastic_Patch}} or \code{\link{oneDay_admSurvival_deterministic_Patch}}
#'  * oneDay_admPupating: see \code{\link{oneDay_admPupating_stochastic_Patch}} or \code{\link{oneDay_admPupating_deterministic_Patch}}
#'  * oneDay_af1Survival: see \code{\link{oneDay_af1Survival_stochastic_Patch}} or \code{\link{oneDay_af1Survival_deterministic_Patch}}
#'  * oneDay_af1Pupation: see \code{\link{oneDay_af1Pupation_stochastic_Patch}} or \code{\link{oneDay_af1Pupation_deterministic_Patch}}
#'  * oneDay_af1Mating: see \code{\link{oneDay_af1Mating_stochastic_Patch}} or \code{\link{oneDay_af1Mating_deterministic_Patch}}
#'
#' @section **Fields**:
#'  * patchID: integer ID of this patch
#'  * EGGt0: vector of initial egg stage population
#'  * LARt0: vector of initial larval stage population
#'  * PUPt0: vector of initial pupae stage population
#'  * ADMt0: vector of initial adult male stage population
#'  * AF1t0: matrix of initial adult female stage population
#'  * EGG: egg stage population
#'  * LAR: larvae stage population
#'  * PUP: pupae stage population
#'  * ADM: adult male stage population
#'  * AF1: adult female stage population
#'  * EGGdly: delay egg stage population
#'  * LARdly: delay larvae stage population
#'  * PUPdly: delay pupae stage population
#'  * ADMdly: delay adult male stage population
#'  * AF1dly: delay adult female stage population
#'  * LARnew: new larval population after difference equations; needed to store population prior to migration exchange
#'  * ADMnew: new adult male population after difference equations; needed to store population prior to migration exchange
#'  * AF1new: new adult female population after difference equations; needed to store population prior to migration exchange
#'  * maleMigration: matrix of outbound migrating males of dimension nGenotypes X nPatch
#'  * femaleMigration: array of outbound migrating females of dimension nGenotypes X nGenotypes X nPatch
#'  * NetworkPointer: a reference to enclosing \code{\link{Network}}
#'  * ovipositG1: new eggs after oviposition by mated female mosquitoes
#'  * larSurvival: surviving larvae
#'  * hatchingFract: fraction of larvae that hatch
#'  * larPupating: fraction of larvae that undergo pupation
#'  * numMaleFemale: number of male vs. female emerging imago stage adults
#'  * admSurvival: number of surviving adult males
#'  * admPupating: number of pupating imago stage adults that become males
#'  * af1Survival: number of surviving adult females
#'  * af1Pupation: number of pupating imago stage adults that become females
#'  * maleMatrix: row of male migration matrix corresponding to migration from this patch
#'  * femaleMatrix: row of female migration matrix corresponding to migration from this patch
#'  * larDDMortal: larval mortality
#'  * f: density dependent factor in larval mortality
#'
Patch <- R6::R6Class(classname = "Patch",
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

                initialize = function(patchID, genotypesID, simTime, windowSize,
                                      EGGt0, LARt0, PUPt0, ADMt0, AF1t0,
                                      maleReleases = NULL, femaleReleases = NULL,
                                      eggReleases = NULL, numPatches){

                  # ID of this patch
                  private$patchID = patchID

                  # Population Array
                  private$EGG = initPopVectorArray(genotypesID,simTime)
                  private$LAR = initPopVectorArray(genotypesID,simTime)
                  private$PUP = initPopVectorArray(genotypesID,simTime)
                  private$ADM = initPopVectorArray(genotypesID,simTime)
                  private$AF1 = initPopMatrixArray(genotypesID,simTime)

                  #initialize objects for simulation. This way, they have dimensions and names done once.
                  private$ovipositG1 = setNames(object = numeric(length = length(genotypesID)), nm = genotypesID)
                  private$hatchingFract = setNames(object = numeric(length = length(genotypesID)), nm = genotypesID)
                  private$eggsFract2 = setNames(object = numeric(length = length(genotypesID)), nm = genotypesID)
                  private$numMaleFemale = matrix(data = 0, nrow = 2, ncol = length(genotypesID), dimnames = list(c("M","F"), genotypesID))
                  private$maleMigration = matrix(data=0,nrow=length(genotypesID),ncol=numPatches)
                  private$femaleMigration = array(data = 0,dim=c(length(genotypesID),length(genotypesID),numPatches))

                  # Population Initial Sizes
                  private$EGGt0 = EGGt0
                  private$LARt0 = LARt0
                  private$PUPt0 = PUPt0
                  private$ADMt0 = ADMt0
                  private$AF1t0 = AF1t0

                  private$EGG[[1]] = EGGt0
                  private$LAR[[1]] = LARt0
                  private$PUP[[1]] = PUPt0
                  private$ADM[[1]] = ADMt0
                  private$AF1[[1]] = AF1t0

                  # Windowed Population Tensors
                  private$EGGdly =  primePopVectorArray(private$EGG[[1]],windowSize)
                  private$LARdly =  primePopVectorArray(private$LAR[[1]],windowSize)
                  private$PUPdly =  primePopVectorArray(private$PUP[[1]],windowSize)
                  private$ADMdly =  primePopVectorArray(private$ADM[[1]],windowSize)
                  private$AF1dly =  primePopMatrixArray(private$AF1[[1]],windowSize)

                  # Mosquito Releases
                  private$maleReleases = maleReleases
                  private$femaleReleases = femaleReleases
                  private$eggReleases = eggReleases

                } # end constructor
              ),

            # private members
            private = list(

              patchID = NULL,

              # initial populations
              EGGt0 = NULL,
              LARt0 = NULL,
              PUPt0 = NULL,
              ADMt0 = NULL,
              AF1t0 = NULL,

              # populations
              EGG = NULL,
              LAR = NULL,
              PUP = NULL,
              ADM = NULL,
              AF1 = NULL,

              # delay populations
              EGGdly = NULL,
              LARdly = NULL,
              PUPdly = NULL,
              ADMdly = NULL,
              AF1dly = NULL,

              # temporary populations
              LARnew = NULL, # new populations after difference equations; needed to store pops prior to migration exchange
              ADMnew = NULL,
              AF1new = NULL,

              # migration
              maleMigration = NULL, # nGenotypes X nPatch matrix
              femaleMigration = NULL, # nGenotypes X nGenotypes X nPatch array

              # pointers
              NetworkPointer = NULL,

              # calculated during simulation
              ovipositG1 = NULL,
              larSurvival = NULL,
              hatchingFract = NULL,
              larPupating = NULL,
              numMaleFemale = NULL,
              admSurvival = NULL,
              admPupating = NULL,
              af1Survival = NULL,
              af1Pupation = NULL,
              maleMatrix = NULL, # migration probabilities
              femaleMatrix = NULL, # migration probabilities
              larDDMortal = numeric(1),
              f = numeric(1)
            )
)


########################################################################
# Getters & Setters
########################################################################

#' Get patchID
#'
#' Return the ID of this patch
#'
get_patchID_Patch <- function(){return(private$patchID)}

Patch$set(which = "public",name = "get_patchID",
  value = get_patchID_Patch,overwrite = TRUE
)

# populations

#' Get AF1new
#'
#' Return the new AF1 females
#'
get_AF1new_Patch <- function(){return(private$AF1new)}

Patch$set(which = "public",name = "get_AF1new",
  value = get_AF1new_Patch,overwrite = TRUE
)

#' Set AF1new
#'
#' Set an element in new AF1 females
#'
#' @param count number of mosquitoes
#' @param genotype_F genotype of female (row of matrix)
#' @param genotype_M genotype of male (column of matrix)
#'
set_AF1new_Patch <- function(count, genotype_F, genotype_M){
  private$AF1new[genotype_F,genotype_M] = count
}

Patch$set(which = "public",name = "set_AF1new",
  value = set_AF1new_Patch,overwrite = TRUE
)

#' Get ADMnew
#'
#' Return the new ADM females
#'
get_ADMnew_Patch <- function(){return(private$ADMnew)}

Patch$set(which = "public",name = "get_ADMnew",
  value = get_ADMnew_Patch,overwrite = TRUE
)

#' Set ADMnew
#'
#' Set an element in new ADM males
#'
#' @param count number of mosquitoes
#' @param genotype_M genotype of male
#'
set_ADMnew_Patch <- function(count, genotype_M){
  private$ADMnew[genotype_M] = count
}

Patch$set(which = "public",name = "set_ADMnew",
  value = set_ADMnew_Patch,overwrite = TRUE
)

#' Accumulate ADMnew
#'
#' Accumulate new ADM males
#'
#' @param count vector of new ADM males
#'
accumulate_ADMnew_Patch <- function(count){
  if(length(count)!=length(private$ADMnew)){
    stop("error in accumulate_ADMnew_Patch: lengths of new population and existing population must match")
  }
  private$ADMnew = private$ADMnew + count
}

Patch$set(which = "public",name = "accumulate_ADMnew",
  value = accumulate_ADMnew_Patch,overwrite = TRUE
)

#' Get EGG
#'
#' Return egg stage population
#'
get_EGG_Patch <- function(){return(private$EGG)}

Patch$set(which = "public",name = "get_EGG",
  value = get_EGG_Patch,overwrite = TRUE
)

#' Get LAR
#'
#' Return larval stage population
#'
get_LAR_Patch <- function(){return(private$LAR)}

Patch$set(which = "public",name = "get_LAR",
  value = get_LAR_Patch,overwrite = TRUE
)

#' Get PUP
#'
#' Return pupae stage population
#'
get_PUP_Patch <- function(){return(private$PUP)}

Patch$set(which = "public",name = "get_PUP",
  value = get_PUP_Patch,overwrite = TRUE
)

#' Get ADM
#'
#' Return adult male stage population
#'
get_M_Patch <- function(){return(private$ADM)}

Patch$set(which = "public",name = "get_ADM",
  value = get_M_Patch,overwrite = TRUE
)

#' Get AF1
#'
#' Return adult female stage population
#'
get_F_Patch <- function(){return(private$AF1)}

Patch$set(which = "public",name = "get_AF1",
  value = get_F_Patch,overwrite = TRUE
)

# population delay windows

#' Get EGGdly
#'
#' Return egg stage delay window population
#'
get_EGGdly_Patch <- function(){return(private$EGGdly)}

Patch$set(which = "public",name = "get_EGGdly",
  value = get_EGGdly_Patch,overwrite = TRUE
)

#' Get LARdly
#'
#' Return larval stage delay window population
#'
get_LARdly_Patch <- function(){return(private$LARdly)}

Patch$set(which = "public",name = "get_LARdly",
  value = get_LARdly_Patch,overwrite = TRUE
)

#' Get PUP
#'
#' Return pupae stage delay window population
#'
get_PUPdly_Patch <- function(){return(private$PUPdly)}

Patch$set(which = "public",name = "get_PUPdly",
  value = get_PUPdly_Patch,overwrite = TRUE
)

#' Get ADMdly
#'
#' Return adult male stage delay window population
#'
get_ADMdly_Patch <- function(){return(private$ADMdly)}

Patch$set(which = "public",name = "get_ADMdly",
  value = get_ADMdly_Patch,overwrite = TRUE
)

#' Get AF1dly
#'
#' Return adult female stage delay window population
#'
get_AF1dly_Patch <- function(){return(private$AF1dly)}

Patch$set(which = "public",name = "get_AF1dly",
  value = get_AF1dly_Patch,overwrite = TRUE
)

# migration

#' Get maleMigration
#'
#' Return outbound males (nGenotypes X nPatch integer matrix)
#'
get_maleMigration_Patch <- function(){return(private$maleMigration)}

Patch$set(which = "public",name = "get_maleMigration",
  value = get_maleMigration_Patch,overwrite = TRUE
)

#' Get maleMigration
#'
#' Return outbound males (nGenotypes X nGenotypes X nPatch array)
#'
get_femaleMigration_Patch <- function(){return(private$femaleMigration)}

Patch$set(which = "public",name = "get_femaleMigration",
  value = get_femaleMigration_Patch,overwrite = TRUE
)

# pointers

#' Set Network Pointer
#'
#' Set a reference to the enclosing \code{\link{Network}} object
#'
#' @param NetworkPointer a \code{\link{Network}} object
#'
set_NetworkPointer_Patch <- function(NetworkPointer){private$NetworkPointer = NetworkPointer}

Patch$set(which = "public",name = "set_NetworkPointer",
  value = set_NetworkPointer_Patch,overwrite = TRUE
)

#' Get Network Pointer
#'
#' Return a reference to the enclosing \code{\link{Network}} object
#'
get_NetworkPointer_Patch <- function(){return(private$NetworkPointer)}

Patch$set(which = "public",name = "get_NetworkPointer",
  value = get_NetworkPointer_Patch,overwrite = TRUE
)
