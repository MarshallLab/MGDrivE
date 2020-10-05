###############################################################################
#        ____        __       __
#       / __ \____ _/ /______/ /_
#      / /_/ / __ `/ __/ ___/ __ \
#     / ____/ /_/ / /_/ /__/ / / /
#    /_/    \__,_/\__/\___/_/ /_/
#
#   Patch Class Definition
#   Marshall Lab
#   jared_bennett@berkeley.edu
#   December 2019
#
###############################################################################

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
#'  * timeAq: integer vector of length 3 specifying the length of each aquatic stage
#'  * numPatches: integer, total number of patches in this simulation
#'  * adultEQ: integer, total adult population in this patch for the duration of the simulation
#'  * larvalEQ: integer, total larval population in this patch for the duration of the simulation
#'  * muAq: double vector, length 3, daily death rate for each aquatic stage
#'  * alpha: double, density-dependent centering parameter, see \code{\link{parameterizeMGDrivE}}
#'  * adultRatioF: named double vector, distribution of adult female genotypes, see \code{\link{parameterizeMGDrivE}}
#'  * adultRatioM: named double vector, distribution of adult male genotypes, see \code{\link{parameterizeMGDrivE}}
#'  * larvalRatio: named double vector, distribution of all aquatic genotypes, see \code{\link{parameterizeMGDrivE}}
#'  * eggReleases: egg release schedule for this patch, see \code{\link{basicRepeatedReleases}}
#'  * maleReleases: male release schedule for this patch, see \code{\link{basicRepeatedReleases}}
#'  * femaleReleases: female release schedule for this patch, see \code{\link{basicRepeatedReleases}}
#'  * matedFemaleReleases: mated females release schedule for this patch, see \code{\link{basicRepeatedReleases}}
#'
#' @section **Methods**:
#'  * set_NetworkPointer: see \code{\link{set_NetworkPointer_Patch}}
#'  * get_femalePopulation: see \code{\link{get_femalePop_Patch}}
#'  * get_malePopulation: see \code{\link{get_malePop_Patch}}
#'  * initialPopulation: see \code{\link{set_initialPopulation_Patch}}
#'  * setPopulation: see \code{\link{set_population_deterministic_Patch}} or \code{\link{set_population_stochastic_Patch}}
#'  * reset: see \code{\link{reset_Patch}}
#'  * oneDay_initOutput: see \code{\link{oneDay_initOutput_Patch}}
#'  * oneDay_writeOutput: see \code{\link{oneDay_writeOutput_Patch}}
#'  * oneDay_migrationIn: see \code{\link{oneDay_migrationIn_Patch}}
#'  * oneDay_PopDynamics: see \code{\link{oneDay_PopDynamics_Patch}}
#'  * oneDay_adultD: see \code{\link{oneDay_adultDeath_deterministic_Patch}} or \code{\link{oneDay_adultDeath_stochastic_Patch}}
#'  * oneDay_pupaDM: see \code{\link{oneDay_pupaDM_deterministic_Patch}} or \code{\link{oneDay_pupaDM_stochastic_Patch}}
#'  * oneDay_larvaDM: see \code{\link{oneDay_larvaDM_deterministic_Patch}} or \code{\link{oneDay_larvaDM_stochastic_Patch}}
#'  * oneDay_eggDM: see \code{\link{oneDay_eggDM_deterministic_Patch}} or \code{\link{oneDay_eggDM_stochastic_Patch}}
#'  * oneDay_pupation: see \code{\link{oneDay_pupation_deterministic_Patch}} or \code{\link{oneDay_pupation_stochastic_Patch}}
#'  * oneDay_releases: see \code{\link{oneDay_releases_Patch}}
#'  * oneDay_releaseEggs: see \code{\link{oneDay_eggReleases_Patch}}
#'  * oneDay_mating: see \code{\link{oneDay_mating_deterministic_Patch}} or \code{\link{oneDay_mating_stochastic_Patch}}
#'  * oneDay_layEggs: see \code{\link{oneDay_oviposit_deterministic_Patch}} or \code{\link{oneDay_oviposit_stochastic_Patch}}
#'
#' @section **Fields**:
#'  * patchID: integer ID of this patch
#'  * popAquatic: matrix, nGenotype x sum(timeAquatic), holding all eggs, larva, and pupa
#'  * popMale: vector, nGenotype x 1, holds adult males
#'  * popFemale: matrix, nGenotype x nGenotype, holds mated adult females
#'  * popHolder: vector, nGenotype x 1, temporary population storage
#'  * popPupSex: vector, nGenotype x 1, used in stochastic pupation as another temporary population
#'  * popUnmated: vector, nGenotype x 1, holds unmated females
#'  * popAquatict0: matrix, nGenotype x sum(timeAquatic), holding all eggs, larva, and pupa for reset, see \code{\link{reset_Patch}}
#'  * popMalet0: vector, nGenotype x 1, holds adult males for reset see \code{\link{reset_Patch}}
#'  * popFemalet0: matrix, nGenotype x nGenotype, holds mated adult females for reset see \code{\link{reset_Patch}}
#'  * eggReleases: list of egg releases for this patch. See \code{\link{oneDay_eggReleases_Patch}}
#'  * maleReleases: list of adult male releases for this patch. See \code{\link{oneDay_releases_Patch}}
#'  * femaleReleases: list of adult female releases for this patch. See \code{\link{oneDay_releases_Patch}}
#'  * matedFemaleReleases: list of mated adult female releases for this patch. See \code{\link{oneDay_releases_Patch}}
#'  * NetworkPointer: a reference to enclosing \code{\link{Network}}
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

                initialize = function(patchID, genotypesID, timeAq, numPatches,
                                      adultEQ, larvalEQ, muAq, alpha,
                                      adultRatioF, adultRatioM, larvalRatio,
                                      eggReleases = NULL,
                                      maleReleases = NULL,
                                      femaleReleases = NULL,
                                      matedFemaleReleases = NULL
                                      ){

                  # ID of this patch
                  private$patchID = patchID

                  # initialize objects for simulation. This way, they have dimensions and names
                  nGeno = length(genotypesID)
                  private$popAquatic = matrix(data = 0,  nrow = nGeno, ncol = sum(timeAq),
                                              dimnames = list(genotypesID, NULL))
                  private$popMale = setNames(object = numeric(length = nGeno), nm = genotypesID)
                  private$popFemale = matrix(data = 0, nrow = nGeno, ncol = nGeno,
                                             dimnames = list(genotypesID, genotypesID))
                  private$popHolder = setNames(object = numeric(length = nGeno), nm = genotypesID)
                  private$popPupSex = setNames(object = numeric(length = nGeno), nm = genotypesID)
                  private$popUnmated = setNames(object = numeric(length = nGeno), nm = genotypesID)

                  # set initial population
                  self$setPopulation(adultEQ = adultEQ, larvalEQ = larvalEQ,
                                     adultRatioF = adultRatioF,
                                     adultRatioM = adultRatioM,
                                     larvalRatio = larvalRatio,
                                     timeAq = timeAq, muAq = muAq, alpha = alpha)

                  # store reset populations
                  private$popAquatict0 = private$popAquatic
                  private$popMalet0 = private$popMale
                  private$popFemalet0 = private$popFemale

                  # Mosquito Releases
                  private$eggReleases = eggReleases
                  private$maleReleases = maleReleases
                  private$femaleReleases = femaleReleases
                  private$matedFemaleReleases = matedFemaleReleases

                } # end constructor
              ),

            # private members
            private = list(

              patchID = NULL,

              # temporary populations
              popAquatic = NULL,
              popMale = NULL,
              popFemale = NULL,
              popHolder = NULL,
              popPupSex = NULL, # only used in stochastic pupation function
              popUnmated = NULL,

              # reset populations
              popAquatict0 = NULL,
              popMalet0 = NULL,
              popFemalet0 = NULL,

              # releases
              eggReleases = NULL,
              maleReleases = NULL,
              femaleReleases = NULL,
              matedFemaleReleases = NULL,

              # pointers
              NetworkPointer = NULL

            ) # end private list
)


###############################################################################
# Getters & Setters
###############################################################################

#' Set Network Pointer
#'
#' Set a reference to the enclosing \code{\link{Network}} object
#'
#' @param NetworkPointer A \code{\link{Network}} object
#'
set_NetworkPointer_Patch <- function(NetworkPointer){private$NetworkPointer = NetworkPointer}

Patch$set(which = "public",name = "set_NetworkPointer",
          value = set_NetworkPointer_Patch,overwrite = TRUE
)

#' Get male Population
#'
#' Return males (nGenotypes vector)
#'
get_malePop_Patch <- function(){return(private$popMale)}

Patch$set(which = "public",name = "get_malePopulation",
          value = get_malePop_Patch,overwrite = TRUE
)

#' Get female Population
#'
#' Return  females (nGenotypes X nGenotypes matrix)
#'
get_femalePop_Patch <- function(){return(private$popFemale)}

Patch$set(which = "public",name = "get_femalePopulation",
          value = get_femalePop_Patch,overwrite = TRUE
)

