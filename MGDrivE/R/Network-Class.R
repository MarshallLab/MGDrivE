###############################################################################
#       _   __     __                      __
#      / | / /__  / /__      ______  _____/ /__
#     /  |/ / _ \/ __/ | /| / / __ \/ ___/ //_/
#    / /|  /  __/ /_ | |/ |/ / /_/ / /  / ,<
#   /_/ |_/\___/\__/ |__/|__/\____/_/  /_/|_|
#
#   Network Class Definition
#   Marshall Lab
#   jared_bennett@berkeley.edu
#   December 2019
#
###############################################################################
###############################################################################
# Class Definition
###############################################################################

#' Network Class Definition
#'
#' A \code{Network} class object stores all the information for a simulation on
#' a defined landscape.
#'
#' @docType class
#' @format An \code{\link[R6]{R6Class}} generator object
#' @keywords R6 class
#'
#' @section **Constructor**:
#'  * params: see \code{\link{parameterizeMGDrivE}}
#'  * driveCube: an inheritance cube
#'  * patchReleases: see \code{\link{basicRepeatedReleases}} for examples on how to set up release schedules
#'  * migrationMale: a stochastic matrix whose dimensions conform to the number of patches
#'  * migrationFemale: a stochastic matrix whose dimensions conform to the number of patches
#'  * migrationBatch: a list of batch migration parameters. See\code{\link{basicBatchMigration}}
#'  * directory: character string of output directory
#'  * verbose: Chatty? Default is TRUE
#'
#' @section **Methods**:
#'  * get_timeAq: see \code{\link{get_timeAq_Network}}
#'  * get_beta: see \code{\link{get_beta_Network}}
#'  * get_muAd: see \code{\link{get_muAd_Network}}
#'  * get_muAq: see \code{\link{get_muAq_Network}}
#'  * get_alpha: see \code{\link{get_alpha_Network}}
#'  * get_drivecubeindex: see \code{\link{get_drivecubeindex_Network}}
#'  * get_tau: see \code{\link{get_tau_Network}}
#'  * get_genotypesID: see \code{\link{get_genotypesID_Network}}
#'  * get_genotypesN: see \code{\link{get_genotypesN_Network}}
#'  * get_eta: see \code{\link{get_eta_Network}}
#'  * get_phi: see \code{\link{get_phi_Network}}
#'  * get_omega: see \code{\link{get_omega_Network}}
#'  * get_xiF: see \code{\link{get_xiF_Network}}
#'  * get_xiM: see \code{\link{get_xiM_Network}}
#'  * get_s: see \code{\link{get_s_Network}}
#'  * get_nPatch: see \code{\link{get_nPatch_Network}}
#'  * get_conADM: see \code{\link{get_conM_Network}}
#'  * get_conADF: see \code{\link{get_conF_Network}}
#'  * get_tNow: see \code{\link{get_tNow_Network}}
#'  * get_patchReleases: see \code{\link{get_patchReleases_Network}}
#'  * oneDay_Migration: see \code{\link{oneDay_Migration_Deterministic_Network}} or see \code{\link{oneDay_Migration_Stochastic_Network}}
#'  * reset: see \code{\link{reset_Network}}
#'  * oneDay: see \code{\link{oneDay_Network}}
#'  * oneRun: see \code{\link{oneRun_Network}}
#'  * multRun: see \code{\link{multRun_Network}}
#'
#' @section **Fields**:
#'  * parameters: see \code{\link{parameterizeMGDrivE}}
#'  * patches: a list of \code{\link{Patch}} objects
#'  * nPatch: number of patches
#'  * simTime: maximum time of simulation
#'  * sampTime: how often to write output, tNow %% sampTime
#'  * driveCube: an inheritance cube
#'  * tNow: current time of simulation (time starts at 2 because time 1 is the initial equilibrium state)
#'  * runID: an identifier for the current simulation run, useful for Monte Carlo simulation
#'  * directory: a character string of where to store output
#'  * conADM: a \code{\link[base]{connection}} to write male population dynamics out to
#'  * conADF: a \code{\link[base]{connection}} to write female population dynamics out to
#'  * migrationMale: a stochastic matrix whose dimensions conform to the number of patches
#'  * migrationFemale: a stochastic matrix whose dimensions conform to the number of patches
#'  * migrationBatch: list of items for batch migration in stochastic sim.
#'  * mMoveMat: holder object for male migration
#'  * fMoveArray: holder object for female migration
#'  * patchReleases: a list of release schedules for each patch
#'
#' @examples
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

                  # safety check
                  if(length(patchReleases) != params$nPatch){
                    stop("length of patchReleases must equal number of patches in params!")
                  }

                  # set private parameters from input
                  private$parameters = params
                  private$nPatch = params$nPatch
                  private$patches = vector(mode="list",length=private$nPatch)
                  private$simTime = params$simTime
                  private$sampTime = params$sampTime
                  private$driveCube = driveCube
                  private$directory = directory
                  private$runID = params$runID

                  # daily migration objects
                  private$migrationMale = migrationMale
                  private$migrationFemale = migrationFemale
                  private$migrationBatch = migrationBatch

                  # network-level migration holders
                  private$mMoveMat = matrix(data = 0, nrow = private$driveCube$genotypesN,
                                            ncol = private$nPatch,
                                            dimnames = list(private$driveCube$genotypesID, NULL))
                  private$fMoveArray = array(data = 0, dim = c(private$driveCube$genotypesN,
                                                               private$driveCube$genotypesN,
                                                               private$nPatch),
                                             dimnames = list(private$driveCube$genotypesID,
                                                             private$driveCube$genotypesID,
                                                             NULL))

                  private$patchReleases = patchReleases

                  # initialize patches
                  for(i in 1:private$nPatch){

                    # initialize patch i
                    if(verbose){cat("initializing patch: ", i, " of ",
                                    private$nPatch, "\n")}

                   # initialize patch
                    private$patches[[i]] = Patch$new(patchID = i,
                                                     genotypesID = driveCube$genotypesID,
                                                     timeAq = params$timeAq,
                                                     numPatches = private$nPatch,
                                                     adultEQ = params$AdPopEQ[i],
                                                     larvalEQ = params$Leq[i],
                                                     muAq = params$muAq,
                                                     alpha = params$alpha[i],
                                                     adultRatioF = params$AdPopRatio_F[i, ],
                                                     adultRatioM = params$AdPopRatio_M[i, ],
                                                     larvalRatio = params$LarPopRatio[i, ],
                                                     eggReleases = patchReleases[[i]]$eggReleases,
                                                     maleReleases = patchReleases[[i]]$maleReleases,
                                                     femaleReleases = patchReleases[[i]]$femaleReleases,
                                                     matedFemaleReleases = patchReleases[[i]]$matedFemaleReleases
                                                   )

                    # set pointers
                    private$patches[[i]]$set_NetworkPointer(self)
                  }

                  # Output
                  if(!all(dir.exists(directory))){
                    for(f in directory){suppressWarnings(dir.create(f))}
                  } else {
                    dirFiles = list.files(path = directory)
                    if(length(dirFiles)>0){
                      if(verbose){
                        cat("warning: ", length(dirFiles), " files found in the output directory;
                            please move files to avoid being overwritten\n", sep="")
                      }
                    }
                  } # end output if

                } # end constructor

              ),

            # private members
            private = list(

                parameters = NULL,
                patches = NULL, # list of patches
                nPatch = NULL, # number of patches
                simTime = NULL, # max time of sim
                sampTime = NULL, # tNow %% sampTime to write output
                driveCube = NULL, # number of genotypes in simulation
                tNow = 2L, # time starts at 2 because time 1 is the initial condition
                runID = numeric(1),

                # output
                directory = NULL, # directory to store all patch output
                conADM = NULL,
                conADF = NULL,

                # inter-patch migration
                migrationMale = NULL,
                migrationFemale = NULL,
                migrationBatch = NULL,

                # migration
                mMoveMat = NULL,
                fMoveArray = NULL,

                # release schedule
                patchReleases = NULL

              ) # end private members
)

###############################################################################
# Getters & Setters: Parameters
###############################################################################

#' Get timeAq
#'
#' Return duration of aquatic stages.
#'
#' @param stage Character in 'E', 'L', 'P'; if \code{NULL} return total duration
#'
#' @keywords internal
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

#' Get beta
#'
#' Return size of wild-type egg batch
#'
#' @keywords internal
#'
get_beta_Network <- function(){return(private$parameters$beta)}

Network$set(which = "public",name = "get_beta",
  value = get_beta_Network,overwrite = TRUE
)

#' Get muAd
#'
#' Return adult mortality
#'
#' @keywords internal
#'
get_muAd_Network <- function(){return(private$parameters$muAd)}

Network$set(which = "public",name = "get_muAd",
  value = get_muAd_Network,overwrite = TRUE
)

#' Get muAq
#'
#' Return larval mortality, see \code{\link{calcLarvalStageMortalityRate}}
#'
#' @keywords internal
#'
get_muAq_Network <- function(){return(private$parameters$muAq)}

Network$set(which = "public",name = "get_muAq",
  value = get_muAq_Network,overwrite = TRUE
)

#' Get alpha
#'
#' Return density dependent mortality, see \code{\link{calcDensityDependentDeathRate}}
#'
#' @param ix Index of patch
#'
#' @keywords internal
#'
get_alpha_Network <- function(ix){return(private$parameters$alpha[ix])}

Network$set(which = "public",name = "get_alpha",
  value = get_alpha_Network,overwrite = TRUE
)

###############################################################################
# Getters & Setters: Drive Cube
###############################################################################

#' Get Element(s) of Drive Cube by Index
#'
#' Return elements or slices of drive cube. If all \code{NULL} return entire cube.
#'
#' @param fG Female genotype index
#' @param mG Male genotype index
#' @param oG Offspring genotype index
#'
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
#'
get_genotypesID_Network <- function(){return(private$driveCube$genotypesID)}

Network$set(which = "public",name = "get_genotypesID",
  value = get_genotypesID_Network,overwrite = TRUE
)

#' Get genotypesN
#'
#' Return number of possible genotypes
#'
#' @keywords internal
#'
get_genotypesN_Network <- function(){return(private$driveCube$genotypesN)}

Network$set(which = "public",name = "get_genotypesN",
  value = get_genotypesN_Network,overwrite = TRUE
)

#' Get eta
#'
#' @param fIdx Index of female genotype to pull
#'
#' @keywords internal
#'
#' Return genotype-specific mating fitness
#'
get_eta_Network <- function(fIdx){return(private$driveCube$eta[fIdx, ])}

Network$set(which = "public",name = "get_eta",
  value = get_eta_Network,overwrite = TRUE
)

#' Get phi
#'
#' Return genotype-specific sex ratio at emergence
#'
#' @keywords internal
#'
get_phi_Network <- function(){return(private$driveCube$phi)}

Network$set(which = "public",name = "get_phi",
  value = get_phi_Network,overwrite = TRUE
)

#' Get omega
#'
#' Return genotype-specific multiplicative modifier of adult mortality
#'
#' @keywords internal
#'
get_omega_Network <- function(){return(private$driveCube$omega)}

Network$set(which = "public",name = "get_omega",
  value = get_omega_Network,overwrite = TRUE
)

#' Get xiF
#'
#' Return genotype-specific female pupatory success
#'
#' @keywords internal
#'
get_xiF_Network <- function(){return(private$driveCube$xiF)}

Network$set(which = "public",name = "get_xiF",
  value = get_xiF_Network,overwrite = TRUE
)

#' Get xiM
#'
#' Return genotype-specific male pupatory success
#'
#' @keywords internal
#'
get_xiM_Network <- function(){return(private$driveCube$xiM)}

Network$set(which = "public",name = "get_xiM",
  value = get_xiM_Network,overwrite = TRUE
)

#' Get s
#'
#' Return genotype-specific fractional reduction(increase) in fertility
#'
#' @keywords internal
#'
get_s_Network <- function(){return(private$driveCube$s)}

Network$set(which = "public",name = "get_s",
  value = get_s_Network,overwrite = TRUE
)

###############################################################################
# Getters & Setters: Other
###############################################################################

#' Get nPatch
#'
#' Return number of patches
#'
#' @keywords internal
#'
get_nPatch_Network <- function(){return(private$nPatch)}

Network$set(which = "public",name = "get_nPatch",
  value = get_nPatch_Network,overwrite = TRUE
)

#' Get conADM
#'
#' Return \code{\link[base]{connection}} where adult male dynamics are written to
#'
#' @keywords internal
#'
get_conM_Network <- function(){return(private$conADM)}

Network$set(which = "public",name = "get_conADM",
  value = get_conM_Network,overwrite = TRUE
)

#' Get conADF
#'
#' Return \code{\link[base]{connection}} where adult female dynamics are written to
#'
#' @keywords internal
#'
get_conF_Network <- function(){return(private$conADF)}

Network$set(which = "public",name = "get_conADF",
  value = get_conF_Network,overwrite = TRUE
)

#' Get tNow
#'
#' Return current simulation time
#'
#' @keywords internal
#'
get_tNow_Network <- function(){return(private$tNow)}

Network$set(which = "public",name = "get_tNow",
  value = get_tNow_Network,overwrite = TRUE
)

#' Get Patch Release Schedule
#'
#' Return the release schedule for a patch for male or female
#'
#' @param patch Index of patch
#' @param sex Character in 'M', 'F', 'Egg', 'mF'
#'
#' @keywords internal
#'
get_patchReleases_Network <- function(patch, sex = "M"){
  switch(sex,
    M = {return(private$patchReleases[[patch]]$maleReleases)},
    F = {return(private$patchReleases[[patch]]$femaleReleases)},
    Egg = {return(private$patchReleases[[patch]]$eggReleases)}, # use Egg because of R CRAN CHECK partial match issue
    mF = {return(private$patchReleases[[patch]]$matedFemaleReleases)}
  )
}

Network$set(which = "public",name = "get_patchReleases",
  value = get_patchReleases_Network,overwrite = TRUE
)

