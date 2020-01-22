###############################################################################
#       __  _____________       _       ______
#      /  |/  / ____/ __ \_____(_)   __/ ____/
#     / /|_/ / / __/ / / / ___/ / | / / __/
#    / /  / / /_/ / /_/ / /  / /| |/ / /___
#   /_/  /_/\____/_____/_/  /_/ |___/_____/
#
#   setupMGDrivE
#   Marshall Lab
#   jared_bennett@berkeley.edu
#   December 2019
#
###############################################################################

#' Setup MGDrivE
#'
#' Initialize methods in \code{\link{Patch}} to run deterministic or stochastic simulations.
#' This sets internal function definitions so that \code{\link{oneRun_Network}}
#' and \code{\link{multRun_Network}} run either deterministic or stochastic functions.
#'
#' @param stochasticityON Enable/disable stochastic simulation. Default is FALSE, implying deterministic simulation
#' @param verbose Chatty? Default is TRUE
#'
#' @examples
#' # run deterministic MGDrivE
#' setupMGDrivE(stochasticityON = FALSE)
#'
#' # run stochastic MGDrivE
#' setupMGDrivE(stochasticityON = TRUE)
#'
#' @export
setupMGDrivE <- function(stochasticityON = FALSE, verbose = TRUE){

  overwrite=TRUE
  if(verbose){cat("initializing MGDrivE\n",sep="")}

  ##########
  # Things that change
  ##########
  if(stochasticityON){
    # stochastic option
    Patch$set(which = "public", name = "setPopulation",
              value = set_population_stochastic_Patch, overwrite = overwrite)

    Patch$set(which = "public",name = "oneDay_adultD",
              value = oneDay_adultDeath_stochastic_Patch, overwrite = overwrite)

    Patch$set(which = "public",name = "oneDay_pupaDM",
              value = oneDay_pupaDM_stochastic_Patch, overwrite = overwrite)

    Patch$set(which = "public",name = "oneDay_larvaDM",
              value = oneDay_larvaDM_stochastic_Patch, overwrite = overwrite)

    Patch$set(which = "public",name = "oneDay_eggDM",
              value = oneDay_eggDM_stochastic_Patch, overwrite = overwrite)

    Patch$set(which = "public",name = "oneDay_pupation",
              value = oneDay_pupation_stochastic_Patch, overwrite = overwrite)

    Patch$set(which = "public",name = "oneDay_mating",
              value = oneDay_mating_stochastic_Patch, overwrite = overwrite)

    Patch$set(which = "public",name = "oneDay_layEggs",
              value = oneDay_oviposit_stochastic_Patch, overwrite = overwrite)

    Patch$set(which = "public",name = "oneDay_migrationOut",
          value = oneDay_migrationOut_stochastic_Patch, overwrite = overwrite)

  } else {
    # deterministic option
    Patch$set(which = "public", name = "setPopulation",
              value = set_population_deterministic_Patch, overwrite = overwrite)

    Patch$set(which = "public",name = "oneDay_adultD",
              value = oneDay_adultDeath_deterministic_Patch, overwrite = overwrite)

    Patch$set(which = "public",name = "oneDay_pupaDM",
              value = oneDay_pupaDM_deterministic_Patch, overwrite = overwrite)

    Patch$set(which = "public",name = "oneDay_larvaDM",
              value = oneDay_larvaDM_deterministic_Patch, overwrite = overwrite)

    Patch$set(which = "public",name = "oneDay_eggDM",
              value = oneDay_eggDM_deterministic_Patch, overwrite = overwrite)

    Patch$set(which = "public",name = "oneDay_pupation",
              value = oneDay_pupation_deterministic_Patch, overwrite = overwrite)

    Patch$set(which = "public",name = "oneDay_mating",
              value = oneDay_mating_deterministic_Patch, overwrite = overwrite)

    Patch$set(which = "public",name = "oneDay_layEggs",
              value = oneDay_oviposit_deterministic_Patch, overwrite = overwrite)

    Patch$set(which = "public",name = "oneDay_migrationOut",
              value = oneDay_migrationOut_deterministic_Patch, overwrite = overwrite)

  } # end if

  ##########
  # Things that stay the same
  ##########
  Patch$set(which = "public", name = "initialPopulation",
          value = set_initialPopulation_Patch, overwrite = overwrite)

  Patch$set(which = "public", name = "reset",
          value = reset_Patch, overwrite = overwrite)

  Patch$set(which = "public",name = "oneDay_initOutput",
          value = oneDay_initOutput_Patch, overwrite = overwrite)

  Patch$set(which = "public",name = "oneDay_writeOutput",
          value = oneDay_writeOutput_Patch, overwrite = overwrite)

  Patch$set(which = "public",name = "oneDay_releases",
          value = oneDay_releases_Patch, overwrite = overwrite)

  Patch$set(which = "public",name = "oneDay_releaseEggs",
          value = oneDay_eggReleases_Patch, overwrite = overwrite)

  Patch$set(which = "public",name = "oneDay_migrationIn",
          value = oneDay_migrationIn_Patch, overwrite = overwrite)

  Patch$set(which = "public",name = "oneDay_PopDynamics",
          value = oneDay_PopDynamics_Patch, overwrite = overwrite)

}
