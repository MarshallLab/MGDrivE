########################################################################
#       __  _____________       _       ______
#      /  |/  / ____/ __ \_____(_)   __/ ____/
#     / /|_/ / / __/ / / / ___/ / | / / __/
#    / /  / / /_/ / /_/ / /  / /| |/ / /___
#   /_/  /_/\____/_____/_/  /_/ |___/_____/
#
#   MGDrivE Setup
#   Marshall Lab
#   November 2017
#
########################################################################

#' Setup MGDrivE
#'
#' Initialize methods in \code{\link{Patch}} to run deterministic or stochastic simulations. This function must be called prior to any objects being created.
#'
#' @param stochasticityON enable/disable stochastic simulation.
#'
#' @export
MGDrivE.Setup <- function(
  stochasticityON=FALSE
){
  overwrite=TRUE
  cat("initializing MGDrivE\n",sep="")

  stoBools=turnStochasticityOnOrOff(stochasticityON)
  stochasticMove=stoBools[["stochasticMove"]]
  layEggs=stoBools[["layEggs"]]
  larSurvival=stoBools[["larSurvival"]]
  hatchingFract=stoBools[["hatchingFract"]]
  larHatching=stoBools[["larHatching"]]
  eggsFract2=stoBools[["eggsFract2"]]
  larPupating=stoBools[["larPupating"]]
  numMaleFemale=stoBools[["numMaleFemale"]]
  admSurvival=stoBools[["admSurvival"]]
  admPupating=stoBools[["admPupating"]]
  af1Survival=stoBools[["af1Survival"]]
  af1Pupation=stoBools[["af1Pupation"]]
  af1Mating=stoBools[["af1Mating"]]


  if(stochasticMove){
    Patch$set(which = "public",name = "oneDay_migrationOut",
              value = oneDay_migrationOut_stochastic_Patch, overwrite = overwrite
    )
  } else {
    Patch$set(which = "public",name = "oneDay_migrationOut",
              value = oneDay_migrationOut_deterministic_Patch, overwrite = overwrite
    )
  }

  if(layEggs){
    Patch$set(which = "public",name = "oneDay_ovipositG1",
              value = oneDay_ovipositG1_stochastic_Patch, overwrite = overwrite
    )
  } else {
    Patch$set(which = "public",name = "oneDay_ovipositG1",
              value = oneDay_ovipositG1_deterministic_Patch, overwrite = overwrite
    )
  }

  if(larSurvival){
    Patch$set(which = "public",name = "oneDay_larSurvival",
              value = oneDay_larSurvival_stochastic_Patch, overwrite = overwrite
    )
  } else {
    Patch$set(which = "public",name = "oneDay_larSurvival",
              value = oneDay_larSurvival_deterministic_Patch, overwrite = overwrite
    )
  }

  if(hatchingFract){
    Patch$set(which = "public",name = "oneDay_hatchingFract",
              value = oneDay_hatchingFract_stochastic_Patch, overwrite = overwrite
    )
  } else {
    Patch$set(which = "public",name = "oneDay_hatchingFract",
              value = oneDay_hatchingFract_deterministic_Patch, overwrite = overwrite
    )
  }

  if(larHatching){
    Patch$set(which = "public",name = "oneDay_larHatching",
              value = oneDay_larHatching_stochastic_Patch, overwrite = overwrite
    )
  } else {
    Patch$set(which = "public",name = "oneDay_larHatching",
              value = oneDay_larHatching_deterministic_Patch, overwrite = overwrite
    )
  }

  if(eggsFract2){
    Patch$set(which = "public",name = "oneDay_eggsFract2",
              value = oneDay_eggsFract2_stochastic_Patch, overwrite = overwrite
    )
  } else {
    Patch$set(which = "public",name = "oneDay_eggsFract2",
              value = oneDay_eggsFract2_deterministic_Patch, overwrite = overwrite
    )
  }

  if(larPupating){
    Patch$set(which = "public",name = "oneDay_larPupating",
              value = oneDay_larPupating_stochastic_Patch, overwrite = overwrite
    )
  } else {
    Patch$set(which = "public",name = "oneDay_larPupating",
              value = oneDay_larPupating_deterministic_Patch, overwrite = overwrite
    )
  }

  if(numMaleFemale){
    Patch$set(which = "public",name = "oneDay_numMaleFemale",
              value = oneDay_numMaleFemale_stochastic_Patch, overwrite = overwrite
    )
  } else {
    Patch$set(which = "public",name = "oneDay_numMaleFemale",
              value = oneDay_numMaleFemale_deterministic_Patch, overwrite = overwrite
    )
  }

  if(admSurvival){
    Patch$set(which = "public",name = "oneDay_admSurvival",
              value = oneDay_admSurvival_stochastic_Patch, overwrite = overwrite
    )
  } else {
    Patch$set(which = "public",name = "oneDay_admSurvival",
              value = oneDay_admSurvival_deterministic_Patch, overwrite = overwrite
    )
  }

  if(admPupating){
    Patch$set(which = "public",name = "oneDay_admPupating",
              value = oneDay_admPupating_stochastic_Patch, overwrite = overwrite
    )
  } else {
    Patch$set(which = "public",name = "oneDay_admPupating",
              value = oneDay_admPupating_deterministic_Patch, overwrite = overwrite
    )
  }

  if(af1Survival){
    Patch$set(which = "public",name = "oneDay_af1Survival",
              value = oneDay_af1Survival_stochastic_Patch, overwrite = overwrite
    )
  } else {
    Patch$set(which = "public",name = "oneDay_af1Survival",
              value = oneDay_af1Survival_deterministic_Patch, overwrite = overwrite
    )
  }

  if(af1Pupation){
    Patch$set(which = "public",name = "oneDay_af1Pupation",
              value = oneDay_af1Pupation_stochastic_Patch, overwrite = overwrite
    )
  } else {
    Patch$set(which = "public",name = "oneDay_af1Pupation",
              value = oneDay_af1Pupation_deterministic_Patch, overwrite = overwrite
    )
  }

  if(af1Mating){
    Patch$set(which = "public",name = "oneDay_af1Mating",
              value = oneDay_af1Mating_stochastic_Patch, overwrite = overwrite
    )
  } else {
    Patch$set(which = "public",name = "oneDay_af1Mating",
              value = oneDay_af1Mating_deterministic_Patch, overwrite = overwrite
    )
  }

}


#' Enable or Disable Stochastic Model
#'
#' Set switches for deterministic or stochastic model
#'
#' @param on enable/disable stochastic behaviour
#'
turnStochasticityOnOrOff=function(on=TRUE){
  returnList=c(
    larSurvival=FALSE,larHatching=FALSE,larPupating=FALSE,sexDet=FALSE,
    admSurvival=FALSE,admPupating=FALSE,af1Survival=FALSE,af1Pupation=FALSE,
    af1Mating=FALSE, layEggs=FALSE, hatchingFract=FALSE, eggsFract2=FALSE,
    stochasticMove=FALSE,numMaleFemale=FALSE
  )
  if(on){
    returnList=c(
      larSurvival=TRUE,larHatching=TRUE,larPupating=TRUE,sexDet=TRUE,
      admSurvival=TRUE,admPupating=TRUE,af1Survival=TRUE,af1Pupation=TRUE,
      af1Mating=TRUE, layEggs=TRUE, hatchingFract=TRUE, eggsFract2=TRUE,
      stochasticMove=TRUE,numMaleFemale=TRUE
    )
  }
  return(returnList)
}
