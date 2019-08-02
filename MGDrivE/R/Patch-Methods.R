########################################################################
#        ____        __       __
#       / __ \____ _/ /______/ /_
#      / /_/ / __ `/ __/ ___/ __ \
#     / ____/ /_/ / /_/ /__/ / / /
#    /_/    \__,_/\__/\___/_/ /_/#
#
#   Patch Class Implementation
#   Marshall Lab
#   November 2017
#
########################################################################

#' Reset Patch to Initial Conditions
#'
#' Resets a patch to its initial configuration so that a new one does not have
#' to be created and allocated in the network (for Monte Carlo simulation).
#'
#' @param verbose Chatty? Default is TRUE
#'
reset_Patch <- function(verbose = TRUE){

  if(verbose){cat("reset patch ",private$patchID,"\n",sep="")}

  # reset initial population size
  private$EGG[[1]] = private$EGGt0
  private$LAR[[1]] = private$LARt0
  private$PUP[[1]] = private$PUPt0
  private$ADM[[1]] = private$ADMt0
  private$AF1[[1]] = private$AF1t0

  # Population Array
  private$EGG[2:private$NetworkPointer$get_simTime()] =  initPopVectorArray(private$NetworkPointer$get_genotypesID(),private$NetworkPointer$get_simTime()-1)
  private$LAR[2:private$NetworkPointer$get_simTime()] =  initPopVectorArray(private$NetworkPointer$get_genotypesID(),private$NetworkPointer$get_simTime()-1)
  private$PUP[2:private$NetworkPointer$get_simTime()] =  initPopVectorArray(private$NetworkPointer$get_genotypesID(),private$NetworkPointer$get_simTime()-1)
  private$ADM[2:private$NetworkPointer$get_simTime()] =  initPopVectorArray(private$NetworkPointer$get_genotypesID(),private$NetworkPointer$get_simTime()-1)
  private$AF1[2:private$NetworkPointer$get_simTime()] =  initPopMatrixArray(private$NetworkPointer$get_genotypesID(),private$NetworkPointer$get_simTime()-1)

  # Windowed Population Tensors
  private$EGGdly =  primePopVectorArray(private$EGG[[1]],private$NetworkPointer$get_windowSize())
  private$LARdly =  primePopVectorArray(private$LAR[[1]],private$NetworkPointer$get_windowSize())
  private$PUPdly =  primePopVectorArray(private$PUP[[1]],private$NetworkPointer$get_windowSize())
  private$ADMdly =  primePopVectorArray(private$ADM[[1]],private$NetworkPointer$get_windowSize())
  private$AF1dly =  primePopMatrixArray(private$AF1[[1]],private$NetworkPointer$get_windowSize())

  # Reset Mosquito Releases
  private$maleReleases = private$NetworkPointer$get_patchReleases(private$patchID,"M")
  private$femaleReleases = private$NetworkPointer$get_patchReleases(private$patchID,"F")

}

Patch$set(which = "public",name = "reset",
          value = reset_Patch, overwrite = TRUE
)


#' Initialize Output from Focal Patch
#'
#' Writes output to the text connections specified in the enclosing \code{\link{Network}}
#'
oneDay_initOutput_Patch <- function(){

  # write males
  ADMout = paste0(c(1,private$patchID,private$ADM[[1]]),collapse = ",")
  writeLines(text = ADMout,con = private$NetworkPointer$get_conADM(),sep = "\n")

  # write females
  AF1out = paste0(c(1,private$patchID,c(t(private$AF1[[1]]))),collapse = ",")
  writeLines(text = AF1out,con = private$NetworkPointer$get_conAF1(),sep = "\n")

}

Patch$set(which = "public",name = "oneDay_initOutput",
          value = oneDay_initOutput_Patch, overwrite = TRUE
)


#' Write Output from Focal Patch
#'
#' Writes output to the text connections specified in the enclosing \code{\link{Network}}
#'
oneDay_writeOutput_Patch <- function(){

  tNow = private$NetworkPointer$get_tNow()

  # write males
  ADMout = paste0(c(tNow,private$patchID,private$ADM[[tNow]]),collapse = ",")
  writeLines(text = ADMout,con = private$NetworkPointer$get_conADM(),sep = "\n")

  # write females
  AF1out = paste0(c(tNow,private$patchID,c(t(private$AF1[[tNow]]))),collapse = ",")
  writeLines(text = AF1out,con = private$NetworkPointer$get_conAF1(),sep = "\n")

}

Patch$set(which = "public",name = "oneDay_writeOutput",
          value = oneDay_writeOutput_Patch, overwrite = TRUE
)
