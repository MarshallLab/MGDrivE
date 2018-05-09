########################################################################
#        ____        __       __
#       / __ \____ _/ /______/ /_
#      / /_/ / __ `/ __/ ___/ __ \
#     / ____/ /_/ / /_/ /__/ / / /
#    /_/    \__,_/\__/\___/_/ /_/
#
#   Patch Migration
#   Marshall Lab
#   November 2017
#
########################################################################

########################################################################
# Outbound Migration
########################################################################

#' Deterministic Oubound Migration from a Patch
#'
#' Deterministc model of outbound migration of \code{AF1new} females from this patch, fills up the \code{femaleMigration} array.
#'
oneDay_migrationOut_deterministic_Patch <- function(){
  genotypesID = private$NetworkPointer$get_genotypesID()
  genotypesN = private$NetworkPointer$get_genotypesN()

  # male migration
  maleMatrix = private$NetworkPointer$get_migrationMaleRow(private$patchID)
  # maleMigration is the matrix of migration to be exposed to the network
  private$maleMigration = matrix(data=private$ADMnew,nrow=genotypesN,dimnames=list(genotypesID,NULL)) %*% maleMatrix


  # female migration
  femaleMatrix = private$NetworkPointer$get_migrationFemaleRow(private$patchID)
  # femaleMigration is the matrix of migration to be exposed to the network
  private$femaleMigration = array(data = rep(0,private$NetworkPointer$get_nPatch()),dim=c(genotypesN,genotypesN,private$NetworkPointer$get_nPatch()),
                                      dimnames = list(genotypesID,genotypesID,NULL))
  for(ix in 1:private$NetworkPointer$get_nPatch()){
    private$femaleMigration[,,ix] = private$AF1new %x% femaleMatrix[1,ix]
  }

}

#' Stochastic Oubound Migration
#'
#' Stochastic model of migration of \code{AF1new} females from this patch, fills up the \code{femaleMigration} array.
#' Migration is modeled as a Dirichlet-Multinomial process parameterized by \code{moveVar} multiplied by the row corresponding to this
#' patch from the stochastic matrix. A Dirichlet distributed random variate is sampled from \code{\link[MCMCpack]{rdirichlet}} according to that
#' parameter vector and then movement is sampled from \code{\link[stats]{rmultinom}}.
#'
oneDay_migrationOut_stochastic_Patch <- function(){

  genotypesID = private$NetworkPointer$get_genotypesID()
  genotypesN = private$NetworkPointer$get_genotypesN()

  # Males
  # Dirchlet sampling of probabilities
  maleMatrix = private$NetworkPointer$get_migrationMaleRow(private$patchID)

  maleProb = rDirichlet(migrationPoint = maleMatrix*private$NetworkPointer$get_moveVar())
  private$maleMigration = matrix(data=0,nrow=genotypesN,ncol=private$NetworkPointer$get_nPatch())
  for(i in 1:nrow(private$maleMigration)){
    if(private$ADMnew[i]<.Machine$double.eps){
      next()
    } else {
      private$maleMigration[i,] = t(rmultinom(n = 1,size = private$ADMnew[i],prob = maleProb))
    }
  }

  # Females
  # Dirchlet sampling of probabilities
  femaleMatrix = private$NetworkPointer$get_migrationFemaleRow(private$patchID)

  private$femaleMigration = array(data = rep(0,private$NetworkPointer$get_nPatch()),dim=c(genotypesN,genotypesN,private$NetworkPointer$get_nPatch()),
                                      dimnames = list(genotypesID,genotypesID,NULL))
  for(i in 1:genotypesN){
    for(j in 1:genotypesN){
      femaleProb = rDirichlet(migrationPoint = as.vector(femaleMatrix)*private$NetworkPointer$get_moveVar())
      private$femaleMigration[i,j,] = t(rmultinom(n = 1,size = private$AF1new[i,j],prob = femaleMatrix))
    }
  }

}


########################################################################
# Inbound Migration
########################################################################

#' Inbound Migration
#'
#' Accumulate all inbound migration to this patch.
#'
#' @param maleIn vector of inbound migration
#' @param femaleIn matrix of inbound migration
#'
oneDay_migrationIn_Patch <- function(maleIn, femaleIn){
  private$ADMnew = maleIn
  private$AF1new = femaleIn
}

Patch$set(which = "public",name = "oneDay_migrationIn",
          value = oneDay_migrationIn_Patch, overwrite = TRUE
)
