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

#' Deterministic Outbound Migration from a Patch
#'
#' Deterministc model of outbound migration of \code{AF1new} females from this patch, fills up the \code{femaleMigration} array.
#'
oneDay_migrationOut_deterministic_Patch <- function(){

  # male migration
  maleMatrix = private$NetworkPointer$get_migrationMaleRow(private$patchID)
  # maleMigration is the matrix of migration to be exposed to the network
  private$maleMigration[] = private$ADMnew %*% maleMatrix


  # female migration
  femaleMatrix = private$NetworkPointer$get_migrationFemaleRow(private$patchID)
  # femaleMigration is the matrix of migration to be exposed to the network
  for(ix in 1:private$NetworkPointer$get_nPatch()){
    private$femaleMigration[,,ix] = private$AF1new %x% femaleMatrix[1,ix]
  }

}

#' Stochastic Outbound Migration
#'
#' Stochastic model of migration of \code{AF1new} females from this patch, fills up the \code{femaleMigration} array.
#' Migration is modeled as a Dirichlet-Multinomial process parameterized by \code{moveVar} multiplied by the row corresponding to this
#' patch from the stochastic matrix. A Dirichlet distributed random variate is sampled from \code{rDirichlet} according to that
#' parameter vector and then movement is sampled from \code{\link[stats]{rmultinom}}.
#'
oneDay_migrationOut_stochastic_Patch <- function(){

  genotypesN = private$NetworkPointer$get_genotypesN()

  # Dirchlet sampling of probabilities
  maleProb = rDirichlet(migrationPoint = private$NetworkPointer$get_migrationMaleRow(private$patchID) *
                          private$NetworkPointer$get_moveVar())
  femaleProb = rDirichlet(migrationPoint = private$NetworkPointer$get_migrationFemaleRow(private$patchID) *
                            private$NetworkPointer$get_moveVar())

  # for each genotype, multinomial over which patch to migrate to
  for(i in 1:genotypesN){
    for(j in 1:genotypesN){
      # pull females
      if(private$AF1new[i,j] < .Machine$double.eps){
        private$femaleMigration[i,j,] = 0
      } else {
        private$femaleMigration[i,j,] = rmultinom(n = 1,size = private$AF1new[i,j],prob = femaleProb)
      }
    }
    # pull males
    if(private$ADMnew[i] < .Machine$double.eps){
      private$maleMigration[i,] = 0
    } else {
      private$maleMigration[i,] = rmultinom(n = 1,size = private$ADMnew[i],prob = maleProb)
    }
  }# end migration loops


  # Batch migration
  if(rbinom(n = 1, size = 1, prob = private$NetworkPointer$get_batchProbs(private$patchID))){

    # probs vecs
    mProb <- private$NetworkPointer$get_batchSex(private$patchID, 1)

    fProb <- private$NetworkPointer$get_batchSex(private$patchID, 2)


    # hold structures
    maleHold <- matrix(data = 0L, nrow = genotypesN, ncol = 2L)
    femaleHold <- array(data = 0L, dim = c(genotypesN,genotypesN,2L))

    # pull over males/females
    for(row in 1:genotypesN){
      for(col in 1:genotypesN){
        # pull over females
        if(private$femaleMigration[row,col,private$patchID] < .Machine$double.eps){
          femaleHold[row,col, ] = 0
        } else {
          femaleHold[row,col, ] <- rmultinom(n = 1L,
                                             size = private$femaleMigration[row,col,private$patchID],
                                             prob = fProb)
        }
      }
      # pull over males
      if(private$maleMigration[row,private$patchID] < .Machine$double.eps){
        maleHold[row, ] = 0
      } else {
        maleHold[row, ] <- rmultinom(n = 1L,
                                     size = private$maleMigration[row,private$patchID],
                                     prob = mProb)
      }
    }# end sampling loops

    # sample over patches, given their relative probabilities
    whichPatch <- sample(x = 1:private$NetworkPointer$get_nPatch(),
                         size = 1, replace = FALSE,
                         prob = private$NetworkPointer$get_batchLocations(private$patchID))

    # combine with regular migration
    private$maleMigration[ ,private$patchID] <- maleHold[ ,1L]
    private$maleMigration[, whichPatch] <- private$maleMigration[, whichPatch] + maleHold[ ,2L]

    private$femaleMigration[ , ,private$patchID] <- femaleHold[ , ,1L]
    private$femaleMigration[ , ,whichPatch] <- private$femaleMigration[ , ,whichPatch] + femaleHold[ , ,2L]

  }# end batch migration

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
