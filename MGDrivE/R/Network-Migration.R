###############################################################################
#       _   __     __                      __
#      / | / /__  / /__      ______  _____/ /__
#     /  |/ / _ \/ __/ | /| / / __ \/ ___/ //_/
#    / /|  /  __/ /_ | |/ |/ / /_/ / /  / ,<
#   /_/ |_/\___/\__/ |__/|__/\____/_/  /_/|_|
#
#   Network Migration
#   Marshall Lab
#   jared_bennett@berkeley.edu
#   June 2020
#
###############################################################################
###############################################################################
# Outbound Migration
###############################################################################

#' Deterministic Inter-Patch Migration
#'
#' Deterministic model of interpatch migration from each patch.
#' \code{popFemale}/\code{popMale} is retrieved from each patch using
#' \code{\link{get_femalePop_Patch}}/\code{\link{get_malePop_Patch}}. Migration
#' location is determined from the supplied matrices, \code{private$migrationFemale}
#' or \code{private$migrationMale}. Final migration status is held in
#' \code{private$fMoveArray} or \code{private$mMoveMat}. \cr
#' Batch migration is not used in the deterministic model.
#'
#' This function handles outbound and inbound migration. See \code{\link{MGDrivE-Model}},
#' 'Migration' section for more details on how inter-patch migration is handled.
#'
oneDay_Migration_Deterministic_Network <- function(){

  ######################################
  # Clear Holder Objects
  ######################################
  private$mMoveMat[] <- 0
  private$fMoveArray[] <- 0


  ######################################
  # Migration Out
  ######################################
  # loop over patches, migrate things
  for(i in 1:private$nPatch){
    # grab moving males
    # outputs matrix
    #  outer() would return an array.
    private$mMoveMat[] <- private$mMoveMat + private$patches[[i]]$get_malePopulation() %*%
                                              private$migrationMale[i, ,drop = FALSE]

    # grab moving females
    # outputs 3d array
    probHolder <- private$migrationFemale[i, ,drop = FALSE]
    for(i in 1:private$nPatch){
      private$fMoveArray[ , ,i] <- private$fMoveArray[ , ,i] + private$patches[[i]]$get_femalePopulation() %x%
                                                          probHolder[1,i]
    } # end female loop

  } # end patch loop

  ######################################
  # Migration In
  ######################################
  for(i in 1:private$nPatch){
    private$patches[[i]]$oneDay_migrationIn(maleIn = private$mMoveMat[ ,i], femaleIn = private$fMoveArray[ , ,i])
  }

}


#' Stochastic Inter-Patch Migration
#'
#' Stochastic model of interpatch migration from each patch.
#' \code{popFemale}/\code{popMale} is retrieved from each patch using
#' \code{\link{get_femalePop_Patch}}/\code{\link{get_malePop_Patch}}. Migration
#' location is determined from the supplied matrices, \code{private$migrationFemale}
#' or \code{private$migrationMale}. Migration is modeled as a Multinomial
#' process parameterized by migration location probabilities corresponding to each
#' patch . Movement is sampled from \code{\link[stats]{rmultinom}}. \cr
#' Batch migration begins as a \code{\link[stats]{rbinom}} sampled from
#' \code{private$migrationBatch$batchProbs}.If there is batch migration, the
#' location of migration is sampled uniformly (see \code{\link[base]{sample}}),
#' parameterized by \code{private$migrationBatch$moveProbs}. The amount of each sex
#' that migrations is sampled from \code{\link[stats]{rbinom}}, parameterized by
#' \code{private$migrationBatch$sexProbs}.
#'
#' This function handles outbound and inbound migration. See \code{\link{MGDrivE-Model}},
#' 'Migration' section for more details on how inter-patch migration is handled.
#'
oneDay_Migration_Stochastic_Network <- function(){

  ######################################
  # Clear holder objects
  ######################################
  private$mMoveMat[] <- 0
  private$fMoveArray[] <- 0

  # things to reuse
  nGeno <- private$driveCube$genotypesN


  ######################################
  # Migration Out
  ######################################
  # loop over patches, migrate things
  for(i in 1:private$nPatch){

    ##########
    # Batch Check
    ##########
    # check for batch migration
    #  if this patch has it, store current state of population for use later
    bCheck <- rbinom(n = 1, size = 1, prob = private$migrationBatch$batchProbs[i])

    if(bCheck){
      # vector, current state of males in patch
      mState <- private$mMoveMat[ ,i]

      # matrix, current state of females in patch
      fState <- private$fMoveArray[ , ,i]
    }

    ##########
    # Migration
    ##########
    # things to reuse
    maleProb <- private$migrationMale[i, ,drop = FALSE]
    femaleProb <- private$migrationFemale[i, ,drop = FALSE]

    # for each genotype, multinomial over which patch to migrate to
    for(k in 1:nGeno){
      ##########
      # females
      ##########
      for(j in 1:nGeno){
        # if no females, skip
        if(private$patches[[i]]$get_femalePopulation()[k,j] > 0){
          # draw multinomial over patches
          private$fMoveArray[k,j, ] = private$fMoveArray[k,j, ] + rmultinom(n = 1,
                                                                            size = private$patches[[i]]$get_femalePopulation()[k,j],
                                                                            prob = femaleProb)
        }
      } # end second female loop

      ##########
      # males
      ##########
      # if no males, skip
      if(private$patches[[i]]$get_malePopulation()[k] > 0){
        # draw multinomial over patches
        private$mMoveMat[k, ] = private$mMoveMat[k, ] + rmultinom(n = 1,
                                                                  size = private$patches[[i]]$get_malePopulation()[k],
                                                                  prob = maleProb)
      }
    } # end migration loops


    ##########
    # Batch Migration
    ##########
    if(bCheck){

      # sample over patches, given their relative probabilities
      whichPatch <- sample(x = 1:private$nPatch,
                           size = 1, replace = FALSE,
                           prob = private$migrationBatch$moveProbs[i, ,drop = FALSE])

      ##########
      # sample batch out
      ##########
      # the population available to move is the current population remaining in the
      #  patch after regular migration has taken place.
      # Because the migration object has everything for the network, we save the
      #  state of the network, then do the current patch migration.
      #  The difference between the updated and the saved state is the population
      #  available for batch migratio, ie, the population of the patch that didn't move.
      mBatch <- rbinom(n = nGeno, size = private$mMoveMat[ ,i] - mState,
                       prob = private$migrationBatch$sexProbs[i,1])

      fBatch <- rbinom(n = nGeno*nGeno, size = private$fMoveArray[ , ,i] - fState,
                       prob = private$migrationBatch$sexProbs[i,2])

      ##########
      # Combine
      ##########
      private$mMoveMat[ ,i] <- private$mMoveMat[ ,i] - mBatch
      private$mMoveMat[ ,whichPatch] <- private$mMoveMat[ ,whichPatch] + mBatch

      private$fMoveArray[ , ,i] <- private$fMoveArray[ , ,i] - fBatch
      private$fMoveArray[ , ,whichPatch] <- private$fMoveArray[ , ,whichPatch] + fBatch

    } # end batch migration

  } # end patch loop


  ######################################
  # Migration In
  ######################################
  for(i in 1:private$nPatch){
    private$patches[[i]]$oneDay_migrationIn(maleIn = private$mMoveMat[ ,i], femaleIn = private$fMoveArray[ , ,i])
  }

}

