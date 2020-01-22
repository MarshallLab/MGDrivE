###############################################################################
#        ____        __       __
#       / __ \____ _/ /______/ /_
#      / /_/ / __ `/ __/ ___/ __ \
#     / ____/ /_/ / /_/ /__/ / / /
#    /_/    \__,_/\__/\___/_/ /_/
#
#   Patch Migration
#   Marshall Lab
#   jared_bennett@berkeley.edu
#   December 2019
#
###############################################################################
###############################################################################
# Outbound Migration
###############################################################################

#' Deterministic Outbound Migration from a Patch
#'
#' Deterministc model of outbound migration from each patch for \code{\link{oneDay_Migration_Network}}.
#' \code{popFemale} from each patch is stored in the \code{fMig} array, and \code{popMale}
#' are stored in the \code{mMig} matrix. Migration location is determined from
#' \code{\link{get_migrationFemaleRow_Network}} and \code{\link{get_migrationMaleRow_Network}}. \cr
#' Batch migration is not used in the deterministic model.
#'
oneDay_migrationOut_deterministic_Patch <- function(){

  # male migration
  # mMig is the matrix of migration to be exposed to the network
  private$mMig[] = private$popMale %*% private$NetworkPointer$get_migrationMaleRow(private$patchID)

  # female migration
  probHolder = private$NetworkPointer$get_migrationFemaleRow(private$patchID)
  # fMig is the matrix of migration to be exposed to the network
  for(i in 1:private$NetworkPointer$get_nPatch()){
    private$fMig[ , ,i] = private$popFemale %x% probHolder[1,i]
  }

}

#' Stochastic Outbound Migration
#'
#' Stochastic model of outbound migration from each patch for \code{\link{oneDay_Migration_Network}}.
#' \code{popFemale} from each patch is stored in the \code{fMig} array, and \code{popMale}
#' are stored in the \code{mMig} matrix. Migration is modeled as a Dirichlet-Multinomial
#' process parameterized by \code{moveVar} multiplied by the migration location
#' probabilities corresponding to this patch (\code{\link{get_migrationFemaleRow_Network}}
#' and \code{\link{get_migrationMaleRow_Network}}). A Dirichlet distributed random
#' variate is sampled from \code{rDirichlet} according to that parameter vector
#' and then movement is sampled from \code{\link[stats]{rmultinom}}. \cr
#' Batch migration begins as a \code{\link[stats]{rbinom}} sampled from \code{\link{get_batchProbs_Network}}.
#' If there is batch migration, the location of migration is sampled uniformly (see \code{\link[base]{sample}}),
#' parameterized by \code{\link{get_batchLocRow_Network}}. The amount of each sex
#' that migrations is sampled from \code{\link[stats]{rbinom}}, parameterized by
#' \code{\link{get_batchSex_Network}}.
#'
oneDay_migrationOut_stochastic_Patch <- function(){

  # things to reuse
  nGeno = private$NetworkPointer$get_genotypesN()

  # Dirchlet sampling of probabilities
  maleProb = rDirichlet(migrationPoint = private$NetworkPointer$get_migrationMaleRow(private$patchID) *
                          private$NetworkPointer$get_moveVar())
  femaleProb = rDirichlet(migrationPoint = private$NetworkPointer$get_migrationFemaleRow(private$patchID) *
                            private$NetworkPointer$get_moveVar())

  # clear migration matrices
  private$mMig[] = 0
  private$fMig[] = 0

  ##########
  # Male/Female
  ##########
  # for each genotype, multinomial over which patch to migrate to
  for(i in 1:nGeno){
    ##########
    # females
    ##########
    for(j in 1:nGeno){
      # if no females, skip
      if(private$popFemale[i,j] > 0){
        # draw multinomial over patches
        private$fMig[i,j, ] = rmultinom(n = 1, size = private$popFemale[i,j], prob = femaleProb)
      }
    } # end second female loop

    ##########
    # males
    ##########
    # if no males, skip
    if(private$popMale[i] > 0){
      # draw multinomial over patches
      private$mMig[i, ] = rmultinom(n = 1, size = private$popMale[i], prob = maleProb)
    }
  } # end migration loops


  ##########
  # Batch
  ##########
  # check for batch migration
  if(rbinom(n = 1, size = 1, prob = private$NetworkPointer$get_batchProbs(private$patchID))){

    # sample over patches, given their relative probabilities
    whichPatch <- sample(x = 1:private$NetworkPointer$get_nPatch(),
                         size = 1, replace = FALSE,
                         prob = private$NetworkPointer$get_batchLocations(private$patchID))

    ##########
    # sample batch out
    ##########
    mBatch <- rbinom(n = nGeno, size = private$mMig[ ,private$patchID],
                     prob = private$NetworkPointer$get_batchSex(private$patchID, 1))
    fBatch <- rbinom(n = nGeno*nGeno, size = private$fMig[ , ,private$patchID],
                     prob = private$NetworkPointer$get_batchSex(private$patchID, 2))

    ##########
    # Combine
    ##########
    private$mMig[ ,private$patchID] <- private$mMig[ ,private$patchID] - mBatch
    private$mMig[ ,whichPatch] <- private$mMig[ ,whichPatch] + mBatch

    private$fMig[ ,private$patchID] <- private$fMig[ ,private$patchID] - fBatch
    private$fMig[ ,whichPatch] <- private$fMig[ ,whichPatch] + fBatch

  } # end batch migration

}


###############################################################################
# Inbound Migration
###############################################################################

#' Inbound Migration
#'
#' Accumulate all inbound migration to this patch.
#'
#' @param maleIn Vector of inbound migration
#' @param femaleIn Matrix of inbound migration
#'
oneDay_migrationIn_Patch <- function(maleIn, femaleIn){
  private$popMale[] = maleIn
  private$popFemale[] = femaleIn
}
