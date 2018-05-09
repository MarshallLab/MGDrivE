########################################################################
#       _   __     __                      __
#      / | / /__  / /__      ______  _____/ /__
#     /  |/ / _ \/ __/ | /| / / __ \/ ___/ //_/
#    / /|  /  __/ /_ | |/ |/ / /_/ / /  / ,<
#   /_/ |_/\___/\__/ |__/|__/\____/_/  /_/|_|
#
#   Network Migration
#   Marshall Lab
#   November 2017
#
########################################################################

#' Inter-Patch Migration
#'
#' Simulate migration between patches. See \code{\link{MGDrivE-Model}}, 'Migration' section for more details on how inter-patch migration is handled.
#'
oneDay_Migration_Network <- function(){

  ######################################
  # migration out
  ######################################

  for(i in 1:private$nPatch){
    private$patches[[i]]$oneDay_migrationOut()
  }

  ######################################
  # male migration
  ######################################

  # grab moving males
  maleMoveOut = vector(mode="list",length=private$nPatch)
  for(ix in 1:private$nPatch){
    maleMoveOut[[ix]] = private$patches[[ix]]$get_maleMigration()
  }

  # sum over patches
  maleMoveIn = vector(mode="list",length=private$nPatch)
  for(ix in 1:private$nPatch){
    maleMoveIn[[ix]] = rowSums(vapply(X = maleMoveOut,FUN = function(x,ix){x[,ix]},FUN.VALUE = numeric(self$get_genotypesN()),ix=ix))
  }

  ######################################
  # female migration
  ######################################

  # same thing as above for females
  femaleMoveOut = vector(mode="list",length=private$nPatch)
  for(ix in 1:private$nPatch){
    femaleMoveOut[[ix]] = private$patches[[ix]]$get_femaleMigration()
  }

  # sum over patches
  femaleMoveIn = vector(mode="list",length=private$nPatch)
  for(ix in 1:private$nPatch){
    femaleMoveIn[[ix]] = Reduce(f = "+",x = lapply(X = femaleMoveOut,FUN = function(x,ix){
      x[,,ix]
    },ix=ix))
  }

  ######################################
  # migration in
  ######################################

  for(ix in 1:private$nPatch){
    private$patches[[ix]]$oneDay_migrationIn(maleIn = maleMoveIn[[ix]], femaleIn = femaleMoveIn[[ix]])
  }

}

Network$set(which = "public",name = "oneDay_Migration",
            value = oneDay_Migration_Network, overwrite = TRUE
)
