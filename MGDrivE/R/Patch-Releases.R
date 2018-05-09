########################################################################
#        ____        __       __
#       / __ \____ _/ /______/ /_
#      / /_/ / __ `/ __/ ___/ __ \
#     / ____/ /_/ / /_/ /__/ / / /
#    /_/    \__,_/\__/\___/_/ /_/#
#
#   Patch Mosquito Releases
#   Marshall Lab
#   November 2017
#
########################################################################

#' Release Male Mosquitoes in a Patch
#'
#' Based on this patch's release schedule, handle daily releases.
#'
oneDay_maleReleases_Patch <- function(){
  # male releases
  if(length(private$maleReleases) > 0){
    if(private$maleReleases[[1]]$tRelease <= private$NetworkPointer$get_tNow()){
      private$ADMnew = private$admSurvival + private$admPupating + private$maleReleases[[1]]$nuM
      private$maleReleases[[1]] = NULL
    } else {
      private$ADMnew = private$admSurvival + private$admPupating
    }

  } else {
    private$ADMnew = private$admSurvival + private$admPupating
  }
}

Patch$set(which = "public",name = "oneDay_maleReleases",
          value = oneDay_maleReleases_Patch, overwrite = TRUE
)

#' Release Female Mosquitoes in a Patch
#'
#' Based on this patch's release schedule, handle daily releases.
#'
oneDay_femaleReleases_Patch <- function(){
  # female releases
  if(length(private$femaleReleases) > 0){
    if(private$femaleReleases[[1]]$tRelease <= private$NetworkPointer$get_tNow()){
      private$af1Pupation = private$af1Pupation + private$femaleReleases[[1]]$nuF
      private$femaleReleases[[1]] = NULL
    }

  }

}

Patch$set(which = "public",name = "oneDay_femaleReleases",
          value = oneDay_femaleReleases_Patch, overwrite = TRUE
)
