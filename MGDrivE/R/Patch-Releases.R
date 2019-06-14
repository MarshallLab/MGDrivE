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
    # browser()
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
    # browser()
    if(private$femaleReleases[[1]]$tRelease <= private$NetworkPointer$get_tNow()){
      private$af1Pupation = private$af1Pupation + private$femaleReleases[[1]]$nuF
      private$femaleReleases[[1]] = NULL
    }
  }

}

Patch$set(which = "public",name = "oneDay_femaleReleases",
          value = oneDay_femaleReleases_Patch, overwrite = TRUE
)

#' Release Eggs into hatching eggs
#'
#' This function performs egg releases into hatching eggs. This allows them
#' to under-go density-dependence
#'
oneDay_eggReleases2eggs_Patch <- function(){

#  if(private$NetworkPointer$get_tNow() > 19) browser()

  if(length(private$eggReleases)>0){
    if(private$eggReleases[[1]]$tRelease <= private$NetworkPointer$get_tNow() && !private$eggReleases[[1]]$releaseFlag){
      private$hatchingFract = private$hatchingFract + private$eggReleases[[1]]$nuE
      private$eggReleases[[1]]$releaseFlag = TRUE # set this to do larval windowed release later
    }
  }

}

Patch$set(which = "public",name = "oneDay_eggReleases2eggs",
          value = oneDay_eggReleases2eggs_Patch, overwrite = TRUE
)

#' Release Eggs into Larval Window
#'
#' Get eggs that were released into larval window t - aquaStageDuration days
#' ago and put them into the oviposit vector
#'
oneDay_eggReleases2adults_Patch <- function(){

  # check that there are releases to be done
  if(length(private$eggReleases)>0){
    # get constants we need

    tNow = private$NetworkPointer$get_tNow()
    period = private$NetworkPointer$get_timeAq("L") + private$NetworkPointer$get_timeAq("P")

    if((private$eggReleases[[1]]$tRelease + period <= tNow) && private$eggReleases[[1]]$releaseFlag ){
      private$ovipositG1 = private$ovipositG1 + private$eggReleases[[1]]$nuE
      private$eggReleases[[1]] = NULL
    }
  }

}

Patch$set(which = "public",name = "oneDay_eggReleases2adults",
          value = oneDay_eggReleases2adults_Patch, overwrite = TRUE
)
