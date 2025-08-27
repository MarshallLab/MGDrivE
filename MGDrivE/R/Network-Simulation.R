###############################################################################
#       _   __     __                      __
#      / | / /__  / /__      ______  _____/ /__
#     /  |/ / _ \/ __/ | /| / / __ \/ ___/ //_/
#    / /|  /  __/ /_ | |/ |/ / /_/ / /  / ,<
#   /_/ |_/\___/\__/ |__/|__/\____/_/  /_/|_|
#
#   Network Class Mosquito Population Simulation
#   Marshall Lab
#   jared_bennett@berkeley.edu
#   December 2019
#
###############################################################################
###############################################################################
# Main Run Functions
###############################################################################

#' Run Simulation
#'
#' Run a single simulation on this network.
#'
#' @param verbose Chatty? Default is TRUE
#'
oneRun_Network <- function(verbose = TRUE){

  ####################
  # open connections
  ####################
  # can't use on.exit() to close files because it closes onces this function goes out
  #  of scope, which happens before the simulation even runs!
  private$conADM = file(description = paste0(private$directory[1],
                                             .Platform$file.sep,
                                             "M_Run",
                                             formatC(x = private$runID, width = 3, format = "d", flag = "0"),
                                             ".csv"),
                        open = "wt")
  private$conADF = file(description = paste0(private$directory[1],
                                             .Platform$file.sep,
                                             "F_Run",
                                             formatC(x = private$runID, width = 3, format = "d", flag = "0"),
                                             ".csv"),
                        open = "wt")

  ####################
  # begin run
  ####################
  if(verbose){cat("begin run ",private$runID,"\n",sep="")}

  ####################
  # setup output
  ####################
  for(i in 1:private$nPatch){
    private$patches[[i]]$oneDay_initOutput()
  }

  if(verbose){pb = txtProgressBar(min = 0,max = private$simTime,style = 3)}

  # run simulation
  while(private$simTime >= private$tNow){

    self$oneDay()

    private$tNow = private$tNow + 1L
    if(verbose){setTxtProgressBar(pb,value = private$tNow)}
  }

  ####################
  # close connections
  ####################
  close(private$conADM)
  close(private$conADF)

  if(verbose){cat("run ",private$runID," over\n",sep="")}
}

Network$set(which = "public",name = "oneRun",
          value = oneRun_Network, overwrite = TRUE
)

#' Run Simulation
#'
#' Run multiple simulations on this network
#'
#' @param verbose Chatty? Default is TRUE
#'
multRun_Network <- function(verbose = TRUE){

  ####################
  # setup loop over number of runs
  ####################
  for(run in 1:length(private$directory)){
    ####################
    # open connections
    ####################
    private$conADM = file(description = paste0(private$directory[run],
                                               .Platform$file.sep,
                                               "M_Run",
                                               formatC(x = private$runID, width = 3, format = "d", flag = "0"),
                                               ".csv"),
                          open = "wt")
    private$conADF = file(description = paste0(private$directory[run],
                                               .Platform$file.sep,
                                               "F_Run",
                                               formatC(x = private$runID, width = 3, format = "d", flag = "0"),
                                               ".csv"),
                          open = "wt")

    ####################
    # begin runs
    ####################
    if(verbose){cat("begin run ",private$runID,"\n",sep="")}

    ####################
    # setup output
    ####################
    for(i in 1:private$nPatch){
      private$patches[[i]]$oneDay_initOutput()
    }

    if(verbose){pb = txtProgressBar(min = 0,max = private$simTime,style = 3)}

    # run simulation
    while(private$simTime >= private$tNow){

      self$oneDay()

      private$tNow = private$tNow + 1L
      if(verbose){setTxtProgressBar(pb,value = private$tNow)}
    }# end rest of sim

    ####################
    # close connections
    ####################
    close(private$conADM)
    close(private$conADF)

    if(verbose){cat("run ",private$runID," over\n",sep="")}

    ####################
    # reset everything
    ####################
    for(i in 1:private$nPatch){
      private$patches[[i]]$reset(verbose = verbose)
    }

    private$tNow = 2L
    private$runID = private$runID + 1L

  }# end repetition loop

}# end function

Network$set(which = "public",name = "multRun",
          value = multRun_Network, overwrite = TRUE
)

###############################################################################
# Auxiliary Run Functions
###############################################################################

#' Run a Single Day on a Network
#'
#' Runs a single day of simulation on a \code{\link{Network}} object, handling
#' population dynamics, migration, population update, and output.
#'
#' @keywords internal
#'
oneDay_Network <- function(){

  # intra-patch population dynamics
  for(i in 1:private$nPatch){
    private$patches[[i]]$oneDay_PopDynamics()
  }

  # inter-patch migration
  self$oneDay_Migration()

  # log output
  if(!(private$tNow %% private$sampTime)){
    for(i in 1:private$nPatch){
      private$patches[[i]]$oneDay_writeOutput()
    }
  }

}

Network$set(which = "public",name = "oneDay",
          value = oneDay_Network, overwrite = TRUE
)

#' Reset Network
#'
#' Reset a \code{\link{Network}} between runs, useful for Monte Carlo simulation.
#' This calls \code{\link{reset_Patch}} on each patch
#' and resets \code{tNow = 2} and increments the \code{runID}.
#'
#' @param verbose Chatty? Default is TRUE
#'
reset_Network <- function(verbose = TRUE){

  if(verbose){cat("reset network\n",sep="")}

  for(i in 1:private$nPatch){
    private$patches[[i]]$reset()
  }

  private$tNow = 2L
  private$runID = private$runID + 1L

}

Network$set(which = "public",name = "reset",
          value = reset_Network, overwrite = TRUE
)

