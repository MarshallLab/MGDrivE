########################################################################
#       _   __     __                      __
#      / | / /__  / /__      ______  _____/ /__
#     /  |/ / _ \/ __/ | /| / / __ \/ ___/ //_/
#    / /|  /  __/ /_ | |/ |/ / /_/ / /  / ,<
#   /_/ |_/\___/\__/ |__/|__/\____/_/  /_/|_|
#
#   Network Class Mosquito Population Simulation
#   Marshall Lab
#   November 2017
#
########################################################################

#' Reset Network
#'
#' Reset a \code{\link{Network}} between runs, useful for Monte Carlo simulation. This calls \code{\link{reset_Patch}} on each patch
#' and resets \code{tNow = 2} and increments the \code{runID}.
#'
reset_Network <- function(){

  cat("reset network\n",sep="")

  for(i in 1:private$nPatch){
    private$patches[[i]]$reset()
  }

  private$tNow = 2L
  private$runID = private$runID + 1L

}

Network$set(which = "public",name = "reset",
          value = reset_Network, overwrite = TRUE
)

#' Run Simulation
#'
#' Run a single simulation on this network.
#'
#' @param conADM an optional \code{\link[base]{connection}} to write male population dynamics to, if \code{NULL} use the directory specified in the constructor of \code{\link{Network}} with the current runID appended to the file.
#' @param conAF1 an optional \code{\link[base]{connection}} to write female population dynamics to, if \code{NULL} use the directory specified in the constructor of \code{\link{Network}} with the current runID appended to the file.
#'
oneRun_Network <- function(conADM = NULL, conAF1 = NULL){

  # open connections & write headers
  # parallel
  if(private$parameters$parallel){

    pid = Sys.getpid()
    if(is.null(conADM)){
      private$conADM = file(description = paste0(private$directory,
                                                 .Platform$file.sep,
                                                 "ADM_pid_",
                                                 pid,
                                                 "_Run",
                                                 formatC(x = private$runID, width = 6, format = "d", flag = "0"),
                                                 ".csv"),
                            open = "wt")
    } else {
      private$conADM = file(description = file.path(private$directory, conADM),open = "wt")
    }

    if(is.null(conAF1)){
      private$conAF1 = file(description = paste0(private$directory,
                                                 .Platform$file.sep,
                                                 "AF1_pid_",
                                                 pid,
                                                 "_Run",
                                                 formatC(x = private$runID, width = 6, format = "d", flag = "0"),
                                                 ".csv"),
                            open = "wt")
    } else {
      private$conAF1 = file(description = file.path(private$directory, conAF1),open = "wt")
    }

  # serial
  } else {
    private$conADM = file(description = paste0(private$directory,
                                               .Platform$file.sep,
                                               "ADM_Run",
                                               formatC(x = private$runID, width = 6, format = "d", flag = "0"),
                                               ".csv"),
                          open = "wt")
    private$conAF1 = file(description = paste0(private$directory,
                                               .Platform$file.sep,
                                               "AF1_Run",
                                               formatC(x = private$runID, width = 6, format = "d", flag = "0"),
                                               ".csv"),
                          open = "wt")
  }

  # males
  writeLines(text = paste0(c("Time","Patch",self$get_genotypesID()),collapse = ","),con = private$conADM,sep = "\n")

  # females
  femaleCrosses = c(t(outer(self$get_genotypesID(),self$get_genotypesID(),FUN = paste0)))
  writeLines(text = paste0(c("Time","Patch",femaleCrosses),collapse = ","),con = private$conAF1,sep = "\n")

  cat("begin run ",private$runID,"\n",sep="")

  # setup output
  for(i in 1:private$nPatch){
    private$patches[[i]]$oneDay_initOutput()
  }

  pb = txtProgressBar(min = 0,max = private$simTime,style = 3)

  while(private$simTime >= private$tNow){

    self$oneDay()

    private$tNow = private$tNow + 1L
    setTxtProgressBar(pb,value = private$tNow)
  }

  close(private$conADM)
  close(private$conAF1)

  cat("run ",private$runID," over\n",sep="")

}

Network$set(which = "public",name = "oneRun",
          value = oneRun_Network, overwrite = TRUE
)


#' Run a Single Day on a Network
#'
#' Runs a single day of simulation on a \code{\link{Network}} object, handling population dynamics, migration, population update, and output.
#'
oneDay_Network <- function(){

  # intra-patch population dynamics
  for(i in 1:private$nPatch){
    private$patches[[i]]$oneDay_PopDynamics()
  }

  # inter-patch migration
  self$oneDay_Migration()

  for(i in 1:private$nPatch){
    private$patches[[i]]$oneDay_updatePopulation()
  }

  for(i in 1:private$nPatch){
    private$patches[[i]]$oneDay_writeOutput()
  }

}

Network$set(which = "public",name = "oneDay",
          value = oneDay_Network, overwrite = TRUE
)
