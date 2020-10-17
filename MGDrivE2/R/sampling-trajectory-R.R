################################################################################
#
#   MGDrivE2: trajectory interface
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   October 2019
#
################################################################################

################################################################################
# User-facing simulation trajectory
################################################################################

#' Simulate Trajectory From a SPN Model
#'
#' This function provides a unified interface to the various simulation algorithms
#' for SPN, returning output sampled at a lattice of time points to the user, and
#' handling various exogenous events that may occur during the simulation
#' (such as release of adult mosquitoes).
#'
#' \code{dt_stoch} is used by the Poisson Time-Step (\code{\link{step_PTS}}) and
#' Chemical Langevin (\code{\link{step_CLE}}) methods to approximate the hazards.
#' A smaller \code{dt_stoch} provides a better approximation, but will take longer
#' to run.
#'
#' The stoichiometry matrix (\code{S}) is generated in \code{\link{spn_S}}.
#'
#' The list of hazards (\code{hazards}) come from \code{\link{spn_hazards}}.
#'
#' Several samplers are provided. The default is a Poisson Time-Step
#' (\code{\link{step_PTS}}) method. Other options are Gillespie's Direct Method
#' (\code{\link{step_DM}}) and a Chemical Langevin sampler (\code{\link{step_CLE}}).
#' Additionally, for convenience, an ODE "sampler" (\code{\link{step_ODE}}) is
#' provided for compatibility with other samplers. This function uses methods from
#' \code{deSolve}.
#'
#' If using the \code{ode} sampler, several \code{methods} are provided in the \code{deSolve}
#' package, see \code{\link[deSolve]{ode}}. For inhomogeneous systems, consider
#' using the "rk4" method to avoid excessive integration times.
#'
#' Additionally, \code{events} objects must follow the format required by
#' \code{deSolve}. This was done for consistency, see \code{\link[deSolve]{events}}
#' for more information.
#'
#' This function tracks state variables by default; an optional argument \code{Sout}
#' can be provided to track number of event firings each time step (for discrete stochastic simulations),
#' or cumulative intensity (for continuous stochastic simulations), or the rate function of
#' particular events for ODE simulation. The matrix must have number of columns equal to
#' number of events in the system (the number of hazard functions), and a row for each tracking
#' variable. The function \code{\link{track_hinf}} is provided, which builds a matrix to track
#' human infection events.
#'
#' To save output as .csv files, see \code{\link{sim_trajectory_CSV}}.
#'
#'
#' @param x0 the initial marking of the SPN (initial state, M0)
#' @param t0 initial time to begin simulation
#' @param tt the final time to end simulation
#' @param dt the time-step at which to return output (\strong{not} the time-step of the sampling algorithm)
#' @param dt_stoch time-step used for approximation of hazards
#' @param num_reps number of repetitions to run, default is 1.
#' @param S a stoichiometry \code{\link[Matrix]{Matrix-class}} object
#' @param hazards list of hazard functions
#' @param Sout an optional matrix to track event firings
#' @param sampler determines sampling algorithm, one of; "ode", "tau", "cle", or "dm"
#' @param method if \code{sampler} is "ode", the solver to use, from \code{deSolve}
#' @param events a \code{data.frame} of events
#' @param verbose print a progress bar?
#' @param ... further named arguments passed to the step function
#'
#' @return a list with 2 elements: "state" is the array of returned state values, and "events" will
#'        return events tracked with \code{Sout} if provided, otherwise is \code{NULL}
#'
#' @export
sim_trajectory_R <- function(x0, t0=0, tt=100, dt=1, dt_stoch = 0.1, num_reps=1,
                             S, hazards, Sout = NULL, sampler = "tau", method = "lsoda",
                             events = NULL, verbose = TRUE,...){

  # check/setup step function
  stepFun <- base_stepFunc(sampler = sampler,S = S,hazards = hazards,Sout = Sout,
                           dt_stoch = dt_stoch,method = method,...)

  # check/setup simulation time
  simTimes <- base_time(t0 = t0,tt = tt,dt = dt)

  # check/organize events and x0
  events <- base_events(x0 = x0, events = events, dt = dt)

  # pass everything down to base function
  #  base function does return
  sim_trajectory_base_R(x0 = switch(EXPR = sampler, "ode" = x0, round(x0)),
                        times = simTimes,dt = dt, num_reps = num_reps, stepFun = stepFun,
                        events = events, Sout = Sout, verbose = verbose)
}


################################################################################
# Base Functions
################################################################################
#######################################
# Step Funtion
#######################################

# This sets up the step function for sampling
# It's the exact same in sim_trajectory_CSV, so using function
# sample = string denoting type of sampler
# S = stoichiometry matrix
# hazards = list of hazard functions and a flag
# dt_stoch = time-step for approximate hazards
# method = sampler for ODE solver
#
base_stepFunc <- function(sampler,S,hazards,Sout = NULL,dt_stoch,method,...){

  ##########
  # checks for step function!
  ##########
  #  sampler
  if(!(sampler %in% c("ode","tau","cle","dm")) ){
    stop("Sampler must be one of: ode, tau, cle, or dm")
  }

  #  check hazards, depending on sampler
  if((sampler %in% c("ode","cle")) && hazards$flag){
    warning(paste0("For numerical stability, it is strongly recommended that ",
                   "ode and cle samplers use approximate hazards."))
  }

  #  warning to use sparse matrix
  if(attr(class(S),"package") != "Matrix"){
    warning("the stoichiometry 'S' matrix is not sparse.\n significant speed issues possible")
  }


  ##########
  # setup step function!
  #########
  #  do the deterministic, then stochastic
  if(sampler == "ode"){
    stepFun <- step_ODE(S=S, Sout = Sout, haz = hazards$hazards, method = method, ...)

  } else {
    # setup hazard function for stochastic things, deterministic one has it's own thing
    hazFunc <- function(M,t){vapply(X = hazards$hazards,FUN = function(h){h(t=t,M=M)},
                                    FUN.VALUE = numeric(1),USE.NAMES = FALSE) }

    # generate stepFun
    if(sampler == "tau"){
      stepFun <- step_PTS(S=S, Sout = Sout, haz=hazFunc, dt = dt_stoch, ...)
    } else if(sampler == "cle"){
      stepFun <- step_CLE(S=S, Sout = Sout, haz=hazFunc, dt = dt_stoch, ...)
    } else if(sampler == "dm"){
      stepFun <- step_DM(S=S, Sout = Sout, haz=hazFunc, ...)
    } else {
      stop("option 'sampler' must be one of 'tau', 'cle', 'dm', 'ode'")
    }

  } # end stepFun setup

  # return step function
  return(stepFun)
}


#######################################
# Time Funtion
#######################################

# check time, and setup sampling times
# t0 = initial time to begin simulation
# tt = the final time to end simulation
# dt = the time-step at which to return output
#
base_time <- function(t0,tt,dt){

  # number of steps we need to take
  n <- (tt-t0) %/% dt + 1
  times <- seq(from=t0,to=tt,by=dt)

  if(length(times) != n){
    stop("error in sequence of times; make sure tt is evenly divisible by dt")
  }

  # return
  return(times)
}


#######################################
# Events Funtion
#######################################

# check and organize events
#  this function also checks the names on x0, for lack of a better place
# x0 = the initial marking of the SPN (initial state)
# events = a \code{data.frame} of events
# dt = the time-step at which to return output
#
base_events <- function(x0, events, dt){

  if(is.null(events)){
    return(NULL)
  }

  # check x0 names
  xNames <- names(x0)
  if(any(is.null(xNames))){
    stop("Marking(s) on x0 is null, please check and try again.")
  }

  #  check events naming
  if(!all(events$var %in% xNames)){
    stop("Events ", paste(events$var[!(events$var %in% xNames)],collapse = ", "),
         " are not valid places in the network.\nPlease check the naming.")
  }


  # if there are events:
  #  check event timing
  #  sort events by time
  if(nrow(events) > 0){

    # make sure events occur at times that will actually be sampled
    if(any(events$time %% dt != 0)){
      warning("event times do not correspond to a multiple of dt.\n",
      "event times will be rounded up to the nearest time-step!")
      events$time = events$time + (events$time %% dt)
    }

    # sort it
    events <- events[order(events$time),]

  } # end check/sort

  # use integer position instead of string matching
  events$var_id <- match(x = events$var, table = names(x0))


  # return updated events
  return(events)
}

################################################################################
# Base simulation trajectory
################################################################################

#' Simulate Trajectory From one SPN Model
#'
#' This is an internal function to \code{\link{sim_trajectory_R}}. It does the
#' actual sampling once all of the functions have been checked and setup.
#'
#'
#' @param x0 the initial marking of the SPN (initial state)
#' @param times sequence of sampling times
#' @param dt the time-step at which to return output (\strong{not} the time-step of the sampling algorithm)
#' @param num_reps number of repetitions to run
#' @param stepFun a sampling function
#' @param events a \code{data.frame} of events (uses the same format as required in package \code{deSolve} for consistency, see \code{\link[deSolve]{events}} for more information)
#' @param Sout an optional matrix to track event firings
#' @param verbose print a progress bar?
#'
#' @return matrix of sampled values
sim_trajectory_base_R <- function(x0, times, dt=1, num_reps, stepFun, events=NULL, Sout = NULL, verbose=TRUE){

  # setup return array
  nTime <- length(times)
  retArray <- array(data = 0, dim = c(nTime, length(x0)+1, num_reps),
                    dimnames = list(times, c("time",names(x0)), 1:num_reps))
  retArray[ ,1, ] <- times
  retArray[1,-1, ] <- x0

  # set up event tracking (tracks differences, so one less time output than state)
  if(!is.null(Sout)){
    track <- TRUE
    ret_events <- array(
      data = 0, dim = c(nTime-1, nrow(Sout)+1, num_reps), dimnames = list(times[-1], c("time",rownames(Sout)), 1:num_reps)
    )
    ret_events[,1,] <- times[-1]
  } else {
    track <- FALSE
    ret_events <- NULL
  }

  # loop over num_reps, calling base function
  for(r in 1:num_reps){

    # set things for reseting later
    state <- list("x"=NULL,"o"=NULL)
    state$x <- x0
    repEvents <- events

    # progress bar
    if(verbose){
      pbar <- txtProgressBar(min = 1,max = nTime,style = 3)
      cat(" --- begin simulation --- \n")
    }

    # main simulation loop
    for(i in 2:nTime){

      # iterate the step function until this delta t is over
      tNow <- times[i]
      state <- stepFun(state$x,tNow,dt)

      if(all(fequal(state$x,0))){
        if(verbose){close(pbar)}
        warning(" --- marking of net is zero; terminating simulation early --- \n")
        break
      }

      # add the event to the state vector
      #  would it be better to get times of events, then add them in at times
      #  instead of checking every day if there are events, and then if it isthe day
      if(!is.null(repEvents)){
        while((nrow(repEvents) > 0) && fequal(repEvents[1,"time"],tNow)){
          state$x[repEvents[1,"var_id"]] <- state$x[repEvents[1,"var_id"]] + repEvents[1,"value"]
          repEvents <- repEvents[-1,]
        }
      }

      # record output
      retArray[i,-1,r] <- state$x
      if(track){
        ret_events[i-1,-1,r] <- state$o
      }

      # progress bar
      if(verbose){setTxtProgressBar(pb = pbar,value = i)}
    } # end sim loop

    if(verbose){
      close(pbar)
      cat(" --- end simulation --- \n")
    }

  } # end repetition loop

  # return stuff
  return(list("state"=retArray,"events"=ret_events))
}
