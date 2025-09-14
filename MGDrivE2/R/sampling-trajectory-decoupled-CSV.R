################################################################################
#
#   MGDrivE2: trajectory interface for decoupled sampling
#   Marshall Lab
#   Agastya Mondal (agastya_mondal@berkeley.edu)
#   September 2021
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
#' (such as release of adult mosquitoes). This function is used in decoupled sampling, where
#' the mosquito and human states are separated. This is used primarily when using the
#' Imperial model of malaria transmission.
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
#' This function writes all output to .csv files. Each simulation is written to
#' a \code{folder} element - the number of repetitions is the number of folders
#' provided. For now, only adult mosquito states, human states, clinical incidence, and pathogen prevalence are
#' written to CSVs.
#'
#' This function tracks state variables specified by argument \code{stage} by default; an optional argument \code{Sout}
#' can be provided to track number of event firings each time step (for discrete stochastic simulations),
#' or cumulative intensity (for continuous stochastic simulations), or the rate function of
#' particular events for ODE simulation. The matrix must have number of columns equal to
#' number of events in the system (the number of hazard functions), and a row for each tracking
#' variable. If \code{Sout} is provided, it outputs an additional csv, "Tracking.csv".
#' The function \code{\link{track_hinf}} is provided, which builds a matrix to track
#' human infection events.
#'
#' To return simulations to R for further processing, see \code{\link{sim_trajectory_R}}.
#'
#'
#' @param x0 the initial marking of the SPN (initial state, M0)
#' @param h0 the initial human state distribution
#' @param inf_labels labels corresponding to female mosquito infection hazard
#' @param tmax the final time to end simulation
#' @param dt the time-step at which to return output (\strong{not} the time-step of the sampling algorithm)
#' @param dt_stoch time-step used for approximation of hazards
#' @param folders vector of folders to write output
#' @param S a stoichiometry \code{\link[Matrix]{Matrix-class}} object
#' @param hazards list of hazard functions
#' @param Sout an optional matrix to track event firings
#' @param sampler determines sampling algorithm, one of; "ode", "tau", "cle", or "dm"
#' @param method if \code{sampler} is "ode", the solver to use, from \code{deSolve}
#' @param events a \code{data.frame} of events
#' @param batch a \code{list} of batch migration events, created from \code{\link[MGDrivE2]{batch_migration}}, may be set to \code{NULL} if not used
#' @param verbose print a progress bar?
#' @param SPN_P stochastic petri net places
#' @param theta parameters
#' @param human_ode human ode function
#' @param cube inheritance cube
#' @param ... further named arguments passed to the step function
#'
#' @return NULL - prints output to .csv files
#'
#' @importFrom stats rbinom
#'
#' @export
sim_trajectory_CSV_decoupled <- function(x0,
                                         h0,
                                         inf_labels,
                                         tmax,
                                         dt = 1,
                                         dt_stoch = 0.1,
                                         folders = "./",
                                         S,
                                         hazards,
                                         SPN_P,
                                         theta,
                                         Sout = NULL,
                                         sampler = "tau",
                                         method = "lsoda",
                                         events = NULL,
                                         batch = NULL,
                                         verbose = TRUE,
                                         human_ode = "Imperial",
                                         cube = cube,
                                         ...) {
  if (sampler == "ode-decoupled" & !is.null(batch)) {
    stop("batch migration is incompatible with deterministic simulations")
  }

  if(human_ode != "Imperial") {
    stop("CSV writing is only supported for decoupled Imperial sampling. Check human_ode parameter.")
  }

  stepFun <- base_stepFunc_decoupled(
    sampler = sampler,
    S = S,
    hazards = hazards,
    Sout = Sout,
    dt_stoch = dt_stoch,
    method = method,
    inf_labels = inf_labels,
    human_ode=human_ode,
    ...
  )

  # check/setup simulation time
  simTimes <- base_time_decoupled(tt = tmax, dt = dt)

  # check/organize events and x0
  events <- base_events_decoupled(x0 = x0, events = events, dt = dt)

  # pass everything down to base function - writes CSVs
  sim_trajectory_base_CSV_decoupled(
    x0 = switch(EXPR = sampler, "ode-decoupled" = x0, round(x0)),
    h0 = h0,
    SPN_P = SPN_P,
    theta = theta,
    times = simTimes,
    stepFun = stepFun,
    events0 = events,
    batch = batch,
    Sout = Sout,
    verbose = verbose,
    human_ode = human_ode,
    cube = cube,
    folders = folders
  )
}


################################################################################
# Base simulation trajectory
################################################################################

#' Simulate Trajectory From one  SPN Model
#'
#' This is an internal function to \code{\link{sim_trajectory_CSV}}. It does the
#' actual sampling once all of the functions have been checked and setup.
#'
#' @param x0 the initial marking of the SPN (initial state)
#' @param h0 initial human state distribution
#' @param SPN_P stochastic petri net, places
#' @param theta parameters
#' @param human_ode ode function used for human states
#' @param cube inheritance cube
#' @param times sequence of sampling times
#' @param stepFun a sampling function
#' @param folders vector of folders to write output
#' @param events0 a \code{data.frame} of events (uses the same format as required
#' in package \code{deSolve} for consistency, see \code{\link[deSolve]{events}}
#' for more information)
#' @param batch a \code{list} of batch migration events, created from \code{\link[MGDrivE2]{batch_migration}}, may be set to \code{NULL} if not used
#' @param Sout an optional matrix to track event firings
#' @param verbose print a progress bar?
#'
#' @return no return, prints .csv files into provided folders
sim_trajectory_base_CSV_decoupled <- function(x0,
                                              h0,
                                              SPN_P,
                                              theta,
                                              times,
                                              stepFun,
                                              events0 = NULL,
                                              batch = NULL,
                                              Sout = NULL,
                                              verbose = TRUE,
                                              human_ode = "Imperial",
                                              cube = NULL,
                                              folders = folders) {
  # stuff for later
  xNames <- names(x0)
  nTime <- length(times)
  na <- theta$na #number of age compartments
  clin_inc <- theta$clin_inc
  mort <- theta$mort_eq

  human_state_labels <- generate_Imperial_human_state_labels(na)

  # make index list for all possible output of mosquito epi states
  # humans will be handled differently
  pStuff <- base_print_stuff_decoupled_mosy(state_names = xNames)
  pList <- pStuff$pList
  stage <- pStuff$stage
  lenPS <- length(stage)
  # track rates/events?
  if (!is.null(Sout)) {
    track <- TRUE
  } else {
    track <- FALSE
  }

  # loop over repetitions
  for (num_reps in 1:length(folders)) {
    # progress bar
    if (verbose) {
      pbar <- txtProgressBar(min = 1,
                             max = nTime,
                             style = 3)
      cat(" --- begin simulation --- \n")
    }

    if (human_ode == "Imperial") {
      # initialize human trace list
      names <- as.character(seq(1, 15))
      human_trace <- vector("list", length(names))
      names(human_trace) <- names

      human_trace[["1"]] <- h0
    }

    # reset at the beginning of every simulation
    events <- events0
    repBatch <- batch

    state <- list("x" = NULL, "o" = NULL, "h" = NULL)
    state$x <- x0
    state$h <- h0

    # Initialize files
    fileCons <-
      setNames(object = vector(mode = "list", length = lenPS + 1),
               nm = c(stage, "H"))

    # mosquito files
    for (curS in 1:lenPS) {
      # open file connection
      fileCons[[curS]] <-
        file(
          description = paste0(folders[num_reps], .Platform$file.sep, stage[curS], ".csv"),
          open = "wt"
        )
      # write header
      writeLines(
        text = paste0(c("Time", xNames[pList[[curS]]]), collapse = ","),
        con = fileCons[[curS]],
        sep = "\n"
      )
      # write T0
      writeLines(
        text = paste0(c(0, x0[pList[[curS]]]), collapse = ","),
        con = fileCons[[curS]],
        sep = "\n"
      )
    }


    # human file
    fileCons[["H"]] <- file(
      description = paste0(folders[num_reps], .Platform$file.sep, "H", ".csv"),
      open = 'wt'
    )
    writeLines(text = paste0(c("Time", human_state_labels), collapse = ","),
               con = fileCons[["H"]],
               sep = "\n")
    writeLines(text = paste0(c(0, as.vector(h0), clin_inc, mort), collapse = ","),
               con = fileCons[["H"]],
               sep = "\n")


    # tracking rates/event firing
    if (track) {
      event_con <- file(
        description = paste0(folders[num_reps], .Platform$file.sep, "Tracking.csv"),
        open = "wt"
      )
      writeLines(text = paste0(c("Time", rownames(Sout)), collapse = ","),
                 con = event_con,
                 sep = "\n")
    }

    # main simulation loop
    for (i in 2:nTime) {
      # iterate the step function until this delta t is over
      t0 <- times[i - 1]
      t1 <- times[i]
      dt <- t1 - t0

      state <- stepFun(state$x, state$h, t0, dt, human_trace)
      mosquito_counts <- aggregate_female_SEI(SPN_P, state$x)

      if (human_ode != "Imperial") {
        stop("CSV output is only supported for human_ode='Imperial'")
      }

      func <- human_Imperial_ODE
      # attach distribution of infectious mosquitoes across genotypes
      # and total number of adult mosquitoes
      theta$I_V <- mosquito_counts$I_V
      theta$N_V <- sum(mosquito_counts$N_V)

      out <- ode(
        y = state$h,
        times = c(t0, t1),
        func = func,
        parms = theta
      )
      # update human state matrix
      y <- tail(out, 1)
      num_states <- 12
      # keep clinical incidence separate, it doesn't contribute to ODE dynamics
      # we're just using it for presentation/viz later on
      human_state_matrix <-
        matrix(y[2:length(y)], ncol = num_states)
      clin_inc_idx <- 11
      mort_idx <- 12

      clin_inc <- human_state_matrix[, clin_inc_idx]
      mort <- human_state_matrix[, mort_idx]
      human_state_matrix <- human_state_matrix[, -c(clin_inc_idx, mort_idx)]
      colnames(human_state_matrix) <-
        c("S", "T", "D", "A", "U", "P", "ICA", "IB", "ID", "IVA")

      state$h <- human_state_matrix
      idx <- as.character(t1 + 1)
      human_trace[[idx]] <- human_state_matrix

      # only keep the last 15 values in memory, otherwise program becomes too slow
      max_entries <- 15
      if (length(human_trace) > max_entries) {
        num_entries <- length(human_trace)
        human_trace <-
          human_trace[(num_entries - 15):num_entries]
      }


      if (all(fequal(state$x, 0))) {
        if (verbose) {
          close(pbar)
        }

        # fill rest of file with 0s
        #  this makes sure all files are the same size/shape for analysis
        for (nT in i:nTime) {
          for (curS in 1:lenPS) {
            writeLines(
              text = paste0(c(times[nT], state$x[pList[[curS]]]), collapse = ","),
              con = fileCons[[curS]],
              sep = "\n"
            )
          }
          # human output
          writeLines(
            text = paste0(c(times[nT], rep(
              0, length(human_state_labels)
            )), collapse = ","),
            con = fileCons[["H"]],
            sep = ","
          )
          # tracking rates/event firing
          if (track) {
            writeLines(
              text = paste0(c(times[nT], rep(
                0, nrow(Sout)
              )), collapse = ","),
              con = event_con,
              sep = "\n"
            )
          }
        } # end print loop

        warning(" --- marking of net is zero; terminating simulation early --- \n")
        break
      }

      # add the event to the state vector
      if (!is.null(events)) {
        while ((nrow(events) > 0) && events[1, "time"] <= t1) {
          state$x[events[1, "var_id"]] <-
            state$x[events[1, "var_id"]] + events[1, "value"]
          events <- events[-1, ]
        }
      }

      # batch migration?
      if (!is.null(repBatch)) {
        while (length(repBatch) > 0 && repBatch[[1]]$time <= t1) {
          # how many go?
          moving <-
            stats::rbinom(
              n = length(repBatch[[1]]$from),
              size = as.integer(state$x[repBatch[[1]]$from]),
              prob = repBatch[[1]]$prob
            )
          # update the state
          state$x[repBatch[[1]]$from] <-
            state$x[repBatch[[1]]$from] - moving
          state$x[repBatch[[1]]$to] <-
            state$x[repBatch[[1]]$to] + moving
          repBatch <- repBatch[-1]
        }
      }

      # record output
      for (curS in 1:lenPS) {
        writeLines(
          text = paste0(c(t1, state$x[pList[[curS]]]), collapse = ","),
          con = fileCons[[curS]],
          sep = "\n"
        )
      }
      writeLines(
        text = paste0(c(t1, as.vector(y[2:length(y)])), collapse = ','),
        con = fileCons[["H"]],
        sep = "\n"
      )
      # tracking rates/event firing
      if (track) {
        writeLines(
          text = paste0(c(t1, state$o), collapse = ","),
          con = event_con,
          sep = "\n"
        )
      }

      # progress bar
      if (verbose) {
        setTxtProgressBar(pb = pbar, value = i)
      }
    } # end sim loop

    # close connections
    for (curS in 1:lenPS) {
      close(con = fileCons[[curS]])
    }
    close(con = fileCons[["H"]])
    if (track) {
      close(con = event_con)
    }


    if (verbose) {
      close(pbar)
      cat(" --- end simulation --- \n")
    }
  } # end repetitions loop

  # no return
}


################################################################################
# Base print
################################################################################

# there was a lot to make sure that everything gets printed properly,
#  so moving it down here
base_print_stuff_decoupled_mosy <- function(state_names) {
  # make index list for all possible output
  # for now we pre-select the adult mosquito states and all human states

  # First get mosquito states
  stage <- c("F", "M")
  pList <-
    lapply(
      X = stage,
      FUN = function(x) {
        grep(pattern = paste0("^", x),
             x = state_names)
      }
    )

  # check for silly user request
  nullIDX <- lengths(pList) == 0
  # warn them about it
  if (any(nullIDX)) {
    warning(paste0(
      stage[nullIDX],
      " is not present in the sim, removing from output request.\n"
    ))
  }

  # remove it
  stage <- stage[!nullIDX]
  pList[nullIDX] <- NULL
  # split female into SEI epi states
  # get index of female
  fIDX <- match(x = "F", table = stage)

  # list of split female names
  fSplitList <- strsplit(x = state_names[pList[[fIDX]]],
                         split = "_",
                         fixed = TRUE)

  # mosquito-only case
  if (is.na(fSplitList[[1]][4]) || (fSplitList[[1]][4] != "S")) {
    # just change female name to match SEI output later
    stage[stage == "F"] <- "FS"

  } else {
    # get index of female, remove from stage
    #  hold on to pList so I can use it
    stage <- stage[-fIDX]


    # get list of female SEI indices
    #  unlist states, then grep them (can't match because E states)
    #  then return index as list
    infStat <-
      vapply(X = fSplitList,
             FUN = '[[',
             4,
             FUN.VALUE = character(length = 1))
    fList <- lapply(
      X = c("S", "E", "I"),
      FUN = function(x) {
        pList[[fIDX]][grep(pattern = x,
                           x = infStat,
                           fixed = TRUE)]
      }
    )

    # now remove old female indices from pList
    pList[fIDX] <- NULL

    # append things and return updated pStages and pList
    stage <- c(stage, "FS", "FE", "FI")
    pList <- c(pList, fList)

  } # end update for F_S/E/I

  # return 2 things in list
  return(list("stage" = stage,
              "pList" = pList))
}
