################################################################################
#
#   MGDrivE2: Poisson time step sampler
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   October 2019
#
################################################################################

#' Make Poisson Time-Step (PTS) Sampler for a SPN Model
#'
#' Make a function closure to implement a Poisson time-step (tau-leaping with fixed tau)
#' sampler for a SPN.
#'
#' This sampling algorithm is based on representing a SPN as a set of competing
#' Poisson processes; it thus uses an integer valued state space but approximates
#' the number of events over \code{dt}.
#'
#' The design of \code{step_PTS} is from: Wilkinson, D. J. (2011). Stochastic
#' modeling for systems biology. CRC press
#'
#' Elements of the \code{N} list come from two places: The stoichiometry matrix
#' (\code{S}) is generated in \code{\link{spn_S}} and the hazards (\code{h}) come
#' from \code{\link{spn_hazards}}.
#'
#' For other samplers, see: \code{\link{step_CLE}}, \code{\link{step_DM}}, \code{\link{step_ODE}}
#'
#'
#' @param S a stoichiometry \code{\link[Matrix]{Matrix-class}} object
#' @param Sout an optional matrix to track of event firings
#' @param haz a list of hazard functions
#' @param dt time-step for tau-leap method
#' @param maxhaz maximum allowable hazard
#' @param human_ode ode used for human states
#'
#' @return function closure for use in \code{\link{sim_trajectory_R}} or \code{\link{sim_trajectory_CSV}}
#' @export
#' @importFrom stats rpois
#' @importFrom deSolve ode
step_PTS_decoupled <-
  function(S,
           Sout,
           haz,
           dt = 0.01,
           maxhaz = 1e6,
           human_ode = "SIS") {
    v = ncol(S)

    # if we are tracking things, this is dimension of tracking vector
    track <- FALSE
    if (!is.null(Sout)) {
      if (ncol(Sout) != v) {
        stop(
          "if providing output tracking matrix 'Sout' it must have same number of columns as stoichiometry matrix S"
        )
      }
      o <- nrow(Sout)
      track <- TRUE
    }

    if (human_ode == "SIS") {
      return(function(x0, h0, t0, deltat) {
        x <- x0 # state vector at t=t0
        human_states <- h0 #human states at t=0
        tNow <- t0
        termt <- t0 + deltat

        # tracking events
        if (track) {
          ovec <- rep(0, o) # output vetor at t=t0
        } else {
          ovec <- NULL
        }

        # sim loop
        repeat {
          # eval hazards
          h <- haz(x, tNow, human_states)
          if (any(h > maxhaz)) {
            stop("hazard too large, terminating simulation.\n\ttry reducing dt")
          }

          # sample event firings
          r <- rpois(n = v, lambda = h * dt)

          # update state and event tracking
          x <- x + as.vector(S %*% r)
          if (track) {
            ovec <- ovec + as.vector(Sout %*% r)
          }

          x[x < 0] <- 0 # absorption at zero
          tNow <- tNow + dt # update time

          # return condition
          if (tNow > termt) {
            return(list(
              "x" = x,
              "o" = ovec,
              "h" = human_states
            ))
          }
        } # end loop
      }
      )
    }
    return(function(x0, h0, t0, deltat, human_trace) {
      x <- x0 # state vector at t=t0
      human_states <- h0 #human states at t=0
      tNow <- t0
      termt <- t0 + deltat

      # tracking events
      if (track) {
        ovec <- rep(0, o) # output vetor at t=t0
      } else {
        ovec <- NULL
      }

      # sim loop
      repeat {
        # eval hazards
        h <- haz(x, tNow, human_states, human_trace)
        if (any(h > maxhaz)) {
          stop("hazard too large, terminating simulation.\n\ttry reducing dt")
        }

        # sample event firings
        r <- rpois(n = v, lambda = h * dt)

        # update state and event tracking
        x <- x + as.vector(S %*% r)
        if (track) {
          ovec <- ovec + as.vector(Sout %*% r)
        }

        x[x < 0] <- 0 # absorption at zero
        tNow <- tNow + dt # update time

        # return condition
        if (tNow > termt) {
          return(list(
            "x" = x,
            "o" = ovec,
            "h" = human_states
          ))
        }
      } # end loop

    }) # end function)
  }
