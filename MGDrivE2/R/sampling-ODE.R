################################################################################
#
#   MGDrivE2: ODE numerical integrator
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   October 2019
#
################################################################################


#' Make Mean-field Approximation (ODE) Numerical Integrator for a SPN Model
#'
#' Make a function closure to implement a first order mean-field ODE approximation
#' for a SPN.
#'
#' This method is equivalent to considering the ODEs describing the time
#' evolution of the mean trajectory (first moment) and setting all higher order
#' moments which appear on the right hand side to zero.
#'
#' The solvers used within can be found in the \code{deSolve} package, see
#' \code{\link[deSolve]{ode}}. For inhomogeneous systems, consider using the "rk4"
#' method to avoid excessive integration times.
#'
#' The stoichiometry matrix (\code{S}) is generated in \code{\link{spn_S}}.
#'
#' The list of hazards (\code{haz}) come from \code{\link{spn_hazards}}.
#'
#' For other samplers, see: \code{\link{step_CLE}}, \code{\link{step_PTS}}, \code{\link{step_DM}}
#'
#'
#' @param S a stoichiometry \code{\link[Matrix]{Matrix-class}} object
#' @param Sout an optional matrix to track of event firings. In the deterministic case it will return
#'        the rate of that event at the end of the time step
#' @param haz a list of hazard functions
#' @param method a character giving the type of numerical integrator used, the default is "lsoda"
#'
#' @return function closure for use in \code{\link{sim_trajectory_R}} or \code{\link{sim_trajectory_CSV}}
#' @export
#' @importFrom deSolve ode
step_ODE <- function(S, Sout, haz, method = "lsoda") {
  # assign to local environment
  S <- S
  v <- ncol(S)
  u <- nrow(S)
  haz <- haz
  method <- method

  # dx/dt vector changes based on if tracking or not
  if (!is.null(Sout)) {
    if (ncol(Sout) != v) {
      stop(
        "if providing output tracking matrix 'Sout' it must have same number of columns as stoichiometry matrix S"
      )
    }
    track <- TRUE

    dxdt <- function(t, state, par = NULL) {
      h <- haz(t = t, M = state)
      list((S %*% h)[, 1],
           (Sout %*% h)[, 1])
    }

    nout <- nrow(Sout)
    xout <- 2:(u + 1)
    oout <- (u + 2):(u + 2 + nout - 1)

  } else {
    track <- FALSE

    dxdt <- function(t, state, par = NULL) {
      h <- haz(t = t, M = state)
      list((S %*% h)[, 1])
    }

  }

  return(function(x0, t0, deltat) {
    # solve ODEs over the step
    X <-
      ode(
        y = x0,
        times = c(t0, t0 + deltat),
        func = dxdt,
        parms = NULL,
        method = method
      )

    if (track) {
      return(list("x" = X[2, xout], "o" = X[2, oout]))
    } else {
      return(list("x" = X[2, -1], "o" = NULL))
    }


  })
}



#' Make Mean-field Approximation (ODE) Numerical Integrator for a SPN Model for Decoupled Epi Dynamics
#'
#' Make a function closure to implement a first order mean-field ODE approximation
#' for a SPN.
#'
#' This method is equivalent to considering the ODEs describing the time
#' evolution of the mean trajectory (first moment) and setting all higher order
#' moments which appear on the right hand side to zero.
#'
#' The solvers used within can be found in the \code{deSolve} package, see
#' \code{\link[deSolve]{ode}}. For inhomogeneous systems, consider using the "rk4"
#' method to avoid excessive integration times.
#'
#' The stoichiometry matrix (\code{S}) is generated in \code{\link{spn_S}}.
#'
#' The list of hazards (\code{haz}) come from \code{\link{spn_hazards}}.
#'
#' For other samplers, see: \code{\link{step_CLE}}, \code{\link{step_PTS}}, \code{\link{step_DM}}
#'
#'
#' @param S a stoichiometry \code{\link[Matrix]{Matrix-class}} object
#' @param Sout an optional matrix to track of event firings. In the deterministic case it will return
#'        the rate of that event at the end of the time step
#' @param haz a list of hazard functions
#' @param method a character giving the type of numerical integrator used, the default is "lsoda"
#' @param human_ode ODE function used for human states
#'
#' @return function closure for use in \code{\link{sim_trajectory_R}} or \code{\link{sim_trajectory_CSV}}
#' @export
#' @importFrom deSolve ode
step_ODE_decoupled <-
  function(S,
           Sout,
           haz,
           method = "lsoda",
           human_ode = "SIS") {
    # assign to local environment
    S <- S
    v <- ncol(S)
    u <- nrow(S)
    haz <- haz
    method <- method
    # dx/dt vector changes based on if tracking or not
    if (!is.null(Sout)) {
      if (ncol(Sout) != v) {
        stop(
          "if providing output tracking matrix 'Sout' it must have same number of columns as stoichiometry matrix S"
        )
      }
      track <- TRUE

      if (human_ode == "SIS") {
        dxdt <- function(t, state, par = NULL) {
          h <- haz(t = t, M = state, h = par)
          list((S %*% h)[, 1],
               (Sout %*% h)[, 1])
        }
      } else {
        dxdt <- function(t, state, par = NULL) {
          h <- haz(
            t = t,
            M = state,
            h = par$h0,
            human_trace = par$human_trace
          )
          list((S %*% h)[, 1],
               (Sout %*% h)[, 1])
        }
      }

      nout <- nrow(Sout)
      xout <- 2:(u + 1)
      oout <- (u + 2):(u + 2 + nout - 1)

    } else {
      track <- FALSE

      if (human_ode == "SIS") {
        dxdt <- function(t, state, par = NULL) {
          h <- haz(t = t, M = state, h = par)
          list((S %*% h)[, 1])
        }
      } else {
        dxdt <- function(t, state, par = NULL) {
          h <- haz(
            t = t,
            M = state,
            h = par$h0,
            human_trace = par$human_trace
          )
          list((S %*% h)[, 1])
        }
      }



    }

    return(function(x0, h0, t0, deltat, human_trace = NULL) {
      # solve ODEs over the step
      X <-
        ode(
          y = x0,
          times = c(t0, t0 + deltat),
          func = dxdt,
          parms = list(h0 = h0, human_trace = human_trace) ,
          method = method
        )

      if (track) {
        return(list(
          "x" = X[2, xout],
          "o" = X[2, oout],
          "h" = h0
        ))
      } else {
        return(list(
          "x" = X[2, -1],
          "o" = NULL,
          "h" = h0
        ))
      }


    })
  }
