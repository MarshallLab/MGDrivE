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
#'
#' @importFrom deSolve ode
step_ODE <- function(S,Sout,haz,method="lsoda"){

  # assign to local environment
  S <- S
  v <- ncol(S)
  u <- nrow(S)
  haz <- haz
  method <- method

  # dx/dt vector changes based on if tracking or not
  if(!is.null(Sout)){
    if(ncol(Sout) != v){
      stop(
        "if providing output tracking matrix 'Sout' it must have same number of columns as stoichiometry matrix S"
      )
    }
    track <- TRUE

    dxdt <- function(t,state,par=NULL){
      h <- vapply(X = haz,FUN = function(h){h(t=t,M=state)},FUN.VALUE = numeric(1),USE.NAMES = FALSE)
      list(
        (S %*% h)[,1],
        (Sout %*% h)[,1]
      )
    }
    nout <- nrow(Sout)
    xout <- 2:(u+1)
    oout <- (u+2):(u+2+nout-1)

  } else {
    track <- FALSE

    dxdt <- function(t,state,par=NULL){
      list((S %*% vapply(X = haz,FUN = function(h){h(t=t,M=state)},FUN.VALUE = numeric(1),USE.NAMES = FALSE))[,1])
    }
  }

  return(
         function(x0, t0, deltat){

           # solve ODEs over the step
           X <- ode(y = x0,times = c(t0,t0+deltat),func = dxdt,parms = NULL,method=method)

           if(track){
             return(
               list("x" = X[2,xout], "o" = X[2,oout])
             )
           } else {
             return(
               list("x" = X[2,-1], "o" = NULL)
             )
           }


         }
       )
}

# step_ODE <- function(S,haz,method="lsoda"){
#
#   # assign to local environment
#   S <- S
#   haz <- haz
#   method <- method
#
#   dxdt <- function(t,state,par=NULL){
#     list((S %*% vapply(X = haz,FUN = function(h){h(t=t,M=state)},FUN.VALUE = numeric(1),USE.NAMES = FALSE))[,1])
#   }
#
#   return(
#          function(x0, t0, deltat){
#            # solve ODEs over the step
#            X <- ode(y = x0,times = c(t0,t0+deltat),func = dxdt,parms = NULL,method=method)
#            return(
#              X[2,-1]
#            )
#          }
#        )
# }
