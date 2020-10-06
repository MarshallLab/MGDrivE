################################################################################
#
#   MGDrivE2: Direct-method sampler
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   October 2019
#
################################################################################

#' Make Gillespie's Direct Method (DM) Sampler for a SPN model
#'
#' Make a function closure to implement Gillespie's Direct Method sampler for a SPN.
#'
#' The direct method is an exact sampling algorithm; it simulates each event individually.
#' Because of this it may be extremely slow for non-trivial population sizes, and
#' thus should be used to debug and test rather than for serious Monte Carlo simulation.
#'
#' The design of \code{step_DM} is from: Wilkinson, D. J. (2011). Stochastic
#' modeling for systems biology. CRC press
#'
#' Elements of the \code{N} list come from two places: The stoichiometry matrix
#' (\code{S}) is generated in \code{\link{spn_S}} and the hazards (\code{h}) come
#' from \code{\link{spn_hazards}}.
#'
#' For other samplers, see: \code{\link{step_CLE}}, \code{\link{step_PTS}}, \code{\link{step_ODE}}
#'
#'
#' @param S a stoichiometry \code{\link[Matrix]{Matrix-class}} object
#' @param Sout an optional matrix to track of event firings
#' @param haz a list of hazard functions
#' @param maxhaz maximum allowable hazard
#'
#' @return function closure for use in \code{\link{sim_trajectory_R}} or \code{\link{sim_trajectory_CSV}}
#'
#' @importFrom stats rexp
step_DM <- function(S,Sout,haz,maxhaz=1e6){

  v = ncol(S)

  # if we are tracking things, this is dimension of tracking vector
  track <- FALSE
  if(!is.null(Sout)){
    if(ncol(Sout) != v){
      stop(
        "if providing output tracking matrix 'Sout' it must have same number of columns as stoichiometry matrix S"
      )
    }
    o <- nrow(Sout)
    track <- TRUE
  }

  return(
    function(x0, t0, deltat){

      # initial state and time
      x <- x0
      tNow <- t0
      termt <- t0+deltat

      # tracking event firings
      if(track){
        ovec <- rep(0,o) # output vetor at t=t0
      } else {
        ovec <- NULL
      }

      repeat {
        # evaluate hazards
        h <- haz(x,tNow)
        h0 <- sum(h)
        if(h0 < 1e-10){
          return(x)
        } else if(h0 > maxhaz){
          stop("hazard too large, terminating simulation.")
        } else {
          tNow <- tNow + rexp(n = 1,rate = h0)
        }

        # return condition
        if(tNow >= termt){
          return(list("x"=x,"o"=ovec))
        }

        j <- sample(x = v,size = 1,prob=h)
        x <- x + S[,j]
        if(track){
          ovec <- ovec + Sout[,j]
        }
      } # end loop
    } # end function
  )
}
