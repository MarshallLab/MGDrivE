################################################################################
#
#   MGDrivE2: CLE sampler
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   October 2019
#
################################################################################

#' Make Chemical Langevin (CLE) Sampler for a SPN model
#'
#' Make a function closure to implement a chemical Langevin (continuous-state)
#' approximation for a SPN.
#'
#' The chemical Langevin approximation is a numerical simulation of a Fokker-Planck
#' approximation to the Master equations (Kolmogorov Forwards Equations) governing
#' the stochastic model; the CLE approximation is a second-order approximation
#' that will get the correct mean and variance but higher order moments will be
#' incorrect.
#'
#' The design of \code{step_CLE} is from: Wilkinson, D. J. (2011). Stochastic
#' modeling for systems biology. CRC press
#'
#' Elements of the \code{N} list come from two places: The stoichiometry matrix
#' (\code{S}) is generated in \code{\link{spn_S}} and the hazards (\code{h}) come
#' from \code{\link{spn_hazards}}.
#'
#' For other samplers, see: \code{\link{step_PTS}}, \code{\link{step_DM}}, \code{\link{step_ODE}}
#'
#'
#' @param S a stoichiometry \code{\link[Matrix]{Matrix-class}} object
#' @param Sout an optional matrix to track of event firings. In the continuous stochastic model this will
#'        be the approximate cumulative intensity of each event.
#' @param haz a list of hazard functions
#' @param dt time-step for Euler-Maruyama method used to solve the SDE system
#' @param maxhaz maximum allowable hazard
#'
#' @return function closure for use in \code{\link{sim_trajectory_R}} or \code{\link{sim_trajectory_CSV}}
#'
#' @importFrom stats rnorm
step_CLE <- function(S,Sout,haz,dt=0.01,maxhaz=1e6){

  v = ncol(S)
  sdt = sqrt(dt)

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

           # initial state
           x <- x0
           tNow <- t0
           termt <- t0+deltat

           # tracking events
           if(track){
             ovec <- rep(0,o) # output vetor at t=t0
           } else {
             ovec <- NULL
           }

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

           # simulation loop
           repeat {

             h <- haz(x, tNow)
             if(any(h > maxhaz)){
               stop("hazard too large, terminating simulation.\n\ttry reducing dt")
             }

             # sample weiner process
             dw <- rnorm(n = v,mean = 0,sd = sdt)

             # update state and event tracking
             dx <- S %*% (h*dt + sqrt(h)*dw)
             x <- x + as.vector(dx)
             if(track){
               ovec <- ovec + as.vector(Sout %*% (h*dt + sqrt(h)*dw))
             }

             x[x<0] <- 0 # "absorption" at 0
             tNow <- tNow+dt

             # return condition
             if(tNow > termt){
               return(list("x"=x,"o"=ovec))
             }
           } # end loop
         } #end fucntion
  )
}
