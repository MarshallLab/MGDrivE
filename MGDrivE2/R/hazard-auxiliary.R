################################################################################
#
#   MGDrivE2: auxiliary functions to calculate hazards, and parameter conversions
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   November 2019
#
################################################################################

################################################################################
# auxiliary fn's
################################################################################

#' Calculate Outbound Movement Rate
#'
#' Given \code{P}, the cumulative probability of moving before dying, and \code{mu},
#' the daily mortality rate, calculate the movement rate \code{gamma} to get \code{P}.
#' The equation comes from integrating the competing risks and solving for \code{gamma}.
#'
#' @param mu daily mortality rate
#' @param P cumulative probability to move before dying
#'
#' @return numeric probability of movement
#'
#' @examples
#'   # parameters, see vignette MGDrivE2: One Node Lifecycle Dynamics
#'   theta <- list(qE = 1/4, nE = 2, qL = 1/3, nL = 3, qP = 1/6, nP = 2,
#'                 muE = 0.05, muL = 0.15, muP = 0.05, muF = 0.09, muM = 0.09,
#'                 beta = 16, nu = 1/(4/24) )
#'
#'   # lets say a 70% chance to move over the entire lifespan
#'   rMoveRate <- calc_move_rate(mu = theta$muF, P = 0.70)
#'
#' @export
calc_move_rate <- function(mu,P){ return((P*mu) / (1 - P)) }


################################################################################
# convert tau, the stochastic movement matrix to a rate matrix
################################################################################

# prob to rate
stay_prob2rate <- function(tau){ return(-log(tau)) }

# take stochastic matrix tau and return a rate matrix Q but split into 2 parts
# gamma: rate to leave nodes (negative diagonal of the Q matrix)
# mat: the conditional probability matrix given something happens, what?

#' Convert Stochastic Matrix to Rate Matrix
#'
#' Given a stochastic matrix, return the rate matrix (infinitesimal generator)
#' that would generate it when exponentiated over the interval of unit time.
#'
#' Warning: if the matrix provided has diagonal-only rows (i.e., the location is
#' independent), the rate matrix will return 0 in that row, as there is no movement
#' rate that can generate that scenario.
#'
#' @param tau a row normalized stochastic matrix
#'
#' @return a list with two elements: \code{gamma} negative diagonal of the rate
#' matrix, \code{mat} matrix of row normalized off-diagonal elements
#'
#' @examples
#'   # generate random matrix for example
#'   #  This represents a 3-node landscape, with random movement between nodes
#'   moveMat <- matrix(data = runif(n = 9), nrow = 3, ncol = 3)
#'   moveMat <- moveMat/rowSums(moveMat)
#'
#'   moveRate <- movement_prob2rate(tau = moveMat)
#'
#' @export
movement_prob2rate <- function(tau){

  # check normalized
  if(any(rowSums(tau) != 1)){
    stop("'tau' has some rows that do not sum to 1, please normalize")
  }

  # calc gamma, leave rate
  gamma <- stay_prob2rate(diag(tau))

  # set diagonal to zero and renormalize
  diag(tau) <- 0
  tau <- tau/rowSums(tau)

  # protect if provided a diagonal matrix
  #  the "rate" to generate that matrix doesn't exist, but returning as 0
  tau[is.nan(tau)] <- 0

  # return list of 2 elements
  return(list("gamma"=gamma,
              "mat"=tau))
}
