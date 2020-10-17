################################################################################
#
#   MGDrivE2: parameters & equilibrium-related functions (some user-facing)
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   December 2019
#
################################################################################

# not exported
# solve for aquatic mortality
# x: muA
fn_muA <- function(x,beta,mu,phi,qE,nE,qL,nL,qP,nP,Rm){
  return(Rm - ( (beta/mu) * (((qE*nE)/((qE*nE) + x))^nE) * (((qL*nL)/((qL*nL) + x))^nL) *
                  (((qP*nP)/((qP*nP) + x))^nP) * phi ) )
}

#' Solve for Constant Aquatic Mortality
#'
#' In \code{MGDrivE}, the model was typically solved at equilibrium assuming the
#' density-independent mortality was constant over aquatic stages (eggs, larvae, pupae),
#' given a daily growth rate, \eqn{r_{M}}. Given that growth rate, it solved for
#' that mortality \eqn{\mu_{Aqua}} by relating it with \eqn{R_{M}}, the per-generation
#' growth rate of the population, calculable from \eqn{r_{M}} and the mean
#' duration of life stages. This function uses \code{\link[stats]{uniroot}} to
#' solve for \eqn{mu_{Aqua}}.
#'
#' This function needs the following parameters in \code{params}:
#'  * \code{muF}: adult female mortality
#'  * \code{beta}: rate of egg laying
#'  * \code{phi}: sex ratio at emergence
#'  * \code{qE}: inverse of mean duration of egg stage
#'  * \code{nE}: shape parameter of Erlang-distributed egg stage
#'  * \code{qL}: inverse of mean duration of larval stage
#'  * \code{nL}: shape parameter of Erlang-distributed larval stage
#'  * \code{qP}: inverse of mean duration of pupal stage
#'  * \code{nP}: shape parameter of Erlang-distributed pupal stage
#'
#' @param params a named list of parameters
#' @param rm the daily growth rate
#'
#' @return location of the root, as provided from uniroot
#'
#' @importFrom stats uniroot
#'
#' @examples
#' theta <- list(qE = 1/4, nE = 2, qL = 1/5, nL = 3, qP = 1/6, nP = 2, muF = 1/12,
#'              beta = 32, phi = 0.5);
#' muAqatic <- solve_muAqua(params = theta, rm = 1.096)
#'
#' @export
solve_muAqua <- function(params,rm){

  # check
  stopifnot(all(c("muF","beta","phi","qE","nE","qL","nL","qP","nP") %in% names(params)))

  # per-generation growth rate
  g <- (1/params$qE)+(1/params$qL)+(1/params$qP)+(1/params$muF)
  Rm <- exp(log(rm)*g)

  # solve
  out <- uniroot(
    f=fn_muA,
    beta=params$beta,
    mu=params$muF,
    phi=params$phi,
    qE=params$qE,
    nE=params$nE,
    qL=params$qL,
    nL=params$nL,
    qP=params$qP,
    nP=params$nP,
    Rm=Rm,
    lower=.Machine$double.eps,
    upper=10
  )

  return(out$root)
}
