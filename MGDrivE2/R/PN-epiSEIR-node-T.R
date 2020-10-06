################################################################################
#
#   MGDrivE2: SPN structure for a single node (SEI-SEIR epi)
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   March 2020
#
################################################################################

################################################################################
# make the transitions (T) of the SPN
################################################################################

#' Make Transitions (T) For a Node (SEI Mosquitoes - SEIR Humans)
#'
#' This function makes the set of transitions (T) for a SPN. It is used alone
#' if our model is a single-node metapopulation of mosquito and human dynamics;
#' otherwise it is used as part of other functions to make SPN models with larger
#' state spaces (metapopulation models, see \code{\link{spn_T_epiSEIR_network}}).
#'
#' This function takes the places produced from \code{\link{spn_P_epiSEIR_node}}
#' and builds all possible transitions between subsets of those places.
#'
#' The \code{params} argument supplies all of the ecological parameters necessary
#' to calculate equilibrium values. This function requires the \code{nE},
#' \code{nL}, \code{nP}, and \code{nEIP} parameters to be specified. For more details, see
#' \code{\link{equilibrium_SEI_SEIR}}
#'
#' While this function produces all structural information related to transitions,
#' hazards are produced by a separate function, \code{\link{spn_hazards}}.
#'
#' @param spn_P set of places produced by \code{\link{spn_P_epiSEIR_node}}
#' @param params a named list of parameters (see details)
#' @param cube an inheritance cube from the \code{MGDrivE} package (e.g. \code{\link[MGDrivE]{cubeMendelian}})
#'
#' @return a list with two elements: \code{T} contains transitions packets as lists,
#' \code{v} is the character vector of transitions (T)
#'
#' @examples
#' \dontrun{
#'   # parameters, see vignette MGDrivE2 SEIR Epidemiological Dynamics
#'   theta <- list(qE = 1/4, nE = 2, qL = 1/3, nL = 3, qP = 1/6, nP = 2,
#'                 muE = 0.05, muL = 0.15, muP = 0.05, muF = 0.09, muM = 0.09,
#'                 beta = 16, nu = 1/(4/24), NH = 250, X = c(1,0,0,0), NFX = 500,
#'                 f = 1/3, Q = 0.9, b = 0.55, c = 0.15, delta = 1/5, r = 1/14,
#'                 muH = 1/(62*365), qEIP = 1/11, nEIP = 3 )
#'
#'   # spn_P needs setup elsewhere, see vignettes
#'
#'   spn_T <- spn_T_epiSEIR_node(spn_P = spn_P, params = theta,
#'                               cube = MGDrivE::cubeMendelian() )
#' }
#'
#' @export
spn_T_epiSEIR_node <- function(spn_P,params,cube){

  # set of places
  u <- spn_P$u

  # time index
  T_index <- 1

  # return list of transitions and the vector
  retList <- spn_T_both_epiSEIR(u = u,nE = params$nE,nL = params$nL,nP = params$nP,
                                nEIP = params$nEIP,cube = cube,node_id = NULL,
                                T_index = T_index)


  # check the set for errors
  labels <- unlist(x = lapply(X = retList$T, FUN = '[[', 'label'), use.names = FALSE)
  indices <- unlist(x = lapply(X = retList$T, FUN = '[[', 'vix'), use.names = FALSE)
  if(any(labels != retList$v[indices])){
    stop(paste0("error in set of transitions T and transition vector at transition(s): ",
                paste0(retList$v[(labels != retList$v[indices])], collapse = ", ")))
  }


  # return the set of transitions and the vector (v)
  return(retList)
}


################################################################################
# TRANSITION FUNCTIONS
# helper functions that make the t in T (each transition in the full set)
# NOTE:
# these don't make the hazards yet; we have to make {P,T,Pre,Post}
# and then construct the hazards last before we are ready to simulate
# these encode the full structural information, {Pre,Post} are for easy
# analyepiSIS and computation.
#
#
#
################################################################################


################################################################################
# HUMANS
################################################################################

## births ---------------------------------------------------------------------

make_transition_human_birth_E_epiSEIR <- function(T_index,u,node=NULL){

  # tokens required
  input <- paste0(c("H","E",node),collapse = "_")

  # produces a new susceptible token
  output <- paste0(c("H","S",node),collapse = "_")

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(input,"->",output) # name of this t (corresponds to v)

  # requires a susceptible token
  t$s <- match(x = input, table = u)
  t$s_w <- 1

  # produces two tokens
  t$o <- match(x = c(output,input), table = u)
  t$o_w <- c(1,1)

  # class of the transition
  t$class <- "H_birth"

  # return the transition
  return(t)
}

make_transition_human_birth_R_epiSEIR <- function(T_index,u,node=NULL){

  # tokens required
  input <- paste0(c("H","R",node),collapse = "_")

  # produces a new susceptible token
  output <- paste0(c("H","S",node),collapse = "_")

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(input,"->",output) # name of this t (corresponds to v)

  # requires a susceptible token
  t$s <- match(x = input, table = u)
  t$s_w <- 1

  # produces two tokens
  t$o <- match(x = c(output,input),table = u)
  t$o_w <- c(1,1)

  # class of the transition
  t$class <- "H_birth"

  # return the transition
  return(t)
}

## deaths ---------------------------------------------------------------------

make_transition_human_death_E_epiSEIR <- function(T_index,u,node=NULL){

  # tokens required
  input <- paste0(c("H","E",node),collapse = "_")

  # no tokens are produced

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(input,"->D") # name of this t (corresponds to v)

  # requires a infected token
  t$s <- match(x = input, table = u)
  t$s_w <- 1

  # no tokens produced
  t$o <- NaN
  t$o_w <- NaN

  # class of the transition
  t$class <- "H_mort"

  # return the transition
  return(t)
}

make_transition_human_death_R_epiSEIR <- function(T_index,u,node=NULL){

  # tokens required
  input <- paste0(c("H","R",node),collapse = "_")

  # no tokens are produced

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(input,"->D") # name of this t (corresponds to v)

  # requires a infected token
  t$s <- match(x = input, table = u)
  t$s_w <- 1

  # no tokens produced
  t$o <- NaN
  t$o_w <- NaN

  # class of the transition
  t$class <- "H_mort"

  # return the transition
  return(t)
}

## infection ------------------------------------------------------------------

make_transition_human_inf_epiSEIR <- function(T_index,u,f_gen,m_gen,node=NULL){

  # tokens required: {1 infectious female, 1 susceptible human}
  ftoken_inf <- paste0(c("F",f_gen,m_gen,"I",node),collapse = "_")
  sus_h <- paste0(c("H","S",node),collapse = "_")

  # tokens produced: {1 infectious female, 1 infected (not yet infectious) human}
  exp_h <- paste0(c("H","E",node),collapse = "_")

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(sus_h,"->",exp_h) # name of this t (corresponds to v)

  # input arcs & multiplicity
  t$s <- match(x = c(ftoken_inf,sus_h),table = u)
  t$s_w <- c(1,1)

  # output arcs & multiplicity
  t$o <- match(x = c(ftoken_inf,exp_h),table = u)
  t$o_w <- c(1,1)

  # class of the transition
  t$class <- "H_infection"

  # return the transition
  return(t)
}

make_transition_human_latent_epiSEIR <- function(T_index,u,node=NULL){

  # tokens required
  input <- paste0(c("H","E",node),collapse = "_")

  # tokens produced
  output <- paste0(c("H","I",node),collapse = "_")

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(input,"->",output) # name of this t (corresponds to v)

  # requires a susceptible token
  t$s <- match(x = input, table = u)
  t$s_w <- 1

  # no tokens produced
  t$o <- match(x = output, table = u)
  t$o_w <- 1

  # class of the transition
  t$class <- "H_latent"

  # return the transition
  return(t)
}

make_transition_human_rec_epiSEIR <- function(T_index,u,node=NULL){

  # tokens required
  input <- paste0(c("H","I",node),collapse = "_")

  # tokens produced
  output <- paste0(c("H","R",node),collapse = "_")

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(input,"->",output) # name of this t (corresponds to v)

  # requires a susceptible token
  t$s <- match(x = input, table = u)
  t$s_w <- 1

  # no tokens produced
  t$o <- match(x = output, table = u)
  t$o_w <- 1

  # class of the transition
  t$class <- "H_recovery"

  # return the transition
  return(t)
}
