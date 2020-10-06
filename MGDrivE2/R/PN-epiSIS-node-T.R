################################################################################
#
#   MGDrivE2: SPN structure for a single node (SEI-SEIR epi)
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   October 2019
#
################################################################################

################################################################################
# make the transitions (T) of the SPN
################################################################################

#' Make Transitions (T) For a Node (SEI Mosquitoes - SIS Humans)
#'
#' This function makes the set of transitions (T) for a SPN. It is used alone if
#' our model is a single-node metapopulation of mosquito and human dynamics; otherwise
#' it is used as part of other functions to make SPN models with larger state
#' spaces (metapopulation models, see \code{\link{spn_T_epiSIS_network}}).
#'
#' This function takes the places produced from \code{\link{spn_P_epiSIS_node}}
#' and builds all possible transitions between subsets of those places.
#'
#' The \code{params} argument supplies all of the ecological parameters necessary
#' to calculate equilibrium values. This function requires the \code{nE},
#' \code{nL}, \code{nP}, and \code{nEIP} parameters to be specified. For more details, see
#' \code{\link{equilibrium_SEI_SIS}}
#'
#' While this function produces all structural information related to transitions,
#' hazards are produced by a separate function, \code{\link{spn_hazards}}.
#'
#' @param spn_P set of places produced by \code{\link{spn_P_epiSIS_node}}
#' @param params a named list of parameters (see details)
#' @param cube an inheritance cube from the \code{MGDrivE} package (e.g. \code{\link[MGDrivE]{cubeMendelian}})
#'
#' @return a list with two elements: \code{T} contains transitions packets as lists,
#' \code{v} is the character vector of transitions (T)
#'
#' @examples
#' \dontrun{
#'   # parameters, see vignette MGDrivE2: One Node Epidemiological Dynamics
#'   theta <- list(qE = 1/4, nE = 2, qL = 1/3, nL = 3, qP = 1/6, nP = 2,
#'                 muE = 0.05, muL = 0.15, muP = 0.05, muF = 0.09, muM = 0.09,
#'                 beta = 16, nu = 1/(4/24), NH = 1e3, X = 0.25, f = 1/3,
#'                 Q = 0.9, b = 0.55, c = 0.15, r = 1/200, muH = 1/(62*365),
#'                 qEIP = 1/11, nEIP = 6)
#'
#'   # spn_P needs setup elsewhere, see vignettes
#'
#'   spn_T <- spn_T_epiSIS_node(spn_P = spn_P, params = theta,
#'                              cube = MGDrivE::cubeMendelian() )
#' }
#'
#' @export
spn_T_epiSIS_node <- function(spn_P,params,cube){

  # set of places
  u <- spn_P$u

  # time index
  T_index <- 1

  # return list of places
  retList <- spn_T_both_epi(u = u,nE = params$nE,nL = params$nL,nP = params$nP,
                            nEIP = params$nEIP,cube = cube,node_id = NULL,
                            T_index = T_index)


  # check the set for errors
  labels <- unlist(x = lapply(X = retList$T, FUN = '[[','label'),use.names = FALSE)
  indices <- unlist(x = lapply(X = retList$T, FUN = '[[','vix'),use.names = FALSE)
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
# analysis and computation.
#
#
#
################################################################################


################################################################################
# PUPAE
################################################################################

# emergence to females event for pupae
# p_gen: genotype of the pupae
# m_gen: genotype of the male mate
# nP: number of pupae compartments
make_transition_pupae_emerge_f_epi <- function(T_index,u,p_gen,m_gen,nP,node=NULL){

  # tokens required
  ptoken <- paste0(c(paste0("P",nP),p_gen,node),collapse = "_")
  mtoken <- paste0(c("M",m_gen,node),collapse = "_")

  # produces a female token (with her genotype + the genotype of her mate + susceptible)
  # also return the male token.
  ftoken <- paste0(c("F",p_gen,m_gen,"S",node),collapse = "_")

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(ptoken,"->",ftoken) # name of this t (corresponds to v)

  # requires a pupae token and a male token
  t$s <- match(x = c(ptoken,mtoken), table = u)
  t$s_w <- c(1,1)

  # produces a female and a male token
  t$o <- match(x = c(ftoken,mtoken), table = u)
  t$o_w <- c(1,1)

  # class of the transition
  t$class <- "pupae_2f"

  # return the transition
  return(t)
}


################################################################################
# OVIPOSITION (NEED TO OVERWRITE BECAUSE FEMALES HAVE ADDITIONAL INFECTION INDEXING)
################################################################################

################################################################################
# FEMALES
################################################################################

# infect a female token
# inf: infection
make_transition_female_inf_epi <- function(T_index,u,f_gen,m_gen,node=NULL){

  # tokens required: 1 susceptible female, 1 infected human
  inf_h <- paste0(c("H_I",node),collapse = "_")
  ftoken <- paste0(c("F",f_gen,m_gen,"S",node),collapse = "_")

  # produces an incubating female token
  ftoken_inf <- paste0(c("F",f_gen,m_gen,"E1",node),collapse = "_")

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(ftoken,"->",ftoken_inf) # name of this t (corresponds to v)

  # requires a female token
  t$s <- match(x = c(ftoken,inf_h), table = u)
  t$s_w <- c(1,1)

  # produces one infected token
  t$o <- match(x = c(ftoken_inf,inf_h), table = u)
  t$o_w <- c(1,1)

  # class of the transition
  t$class <- "female_inf"

  # return the transition
  return(t)
}

# a female token advances incubation
# inc: incubation
# inc1: stage 1
# inc2: stage 2
make_transition_female_eip_epi <- function(T_index,u,f_gen,m_gen,inc1,inc2,node=NULL){

  # tokens required
  ftoken <- paste0(c("F",f_gen,m_gen,paste0("E",inc1),node),collapse = "_")

  # produces a token
  ftoken_inf <- paste0(c("F",f_gen,m_gen,paste0("E",inc2),node),collapse = "_")

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(ftoken,"->",ftoken_inf) # name of this t (corresponds to v)

  # requires a female token
  t$s <- match(x = ftoken, table = u)
  t$s_w <- 1

  # produces one infected token
  t$o <- match(x = ftoken_inf, table = u)
  t$o_w <- 1

  # class of the transition
  t$class <- "female_eip"

  # return the transition
  return(t)
}

# a female token finishes incubation and is now infectious
# inc: incubation
# nEIP: max bin of EIP
make_transition_female_inc_epi <- function(T_index,u,f_gen,m_gen,nEIP,node=NULL){

  # tokens required
  ftoken <- paste0(c("F",f_gen,m_gen,paste0("E",nEIP),node),collapse = "_")

  # produces an infected female token
  ftoken_inf <- paste0(c("F",f_gen,m_gen,"I",node),collapse = "_")

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(ftoken,"->",ftoken_inf) # name of this t (corresponds to v)

  # requires a female token
  t$s <- match(x = ftoken, table = u)
  t$s_w <- 1

  # produces one infected token
  t$o <- match(x = ftoken_inf, table = u)
  t$o_w <- 1

  # class of the transition
  t$class <- "female_inc"

  # return the transition
  return(t)
}


################################################################################
# HUMANS
################################################################################

## births ---------------------------------------------------------------------

make_transition_human_birth_S_epiSIS <- function(T_index,u,node=NULL){

  # tokens required
  #  produced token is the same, as it produces 2 susceptibles
  sus <- paste0(c("H","S",node),collapse = "_")

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(sus,"->",sus) # name of this t (corresponds to v)

  # requires a susceptible token
  t$s <- match(x = sus, table = u)
  t$s_w <- 1

  # produces two suceptible tokens
  t$o <- match(x = sus, table = u)
  t$o_w <- 2

  # class of the transition
  t$class <- "H_birth"

  # return the transition
  return(t)
}

# we assume no vertical transmission of diseases we are considering, so I's produce S's
make_transition_human_birth_I_epiSIS <- function(T_index,u,node=NULL){

  # tokens required
  inf <- paste0(c("H","I",node),collapse = "_")

  # produces a new susceptible token
  new_sus <-paste0(c("H","S",node),collapse = "_")

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(inf,"->",new_sus) # name of this t (corresponds to v)

  # requires a susceptible token
  t$s <- match(x = inf, table = u)
  t$s_w <- 1

  # produces two suceptible tokens
  t$o <- match(x = c(new_sus,inf), table = u)
  t$o_w <- c(1,1)

  # class of the transition
  t$class <- "H_birth"

  # return the transition
  return(t)
}

## deaths ---------------------------------------------------------------------

make_transition_human_death_S_epiSIS <- function(T_index,u,node=NULL){

  # tokens required
  sus <- paste0(c("H","S",node),collapse = "_")

  # no tokens are produced

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(sus,"->D") # name of this t (corresponds to v)

  # requires a susceptible token
  t$s <- match(x = sus, table = u)
  t$s_w <- 1

  # no tokens produced
  t$o <- NaN
  t$o_w <- NaN

  # class of the transition
  t$class <- "H_mort"

  # return the transition
  return(t)
}

make_transition_human_death_I_epiSIS <- function(T_index,u,node=NULL){

  # tokens required
  inf <- paste0(c("H","I",node),collapse = "_")

  # no tokens are produced

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(inf,"->D") # name of this t (corresponds to v)

  # requires a infected token
  t$s <- match(x = inf, table = u)
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

make_transition_human_inf_epiSIS <- function(T_index,u,f_gen,m_gen,node=NULL){

  # tokens required: {1 infectious female, 1 susceptible human}
  ftoken_inf <- paste0(c("F",f_gen,m_gen,"I",node),collapse = "_")
  sus_h <- paste0(c("H","S",node),collapse = "_")

  # tokens produced: {1 infectious female, 1 infected human}
  inf_h <- paste0(c("H","I",node),collapse = "_")

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(sus_h,"->",inf_h) # name of this t (corresponds to v)

  # input arcs & multiplicity
  t$s <- match(x = c(ftoken_inf,sus_h), table = u)
  t$s_w <- c(1,1)

  # output arcs & multiplicity
  t$o <- match(x = c(ftoken_inf,inf_h), table = u)
  t$o_w <- c(1,1)

  # class of the transition
  t$class <- "H_infection"

  # return the transition
  return(t)
}

make_transition_human_rec_epiSIS <- function(T_index,u,node=NULL){

  # tokens required
  inf <- paste0(c("H","I",node),collapse = "_")

  # tokens produced
  sus <- paste0(c("H","S",node),collapse = "_")

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(inf,"->",sus) # name of this t (corresponds to v)

  # input arcs & multiplicity
  t$s <- match(x = inf, table = u)
  t$s_w <- 1

  # output arcs & multiplicity
  t$o <- match(x = sus, table = u)
  t$o_w <- 1

  # class of the transition
  t$class <- "H_recovery"

  # return the transition
  return(t)
}
