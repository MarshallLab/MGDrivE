################################################################################
#
#   MGDrivE2: SPN structure for a single node (SEI-SIS epi) in Decoupled sampling framework
#   Marshall Lab
#   Agastya Mondal (agastya_mondal@berkeley.edu)
#   February 2021
#
################################################################################

################################################################################
# make the transitions (T) of the SPN
################################################################################

#' Make Transitions (T) For a Node (SEI Mosquitoes)
#'
#' This function makes the set of transitions (T) for a SPN. It is used alone if
#' our model is a single-node metapopulation of mosquito; otherwise
#' it is used as part of other functions to make SPN models with larger state
#' spaces (metapopulation models, see \code{\link{spn_T_epiSIS_network}}).
#'
#' This function takes the places produced from \code{\link{spn_P_epi_decoupled_node}}
#' and builds all possible transitions between subsets of those places.
#'
#' The \code{params} argument supplies all of the ecological parameters necessary
#' to calculate equilibrium values. This function requires the \code{nE},
#' \code{nL}, \code{nP}, and \code{nEIP} parameters to be specified. For more details, see
#' \code{\link{equilibrium_SEI_SIS}}
#'
#' While this function produces all structural information related to transitions,
#' hazards are produced by a separate function, \code{\link{spn_hazards}}.
#' This is used by both decoupled SIS and Imperial transmission model sampling.
#' For examples of using this function, see:
#' \code{vignette("epi-node-decoupled", package = "MGDrivE2")}
#'
#' @param spn_P set of places produced by \code{\link{spn_P_epi_decoupled_node}} function
#' @param params a named list of parameters (see details)
#' @param cube an inheritance cube from the \code{MGDrivE} package (e.g. \code{\link[MGDrivE]{cubeMendelian}})
#'
#' @return a list with two elements: \code{T} contains transitions packets as lists,
#' \code{v} is the character vector of transitions (T)
#'
#' @export
spn_T_epi_decoupled_node <- function(spn_P,params,cube){

  # set of places
  u <- spn_P$u

  # time index
  T_index <- 1

  # return list of places
  retList <- spn_T_both_epi_decoupled(u = u,nE = params$nE,nL = params$nL,nP = params$nP,
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

spn_T_both_epi_decoupled <- function(u,nE,nL,nP,nEIP,cube,node_id,T_index){

  # genetic states
  g <- cube$genotypesID
  nG <- cube$genotypesN

  epi_stages <- c("S",paste0("E",as.character(1:nEIP)),"I")

  # MOSQUITO TRANSITIONS
  # base things in a list
  base_mos <- base_T_mosy_epi(u = u, nE = nE, nL = nL, nP = nP, nEIP = nEIP,
                               cube = cube, node_id = node_id, T_index = T_index,
                               epi_stages = epi_stages)

  # make female infection
  female_inf_tt <- vector("list",nG^2)
  vv <- 1

  for(j_f in 1:nG){
    for(j_m in 1:nG){
      female_inf_tt[[vv]] <- make_transition_female_inf_epi_decoupled(T_index,u=u,f_gen=g[j_f],
                                                            m_gen=g[j_m],node=node_id)
      T_index <- T_index + 1
      vv <- vv + 1
    }
  }

  # the set of transitions
  t <- c(base_mos,
         "female_inf" = list(female_inf_tt))

  # keep track of the female infection transitions
  # since these rely on human states, we will need them later
  # when evaluating hazards

  inf_labels <- unlist(sapply(female_inf_tt, function(inf) { return(inf$label) }))

  # transitions (v)
  v <- unlist(x = lapply(X = t, FUN = lapply, '[[', 'label'), use.names = FALSE)

  # one long vector
  t <- unlist(x = t, recursive = FALSE, use.names = FALSE)


  # set T_index in parent environment
  #  ie, update the counter
  assign(x = "T_index", value = T_index, pos = parent.frame())

  # return the set of transitions and the vector (v)
  return(list("T" = t,
              "v" = v,
              "inf_labels" = inf_labels) )
}



################################################################################
# FEMALES
################################################################################

# infect a female token
# inf: infection
make_transition_female_inf_epi_decoupled <- function(T_index,u,f_gen,m_gen,node=NULL){

  # tokens required: 1 susceptible female, human tokens will be handled in
  # corresponding hazard function
  ftoken <- paste0(c("F",f_gen,m_gen,"S",node),collapse = "_")

  # produces an incubating female token
  ftoken_inf <- paste0(c("F",f_gen,m_gen,"E1",node),collapse = "_")

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
  t$class <- "female_inf"

  # return the transition
  return(t)
}

