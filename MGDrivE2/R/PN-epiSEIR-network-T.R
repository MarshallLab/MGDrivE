################################################################################
#
#   MGDrivE2: SPN structure for a metapopulation network (SEI-SIS epi)
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   November 2019
#
################################################################################

################################################################################
#   Compose transitions to create the metapopulation transition structure
################################################################################

#' Make Transitions (T) For a Network (SEI Mosquitoes - SEIR Humans)
#'
#' This function makes the set of transitions (T) for a SPN model of a
#' metapopulation network for simulation of coupled SEI-SEIR dynamics. It is the
#' network version of \code{\link{spn_T_epiSEIR_node}}.
#'
#' This function takes the places produced from \code{\link{spn_P_epiSEIR_network}}
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
#' For larger networks, this function may take some time to return, please be patient;
#' the Petri Net modeling formalism trades additional computation time at model
#' initialization for faster sampling of trajectories within a simulation.
#'
#' Please note, the movement matrices (\code{h_move} and \code{m_move}) are NOT
#' stochastic matrices, just binary matrices that say if i,j can exchange population.
#' Diagonal elements must be \code{FALSE}, and both matrices are checked for validity; the
#' function will stop with errors if the adjacency matrix specifies illegal movement
#' rules (e.g.; mosquito movement from a "h" node to a "b" node)
#'
#' For examples of using this function, see:
#' \code{vignette("seir-dynamics", package = "MGDrivE2")}
#'
#' @param node_list a character vector specifying what type of nodes to create;
#' (m = a node with only mosquitoes, h = a node with only humans, b = a node with both humans and mosquitoes)
#' @param spn_P set of places produced by \code{\link{spn_P_epiSEIR_network}}
#' @param params a named list of parameters (see details)
#' @param cube an inheritance cube from the \code{MGDrivE} package (e.g. \code{\link[MGDrivE]{cubeMendelian}})
#' @param h_move binary adjacency matrix indicating if movement of humans between nodes is possible or not
#' @param m_move binary adjacency matrix indicating if movement of mosquitoes between nodes is possible or not
#'
#' @return a list with two elements: \code{T} contains transitions packets as lists,
#' \code{v} is the character vector of transitions (T)
#'
#' @export
spn_T_epiSEIR_network <- function(node_list,spn_P,params,cube,h_move,m_move){

  # set of places
  u <- spn_P$u

  # genotypes
  nG <- cube$genotypesN
  g <- cube$genotypesID

  # infection stages
  epi_stages <- c("S",paste0("E",as.character(1:params$nEIP)),"I")

  # within node transitions
  T_meta <- vector("list",length(node_list))
  T_index <- 1

  # make all of the transitions in each node
  for(t in 1:length(node_list)){
    if(node_list[t]=="b"){
      T_meta[[t]] <- spn_T_both_epiSEIR(u = u,nE = params$nE,nL = params$nL,nP = params$nP,
                                        nEIP = params$nEIP,cube = cube,node_id = t,
                                        T_index = T_index)
    } else if(node_list[t]=="h"){
      T_meta[[t]] <- spn_T_humans_epiSEIR(u = u,node_id = t,T_index = T_index)
    } else if(node_list[t]=="m"){
      T_meta[[t]] <- spn_T_mosy_epi(u = u,nE = params$nE,nL = params$nL,nP = params$nP,
                                    nEIP = params$nEIP,cube = cube,node_id = t,
                                    T_index = T_index)
    } else {
      stop(paste0("error: bad entry in node_list, ",node_list[t]))
    }
  }


  # make mosquito movement events
  # This checks for valid movement, then sets lists for male/female mosquitoes,
  #  and human movement
  if(any(Matrix::diag(h_move))){
    stop("adjacency matrix 'h_move' must have all FALSE elements along its diagonal")
  }
  if(any(Matrix::diag(m_move))){
    stop("adjacency matrix 'm_move' must have all FALSE elements along its diagonal")
  }
  moveList <- spn_T_move_epi(node_list = node_list,u = u,m_move = m_move,
                             h_move = h_move,nG = nG,g = g,nEIP = params$nEIP,
                             epi_stages = epi_stages,
                             h_state = c("S","E","I","R"),T_index = T_index)


  # set of transitions with movement appended to it
  v <- c(unlist(x = lapply(X = T_meta, FUN = '[[', 'v'), use.names = FALSE),
         unlist(x = lapply(X = moveList, FUN = lapply, '[[', 'label'),
                use.names = FALSE))

  # all the transitions (and flatten them)
  T_all <- c(list("T_win" = do.call(c,lapply(X = T_meta,FUN = '[[', 'T'))),
            moveList)

  T_all <- unlist(x = T_all, recursive = FALSE, use.names = FALSE)


  # check the set for errors
  labels <- unlist(x = lapply(X = T_all, FUN = '[[', 'label'), use.names = FALSE)
  indices <- unlist(x = lapply(X = T_all, FUN = '[[', 'vix'), use.names = FALSE)
  if(any(labels != v[indices])){
    stop(paste0("error in set of transitions T and transition vector at transition(s): ",
                paste0(v[(labels != v[indices])], collapse = ", ")))
  }


  # return the complete set of transitions {T}
  return(list("T" = T_all,
              "v" = v) )
}


################################################################################
# WITHIN NODE INTERNAL TRANSITION FUNCTIONS
#
# NOTE:
# these functions make the within node transitions for the T set.
# they are not exported from the package for users
# as the package provides support only to generate Petri Nets by
# returning complete sets (P,T); if the below functions were exported
# it would be possible for a user to return an incomplete set of places (P)
# or transitions (T)
#
################################################################################

################################################################################
# make internal transitions (T) for human-only node
################################################################################
# base SEIR human function
#  used in spn_T_humans_epiSEIR and spn_T_both_epiSEIR
# only ~40% similar to base_T_humans_epi
base_T_humans_epiSEIR <- function(u,node_id,T_index){

  # HUMAN TRANSITIONS
  human_tt <- vector("list",10)
  vv <- 1

  # births
  human_tt[[vv]] <- make_transition_human_birth_S_epiSIS(T_index,u=u,node=node_id)
  T_index <- T_index + 1
  vv <- vv + 1

  human_tt[[vv]] <- make_transition_human_birth_E_epiSEIR(T_index,u=u,node=node_id)
  T_index <- T_index + 1
  vv <- vv + 1

  human_tt[[vv]] <- make_transition_human_birth_I_epiSIS(T_index,u=u,node=node_id)
  T_index <- T_index + 1
  vv <- vv + 1

  human_tt[[vv]] <- make_transition_human_birth_R_epiSEIR(T_index,u=u,node=node_id)
  T_index <- T_index + 1
  vv <- vv + 1

  # deaths
  human_tt[[vv]] <- make_transition_human_death_S_epiSIS(T_index,u=u,node=node_id)
  T_index <- T_index + 1
  vv <- vv + 1

  human_tt[[vv]] <- make_transition_human_death_E_epiSEIR(T_index,u=u,node=node_id)
  T_index <- T_index + 1
  vv <- vv + 1

  human_tt[[vv]] <- make_transition_human_death_I_epiSIS(T_index,u=u,node=node_id)
  T_index <- T_index + 1
  vv <- vv + 1

  human_tt[[vv]] <- make_transition_human_death_R_epiSEIR(T_index,u=u,node=node_id)
  T_index <- T_index + 1
  vv <- vv + 1

  # latency and recovery
  human_tt[[vv]] <- make_transition_human_latent_epiSEIR(T_index,u=u,node=node_id)
  T_index <- T_index + 1
  vv <- vv + 1

  human_tt[[vv]] <- make_transition_human_rec_epiSEIR(T_index,u=u,node=node_id)
  T_index <- T_index + 1
  vv <- vv + 1


  # set T_index in parent environment
  #  ie, update the counter
  assign(x = "T_index", value = T_index, pos = parent.frame())

  return(human_tt)
}


# u: set of places
# if this is for node 1, T_index = 1, otherwise its max(node[i-1].T_index)+1
# node_id: id of the node these transitions are for
spn_T_humans_epiSEIR <- function(u,node_id,T_index){

  # HUMAN TRANSITIONS
  human_tt <- base_T_humans_epiSEIR(u = u,node_id = node_id,T_index = T_index)

  # transitions (v)
  v <- unlist(x = lapply(X = human_tt, FUN = '[[', 'label'), use.names = FALSE)

  # set T_index in parent environment
  #  ie, update the counter
  assign(x = "T_index", value = T_index, pos = parent.frame())

  # return the set of transitions and the vector (v)
  return(list("T" = human_tt,
              "v" = v) )
}


################################################################################
# make internal transitions (T) for human & mosquito node
################################################################################

# u: set of places
# if this is for node 1, T_index = 1, otherwise its max(node[i-1].T_index)+1
spn_T_both_epiSEIR <- function(u,nE,nL,nP,nEIP,cube,node_id,T_index){

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
      female_inf_tt[[vv]] <- make_transition_female_inf_epi(T_index,u=u,f_gen=g[j_f],m_gen=g[j_m],node=node_id)
      T_index <- T_index + 1
      vv <- vv + 1
    }
  }


  # HUMAN TRANSITIONS
  human_tt <- vector("list",nG^2)
  vv <- 1

  # add epi transitions not in base
  for(j_f in 1:nG){
    for(j_m in 1:nG){
      human_tt[[vv]] <- make_transition_human_inf_epiSEIR(T_index,u=u,f_gen=g[j_f],m_gen=g[j_m],node=node_id)
      T_index <- T_index + 1
      vv <- vv + 1
    }
  }

  # combine with base transitions
  human_tt <- c(human_tt,
                base_T_humans_epiSEIR(u = u, node_id = node_id, T_index = T_index))


  # the set of transitions
  t <- c(base_mos,
         "female_inf" = list(female_inf_tt),
         "human" = list(human_tt))

  # transitions (v)
  v <- unlist(x = lapply(X = t, FUN = lapply, '[[', 'label'), use.names = FALSE)

  # one long vector
  t <- unlist(x = t, recursive = FALSE, use.names = FALSE)


  # set T_index in parent environment
  #  ie, update the counter
  assign(x = "T_index", value = T_index, pos = parent.frame())


  # return the set of transitions and the vector (v)
  return(list("T" = t,
              "v" = v) )
}
