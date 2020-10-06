################################################################################
#
#   MGDrivE2: SPN structure for a network (lifecycle only)
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   October 2019
#
################################################################################

################################################################################
# make the transitions (T) of the SPN
################################################################################

#' Make Transitions (T) For a Network (Mosquitoes only)
#'
#' This function makes the set of transitions (T) for a SPN model of a
#' metapopulation network. It is the network version of \code{\link{spn_T_lifecycle_node}}.
#'
#' This function takes the places produced from \code{\link{spn_P_lifecycle_network}}
#' and builds all possible transitions between subsets of those places.
#'
#' The \code{params} argument supplies all of the ecological parameters necessary
#' to calculate equilibrium values. This function requires the \code{nE},
#' \code{nL}, and \code{nP} parameters to be specified. For more details, see
#' \code{\link{equilibrium_lifeycle}}
#'
#' While this function produces all structural information related to transitions,
#' hazards are produced by a separate function, \code{\link{spn_hazards}}.
#'
#' For larger networks, this function may take some time to return, please be patient;
#' the Petri Net modeling formalism trades additional computation time at model
#' initialization for faster sampling of trajectories within a simulation.
#'
#' Please note, the movement matrix (\code{m_move}) is NOT a
#' stochastic matrices, just a binary matrix that say if i,j can exchange population.
#' Diagonal elements must be \code{FALSE}.
#'
#' @param spn_P set of places produced by \code{\link{spn_P_lifecycle_network}}
#' @param params a named list of parameters (see details)
#' @param cube an inheritance cube from the \code{MGDrivE} package (e.g. \code{\link[MGDrivE]{cubeMendelian}})
#' @param m_move binary adjacency matrix indicating if movement of mosquitoes between nodes is possible or not
#'
#' @return a list with two elements: \code{T} contains transitions packets as lists,
#' \code{v} is the character vector of transitions (T)
#'
#' @importFrom Matrix diag which
#' @importFrom utils tail
#'
#' @examples
#' \dontrun{
#'   # parameters, see vignette MGDrivE2: One Node Lifecycle Dynamics
#'   theta <- list(qE = 1/4, nE = 2, qL = 1/3, nL = 3, qP = 1/6, nP = 2,
#'                 muE = 0.05, muL = 0.15, muP = 0.05, muF = 0.09, muM = 0.09,
#'                 beta = 16, nu = 1/(4/24) )
#'
#'   # spn_P needs setup elsewhere, see vignettes
#'
#'   # setup 2-node allowed migrations
#'   #  let node 1 migrate to node 2, but not vice-versa
#'   m_move <- matrix(data = FALSE, nrow = 2, ncol = 2)
#'   m_move[1,2] <- TRUE
#'
#'   spn_T <- spn_T_lifecycle_network(spn_P = spn_P, params = theta,
#'                                    cube = MGDrivE::cubeMendelian(),
#'                                    m_move = m_move )
#' }
#'
#' @export
spn_T_lifecycle_network <- function(spn_P,params,cube,m_move){

  # set of places
  u <- spn_P$u

  # genotypes
  nG <- cube$genotypesN
  g <- cube$genotypesID

  # within node transitions
  T_meta <- vector("list",nrow(m_move))
  T_index <- 1

  # make all of the transitions in each node
  for(id in 1:nrow(m_move)){
    # only 1 type of nodes
    T_meta[[id]] <- spn_T_mosy_lifecycle(u = u, nE = params$nE, nL = params$nL,
                                         nP = params$nP, cube = cube, node_id = id,
                                         T_index = T_index)
  } # end loop over nodes


  # make mosquito movement events
  # checks first
  if(any(Matrix::diag(m_move))){
    stop("adjacency matrix 'm_move' must have all FALSE elements along its diagonal")
  }

  # make transition events
  # edges: need to make mosquito movement for these edges
  edges <- arrayInd(Matrix::which(m_move),.dim=dim(m_move))
  edges <- edges[order(edges[,1]),]

  # make female movement
  female_move <- vector("list",nrow(edges)*(nG^2))
  vv <- 1

  # iterate over edges
  for(e in 1:nrow(edges)){
    # ... and over female genotypes
    for(f in 1:nG){
      # ... and over male mate genotypes
      for(m in 1:nG){
        female_move[[vv]] <- make_transition_move_female(T_index = T_index,u = u,
                                                         f_gen = g[f],m_gen = g[m],
                                                         origin = edges[e,1],
                                                         dest = edges[e,2])
        T_index <- T_index + 1
        vv <- vv + 1
      }
    }
  }

  # make male movement
  male_move <- vector("list",nrow(edges)*nG)
  vv <- 1

  # iterate over edges
  for(e in 1:nrow(edges)){
    # ... and over male genotypes
    for(m in 1:nG){
      male_move[[vv]] <- make_transition_move_male(T_index = T_index,u = u,
                                                   m_gen = g[m],origin = edges[e,1],
                                                   dest = edges[e,2])
      T_index <- T_index + 1
      vv <- vv + 1
    }
  }


  # set of transitions
  v <- c(unlist(x = lapply(X = T_meta, FUN = '[[', 'v'), use.names = FALSE),
         unlist(x = lapply(X = c(female_move, male_move), FUN = '[[', 'label'),
                use.names = FALSE))


  # all the transitions
  T_all <- list("T_win" = do.call(c,lapply(X = T_meta,FUN = '[[', 'T')),
                "female_move" = female_move,
                "male_move" = male_move)

  # # flatten it for testing (maybe use this later?)
  T_all <- unlist(x = T_all, recursive = FALSE, use.names = FALSE)


  # check the set for errors
  labels <- unlist(x = lapply(X = T_all, FUN = '[[','label'), use.names = FALSE)
  indices <- unlist(x = lapply(X = T_all, FUN = '[[','vix'), use.names = FALSE)
  if(any(labels != v[indices])){
    stop(paste0("error in set of transitions T and transition vector at transition(s): ",
                paste0(v[(labels != v[indices])], collapse = ", ")))
  }


  return(list("T" = T_all,
              "v" = v) )
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
#   some transition functions for movement
################################################################################

# move one male
make_transition_move_male <- function(T_index,u,m_gen,origin,dest){

  # tokens required
  mtoken1 <- file.path("M_",m_gen,"_",origin,fsep = "")

  # tokens produced
  mtoken2 <- file.path("M_",m_gen,"_",dest,fsep = "")

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- file.path(mtoken1,"->",mtoken2,fsep = "") # name of this t (corresponds to v)

  # requires one male token
  t$s <- match(x = mtoken1,table = u)
  t$s_w <- 1

  # produces one male token
  t$o <- match(x = mtoken2,table = u)
  t$o_w <- 1

  # class of the transition
  t$class <- "move_male"

  # return the transition
  return(t)
}

# move one female
make_transition_move_female <- function(T_index,u,f_gen,m_gen,origin,dest,inf=NULL){

  # tokens required/produced
  ftoken1_2 <- file.path(paste0(c("F",f_gen,m_gen,inf),collapse = "_"),
                         c(origin,dest), fsep = "_")

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(ftoken1_2, collapse = "->") # name of this t (corresponds to v)

  # requires one female token
  t$s <- match(x = ftoken1_2[1], table = u)
  t$s_w <- 1

  # produces one female token
  t$o <- match(x = ftoken1_2[2], table = u)
  t$o_w <- 1

  # class of the transition
  t$class <- "move_female"

  # return the transition
  return(t)
}


################################################################################
# make base internal transitions (T) for other functions below
################################################################################

# base mosy function
#  used in spn_T_mosy, and spn_T_mosy_epi
base_T_mosy <- function(u,nE,nL,nP,nG,g,node_id = NULL,T_index){

  # EGG TRANSITIONS

  # make egg mortality transitions
  egg_mort_tt <- vector("list",nE*nG)
  vv <- 1

  for(i in 1:nE){
    for(j in 1:nG){
      egg_mort_tt[[vv]] <- make_transition_egg_mort(T_index,u=u,e_gen=g[j],stage=i,node=node_id)
      T_index <- T_index + 1
      vv <- vv + 1
    }
  }

  # make egg advancement transitions
  egg_adv_tt <- vector("list",nE*nG)
  vv <- 1

  for(i in 1:nE){
    for(j in 1:nG){
      if(i == nE){
        egg_adv_tt[[vv]] <- make_transition_egg_adv(T_index,u=u,e_gen=g[j],stage1=i,
                                                    stage2=NULL,node=node_id)
      } else {
        egg_adv_tt[[vv]] <- make_transition_egg_adv(T_index,u=u,e_gen=g[j],stage1=i,
                                                    stage2=i+1,node=node_id)
      }
      T_index <- T_index + 1
      vv <- vv + 1
    }
  }

  # LARVAE TRANSITIONS

  # make larvae mortality transitions
  larvae_mort_tt <- vector("list",nL*nG)
  vv <- 1

  for(i in 1:nL){
    for(j in 1:nG){
      larvae_mort_tt[[vv]] <- make_transition_larvae_mort(T_index,u=u,l_gen=g[j],
                                                          stage=i,node=node_id)
      T_index <- T_index + 1
      vv <- vv + 1
    }
  }

  # make larvae advancement transitions
  larvae_adv_tt <- vector("list",nL*nG)
  vv <- 1

  for(i in 1:nL){
    for(j in 1:nG){
      if(i == nL){
        larvae_adv_tt[[vv]] <- make_transition_larvae_adv(T_index,u=u,l_gen=g[j],
                                                          stage1=i,stage2=NULL,
                                                          node=node_id)
      } else {
        larvae_adv_tt[[vv]] <- make_transition_larvae_adv(T_index,u=u,l_gen=g[j],
                                                          stage1=i,stage2=i+1,
                                                          node=node_id)
      }
      T_index <- T_index + 1
      vv <- vv + 1
    }
  }

  # PUPAE TRANSITIONS

  # make pupae mortality transitions
  pupae_mort_tt <- vector("list",nP*nG)
  vv <- 1

  for(i in 1:nP){
    for(j in 1:nG){
      pupae_mort_tt[[vv]] <- make_transition_pupae_mort(T_index,u=u,p_gen=g[j],
                                                        stage=i,node=node_id)
      T_index <- T_index + 1
      vv <- vv + 1
    }
  }

  # make pupae advancement transitions
  pupae_adv_tt <- vector("list",0)

  if(nP > 1){
    pupae_adv_tt <- vector("list",(nP-1)*nG)
    vv <- 1

    for(i in 1:(nP-1)){
      for(j in 1:nG){
        pupae_adv_tt[[vv]] <- make_transition_pupae_adv(T_index,u=u,p_gen=g[j],
                                                        stage1=i,stage2=i+1,
                                                        node=node_id)
        T_index <- T_index + 1
        vv <- vv + 1
      }
    }
  }

  # make pupae -> male emergence
  pupae_2male_tt <- vector("list",nG)
  vv <- 1

  for(j in 1:nG){
    pupae_2male_tt[[vv]] <- make_transition_pupae_emerge_m(T_index,u=u,p_gen=g[j],
                                                           nP=nP,node=node_id)
    T_index <- T_index + 1
    vv <- vv + 1
  }

  # UNMATED FEMALE TRANSITIONS

  # make pupae -> unmated female emergence
  pupae_2unmated_tt <- vector("list",nG)
  vv <- 1

  for(j in 1:nG){
    pupae_2unmated_tt[[vv]] <- make_transition_pupae_emerge_unmated(T_index,u=u,
                                                                    p_gen=g[j],
                                                                    nP=nP,node=node_id)
    T_index <- T_index + 1
    vv <- vv + 1
  }

  # make unmated female mortality
  unmated_mort_tt <- vector("list",nG)
  vv <- 1

  for(j in 1:nG){
    unmated_mort_tt[[vv]] <- make_transition_female_unmated_mort(T_index,u=u,f_gen=g[j],
                                                                 node=node_id)
    T_index <- T_index + 1
    vv <- vv + 1
  }

  # MALE TRANSITIONS

  # make male mortality
  male_mort_tt <- vector("list",nG)
  vv <- 1

  for(j in 1:nG){
    male_mort_tt[[vv]] <- make_transition_male_mort(T_index,u=u,m_gen=g[j],node=node_id)
    T_index <- T_index + 1
    vv <- vv + 1
  }


  # push T_index back up
  assign(x = "T_index", value = T_index, pos = parent.frame())

  # return as a list for further processing
  return(list("egg_mort" = egg_mort_tt,
              "egg_adv" = egg_adv_tt,
              "larvae_mort" = larvae_mort_tt,
              "larvae_adv" = larvae_adv_tt,
              "pupae_mort" = pupae_mort_tt,
              "pupae_adv" = pupae_adv_tt,
              "pupae_2male" = pupae_2male_tt,
              "pupae_2unmated" = pupae_2unmated_tt,
              "unmated_mort" = unmated_mort_tt,
              "male_mort" = male_mort_tt) )
} # end base mosy func

spn_T_mosy_lifecycle <- function(u,nE,nL,nP,cube,node_id = NULL,T_index){

  # genetic states
  g <- cube$genotypesID
  nG <- cube$genotypesN


  # get most transitions
  trans <- base_T_mosy(u=u,nE=nE,nL=nL,nP=nP,nG=nG,g=g,node_id=node_id,T_index=T_index)


  # empty list to put the transitions in (X_tt is the set of transitions in this subset of the total T)
  ovi_tt <- vector("list",sum((cube$tau * cube$ih) > 0))
  vv <- 1

  # OVIPOSITION
  ovi_dims <- dim(cube$ih)
  for(i in 1:ovi_dims[1]){
    for(j in 1:ovi_dims[2]){
      for(k in 1:ovi_dims[3]){
        # only make valid events (based on tau and ih)
        if(!fequal(cube$tau[i,j,k],0) && !fequal(cube$ih[i,j,k],0)){
          ovi_tt[[vv]] <- make_transition_ovi_epi(T_index,u=u,f_gen=g[i],m_gen=g[j],
                                                  o_gen=g[k],node=node_id)
          T_index <- T_index + 1
          vv <- vv + 1
        }
      }
    }
  }

  # make pupae -> female emergence
  pupae_2female_tt <- vector("list",nG^2)
  vv <- 1

  for(j_f in 1:nG){
    for(j_m in 1:nG){
      pupae_2female_tt[[vv]] <- make_transition_pupae_emerge_f(T_index,u=u,p_gen=g[j_f],
                                                               m_gen=g[j_m],nP=nP,
                                                               node=node_id)
      T_index <- T_index + 1
      vv <- vv + 1
    }
  }


  # FEMALE MORTALITY
  female_mort_tt <- vector("list",nG^2)
  vv <- 1

  for(j_f in 1:nG){
    for(j_m in 1:nG){
      female_mort_tt[[vv]] <- make_transition_female_mort_epi(T_index,u=u,f_gen=g[j_f],
                                                              m_gen=g[j_m],node=node_id)
      T_index <- T_index + 1
      vv <- vv + 1
    }
  }

  # UNMATED FEMALE MATING
  unmated_mate_tt <- vector("list",nG)
  vv <- 1

  for(j_f in 1:nG){
    for(j_m in 1:nG){
      unmated_mate_tt[[vv]] <- make_transition_female_unmated_mate(T_index,u=u,
                                                                   f_gen=g[j_f],
                                                                   m_gen=g[j_m],
                                                                   node=node_id,
                                                                   epi=FALSE)
      T_index <- T_index + 1
      vv <- vv + 1
    }
  }

  # the set of transitions
  tot_trans <- c(trans,
                 "oviposit" = list(ovi_tt),
                 "pupae_2female" = list(pupae_2female_tt),
                 "female_mort" = list(female_mort_tt),
                 "unmated_mate" = list(unmated_mate_tt)
               )

  # transitions (v)
  v <- unlist(x = lapply(X = tot_trans, FUN = lapply, '[[', 'label'), use.names = FALSE)

  # one long vector
  trans_vec <- unlist(x = tot_trans, recursive = FALSE, use.names = FALSE)


  # set T_index in parent environment
  #  ie, update the counter
  assign(x = "T_index", value = T_index, pos = parent.frame())

  # return the set of transitions and the vector (v)
  return(list("T" = trans_vec,
              "v" = v) )
}
