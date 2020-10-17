################################################################################
#
#   MGDrivE2: SPN structure for a single node (lifecycle only)
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   October 2019
#
################################################################################

################################################################################
# make the transitions (T) of the SPN
################################################################################

#' Make Transitions (T) For a Node (Mosquitoes only)
#'
#' This function makes the set of transitions (T) for a SPN. It is used alone if
#' our model is a single-node metapopulation for mosquito dynamics only; otherwise
#' it is used as part of other functions to make SPN models with larger state
#' spaces (metapopulation models, see \code{\link{spn_T_lifecycle_network}}).
#'
#' This function takes the places produced from \code{\link{spn_P_lifecycle_node}}
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
#' For examples of using this function, see:
#' \code{vignette("lifecycle-node", package = "MGDrivE2")}
#'
#' @param spn_P set of places produced by \code{\link{spn_P_lifecycle_node}}
#' @param params a named list of parameters (see details)
#' @param cube an inheritance cube from the \code{MGDrivE} package (e.g. \code{\link[MGDrivE]{cubeMendelian}})
#'
#' @return a list with two elements: \code{T} contains transitions packets as lists,
#' \code{v} is the character vector of transitions (T)
#'
#' @export
spn_T_lifecycle_node <- function(spn_P,params,cube){

  # set of places
  u <- spn_P$u

  # time index
  T_index <- 1

  # list of places
  retList <- spn_T_mosy_lifecycle(u = u,nE = params$nE,nL = params$nL,nP = params$nP,
                                  cube = cube,node_id = NULL,T_index = T_index)


  # check the set for errors
  labels <- unlist(x = lapply(X = retList$T, FUN = '[[', 'label'), use.names = FALSE)
  indices <- unlist(x = lapply(X = retList$T, FUN = '[[', 'vix'), use.names = FALSE)
  if(any(labels != retList$v[indices])){
    stop(paste0("error in set of transitions T and transition vector v at transition(s): ",
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
# OVIPOSITION
################################################################################

# transition event for oviposition
make_transition_ovi_epi <- function(T_index,u,f_gen,m_gen,o_gen,inf=NULL,node=NULL){

  # tokens required
  ftoken <- paste0(c("F",f_gen,m_gen,inf,node),collapse = "_")

  # tokens produced
  etoken <- paste0(c("E1",o_gen,node),collapse = "_")

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(ftoken,"->",etoken) # name of this t (corresponds to v)

  # requires a female token
  t$s <- match(x = ftoken,table = u)
  t$s_w <- 1

  # outputs a female token and an egg token
  t$o <- c(t$s,match(x = etoken,table = u))
  t$o_w <- c(1,1)

  # class of the transition
  t$class <- "oviposit"

  # return the transition
  return(t)
}

################################################################################
# EGGS
################################################################################

# transition event for egg advancement
# a NULL stage2 means they advance to larvae
make_transition_egg_adv <- function(T_index,u,e_gen,stage1,stage2=NULL,node=NULL){

  # tokens required
  etoken <- paste0(c(paste0("E",stage1),e_gen,node),collapse = "_")

  # tokens produced (otoken = output token, maybe nomenclature could be better)
  if(is.null(stage2)){
    otoken <- paste0(c("L1",e_gen,node),collapse = "_")
  } else {
    otoken <- paste0(c(paste0("E",stage2),e_gen,node),collapse = "_")
  }

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(etoken,"->",otoken) # name of this t (corresponds to v)

  # requires a egg token
  t$s <- match(x = etoken,table = u)
  t$s_w <- 1

  # outputs a new egg token of the next stage; or a larvae of stage 1
  t$o <- match(x = otoken,table = u)
  t$o_w <- 1

  # class of the transition
  t$class <- "egg_adv"

  # return the transition
  return(t)
}

# death event for eggs
make_transition_egg_mort <- function(T_index,u,e_gen,stage,node=NULL){

  # tokens required
  etoken <- paste0(c(paste0("E",stage),e_gen,node),collapse = "_")

  # no tokens are produced

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(etoken,"->D") # name of this t (corresponds to v)

  # requires a egg token
  t$s <- match(x = etoken,table = u)
  t$s_w <- 1

  # no tokens produced
  t$o <- NaN
  t$o_w <- NaN

  # class of the transition
  t$class <- "egg_mort"

  # return the transition
  return(t)
}

################################################################################
# LARVAE
################################################################################

# larvae advancement
# NULL stage2 indicates advancement to pupae 1
make_transition_larvae_adv <- function(T_index,u,l_gen,stage1,stage2=NULL,node=NULL){

  # tokens required
  ltoken <- paste0(c(paste0("L",stage1),l_gen,node),collapse = "_")

  # tokens produced (otoken = output token, maybe nomenclature could be better)
  if(is.null(stage2)){
    otoken <- paste0(c("P1",l_gen,node),collapse = "_")
  } else {
    otoken <- paste0(c(paste0("L",stage2),l_gen,node),collapse = "_")
  }

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(ltoken,"->",otoken) # name of this t (corresponds to v)

  # requires a egg token
  t$s <- match(x = ltoken,table = u)
  t$s_w <- 1

  # outputs a new egg token of the next stage; or a larvae of stage 1
  t$o <- match(x = otoken,table = u)
  t$o_w <- 1

  # class of the transition
  t$class <- "larvae_adv"

  # return the transition
  return(t)
}

# death event for larvae
make_transition_larvae_mort <- function(T_index,u,l_gen,stage,node=NULL){

  # tokens required
  ltoken <- paste0(c(paste0("L",stage),l_gen,node),collapse = "_")

  # no tokens are produced

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(ltoken,"->D") # name of this t (corresponds to v)

  # requires a larvae token
  t$s <- match(x = ltoken,table = u)
  t$s_w <- 1

  # no tokens produced
  t$o <- NaN
  t$o_w <- NaN

  # class of the transition
  t$class <- "larvae_mort"

  # return the transition
  return(t)
}


################################################################################
# PUPAE
################################################################################

# inter-pupal stages transitions (handle emergence "delicately" ... )
# so only from 1,...,nP-1
make_transition_pupae_adv <- function(T_index,u,p_gen,stage1,stage2,node=NULL){

  # tokens required
  ptoken <- paste0(c(paste0("P",stage1),p_gen,node),collapse = "_")

  # tokens produced (otoken = output token, maybe nomenclature could be better)
  otoken <- paste0(c(paste0("P",stage2),p_gen,node),collapse = "_")

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(ptoken,"->",otoken) # name of this t (corresponds to v)

  # requires a egg token
  t$s <- match(x = ptoken,table = u)
  t$s_w <- 1

  # outputs a new egg token of the next stage; or a pupae of stage 1
  t$o <- match(x = otoken,table = u)
  t$o_w <- 1

  # class of the transition
  t$class <- "pupae_adv"

  # return the transition
  return(t)
}

# death event for pupae
make_transition_pupae_mort <- function(T_index,u,p_gen,stage,node=NULL){

  # tokens required
  ptoken <- paste0(c(paste0("P",stage),p_gen,node),collapse = "_")

  # no tokens are produced

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(ptoken,"->D") # name of this t (corresponds to v)

  # requires a pupae token
  t$s <- match(x = ptoken,table = u)
  t$s_w <- 1

  # no tokens produced
  t$o <- NaN
  t$o_w <- NaN

  # class of the transition
  t$class <- "pupae_mort"

  # return the transition
  return(t)
}

# emergence to males event for pupae
make_transition_pupae_emerge_m <- function(T_index,u,p_gen,nP,node=NULL){

  # tokens required
  ptoken <- paste0(c(paste0("P",nP),p_gen,node),collapse = "_")

  # produces a male token
  mtoken <- paste0(c("M",p_gen,node),collapse = "_")

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(ptoken,"->",mtoken) # name of this t (corresponds to v)

  # requires a pupae token
  t$s <- match(x = ptoken,table = u)
  t$s_w <- 1

  # male token produced
  t$o <- match(x = mtoken,table = u)
  t$o_w <- 1

  # class of the transition
  t$class <- "pupae_2m"

  # return the transition
  return(t)
}

# emergence to females event for pupae
# p_gen: genotype of the pupae
# m_gen: genotype of the male mate
# nP: number of pupae compartments
make_transition_pupae_emerge_f <- function(T_index,u,p_gen,m_gen,nP,node=NULL){

  # tokens required
  ptoken <- paste0(c(paste0("P",nP),p_gen,node),collapse = "_")
  mtoken <- paste0(c("M",m_gen,node),collapse = "_")

  # produces a female token (with her genotype + the genotype of her mate)
  # also return the male token.
  ftoken <- paste0(c("F",p_gen,m_gen,node),collapse = "_")

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(ptoken,"->",ftoken) # name of this t (corresponds to v)

  # requires a pupae token and a male token
  t$s <- match(x = c(ptoken,mtoken),table = u)
  t$s_w <- c(1,1)

  # produces a female and a male token
  t$o <- match(x = c(ftoken,mtoken),table = u)
  t$o_w <- c(1,1)

  # class of the transition
  t$class <- "pupae_2f"

  # return the transition
  return(t)
}

# emergence to females (unmated) event for pupae
# p_gen: genotype of the pupae
# m_gen: genotype of the male mate
# nP: number of pupae compartments
make_transition_pupae_emerge_unmated <- function(T_index,u,p_gen,nP,node=NULL){

  # tokens required
  ptoken <- paste0(c(paste0("P",nP),p_gen,node),collapse = "_")

  # produces a female token (with her genotype + the genotype of her mate)
  # also return the male token.
  ftoken <- paste0(c("U",p_gen,node),collapse = "_")

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(ptoken,"->",ftoken) # name of this t (corresponds to v)

  # requires a pupae token and a male token
  t$s <- match(x = ptoken,table = u)
  t$s_w <- 1L

  # produces a female and a male token
  t$o <- match(x = ftoken,table = u)
  t$o_w <- 1L

  # class of the transition
  t$class <- "pupae_2unmated"

  # return the transition
  return(t)
}


################################################################################
# MALES
################################################################################

# kill a male token
make_transition_male_mort <- function(T_index,u,m_gen,node=NULL){

  # tokens required
  mtoken <- paste0(c("M",m_gen,node),collapse = "_")

  # no tokens are produced

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(mtoken,"->D") # name of this t (corresponds to v)

  # requires a male token
  t$s <- match(x = mtoken,table = u)
  t$s_w <- 1

  # no tokens produced
  t$o <- NaN
  t$o_w <- NaN

  # class of the transition
  t$class <- "male_mort"

  # return the transition
  return(t)
}


################################################################################
# FEMALES
################################################################################

# kill a female token
# inf: a character in ("S","E1",...,"En","I")
# where n is number of bins in EIP
make_transition_female_mort_epi <- function(T_index,u,f_gen,m_gen,inf=NULL,node=NULL){

  # tokens required
  ftoken <- paste0(c("F",f_gen,m_gen,inf,node),collapse = "_")

  # no tokens are produced

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(ftoken,"->D") # name of this t (corresponds to v)

  # requires a female token
  t$s <- match(x = ftoken,table = u)
  t$s_w <- 1

  # no tokens produced
  t$o <- NaN
  t$o_w <- NaN

  # class of the transition
  t$class <- "female_mort"

  # return the transition
  return(t)
}

# kill an unmated female token
make_transition_female_unmated_mort <- function(T_index,u,f_gen,node=NULL){

  # tokens required
  ftoken <- paste0(c("U",f_gen,node),collapse = "_")

  # no tokens are produced

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(ftoken,"->D") # name of this t (corresponds to v)

  # requires a female token
  t$s <- match(x = ftoken,table = u)
  t$s_w <- 1

  # no tokens produced
  t$o <- NaN
  t$o_w <- NaN

  # class of the transition
  t$class <- "female_unmated_mort"

  # return the transition
  return(t)
}

# mate an unmated female token
make_transition_female_unmated_mate <- function(T_index,u,f_gen,m_gen,node=NULL,epi=FALSE){

  # tokens required
  ftoken_u <- paste0(c("U",f_gen,node),collapse = "_")
  mtoken <- paste0(c("M",m_gen,node),collapse = "_")

  # produces a female token (with her genotype + the genotype of her mate)
  # also return the male token.
  if(epi){
    ftoken_m <- paste0(c("F",f_gen,m_gen,"S",node),collapse = "_")
  } else {
    ftoken_m <- paste0(c("F",f_gen,m_gen,node),collapse = "_")
  }

  # t: {index into v, label, input arcs/weights, output arcs/weights}
  t <- list()
  t$vix <- T_index # where we can find this t in v
  t$label <- paste0(ftoken_u,"->",ftoken_m) # name of this t (corresponds to v)

  # requires a pupae token and a male token
  t$s <- match(x = c(ftoken_u,mtoken),table = u)
  t$s_w <- c(1,1)

  # produces a female and a male token
  t$o <- match(x = c(ftoken_m,mtoken),table = u)
  t$o_w <- c(1,1)

  # class of the transition
  t$class <- "female_unmated_mate"

  # return the transition
  return(t)
}
