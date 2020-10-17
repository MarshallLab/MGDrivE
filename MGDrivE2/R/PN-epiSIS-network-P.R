################################################################################
#
#   MGDrivE2: SPN structure for a metapopulation network (SEI-SIS epi)
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   November 2019
#
################################################################################

################################################################################
#   stitch together places {P} for the metapop
################################################################################

#' Make Places (P) For a Network (SEI Mosquitoes - SIS Humans)
#'
#' This function makes the set of places (P) for a SPN model of a metapopulation
#' network for simulation of coupled SEI-SIS dynamics. It is the network version
#' of \code{\link{spn_P_epiSIS_node}}.
#'
#' The \code{params} argument supplies all of the ecological parameters necessary
#' to calculate equilibrium values. This function requires the \code{nE},
#' \code{nL}, \code{nP}, and \code{nEIP} parameters to be specified. For more details, see
#' \code{\link{equilibrium_SEI_SIS}}
#'
#' For examples of using this function, see:
#' \code{vignette("epi-network", package = "MGDrivE2")}
#'
#' @param node_list a character vector specifying what type of nodes to create;
#' (m = a node_id with only mosquitoes, h = a node_id with only humans, b = a node_id with both humans and mosquitoes)
#' @param params a named list of parameters (see details)
#' @param cube an inheritance cube from the \code{MGDrivE} package (e.g. \code{\link[MGDrivE]{cubeMendelian}})
#'
#' @return a list with two elements: \code{ix} contains labeled indices of the
#' places by life stage and node_id, \code{u} is the character vector of places (P)
#'
#' @export
spn_P_epiSIS_network <- function(node_list,params,cube){

  # checks
  nE <- params$nE
  nL <- params$nL
  nP <- params$nP
  nEIP <- params$nEIP
  stopifnot(node_list %in% c("b","m","h"))
  if(nE < 2 || nL < 2 || nP < 2 || nEIP < 2){
    warning(paste0("A shape parameter ('nE', 'nL', 'nP', 'nEIP') of 1 implies ",
                   "exponentially distributed dwell times in that compartment."))
  }

  # genetic information
  nG <- cube$genotypesN
  g <- cube$genotypesID

  # within node_id places
  P_meta <- vector(mode = "list",length = length(node_list))
  P_offset <- 0

  # loop through all nodes
  for(id in 1:length(node_list)){
    if(node_list[id] == "m"){
      P_meta[[id]] <- spn_P_mosy_epi(nE = nE,nL = nL,nP = nP,nEIP = nEIP,
                                         nG = nG, g = g, node_id = id,
                                         P_offset = P_offset)

    } else if(node_list[id] == "h"){
      P_meta[[id]] <- spn_P_humans_epiSIS(node_id = id,P_offset = P_offset)

    } else if(node_list[id] == "b"){
      P_meta[[id]] <- spn_P_both_epiSIS(nE = nE,nL = nL,nP = nP,nEIP = nEIP,
                                        nG = nG, g = g, node_id = id,
                                        P_offset = P_offset)

    } else {
      stop("warning: unrecognized character in 'node_list' object")
    }

  } # end node_id loop


  # return the stitched-together set of places P
  #  pull out the indexing (ix), store as one-depth list "ix"
  #  pull out places (u), store as vector "u"
  return(list("ix" = lapply(X = P_meta, FUN = '[[', 'ix'),
              "u" = unlist(lapply(X = P_meta, FUN = '[[', 'u')) )
         )
}


################################################################################
# NOTE: the below functions are not exported from the package for users
# as the package provides support only to generate Petri Nets by
# returning complete sets (P,T); if the below functions were exported
# it would be possible for a user to return an incomplete set of places (P)
# or transitions (T)
################################################################################

################################################################################
# make the places (P) of the SPN: mosquito
################################################################################

# node_id with only mosquitoes
spn_P_mosy_epi <- function(nE,nL,nP,nEIP,nG,g,node_id,P_offset){

  # setup place names
  eggs <- file.path("E",1:nE,"_",rep(g, each = nE),"_",node_id, fsep = "")

  larvae <- file.path("L",1:nL,"_",rep(g, each = nL),"_",node_id, fsep = "")

  pupae <- file.path("P",1:nP,"_",rep(g, each = nP),"_",node_id, fsep = "")

  females_unmated = file.path("U_",g,"_",node_id, fsep = "")

  stages <- c("S",paste0("E",1:nEIP),"I")
  females <- file.path("F", rep(x = rep(x = g, each = nG), times = nEIP+2),
                       rep(x = g, times = nG*(nEIP+2)),
                       rep(x = stages, each = nG^2), node_id, fsep = "_")

  males <- file.path("M",g,node_id,fsep = "_")


  # indices of states
  ix <- list()
  ix$egg <- matrix(data = seq_along(eggs) + P_offset,nrow = nE,byrow = FALSE,dimnames = list(1:nE,g))
  ix$larvae <- matrix(data = seq_along(larvae) + nG*nE + P_offset,nrow = nL,byrow = FALSE,dimnames = list(1:nL,g))
  ix$pupae <- matrix(data = seq_along(pupae) + nG*(nE + nL) + P_offset,nrow = nP,byrow = FALSE,dimnames = list(1:nP,g))
  ix$females_unmated <- setNames(object = seq_along(females_unmated) + nG*(nE + nL + nP) + P_offset, nm = g)

  # For the females, we have to permute because the array function doesn't have
  # something like "byrow" for matrices. This makes sure the "rows" are the female
  # genotype, cols are male mate, and z-axis is stage of incubation.
  # "resize=TRUE" is a cheat to keep dim names, doesn't cause issues in this instance
  ix$females <- aperm(a = array(data = seq_along(females) + nG*(nE + nL + nP + 1) + P_offset,dim = c(nG,nG,nEIP+2),
                                dimnames = list(g,g,stages)),
                      perm = c(2,1,3),resize = TRUE)

  ix$males <- setNames(object = seq_along(males) + nG*(nE+nL+nP+nG*(nEIP+2) + 1) + P_offset, nm = g)


  # places (u)
  u <- c(eggs,larvae,pupae,females_unmated,females,males)

  # set P_offset in parent environment
  #  ie, update the counter
  assign(x = "P_offset", value = P_offset + nG*(nE+nL+nP+nG*(nEIP+2)+2), pos = parent.frame())

  return(list("ix" = ix,
              "u" = u) )
}


################################################################################
# make the places (P) of the SPN: humans
################################################################################

# node_id with only humans
spn_P_humans_epiSIS <- function(node_id,P_offset){

  # human names
  hNames <- paste0(c("H_S","H_I"),"_",node_id)

  # indices of states
  ix <- list("humans" = setNames(object = c(1,2) + P_offset, nm = hNames))

  # places (u)
  u <- hNames

  # set P_offset in parent environment
  #  ie, update the counter
  assign(x = "P_offset", value = P_offset + 2, pos = parent.frame())

  return(list("ix" = ix,
              "u" = u) )
}

################################################################################
# make the places (P) of the SPN: mosquitoes and humans
################################################################################

spn_P_both_epiSIS <- function(nE,nL,nP,nEIP,nG,g,node_id,P_offset){

  # mosquitoes
  mList <- spn_P_mosy_epi(nE = nE,nL = nL,nP = nP,nEIP = nEIP,
                              nG = nG, g = g, node_id = node_id,
                              P_offset = P_offset)

  # humans
  hList <- spn_P_humans_epiSIS(node_id = node_id,P_offset = P_offset)


  # set P_offset in parent environment
  #  ie, update the counter
  assign(x = "P_offset", value = P_offset, pos = parent.frame())


  return(list("ix"=c(mList$ix,hList$ix),
              "u"=c(mList$u,hList$u)) )
}
