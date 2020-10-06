################################################################################
#
#   MGDrivE2: SPN structure for a metapopulation network (SEI-SEIR epi)
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   April 2020
#
################################################################################

################################################################################
#   stitch together places {P} for the metapop
################################################################################

#' Make Places (P) For a Network (SEI Mosquitoes - SEIR Humans)
#'
#' This function makes the set of places (P) for a SPN model of a metapopulation
#' network for simulation of coupled SEI-SEIR dynamics. It is the network version
#' of \code{\link{spn_P_epiSEIR_node}}.
#'
#' The \code{params} argument supplies all of the ecological parameters necessary
#' to calculate equilibrium values. This function requires the \code{nE},
#' \code{nL}, \code{nP}, and \code{nEIP} parameters to be specified. For more details, see
#' \code{\link{equilibrium_SEI_SEIR}}
#'
#' @param node_list a character vector specifying what type of nodes to create;
#' (m = a node with only mosquitoes, h = a node with only humans, b = a node with both humans and mosquitoes)
#' @param params a named list of parameters (see details)
#' @param cube an inheritance cube from the \code{MGDrivE} package (e.g. \code{\link[MGDrivE]{cubeMendelian}})
#'
#' @return a list with two elements: \code{ix} contains labeled indices of the places
#' by life stage and node, \code{u} is the character vector of places (P)
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
#'   # setup a 3-node network
#'   node_list <- c("m","b", "h")
#'
#'   spn_P <- spn_P_epiSEIR_network(node_list = node_list, params = theta,
#'                                  cube = MGDrivE::cubeMendelian() )
#' }
#'
#' @export
spn_P_epiSEIR_network <- function(node_list,params,cube){

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

  # within node places
  P_meta <- vector(mode = "list",length = length(node_list))
  P_offset <- 0

  # loop through all nodes
  for(id in 1:length(node_list)){
    if(node_list[id] == "m"){
      P_meta[[id]] <- spn_P_mosy_epi(nE = nE,nL = nL,nP = nP,nEIP = nEIP,
                                         nG = nG, g = g, node_id = id,
                                         P_offset = P_offset)

    } else if(node_list[id] == "h"){
      P_meta[[id]] <- spn_P_humans_epiSEIR(node_id = id,P_offset = P_offset)

    } else if(node_list[id] == "b"){
      P_meta[[id]] <- spn_P_both_epiSEIR(nE = nE,nL = nL,nP = nP,nEIP = nEIP,
                                        nG = nG, g = g, node_id = id,
                                        P_offset = P_offset)

    } else {
      stop("warning: unrecognized character in 'node_list' object")
    }

  } # end node loop


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
# make the places (P) of the SPN: humans
################################################################################

# node with only humans
spn_P_humans_epiSEIR <- function(node_id,P_offset){

  # human places
  hNames <- paste0(c("H_S","H_E","H_I","H_R"),"_",node_id)

  # indices of states
  ix <- list("humans"= setNames(object = c(1,2,3,4) + P_offset,nm = hNames))

  # places (u)
  u <- hNames

  # set P_offset in parent environment
  #  ie, update the counter
  assign(x = "P_offset", value = P_offset + 4, pos = parent.frame())

  return(list("ix" = ix,
              "u" = u) )
}

################################################################################
# make the places (P) of the SPN: mosquitoes and humans
################################################################################

spn_P_both_epiSEIR <- function(nE,nL,nP,nEIP,nG,g,node_id,P_offset){

  # mosquitoes
  mList <- spn_P_mosy_epi(nE = nE,nL = nL,nP = nP,nEIP = nEIP,
                              nG = nG, g = g, node_id = node_id,
                              P_offset = P_offset)

  # humans
  hList <- spn_P_humans_epiSEIR(node_id = node_id,P_offset = P_offset)


  # set P_offset in parent environment
  #  ie, update the counter
  assign(x = "P_offset", value = P_offset, pos = parent.frame())

  return(list("ix"=c(mList$ix,hList$ix),
              "u"=c(mList$u,hList$u)) )
}
