################################################################################
#
#   MGDrivE2: SPN structure for a network (lifecycle only)
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   October 2019
#
################################################################################

################################################################################
# make the places (P) of the SPN
################################################################################

#' Make Places (P) For a Network (Mosquitoes only)
#'
#' This function makes the set of places (P) for a SPN model of a metapopulation
#' network. It is the network version of \code{\link{spn_P_lifecycle_node}}.
#'
#' The \code{params} argument supplies all of the ecological parameters necessary
#' to calculate equilibrium values. This function requires the \code{nE},
#' \code{nL}, and \code{nP} parameters to be specified. For more details, see
#' \code{\link{equilibrium_lifeycle}}
#'
#' @param num_nodes number of nodes in the network
#' @param params a named list of parameters (see details)
#' @param cube an inheritance cube from the \code{MGDrivE} package (e.g. \code{\link[MGDrivE]{cubeMendelian}})
#'
#' @return a list with two elements: \code{ix} contains labeled indices of the
#' places by life stage and node_id, \code{u} is the character vector of places (P)
#'
#' @examples
#' \dontrun{
#'   # parameters, see vignette MGDrivE2: One Node Lifecycle Dynamics
#'   theta <- list(qE = 1/4, nE = 2, qL = 1/3, nL = 3, qP = 1/6, nP = 2,
#'                 muE = 0.05, muL = 0.15, muP = 0.05, muF = 0.09, muM = 0.09,
#'                 beta = 16, nu = 1/(4/24) )
#'
#'   # setup a 2-node network
#'   spn_P <- spn_P_lifecycle_network(num_nodes = 2, params = theta,
#'                                    cube = MGDrivE::cubeMendelian() )
#' }
#'
#' @export
spn_P_lifecycle_network <- function(num_nodes,params,cube){

  # checks
  nE <- params$nE
  nL <- params$nL
  nP <- params$nP
  if(nE < 2 || nL < 2 || nP < 2 ){
    warning(paste0("A shape parameter ('nE', 'nL', or 'nP') of 1 implies ",
                   "exponentially distributed dwell times in that compartment."))
  }

  # genetic information
  nG <- cube$genotypesN
  g <- cube$genotypesID

  # within node_id places
  P_meta <- vector(mode = "list",length = num_nodes)
  P_offset <- 0

  # loop through all nodes
  for(id in 1:num_nodes){
    # only 1 type of node in this network
    P_meta[[id]] <- spn_P_mosy_lifecycle(nE = nE,nL = nL,nP = nP,
                                             nG = nG, g = g, node_id = id,
                                             P_offset = P_offset)

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
#
# In this case, the function is almost an exact copy of the 1 node places function
# I needed a counter, and a way to handle node IDs
# I couldn't figure out how to append node IDs without breaking stuff :-(
################################################################################

spn_P_mosy_lifecycle <- function(nE,nL,nP,nG,g,node_id,P_offset){

  # setup place names
  eggs <- file.path("E",1:nE,"_",rep(g, each = nE),"_",node_id, fsep = "")

  larvae <- file.path("L",1:nL,"_",rep(g, each = nL),"_",node_id, fsep = "")

  pupae <- file.path("P",1:nP,"_",rep(g, each = nP),"_",node_id, fsep = "")

  females_unmated <- file.path("U_", g, "_", node_id, fsep = "")

  females <- file.path("F_",rep(g, each = nG),"_",g,"_",node_id, fsep = "")

  males <- file.path("M",g,node_id,fsep = "_")


  # indices of states
  ix <- list()
  ix$egg <- matrix(seq_along(eggs) + P_offset,nrow = nE,byrow = FALSE,dimnames = list(1:nE,g))
  ix$larvae <- matrix(data = seq_along(larvae) + nG*nE + P_offset,
                     nrow = nL,byrow = FALSE,dimnames = list(1:nL,g))
  ix$pupae <- matrix(data = seq_along(pupae) + nG*(nE + nL) + P_offset,
                    nrow = nP,byrow = FALSE,dimnames = list(1:nP,g))
  ix$females_unmated <- setNames(object = seq_along(females_unmated) + nG*(nE + nL + nP) + P_offset, nm = g)
  ix$females <- aperm(a = array(data = seq_along(females) + nG*(nE + nL + nP + 1) + P_offset,
                               dim = c(nG,nG,1), dimnames = list(g,g,"S")),
                     perm = c(2,1,3),resize = TRUE)
  ix$males <- setNames(object = seq_along(males) + nG*(nE+nL+nP+nG+1) + P_offset, nm = g)


  # places (u)
  u <- c(eggs,larvae,pupae,females_unmated,females,males)

  # set P_offset in parent environment
  #  ie, update the counter
  assign(x = "P_offset", value = P_offset + nG*(nE+nL+nP+nG+2), pos = parent.frame())

  return(list("ix" = ix,
              "u" = u) )
}
