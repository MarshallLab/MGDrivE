################################################################################
#
#   MGDrivE2: SPN structure for a single node (lifecycle only)
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   October 2019
#
################################################################################

################################################################################
# make the places (P) of the SPN
################################################################################

#' Make Places (P) For a Node (Mosquitoes only)
#'
#' This function makes the set of places (P) for a SPN. It is used alone if our model
#' is a single-node metapopulation for mosquito dynamics only; otherwise it is used
#' as part of other functions to make SPN models with larger state spaces
#' (metapopulation models, see \code{\link{spn_P_lifecycle_network}}).
#'
#' The \code{params} argument supplies all of the ecological parameters necessary
#' to calculate equilibrium values. This function requires the \code{nE},
#' \code{nL}, and \code{nP} parameters to be specified. For more details, see
#' \code{\link{equilibrium_lifeycle}}
#'
#' @param params a named list of parameters (see details)
#' @param cube an inheritance cube from the \code{MGDrivE} package (e.g. \code{\link[MGDrivE]{cubeMendelian}})
#'
#' @return a list with two elements: \code{ix} contains labeled indices of the places
#' by life stage, \code{u} is the character vector of places (P)
#'
#' @importFrom stats setNames
#'
#' @examples
#' \dontrun{
#'   # parameters, see vignette MGDrivE2: One Node Lifecycle Dynamics
#'   theta <- list(qE = 1/4, nE = 2, qL = 1/3, nL = 3, qP = 1/6, nP = 2,
#'                 muE = 0.05, muL = 0.15, muP = 0.05, muF = 0.09, muM = 0.09,
#'                 beta = 16, nu = 1/(4/24) )
#'
#'   spn_P <- spn_P_lifecycle_node(params = theta, cube = MGDrivE::cubeMendelian() )
#' }
#'
#' @export
spn_P_lifecycle_node <- function(params,cube){

  # checks
  nE <- params$nE
  nL <- params$nL
  nP <- params$nP
  if(nE < 2 || nL < 2 || nP < 2){
    warning(paste0("A shape parameter ('nE', 'nL', or 'nP') of 1 implies ",
                   "exponentially distributed dwell times in that compartment."))
  }

  # genetic information
  nG <- cube$genotypesN
  g <- cube$genotypesID

  # setup place names
  eggs <- file.path("E",1:nE,"_",rep(g, each = nE), fsep = "")

  larvae <- file.path("L",1:nL,"_",rep(g, each = nL), fsep = "")

  pupae <- file.path("P",1:nP,"_",rep(g, each = nP), fsep = "")

  females_unmated <- file.path("U",g, fsep = "_")

  females <- file.path("F_",rep(g, each = nG),"_",g, fsep = "")

  males <- file.path("M",g,fsep = "_")

  # indices of states
  ix <- list()
  ix$egg <- matrix(seq_along(eggs),nrow = nE,byrow = FALSE,dimnames = list(1:nE,g))
  ix$larvae <- matrix(data = seq_along(larvae) + nG*nE,nrow = nL,byrow = FALSE,dimnames = list(1:nL,g))
  ix$pupae <- matrix(data = seq_along(pupae) + nG*(nE + nL),nrow = nP,byrow = FALSE,dimnames = list(1:nP,g))
  ix$females_unmated <- setNames(object = seq_along(females_unmated) + nG*(nE + nL + nP), nm = g)
  ix$females <- aperm(a = array(data = seq_along(females) + nG*(nE + nL + nP + 1),dim = c(nG,nG,1),
                               dimnames = list(g,g,"S")),
                     perm = c(2,1,3),resize = TRUE)
  ix$males <- setNames(object = seq_along(males) + nG*(nE+nL+nP+nG+1), nm = g)


  # places (u)
  u <- c(eggs,larvae,pupae,females_unmated,females,males)

  # return list of places
  #  make ix a list to match network version
  return(list("ix" = list(ix),
              "u" = u) )
}
