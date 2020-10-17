################################################################################
#
#   MGDrivE2: generate sparse matrices from PN
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   October 2019
#
################################################################################

################################################################################
# make the Pre matrix (v by u)
################################################################################

#' Make Pre Matrix For a Petri Net
#'
#' Generate the Pre (|v| by |u|) matrix for the SPN. This gives the edges from P
#' to T (input arcs) in the bipartite network.
#'
#' The places (\code{spn_P}) object is generated from one of the following:
#' \code{\link{spn_P_lifecycle_node}}, \code{\link{spn_P_lifecycle_network}},
#' \code{\link{spn_P_epiSIS_node}}, \code{\link{spn_P_epiSIS_network}},
#' \code{\link{spn_P_epiSEIR_node}}, or \code{\link{spn_P_epiSEIR_network}}.
#'
#' The set of transitions (\code{spn_T}) is generated from one of the following:
#' \code{\link{spn_T_lifecycle_node}}, \code{\link{spn_T_lifecycle_network}},
#' \code{\link{spn_T_epiSIS_node}}, \code{\link{spn_T_epiSIS_network}},
#' \code{\link{spn_T_epiSEIR_node}}, \code{\link{spn_T_epiSEIR_network}}.
#'
#'
#' @param spn_P set of places (P) (see details)
#' @param spn_T set of transitions (T) (see details)
#'
#' @return a matrix of type \code{\link[Matrix]{dgCMatrix-class}}
#'
#' @importFrom Matrix sparseMatrix
#'
spn_Pre <- function(spn_P,spn_T){

  u <- spn_P$u # dimension of the places
  v <- spn_T$v # dimension of the transitions

  Pre <- sparseMatrix(i = {},j = {},x = c(0L),dims=c(length(v),length(u)),dimnames=list(v,u))

  # fill in the Pre matrix
  for(i in 1:length(spn_T$T)){

    # index into v (rows)
    ix <- spn_T$T[[i]]$vix

    # input arcs and weights
    s <- spn_T$T[[i]]$s
    s_w <- spn_T$T[[i]]$s_w

    # skip null stuff (in this case, transitions that are always on)
    if(all(is.nan(s)) || all(is.nan(s_w))){
      next
    }

    # input arcs go into the matrix
    Pre[ix,s] <- as.integer(s_w)
  }


  return(Pre)
}


################################################################################
# make the Post matrix (v by u)
################################################################################

#' Make Post Matrix For a Petri Net
#'
#' Generate the Post (|v| by |u|) matrix for the SPN. This gives the edges from
#' T to P (output arcs) in the bipartite network.
#'
#' The places (\code{spn_P}) object is generated from one of the following:
#' \code{\link{spn_P_lifecycle_node}}, \code{\link{spn_P_lifecycle_network}},
#' \code{\link{spn_P_epiSIS_node}}, \code{\link{spn_P_epiSIS_network}},
#' \code{\link{spn_P_epiSEIR_node}}, or \code{\link{spn_P_epiSEIR_network}}.
#'
#' The set of transitions (\code{spn_T}) is generated from one of the following:
#' \code{\link{spn_T_lifecycle_node}}, \code{\link{spn_T_lifecycle_network}},
#' \code{\link{spn_T_epiSIS_node}}, \code{\link{spn_T_epiSIS_network}},
#' \code{\link{spn_T_epiSEIR_node}}, \code{\link{spn_T_epiSEIR_network}}.
#'
#'
#' @param spn_P set of places (P) (see details)
#' @param spn_T set of transitions (T) (see details)
#'
#' @return a matrix of type \code{\link[Matrix]{dgCMatrix-class}}
#'
spn_Post <- function(spn_P,spn_T){

  u <- spn_P$u # dimension of the places
  v <- spn_T$v # dimension of the transitions

  Post <- sparseMatrix(i = {},j = {},x = c(0L),dims=c(length(v),length(u)),dimnames=list(v,u))

  # fill in the Post matrix
  for(i in 1:length(spn_T$T)){

    # index into v (rows)
    ix <- spn_T$T[[i]]$vix

    # output arcs and weights
    o <- spn_T$T[[i]]$o
    o_w <- spn_T$T[[i]]$o_w

    # skip null stuff (in this case, transitions that only consume tokens)
    if(all(is.nan(o)) || all(is.nan(o_w))){
      next
    }

    # output arcs go into the matrix
    Post[ix,o] <- as.integer(o_w)
  }

  return(Post)
}


################################################################################
# make the other matrices that are useful
# A: reaction matrix; v by u
# S: stoichiometry matrix; u by v
################################################################################

#' Make stoichiometry Matrix For a Petri Net
#'
#' Generate the stoichiometry (|u| by |v|) matrix for the SPN.
#' Each column gives the net effect of that transition firing upon the state
#' space of the model. Internally, this creates a Pre (\code{\link{spn_Pre}}) and
#' Post (\code{\link{spn_Post}}) matrix, and then calculates the final stoichiometry.
#'
#' The places (\code{spn_P}) object is generated from one of the following:
#' \code{\link{spn_P_lifecycle_node}}, \code{\link{spn_P_lifecycle_network}},
#' \code{\link{spn_P_epiSIS_node}}, \code{\link{spn_P_epiSIS_network}},
#' \code{\link{spn_P_epiSEIR_node}}, or \code{\link{spn_P_epiSEIR_network}}.
#'
#' The set of transitions (\code{spn_T}) is generated from one of the following:
#' \code{\link{spn_T_lifecycle_node}}, \code{\link{spn_T_lifecycle_network}},
#' \code{\link{spn_T_epiSIS_node}}, \code{\link{spn_T_epiSIS_network}},
#' \code{\link{spn_T_epiSEIR_node}}, \code{\link{spn_T_epiSEIR_network}}.
#'
#'
#' @param spn_P set of places (P) (see details)
#' @param spn_T set of transitions (T) (see details)
#'
#' @importFrom Matrix drop0 t
#'
#' @export
spn_S <- function(spn_P,spn_T){

  # create pre matrix
  Pre <- spn_Pre(spn_P = spn_P,spn_T = spn_T)

  # create post matrix
  Post <- spn_Post(spn_P = spn_P,spn_T = spn_T)

  # calculate difference stoichiometry
  # A matrix
  A <- Post - Pre

  # A has 0 being stored; get rid of it
  A_new <- drop0(x = A)

  # return (u X v) matrix
  return(t(A_new))
}
