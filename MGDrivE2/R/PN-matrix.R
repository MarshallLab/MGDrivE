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
#'
spn_Pre <- function(spn_P, spn_T){

  u <- spn_P$u # dimension of the places
  v <- spn_T$v # dimension of the transitions

  ix <- lapply(X = spn_T$T, FUN = "[[", "vix")
  s <- lapply(X = spn_T$T, FUN = "[[", "s")
  s_w <- lapply(X = spn_T$T, FUN = "[[", "s_w")

  # for transitions with multiple input arcs
  s_len <- vapply(X = s, FUN = length, FUN.VALUE = integer(1))

  # replicate the elements of the ix (transition index) that number of times
  ix[which(s_len > 1)] <- mapply(FUN = function(ixv, num){
    rep(x = ixv, times = num)
  }, ixv = ix[which(s_len > 1)], num = s_len[which(s_len > 1)],
  SIMPLIFY = FALSE, USE.NAMES = FALSE
  )

  # ignore transitions that are always on (if any)
  null_idx <- vapply(X = spn_T$T, FUN = function(v){
    all(is.nan(v$s)) | all(is.nan(v$s_w))
  }, FUN.VALUE = logical(1))

  ix <- ix[!null_idx]
  s <- s[!null_idx]
  s_w <- s_w[!null_idx]

  ix <- unlist(ix)
  s <- unlist(s)
  s_w <- unlist(s_w)

  Pre <- sparseMatrix(
    i = as.integer(ix),
    j = as.integer(s),
    x = as.integer(s_w),
    dims = c(length(v), length(u)),
    dimnames = list(v, u)
  )

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
#' @importFrom Matrix sparseMatrix
#'
#' @return a matrix of type \code{\link[Matrix]{dgCMatrix-class}}
#'
spn_Post <- function(spn_P, spn_T){

  u <- spn_P$u # dimension of the places
  v <- spn_T$v # dimension of the transitions

  ix <- lapply(X = spn_T$T, FUN = "[[", "vix")
  o <- lapply(X = spn_T$T, FUN = "[[", "o")
  o_w <- lapply(X = spn_T$T, FUN = "[[", "o_w")

  # for transitions with multiple output arcs
  o_len <- vapply(X = o, FUN = length, FUN.VALUE = integer(1))

  # replicate the elements of the ix (transition index) that number of times
  ix[which(o_len > 1)] <- mapply(FUN = function(ixv, num){
    rep(x = ixv, times = num)
  }, ixv = ix[which(o_len > 1)], num = o_len[which(o_len > 1)],
  SIMPLIFY = FALSE, USE.NAMES = FALSE
  )

  # ignore transitions that do not produce tokens
  null_idx <- vapply(X = spn_T$T, FUN = function(v){
    all(is.nan(v$o)) | all(is.nan(v$o_w))
  }, FUN.VALUE = logical(1))

  ix <- ix[!null_idx]
  o <- o[!null_idx]
  o_w <- o_w[!null_idx]

  ix <- unlist(ix)
  o <- unlist(o)
  o_w <- unlist(o_w)

  Post <- sparseMatrix(
    i = as.integer(ix),
    j = as.integer(o),
    x = as.integer(o_w),
    dims = c(length(v), length(u)),
    dimnames = list(v, u)
  )

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
#' @return stoichiometry matrix representing the net effect of a transition in the SPN state model.
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
