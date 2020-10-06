# --------------------------------------------------------------------------------
#
#   MGDrivE2: create tracking matrices
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   October 2020
#
# --------------------------------------------------------------------------------

#' Make tracking matrix for human infection events
#'
#' Create a matrix object for tracking incidence in human population
#' to be passed to either \code{\link{sim_trajectory_CSV}} or \code{\link{sim_trajectory_R}}.
#'
#' The returned matrix can be passed to the \code{Sout} argument of \code{\link{sim_trajectory_CSV}} or \code{\link{sim_trajectory_R}}.
#'
#'
#' @param spn_T set of transitions
#' @param S stoichiometry matrix
#'
#' @return a \code{\link[Matrix]{sparseMatrix}} object
#'
#' @export
track_hinf <- function(spn_T,S){

  # which events correspnd to infection?
  inf_events <- which(sapply(spn_T$T,function(x){x$class}) == "H_infection")

  # what nodes are those happening in?
  inf_nodes <- as.integer(
    sapply(
      X = strsplit(x = sapply(X = spn_T$T[inf_events],FUN = function(x){x$label}),split = "->"),
      FUN = function(xx){
        strsplit(x = xx[1],split = "_")[[1]][3]
      }
    )
  )
  # if only one node, it will return all NA
  if(all(is.na(inf_nodes))){
    inf_nodes <- rep(1,length(inf_nodes))
  }

  # unique nodes with infection
  u_nodes <- unique(inf_nodes)

  # rle for rows of Sout
  node_rle <- rle(inf_nodes)
  node_rle <- rep(x = seq_along(node_rle$values),times = node_rle$lengths)

  # matrix object to track output
  Sout <- Matrix::sparseMatrix(
    i = node_rle,j = inf_events,
    x = 1, dimnames = list(paste0("inc_",u_nodes),colnames(S)),dims = c(length(u_nodes),ncol(S))
  )

  return(Sout)
}
