# --------------------------------------------------------------------------------
#
#   MGDrivE2: batch mgiration
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   November 2020
#
# --------------------------------------------------------------------------------

#' Sample Batch Migration Events
#'
#' Sample batch migration events for simulation given rates of occurance and probability of destination for each patch.
#' Batch migration can be simulated for the aquatic life stages (eggs, larvae, pupae), adult females, and/or adult males.
#' To simulate batch migration, each life stage needs all 3 of its arguments specified. If any arguments are left
#' unspecified (\code{NULL}), batch migration for that life stage will not be sampled.
#' The output of this function should be passed to \code{\link[MGDrivE2]{sim_trajectory_R}} or \code{\link[MGDrivE2]{sim_trajectory_CSV}}
#' as the argument \code{batch}.
#' Calls the internal function \code{\link[MGDrivE2]{batch_migration_stage}}.
#'
#' @param SPN_P places of the SPN
#' @param tmax maximum time of the simulation
#' @param ELPrates rate at which aquatic stage batch migration occurs for each node (nodes without mosquitoes should be set to \code{NaN} or \code{NA})
#' @param ELPmove movement matrix for destinations of aquatic stage batch migration events (diagonal will be set to zero and off-diagonal elements normalized)
#' @param ELPprob probability for each individual to be chosen for aquatic stage batch migration events (must be same length as \code{ELPrates})
#' @param Frates rate at which adult female batch migration occurs for each node (nodes without mosquitoes should be set to \code{NaN} or \code{NA})
#' @param Fmove movement matrix for destinations of adult female batch migration events (diagonal will be set to zero and off-diagonal elements normalized)
#' @param Fprob probability for each individual to be chosen for adult female batch migration events (must be same length as \code{Frates})
#' @param Mrates rate at which adult male batch migration occurs for each node (nodes without mosquitoes should be set to \code{NaN} or \code{NA})
#' @param Mmove movement matrix for destinations of adult male batch migration events (diagonal will be set to zero and off-diagonal elements normalized)
#' @param Mprob probability for each individual to be chosen for adult male batch migration events (must be same length as \code{Mrates})
#' @param stage either \code{NULL} or "E", "L", or "P". If not \code{NULL} and migration for aquatic stages is specified by \code{ELPrates}, only the aquatic stage specified here will move
#'
#' @return vector of lists describing all batch migration events, segmented by life stage.
#' @importFrom stats runif
#'
#' @export
batch_migration <- function(
  SPN_P, tmax,
  ELPrates = NULL, ELPmove = NULL, ELPprob = NULL,
  Frates = NULL, Fmove = NULL, Fprob = NULL,
  Mrates = NULL, Mmove = NULL, Mprob = NULL,
  stage = NULL
){

  if(length(SPN_P$ix)<2){
    stop("cannot set up batch migration for a 1-node system")
  }

  # only set up ELP events if all of the parameters are not NULL
  if( all(!is.null(ELPrates), !is.null(ELPmove), !is.null(ELPprob)) ){

    if(!all(dim(ELPmove) == length(SPN_P$ix))){
      stop("dimension of 'ELPmove' must be the same as the number of nodes in the Petri net")
    }
    if(length(ELPrates) != length(ELPprob)){
      stop("length of 'ELPrates' must be same as 'ELPprob'")
    }
    diag(ELPmove) <- 0

    if(!is.null(stage)){
      stopifnot(stage %in% c("E","L","P"))
    }

    aqua_stage <- ifelse(test = is.null(stage),yes = "ELP",no = stage)
    ELP_events <- batch_migration_stage(SPN_P = SPN_P, rates = ELPrates, move = ELPmove, prob = ELPprob, stage = aqua_stage, tmax = tmax)

  } else {

    ELP_events <- NULL

  }

  # adult female movement
  if( all(!is.null(Frates), !is.null(Fmove), !is.null(Fprob))){

    if(!all(dim(Fmove) == length(SPN_P$ix))){
      stop("dimension of 'Fmove' must be the same as the number of nodes in the Petri net")
    }
    if(length(Frates) != length(Fprob)){
      stop("length of 'Frates' must be same as 'Fprob'")
    }
    diag(Fmove) <- 0

    F_events <- batch_migration_stage(SPN_P = SPN_P, rates = Frates, move = Fmove, prob = Fprob, stage = "F", tmax = tmax)

  } else {

    F_events <- NULL

  }

  # adult male movement
  if( all(!is.null(Mrates), !is.null(Mmove), !is.null(Mprob))){

    if(!all(dim(Mmove) == length(SPN_P$ix))){
      stop("dimension of 'Mmove' must be the same as the number of nodes in the Petri net")
    }
    if(length(Mrates) != length(Mprob)){
      stop("length of 'Mrates' must be same as 'Mprob'")
    }
    diag(Mmove) <- 0

    M_events <- batch_migration_stage(SPN_P = SPN_P, rates = Mrates, move = Mmove, prob = Mprob, stage = "M", tmax = tmax)

  } else {

    M_events <- NULL

  }

  batch_list <- c(ELP_events, F_events, M_events)

  # sort by time of event
  batch_list <- batch_list[order(sapply(batch_list,function(x){x$time}))]
  return(batch_list)
}


#' Internal function to sample and set up data structure for batch migration
#'
#' @param SPN_P a set of Petri net places
#' @param rates a vector of rates for each node
#' @param move a movement matrix (where do the batches go?)
#' @param prob the probability vector for each individual moving in each batch
#' @param stage the life stage (one of 'ELP', 'F', 'M')
#' @param tmax maximum simulation time
#' @keywords internal
batch_migration_stage <- function(SPN_P, rates, move, prob, stage, tmax){

  if(!(stage %in% c("ELP","F","M", "E", "L", "P"))){
    stop("called from 'batch_migration_stage', argument 'stage' must be one of 'ELP', 'E', 'L', 'P', 'F', 'M'")
  }

  rates[which(!is.finite(rates))] <- 0
  events_num <- rpois(n = length(rates), lambda = rates * tmax)
  events_times <- lapply(X = events_num, FUN = function(nn){
    if(nn > 0){
      # sample from times points that are guaranteed to be hit by the sampler
      return(unique(round(stats::runif(n = nn, min = 1, max = tmax-1))))
    } else {
      return(NULL)
    }
  })

  events_tofrom <- mapply(FUN = function(id, tt){
    if(is.null(tt)){
      return(NULL)
    } else {
      # where did they go?
      dest <- sample.int(n = nrow(move),size = length(tt),replace = TRUE,prob = move[id, ])
    }
  },
    id = 1:nrow(move), tt = events_times, SIMPLIFY = FALSE,USE.NAMES = FALSE
  )

  events_idx <- mapply(FUN = function(id, event){
      if(is.null(event)){
        return(NULL)
      } else {
        out <- list("from"=NULL,"to"=NULL)
        if(stage=="ELP"){
          out$from <- c(SPN_P$ix[[id]]$egg,SPN_P$ix[[id]]$larvae,SPN_P$ix[[id]]$pupae)
        }
        if(stage == "E"){
          out$from <- c(SPN_P$ix[[id]]$egg)
        }
        if(stage == "L"){
          out$from <- c(SPN_P$ix[[id]]$larvae)
        }
        if(stage == "P"){
          out$from <- c(SPN_P$ix[[id]]$pupae)
        }
        if(stage == "F"){
          out$from <- c(SPN_P$ix[[id]]$females_unmated,as.vector(SPN_P$ix[[id]]$females))
        }
        if(stage == "M"){
          out$from <- as.vector(SPN_P$ix[[id]]$males)
        }
        out$from <- unname(out$from)
        if(is.null(out$from)){stop("origin node has null indices, perhaps the rate vector is misspecified?")}
        out$to <- matrix(data=0,ncol=length(event),nrow=length(out$from))
        # vector 'event' tells me destinations for these movements
        for(i in 1:length(event)){
          if(stage=="ELP"){
            out$to[,i] <- c(SPN_P$ix[[event[i]]]$egg,SPN_P$ix[[event[i]]]$larvae,SPN_P$ix[[event[i]]]$pupae)
          }
          if(stage == "E"){
            out$to[,i] <- c(SPN_P$ix[[event[i]]]$egg)
          }
          if(stage == "L"){
            out$to[,i] <- c(SPN_P$ix[[event[i]]]$larvae)
          }
          if(stage == "P"){
            out$to[,i] <- c(SPN_P$ix[[event[i]]]$pupae)
          }
          if(stage == "F"){
            out$to[,i] <- c(SPN_P$ix[[event[i]]]$females_unmated,as.vector(SPN_P$ix[[event[i]]]$females))
          }
          if(stage == "M"){
            out$to[,i] <- as.vector(SPN_P$ix[[event[i]]]$males)
          }
          if(is.null(out$to[,i])){stop("destination node has null indices, perhaps a movement matrix is misspecified?")}
        }
        return(out)
      }
    }, id = 1:length(events_tofrom), event = events_tofrom,
    SIMPLIFY = FALSE, USE.NAMES = FALSE
  )

  batch_events <- vector(mode="list",length=sum(sapply(events_times,length)))
  i <- 1

  for(j in 1:length(events_times)){
      if(is.null(events_times[[j]])){
          next
      } else {
          for(k in 1:length(events_times[[j]])){
              batch_events[[i]]$time <- events_times[[j]][k]
              batch_events[[i]]$from <- events_idx[[j]]$from
              batch_events[[i]]$to <- as.vector(events_idx[[j]]$to[,k])
              batch_events[[i]]$prob <- prob[j]
              i <- i + 1
          }
      }
  }

  return(batch_events)
}
