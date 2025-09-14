################################################################################
#
#   MGDrivE2: Auxiliary & internal functions
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   October 2019
#
################################################################################

################################################################################
# Internal Functions
################################################################################

# float compare
fequal <- function(x,y,tol=1.490116e-08){
  # tol = sqrt(.Machine$double.eps)
  return(abs(x-y) <= tol)
}

# used to check exact vs approximate hazards
#  just throws a warning if it could be an issue
check_approx <- function(tol, exact){
  if((tol > 1e-6) && (!exact)){
    warning(paste0("hazard function tolerance ",tol,
                   " is large.\n\tconsider tolerance < 1e-6 for sufficient accuracy"))
  }
}

# something that will throw errors if you pass it anything except a well defined float
check_double <- function(dbl){
  return(is.null(dbl) || any(is.na(dbl)) || any(is.nan(dbl)) || !length(dbl))
}

# safe normalization
normalize <- function(myVec,tol = 1.490116e-08){
  # tol = sqrt(.Machine$double.eps)

  # get sum
  vecSum <- sum(myVec)

  # if sum is too small, return identical sized vector of zero
  #  else, return normalized vector
  if(vecSum <= tol){
    return(numeric(length = length(myVec)))
  } else {
    return(myVec / vecSum)
  }
}

# check density dependence choice and parameterization
#  log_dd: TRUE for logistic, FALSE for Lotka-Volterra
#  params: parameters for hazards
# Used directly in: spn_hazards_lifecycle_node,spn_hazards_lifecycle_network,spn_hazards_epiSIS_node,
check_dd <- function(log_dd, params){

  # check type first
  if(log_dd) {
    # make sure K is specified
    if(is.null(params$K)){
      stop("if using logistic (carrying capacity) based density-dependent",
           "larval mortality, please specify parameter 'K' in params")
    }
    # make sure it is specified properly
    if(check_double(params$K)){
      stop("K is improperly specified")
    }
  } else {
    if(is.null(params$gamma)){
      stop("if using Lotka-Volterra based density-dependent larval mortality,",
           "please specify parameter 'gamma' in params")
    }
    # make sure it is specified properly
    if(check_double(params$gamma)){
      stop("gamma is improperly specified")
    }
  }

}

# check parameters for correctness, throw error if not
#  used directly in spn_hazards
#  used indirectly below
check_params_life <- function(params){
  #  pull generic ones from the params
  #  check if doubles
  #  throw error if false
  boolVec = sapply(X = params[c("qE","nE","qL","nL","qP","nP","muE",
                                 "muL","muP","muF","muM","beta","nu")], FUN = check_double)
  if(all(boolVec)){
    stop(paste0(names(boolVec)[boolVec], collapse = ", "), " are improperly specified")
  }

}

#  used directly in spn_hazards
#  used indirectly below
check_params_sis <- function(params){

  # check basic life params
  check_params_life(params = params)

  # check new things for epi stuff
  boolVec = sapply(X = params[c("a","qEIP","nEIP","muH","r")], FUN = check_double)
  if(all(boolVec)){
    stop(paste0(names(boolVec)[boolVec], collapse = ", "), " are improperly specified")
  }

}

#  used directly in spn_hazards
check_params_SEIR <- function(params){

  # check basic sis params
  check_params_sis(params = params)

  # check extra SEIR stuff
  if(check_double(params$nu)){
    stop("nu is improperly specified in params")
  }

}

# used in decoupled sampling to aggregate the number
# of infectious mosquitoes per genotype and total number of SEI
# female mosquitoes to pass to the human model.
# similar to summarize_females_epi but this function
# is called at every time step to send to the human model.
# can probably consolidate the two functions later
# also only works with one node for now
aggregate_female_SEI <- function(spn_P, state) {
  f_idx <- lapply(X = spn_P$ix, '[[', 'females')
  f_idx[vapply(X = f_idx, FUN = is.null, FUN.VALUE = logical(length = 1))] <- NULL
  # get genotype names
  genotypes <- rownames(f_idx[[1]][,,"I"])

  # sum over mated genotypes to obtain aggregated number of females per genotype
  infectious_mosquitoes <- rep(0, length(genotypes))
  for(i in 1:length(genotypes)) {
    infectious_mosquitoes[i] <- sum(state[f_idx[[1]][,,"I"][genotypes[i],]])
  }
  
  names(infectious_mosquitoes) <- genotypes
  total_mosquitoes <- sum(state[f_idx[[1]]])

  return(list(I_V = infectious_mosquitoes, N_V = total_mosquitoes))
}


# internal formatting of CSV output labels for simulation
# pipeline analysis
generate_Imperial_human_state_labels <- function(na) { 
  # get 0-indexed indices
  indices <- 0:(na-1)
  # get labels in the form of <lower index>_<higher_index>
  labels <- paste(
    formatC(indices, width=2, flag="0"), 
    formatC((indices+1), width=2, flag="0"),
    sep="_"
  )
  
  # add state labels
  human_state_labels <- c(
    paste0("S", labels),
    paste0("T", labels),
    paste0("D", labels),
    paste0("A", labels),
    paste0("U", labels),
    paste0("P", labels),
    paste0("ICA", labels),
    paste0("IB", labels),
    paste0("ID", labels),
    paste0("IVA", labels),
    paste0("clin_inc", labels),
    paste0("mort", labels)
  )

  return(human_state_labels)
}