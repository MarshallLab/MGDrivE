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
