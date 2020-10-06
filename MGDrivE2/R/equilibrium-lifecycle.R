################################################################################
#
#   MGDrivE2: equilibrium calculations for mosquito lifecycle model
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   October 2019
#
################################################################################

################################################################################
# equilibrium for logistic density dependence
################################################################################

#' Calculate Equilibrium for Lifecycle Model (Logistic or Lotka-Volterra)
#'
#' This function calculates deterministic equilibria for the mosquito lifecycle
#' model.
#'
#' Equilibrium can be calculated using one of two models: classic logistic dynamics
#' or following the Lotka-Volterra competition model. This is determined by the
#' parameter \code{log_dd}, and it changes elements of the return list: \code{K} is
#' returned for logistic dynamics, or \code{gamma} is returned for Lotka-Volterra
#' dynamics.
#'
#' The places (\code{spn_P}) object is generated from one of the following:
#' \code{\link{spn_P_lifecycle_node}}, \code{\link{spn_P_lifecycle_network}},
#' \code{\link{spn_P_epiSIS_node}}, \code{\link{spn_P_epiSIS_network}},
#' \code{\link{spn_P_epiSEIR_node}}, or \code{\link{spn_P_epiSEIR_network}}.
#'
#' The initial population genotype ratios are set by supplying the \code{pop_ratio_Aq},
#' \code{pop_ratio_F}, and \code{pop_ratio_M} values. The default value is NULL,
#' and the function will use the wild-type alleles provided in the \code{cube}
#' object. However, one can supply
#' several different objects to set the initial genotype ratios. All genotypes provided
#' must exist in the \code{cube} (this is checked by the function). If a single, named vector
#' is provided, then all patches will be initialized with the same ratios. If a
#' matrix is provided, with the number of columns (and column names) giving the
#' initial genotypes, and a row for each patch, each patch can be set to a different
#' initial ratio. The three parameters do not need to match each other.
#'
#' The \code{params} argument supplies all of the ecological parameters necessary
#' to calculate equilibrium values. This is used to set the initial population
#' distribution and during the simulation to maintain equilibrium. \code{params}
#' must include the following named parameters:
#'  * \code{qE}: inverse of mean duration of egg stage
#'  * \code{nE}: shape parameter of Erlang-distributed egg stage
#'  * \code{qL}: inverse of mean duration of larval stage
#'  * \code{nL}: shape parameter of Erlang-distributed larval stage
#'  * \code{qP}: inverse of mean duration of pupal stage
#'  * \code{nP}: shape parameter of Erlang-distributed pupal stage
#'  * \code{muE}: egg mortality
#'  * \code{muL}: density-independent larvae mortality
#'  * \code{muP}: pupae mortality
#'  * \code{muF}: adult female mortality
#'  * \code{muM}: adult male mortality
#'  * \code{beta}: egg-laying rate, daily
#'  * \code{nu}: mating rate of unmated females
#'
#' The return list contains all of the \code{params} parameters, along with the
#' density-dependent parameter, either \code{K} or \code{gamma}. These are the
#' parameters necessary later in the simulations. This was done for compatibility
#' with \code{\link{equilibrium_SEI_SIS}}, which requires several extra parameters
#' not required further in the simulations.
#'
#' For equilibrium with epidemiological parameters, see \code{\link{equilibrium_SEI_SIS}}.
#' For equilibrium with latent humans (SEIR dynamics), see \code{\link{equilibrium_SEI_SEIR}}.
#'
#'
#' @param params a named list of parameters (see details)
#' @param NF vector of female mosquitoes at equilibrium, for every population in the environment
#' @param phi sex ratio of mosquitoes at emergence
#' @param log_dd Boolean: TRUE implies logistic density dependence, FALSE implies Lotka-Volterra model
#' @param spn_P the set of places (P) (see details)
#' @param pop_ratio_Aq May be empty; if not, a named vector or matrix. (see details)
#' @param pop_ratio_F May be empty; if not, a named vector or matrix. (see details)
#' @param pop_ratio_M May be empty; if not, a named vector or matrix. (see details)
#' @param cube an inheritance cube from the \code{MGDrivE} package (e.g. \code{\link[MGDrivE]{cubeMendelian}})
#'
#' @return a list with 3 elements: \code{init} a matrix of equilibrium values for every life-cycle stage,
#' \code{params} a list of parameters for the simulation, \code{M0} a vector of initial conditions
#'
#' @examples
#' \dontrun{
#'   # parameters, see vignette MGDrivE2: One Node Lifecycle Dynamics
#'   theta <- list(qE = 1/4, nE = 2, qL = 1/3, nL = 3, qP = 1/6, nP = 2,
#'                 muE = 0.05, muL = 0.15, muP = 0.05, muF = 0.09, muM = 0.09,
#'                 beta = 16, nu = 1/(4/24) )
#'
#'   # spn_P needs to be setup elsewhere, see vignettes
#'
#'   eqLife <- equilibrium_lifeycle(params = theta, NF = 500, phi = 0.5,
#'                                  log_dd = TRUE, spn_P = spn_P,
#'                                  cube = MGDrivE::cubeMendelian() )
#' }
#'
#' @export
equilibrium_lifeycle <- function(params,NF,phi=0.5,log_dd=TRUE, spn_P,
                                 pop_ratio_Aq=NULL, pop_ratio_F=NULL,
                                 pop_ratio_M=NULL, cube){

  # checks
  #  length of adult vec matches numberof patches
  num_patch <- length(spn_P$ix)
  if(length(NF) != num_patch){
    stop("NF must be the same length as the number of patches")
  }

  # check params
  tCheck <- vapply(X = unlist(x = params,recursive = TRUE, use.names = FALSE),
                   FUN = check_double, FUN.VALUE = logical(length = 1))
  if(any( tCheck )){
    stop("A member of params is not a standard number.\n\tPlease check and try again.")
  }


  # setup parameters and eq distribution
  mosyList <- basic_eq_life(params=params,NF=NF,phi=phi,log_dd=log_dd)

  # add rest of parameters
  mosyList$params <- c(params,mosyList$params)


  # setup genotype breakdowns
  popAq <- calc_Ratio_Aq(pop_ratio_Aq = pop_ratio_Aq, cube = cube, num_patch = num_patch)
  popF <- calc_Ratio_F(pop_ratio_F = pop_ratio_F, cube = cube, num_patch = num_patch)
  popM <- calc_Ratio_M(pop_ratio_M = pop_ratio_M,cube = cube,num_patch = num_patch)

  # setup M0 population
  M0 <- setNames(object = numeric(length = length(spn_P$u)), nm = spn_P$u)
  nE <- params$nE
  nL <- params$nL
  nP <- params$nP

  # fill population
  for(node in 1:num_patch){
    M0[spn_P$ix[[node]]$egg[ ,colnames(popAq)] ] <- outer(X = mosyList$init[node,1:nE],Y = popAq[node, ], FUN = "*")
    M0[spn_P$ix[[node]]$larvae[ ,colnames(popAq)] ] <- outer(X = mosyList$init[node,1:nL+nE],Y = popAq[node, ], FUN = "*")
    M0[spn_P$ix[[node]]$pupae[ ,colnames(popAq)] ] <- outer(X = mosyList$init[node,1:nP+nE+nL],Y = popAq[node, ], FUN = "*")
    M0[spn_P$ix[[node]]$males[colnames(popM)] ] <- outer(X = mosyList$init[node,nE+nL+nP+1],Y = popM[node, ], FUN = "*")
    M0[spn_P$ix[[node]]$females[colnames(popF),colnames(popM),1] ] <- outer(X = mosyList$init[node,nE+nL+nP+2],
                                                                          Y = as.vector(outer(popF[node, ], popM[node, ])),
                                                                          FUN = "*")
  }

  # return everythingc
  return(c(mosyList,list("M0"=M0)))
}


################################################################################
# equilibrium base func
################################################################################
# base function to calculate just the equilibrium
# vectorized over NF
#  returns a matrix of initial conditions
# used in: equilibrium_lifeycle and equilibrium_SEI_SIS
basic_eq_life <- function(params,NF,phi,log_dd){

  # setup "global" return list
  eqList <- list("params"=NULL)

  # calculate equilibrium distribution
  with(params,{
    # get number of nodes to setup objects
    nNode <- length(NF)

    # setup prelim objects
    P_eq <- matrix(data = 0, nrow = nNode, ncol = nP)
    L_eq <- matrix(data = 0, nrow = nNode, ncol = nL)
    E_eq <- matrix(data = 0, nrow = nNode, ncol = nE)


    # calc eq for P_np and NM
    P_eq[ ,nP] <- (NF*muF) / (nP*qP*phi)
    NM_eq <- ((1-phi)*qP*nP*P_eq[ ,nP]) / muM

    # calc eq for rest of P_i's
    if(nP > 1){
      for(i in (nP-1):1){
        P_eq[ ,i] <- ((muP + (qP*nP)) / (qP*nP)) * P_eq[ ,i+1]
      }
    }

    # calc eq for E_i's
    E_eq[ ,1] <- (beta*NF) / (muE + (qE*nE))
    if(nE > 1){
      for(i in 2:nE){
        E_eq[ ,i] <- (qE*nE*E_eq[ ,i-1]) / (muE + (qE*nE))
      }
    }

    # calc eq for L_nl
    L_eq[ ,nL] <- ((muP + (qP*nP)) / (qL*nL)) * P_eq[ ,1]


    # check density dependence type
    #  calculate rest of larvae at EQ
    #  return proper EQ parameterization
    if(log_dd){
      # LOGISTIC DENSITY DEPENDENCE
      # calc eq for rest of L_i's
      L_nL <- L_eq[ ,nL]
      E_nE <- E_eq[ ,nE]
      if(nL > 1){
        for(i in 1:(nL-1)){
          L_eq[ ,i] <- ((L_nL^(i/nL)) * (E_nE^((nL-i)/nL)) * (nE^((nL-i)/nL)) * (qE^((nL-i)/nL))) /
                         ((nL^((nL-i)/nL)) * (qL^((nL-i)/nL)))
        }
      }

      # calc eq for K
      eqList$params$K <<- rowSums(L_eq) / (((qE*nE*E_nE)/(muL*L_eq[ ,1])) - ((qL*nL)/muL) - 1)

    } else {
        # LOTKA-VOLTERRA DENSITY DEPENDENCE
        # calc eq for rest of L_i's
        if(nL > 1){
          L_eq[ ,1] <- ((E_eq[ ,nE]^((nL-1)/nL)) * (L_eq[ ,nL]^(1/nL)) * (nE^((nL-1)/nL)) * (qE^((nL-1)/nL))) /
            ((nL^((nL-1)/nL)) * (qL^((nL-1)/nL)))
          for(i in 2:(nL-1)){
            L_eq[ ,i] <- (L_eq[ ,1]^i) * (((nL*qL)/(E_eq[ ,nE]*nE*qE))^(i-1))
          }
        }

        # calc eq for gamma
        eqList$params$gamma <<- (((E_eq[ ,nE]*nE*qE) / L_eq[ ,1]) - (muL + (nL*qL))) / rowSums(L_eq)

      }

    # initial state for all nodes
    eqList$init <<- matrix(data = cbind(E_eq,L_eq,P_eq,NM_eq,NF), nrow = nNode, ncol = nE+nL+nP+2,
                           dimnames = list(1:nNode, c(paste0("E",1:nE),paste0("L",1:nL),paste0("P",1:nP),"NM","NF")))

  })


  # return list with initial conditions and density-dependent parameters
  return(eqList)
}


################################################################################
# genotype ratios
################################################################################
# these are stolen directly from MGDrivE
# handles: any number of initial conditions,none at all, different ones for each patch
#  defaults to the wild-type of the cube, including X/Y cubes
# used in: equilibrium_lifeycle and equilibrium_SEI_SIS
calc_Ratio_Aq <- function(pop_ratio_Aq,cube,num_patch){
  # setup larval initial pop ratio
  if(is.null(pop_ratio_Aq)){
    # default behaviour - this way nothing needs to be specified
    # setup vector of the correct length with all weights equal
    larRatio <- rep.int(x = 1, times = length(cube$wildType))
    # normalize to 1
    larRatio <- larRatio/sum(larRatio)
    # set all patches to the same thing
    pop_ratio_Aq <- matrix(data = larRatio, nrow = num_patch, ncol = length(cube$wildType),
                                 byrow = TRUE, dimnames = list(NULL,cube$wildType))

  } else if(is.null(dim(pop_ratio_Aq)) ){
    # behaviour of user supplied 1 patch worth of weights.
    # has to be supplied as a vector, the matrix stuff makes this way too difficult
    if(abs(sum(pop_ratio_Aq) - 1) > sqrt(.Machine$double.eps) ) stop('pop_ratio_Aq must sum to 1')
    # check that columns names are in inheritance cube
    if(is.null(names(pop_ratio_Aq)) || !all(names(pop_ratio_Aq) %in% cube$genotypesID)) {
      stop("Names for pop_ratio_Aq must be specified as one of the genotypesID in the inheritance cube.")
    }
    # set all patches equal
    pop_ratio_Aq <- matrix(data = pop_ratio_Aq, nrow = num_patch, ncol = length(pop_ratio_Aq),
                                 byrow = TRUE, dimnames = list(NULL,names(pop_ratio_Aq)) )

  } else if(dim(pop_ratio_Aq)[1] == num_patch){
    # behaviour if user supplies a matrix of probabilities.
    # each row is a different patch
    # check that all patches sum to 1
    if(any(abs(rowSums(pop_ratio_Aq) - 1) > sqrt(.Machine$double.eps)) ) stop('Each row of pop_ratio_Aq must sum to 1')
    # check that columns names are in inheritance cube
    if(is.null(colnames(pop_ratio_Aq)) || !all(colnames(pop_ratio_Aq) %in% cube$genotypesID)) {
      stop("Column names for pop_ratio_Aq must be specified as one of the genotypesID in the inheritance cube.")
    }
    # store
    pop_ratio_Aq <- pop_ratio_Aq

  } else {
    stop("pop_ratio_Aq has been miss specified.\n
         Left blank - default, all populations are the same and begin as wild-type individuals\n
         Vector - a named vector that sums to one. All populations will be the same\n
         Matrix - an num_patch by nGenotype matrix with column names and all rows sum to 1.
         Specifies each population individually.")
  }

  # return aquatic genotypes
  return(pop_ratio_Aq)
}

calc_Ratio_F <- function(pop_ratio_F,cube,num_patch){
  # setup female initial pop ratio
  if(is.null(pop_ratio_F)){
    # default behaviour - this way nothing needs to be specified

    # now, two options. This is not sex based, so only 1 wild-type
    #  or, sex based cube, we have an X and a Y, females are not Y
    if(length(cube$wildType) == 1){
      # only 1 wild-type, so no problem
      pop_ratio_F <- matrix(data = 1, nrow = num_patch, ncol = length(cube$wildType),
                                  dimnames = list(NULL,cube$wildType))
    } else if(length(cube$wildType) == 2){
      # one female and one male wild-type.
      # get index of female wild-type
      whichGeno <- grep(pattern = "Y", x = cube$wildType, fixed = TRUE, invert = TRUE)
      # setup default matrix
      pop_ratio_F <- matrix(data = 0, nrow = num_patch, ncol = length(cube$wildType),
                                  dimnames = list(NULL,cube$wildType))
      # set all to female genotype
      pop_ratio_F[ ,whichGeno] <- 1

    } else {
      stop("Default pop_ratio_F only handles 1 or 2 wild-type genotypes.\n
           Please provide a matrix specifying female initial genotypes for every patch.")
    }

  } else if(is.null(dim(pop_ratio_F)) ){
    # behaviour of user supplied 1 patch worth of weights.
    # has to be supplied as a vector, the matrix stuff makes this way too difficult
    if(abs(sum(pop_ratio_F) - 1) > sqrt(.Machine$double.eps) ) stop('pop_ratio_F must sum to 1')
    # check that columns names are in inheritance cube
    if(is.null(names(pop_ratio_F)) || !all(names(pop_ratio_F) %in% cube$genotypesID)) {
      stop("Names for pop_ratio_F must be specified as one of the genotypesID in the inheritance cube.")
    }
    # set all patches equal
    pop_ratio_F <- matrix(data = pop_ratio_F, nrow = num_patch, ncol = length(pop_ratio_F),
                               byrow = TRUE, dimnames = list(NULL,names(pop_ratio_F)) )

  } else if(dim(pop_ratio_F)[1] == num_patch){
    # behaviour if user supplies a matrix of probabilities.
    # each row is a different patch
    # check that all patches sum to 1
    if(any(abs(rowSums(pop_ratio_F) - 1) > sqrt(.Machine$double.eps)) ) stop('Each row of pop_ratio_F must sum to 1')
    # check that columns names are in inheritance cube
    if(is.null(colnames(pop_ratio_F)) || !all(colnames(pop_ratio_F) %in% cube$genotypesID)) {
      stop("Column names for pop_ratio_F must be specified as one of the genotypesID in the inheritance cube.")
    }
    # store
    pop_ratio_F <- pop_ratio_F

  } else {
    stop("pop_ratio_F has been miss specified.\n
         Left blank - default, all populations are the same and begin as wild-type individuals\n
         Vector - a named vector that sums to one. All populations will be the same\n
         Matrix - an num_patch by nGenotype matrix with column names and all rows sum to 1.
         Specifies each population individually.")
  }

  # return female locations for each patch
  return(pop_ratio_F)
}

calc_Ratio_M <- function(pop_ratio_M,cube,num_patch){
  # setup male initial pop ratio
  if(is.null(pop_ratio_M)){
    # default behaviour - this way nothing needs to be specified

    # now, two options. This is not sex based, so only 1 wild-type
    #  or, sex based cube, we have an X and a Y, females are not Y
    if(length(cube$wildType) == 1){
      # only 1 wild-type, so no problem
      pop_ratio_M <- matrix(data = 1, nrow = num_patch, ncol = length(cube$wildType),
                                  dimnames = list(NULL,cube$wildType))
    } else if(length(cube$wildType) == 2){
      # one female and one male wild-type.
      # get index of female wild-type
      whichGeno <- grep(pattern = "Y", x = cube$wildType, fixed = TRUE)
      # setup default matrix
      pop_ratio_M <- matrix(data = 0, nrow = num_patch, ncol = length(cube$wildType),
                                  dimnames = list(NULL,cube$wildType))
      # set all to female genotype
      pop_ratio_M[ ,whichGeno] <- 1

    } else {
      stop("Default pop_ratio_M only handles 1 or 2 wild-type genotypes.\n
           Please provide a matrix specifying female initial genotypes for every patch.")
    }

  } else if(is.null(dim(pop_ratio_M)) ){
    # behaviour of user supplied 1 patch worth of weights.
    # has to be supplied as a vector, the matrix stuff makes this way too difficult
    if(abs(sum(pop_ratio_M) - 1) > sqrt(.Machine$double.eps) ) stop('pop_ratio_M must sum to 1')
    # check that columns names are in inheritance cube
    if(is.null(names(pop_ratio_M)) || !all(names(pop_ratio_M) %in% cube$genotypesID)) {
      stop("Names for pop_ratio_M must be specified as one of the genotypesID in the inheritance cube.")
    }
    # set all patches equal
    pop_ratio_M <- matrix(data = pop_ratio_M, nrow = num_patch, ncol = length(pop_ratio_M),
                                byrow = TRUE, dimnames = list(NULL,names(pop_ratio_M)) )

  } else if(dim(pop_ratio_M)[1] == num_patch){
    # behaviour if user supplies a matrix of probabilities.
    # each row is a different patch
    # check that all patches sum to 1
    if(any(abs(rowSums(pop_ratio_M) - 1) > sqrt(.Machine$double.eps)) ) stop('Each row of pop_ratio_M must sum to 1')
    # check that columns names are in inheritance cube
    if(is.null(colnames(pop_ratio_M)) || !all(colnames(pop_ratio_M) %in% cube$genotypesID)) {
      stop("Column names for pop_ratio_M must be specified as one of the genotypesID in the inheritance cube.")
    }
    # store
    pop_ratio_M <- pop_ratio_M

  } else {
    stop("pop_ratio_M has been miss specified.\n
         Left blank - default, all populations are the same and begin as wild-type individuals\n
         Vector - a named vector that sums to one. All populations will be the same\n
         Matrix - an num_patch by nGenotype matrix with column names and all rows sum to 1.
         Specifies each population individually.")
  }

  # return male genotypes for each patch
  return(pop_ratio_M)
}
