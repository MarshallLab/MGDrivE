################################################################################
#
#   MGDrivE2: equilibrium calculations for epidemiological extension
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   October 2019
#
################################################################################

###############################################################################
#   equilibrium for SIS
###############################################################################
# This is extremely similar to the SEIR equilibrium function, as much as possible
#  has been turned into a function and reused. However, because of the different
#  shapes for defining humans, several small things had to change, and lots of this
#  ends up straight copied. If there is a better way to do this, it would be great
#  to find.

#' Calculate Equilibrium for Mosquito SEI - Human SIS Model
#'
#' Given prevalence of disease in humans (modeled as an SIS: Susceptible-Infected-Susceptible
#' process with birth and death) and entomological parameters of transmission, this
#' function calculates the quasi-stationary distribution of adult female mosquitoes
#' across SEI (Susceptible-Exposed-Infectious) stages, allowing for Erlang distributed E stage.
#'
#' This function handles 3 types of nodes: Human only, mosquito only, and nodes
#' with both. These nodes are set using the \code{node_list} parameter.
#' Mosquito-only node equilibrium calls \code{\link{equilibrium_lifeycle}}, which
#' follows one of two models: classic logistic dynamics or the Lotka-Volterra
#' competition model. This is determined by the parameter \code{log_dd}, and it
#' changes elements of the return list: \code{K} is returned for logistic dynamics,
#' or \code{gamma} is returned for Lotka-Volterra dynamics. This
#' is parameterized with the \code{NF} parameter to define the adult female numbers.
#' This parameter only needs to be supplied if there are mosquito-only nodes.
#'
#' Human-only nodes don't require any equilibrium calculations. These nodes use
#' the \code{NH} and \code{pop_ratio_H} to set adult human populations and
#' infection rates in nodes. These two parameters only need to be supplied
#' if there are human-only nodes.
#'
#' For human and mosquito nodes, this function calls \code{\link{make_Q_SEI}} to construct the
#' infinitesimal generator matrix which is used to solve for the quasi-stationary
#' (stochastic) or equilibrium (deterministic) distribution of mosquitoes over stages.
#' Parameters are provided by \code{params}.
#'
#' For information on the method used to solve this distribution, see section
#' "3.1.3 Nonsingularity of the Subintensity Matrix" of:
#'    * Bladt, Mogens, and Bo Friis Nielsen. Matrix-exponential distributions in
#'    applied probability. Vol. 81. New York: Springer, 2017.
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
#' The \code{params} argument supplies all of the ecological and epidemiological
#' parameters necessary to calculate equilibrium values. This is used to set the
#' initial population distribution and during the simulation to maintain equilibrium.
#' This \code{params} must include the following named parameters, noted as being
#' the same as lifecycle parameters, or new for the epidemiological equilibrium
#'  * \strong{(Lifecycle parameters)}
#'    * \code{qE}: inverse of mean duration of egg stage
#'    * \code{nE}: shape parameter of Erlang-distributed egg stage
#'    * \code{qL}: inverse of mean duration of larval stage
#'    * \code{nL}: shape parameter of Erlang-distributed larval stage
#'    * \code{qP}: inverse of mean duration of pupal stage
#'    * \code{nP}: shape parameter of Erlang-distributed pupal stage
#'    * \code{muE}: egg mortality
#'    * \code{muL}: density-independent larvae mortality
#'    * \code{muP}: pupae mortality
#'    * \code{muF}: adult female mortality
#'    * \code{muM}: adult male mortality
#'    * \code{beta}: egg-laying rate, daily
#'    * \code{nu}: mating rate of unmated females
#'  * \strong{(Epidemiological parameters)}
#'    * \code{NH}: number of humans, can be a vector
#'    * \code{X}: prevalence in humans, can be a vector
#'    * \code{NFX}: number of female mosquitoes, only required if any prevalence (X) is zero
#'    * \code{b}: mosquito to human transmission efficiency, can be a vector
#'    * \code{c}: human to mosquito transmission efficiency, can be a vector
#'    * \code{r}: rate of recovery in humans (1/duration of infectiousness)
#'    * \code{muH}: death rate of humans (1/avg lifespan)
#'    * \code{f}: rate of blood feeding
#'    * \code{Q}: human blood index
#'    * \code{qEIP}: related to scale parameter of Gamma distributed EIP (1/qEIP is mean length of EIP)
#'    * \code{nEIP}: shape parameter of Gamma distributed EIP
#'
#' The return list contains all of the parameters necessary later in the simulations.
#'
#' For equilibrium without epidemiological parameters, see \code{\link{equilibrium_lifeycle}}.
#' For equilibrium with latent humans (SEIR dynamics), see \code{\link{equilibrium_SEI_SEIR}}.
#'
#' For examples of using this function, see:
#' \code{vignette("lifecycle-node", package = "MGDrivE2")}
#'
#' @param params a named list of parameters (see details)
#' @param node_list a character vector specifying what type of nodes to create;
#' (m = a node with only mosquitoes, h = a node with only humans, b = a node with both humans and mosquitoes)
#' @param NF vector of female mosquitoes at equilibrium, for mosquito-only nodes
#' @param phi sex ratio of mosquitoes at emergence
#' @param NH vector of humans at equilibrium, for human-only nodes
#' @param log_dd Boolean: TRUE implies logistic density dependence, FALSE implies Lotka-Volterra model
#' @param spn_P the set of places (P) (see details)
#' @param pop_ratio_Aq May be empty; if not, a named vector or matrix. (see details)
#' @param pop_ratio_F May be empty; if not, a named vector or matrix. (see details)
#' @param pop_ratio_M May be empty; if not, a named vector or matrix. (see details)
#' @param pop_ratio_H Prevalence in human-only nodes
#' @param cube an inheritance cube from the \code{MGDrivE} package (e.g. \code{\link[MGDrivE]{cubeMendelian}})
#'
#'
#' @return a vector of the equilibrium number of females in each SEI stage
#'
#' @importFrom Matrix solve
#'
#' @export
equilibrium_SEI_SIS <- function(params, node_list="b",NF=NULL,phi=0.5, NH=NULL,log_dd=TRUE,
                                spn_P, pop_ratio_Aq=NULL, pop_ratio_F=NULL,
                                pop_ratio_M=NULL, pop_ratio_H=1, cube){

  # check things first
  # check params
  tCheck <- vapply(X = unlist(x = params,recursive = TRUE, use.names = FALSE),
                   FUN = check_double, FUN.VALUE = logical(length = 1))
  if(any( tCheck )){
    stop("A member of params is not a standard number.\n\tPlease check and try again.")
  }

  # checks to perform if there are mixed nodes
  numBoth <- sum(node_list=="b")
  if(numBoth){
    if(any(numBoth != c(length(params$NH),length(params$b),length(params$c),length(params$X)) ) ){
      stop("Parameters NH, b, c, and X must be the same length as the number of mosquito/human patches")
    }
    if(any(params$X == 0)){
      if(length(params$X) != length(params$NFX)){
        stop(paste0("Some prevalence (params$X) in humans is 0, please specify NFX",
                    " to set female population. NFX must be length(number of 'b' nodes)."))
      }
    }
  } # end checks for mixed nodes

  # checks for human-only nodes
  numH <- sum(node_list=="h")
  if(numH){
    # total and prevalence in humans
    if(all(numH != c(length(NH),length(pop_ratio_H)) ) ){
      stop("NH and pop_ratio_H must be the same length as the number of human-only patches")
    }
  }

  if(sum(node_list=="m") != length(NF)){
    stop("NF must be the same length as the number of mosquito-only patches")
  }

  num_nodes <- length(node_list)
  eMsg <- paste0("population ratios (popRatio_{Aq,F,M} must be one of three values:\n",
           "    NULL (defult), where genotypes are taken from the cube.\n",
           "    a named vector of genotype ratios, matching genotypes in the cube.\n",
           "    a matrix, with nRows == num_patches, and column names matching genotypes in the cube.")
  if(!any(is.null(pop_ratio_Aq) || is.null(dim(pop_ratio_Aq)) || dim(pop_ratio_Aq)[1]==num_nodes)){
    stop(eMsg)
  }
  if(!any(is.null(pop_ratio_F) || is.null(dim(pop_ratio_F)) || dim(pop_ratio_F)[1]==num_nodes)){
    stop(eMsg)
  }
  if(!any(is.null(pop_ratio_M) || is.null(dim(pop_ratio_M)) || dim(pop_ratio_M)[1]==num_nodes)){
    stop(eMsg)
  }


  # create population objects depending on what nodes are here
  #  This function returns objects directly to this environment.
  #  "fSEI", "mosyHList", "hPopAq", "hPopF", "hPopM", "mosyList", "mPopAq", "mPopF", "mPopM"
  set_populations(node_list = node_list,params = params,phi = phi,log_dd = log_dd,
                  NF = NF,pop_ratio_Aq = pop_ratio_Aq,pop_ratio_F = pop_ratio_F,
                  pop_ratio_M = pop_ratio_M,cube = cube,fem_func = base_female_SIS)


  # begin setting up return object
  #  list for all the things
  retList <- list()

  #  start filling parameters
  retList$params <- c(params[c("qE","nE","qL","nL","qP","nP","muE","muL","muP",
                              "muF","muM","beta","muH","qEIP","nEIP","r","nu")],
                      list("a"=params$f*params$Q))

  # make initial conditions matrix
  #  columns are egg, larva, pupa, male, female infections, human S/I
  nE <- params$nE
  nL <- params$nL
  nP <- params$nP
  nELP <- nE+nL+nP
  nEIP <- params$nEIP

  M0 <- setNames(object = numeric(length = length(spn_P$u)), nm = spn_P$u)
  dPar <- numeric(length = num_nodes)
  initMat <- matrix(data = NA, nrow = num_nodes,
                     ncol = (nELP + 1) + (nEIP+2) + 2,
                     dimnames = list(1:num_nodes,c(paste0("E",1:nE),paste0("L",1:nL),
                                                  paste0("P",1:nP),"M",
                                                  c("F_S",paste0("F_E",1:nEIP),"F_I"),
                                                  "H_S","H_I")) )


  # loop over node_list, fill M0, inits mat, and density parameter
  idxB <- 1
  idxM <- 1
  idxH <- 1
  for(node in 1:num_nodes){
    if(node_list[node] == "b"){
      # S/I humans
      # calc breakdown
      hRatio <- params$NH[idxB] * c(1-params$X[idxB],params$X[idxB])

      # fill S/I humans
      M0[spn_P$ix[[node]]$humans] <- hRatio

      # fill mosquitoes
      M0[spn_P$ix[[node]]$egg[ ,colnames(hPopAq)] ] <- outer(X = mosyHList$init[idxB,1:nE],Y = hPopAq[idxB, ], FUN = "*")
      M0[spn_P$ix[[node]]$larvae[ ,colnames(hPopAq)] ] <- outer(X = mosyHList$init[idxB,1:nL+nE],
                                                               Y = hPopAq[idxB, ], FUN = "*")
      M0[spn_P$ix[[node]]$pupae[ ,colnames(hPopAq)] ] <- outer(X = mosyHList$init[idxB,1:nP+nE+nL],
                                                              Y = hPopAq[idxB, ], FUN = "*")
      M0[spn_P$ix[[node]]$males[colnames(hPopM)] ] <- outer(X = mosyHList$init[idxB,nELP+1],
                                                            Y = hPopM[idxB, ], FUN = "*")
      # females have SEI distribution
      for(SEI in 1:(nEIP+2)){
        M0[spn_P$ix[[node]]$females[colnames(hPopF),colnames(hPopM),SEI] ] <- fSEI[idxB,SEI] *
                                          as.vector(outer(X = hPopF[idxB, ],Y = hPopM[idxB, ]))
      }

      # fill inits conditions
      initMat[node,1:(nELP+1)] <- mosyHList$init[idxB, -(nELP+2)]
      initMat[node,1:(nEIP+2) + (nELP+1)] <- fSEI[idxB, ]
      initMat[node,1:2 + (nELP+nEIP+3)] <- hRatio

      # fill density parameter
      #  there is only 1 element in params list, K or gamma, so [[1]] gets it right
      dPar[node] <- mosyHList$params[[1]][idxB]

      # increment counter
      idxB <- idxB + 1

    } else if(node_list[node] == "h"){
      # calc breakdown
      hRatio <- NH[idxH] * c(1-pop_ratio_H[idxH],pop_ratio_H[idxH])

      # fill S/I humans
      M0[spn_P$ix[[node]]$humans] <- hRatio

      # fill inits conditions
      initMat[node,1:2 + (nELP+nEIP+3)] <- hRatio

      # increment counter
      idxH <- idxH+1

    } else if(node_list[node] == "m"){
      # fill M0 conditions
      M0[spn_P$ix[[node]]$egg[ ,colnames(mPopAq)] ] <- outer(X = mosyList$init[idxM,1:nE],Y = mPopAq[idxM, ], FUN = "*")
      M0[spn_P$ix[[node]]$larvae[ ,colnames(mPopAq)] ] <- outer(X = mosyList$init[idxM,1:nL+nE],
                                                               Y = mPopAq[idxM, ], FUN = "*")
      M0[spn_P$ix[[node]]$pupae[ ,colnames(mPopAq)] ] <- outer(X = mosyList$init[idxM,1:nP+nE+nL],
                                                              Y = mPopAq[idxM, ], FUN = "*")
      M0[spn_P$ix[[node]]$males[colnames(mPopM)] ] <- outer(X = mosyList$init[idxM,nELP+1],
                                                            Y = mPopM[idxM, ], FUN = "*")
      # Assume all females are uninfected
      M0[spn_P$ix[[node]]$females[colnames(mPopF),colnames(mPopM),1] ] <- outer(X = mosyList$init[idxM,nELP+2],
                                                                            Y = as.vector(outer(X = mPopF[idxM, ],
                                                                                                Y = mPopM[idxM, ])),
                                                                            FUN = "*")

      # fill inits conditions
      initMat[node,1:(nELP+2)] <- mosyList$init[idxM, ]

      # fill density parameter
      #  there is only 1 element in params list, K or gamma, so [[1]] gets it right
      dPar[node] <- mosyList$params[[1]][idxM]

      # increment counter
      idxM <- idxM + 1

    } else {
      stop(paste0("error: bad entry in node_list, ",node_list[node]))
    }

  } # end loop over nodes


  # put things into return list
  if(log_dd){
    retList$params <- c(retList$params,"K"=list(dPar))
  } else {
    retList$params <- c(retList$params,"gamma"=list(dPar))
  }
  retList$init <- initMat
  retList$M0 <- M0

  return(retList)

}


###############################################################################
# Populations base func
###############################################################################
set_populations <- function(node_list,params,phi,log_dd,NF,
                            pop_ratio_Aq,pop_ratio_F,pop_ratio_M,cube,
                            fem_func){

    # create objects depending on what nodes are here
  if(any(node_list %in% "b")){
    # female SEI equilibrium
    #  matrix over all human and mosquito habitats
    fSEI <- fem_func(params = params)

    # mosquito distribution for human/mosquito habitats
    mosyHList <- basic_eq_life(params = params,NF = rowSums(fSEI),phi = phi,log_dd = log_dd)

    # generate population distribution
    if(!any(is.null(pop_ratio_Aq) || length(pop_ratio_Aq)==1)){
      hAq <- pop_ratio_Aq[which(node_list == "b")]
    } else {
      hAq <- pop_ratio_Aq
    }
    if(!any(is.null(pop_ratio_F) || length(pop_ratio_F)==1)){
      hF <- pop_ratio_F[which(node_list == "b")]
    } else {
      hF <- pop_ratio_F
    }
    if(!any(is.null(pop_ratio_M) || length(pop_ratio_M)==1)){
      hM <- pop_ratio_M[which(node_list == "b")]
    } else {
      hM <- pop_ratio_M
    }

    # pass things back to parent environment
    assign(x = "fSEI", value = fSEI, pos = parent.frame())
    assign(x = "mosyHList", value = mosyHList, pos = parent.frame())
    assign(x = "hPopAq", value = calc_Ratio_Aq(pop_ratio_Aq = hAq, cube = cube, num_patch = nrow(fSEI)),
           pos = parent.frame())
    assign(x = "hPopF", value = calc_Ratio_F(pop_ratio_F = hF, cube = cube, num_patch = nrow(fSEI)),
           pos = parent.frame())
    assign(x = "hPopM", value = calc_Ratio_M(pop_ratio_M = hM,cube = cube,num_patch = nrow(fSEI)),
           pos = parent.frame())

  } # end human/mosquito node calcs
  if(any(node_list %in% "m")){
    # setup mosquito only equilibrium
    #  list of initial distributions
    #  "params" element with gamma or K
    mosyList <- basic_eq_life(params = params,NF = NF,phi = phi,log_dd = log_dd)

    # generate population distribution
    if(!any(is.null(pop_ratio_Aq) || length(pop_ratio_Aq)==1)){
      mAq <- pop_ratio_Aq[which(node_list == "m")]
    } else {
      mAq <- pop_ratio_Aq
    }
    if(!any(is.null(pop_ratio_F) || length(pop_ratio_F)==1)){
      mF <- pop_ratio_F[which(node_list == "m")]
    } else {
      mF <- pop_ratio_F
    }
    if(!any(is.null(pop_ratio_M) || length(pop_ratio_M)==1)){
      mM <- pop_ratio_M[which(node_list == "m")]
    } else {
      mM <- pop_ratio_M
    }

    # pass things back to parent environment
    assign(x = "mosyList", value = mosyList, pos = parent.frame())
    assign(x = "mPopAq", value = calc_Ratio_Aq(pop_ratio_Aq = mAq, cube = cube, num_patch = nrow(mosyList$init)),
           pos = parent.frame())
    assign(x = "mPopF", value = calc_Ratio_F(pop_ratio_F = mF, cube = cube, num_patch = nrow(mosyList$init)),
           pos = parent.frame())
    assign(x = "mPopM", value = calc_Ratio_M(pop_ratio_M = mM,cube = cube,num_patch = nrow(mosyList$init)),
           pos = parent.frame())
  } # end mosquito node calcs
}



###############################################################################
# equilibrium base func
###############################################################################
# this calculates equilibrium females, SEI, based on humans and epi info
# vectorized over: NH, X, b, and c
# optional parameter, NFX, if X = 0
# only used in equilibrium_SEI_SIS
base_female_SIS <- function(params){

  with(params,{

    # compund values
    a <- f*Q # HBR

    # human pop
    SH <- NH*(1-X)
    IH <- NH*(X)

    # number of infectious mosquitos required to sustain transmission
    IV <- (IH*NH*(r + muH)) / (a*b*SH)

    # setup female pops for return
    nNode <- length(NH)
    femPop <- matrix(data = 0, nrow = nNode, ncol = nEIP+2,
                     dimnames = list(1:nNode,c("S",paste0("E",1:nEIP),"I")))

    # loop over all nodes and compute!
    for(node in 1:nNode){
      # check if there are any infected humans
      #  the eq calcs fail if there aren't
      if(X[node] != 0){
        # make the generator matrix
        Qmat <- make_Q_SEI(q=qEIP,n=nEIP,mu=muF,c=c[node],a=a,x=X[node])

        # compute Green matrix and condition on starting in S compartment
        D0 <- Qmat[1:(nrow(Qmat)-1),1:(ncol(Qmat)-1)]
        # D0inv <- solve(-D0)
        # ttd <- D0inv[1,]
        ttd <- solve(-D0)[1,]
        ttd <- ttd/sum(ttd)

        # total number of mosquitos
        NV <- IV[node] / tail(ttd,1)

        # set female mosquito distribution
        femPop[node, ] <- NV*ttd
      } else {
        # set female mosquito distribution
        femPop[node,"S"] <- NFX[node]
      }

    } # end node loop

    # return matrix of female mosquitoes
    return(femPop)

  }) # end with
}


###############################################################################
# rate matrix
###############################################################################

#' Rate Matrix (Q) for Adult Mosquito SEI Dynamics
#'
#' Construct the infinitesimal generator matrix for (individual) adult female
#' infection dynamics. Adult females follow SEI (Susceptible-Exposed-Infectious)
#' style dynamics with a Gamma distributed EIP, with a mean duration 1/q and
#' variance 1/nq^2 (following shape-scale parameterization, EIP ~ Gamma(n,1/nq)).
#' This function only constructs the rate matrix for either a single mosquito or
#' cohort that all emerged at the same time (the rate matrix for a population
#' with emergence is infinite in dimension).
#'
#' @param q related to scale parameter of Gamma distributed EIP (1/q is mean length of EIP)
#' @param n shape parameter of Gamma distributed EIP
#' @param mu mosquito mortality rate
#' @param c human to mosquito transmission efficiency
#' @param a human biting rate
#' @param x prevalence of disease in humans
#'
#' @return rate matrix for a single (emergence) cohort of SEI mosquito
#'
#' @importFrom Matrix Matrix
#'
make_Q_SEI <- function(q, n, mu, c, a, x){

  # check nEIP is at least 1
  stopifnot(n >= 1)

  states <- c("S",paste0("E",1:n),"I","D")
  D <- (n+3)
  Q <- Matrix::Matrix(data = 0.0,nrow=D,ncol=D,dimnames=list(states,states),
                      sparse = FALSE,doDiag = FALSE)

  # FOI
  lambda <- a*c*x

  # Sv
  Q[1,2] <- lambda
  Q[1,D] <- mu
  Q[1,1] <- -sum(Q[1,]) # because of floating point inaccuracies

  # Ev(1)
  Q[2,3] <- q*n
  Q[2,D] <- mu
  Q[2,2] <- -sum(Q[2,])
  i <- 3

  # EIP
  if(n > 1){

    # Ev(2,...,EIP)
    for(i in (2:n)+1){
      Q[i,i+1] <- q*n
      Q[i,D] <- mu
      Q[i,i] <- -sum(Q[i,])
    }
    i <- i + 1

  }

  # Iv
  Q[i,i] <- -mu
  Q[i,D] <- mu

  # D
  Q[D,D] <- 0

  return(Q)
}
