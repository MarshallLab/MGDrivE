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

#' Calculate Equilibrium for Decoupled Mosquito SEI model. Human states will be handled separately.
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
equilibrium_SEI_decoupled_mosy <- function(params, node_list="b",NF=NULL,phi=0.5, NH=NULL,log_dd=TRUE,
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
                     ncol = (nELP + 1) + (nEIP+2),
                     dimnames = list(1:num_nodes,c(paste0("E",1:nE),paste0("L",1:nL),
                                                  paste0("P",1:nP),"M",
                                                  c("F_S",paste0("F_E",1:nEIP),"F_I")
                                                  )) )


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

#' This function calculates the equilibrium values for the decoupled SIS human states.
#' Currently this only works in one node.
#' @param params a named list of parameters (see details)
#' @return a vector of the equilibrium number of humans in each SIS stage
#' @export
equilibrium_SEI_decoupled_human <- function(params) {
    idxB <- 1
    hRatio <- params$NH[idxB] * c(1-params$X[idxB],params$X[idxB])
    names(hRatio) <- c("H_S", "H_I")
    return(hRatio)
}

#' This function calculates the human equilibrium values for the decoupled Imperial model.
#' Requires the age structure of the population
#' Currently this only works in one node.
#' @param age_vector age structure of population (see vignette for example)
#' @param ft proportion of population seeking treatment
#' @param EIR desired annual EIR
#' @param model_param_list parameters for the
#' @return a matrix of the equilibrium number of humans in each Imperial stage by age, and immunity
#' @importFrom statmod gauss.quad.prob
equilibrium_Imperial_decoupled_human <- function(age_vector, ft, EIR, model_param_list)
{
  mpl <- model_param_list
  het_brackets <- 1

  ## Check Parameters
  if(!is.numeric(age_vector)) stop("age_vector provided is not numeric")
  if(!is.numeric(ft)) stop("ft provided is not numeric")
  if(!is.numeric(EIR)) stop("EIR provided is not numeric")

  ## population demographics
  age <- age_vector * mpl$DY
  na <- as.integer(length(age))  # number of age groups
  nh <- as.integer(het_brackets)  # number of heterogeneity groups
  h <- statmod::gauss.quad.prob(nh, dist = "normal")
  age0 <- 2
  age1 <- 10
  num_int <- mpl$num_int

  age_rate <- age_width <- age_mid_point <- den <- c()
  for (i in 1:(na-1))
  {
    age_width[i] <- age[i+1] - age[i]
    age_rate[i] <- 1/(age[i + 1] - age[i])  # vector of rates at which people leave each age group (1/age group width)
    age_mid_point[i] <- 0.5 * (age[i] + age[i + 1])  # set age group vector to the midpoint of the group

  }
  age_rate[na] = 0


  den <- 1/(1 + age_rate[1]/mpl$eta)
  for (i in 1:(na-1))
  {
    den[i+1] <- age_rate[i] * den[i]/(age_rate[i+1] + mpl$eta)  # proportion in each age_vector group
  }

  age59 <- which(age_vector * 12 > 59)[1] - 1  # index of age vector before age is >59 months
  age05 <- which(age_vector > 5)[1] - 1  # index of age vector before age is 5 years

  ## force of infection
  foi_age <- c()
  for (i in 1:na)
  {
    foi_age[i] <- 1 - (mpl$rho * exp(-age[i]/mpl$a0))  #force of infection for each age group
  }
  fden <- foi_age * den
  omega <- sum(fden)  #normalising constant

  ## heterogeneity
  het_x <- h$nodes
  het_wt <- h$weights
  den_het <- outer(den, het_wt)
  rel_foi <- exp(-mpl$sigma2/2 + sqrt(mpl$sigma2) * het_x)/sum(het_wt * exp(-mpl$sigma2/2 + sqrt(mpl$sigma2) * het_x))

  ## EIR
  EIRY_eq <- EIR  # initial annual EIR
  EIRd_eq <- EIRY_eq/mpl$DY
  EIR_eq   <- outer(foi_age, rel_foi) * EIRd_eq

  ## Immunity and FOI
  x_I <- den[1]/mpl$eta
  for (i in 2:na)
  {
    x_I[i] <- den[i]/(den[i - 1] * age_rate[i - 1])  #temporary variables
  }
  fd <- 1 - (1 - mpl$fD0)/(1 + (age/mpl$aD)^mpl$gammaD)

  # maternal immunity begins at a level proportional to the clinical
  # immunity of a 20 year old, this code finds that level
  age20i <- rep(0, na)
  for (i in 2:na)
  {
    age20i[i] <- ifelse(age[i] >= (20 * mpl$DY) & age[i - 1] < (20 * mpl$DY),
                        i, age20i[i - 1])
  }
  age20u <- as.integer(age20i[na])
  age20l <- as.integer(age20u - 1)
  age_20_factor <- (20 * mpl$DY - age[age20l] - 0.5 * age_width[age20l]) *
    2/(age_width[age20l] + age_width[age20u])

  # finding initial values for all immunity states
  IB_eq <- matrix(0, na, nh)
  FOI_eq <- matrix(0, na, nh)
  ID_eq <- matrix(0, na, nh)
  ICA_eq <- matrix(0, na, nh)
  IVA_eq <- matrix(0,na,nh)
  ICM_init_eq <- vector(length = nh, mode = "numeric")
  ICM_eq <- matrix(0, na, nh)
  IVM_init_eq <- vector(length = nh, mode = "numeric")
  IVM_eq <- matrix(0, na, nh)
  cA_eq <- matrix(0, na, nh)
  p_det_eq <- matrix(0, na, nh)
  FOIvij_eq <- matrix(0, na, nh)

  for (j in 1:nh)
  {
    for (i in 1:na)
    {
      IB_eq[i, j] <- (ifelse(i == 1, 0, IB_eq[i - 1, j]) +
                        EIR_eq[i,j]/(EIR_eq[i, j] * mpl$uB + 1) * x_I[i])/(1 + x_I[i]/mpl$dB)
      FOI_eq[i, j] <- EIR_eq[i, j] * ifelse(IB_eq[i, j] == 0, mpl$b0,
                                            mpl$b0 * ((1 - mpl$b1)/(1 + (IB_eq[i, j]/mpl$IB0)^mpl$kB) + mpl$b1))
      ID_eq[i, j] <- (ifelse(i == 1, 0, ID_eq[i - 1, j]) +
                        FOI_eq[i, j]/(FOI_eq[i, j] * mpl$uD + 1) * x_I[i])/(1 + x_I[i]/mpl$dID)
      ICA_eq[i, j] <- (ifelse(i == 1, 0, ICA_eq[i - 1, j]) +
                         FOI_eq[i,j]/(FOI_eq[i, j] * mpl$uCA + 1) * x_I[i])/(1 + x_I[i]/mpl$dCA)
      IVA_eq[i, j] <- (ifelse(i == 1, 0, IVA_eq[i - 1, j]) +
                         FOI_eq[i,j]/(FOI_eq[i, j] * mpl$uVA + 1) * x_I[i])/(1 + x_I[i]/mpl$dVA)
      p_det_eq[i, j] <- mpl$d1 + (1 - mpl$d1)/(1 + fd[i] * (ID_eq[i, j]/mpl$ID0)^mpl$kD)
      cA_eq[i, j] <- mpl$cU + (mpl$cD - mpl$cU) * p_det_eq[i, j]^mpl$gamma1
    }
  }

  # needs to be calculated after because it references ICA
  for (j in 1:nh)
  {
    for (i in 1:na)
    {
      ICM_init_eq[j] <- mpl$PM * (ICA_eq[age20l, j] + age_20_factor *
                                    (ICA_eq[age20u, j] - ICA_eq[age20l, j]))
      ICM_eq[i, j] <- ifelse(i == 1,
                             ICM_init_eq[j], ICM_eq[i - 1,j])/(1 + x_I[i]/mpl$dCM)
      IVM_init_eq[j] <- mpl$PVM * (IVA_eq[age20l, j] + age_20_factor *
                                    (IVA_eq[age20u, j] - IVA_eq[age20l, j]))
      IVM_eq[i, j] <- ifelse(i == 1,
                             IVM_init_eq[j], IVM_eq[i - 1,j])/(1 + x_I[i]/mpl$dVM)
    }
  }

  IC_eq <- ICM_eq + ICA_eq
  phi_eq <- mpl$phi0 * ((1 - mpl$phi1)/(1 + (IC_eq/mpl$IC0)^mpl$kC) + mpl$phi1)


  # human states
  gamma <- mpl$eta + c(age_rate[1:(na - 1)], 0)
  delta <- c(mpl$eta, age_rate[1:(na - 1)])

  betaT <- matrix(rep(mpl$rT + gamma, rep(nh, na)), ncol = nh, byrow = TRUE)
  betaD <- matrix(rep(mpl$rD + gamma, rep(nh, na)), ncol = nh, byrow = TRUE)
  betaP <- matrix(rep(mpl$rP + gamma, rep(nh, na)), ncol = nh, byrow = TRUE)

  aT <- FOI_eq * phi_eq * ft/betaT
  aD <- FOI_eq * phi_eq * (1 - ft)/betaD
  aP <- mpl$rT * aT/betaP

  Z_eq <- array(dim = c(na, nh, 4))
  Z_eq[1, , 1] <- den_het[1, ]/(1 + aT[1, ] + aD[1, ] + aP[1, ])
  Z_eq[1, , 2] <- aT[1, ] * Z_eq[1, , 1]
  Z_eq[1, , 3] <- aD[1, ] * Z_eq[1, , 1]
  Z_eq[1, , 4] <- aP[1, ] * Z_eq[1, , 1]

  for (j in 1:nh)
  {
    for (i in 2:na)
      {
      Z_eq[i, j, 1] <- (den_het[i, j] - delta[i] * (Z_eq[i - 1, j, 2]/betaT[i, j] +
                                                      Z_eq[i - 1, j, 3]/betaD[i, j] +
                                                      (mpl$rT *  Z_eq[i - 1, j, 2]/betaT[i, j]
                                                       + Z_eq[i - 1, j, 4])/betaP[i, j]))/(1 + aT[i, j] + aD[i, j] + aP[i, j])
      Z_eq[i, j, 2] <- aT[i, j] * Z_eq[i, j, 1] + delta[i] * Z_eq[i -
                                                                    1, j, 2]/betaT[i, j]
      Z_eq[i, j, 3] <- aD[i, j] * Z_eq[i, j, 1] + delta[i] * Z_eq[i -
                                                                    1, j, 3]/betaD[i, j]
      Z_eq[i, j, 4] <- aP[i, j] * Z_eq[i, j, 1] + delta[i] * (mpl$rT *
                                                                Z_eq[i - 1, j, 2]/betaT[i, j] + Z_eq[i - 1, j, 4])/betaP[i,j]

    }
  }

  Y_eq <- matrix(Z_eq[, , 1], nrow = na, ncol=nh)
  T_eq <- matrix(Z_eq[, , 2], nrow = na, ncol=nh)
  D_eq <- matrix(Z_eq[, , 3], nrow = na, ncol=nh)
  P_eq <- matrix(Z_eq[, , 4], nrow = na, ncol=nh)

  betaS <- apply(FOI_eq, MARGIN = 2, FUN = function(x, y)
  {
    x + y
  }, y = gamma)
  betaA <- apply(FOI_eq * phi_eq + mpl$rA, MARGIN = 2, FUN = function(x, y)
  {
    x + y
  }, y = gamma)
  betaU <- apply(FOI_eq + mpl$rU, MARGIN = 2, FUN = function(x, y)
  {
    x + y
  }, y = gamma)

  A_eq <- matrix(ncol = nh, nrow = na)
  U_eq <- matrix(ncol = nh, nrow = na)
  S_eq <- matrix(ncol = nh, nrow = na)

  for (i in 1:na)
  {
    for (j in 1:nh)
    {
      A_eq[i, j] <- (delta[i] * ifelse(i == 1, 0, A_eq[i - 1, j]) +
                       FOI_eq[i, j] * (1 - phi_eq[i, j]) * Y_eq[i, j] +
                       mpl$rD * D_eq[i,j])/(betaA[i, j] + FOI_eq[i, j] * (1 - phi_eq[i, j]))
      U_eq[i, j] <- (mpl$rA * A_eq[i, j] + delta[i] * ifelse(i == 1,
                                                             0, U_eq[i - 1, j]))/betaU[i, j]
      S_eq[i, j] <- Y_eq[i, j] - A_eq[i, j] - U_eq[i, j]
      FOIvij_eq[i, j] <- foi_age[i] * mpl$av0 * (mpl$cT * T_eq[i, j] + mpl$cD *
                                                   D_eq[i, j] + cA_eq[i, j] * A_eq[i, j] + mpl$cU * U_eq[i, j]) * rel_foi[j]/omega
    }
  }


  # number of infectious mosquitoes states
  FOIv_eq <- sum(FOIvij_eq)
  Iv_eq <- FOIv_eq * mpl$Surv0/(FOIv_eq + mpl$muF)
  Sv_eq <- mpl$muF * Iv_eq/(FOIv_eq * mpl$Surv0)
  Ev_eq <- 1 - Sv_eq - Iv_eq

  # mosquito density needed to give this EIR
  mv0 <- omega * EIRd_eq/(Iv_eq * mpl$av0)

  # pathogen prevalence
  prev_eq <- sum(A_eq) + sum(U_eq) + sum(D_eq) + sum(T_eq)

  fv <- 1 - (1-mpl$fvS)/(1+(age/mpl$av)^mpl$gammaV)
  immunity_multiplier <- ((IVA_eq+IVM_eq)/mpl$iv0)^mpl$kv
  severe_disease <- mpl$theta0*(mpl$theta1 + ((1-mpl$theta1)/(1+fv*immunity_multiplier)))*(A_eq + U_eq + T_eq + D_eq)
  mort_eq <- mpl$pctMort*severe_disease

  ## collate into state variables and parameters needed for ODE model
  state <- c(
    S = S_eq,
    T = T_eq,
    D = D_eq,
    A = A_eq,
    U = U_eq,
    P = P_eq,
    ICA = ICA_eq,
    IB = IB_eq,
    ID = ID_eq,
    IVA = IVA_eq
  )

  params <- list(
    d1 = mpl$d1,
    fd = fd,
    ID0 = mpl$ID0,
    kd = mpl$kD,
    av0 = mpl$av0,
    omega = omega,
    uB = mpl$uB,
    uCA = mpl$uCA,
    ft = ft,
    rP = mpl$rP,
    rU = mpl$rU,
    eta = mpl$eta,
    het_wt = het_wt,
    age_rate = age_rate,
    foi_age = foi_age,
    rel_foi = rel_foi,
    b0 = mpl$b0,
    b1 = mpl$b1,
    IB0 = mpl$IB0,
    kB = mpl$kB,
    phi0 = mpl$phi0,
    phi1 = mpl$phi1,
    kC = mpl$kC,
    rD = mpl$rD,
    recA = mpl$rA,
    rT = mpl$rT,
    dc = mpl$dCA,
    db = mpl$dB,
    dd = mpl$dID,
    dv = mpl$dVA,
    x_I = x_I,
    FOI_eq = FOI_eq,
    na = na,
    ICM=ICM_eq,
    IVM = IVM_eq,
    IC0 = mpl$IC0,
    uD = mpl$uD,
    Iv_eq = Iv_eq,
    Sv_eq = Sv_eq,
    Ev_eq = Ev_eq,
    mv0=mv0,
    Sv_eq=Sv_eq,
    Ev_eq=Ev_eq,
    FOIv_eq=FOIv_eq,
    clin_inc = phi_eq * FOI_eq * c(S_eq+A_eq+U_eq),
    prev_eq = prev_eq,
    mort_eq = mort_eq,
    age_days = age
  )
  return(list(state=state, params = params))
}

#' Generally, pathogen prevalence is a more accesible metric for users, but
#' the Imperial equilibrium function requires an annual EIR.
#' This function converts a given pathogen prevalence to an EIR
#'
#' @param prevalence desired prevalence value
#' @param age_vector age distribution of the population
#' @param ft proportion treated
#' @param params entomological and epidemiological parameters
#' @return a vector of the equilibrium number of humans in each SIS stage
#' @importFrom stats predict
#' @importFrom stats smooth.spline
#' @export
convert_prevalence_to_eir <- function(prevalence, age_vector, ft, params) {
  # first evaluate equilibrium for a range of EIRs
  # pull corresponding prevalence values
  eir_range <- seq(0, 100, by=0.1)
  prevalences <- sapply(eir_range, function(eir) {
    eqm <- equilibrium_Imperial_decoupled_human(age_vector, ft, eir, params)
    return(eqm$params$prev_eq)
  })

  # fit spline to values
  f <- smooth.spline(prevalences, eir_range)

  # return predicted eir value from input prevalence
  eir <- predict(f, prevalence)$y
  return(eir)

}

#' This function calculates the human and mosquito equilibrium values for the decoupled Imperial model.
#' Currently this only works in one node.
#' @param age_vector age structure of population (see vignette for example)
#' @param ft proportion of population seeking treatment
#' @param eir desired annual EIR
#' @param theta parameters
#' @param cube inheritance cube
#' @param spn_P places of the stochastic petri net
#' @return a matrix of the equilibrium number of humans in each Imperial stage by age, and immunity. mosquito equilibrium values, and full theta vector
#' @importFrom statmod gauss.quad.prob
#' @export
equilibrium_Imperial_decoupled <-
  function(age_vector, ft, eir, theta, cube, spn_P) {
    human_eqm <-
      equilibrium_Imperial_decoupled_human(age_vector, ft, eir, theta)

    #  SEI mosquitoes and Imperial humans equilibrium
    #  outputs required parameters in the named list "params"
    #  outputs initial equilibrium for adv users, "init
    #  outputs properly filled initial markings, "M0"

    # augment theta with Imperial params needed to calculate mosquito equilibrium
    # also add derived beta (egg laying rate parameter)
    # calculate the total number of adult mosquitos from eqm distribution
    total_daily_eir <- eir * theta$NH / 365
    lifetime <- human_eqm$params$FOIv_eq * exp(-theta$muF * (1/theta$qEIP)) / (human_eqm$params$FOIv_eq + theta$muF)
    total_M <- total_daily_eir / sum(theta$av0 * theta$Q0 * lifetime)
    theta$FOIv <- human_eqm$params$FOIv_eq
    theta$Iv_eq <- human_eqm$params$Iv_eq
    theta$Sv_eq <- human_eqm$params$Sv_eq
    theta$Ev_eq <- human_eqm$params$Ev_eq
    theta$mv0 <- human_eqm$params$mv0
    theta$total_M <- total_M 

    initialCons <-
      equilibrium_SEI_Imperial(
        params = theta,
        phi = 0.5,
        log_dd = TRUE,
        spn_P = spn_P,
        cube = cube
      )

    # add all imperial params to theta - will be used in ODE simulation
    theta <- c(theta, initialCons$params, human_eqm$params)

    # convert to matrix for use in ODE
    initialCons$H <- matrix(human_eqm$state, ncol = 10)
    colnames(initialCons$H) <-
      c("S", "T", "D", "A", "U", "P", "ICA", "IB", "ID", "IVA")

    # update inheritance cube with transmission blocking parameters
    # augment the cube with RM transmission parameters
    cube$cT <-
      setNames(
        object = rep(x = theta$cT, times = cube$genotypesN),
        nm = cube$genotypesID
      )

    cube$cU <-
      setNames(
        object = rep(x = theta$cU, times = cube$genotypesN),
        nm = cube$genotypesID
      )

    cube$cD <-
      setNames(
        object = rep(x = theta$cD, times = cube$genotypesN),
        nm = cube$genotypesID
      )

    # add b0 to theta, since theta is used in the ODE simulation
    # also add genotype IDs to match b0 values
    theta$b0 <-
      setNames(
        object = rep(x = theta$b0, times = cube$genotypesN),
        nm = cube$genotypesID
      )
    theta$genotypesID <- cube$genotypesID

    return(list(
      theta = theta,
      initialCons = initialCons,
      cube = cube
    ))

  }
