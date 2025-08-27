###############################################################################
#       _   __     __                      __
#      / | / /__  / /__      ______  _____/ /__
#     /  |/ / _ \/ __/ | /| / / __ \/ ___/ //_/
#    / /|  /  __/ /_ | |/ |/ / /_/ / /  / ,<
#   /_/ |_/\___/\__/ |__/|__/\____/_/  /_/|_|
#
#   parameterizeMGDrivE
#   Marshall Lab
#   jared_bennett@berkeley.edu
#   December 2019
#
###############################################################################

#' parameterizeMGDrivE
#'
#' Generate parameters for simulation on a \code{\link{Network}}.
#' Parameters include: average generation time \eqn{g}, population growth rate \eqn{R_{m}},
#' aquatic mortality \eqn{\mu_{Aq}}, and aquatic survival \eqn{\theta_{Aq}}, which
#' are shared between patches and calculated by \code{\link{calcAverageGenerationTime}},
#' \code{\link{calcPopulationGrowthRate}}, and \code{\link{calcLarvalStageMortalityRate}}. \cr
#' Patch-specific parameters \eqn{\alpha} and \eqn{L_{eq}}
#' are calculated for each patch by \code{\link{calcDensityDependentDeathRate}}
#' and \code{\link{calcLarvalPopEquilibrium}}.
#'
#' @param runID Begin counting runs with this set of parameters from this value
#' @param nPatch Number of \code{\link{Patch}}
#' @param simTime Maximum time to run simulation
#' @param sampTime Times to sample, used as tNow %% sampTime, default is every day
#' @param tEgg Length of egg stage
#' @param tLarva Length of larval instar stage
#' @param tPupa Length of pupal stage
#' @param beta Female egg batch size of wild-type
#' @param muAd Wild-type daily adult mortality (1/muAd is average wild-type lifespan)
#' @param popGrowth Daily population growth rate (used to calculate equilibrium)
#' @param AdPopEQ Single number or vector of adult population size at equilibrium
#' (single number implies all patches have the same population)
#' @param LarPopRatio May be empty; if not, a vector gives the wildtype gene frequencies
#' among larval stages at the beginning of simulation or a matrix provides different
#' initial frequencies for each patch (every row is a different patch, must have nrow = nPatch)
#' @param AdPopRatio_F May be empty; if not, a vector gives the wildtype gene frequencies
#' among adult females at the beginning of simulation or a matrix provides different
#' initial frequencies for each patch (every row is a different patch, must have nrow = nPatch)
#' @param AdPopRatio_M May be empty; if not, a vector gives the wildtype gene frequencies
#' among adult males at the beginning of simulation or a matrix provides different
#' initial frequencies for each patch (every row is a different patch, must have nrow = nPatch)
#' @param inheritanceCube Inheritance cube to check/set population ratios at the beginning of the simulation
#'
#' @examples
#' # using default parameters for 2 patches
#' #  using different population sizes for patches
#' simPars <- parameterizeMGDrivE(nPatch = 2, simTime = 365,
#'                                AdPopEQ = c(100,200), inheritanceCube = cubeMendelian())
#'
#' @export
parameterizeMGDrivE <- function(
  runID = 1L,
  nPatch,
  simTime,
  sampTime = 1L,
  tEgg = 1L,
  tLarva = 14L,
  tPupa = 1L,
  beta = 32,
  muAd = 0.123,
  popGrowth = 1.096,
  AdPopEQ,
  LarPopRatio,
  AdPopRatio_F,
  AdPopRatio_M,
  inheritanceCube
){

  # check required parameters
  if(any(missing(nPatch),missing(simTime),missing(AdPopEQ),missing(inheritanceCube))){
    stop("nPatch, simTime, AdPopEQ, and inheritanceCube must be provided by the user.")
  }

  # make empty parameter list
  pars = list()

  # fill list
  pars$nPatch = nPatch
  pars$simTime = simTime
  pars$sampTime = sampTime
  pars$runID = runID

  # biological parameters
  pars$timeAq = c("E"=tEgg, "L"=tLarva, "P"=tPupa)
  pars$beta = beta

  # initial parameters
  pars$muAd = muAd
  pars$dayPopGrowth = popGrowth

  if(length(AdPopEQ) == 1){
    AdPopEQ = rep.int(x = AdPopEQ, times = nPatch)
  } else if(length(AdPopEQ)!=nPatch){
    stop("length of AdPopEQ vector must be 1 or nPatch (number of patches)")
  }
  pars$AdPopEQ = AdPopEQ


  # setup larval initial pop ratio
  if(missing(LarPopRatio)){
    # default behaviour - this way nothing needs to be specified
    # setup vector of the correct length with all weights equal
    larRatio <- rep.int(x = 1, times = length(inheritanceCube$wildType))
    # normalize to 1
    larRatio <- larRatio/sum(larRatio)
    # set all patches to the same thing
    pars$LarPopRatio <- matrix(data = larRatio, nrow = nPatch, ncol = length(inheritanceCube$wildType),
                                 byrow = TRUE, dimnames = list(NULL,inheritanceCube$wildType))

  } else if(is.null(dim(LarPopRatio)) ){
    # behaviour of user supplied 1 patch worth of weights.
    # has to be supplied as a vector, the matrix stuff makes this way too difficult
    if(abs(sum(LarPopRatio) - 1) > sqrt(.Machine$double.eps) ) stop('LarPopRatio must sum to 1')
    # check that columns names are in inheritance cube
    if(is.null(names(LarPopRatio)) || !all(names(LarPopRatio) %in% inheritanceCube$genotypesID)) {
      stop("Names for LarPopRatio must be specified as one of the genotypesID in the inheritance cube.")
    }
    # set all patches equal
    pars$LarPopRatio <- matrix(data = LarPopRatio, nrow = nPatch, ncol = length(LarPopRatio),
                                 byrow = TRUE, dimnames = list(NULL,names(LarPopRatio)) )

  } else if(dim(LarPopRatio)[1] == nPatch){
    # behaviour if user supplies a matrix of probabilities.
    # each row is a different patch
    # check that all patches sum to 1
    if(any(abs(rowSums(LarPopRatio) - 1) > sqrt(.Machine$double.eps)) ) stop('Each row of LarPopRatio must sum to 1')
    # check that columns names are in inheritance cube
    if(is.null(colnames(LarPopRatio)) || !all(colnames(LarPopRatio) %in% inheritanceCube$genotypesID)) {
      stop("Column names for LarPopRatio must be specified as one of the genotypesID in the inheritance cube.")
    }
    # store
    pars$LarPopRatio <- LarPopRatio

  } else {
    stop("LarPopRatio has been miss specified.\n
         Left blank - default, all populations are the same and begin as wild-type individuals\n
         Vector - a named vector that sums to one. All populations will be the same\n
         Matrix - an nPatch by nGenotype matrix with column names and all rows sum to 1.
         Specifies each population individually.")
  }


  # setup female initial pop ratio
  if(missing(AdPopRatio_F)){
    # default behaviour - this way nothing needs to be specified

    # now, two options. This is not sex based, so only 1 wild-type
    #  or, sex based cube, we have an X and a Y, females are not Y
    if(length(inheritanceCube$wildType) == 1){
      # only 1 wild-type, so no problem
      pars$AdPopRatio_F <- matrix(data = 1, nrow = nPatch, ncol = length(inheritanceCube$wildType),
                                  dimnames = list(NULL,inheritanceCube$wildType))
    } else if(length(inheritanceCube$wildType) == 2){
      # one female and one male wild-type.
      # get index of female wild-type
      whichGeno <- grep(pattern = "Y", x = inheritanceCube$wildType, fixed = TRUE, invert = TRUE)
      # setup default matrix
      pars$AdPopRatio_F <- matrix(data = 0, nrow = nPatch, ncol = length(inheritanceCube$wildType),
                                  dimnames = list(NULL,inheritanceCube$wildType))
      # set all to female genotype
      pars$AdPopRatio_F[ ,whichGeno] <- 1

    } else {
      stop("Default AdPopRatio_F only handles 1 or 2 wild-type genotypes.\n
           Please provide a matrix specifying female initial genotypes for every patch.")
    }

  } else if(is.null(dim(AdPopRatio_F)) ){
    # behaviour of user supplied 1 patch worth of weights.
    # has to be supplied as a vector, the matrix stuff makes this way too difficult
    if(abs(sum(AdPopRatio_F) - 1) > sqrt(.Machine$double.eps) ) stop('AdPopRatio_F must sum to 1')
    # check that columns names are in inheritance cube
    if(is.null(names(AdPopRatio_F)) || !all(names(AdPopRatio_F) %in% inheritanceCube$genotypesID)) {
      stop("Names for AdPopRatio_F must be specified as one of the genotypesID in the inheritance cube.")
    }
    # set all patches equal
    pars$AdPopRatio_F <- matrix(data = AdPopRatio_F, nrow = nPatch, ncol = length(AdPopRatio_F),
                               byrow = TRUE, dimnames = list(NULL,names(AdPopRatio_F)) )

  } else if(dim(AdPopRatio_F)[1] == nPatch){
    # behaviour if user supplies a matrix of probabilities.
    # each row is a different patch
    # check that all patches sum to 1
    if(any(abs(rowSums(AdPopRatio_F) - 1) > sqrt(.Machine$double.eps)) ) stop('Each row of AdPopRatio_F must sum to 1')
    # check that columns names are in inheritance cube
    if(is.null(colnames(AdPopRatio_F)) || !all(colnames(AdPopRatio_F) %in% inheritanceCube$genotypesID)) {
      stop("Column names for AdPopRatio_F must be specified as one of the genotypesID in the inheritance cube.")
    }
    # store
    pars$AdPopRatio_F <- AdPopRatio_F

  } else {
    stop("AdPopRatio_F has been miss specified.\n
         Left blank - default, all populations are the same and begin as wild-type individuals\n
         Vector - a named vector that sums to one. All populations will be the same\n
         Matrix - an nPatch by nGenotype matrix with column names and all rows sum to 1.
         Specifies each population individually.")
  }


  # setup male initial pop ratio
  if(missing(AdPopRatio_M)){
    # default behaviour - this way nothing needs to be specified

    # now, two options. This is not sex based, so only 1 wild-type
    #  or, sex based cube, we have an X and a Y, females are not Y
    if(length(inheritanceCube$wildType) == 1){
      # only 1 wild-type, so no problem
      pars$AdPopRatio_M <- matrix(data = 1, nrow = nPatch, ncol = length(inheritanceCube$wildType),
                                  dimnames = list(NULL,inheritanceCube$wildType))
    } else if(length(inheritanceCube$wildType) == 2){
      # one female and one male wild-type.
      # get index of female wild-type
      whichGeno <- grep(pattern = "Y", x = inheritanceCube$wildType, fixed = TRUE)
      # setup default matrix
      pars$AdPopRatio_M <- matrix(data = 0, nrow = nPatch, ncol = length(inheritanceCube$wildType),
                                  dimnames = list(NULL,inheritanceCube$wildType))
      # set all to female genotype
      pars$AdPopRatio_M[ ,whichGeno] <- 1

    } else {
      stop("Default AdPopRatio_M only handles 1 or 2 wild-type genotypes.\n
           Please provide a matrix specifying female initial genotypes for every patch.")
    }

  } else if(is.null(dim(AdPopRatio_M)) ){
    # behaviour of user supplied 1 patch worth of weights.
    # has to be supplied as a vector, the matrix stuff makes this way too difficult
    if(abs(sum(AdPopRatio_M) - 1) > sqrt(.Machine$double.eps) ) stop('AdPopRatio_M must sum to 1')
    # check that columns names are in inheritance cube
    if(is.null(names(AdPopRatio_M)) || !all(names(AdPopRatio_M) %in% inheritanceCube$genotypesID)) {
      stop("Names for AdPopRatio_M must be specified as one of the genotypesID in the inheritance cube.")
    }
    # set all patches equal
    pars$AdPopRatio_M <- matrix(data = AdPopRatio_M, nrow = nPatch, ncol = length(AdPopRatio_M),
                                byrow = TRUE, dimnames = list(NULL,names(AdPopRatio_M)) )

  } else if(dim(AdPopRatio_M)[1] == nPatch){
    # behaviour if user supplies a matrix of probabilities.
    # each row is a different patch
    # check that all patches sum to 1
    if(any(abs(rowSums(AdPopRatio_M) - 1) > sqrt(.Machine$double.eps)) ) stop('Each row of AdPopRatio_M must sum to 1')
    # check that columns names are in inheritance cube
    if(is.null(colnames(AdPopRatio_M)) || !all(colnames(AdPopRatio_M) %in% inheritanceCube$genotypesID)) {
      stop("Column names for AdPopRatio_M must be specified as one of the genotypesID in the inheritance cube.")
    }
    # store
    pars$AdPopRatio_M <- AdPopRatio_M

  } else {
    stop("AdPopRatio_M has been miss specified.\n
         Left blank - default, all populations are the same and begin as wild-type individuals\n
         Vector - a named vector that sums to one. All populations will be the same\n
         Matrix - an nPatch by nGenotype matrix with column names and all rows sum to 1.
         Specifies each population individually.")
  }


  # derived parameters
  pars$g = calcAverageGenerationTime(pars$timeAq,muAd)
  pars$genPopGrowth = calcPopulationGrowthRate(popGrowth,pars$g)
  pars$muAq = calcLarvalStageMortalityRate(pars$genPopGrowth,muAd,beta,pars$timeAq)
  pars$thetaAq = c("E"=calcAquaticStageSurvivalProbability(pars$muAq,tEgg),
                   "L"=calcAquaticStageSurvivalProbability(pars$muAq,tLarva),
                   "P"=calcAquaticStageSurvivalProbability(pars$muAq,tPupa))


  # patch-specific derived parameters
  pars$alpha = calcDensityDependentDeathRate(beta, pars$thetaAq, pars$timeAq,
                                             AdPopEQ, pars$genPopGrowth)
  pars$Leq = calcLarvalPopEquilibrium(pars$alpha,pars$genPopGrowth)


  # check the list
  invisible(Map(f = check, pars))


  # if pass the check, return the parameter vector
  return(pars)
}

########################################################################
# Equations and Equilibrium Parameters for parameterizeMGDrivE()
########################################################################

# check for positive parameter values
check <- function(x){
  if(is.numeric(x)||is.integer(x)){
    if(any(x < 0)){
      stop("only nonnegative parameter values allowed")
    }
  }
}

#' Calculate Density-dependent Larval Mortality
#'
#' Calculate \eqn{\alpha}, the strength of density-dependent mortality during the
#' larval stage, given by: \deqn{\alpha=\Bigg( \frac{1/2 * \beta * \theta_e * Ad_{eq}}{R_m-1} \Bigg) * \Bigg( \frac{1-(\theta_l / R_m)}{1-(\theta_l / R_m)^{1/T_l}} \Bigg)}
#'
#' @param fertility Number of eggs per oviposition for wild-type females, \eqn{\beta}
#' @param thetaAq Vector of density-independent survival probabilities of aquatic stages, \eqn{\theta_{e}, \theta_{l}}
#' @param tAq Vector of lengths of aquatic stages, \eqn{T_{e}, T_{l}, T_{p}}
#' @param adultPopSizeEquilibrium Adult population size at equilibrium, \eqn{Ad_{eq}}
#' @param populationGrowthRate Population growth in absence of density-dependent mortality \eqn{R_{m}}
#'
#' @keywords internal
#'
calcDensityDependentDeathRate <- function(fertility, thetaAq, tAq,
                                          adultPopSizeEquilibrium, populationGrowthRate){

  prodA = (fertility * thetaAq[["E"]] * (adultPopSizeEquilibrium/2)) / (populationGrowthRate-1)
  prodB_numerator = (1 - (thetaAq[["L"]] / populationGrowthRate))
  prodB_denominator = (1 - ((thetaAq[["L"]]/populationGrowthRate)^(1/tAq[["L"]])))

  return(prodA*(prodB_numerator/prodB_denominator))
}

#' Calculate Average Generation Time
#'
#' Calculate \eqn{g}, average generation time, given by: \deqn{g=T_e+T_l+T_p+\frac{1}{\mu_{ad}}}
#'
#' @param stagesDuration Vector of lengths of aquatic stages, \eqn{T_{e}, T_{l}, T_{p}}
#' @param adultMortality Adult mortality rate, \eqn{\mu_{ad}}
#'
#' @keywords internal
#'
calcAverageGenerationTime <- function(stagesDuration, adultMortality){
  return(sum(stagesDuration) + (1.0 / adultMortality))
}

#' Calculate Generational Population Growth Rate
#'
#' Calculate \eqn{R_{m}}, population growth in absence of density-dependent mortality,
#' given by: \deqn{(r_{m})^{g}}
#'
#' @param dailyPopGrowthRate Daily population growth rate, \eqn{r_{m}}
#' @param averageGenerationTime See \code{\link{calcAverageGenerationTime}}
#'
#' @keywords internal
#'
calcPopulationGrowthRate <- function(dailyPopGrowthRate, averageGenerationTime){
  return(dailyPopGrowthRate^averageGenerationTime)
}

#' Calculate Aquatic Stage Survival Probability
#'
#' Calculate \eqn{\theta_{st}}, density-independent survival probability, given
#' by: \deqn{\theta_{st}=(1-\mu_{st})^{T_{st}}}
#'
#' @param mortalityRate Daily mortality probability, \eqn{\mu_{st}}
#' @param stageDuration Duration of aquatic stage, \eqn{T^{st}}
#'
#' @keywords internal
#'
calcAquaticStageSurvivalProbability <- function(mortalityRate, stageDuration){
  return((1-mortalityRate)^stageDuration)
}

#' Calculate Larval Stage Mortality Rate
#'
#' Calculate \eqn{\mu_{l}}, the larval mortality, given by
#' \deqn{\mu_l=1-\Bigg( \frac{R_m * \mu_{ad}}{1/2 * \beta * (1-\mu_m)} \Bigg)^{\frac{1}{T_e+T_l+T_p}}}
#'
#' @param generationPopGrowthRate See \code{\link{calcPopulationGrowthRate}}
#' @param adultMortality Adult mortality rate, \eqn{\mu_{ad}}
#' @param fertility Number of eggs per oviposition for wild-type females, \eqn{\beta}
#' @param aquaticStagesDuration Vector of lengths of aquatic stages, \eqn{T_{e}, T_{l}, T_{p}}
#'
#' @keywords internal
#'
calcLarvalStageMortalityRate <- function(generationPopGrowthRate, adultMortality,
                                         fertility, aquaticStagesDuration){

  a = generationPopGrowthRate*adultMortality
  b = (fertility/2)*(1-adultMortality)
  c = sum(aquaticStagesDuration)

  return(1-(a/b)^(1/c))
}

#' Calculate Equilibrium Larval Population
#'
#' Equilibrium larval population size to sustain adult population.
#'
#' @param alpha See \code{\link{calcDensityDependentDeathRate}}
#' @param Rm See \code{\link{calcPopulationGrowthRate}}
#'
#' @keywords internal
#'
calcLarvalPopEquilibrium <- function(alpha, Rm){
  return(as.integer(round(alpha * (Rm-1))))
}
