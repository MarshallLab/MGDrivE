########################################################################
#       _   __     __                      __
#      / | / /__  / /__      ______  _____/ /__
#     /  |/ / _ \/ __/ | /| / / __ \/ ___/ //_/
#    / /|  /  __/ /_ | |/ |/ / /_/ / /  / ,<
#   /_/ |_/\___/\__/ |__/|__/\____/_/  /_/|_|
#
#   parameterizeMGDrivE
#   Marshall Lab
#   November 2017
#
########################################################################

#' parameterizeMGDrivE
#'
#' Generate parameters for simulation on a \code{\link{Network}}.
#' Parameters average generation time \eqn{g}, population growth rate \eqn{R_{m}},
#' aquatic mortality \eqn{\mu_{Aq}}, and aquatic survival \eqn{\theta_{Aq}}
#' are shared between patches and calculated by \code{\link{calcAverageGenerationTime}},
#' \code{\link{calcPopulationGrowthRate}}, \code{\link{calcLarvalStageMortalityRate}},
#' and \code{\link{calcAquaticStagesSurvivalProbability}}. \cr
#' Patch-specific parameters \eqn{\alpha} and \eqn{L_{eq}}
#' are calculated for each patch by \code{\link{calcDensityDependentDeathRate}}
#' and \code{\link{calcLarvalPopEquilibrium}}.
#'
#' @param runID begin counting runs with this set of parameters from this value
#' @param nPatch number of \code{\link{Patch}}
#' @param simTime maximum time to run simulation
#' @param parallel append process id (see \code{link[base]{Sys.getpid}}) to output files for running in parallel
#' @param moveVar variance of stochastic movement (not used in diffusion model of migration). It affects the concentration of probability in the Dirchlet simplex, small values lead to high variance and large values lead to low variance.
#' @param tEgg length of egg stage
#' @param tLarva length of larval instar stage
#' @param tPupa length of pupal stage
#' @param beta female egg batch size of wild-type
#' @param muAd wild-type daily adult mortality (1/muAd is average wild-type lifespan)
#' @param popGrowth daily population growth rate (used to calculate equilibrium)
#' @param AdPopEQ vector of adult population size at equilibrium
#' @param LarPopRatio may be \code{NULL}; if not, gives the wildtype gene frequencies among larval stages at the beginning of simulation
#' @param AdPopRatio_F may be \code{NULL}; if not, gives the wildtype gene frequencies among adult females at the beginning of simulation
#' @param AdPopRatio_M may be \code{NULL}; if not, gives the wildtype gene frequencies among adult males at the beginning of simulation
#'
#' @examples
#' # using default parameters for 2 patches
#' #  using different population sizes for patches
#' simPars <- parameterizeMGDrivE(nPatch = 2, simTime = 365,
#'                                AdPopEQ = c(100,200))
#'
#'
#' @export
parameterizeMGDrivE <- function(
  runID = 1L,
  nPatch,
  simTime,
  parallel = FALSE,
  moveVar = 1000,
  tEgg = 1L,
  tLarva = 14L,
  tPupa = 1L,
  beta = 32,
  muAd = 0.123,
  popGrowth = 1.096,
  AdPopEQ,
  LarPopRatio = NULL,
  AdPopRatio_F = NULL,
  AdPopRatio_M = NULL
){

  # make empty parameter list
  pars = list()

  # fill list
  pars$nPatch = nPatch
  pars$simTime = simTime
  pars$parallel = parallel
  pars$moveVar = moveVar
  pars$runID = runID

  # biological parameters
  pars$timeAq = initStagesDurations(egg = tEgg,larva = tLarva,pupa = tPupa)
  pars$timeE = pars$timeAq[1]
  pars$timeL = pars$timeAq[2]
  pars$timeP = pars$timeAq[3]
  pars$windowSize = calcMemoryWindow(pars$timeAq)
  pars$beta = beta

  # initial parameters
  pars$muAd = muAd
  pars$dayPopGrowth = popGrowth

  if(length(AdPopEQ)!=nPatch){
    stop("length of AdPopEQ vector must equal nPatch (number of patches)")
  }
  pars$AdPopEQ = AdPopEQ

  # if not \code{NULL} add initial pop ratios
  if(!is.null(LarPopRatio)){
    if(sum(LarPopRatio)!=1){
      stop("if 'LarPopRatio' specified, it must sum to 1!")
    }
    pars$LarPopRatio = LarPopRatio
  } else {
    pars$LarPopRatio <- 1.0
  }

  if(!is.null(AdPopRatio_F)){
    if(sum(AdPopRatio_F)!=1){
      stop("if 'AdPopRatio_F' specified, it must sum to 1!")
    }
    pars$AdPopRatio_F = AdPopRatio_F
  } else {
    pars$AdPopRatio_F <- 1.0
  }

  if(!is.null(AdPopRatio_M)){
    if(sum(AdPopRatio_M)!=1){
      stop("if 'AdPopRatio_M' specified, it must sum to 1!")
    }
    pars$AdPopRatio_M = AdPopRatio_M
  } else {
    pars$AdPopRatio_M <- 1.0
  }

  # derived parameters
  pars$g = calcAverageGenerationTime(pars$timeAq,muAd)
  pars$genPopGrowth = calcPopulationGrowthRate(popGrowth,pars$g)
  pars$muAq = calcLarvalStageMortalityRate(pars$genPopGrowth,muAd,beta,pars$timeAq)
  pars$thetaAq = calcAquaticStagesSurvivalProbability(
                    calcAquaticStageSurvivalProbability(pars$muAq,tEgg),
                    calcAquaticStageSurvivalProbability(pars$muAq,tLarva),
                    calcAquaticStageSurvivalProbability(pars$muAq,tPupa)
  )
  pars$thetaE = pars$thetaAq[1]
  pars$thetaL = pars$thetaAq[2]
  pars$thetaP = pars$thetaAq[3]

  # patch-specific derived parameters
  pars$alpha = rep(0,nPatch)
  pars$Leq = rep(0,nPatch)
  for(i in 1:nPatch){
    pars$alpha[i] = calcDensityDependentDeathRate(beta,pars$thetaAq,
                                                  pars$timeAq,
                                                  AdPopEQ[i],
                                                  pars$genPopGrowth)
    pars$Leq[i] = calcLarvalPopEquilibrium(pars$alpha[i],pars$genPopGrowth)
  }

  # check for positive parameter values
  check <- function(x){
    if(is.numeric(x)||is.integer(x)){
      if(length(x)>1){
        if(any(x < 0)){
          stop("only nonnegative parameter values allowed")
        }
      } else {
        if(x < 0){
          stop("only nonnegative parameter values allowed")
        }
      }
    }
  }
  # check the list
  invisible(Map(
    f = check, pars
  ))

  # if pass the check, return the parameter vector
  return(pars)
}

########################################################################
# Equations and Equilibrium Parameters for parameterizeMGDrivE()
########################################################################

#' Initialize Aquatic Stages Durations
#'
#' Initialises the vector that holds the duration of each aquatic stage
#'
#' @param egg length of egg stage (days)
#' @param larva length of larval stage (days)
#' @param pupa length of pupal stage (days)
#'
initStagesDurations <- function(egg=1,larva=14,pupa=1){
  out = c(egg,larva,pupa)
  names(out) = c("E","L","P")
  return(out)
}

#' Calculate Memory Window
#'
#' Calculates the necessary window of population history required for the model to work
#'
#' @param stagesDuration vector of aquatic stages length (egg/larva/pupa in days)
#'
calcMemoryWindow <- function(stagesDuration){
  return(sum(stagesDuration))
}

#' Calculate Density-dependent Larval Mortality
#'
#' Calculate \eqn{\alpha}, the strength of density-dependent mortality during the larval stage, given by: \deqn{\alpha=\Bigg( \frac{1/2 * \beta * \theta_e * Ad_{eq}}{R_m-1} \Bigg) * \Bigg( \frac{1-(\theta_l / R_m)}{1-(\theta_l / R_m)^{1/T_l}} \Bigg)}
#'
#' @param fertility number of eggs per oviposition for wild-type females, \eqn{\beta}
#' @param thetaAq vector of density-independent survival probabilities of aquatic stages, \eqn{\theta_{e}, \theta_{l}}
#' @param tAq vector of lengths of aquatic stages, \eqn{T_{e}, T_{l}, T_{p}}
#' @param adultPopSizeEquilibrium adult population size at equilbrium, \eqn{Ad_{eq}}
#' @param populationGrowthRate population growth in absence of density-dependent mortality \eqn{R_{m}}
#'
calcDensityDependentDeathRate <- function(fertility, thetaAq, tAq, adultPopSizeEquilibrium, populationGrowthRate){
    prodA = (fertility * thetaAq[["E"]] * (adultPopSizeEquilibrium/2)) / (populationGrowthRate-1)
    prodB_numerator = (1 - (thetaAq[["L"]] / populationGrowthRate))
    prodB_denominator = (1 - ((thetaAq[["L"]]/populationGrowthRate)^(1/tAq[["L"]])))
    return(prodA*(prodB_numerator/prodB_denominator))
}

#' Calculate Average Generation Time
#'
#' Calculate \eqn{g}, average generation time, given by: \deqn{g=T_e+T_l+T_p+\frac{1}{\mu_{ad}}}
#'
#' @param stagesDuration vector of lengths of aquatic stages, \eqn{T_{e}, T_{l}, T_{p}}
#' @param adultMortality adult mortality rate, \eqn{\mu_{ad}}
#'
calcAverageGenerationTime <- function(stagesDuration, adultMortality){
  return(sum(stagesDuration) + (1.0 / adultMortality))
}

#' Calculate Generational Population Growth Rate
#'
#' Calculate \eqn{R_{m}}, population growth in absence of density-dependent mortality, given by: \deqn{(r_{m})^{g}}
#'
#' @param dailyPopGrowthRate daily population growth rate, \eqn{r_{m}}
#' @param averageGenerationTime see \code{\link{calcAverageGenerationTime}}
#'
calcPopulationGrowthRate <- function(dailyPopGrowthRate, averageGenerationTime){
  return(dailyPopGrowthRate^averageGenerationTime)
}

#' Calculate Aquatic Stage Surival Probability
#'
#' Calculate \eqn{\theta_{st}}, density-independent survival probability, given by: \deqn{\theta_{st}=(1-\mu_{st})^{T_{st}}}
#'
#' @param mortalityRate daily mortality probability, \eqn{\mu_{st}}
#' @param stageDuration duration of aquatic stage, \eqn{T^{st}}
#'
calcAquaticStageSurvivalProbability <- function(mortalityRate, stageDuration){
  # return(exp(-mortalityRate)^stageDuration)
  return((1-mortalityRate)^stageDuration)
}

#' Calculate Larval Stage Mortality Rate
#'
#' Calculate \eqn{\mu_{l}}, the larval mortality, given by \deqn{\mu_l=1-\Bigg( \frac{R_m * \mu_{ad}}{1/2 * \beta * (1-\mu_m)} \Bigg)^{\frac{1}{T_e+T_l+T_p}}}
#'
#' @param generationPopGrowthRate see \code{\link{calcPopulationGrowthRate}}
#' @param adultMortality adult mortality rate, \eqn{\mu_{ad}}
#' @param fertility number of eggs per oviposition for wild-type females, \eqn{\beta}
#' @param aquaticStagesDuration vector of lengths of aquatic stages, \eqn{T_{e}, T_{l}, T_{p}}
#'
calcLarvalStageMortalityRate <- function(generationPopGrowthRate, adultMortality, fertility, aquaticStagesDuration){
    a = generationPopGrowthRate*adultMortality
    b = (fertility/2)*(1-adultMortality)
    c = sum(aquaticStagesDuration)
    return(1-(a/b)^(1/c))
}

#' Calculate Survival Probability of entire Aquatic Stage Life-cycle
#'
#' Calculate vector of survival probabilities for each stage of aquatic lifecycle.
#'
#' @param eggSurvivalProbability see \code{\link{calcAquaticStageSurvivalProbability}}
#' @param larvaSurvivalProbability see \code{\link{calcAquaticStageSurvivalProbability}}
#' @param pupaSurvivalProbability see \code{\link{calcAquaticStageSurvivalProbability}}
#'
calcAquaticStagesSurvivalProbability <- function(eggSurvivalProbability, larvaSurvivalProbability, pupaSurvivalProbability){
  out = c(eggSurvivalProbability,larvaSurvivalProbability,pupaSurvivalProbability)
  names(out) = c("E","L","P")
  return(out)
}

#' Calculate Equilibrium Larval Population
#'
#' Equilibrium larval population to sustain population.
#'
#' @param alpha see \code{\link{calcDensityDependentDeathRate}}
#' @param Rm see \code{\link{calcPopulationGrowthRate}}
#'
calcLarvalPopEquilibrium <- function(alpha,Rm){
  return(as.integer(round(alpha * (Rm-1))))
}
