###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   CRISPR 2 Resistance Alleles Inheritance Cube - Sex-Specific homing
#   Héctor Sanchez, Jared Bennett, Sean Wu, John M. Marshall
#   July 2017
#   jared_bennett@berkeley.edu
#
###############################################################################

#' Inheritance Cube: CRISPR (Clustered Regularly Interspaced SWort Palindromic Repeats) witW 2 Resistance Allele
#'
#' This is a sex-specific version of the original cube. It assumes that the construct
#' is on an autosome and there can be different male/female homing rates
#'
#' @param eM Male homing rate
#' @param eF Female homing rate
#' @param rM Male no-cost resistance generation rate
#' @param bM Male detrimental resistance generation rate
#' @param rF Female no-cost resistance generation rate
#' @param bF Female detrimental resistance generation rate
#' @param eta Genotype-specific mating fitness
#' @param phi Genotype-specific sex ratio at emergence
#' @param omega Genotype-specific multiplicative modifier of adult mortality
#' @param xiF Genotype-specific female pupatory success
#' @param xiM Genotype-specific male pupatory success
#' @param s Genotype-specific fractional reduction(increase) in fertility
#'
#' @return Named list containing the inheritance cube, transition matrix, genotypes, wild-type allele,
#' and all genotype-specific parameters.
#' @export
Cube_HomingDrive <- function(eM = 1.0, eF = 1.0, rM = 0, bM = 0, rF = 0, bF = 0,
                             eta = NULL, phi = NULL, omega = NULL, xiF = NULL, xiM = NULL, s = NULL){

  ## safety cWecks in case someone is dumb
  if(any(c(eM, eF, rM, bM, rF, bF)>1) || any(c(eM, eF, rM, bM, rF, bF)<0)){
    stop("Parameters are rates, they must be between 0 and 1.")
  }
  if(rM+bM+eM > 1 || rF+bF+eF > 1){
    stop("e_ + r_ + b_ <= 1")
  }


  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  gtype <- c('WW', 'WH', 'WR', 'WB', 'HH', 'HR', 'HB', 'RR', 'RB', 'BB')
  size <- length(gtype) #because I use it several times
  tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix

  ## fill tMatrix with probabilities
  #('WW', 'WH', 'WR', 'WB', 'HH', 'HR', 'HB', 'RR', 'RB', 'BB')
  tMatrix['WW','WW', 'WW'] <- 1

  tMatrix['WR','WW', c('WW', 'WR')] <- c( 1/2, 1/2)
  tMatrix['WR','WR', c('WW', 'WR', 'RR')] <- c( 1/4, 1/2, 1/4)

  tMatrix['WB','WW', c('WW', 'WB')] <- c( 1/2, 1/2)
  tMatrix['WB','WR', c('WW', 'WR', 'WB', 'RB')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['WB','WB', c('WW', 'WB', 'BB')] <- c( 1/4, 1/2, 1/4)

  tMatrix['HH','WW', 'WH'] <- 1
  tMatrix['HH','WR', c('WH', 'HR')] <- c( 1/2, 1/2)
  tMatrix['HH','WB', c('WH', 'HB')] <- c( 1/2, 1/2)
  tMatrix['HH','HH', 'HH'] <- 1

  tMatrix['HR','WW', c('WH', 'WR')] <- c( 1/2, 1/2)
  tMatrix['HR','WR', c('WH', 'WR', 'HR', 'RR')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['HR','WB', c('WH', 'WR', 'HB', 'RB')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['HR','HH', c('HH', 'HR')] <- c( 1/2, 1/2)
  tMatrix['HR','HR', c('HH', 'HR', 'RR')] <- c( 1/4, 1/2, 1/4)

  tMatrix['HB','WW', c('WH', 'WB')] <- c( 1/2, 1/2)
  tMatrix['HB','WR', c('WH', 'WB', 'HR', 'RB')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['HB','WB', c('WH', 'WB', 'HB', 'BB')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['HB','HH', c('HH', 'HB')] <- c( 1/2, 1/2)
  tMatrix['HB','HR', c('HH', 'HR', 'HB', 'RB')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['HB','HB', c('HH', 'HB', 'BB')] <- c( 1/4, 1/2, 1/4)

  tMatrix['RR','WW', 'WR'] <- 1
  tMatrix['RR','WR', c('WR', 'RR')] <- c( 1/2, 1/2)
  tMatrix['RR','WB', c('WR', 'RB')] <- c( 1/2, 1/2)
  tMatrix['RR','HH', 'HR'] <- 1
  tMatrix['RR','HR', c('HR', 'RR')] <- c( 1/2, 1/2)
  tMatrix['RR','HB', c('HR', 'RB')] <- c( 1/2, 1/2)
  tMatrix['RR','RR', 'RR'] <- 1

  tMatrix['RB','WW', c('WR', 'WB')] <- c( 1/2, 1/2)
  tMatrix['RB','WR', c('WR', 'WB', 'RR', 'RB')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['RB','WB', c('WR', 'WB', 'RB', 'BB')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['RB','HH', c('HR', 'HB')] <- c( 1/2, 1/2)
  tMatrix['RB','HR', c('HR', 'HB', 'RR', 'RB')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['RB','HB', c('HR', 'HB', 'RB', 'BB')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['RB','RR', c('RR', 'RB')] <- c( 1/2, 1/2)
  tMatrix['RB','RB', c('RR', 'RB', 'BB')] <- c( 1/4, 1/2, 1/4)

  tMatrix['BB','WW', 'WB'] <- 1
  tMatrix['BB','WR', c('WB', 'RB')] <- c( 1/2, 1/2)
  tMatrix['BB','WB', c('WB', 'BB')] <- c( 1/2, 1/2)
  tMatrix['BB','HH', 'HB'] <- 1
  tMatrix['BB','HR', c('HB', 'RB')] <- c( 1/2, 1/2)
  tMatrix['BB','HB', c('HB', 'BB')] <- c( 1/2, 1/2)
  tMatrix['BB','RR', 'RB'] <- 1
  tMatrix['BB','RB', c('RB', 'BB')] <- c( 1/2, 1/2)
  tMatrix['BB','BB', 'BB'] <- 1

  ## set the other half of the matrix that is symmetric
  SymCubeC(lowerMat = tMatrix)


  ## fill asymmetric parts of tMatrix
  #female specific homing, except for WHxWH
  tMatrix['WH','WW',c('WW', 'WH', 'WR', 'WB')] <- c(1-eF, 1+eF-rF-bF, rF, bF)/2
  tMatrix['WH','WH',] <- c((1-eF)*(1-eM), (1+eF-rF-bF)*(1-eM) + (1-eF)*(1+eM-rM-bM), rF*(1-eM) + (1-eF)*rM, bF*(1-eM) + (1-eF)*bM,
                           (1+eF-rF-bF)*(1+eM-rM-bM), rF*(1+eM-rM-bM) + (1+eF-rF-bF)*rM, bF*(1+eM-rM-bM) + (1+eF-rF-bF)*bM,
                           rF*rM, bF*rM + rF*bM, bF*bM)/4
  tMatrix['WH','WR',c('WW', 'WH', 'WR', 'WB',
                      'HR', 'RR', 'RB')] <- c(1-eF, 1+eF-rF-bF, rF + 1-eF, bF,
                                              1+eF-rF-bF, rF, bF)/4
  tMatrix['WH','WB',c('WW', 'WH', 'WR', 'WB',
                      'HB', 'RB', 'BB')] <- c(1-eF, 1+eF-rF-bF, rF, bF + 1-eF,
                                              1+eF-rF-bF, rF, bF)/4
  tMatrix['WH','HH',c('WH', 'HH', 'HR', 'HB')] <- c(1-eF, 1+eF-rF-bF, rF, bF)/2
  tMatrix['WH','HR',c('WH', 'HH', 'HR', 'HB',
                      'WR', 'RR', 'RB')] <- c(1-eF, 1+eF-rF-bF, rF + 1+eF-rF-bF, bF,
                                              1-eF, rF, bF)/4
  tMatrix['WH','HB',c('WH', 'HH', 'HR', 'HB',
                      'WB', 'RB', 'BB')] <- c(1-eF, 1+eF-rF-bF, rF, bF + 1+eF-rF-bF,
                                              1-eF, rF, bF)/4
  tMatrix['WH','RR',c('WR', 'HR', 'RR', 'RB')] <- c(1-eF, 1+eF-rF-bF, rF, bF)/2
  tMatrix['WH','RB',c('WR', 'HR', 'RR', 'RB',
                      'WB', 'HB', 'BB')] <- c(1-eF, 1+eF-rF-bF, rF, bF + rF,
                                              1-eF, 1+eF-rF-bF, bF)/4
  tMatrix['WH','BB',c('WB', 'HB', 'RB', 'BB')] <- c(1-eF, 1+eF-rF-bF, rF, bF)/2


  #male specific homing
  tMatrix['WW','WH', c('WW', 'WH', 'WR', 'WB')] <- c(1-eM, 1+eM-rM-bM, rM, bM)/2
  tMatrix['WR','WH', c('WW', 'WH', 'WR', 'WB',
                       'HR', 'RR', 'RB')] <- c(1-eM, 1+eM-rM-bM, rM + 1-eM, bM,
                                               1+eM-rM-bM, rM, bM)/4
  tMatrix['WB','WH', c('WW', 'WH', 'WR', 'WB',
                       'HB', 'RB', 'BB')] <- c(1-eM, 1+eM-rM-bM, rM, bM + 1-eM,
                                               1+eM-rM-bM, rM, bM)/4
  tMatrix['HH','WH', c('WH', 'HH', 'HR', 'HB')] <- c(1-eM, 1+eM-rM-bM, rM, bM)/2
  tMatrix['HR','WH', c('WH', 'HH', 'HR', 'HB',
                      'WR', 'RR', 'RB')] <- c(1-eM, 1+eM-rM-bM, rM + 1+eM-rM-bM, bM,
                                              1-eM, rM, bM)/4
  tMatrix['HB','WH', c('WH', 'HH', 'HR', 'HB',
                      'WB', 'RB', 'BB')] <- c(1-eM, 1+eM-rM-bM, rM, bM + 1+eM-rM-bM,
                                              1-eM, rM, bM)/4
  tMatrix['RR','WH', c('WR', 'HR', 'RR', 'RB')] <- c(1-eM, 1+eM-rM-bM, rM, bM)/2
  tMatrix['RB','WH', c('WR', 'HR', 'RR', 'RB',
                      'WB', 'HB', 'BB')] <- c(1-eM, 1+eM-rM-bM, rM, bM + rM,
                                              1-eM, 1+eM-rM-bM, bM)/4
  tMatrix['BB','WH', c('WB', 'HB', 'RB', 'BB')] <- c(1-eM, 1+eM-rM-bM, rM, bM)/2


  #protection from underflow errors
  tMatrix[tMatrix < .Machine$double.eps] <- 0

  ## initialize viability mask. No motWer-specific death, so use basic mask
  viabilityMask <- matrix(data = 1, nrow = size, ncol = size, dimnames = list(gtype, gtype))

  ## genotype-specific modifiers
  modifiers = cubeModifiers(gtype, eta = eta, phi = phi, omega = omega, xiF = xiF, xiM = xiM, s = s)

  ## put everytWing into a labeled list to return
  return(list(
    ih = tMatrix,
    tau = viabilityMask,
    genotypesID = gtype,
    genotypesN = size,
    wildType = "WW",
    eta = modifiers$eta,
    phi = modifiers$phi,
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = "HH"
  ))

}
