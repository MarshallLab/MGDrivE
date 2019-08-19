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
#   December 2018
#    Modified to reflect new cutting, homing, resistance generation rates
#
###############################################################################

#' Inheritance Cube: CRISPR (Clustered Regularly Interspaced Short Palindromic Repeats) with 2 Resistance Alleles and maternal deposition
#'
#' This is a sex-specific version of the original cube. It assumes that the construct
#' is on an autosome and there can be different male/female homing rates. It also has
#' maternal deposition, ie, when the male provides a W allele to a female with a H allele,
#' some portion are cut during oogenesis.
#' If the maternal deposition parameters are zero (d* parameters), this is a normal
#' CRISPR drive
#'
#' @param cM Male homing rate
#' @param cF Female homing rate
#' @param dF Female deposition homing rate
#' @param chF Female correct homing rate
#' @param crF Female resistance generating rate
#' @param chM Male correct homing rate
#' @param crM Male resistance generating rate
#' @param dhF Female correct deposition rate
#' @param drF Female resistance deposition rate
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
cubeHomingDrive <- function(cM = 1.0, cF = 1.0, dF=0, chM = 0, crM = 0, chF = 0,
                             crF = 0, dhF = 0, drF = 0, eta = NULL, phi = NULL,
                            omega = NULL, xiF = NULL, xiM = NULL, s = NULL){

  ## safety cWecks in case someone is dumb
  if(any(c(cM, cF, dF, chM, crM, chF, crF, dhF, drF)>1) || any(c(cM, cF, dF, chM, crM, chF, crF, dhF, drF)<0)){
    stop("Parameters are rates, they must be between 0 and 1.")
  }


  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  gtype <- c('WW', 'WH', 'WR', 'WB', 'HH', 'HR', 'HB', 'RR', 'RB', 'BB')
  size <- length(gtype) #because I use it several times
  tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix

  ## fill tMatrix with probabilities
  #('WW', 'WH', 'WR', 'WB', 'HH', 'HR', 'HB', 'RR', 'RB', 'BB')
  tMatrix['WW','WW', 'WW'] <- 1

  tMatrix['WR','WW', c('WW', 'WR')] <- c( 1, 1)/2
  tMatrix['WR','WR', c('WW', 'WR', 'RR')] <- c( 1/2, 1, 1/2)/2

  tMatrix['WB','WW', c('WW', 'WB')] <- c( 1, 1)/2
  tMatrix['WB','WR', c('WW', 'WR', 'WB', 'RB')] <- c( 1, 1, 1, 1)/4
  tMatrix['WB','WB', c('WW', 'WB', 'BB')] <- c( 1/2, 1, 1/2)/2

  tMatrix['HH','HH', 'HH'] <- 1

  tMatrix['HR','HH', c('HH', 'HR')] <- c( 1, 1)/2
  tMatrix['HR','HR', c('HH', 'HR', 'RR')] <- c( 1/2, 1, 1/2)/2

  tMatrix['HB','HH', c('HH', 'HB')] <- c( 1, 1)/2
  tMatrix['HB','HR', c('HH', 'HR', 'HB', 'RB')] <- c( 1, 1, 1, 1)/4
  tMatrix['HB','HB', c('HH', 'HB', 'BB')] <- c( 1/2, 1, 1/2)/2

  tMatrix['RR','WW', 'WR'] <- 1
  tMatrix['RR','WR', c('WR', 'RR')] <- c( 1, 1)/2
  tMatrix['RR','WB', c('WR', 'RB')] <- c( 1, 1)/2
  tMatrix['RR','HH', 'HR'] <- 1
  tMatrix['RR','HR', c('HR', 'RR')] <- c( 1, 1)/2
  tMatrix['RR','HB', c('HR', 'RB')] <- c( 1, 1)/2
  tMatrix['RR','RR', 'RR'] <- 1

  tMatrix['RB','WW', c('WR', 'WB')] <- c( 1, 1)/2
  tMatrix['RB','WR', c('WR', 'WB', 'RR', 'RB')] <- c( 1, 1, 1, 1)/4
  tMatrix['RB','WB', c('WR', 'WB', 'RB', 'BB')] <- c( 1, 1, 1, 1)/4
  tMatrix['RB','HH', c('HR', 'HB')] <- c( 1, 1)/2
  tMatrix['RB','HR', c('HR', 'HB', 'RR', 'RB')] <- c( 1, 1, 1, 1)/4
  tMatrix['RB','HB', c('HR', 'HB', 'RB', 'BB')] <- c( 1, 1, 1, 1)/4
  tMatrix['RB','RR', c('RR', 'RB')] <- c( 1, 1)/2
  tMatrix['RB','RB', c('RR', 'RB', 'BB')] <- c( 1/2, 1, 1/2)/2

  tMatrix['BB','WW', 'WB'] <- 1
  tMatrix['BB','WR', c('WB', 'RB')] <- c( 1, 1)/2
  tMatrix['BB','WB', c('WB', 'BB')] <- c( 1, 1)/2
  tMatrix['BB','HH', 'HB'] <- 1
  tMatrix['BB','HR', c('HB', 'RB')] <- c( 1, 1)/2
  tMatrix['BB','HB', c('HB', 'BB')] <- c( 1, 1)/2
  tMatrix['BB','RR', 'RB'] <- 1
  tMatrix['BB','RB', c('RB', 'BB')] <- c( 1, 1)/2
  tMatrix['BB','BB', 'BB'] <- 1

  ## set the other half of the matrix that is symmetric
  # Boolean matrix for subsetting, used several times
  boolMat <- upper.tri(x = tMatrix[ , ,1], diag = FALSE)
  # loop over depth, set upper triangle
  for(z in 1:size){tMatrix[ , ,z][boolMat] <- t(tMatrix[ , ,z])[boolMat]}


  ## fill asymmetric parts of tMatrix
  #female specific homing, except for WHxWH
  tMatrix['WH','WW', ] <- c((1-cF)*(1-dF), (1+cF*chF)*(1-dF), (1-cF)*dF*drF + (cF*(1-chF)*crF)*(1-dF), (1-cF)*dF*(1-drF) + (cF*(1-chF)*(1-crF))*(1-dF),
                            (1+cF*chF)*dF*dhF, (1+cF*chF)*dF*(1-dhF)*drF, (1+cF*chF)*dF*(1-dhF)*(1-drF),
                            (cF*(1-chF)*crF)*dF*drF, (cF*(1-chF)*crF)*dF*(1-drF) + (cF*(1-chF)*(1-crF))*dF*drF, (cF*(1-chF)*(1-crF))*dF*(1-drF))/2

  tMatrix['WH','WH', ] <- c((1-cF)*(1-cM)*(1-dF), (1+cF*chF)*(1-cM)*(1-dF) + (1-cF)*(1+cM*chM),
                            (1-cF)*(1-cM)*dF*drF + cF*(1-chF)*crF*(1-cM)*(1-dF) + (1-cF)*cM*(1-chM)*crM,
                            (1-cF)*(1-cM)*dF*(1-drF) + cF*(1-chF)*(1-crF)*(1-cM)*(1-dF) + (1-cF)*cM*(1-chM)*(1-crM),
                            (1+cF*chF)*(1-cM)*dF*dhF + (1+cF*chF)*(1+cM*chM),
                            (1+cF*chF)*(1-cM)*dF*(1-dhF)*drF + cF*(1-chF)*crF*(1+cM*chM) + (1+cF*chF)*cM*(1-chM)*crM,
                            (1+cF*chF)*(1-cM)*dF*(1-dhF)*(1-drF) + cF*(1-chF)*(1-crF)*(1+cM*chM) + (1+cF*chF)*cM*(1-chM)*(1-crM),
                            cF*(1-chF)*crF*(1-cM)*dF*drF + cF*(1-chF)*crF*cM*(1-chM)*crM,
                            cF*(1-chF)*crF*(1-cM)*dF*(1-drF) + cF*(1-chF)*(1-crF)*(1-cM)*dF*drF + cF*(1-chF)*(1-crF)*cM*(1-chM)*crM + cF*(1-chF)*crF*cM*(1-chM)*(1-crM),
                            cF*(1-chF)*(1-crF)*(1-cM)*dF*(1-drF) + cF*(1-chF)*(1-crF)*cM*(1-chM)*(1-crM))/4

  tMatrix['WH','WR', ] <- c((1-cF)*(1-dF), (1+cF*chF)*(1-dF), (1-cF)*dF*drF + (cF*(1-chF)*crF)*(1-dF) + 1-cF, (1-cF)*dF*(1-drF) + (cF*(1-chF)*(1-crF))*(1-dF),
                            (1+cF*chF)*dF*dhF, (1+cF*chF)*dF*(1-dhF)*drF + 1+cF*chF, (1+cF*chF)*dF*(1-dhF)*(1-drF),
                            (cF*(1-chF)*crF)*dF*drF + cF*(1-chF)*crF, (cF*(1-chF)*crF)*dF*(1-drF) + (cF*(1-chF)*(1-crF))*dF*drF + cF*(1-chF)*(1-crF),
                            (cF*(1-chF)*(1-crF))*dF*(1-drF))/4
  tMatrix['WH','WB', ] <- c((1-cF)*(1-dF), (1+cF*chF)*(1-dF), (1-cF)*dF*drF + (cF*(1-chF)*crF)*(1-dF), (1-cF)*dF*(1-drF) + (cF*(1-chF)*(1-crF))*(1-dF) + 1-cF,
                            (1+cF*chF)*dF*dhF, (1+cF*chF)*dF*(1-dhF)*drF, (1+cF*chF)*dF*(1-dhF)*(1-drF) + 1+cF*chF,
                            (cF*(1-chF)*crF)*dF*drF, (cF*(1-chF)*crF)*dF*(1-drF) + (cF*(1-chF)*(1-crF))*dF*drF + cF*(1-chF)*crF,
                            (cF*(1-chF)*(1-crF))*dF*(1-drF) + cF*(1-chF)*(1-crF))/4

  tMatrix['WH','HH',c('WH', 'HH', 'HR', 'HB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/2
  tMatrix['WH','HR',c('WH', 'HH', 'HR', 'HB',
                      'WR', 'RR', 'RB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF + 1+cF*chF, cF*(1-chF)*(1-crF),
                                              1-cF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/4
  tMatrix['WH','HB',c('WH', 'HH', 'HR', 'HB',
                      'WB', 'RB', 'BB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF) + 1+cF*chF,
                                              1-cF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/4
  tMatrix['WH','RR',c('WR', 'HR', 'RR', 'RB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/2
  tMatrix['WH','RB',c('WR', 'HR', 'RR', 'RB',
                      'WB', 'HB', 'BB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF) + cF*(1-chF)*crF,
                                              1-cF, 1+cF*chF, cF*(1-chF)*(1-crF))/4
  tMatrix['WH','BB',c('WB', 'HB', 'RB', 'BB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/2


  # female deposition things
  tMatrix['HH','WW', c('WH', 'HH', 'HR', 'HB')] <- c( 1-dF, dF*dhF, dF*(1-dhF)*drF, dF*(1-dhF)*(1-drF))

  tMatrix['HH','WR', c('WH', 'HH', 'HR', 'HB')] <- c( 1-dF, dF*dhF, dF*(1-dhF)*drF + 1, dF*(1-dhF)*(1-drF))/2
  tMatrix['HH','WB', c('WH', 'HH', 'HR', 'HB')] <- c( 1-dF, dF*dhF, dF*(1-dhF)*drF, dF*(1-dhF)*(1-drF) + 1)/2

  tMatrix['HR','WW', c('WH', 'WR', 'HH', 'HR',
                       'HB', 'RR', 'RB')] <- c( 1-dF, 1-dF, dF*dhF, dF*(1-dhF)*drF,
                                                dF*(1-dhF)*(1-drF), dF*drF, dF*(1-drF))/2
  tMatrix['HR','WR', c('WH', 'WR', 'HH', 'HR',
                       'HB', 'RR', 'RB')] <- c( 1-dF, 1-dF, dF*dhF, dF*(1-dhF)*drF + 1,
                                                dF*(1-dhF)*(1-drF), dF*drF + 1, dF*(1-drF))/4
  tMatrix['HR','WB', c('WH', 'WR', 'HH', 'HR',
                       'HB', 'RR', 'RB')] <- c( 1-dF, 1-dF, dF*dhF, dF*(1-dhF)*drF,
                                                dF*(1-dhF)*(1-drF) + 1, dF*drF, dF*(1-drF) + 1)/4

  tMatrix['HB','WW', c('WH', 'WB', 'HH', 'HR',
                       'HB', 'RB', 'BB')] <- c( 1-dF, 1-dF, dF*dhF, dF*(1-dhF)*drF,
                                                dF*(1-dhF)*(1-drF), dF*drF, dF*(1-drF))/2
  tMatrix['HB','WR', c('WH', 'WB', 'HH','HR',
                       'HB','RB', 'BB')] <- c( 1-dF, 1-dF, dF*dhF, dF*(1-dhF)*drF + 1,
                                                dF*(1-dhF)*(1-drF), dF*drF + 1, dF*(1-drF))/4
  tMatrix['HB','WB', c('WH', 'WB', 'HH','HR',
                       'HB','RB', 'BB')] <- c( 1-dF, 1-dF, dF*dhF, dF*(1-dhF)*drF,
                                                dF*(1-dhF)*(1-drF) + 1, dF*drF, dF*(1-drF) + 1)/4


  #male specific homing
  tMatrix['WW','WH', c('WW', 'WH', 'WR', 'WB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/2
  tMatrix['WR','WH', c('WW', 'WH', 'WR', 'WB',
                       'HR', 'RR', 'RB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM + 1-cM, cM*(1-chM)*(1-crM),
                                               1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/4
  tMatrix['WB','WH', c('WW', 'WH', 'WR', 'WB',
                       'HB', 'RB', 'BB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + 1-cM,
                                               1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/4
  tMatrix['HH','WH', c('WH', 'HH',
                       'HR', 'HB')] <- c((1-cM)*(1-dF), 1+cM*chM + (1-cM)*dF*dhF,
                                         cM*(1-chM)*crM + (1-cM)*dF*(1-dhF)*drF, cM*(1-chM)*(1-crM) + (1-cM)*dF*(1-dhF)*(1-drF))/2
  tMatrix['HR','WH', c('WH', 'HH', 'HR', 'HB',
                       'WR', 'RR', 'RB')] <- c((1-cM)*(1-dF), 1+cM*chM + (1-cM)*dF*dhF, cM*(1-chM)*crM + 1+cM*chM + (1-cM)*dF*(1-dhF)*drF,
                                               cM*(1-chM)*(1-crM) + (1-cM)*dF*(1-dhF)*(1-drF), (1-cM)*(1-dF),
                                               cM*(1-chM)*crM + (1-cM)*dF*drF, cM*(1-chM)*(1-crM) + (1-cM)*dF*(1-drF))/4
  tMatrix['HB','WH', c('WH', 'HH', 'HR', 'HB',
                      'WB', 'RB', 'BB')] <- c((1-cM)*(1-dF), 1+cM*chM + (1-cM)*dF*dhF, cM*(1-chM)*crM + (1-cM)*dF*(1-dhF)*drF,
                                              cM*(1-chM)*(1-crM) + 1+cM*chM + (1-cM)*dF*(1-dhF)*(1-drF), (1-cM)*(1-dF),
                                              cM*(1-chM)*crM + (1-cM)*dF*drF, cM*(1-chM)*(1-crM) + (1-cM)*dF*(1-drF))/4
  tMatrix['RR','WH', c('WR', 'HR', 'RR', 'RB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/2
  tMatrix['RB','WH', c('WR', 'HR', 'RR', 'RB',
                       'WB', 'HB', 'BB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + cM*(1-chM)*crM,
                                               1-cM, 1+cM*chM, cM*(1-chM)*(1-crM))/4
  tMatrix['BB','WH', c('WB', 'HB', 'RB', 'BB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/2


  #male stuff from female deposition
  tMatrix['WW','HH', 'WH'] <- 1
  tMatrix['WR','HH', c('WH', 'HR')] <- c( 1, 1)/2
  tMatrix['WB','HH', c('WH', 'HB')] <- c( 1, 1)/2

  tMatrix['WW','HR', c('WH', 'WR')] <- c( 1, 1)/2
  tMatrix['WR','HR', c('WH', 'WR', 'HR', 'RR')] <- c( 1, 1, 1, 1)/4
  tMatrix['WB','HR', c('WH', 'WR', 'HB', 'RB')] <- c( 1, 1, 1, 1)/4

  tMatrix['WW','HB', c('WH', 'WB')] <- c( 1, 1)/2
  tMatrix['WR','HB', c('WH', 'WB', 'HR', 'RB')] <- c( 1, 1, 1, 1)/4
  tMatrix['WB','HB', c('WH', 'WB', 'HB', 'BB')] <- c( 1, 1, 1, 1)/4


  #protection from underflow errors
  tMatrix[tMatrix < .Machine$double.eps] <- 0

  ## initialize viability mask. No mother/father-specific death, so use basic mask
  viabilityMask <- array(data = 1L, dim = c(size,size,size), dimnames = list(gtype, gtype, gtype))

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
