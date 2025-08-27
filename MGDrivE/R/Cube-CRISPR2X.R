###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   CRISPR 2 Resistance Alleles Inheritance Cube - Sex-Specific homing
#   Héctor Sanchez, Jared Bennett, Sean Wu, John Marshall
#   July 2017
#   jared_bennett@berkeley.edu
#   December 2018
#    Modified to reflect new cutting, homing, resistance generation rates
#
###############################################################################

#' Inheritance Cube: CRISPR (Clustered Regularly Interspaced Short Palindromic Repeats) X-linked with 2 Resistance Allele and Maternal Deposition
#'
#' This is an X-linked version of the 2 allele cube. It assumes that the construct
#' is on the X chromosome and there is no male homing. It also has maternal deposition,
#' i.e., when the male provides a W allele to a female with an H allele, some portion
#' are cut during oogenesis.
#' If the deposition parameters are zero (*D parameters), this is just
#' an X-linked drive.
#'
#' @param cF Female cutting rate
#' @param chF Female proper homing rate
#' @param crF Female no-cost resistance generation rate
#' @param dF Female deposition cutting rate
#' @param dhF Female deposition proper homing rate
#' @param drF Female deposition no-cost resistance generation rate
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
cubeXHomingDeposition <- function(cF=1.0, chF=0, crF=0, dF=0, dhF=0, drF=0,
                                   eta = NULL, phi = NULL, omega = NULL, xiF = NULL,
                                   xiM = NULL, s = NULL){

  ## safety checks
  if(any(c(cF,chF,crF,dF,dhF,drF)>1) || any(c(cF,chF,crF,dF,dhF,drF)<0)){
    stop("Parameters are rates, they must be between 0 and 1.")
  }


  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  gtype <- c('WW', 'WH', 'WR', 'WB', 'WY',
             'HH', 'HR', 'HB', 'HY',
             'RR', 'RB', 'RY', 'BB', 'BY')
  size <- length(gtype) #because I use it several times
  tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix


  ## fill tMatrix with probabilities
  tMatrix['WW','WY',c('WW','WY')] <- c(1,1)/2
  tMatrix['WW','HY',c('WH','WY')] <- c(1,1)/2
  tMatrix['WW','RY',c('WR','WY')] <- c(1,1)/2
  tMatrix['WW','BY',c('WB','WY')] <- c(1,1)/2

  tMatrix['WH','WY', ] <- c((1-cF)*(1-dF), (1+cF*chF)*(1-dF), (1-cF)*dF*drF+cF*(1-chF)*crF*(1-dF), (1-cF)*dF*(1-drF) + cF*(1-chF)*(1-crF)*(1-dF), 1-cF,
                            (1+cF*chF)*(dF*dhF), (1+cF*chF)*dF*(1-dhF)*drF, (1+cF*chF)*dF*(1-dhF)*(1-drF), (1+cF*chF),
                            cF*(1-chF)*crF*dF*drF, cF*(1-chF)*crF*dF*(1-drF) + cF*(1-chF)*(1-crF)*dF*drF, cF*(1-chF)*crF,
                            cF*(1-chF)*(1-crF)*dF*(1-drF), cF*(1-chF)*(1-crF))/4
  tMatrix['WH','HY',c('WH','HH','HR','HB',
                      'WY','HY','RY','BY')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF),
                                                 1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/4
  tMatrix['WH','RY',c('WR','HR','RR','RB',
                      'WY','HY','RY','BY')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF),
                                                 1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/4
  tMatrix['WH','BY',c('WB','HB','RB','BB',
                      'WY','HY','RY','BY')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF),
                                                 1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/4

  tMatrix['WR','WY',c('WW','WR','WY','RY')] <- c(1,1,1,1)/4
  tMatrix['WR','HY',c('WH','HR','WY','RY')] <- c(1,1,1,1)/4
  tMatrix['WR','RY',c('WR','RR','WY','RY')] <- c(1,1,1,1)/4
  tMatrix['WR','BY',c('WB','RB','WY','RY')] <- c(1,1,1,1)/4

  tMatrix['WB','WY',c('WW','WB','WY','BY')] <- c(1,1,1,1)/4
  tMatrix['WB','HY',c('WH','HB','WY','BY')] <- c(1,1,1,1)/4
  tMatrix['WB','RY',c('WR','RB','WY','BY')] <- c(1,1,1,1)/4
  tMatrix['WB','BY',c('WB','BB','WY','BY')] <- c(1,1,1,1)/4

  tMatrix['HH','WY',c('WH','HH','HR','HB','HY')] <- c(1-dF, dF*dhF, dF*(1-dhF)*drF, dF*(1-dhF)*(1-drF),1)/2
  tMatrix['HH','HY',c('HH','HY')] <- c(1,1)/2
  tMatrix['HH','RY',c('HR','HY')] <- c(1,1)/2
  tMatrix['HH','BY',c('HB','HY')] <- c(1,1)/2

  tMatrix['HR','WY',c('WH','WR','HH','HR',
                      'HB','RR','RB','HY','RY')] <- c(1-dF, 1-dF, dF*dhF, dF*(1-dhF)*drF,
                                                      dF*(1-dhF)*(1-drF), dF*drF, dF*(1-drF),1,1)/4
  tMatrix['HR','HY',c('HH','HR','HY','RY')] <- c(1,1,1,1)/4
  tMatrix['HR','RY',c('HR','RR','HY','RY')] <- c(1,1,1,1)/4
  tMatrix['HR','BY',c('HB','RB','HY','RY')] <- c(1,1,1,1)/4

  tMatrix['HB','WY',c('WH','WB','HH','HR',
                      'HB','RB','BB','HY','BY')] <- c(1-dF, 1-dF, dF*dhF, dF*(1-dhF)*drF,
                                                      dF*(1-dhF)*(1-drF), dF*drF, dF*(1-drF),1,1)/4
  tMatrix['HB','HY',c('HH','HB','HY','BY')] <- c(1,1,1,1)/4
  tMatrix['HB','RY',c('HR','RB','HY','BY')] <- c(1,1,1,1)/4
  tMatrix['HB','BY',c('HB','BB','HY','BY')] <- c(1,1,1,1)/4

  tMatrix['RR','WY',c('WR','RY')] <- c(1,1)/2
  tMatrix['RR','HY',c('HR','RY')] <- c(1,1)/2
  tMatrix['RR','RY',c('RR','RY')] <- c(1,1)/2
  tMatrix['RR','BY',c('RB','RY')] <- c(1,1)/2

  tMatrix['RB','WY',c('WR','WB','RY','BY')] <- c(1,1,1,1)/4
  tMatrix['RB','HY',c('HR','HB','RY','BY')] <- c(1,1,1,1)/4
  tMatrix['RB','RY',c('RR','RB','RY','BY')] <- c(1,1,1,1)/4
  tMatrix['RB','BY',c('RB','BB','RY','BY')] <- c(1,1,1,1)/4

  tMatrix['BB','WY',c('WB','BY')] <- c(1,1)/2
  tMatrix['BB','HY',c('HB','BY')] <- c(1,1)/2
  tMatrix['BB','RY',c('RB','BY')] <- c(1,1)/2
  tMatrix['BB','BY',c('BB','BY')] <- c(1,1)/2


  #protection from underflow errors
  tMatrix[tMatrix < .Machine$double.eps] <- 0

  ## initialize viability mask. No mother/father-specific death, so use basic mask
  viabilityMask <- array(data = 1, dim = c(size,size,size), dimnames = list(gtype, gtype, gtype))

  ## genotype-specific modifiers
  if(!is.null(phi)){
    stop("This cube has a special phi, due to being male/female specific.
         Please edit it after the cube is built, if you are absolutely sure it needs to change.")
  }
  phi = setNames(object = c(1,1,1,1,0,1,1,1,0,1,1,0,1,0), nm = gtype)
  modifiers = cubeModifiers(gtype, eta = eta, phi = phi, omega = omega, xiF = xiF, xiM = xiM, s = s)

  ## put everything into a labeled list to return
  return(list(
    ih = tMatrix,
    tau = viabilityMask,
    genotypesID = gtype,
    genotypesN = size,
    wildType = c("WW","WY"),
    eta = modifiers$eta,
    phi = modifiers$phi,
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = "HY"
  ))
}
