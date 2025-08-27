###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   Small-molecule CRISPR cube with 1 resistant allele - Sex-Specific homing
#   Héctor Sanchez, Jared Bennett, Sean Wu, John Marshall
#   July 2017
#   jared_bennett@berkeley.edu
#   December 2018
#    Modified to reflect new cutting, homing, resistance generation rates
#
###############################################################################

#' Inheritance Cube: CRISPR-SM (Clustered Regularly Interspaced Short Palindromic Repeats) with Small-Molecule Induction and 1 Resistance Allele and Maternal Deposition
#'
#' This is a sex-specific version of CRISPR with small-molecule induced homing.
#' It assumes that the construct is on an autosome and there can be different
#' male/female homing rates. It also has maternal deposition, i.e., when the
#' male provides a W allele to a female with a H allele, some portion are cut
#' during oogenesis. Additionally, this cube is designed for small-molecule
#' induction, i.e., with the SM branch of MGDrivE. It allows the homing (H) allele
#' to be turned off into an O allele, which inherits stably, and so that all
#' offspring of H individuals are O until turned on with the spray.
#' If the maternal deposition parameters are zero (d* parameters), this is a normal
#' CRISPR drive.
#'
#' @param cM Male homing rate
#' @param cF Female homing rate
#' @param dF Female deposition homing rate
#' @param chF Female correct homing rate
#' @param chM Male correct homing rate
#' @param dhF Female correct deposition rate
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
cubeHomingDriveSM <- function(cM = 1.0, cF = 1.0, dF=0, chM = 0, chF = 0,
                            dhF = 0, eta = NULL, phi = NULL,
                            omega = NULL, xiF = NULL, xiM = NULL, s = NULL){

  ## safety checks
  if(any(c(cM, cF, dF, chM, chF, dhF)>1) || any(c(cM, cF, dF, chM, chF, dhF)<0)){
    stop("Parameters are rates, they must be between 0 and 1.")
  }


  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  gtype <- c('WW', 'WH', 'WR', 'WO', 'HH', 'HR', 'RR', 'RO', 'OO')
  size <- length(gtype) #because I use it several times
  tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix

  ## fill tMatrix with probabilities
  #('WW', 'WH', 'WR', 'WO', 'HH', 'HR', 'RR', 'RO', 'OO')
  tMatrix['WW','WW', 'WW'] <- 1

  tMatrix['WR','WW', c('WW', 'WR')] <- c( 1, 1)/2
  tMatrix['WR','WR', c('WW', 'WR', 'RR')] <- c( 1/2, 1, 1/2)/2

  tMatrix['WO','WW', c('WW', 'WO')] <- c( 1, 1)/2
  tMatrix['WO','WR', c('WW', 'WR', 'WO', 'RO')] <- c( 1, 1, 1, 1)/4
  tMatrix['WO','WO', c('WW', 'WO', 'OO')] <- c( 1/2, 1, 1/2)/2

  tMatrix['HH','HH', 'OO'] <- 1

  tMatrix['HR','HH', c('OO', 'RO')] <- c( 1, 1)/2
  tMatrix['HR','HR', c('OO', 'RO', 'RR')] <- c( 1/2, 1, 1/2)/2

  tMatrix['RR','WW', 'WR'] <- 1
  tMatrix['RR','WR', c('WR', 'RR')] <- c( 1, 1)/2
  tMatrix['RR','WO', c('WR', 'RO')] <- c( 1, 1)/2
  tMatrix['RR','HH', 'RO'] <- 1
  tMatrix['RR','HR', c('RO', 'RR')] <- c( 1, 1)/2
  tMatrix['RR','RR', 'RR'] <- 1

  tMatrix['RO','WW', c('WR', 'WO')] <- c( 1, 1)/2
  tMatrix['RO','WR', c('WR', 'WO', 'RR', 'RO')] <- c( 1, 1, 1, 1)/4
  tMatrix['RO','WO', c('WR', 'WO', 'RO', 'OO')] <- c( 1, 1, 1, 1)/4
  tMatrix['RO','HH', c('RO', 'OO')] <- c( 1, 1)/2
  tMatrix['RO','HR', c('OO', 'RR', 'RO')] <- c( 1, 1, 2)/4
  tMatrix['RO','RR', c('RR', 'RO')] <- c( 1, 1)/2
  tMatrix['RO','RO', c('RR', 'RO', 'OO')] <- c( 1/2, 1, 1/2)/2

  tMatrix['OO','WW', 'WO'] <- 1
  tMatrix['OO','WR', c('WO', 'RO')] <- c( 1, 1)/2
  tMatrix['OO','WO', c('WO', 'OO')] <- c( 1, 1)/2
  tMatrix['OO','HH', 'OO'] <- 1
  tMatrix['OO','HR', c('OO', 'RO')] <- c( 1, 1)/2
  tMatrix['OO','RR', 'RO'] <- 1
  tMatrix['OO','RO', c('RO', 'OO')] <- c( 1, 1)/2
  tMatrix['OO','OO', 'OO'] <- 1

  ## set the other half of the matrix that is symmetric
  # Boolean matrix for subsetting, used several times
  boolMat <- upper.tri(x = tMatrix[ , ,1], diag = FALSE)
  # loop over depth, set upper triangle
  for(z in 1:size){tMatrix[ , ,z][boolMat] <- t(tMatrix[ , ,z])[boolMat]}


  ## fill asymmetric parts of tMatrix
  #female specific homing, except for WHxWH
  tMatrix['WH','WW', ] <- c((1-cF)*(1-dF), 0, cF*(1-chF)*(1-dF) + (1-cF)*dF*(1-dhF),
                            (1+cF*chF)*(1-dF) + (1-cF)*dF*dhF, 0, 0,
                            cF*(1-chF)*dF*(1-dhF), cF*(1-chF)*dF*dhF + (1+cF*chF)*dF*(1-dhF),
                            (1+cF*chF)*dF*dhF)/2

  tMatrix['WH','WH', ] <- c((1-cF)*(1-cM)*(1-dF),
                            0,
                            cF*(1-chF)*(1-cM)*(1-dF) + (1-cF)*((1-cM)*dF*(1-dhF) + cM*(1-chM)),
                            (1+cF*chF)*(1-cM)*(1-dF) + (1-cF)*((1-cM)*dF*dhF + 1+cM*chM),
                            0,
                            0,
                            cF*(1-chF)*((1-cM)*dF*(1-dhF) + cM*(1-chM)),
                            cF*(1-chF)*((1-cM)*dF*dhF + 1+cM*chM) + (1+cF*chF)*((1-cM)*dF*(1-dhF) + cM*(1-chM)),
                            (1+cF*chF)*((1-cM)*dF*dhF + 1+cM*chM))/4

  tMatrix['WH','WR', ] <- c((1-cF)*(1-dF), 0, cF*(1-chF)*(1-dF) + (1-cF)*(1+dF*(1-dhF)),
                            (1+cF*chF)*(1-dF) + (1-cF)*dF*dhF, 0, 0,
                            cF*(1-chF)*(1+dF*(1-dhF)), cF*(1-chF)*dF*dhF + (1+cF*chF)*(1+dF*(1-dhF)),
                            (1+cF*chF)*dF*dhF)/4

  tMatrix['WH','WO', ] <- c((1-cF)*(1-dF), 0, cF*(1-chF)*(1-dF) + (1-cF)*dF*(1-dhF),
                            (1+cF*chF)*(1-dF) + (1-cF)*(1+dF*dhF), 0, 0,
                            cF*(1-chF)*dF*(1-dhF), cF*(1-chF)*(1+dF*dhF) + (1+cF*chF)*dF*(1-dhF),
                            (1+cF*chF)*(1+dF*dhF))/4

  tMatrix['WH','HH',c('WO', 'OO', 'RO')] <- c(1-cF, 1+cF*chF, cF*(1-chF))/2
  tMatrix['WH','HR',c('WO', 'OO', 'RO',
                      'WR', 'RR')] <- c(1-cF, 1+cF*chF, cF*(1-chF) + 1+cF*chF,
                                        1-cF, cF*(1-chF))/4

  tMatrix['WH','RR',c('WR', 'RO', 'RR')] <- c(1-cF, 1+cF*chF, cF*(1-chF))/2
  tMatrix['WH','RO',c('WR', 'RO', 'RR',
                      'WO', 'OO')] <- c(1-cF, 1+cF*chF + cF*(1-chF), cF*(1-chF),
                                        1-cF, 1+cF*chF)/4
  tMatrix['WH','OO',c('WO', 'OO', 'RO')] <- c(1-cF, 1+cF*chF, cF*(1-chF))/2


  # female deposition things
  tMatrix['HH','WW', c('WO', 'OO', 'RO')] <- c( 1-dF, dF*dhF, dF*(1-dhF))

  tMatrix['HH','WR', c('WO', 'OO', 'RO')] <- c( 1-dF, dF*dhF, dF*(1-dhF) + 1)/2
  tMatrix['HH','WO', c('WO', 'OO', 'RO')] <- c( 1-dF, 1+dF*dhF, dF*(1-dhF))/2

  tMatrix['HR','WW', c('WO', 'WR', 'OO',
                       'RO', 'RR')] <- c( 1-dF, 1-dF, dF*dhF,
                                          dF*dhF + dF*(1-dhF), dF*(1-dhF))/2
  tMatrix['HR','WR', c('WO', 'WR', 'OO',
                       'RO', 'RR')] <- c( 1-dF, 1-dF, dF*dhF,
                                          1 + dF*dhF + dF*(1-dhF), 1 + dF*(1-dhF))/4
  tMatrix['HR','WO', c('WO', 'WR', 'OO',
                       'RO', 'RR')] <- c( 1-dF, 1-dF, 1 + dF*dhF,
                                          1 + dF*dhF + dF*(1-dhF), dF*(1-dhF))/4


  #male specific homing
  tMatrix['WW','WH', c('WW', 'WO', 'WR')] <- c(1-cM, 1+cM*chM, cM*(1-chM))/2
  tMatrix['WR','WH', c('WW', 'WO', 'WR',
                       'RO', 'RR')] <- c(1-cM, 1+cM*chM, cM*(1-chM) + 1-cM,
                                               1+cM*chM, cM*(1-chM))/4
  tMatrix['WO','WH', c('WW', 'WO', 'WR',
                       'RO', 'OO')] <- c(1-cM, 1+cM*chM + 1-cM, cM*(1-chM),
                                         cM*(1-chM), 1+cM*chM)/4
  tMatrix['HH','WH', c('WO','RO','OO')] <- c((1-cM)*(1-dF), cM*(1-chM) + (1-cM)*dF*(1-dhF), 1+cM*chM + (1-cM)*dF*dhF)/2
  tMatrix['HR','WH', c('WO', 'RO', 'OO',
                       'WR', 'RR')] <- c((1-cM)*(1-dF), cM*(1-chM) + (1-cM)*dF*(1-dhF) + (1-cM)*dF*dhF + 1+cM*chM, 1+cM*chM + (1-cM)*dF*dhF,
                                         (1-cM)*(1-dF), cM*(1-chM) + (1-cM)*dF*(1-dhF))/4
  tMatrix['RR','WH', c('WR', 'RO', 'RR')] <- c(1-cM, 1+cM*chM, cM*(1-chM))/2
  tMatrix['RO','WH', c('WR', 'RO', 'RR',
                       'WO', 'OO')] <- c(1-cM, 1+cM*chM + cM*(1-chM), cM*(1-chM),
                                         1-cM, 1+cM*chM)/4
  tMatrix['OO','WH', c('WO', 'RO', 'OO')] <- c(1-cM, cM*(1-chM), 1+cM*chM)/2


  #male stuff from female deposition
  tMatrix['WW','HH', 'WO'] <- 1
  tMatrix['WR','HH', c('WO', 'RO')] <- c( 1, 1)/2
  tMatrix['WO','HH', c('WO', 'OO')] <- c( 1, 1)/2

  tMatrix['WW','HR', c('WO', 'WR')] <- c( 1, 1)/2
  tMatrix['WR','HR', c('WO', 'WR', 'RO', 'RR')] <- c( 1, 1, 1, 1)/4
  tMatrix['WO','HR', c('WO', 'WR', 'OO', 'RO')] <- c( 1, 1, 1, 1)/4


  #protection from underflow errors
  tMatrix[tMatrix < .Machine$double.eps] <- 0

  ## initialize viability mask. No mother/father-specific death, so use basic mask
  viabilityMask <- array(data = 1, dim = c(size,size,size), dimnames = list(gtype, gtype, gtype))

  ## genotype-specific modifiers
  modifiers = cubeModifiers(gtype, eta = eta, phi = phi, omega = omega, xiF = xiF, xiM = xiM, s = s)

  ## put everything into a labeled list to return
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
    releaseType = "OO"
  ))

}
