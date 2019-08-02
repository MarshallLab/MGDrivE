###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   MEDEA
#   Héctor SancWez, Jared Bennett, Sean Wu, JoWn M. MarsWall
#   August 2017
#   jared_bennett@berkeley.edu
#
###############################################################################

#' Inheritance Cube: MEDEA (Maternal Effect Dominant Embryonic Arrest)
#'
#' This function creates an inheritance cube to model a MEDEA drive system. This
#' system was first discovered in flour beetles. It biases inheritance by expressing
#' a maternal toxin such that offspring die unless they express a zygotic antidote. \cr
#' This drive has 3 alleles at 1 locus:
#'  * W: Wild-type allele
#'  * M: MEDEA allele
#'  * R: Resistance allele
#'
#'
#' @param rM Breakdown of MEDEA allele, no homing/toxin/antidote, M -> R conversion
#' @param rW De novo resistance generation, W -> R conversion
#' @param Teff Efficacy of the toxin
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
cubeMEDEA <- function(rM = 0, rW = 0, Teff = 1.0, eta = NULL, phi = NULL,
                      omega = NULL, xiF = NULL, xiM = NULL, s = NULL){

  #rM=0; rW=0;

  ## safety checks
  if(any(c(rM,rW,Teff)<0) || any(c(rM,rW,Teff)>1)){
    stop("Parameters are rates.
         0 <= x <= 1")
  }


  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  gtype <- c('WW', 'WM', 'WR', 'MM', 'MR', 'RR')
  size <- length(gtype) #because I use it several times
  tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix

  ## fill tMatrix with probabilities
                          #( 'WW', 'WM', 'WR', 'MM', 'MR', 'RR')
  tMatrix['WW','WW', c('WW', 'WR', 'RR')] <- c( (1-rW)^2, 2*(1-rW)*rW, rW^2)

  tMatrix['WM','WW', ] <- c( (1-rW)^2, (1-rW)*(1-rM), (1-rW)*(rM+rW), 0, (1-rM)*rW, (rM+rW)*rW)/2
  tMatrix['WM','WM', ] <- c( (1-rW)^2, 2*(1-rW)*(1-rM), 2*(1-rW)*(rW+rM), (1-rM)^2, (1-rM)*(rW+rM), (rW+rM)^2)/4

  tMatrix['WR','WW', c('WW', 'WR', 'RR')] <- c( (1-rW)^2, (1-rW)*(1+rW) + (1-rW)*rW, rW*(1+rW))/2
  tMatrix['WR','WM', ] <- c( (1-rW)^2, (1-rW)*(1-rM), (1-rW)*(1+rW) + (1-rW)*(rW+rM), 0, (1+rW)*(1-rM), (1+rW)*(rW+rM))/4
  tMatrix['WR','WR', c('WW', 'WR', 'RR')] <- c( (1-rW)^2, 2*(1+rW)*(1-rW), (1+rW)^2)/4

  tMatrix['MM','WW', c('WM', 'WR', 'MR', 'RR')] <- c( (1-rM)*(1-rW), rM*(1-rW), (1-rM)*rW, rM*rW)
  tMatrix['MM','WM', ] <- c( 0, (1-rM)*(1-rW), (1-rW)*rM, (1-rM)^2, rM*(1-rM) + (1-rM)*(rW+rM), rM*(rW+rM))/2
  tMatrix['MM','WR', c('WM', 'WR', 'MR', 'RR')] <- c( (1-rW)*(1-rM), rM*(1-rW), (1+rW)*(1-rM), rM*(1+rW))/2
  tMatrix['MM','MM', c('MM', 'MR', 'RR')] <- c( (1-rM)^2, 2*rM*(1-rM), rM^2)

  tMatrix['MR','WW', c('WM', 'WR', 'MR', 'RR')] <- c( (1-rW)*(1-rM), (1-rW)*(1+rM), rW*(1-rM), rW*(1+rM))/2
  tMatrix['MR','WM', ] <- c( 0, (1-rW)*(1-rM), (1-rW)*(1+rM), (1-rM)^2, (1+rM)*(1-rM) + (1-rM)*(rW+rM), (1+rM)*(rM+rW))/4
  tMatrix['MR','WR', c('WM', 'WR', 'MR', 'RR')] <- c( (1-rW)*(1-rM), (1-rW)*(1+rM), (1+rW)*(1-rM), (1+rW)*(1+rM))/4
  tMatrix['MR','MM', c('MM', 'MR', 'RR')] <- c( (1-rM)^2, (1+rM)*(1-rM) + (1-rM)*rM, rM*(1+rM))/2
  tMatrix['MR','MR', c('MM', 'MR', 'RR')] <- c( (1-rM)^2, 2*(1+rM)*(1-rM), (1+rM)^2)/4

  tMatrix['RR','WW', c('WR', 'RR')] <- c( (1-rW), rW)
  tMatrix['RR','WM', c('WR', 'MR', 'RR')] <- c( (1-rW), (1-rM),  rW+rM )/2
  tMatrix['RR','WR', c('WR', 'RR')] <- c( (1-rW), (1+rW))/2
  tMatrix['RR','MM', c('MR', 'RR')] <- c( 1-rM, rM)
  tMatrix['RR','MR', c('MR', 'RR')] <- c((1-rM), (1+rM))/2
  tMatrix['RR','RR', 'RR'] <- 1

  ## set the other half of the matrix
  # Boolean matrix for subsetting, used several times
  boolMat <- upper.tri(x = tMatrix[ , ,1], diag = FALSE)
  # loop over depth, set upper triangle
  for(z in 1:size){tMatrix[ , ,z][boolMat] <- t(tMatrix[ , ,z])[boolMat]}

  tMatrix[tMatrix < .Machine$double.eps] <- 0 #protection from underflow errors

  ## initialize viability mask.
  viabilityMask <- array(data = 1, dim = c(size,size,size), dimnames = list(gtype, gtype, gtype))

  ## fill mother/offspring specific death, then muliply by efficacy of toxin and antidote
  for(slice in 1:size){
    viabilityMask[c('WM', 'MR'),slice, ] <- matrix( c( 1-Teff, 1, 1-Teff,1, 1, 1), nrow = 2, ncol = size, byrow = TRUE )
    viabilityMask['MM',slice, ] <- c( 1-2*Teff+Teff^2, 1, 1-2*Teff+Teff^2, 1, 1, 1)
  }


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
    releaseType = "MM"
  ))
}
