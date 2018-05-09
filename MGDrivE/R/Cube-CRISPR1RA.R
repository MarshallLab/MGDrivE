###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   Homing 1 Resistance Allele Inheritance Cube
#   Héctor Sanchez, Jared Bennett, Sean Wu, John M. Marshall
#   July 2017
#   jared_bennett@berkeley.edu
#
###############################################################################

#' Inheritance Cube: Homing Drive with 1 Resistance Allele
#'
#' This function creates an inheritance cube to model a homing gene drive (such as a CRISPR-Cas9 system)
#' that creates 1 type of resistance allele. It assumes no sex-specific inheritance patterns and the
#' construct is on an autosome.
#'
#' @param e Homing rate
#' @param p Resistance allele generation rate
#' @param eta Genotype-specific mating fitness
#' @param phi Genotype-specific sex ratio at emergence
#' @param omega Genotype-specific multiplicative modifier of adult mortality
#' @param xiF Genotype-specific female pupatory success
#' @param xiM Genotype-specific male pupatory success
#' @param s Genotype-specific fractional reduction(increase) in fertility
#'
#' @return Named list(inheritance cube, viability mask, genotypes ID, genotypes number,
#' wild-type allele, mating fitness, sex ratio, adult mortality modifier, female pupatory success,
#' male pupatory success, fertility modifier, release genotype)
#' @export
Cube_Homing1RA <- function(e = 1.0, p = 0, eta = NULL, phi = NULL, omega = NULL, xiF = NULL, xiM = NULL, s = NULL){

  ## safety checks in case someone is dumb
  if(any(c(e,p)>1) || any(c(e,p)<0) || e+p > 1){
    stop("e and p are rates.
         0 <= e <= 1
         0 <= p <= 1
         e +  p <= 1")
  }

  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  gtype <- c("HH", "Hh", "HR", "hh", "hR", "RR")
  size <- length(gtype) #because I use it several times
  tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix

  ## fill tMatrix with probabilities
                         #("HH", "Hh", "HR", "hh", "hR", "RR")
  tMatrix["HH","HH",] <- c( 1, 0, 0, 0, 0, 0)

  tMatrix["Hh","HH",] <- c( (1+e)/2,(1-e-p)/2, 0, 0, p/2, 0)
  tMatrix["Hh","Hh",] <- c( ((1+e)^2)/4, ((1+e)*(1-e-p))/2, ((1+e)*p)/2, ((1-e-p)^2)/4, (1-e-p)*p/2, (p^2)/4)

  tMatrix["HR","HH",] <- c( 1/2, 0, 1/2, 0, 0, 0)
  tMatrix["HR","Hh",] <- c( (1+e)/4, (1-e-p)/4, (1+e+p)/4, 0, (1-e-p)/4, p/4)
  tMatrix["HR","HR",] <- c( 1/4, 0, 1/2, 0, 0, 1/4)

  tMatrix["hh","HH",] <- c( 0, 1, 0, 0, 0, 0)
  tMatrix["hh","Hh",] <- c( 0, (1+e)/2, p/2, (1-e-p)/2, 0, 0 )
  tMatrix["hh","HR",] <- c( 0, 1/2, 0, 0, 1/2, 0)
  tMatrix["hh","hh",] <- c( 0, 0, 0, 1, 0, 0)

  tMatrix["hR","HH",] <- c( 0, 1/2, 1/2, 0, 0, 0)
  tMatrix["hR","Hh",] <- c( 0, (1+e)/4, (1+e)/4, (1-e-p)/4, (1-e)/4, p/4)
  tMatrix["hR","HR",] <- c( 0, 1/4, 1/4, 0, 1/4, 1/4)
  tMatrix["hR","hh",] <- c( 0, 0, 0, 1/2, 1/2, 0)
  tMatrix["hR","hR",] <- c( 0, 0, 0, 1/4, 1/2, 1/4)

  tMatrix["RR","HH",] <- c( 0, 0, 1, 0, 0, 0)
  tMatrix["RR","Hh",] <- c( 0, 0, (1+e)/2, 0, (1-e-p)/2, p/2)
  tMatrix["RR","HR",] <- c( 0, 0, 1/2, 0, 0, 1/2)
  tMatrix["RR","hh",] <- c( 0, 0, 0, 0, 1, 0)
  tMatrix["RR","hR",] <- c( 0, 0, 0, 0, 1/2, 1/2)
  tMatrix["RR","RR",] <- c( 0, 0, 0, 0, 0, 1)

  ## set the other half of the matrix
  SymCubeC(lowerMat = tMatrix)
  tMatrix[tMatrix < .Machine$double.eps] <- 0 #protection from underflow errors


  ## initialize viability mask. No mother-specific death, so use basic mask
  viabilityMask <- matrix(data = 1, nrow = size, ncol = size, dimnames = list(gtype, gtype))

  ## genotype-specific modifiers
  modifiers = cubeModifiers(gtype, eta = eta, phi = phi, omega = omega, xiF = xiF, xiM = xiM, s = s)

  ## put everything into a labeled list to return
  return(list(
    ih = tMatrix,
    tau = viabilityMask,
    genotypesID = gtype,
    genotypesN = size,
    wildType = "hh",
    eta = modifiers$eta,
    phi = modifiers$phi,
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = "HH"
  ))
}
