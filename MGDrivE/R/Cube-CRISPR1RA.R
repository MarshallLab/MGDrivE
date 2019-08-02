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
#   December 2018
#    Update to reflect cutting, homing, resistance generation rates
#    Did not add male/female specific stuff
#    If any of that is needed, used the CRISPR2RA cube
#   Jan 2019
#    Did not actually update this.
#
###############################################################################

#' Inheritance Cube: Homing Drive with 1 Resistance Allele
#'
#' This function creates an inheritance cube to model a homing gene drive (such as a CRISPR-Cas9 system)
#' that creates 1 type of resistance allele. It assumes no sex-specific inheritance patterns and the
#' construct is on an autosome.
#'
#' @param c Cutting rate
#' @param ch Successful homing rate rate
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
cubeHoming1RA <- function(c = 1.0, ch = 0, eta = NULL, phi = NULL, omega = NULL,
                           xiF = NULL, xiM = NULL, s = NULL){

  ## safety checks in case someone is dumb
  if(any(c(c,ch)>1) || any(c(c,ch)<0) ){
    stop("e and p are rates.
         0 <= e <= 1
         0 <= p <= 1")
  }

  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  gtype <- c("HH", "HW", "HR", "WW", "WR", "RR")
  size <- length(gtype) #because I use it several times
  tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix

  ## fill tMatrix with probabilities
  #("HH", "HW", "HR", "WW", "WR", "RR")
  tMatrix["HH","HH","HH"] <-1

  tMatrix["HW","HH",c("HH","HW","HR")] <- c( 1+c*ch, 1-c, c*(1-ch))/2
  tMatrix["HW","HW",] <- c( (1+c*ch)^2, 2*(1+c*ch)*(1-c), 2*(1+c*ch)*(c*(1-ch)),
                            (1-c)^2, 2*(1-c)*(c*(1-ch)), (c*(1-ch))^2)/4

  tMatrix["HR","HH",c("HH","HR")] <- c( 1, 1)/2
  tMatrix["HR","HW",] <- c( 1+c*ch, 1-c, 1+c*ch + c*(1-ch), 0, 1-c, c*(1-ch))/4
  tMatrix["HR","HR",c("HH","HR","RR")] <- c( 1/2, 1, 1/2)/2

  tMatrix["WW","HH","HW"] <- 1
  tMatrix["WW","HW",c("HW","WW","WR")] <- c( 1+c*ch, 1-c, c*(1-ch))/2
  tMatrix["WW","HR",c("HW","WR")] <- c( 1, 1)/2
  tMatrix["WW","WW","WW"] <- 1

  tMatrix["WR","HH",c("HW","HR")] <- c( 1, 1)/2
  tMatrix["WR","HW",] <- c( 0, 1+c*ch, 1+c*ch, 1-c, 1-c+c*(1-ch), c*(1-ch))/4
  tMatrix["WR","HR",c("HW","HR","WR","RR")] <- c( 1, 1, 1, 1)/4
  tMatrix["WR","WW",c("WW","WR")] <- c( 1, 1)/2
  tMatrix["WR","WR",c("WW","WR","RR")] <- c( 1/2, 1, 1/2)/2

  tMatrix["RR","HH","HR"] <- 1
  tMatrix["RR","HW",c("HR","WR","RR")] <- c( 1+c*ch, 1-c, c*(1-ch))/2
  tMatrix["RR","HR",c("HR","RR")] <- c(1, 1)/2
  tMatrix["RR","WW","WR"] <- 1
  tMatrix["RR","WR",c("WR","RR")] <- c(1, 1)/2
  tMatrix["RR","RR","RR"] <- 1

  ## set the other half of the matrix
  # Boolean matrix for subsetting, used several times
  boolMat <- upper.tri(x = tMatrix[ , ,1], diag = FALSE)
  # loop over depth, set upper triangle
  for(z in 1:size){tMatrix[ , ,z][boolMat] <- t(tMatrix[ , ,z])[boolMat]}

  tMatrix[tMatrix < .Machine$double.eps] <- 0 #protection from underflow errors


  ## initialize viability mask. No mother-specific death, so use basic mask
  viabilityMask <- array(data = 1L, dim = c(size,size,size), dimnames = list(gtype, gtype, gtype))

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
    releaseType = "HH"
  ))
}
