###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   RIDL
#   Héctor Sanchez, Jared Bennett, Sean Wu, John Marshall
#   August 2017
#   Jared_Bennett@berkeley.edu
#
###############################################################################

#' Inheritance Cube: RIDL (Release of Insects with Dominant Lethality)
#'
#' This function creates a RIDL system.
#' RIDL (Release of Insects with Dominant Lethality), is a form of SIT.
#' Created by Oxitec, this is based on a positive feedback loop using the
#' toxic tTAV gene, controlled under lab conditions by the TetO promoter.
#' This has 2 alleles at 1 locus
#'  * W: Wild-type allele
#'  * R: OX513 RIDL allele
#'
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
cubeRIDL <- function(eta = NULL, phi = NULL, omega = NULL, xiF = NULL, xiM = NULL, s = NULL){

  gtype <- c('WW', 'WR', 'RR')
  size <- length(gtype) #because I use it several times
  tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix

  ## fill tMatrix with probabilities
                          #( 'rr', 'rR', 'RR')
  tMatrix['WW','WW', ] <- c( 1, 0, 0)

  tMatrix['WR','WW', ] <- c( 1/2, 1/2, 0)
  tMatrix['WR','WR', ] <- c( 1/4, 1/2, 1/4)

  tMatrix['RR','WW', ] <- c( 0, 1, 0)
  tMatrix['RR','WR', ] <- c( 0, 1/2, 1/2)
  tMatrix['RR','RR', ] <- c( 0, 0, 1)

  ## set the other half of the matrix
  # Boolean matrix for subsetting, used several times
  boolMat <- upper.tri(x = tMatrix[ , ,1], diag = FALSE)
  # loop over depth, set upper triangle
  for(z in 1:size){tMatrix[ , ,z][boolMat] <- t(tMatrix[ , ,z])[boolMat]}

  ## initialize viability mask. No mother-specific death, so use basic mask
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
    releaseType = "RR"
    ))

}
