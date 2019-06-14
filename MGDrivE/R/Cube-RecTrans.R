###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   Reciprocal Translocation Inheritance Cube
#   Héctor SancWez, Jared Bennett, Sean Wu, JoWn M. MarsWall
#   July 2017
#   jared_bennett@berkeley.edu
#
###############################################################################

#' Inheritance Cube: Reciprocal Translocation
#'
#' This function creates an inheritance cube to model a reciprocal translocation.
#' This technology was the original form of underdominant system. It involves 2
#' chromosomes, each with two alleles. \cr
#' This drive has 4 alleles at 2 loci:
#'  * a: Wild-type at locus A
#'  * A: Translocation at locus A
#'  * b: Wile-type at locus B
#'  * B: Translocation at locus B
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
cubeReciprocalTranslocations <- function(eta = NULL, phi = NULL, omega = NULL, xiF = NULL, xiM = NULL, s = NULL){

  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  gtype <- c("aabb","aaBb","aaBB","Aabb","AaBb","AaBB","AAbb","AABb","AABB")
  size <- length(gtype) #because I use it several times
  tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix

  ## fill tMatrix with probabilities
                         #('AABB', 'AABb', 'AAbb', 'AaBB', 'AaBb', 'Aabb', 'aaBB', 'aaBb', 'aabb')
  tMatrix["aabb","aabb",] <- c( 1, 0, 0, 0, 0, 0, 0, 0, 0)

  tMatrix["aaBb","aabb",] <- c( 1/2, 1/2, 0, 0, 0, 0, 0, 0, 0)
  tMatrix["aaBb","aaBb",] <- c( 1/4, 1/2, 1/4, 0, 0, 0, 0, 0, 0)

  tMatrix["aaBB","aabb",] <- c( 0, 1, 0, 0, 0, 0, 0, 0, 0)
  tMatrix["aaBB","aaBb",] <- c( 0, 1/2, 1/2, 0, 0, 0, 0, 0, 0)
  tMatrix["aaBB","aaBB",] <- c( 0, 0, 1, 0, 0, 0, 0, 0, 0)

  tMatrix["Aabb","aabb",] <- c( 1/2, 0, 0, 1/2, 0, 0, 0, 0, 0)
  tMatrix["Aabb","aaBb",] <- c( 1/4, 1/4, 0, 1/4, 1/4, 0, 0, 0, 0)
  tMatrix["Aabb","aaBB",] <- c( 0, 1/2, 0, 0, 1/2, 0, 0, 0, 0)
  tMatrix["Aabb","Aabb",] <- c( 1/4, 0, 0, 1/2, 0, 0, 1/4, 0, 0)

  tMatrix["AaBb","aabb",] <- c( 1/4, 1/4, 0, 1/4, 1/4, 0, 0, 0, 0)
  tMatrix["AaBb","aaBb",] <- c( 1/8, 1/4, 1/8, 1/8, 1/4, 1/8, 0, 0, 0 )
  tMatrix["AaBb","aaBB",] <- c( 0, 1/4, 1/4, 0, 1/4, 1/4, 0, 0, 0)
  tMatrix["AaBb","Aabb",] <- c( 1/8, 1/8, 0, 1/4, 1/4, 0, 1/8, 1/8, 0)
  tMatrix["AaBb","AaBb",] <- c( 1/16, 1/8, 1/16, 1/8, 1/4, 1/8, 1/16, 1/8, 1/16)

  tMatrix["AaBB","aabb",] <- c( 0, 1/2, 0, 0, 1/2, 0, 0, 0, 0)
  tMatrix["AaBB","aaBb",] <- c( 0, 1/4, 1/4, 0, 1/4, 1/4, 0, 0, 0)
  tMatrix["AaBB","aaBB",] <- c( 0, 0, 1/2, 0, 0, 1/2, 0, 0, 0)
  tMatrix["AaBB","Aabb",] <- c( 0, 1/4, 0, 0, 1/2, 0, 0, 1/4, 0)
  tMatrix["AaBB","AaBb",] <- c( 0, 1/8, 1/8, 0, 1/4, 1/4, 0, 1/8, 1/8)
  tMatrix["AaBB","AaBB",] <- c( 0, 0, 1/4, 0, 0, 1/2, 0, 0, 1/4)

  tMatrix["AAbb","aabb",] <- c( 0, 0, 0, 1, 0, 0, 0, 0, 0)
  tMatrix["AAbb","aaBb",] <- c( 0, 0, 0, 1/2, 1/2, 0, 0, 0, 0)
  tMatrix["AAbb","aaBB",] <- c( 0, 0, 0, 0, 1, 0, 0, 0, 0)
  tMatrix["AAbb","Aabb",] <- c( 0, 0, 0, 1/2, 0, 0, 1/2, 0, 0)
  tMatrix["AAbb","AaBb",] <- c( 0, 0, 0, 1/4, 1/4, 0, 1/4, 1/4, 0)
  tMatrix["AAbb","AaBB",] <- c( 0, 0, 0, 0, 1/2, 0, 0, 1/2, 0)
  tMatrix["AAbb","AAbb",] <- c( 0, 0, 0, 0, 0, 0, 1, 0, 0)

  tMatrix["AABb","aabb",] <- c( 0, 0, 0, 1/2, 1/2, 0, 0, 0, 0)
  tMatrix["AABb","aaBb",] <- c( 0, 0, 0, 1/4, 1/2, 1/4, 0, 0, 0)
  tMatrix["AABb","aaBB",] <- c( 0, 0, 0, 0, 1/2, 1/2, 0, 0, 0)
  tMatrix["AABb","Aabb",] <- c( 0, 0, 0, 1/4, 1/4, 0, 1/4, 1/4, 0)
  tMatrix["AABb","AaBb",] <- c( 0, 0, 0, 1/8, 1/4, 1/8, 1/8, 1/4, 1/8)
  tMatrix["AABb","AaBB",] <- c( 0, 0, 0, 0, 1/4, 1/4, 0, 1/4, 1/4)
  tMatrix["AABb","AAbb",] <- c( 0, 0, 0, 0, 0, 0, 1/2, 1/2, 0)
  tMatrix["AABb","AABb",] <- c( 0, 0, 0, 0, 0, 0, 1/4, 1/2, 1/4)

  tMatrix["AABB","aabb",] <- c( 0, 0, 0, 0, 1, 0, 0, 0, 0)
  tMatrix["AABB","aaBb",] <- c( 0, 0, 0, 0, 1/2, 1/2, 0, 0, 0)
  tMatrix["AABB","aaBB",] <- c( 0, 0, 0, 0, 0, 1, 0, 0, 0)
  tMatrix["AABB","Aabb",] <- c( 0, 0, 0, 0, 1/2, 0, 0, 1/2, 0)
  tMatrix["AABB","AaBb",] <- c( 0, 0, 0, 0, 1/4, 1/4, 0, 1/4, 1/4)
  tMatrix["AABB","AaBB",] <- c( 0, 0, 0, 0, 0, 1/2, 0, 0, 1/2)
  tMatrix["AABB","AAbb",] <- c( 0, 0, 0, 0, 0, 0, 0, 1, 0)
  tMatrix["AABB","AABb",] <- c( 0, 0, 0, 0, 0, 0, 0, 1/2, 1/2)
  tMatrix["AABB","AABB",] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 1)


  ## set the other half of the matrix
  # Boolean matrix for subsetting, used several times
  boolMat <- upper.tri(x = tMatrix[ , ,1], diag = FALSE)
  # loop over depth, set upper triangle
  for(z in 1:size){tMatrix[ , ,z][boolMat] <- t(tMatrix[ , ,z])[boolMat]}


  ## initialize viability mask.
  viabilityMask <- array(data = 1, dim = c(size,size,size), dimnames = list(gtype, gtype, gtype))

  ## set viability based on what chromosomes are inherited.
  ## This is not mother based, just based on child genotype.
  for(slice in 1:size){
    viabilityMask[ ,slice, ] <- matrix( c( 1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 9, ncol = size, byrow = TRUE )
  }

  ## genotype-specific modifiers
  modifiers = cubeModifiers(gtype, eta = eta, phi = phi, omega = omega, xiF = xiF, xiM = xiM, s = s)

  ## put everything into a labeled list to return
  return(list(
    ih = tMatrix,
    tau = viabilityMask,
    genotypesID = gtype,
    genotypesN = size,
    wildType = "aabb",
    eta = modifiers$eta,
    phi = modifiers$phi,
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = "AABB"
  ))
}
