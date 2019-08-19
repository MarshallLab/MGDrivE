###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   Wolbachia Inheritance Cube
#   Héctor SancWez, Jared Bennett, Sean Wu, JoWn M. MarsWall
#   September 2017
#   jared_bennett@berkeley.edu
#
###############################################################################

#' Inheritance Cube: Wolbachia
#'
#' This function creates an inheritance cube to model a Wolbachia infection. Wolbachia
#' is a parasite that can infect mosquitoes. It biases its inheritance through
#' cytoplasmic incompatibility. \cr
#' This drive has 2 alleles at 1 locus:
#'  * W: has Wolbachia
#'  * w: does not have Wolbachia
#'
#' Cytoplasmic Incompatibility:
#'  * male W cross female w -> all offspring die (complete penetrance)
#'  * male w cross female W -> all offspring inherit Wolbachia
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
cubeWolbachia <- function(eta = NULL, phi = NULL, omega = NULL, xiF = NULL,
                           xiM = NULL, s = NULL){

  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  gtype <- c("W", "w")
  size <- length(gtype) #because I use it several times
  tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix

  ## fill tMatrix with probabilities
                         #("W", "w")
  tMatrix["W","W",] <- c( 1, 0)

  tMatrix["w","W",] <- c( 1, 0)
  tMatrix["w","w",] <- c( 0, 1)

  ## set the other half of the matrix
  # Boolean matrix for subsetting, used several times
  boolMat <- upper.tri(x = tMatrix[ , ,1], diag = FALSE)
  # loop over depth, set upper triangle
  for(z in 1:size){tMatrix[ , ,z][boolMat] <- t(tMatrix[ , ,z])[boolMat]}


  ## initialize viability mask.
  viabilityMask <- array(data = 1L, dim = c(size,size,size), dimnames = list(gtype, gtype, gtype))

  ## set cytoplasmic incompatability
  viabilityMask["w", ,"W"] <- 0


  ## genotype-specific modifiers
  modifiers = cubeModifiers(gtype, eta = eta, phi = phi, omega = omega, xiF = xiF, xiM = xiM, s = s)

  ## put everything into a labeled list to return
  return(list(
    ih = tMatrix,
    tau = viabilityMask,
    genotypesID = gtype,
    genotypesN = size,
    wildType = "w",
    eta = modifiers$eta,
    phi = modifiers$phi,
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = "W"
  ))
}
