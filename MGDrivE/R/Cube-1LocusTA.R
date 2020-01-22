###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   1-locus UDmel
#   Héctor Sanchez, Jared Bennett, Sean Wu, John Marshall
#   jared_bennett@berkeley.edu
#   August 2017
#
###############################################################################

#' Inheritance Cube: 1 Locus Maternal-Toxin/Zygotic-Antidote System
#'
#' This function creates a 1 locus maternal-toxin/zygotic-antidote system. This
#' is similar to the construct called UDmel. There is no resistance generation
#' in this model. \cr
#' This drive has 3 alleles at 1 locus:
#'  * A: Maternal-toxin 1, zygotic-antidote 2
#'  * B: Maternal-toxin 2, zygotic-antidote 1
#'  * W: Wild-type allele
#'
#' @param TAEfficacy Maternal toxin A efficacy
#' @param TBEfficacy Maternal toxin B efficacy
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
cubeOneLocusTA <- function(TAEfficacy = 1.0, TBEfficacy = 1.0,
                                    eta = NULL, phi = NULL, omega = NULL,
                                    xiF = NULL, xiM = NULL, s = NULL){

  ## safety checks
  if(any(c(TAEfficacy,TBEfficacy)>1) || any(c(TAEfficacy,TBEfficacy)<0)){
    stop("TAEfficacy, and TBEfficacy are rates.
         0 <= TAEfficacy <= 1
         0 <= TBEfficacy <= 1")
  }

  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  gtype <- c('WW', 'WA', 'WB', 'AA', 'AB', 'BB')
  size <- length(gtype) #because I use it several times
  tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix


  ## fill tMatrix with probabilities
                         #( 'WW', 'WA', 'WB', 'AA', 'AB', 'BB')
  tMatrix['WW','WW', 'WW'] <- 1

  tMatrix['WA','WW', c('WW', 'WA')] <- c( 1/2, 1/2)
  tMatrix['WA','WA', c('WW', 'WA', 'AA')] <- c( 1/4, 1/2, 1/4)

  tMatrix['WB','WW', c('WW', 'WB')] <- c( 1/2, 1/2)
  tMatrix['WB','WA', c('WW', 'WA', 'WB', 'AB')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['WB','WB', c('WW', 'WB', 'BB')] <- c( 1/4, 1/2, 1/4)

  tMatrix['AA', 'WW', 'WA'] <- 1
  tMatrix['AA', 'WA', c('WA', 'AA')] <- c( 1/2, 1/2)
  tMatrix['AA', 'WB', c('WA', 'AB')] <- c( 1/2, 1/2)
  tMatrix['AA', 'AA', 'AA'] <- 1

  tMatrix['AB','WW', c('WA', 'WB')] <- c( 1/2, 1/2)
  tMatrix['AB','WA', c('WA', 'WB', 'AA', 'AB')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['AB','WB', c('WA', 'WB', 'AB', 'BB')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['AB','AA', c('AA', 'AB')] <- c( 1/2, 1/2)
  tMatrix['AB','AB', c('AA', 'AB', 'BB')] <- c( 1/4, 1/2, 1/4)

  tMatrix['BB','WW', 'WB'] <- 1
  tMatrix['BB','WA', c('WB', 'AB')] <- c( 1/2, 1/2)
  tMatrix['BB','WB', c('WB', 'BB')] <- c( 1/2, 1/2)
  tMatrix['BB','AA', 'AB'] <- 1
  tMatrix['BB','AB', c('AB', 'BB')] <- c( 1/2, 1/2)
  tMatrix['BB','BB', 'BB'] <- 1

  ## set the other half of the matrix
  # Boolean matrix for subsetting, used several times
  boolMat <- upper.tri(x = tMatrix[ , ,1], diag = FALSE)
  # loop over depth, set upper triangle
  for(z in 1:size){tMatrix[ , ,z][boolMat] <- t(tMatrix[ , ,z])[boolMat]}


  ## initialize viability mask.
  viabilityMask <- array(data = 1L, dim = c(size,size,size), dimnames = list(gtype, gtype, gtype))

  ## fill mother/offspring specific death, then muliply by efficacy of toxins
  mixed <- (1-TAEfficacy)+(1-TBEfficacy)-(1-TAEfficacy)*(1-TBEfficacy)
  viabilityMask[c('WA','AA'), , ] <- matrix( c( 1-TAEfficacy, 1-TAEfficacy, 1, 1-TAEfficacy, 1, 1), nrow = 2, ncol = size, byrow = TRUE )
  viabilityMask[c('WB','BB'), , ] <- matrix( c( 1-TBEfficacy, 1, 1-TBEfficacy, 1, 1, 1-TBEfficacy), nrow = 2, ncol = size, byrow = TRUE )
  viabilityMask['AB', , ] <- c( mixed, 1-TAEfficacy, 1-TBEfficacy, 1-TAEfficacy, 1, 1-TBEfficacy)

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
    releaseType = "AB"
  ))
}
