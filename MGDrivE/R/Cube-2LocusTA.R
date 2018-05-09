###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   2-locus UDmel
#   Héctor Sanchez, Jared Bennett, Sean Wu, John M. Marshall
#   jared_bennett@berkeley.edu
#   August 2017
#
###############################################################################

#' Inheritance Cube: 2 Locus Maternal-Toxin/Zygotic-Antidote System
#'
#' This function creates a 2 locus maternal-toxin/zygotic-antidote system. This
#' is similar to the construct called UDmel. There is no resistance generation
#' in this model. \cr
#' This drive has 2 unlinked alleles, 1 allele each at 2 loci:
#'  * A: Maternal-toxin 1, zygotic-antidote 2
#'  * a: Wild-type at locus A
#'  * B: Maternal-toxin 2, zygotic-antidote 1
#'  * b: Wild-type at locus B
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
Cube_twoLocusTA <- function(TAEfficacy = 1.0, TBEfficacy = 1.0, eta = NULL,
                               phi = NULL, omega = NULL, xiF = NULL, xiM = NULL, s = NULL){

  ## safety checks in case someone is dumb
  if(any(c(TAEfficacy,TBEfficacy)>1) || any(c(TAEfficacy,TBEfficacy)<0)){
    stop("TAEfficacy, and TBEfficacy are rates.
         0 <= TAEfficacy <= 1
         0 <= TBEfficacy <= 1")
  }

  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  gtype <- c('AABB', 'AABb', 'AAbb', 'AaBB', 'AaBb', 'Aabb', 'aaBB', 'aaBb', 'aabb')
  size <- length(gtype) #because I use it several times
  tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix


  ## fill tMatrix with probabilities
                             #( 'AABB', 'AABb', 'AAbb', 'AaBB', 'AaBb', 'Aabb', aaBB', 'aaBb', 'aabb')
  tMatrix['AABB','AABB', 'AABB'] <- 1

  tMatrix['AABb','AABB', c('AABB', 'AABb')] <- c( 1/2, 1/2)
  tMatrix['AABb','AABb', c('AABB', 'AABb', 'AAbb')] <- c( 1/4, 1/2, 1/4)

  tMatrix['AAbb','AABB', 'AABb'] <- 1
  tMatrix['AAbb','AABb', c('AABb', 'AAbb')] <- c( 1/2, 1/2)
  tMatrix['AAbb','AAbb', 'AAbb'] <- 1

  tMatrix['AaBB','AABB', c('AABB', 'AaBB')] <- c( 1/2, 1/2)
  tMatrix['AaBB','AABb', c('AABB', 'AABb', 'AaBB', 'AaBb')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['AaBB','AAbb', c('AABb', 'AaBb')] <- c( 1/2, 1/2)
  tMatrix['AaBB','AaBB', c('AABB', 'AaBB', 'aaBB')] <- c( 1/4, 1/2, 1/4)

  tMatrix['AaBb','AABB', ] <- c( 1/4, 1/4, 0, 1/4, 1/4, 0, 0, 0, 0)
  tMatrix['AaBb','AABb', ] <- c( 1/8, 1/4, 1/8, 1/8, 1/4, 1/8, 0, 0, 0)
  tMatrix['AaBb','AAbb', ] <- c( 0, 1/4, 1/4, 0, 1/4, 1/4, 0, 0, 0)
  tMatrix['AaBb','AaBB', ] <- c( 1/8, 1/8, 0, 1/4, 1/4, 0, 1/8, 1/8, 0)
  tMatrix['AaBb','AaBb', ] <- c( 1/16, 1/8, 1/16, 1/8, 1/4, 1/8, 1/16, 1/8, 1/16)

  tMatrix['Aabb','AABB', c('AABb', 'AaBb')] <- c( 1/2, 1/2)
  tMatrix['Aabb','AABb', c('AABb', 'AAbb', 'AaBb', 'Aabb')] <- c( 1, 1, 1, 1)/4
  tMatrix['Aabb','AAbb', c('AAbb', 'Aabb')] <- c( 1, 1)/2
  tMatrix['Aabb','AaBB', c('AABb', 'AaBb', 'aaBb')] <- c( 1/4, 1/2, 1/4)
  tMatrix['Aabb','AaBb', ] <- c( 0, 1/8, 1/8, 0, 1/4, 1/4, 0, 1/8, 1/8)
  tMatrix['Aabb','Aabb', c('AAbb', 'Aabb', 'aabb')] <- c( 1/4, 1/2, 1/4)

  tMatrix['aaBB','AABB', 'AaBB'] <- 1
  tMatrix['aaBB','AABb', c('AaBB', 'AaBb')] <- c( 1/2, 1/2)
  tMatrix['aaBB','AAbb', 'AaBb'] <- 1
  tMatrix['aaBB','AaBB', c('AaBB', 'aaBB')] <- c( 1/2, 1/2)
  tMatrix['aaBB','AaBb', c('AaBB', 'AaBb', 'aaBB', 'aaBb')] <- c( 1, 1, 1, 1)/4
  tMatrix['aaBB','Aabb', c('AaBb', 'aaBb')] <- c( 1/2, 1/2)
  tMatrix['aaBB','aaBB', 'aaBB'] <- 1

  tMatrix['aaBb','AABB', c('AaBB', 'AaBb')] <- c( 1/2, 1/2)
  tMatrix['aaBb','AABb', c('AaBB', 'AaBb', 'Aabb')] <- c( 1/4, 1/2, 1/4)
  tMatrix['aaBb','AAbb', c('AaBb', 'Aabb')] <- c( 1/2, 1/2)
  tMatrix['aaBb','AaBB', c('AaBB', 'AaBb', 'aaBB', 'aaBb')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['aaBb','AaBb', ] <- c( 0, 0, 0, 1/2, 1, 1/2, 1/2, 1, 1/2)/4
  tMatrix['aaBb','Aabb', c('AaBb', 'Aabb', 'aaBb', 'aabb')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['aaBb','aaBB', c('aaBB', 'aaBb')] <- c( 1/2, 1/2)
  tMatrix['aaBb','aaBb', c('aaBB', 'aaBb', 'aabb')] <- c( 1/4, 1/2, 1/4)

  tMatrix['aabb','AABB', 'AaBb'] <- 1
  tMatrix['aabb','AABb', c('AaBb', 'Aabb')] <- c( 1/2, 1/2)
  tMatrix['aabb','AAbb', 'Aabb'] <- 1
  tMatrix['aabb','AaBB', c('AaBb', 'aaBb')] <- c( 1/2, 1/2)
  tMatrix['aabb','AaBb', c('AaBb', 'Aabb', 'aaBb', 'aabb')] <- c( 1, 1, 1, 1)/4
  tMatrix['aabb','Aabb', c('Aabb', 'aabb')] <- c( 1/2, 1/2)
  tMatrix['aabb','aaBB', 'aaBb'] <- 1
  tMatrix['aabb','aaBb', c('aaBb', 'aabb')] <- c( 1/2, 1/2)
  tMatrix['aabb','aabb', 'aabb'] <- 1

  ## set the other half of the matrix
  SymCubeC(lowerMat = tMatrix)


  ## initialize viability mask.
  viabilityMask <- matrix(data = 1, nrow = size, ncol = size, dimnames = list(gtype, gtype))

  ## fill mother/offspring specific death, then muliply by efficacy of toxins
  mixed <- (1-TAEfficacy)+(1-TBEfficacy)-(1-TAEfficacy)*(1-TBEfficacy)
  viabilityMask[c('AABB','AABb', 'AaBB', 'AaBb'), ] <- matrix( c( 1, 1, 1-TAEfficacy, 1, 1, 1-TAEfficacy, 1-TBEfficacy, 1-TBEfficacy, mixed), nrow = 4, ncol = size, byrow = TRUE )
  viabilityMask[c('AAbb','Aabb'), ] <- matrix( c( 1, 1, 1-TAEfficacy, 1, 1, 1-TAEfficacy, 1, 1, 1-TAEfficacy), nrow = 2, ncol = size, byrow = TRUE )
  viabilityMask[c('aaBB','aaBb'), ] <- matrix( c( 1, 1, 1, 1, 1, 1, 1-TBEfficacy, 1-TBEfficacy, 1-TBEfficacy), nrow = 2, ncol = size, byrow = TRUE )


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
    releaseType="AABB"
  ))
}
