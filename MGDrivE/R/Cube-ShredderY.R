###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   Y-Linked X-Shredder
#   Héctor Sanchez, Jared Bennett, Sean Wu, John Marshall
#   October 2018
#   jared_bennett@berkeley.edu
#
###############################################################################

#' Inheritance Cube: Y-Linked X-Shredder
#'
#' This function creates an inheritance cube to model a Y-linked X-Shredder construct.
#' This construct resides on the Y chromosome, and chops the X chromosome into
#' many pieces during male spermatogenesis, destroying the X chromosome. Thus, males
#' only produce Y gametes. \cr
#' This drive has 5 alleles at 1 locus:
#'  * X: Wild-type X chromosome
#'  * R: X chromosome resistant to destruction by the shredder construct
#'  * Y: Wild-type Y chromosome
#'  * A: Attacking Y chromosome, a Y chromosome with the shredder construct
#'  * B: Broken Y chromosome, a Y chromosome with a defunct shredder construct
#'
#' @param cX Rate of X shredding (default is 1, complete shredding)
#' @param crX Rate of resistance chromosome generation (default is 0)
#' @param cB Rate of shredder construct breakdown (default is 0)
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
cubeXShredderY <- function(cX = 1, crX = 0, cB = 0,
                           eta = NULL, phi = NULL, omega = NULL,
                           xiF = NULL, xiM = NULL, s = NULL){

  # # test stuff
  # testVec <- runif(n = 3)
  # cX <- testVec[1]; crX <- testVec[2]; cB <- testVec[3]

  ## safety checks
  if(any(c(cX, crX, cB)>1) || any(c(cX, crX, cB)<0)){
    stop("Parameters are rates, they must be between 0 and 1.")
  }

  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  gtype <- c('XX','XR','RR','XY','XA','XB','RY','RA','RB')
  size <- length(gtype) #because I use it several times
  tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix

  ## fill tMatrix with probabilities
                                     #('XX','XR','RR','XY','XA','XB','RY','RA','RB')
  tMatrix['XX','XY',c('XX','XY')] <- c( 1, 1)/2
  tMatrix['XX','XA',c('XX','XR','XA','XB')] <- c( 1-cX, cX*crX, 1-cB, cB)/(2 + cX*(crX-1))
  tMatrix['XX','XB',c('XX','XB')] <- c( 1, 1)/2
  tMatrix['XX','RY',c('XR','XY')] <- c( 1, 1)/2
  tMatrix['XX','RA',c('XR','XA','XB')] <- c( 1, 1-cB, cB)/2
  tMatrix['XX','RB',c('XR','XB')] <- c( 1, 1)/2

  tMatrix['XR','XY',c('XX','XR','XY','RY')] <- c( 1, 1, 1, 1)/4
  tMatrix['XR','XA',c('XX','XR','XA','XB',
                      'RR','RA','RB')] <- c( 1-cX, cX*crX + 1-cX, 1-cB, cB,
                                             cX*crX, 1-cB, cB)/(4 + 2*cX*(crX-1))
  tMatrix['XR','XB',c('XX','XR','XB','RB')] <- c( 1, 1, 1, 1)/4
  tMatrix['XR','RY',c('XR','RR','XY','RY')] <- c( 1, 1, 1, 1)/4
  tMatrix['XR','RA',c('XR','XA','XB',
                      'RR','RA','RB')] <- c( 1, 1-cB, cB,
                                             1, 1-cB, cB)/4
  tMatrix['XR','RB',c('XR','RR','XB','RB')] <- c( 1, 1, 1, 1)/4

  tMatrix['RR','XY',c('XR','RY')] <- c( 1, 1)/2
  tMatrix['RR','XA',c('XR','RR','RA','RB')] <- c( 1-cX, cX*crX, 1-cB, cB)/(2 + cX*(crX-1))
  tMatrix['RR','XB',c('XR','RB')] <- c( 1, 1)/2
  tMatrix['RR','RY',c('RR','RY')] <- c( 1, 1)/2
  tMatrix['RR','RA',c('RR','RA','RB')] <- c( 1, 1-cB, cB)/2
  tMatrix['RR','RB',c('RR','RB')] <- c( 1, 1)/2

  #protection from underflow errors
  tMatrix[tMatrix < .Machine$double.eps] <- 0


  # # test stuff
  # for(fem in c('XX','XR','RR')){
  #   for(mal in c('XY','XA','XB','RY','RA','RB')){
  #     if((1-sum(tMatrix[fem,mal, ])) > sqrt(.Machine$double.eps)){
  #       print("fuck")
  #     }
  #   }
  # }


  # so, males will cut some rate of X alleles before forming gametes,
  # so all gametes are fine (assuem more sperm than necessary)
  # make sure to watch sex ratios at pupation


  ## initialize viability mask.
  viabilityMask <- array(data = 1, dim = c(size,size,size), dimnames = list(gtype, gtype, gtype))


  ## genotype-specific modifiers
  if(!is.null(phi)){
    stop("This cube has a special phi, due to being male/female specific.
         Please edit it after the cube is built, if you are absolutely sure it needs to change.")
  }
  phi = setNames(object = c(1, 1, 1, 0, 0, 0, 0, 0, 0), nm = gtype)
  modifiers = cubeModifiers(gtype, eta = eta, phi = phi, omega = omega, xiF = xiF, xiM = xiM, s = s)

  ## put everything into a labeled list to return
  return(list(
    ih = tMatrix,
    tau = viabilityMask,
    genotypesID = gtype,
    genotypesN = size,
    wildType = c('XX','XY'),
    eta = modifiers$eta,
    phi = modifiers$phi,
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = 'XA'
  ))
}
