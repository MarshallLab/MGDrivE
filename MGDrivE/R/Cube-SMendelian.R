###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   Mendelian
#   Héctor Sanchez, Jared Bennett, Sean Wu, John Marshall
#   jared_bennett@berkeley.edu
#   August 2017
#
###############################################################################

#' Inheritance Cube: Mendelian
#'
#' This function creates a Mendelian Inheritance Cube. It only handles simple,
#' alphabetic genotypes. \cr
#' The default is 3 alleles at 1 locus, this can be extended to however many
#' alleles one is interested in, but only at 1 locus.
#'
#' @param gtype Vector of genotypes, with the wild-type in the first position
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
cubeMendelian <- function(gtype = c("AA", "Aa", "aa"), eta = NULL, phi = NULL,
                                omega = NULL, xiF = NULL, xiM = NULL, s = NULL){

  ## safety check
  if(!all(nchar(gtype[1])==nchar(gtype))){
    stop("All the genotypes are not the same length")
  }


  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  size <- length(gtype) #because I use it several times
  tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix
  testVec <- setNames(object = numeric(size), nm = gtype)  #need later

  ## fill tMatrix with probabilities
  for (i in gtype)  # loop over female genotypes
  {
    for (j in gtype)  # loop over male genotypes
    {

      male <- strsplit(j, split='')[[1]] # male genotype for this cross
      female <- strsplit(i, split='')[[1]] # female genotype for this cross
      offspring <- as.vector( outer(male, female, paste0, sep='') )  # offspring genotypes for this cross

      # reorder all offspring alleles to match allele order in gtype
      offspring <- vapply( strsplit(offspring, split=''),
                           function(x) {paste0(sort(x, method = 'radix'), collapse='')},
                           FUN.VALUE = character(1))

      # count genotypes of offspring, order according to gtype
      for (k in offspring)
      {
        testVec[k] <- testVec[k]+1
      }

      testVec[] <- testVec/sum(testVec) # normalize offspring frequencies to 1

      tMatrix[i,j, ] <- testVec  # store offspring frequences in transition matrix

      testVec[] <- 0  # clear testVec
    }
  }


  ## initialize viability mask. No mother-specific death, so use basic mask
  viabilityMask <- array(data = 1, dim = c(size,size,size),
                         dimnames = list(gtype, gtype, gtype))

  ## genotype-specific modifiers
  modifiers = cubeModifiers(gtype, eta = eta, phi = phi, omega = omega,
                            xiF = xiF, xiM = xiM, s = s)

  ## put everything into a labeled list to return
  return(list(
    ih = tMatrix,
    tau = viabilityMask,
    genotypesID = gtype,
    genotypesN = size,
    wildType = gtype[[1]],
    eta = modifiers$eta,
    phi = modifiers$phi,
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = gtype[size]
  ))
}
