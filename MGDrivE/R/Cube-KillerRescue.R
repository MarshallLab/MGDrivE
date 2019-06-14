###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   Killer-Rescue System
#   Héctor Sanchez, Jared Bennett, Sean Wu, John M. Marshall
#   jared_bennett@berkeley.edu
#   August 2017
#
###############################################################################

#' Inheritance Cube: Killer-Rescue System
#'
#' This function creates an inheritance cube to model a Killer-Rescue system.
#' Killer-Resuce is a 2-locus system: one locus has a toxin and the other locus contains
#' the antidote. The loci are assumed independent and are non-homing. \cr
#' This drive has 3 alleles at locus 1 and 2 alleles and locus 2:
#'  * Locus 1
#'    * T: Wild-type allele
#'    * K: "Killer" toxin allele
#'    * R: Rroken toxin allele
#'  * Locus 2
#'    * W: Wild-type allele
#'    * A: Antidote allele
#'
#' @param eR Conversion of K allele to R allele, a basal mutation rate
#' @param Keff Toxin efficacy
#' @param Aeff Antidote efficacy
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
cubeKillerRescue <- function(eR = 0, Keff = 1.0, Aeff = 1.0, eta = NULL,
                              phi = NULL, omega = NULL, xiF = NULL, xiM = NULL, s = NULL){

  #safety checks
  if(any(c(eR, Keff, Aeff)<0) || any(c(eR, Keff, Aeff)>1)){
    stop("Parameters are rates, they must be between 0 and 1.")
  }

  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  gtype <- c('TTWW', 'TTWA', 'TTAA', 'TKWW', 'TKWA', 'TKAA', 'TRWW', 'TRWA', 'TRAA',
             'KKWW', 'KKWA', 'KKAA', 'KRWW', 'KRWA', 'KRAA', 'RRWW', 'RRWA', 'RRAA')
  size <- length(gtype) #because I use it several times
  tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix

  ## fill tMatrix with probabilities
                             ##('TTWW', 'TTWA', 'TTAA', 'TKWW', 'TKWA', 'TKAA', 'TRWW', 'TRWA', 'TRAA',
                             ## 'KKWW', 'KKWA', 'KKAA', 'KRWW', 'KRWA', 'KRAA', 'RRWW', 'RRWA', 'RRAA')
  tMatrix['TTWW','TTWW', ] <- c( 1, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)

  tMatrix['TTWA','TTWW', ] <- c( 1/2, 1/2, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['TTWA','TTWA', ] <- c( 1/4, 1/2, 1/4, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)

  tMatrix['TTAA','TTWW', ] <- c( 0, 1, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['TTAA','TTWA', ] <- c( 0, 1/2, 1/2, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['TTAA','TTAA', ] <- c( 0, 0, 1, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)

  tMatrix['TKWW','TTWW', ] <- c( 1/2, 0, 0, 1/2-eR, 0, 0, eR, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['TKWW','TTWA', ] <- c( 1/4, 1/4, 0, (1/2-eR)/2, (1/2-eR)/2, 0, eR/2, eR/2, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['TKWW','TTAA', ] <- c( 0, 1/2, 0, 0, 1/2-eR, 0, 0, eR, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['TKWW','TKWW', ] <- c( 1/4, 0, 0, 1/2-eR, 0, 0, eR, 0, 0,
                                 (1/2-eR)^2, 0, 0, 2*eR*(1/2-eR), 0, 0, eR^2, 0, 0)

  tMatrix['TKWA','TTWW', ] <- c( 1/4, 1/4, 0, (1/2-eR)/2, (1/2-eR)/2, 0, eR/2, eR/2, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['TKWA','TTWA', ] <- c( 1/8, 1/4, 1/8, (1/2-eR)/4, (1/2-eR)/2, (1/2-eR)/4, eR/4, eR/2, eR/4,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['TKWA','TTAA', ] <- c( 0, 1/4, 1/4, 0, (1/2-eR)/2, (1/2-eR)/2, 0, eR/2, eR/2,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['TKWA','TKWW', ] <- c( 1/8, 1/8, 0, (1/2-eR)/2, (1/2-eR)/2, 0, eR/2, eR/2, 0,
                                 (1/2-eR)^2/2, (1/2-eR)^2/2, 0, eR*(1/2-eR), eR*(1/2-eR), 0, eR^2/2, eR^2/2, 0)
  tMatrix['TKWA','TKWA', ] <- c( 1/16, 1/8, 1/16, (1/2-eR)/4, (1/2-eR)/2, (1/2-eR)/4, eR/4, eR/2, eR/4,
                                 (1/2-eR)^2/4, (1/2-eR)^2/2, (1/2-eR)^2/4, eR*(1/2-eR)/2, eR*(1/2-eR), eR*(1/2-eR)/2, eR^2/4, eR^2/2, eR^2/4)

  tMatrix['TKAA','TTWW', ] <- c( 0, 1/2, 0, 0, (1/2-eR), 0, 0, eR, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['TKAA','TTWA', ] <- c( 0, 1/4, 1/4, 0, (1/2-eR)/2, (1/2-eR)/2, 0, eR/2, eR/2,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['TKAA','TTAA', ] <- c( 0, 0, 1/2, 0, 0, (1/2-eR), 0, 0, eR,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['TKAA','TKWW', ] <- c( 0, 1/4, 0, 0, (1/2-eR), 0, 0, eR, 0,
                                 0, (1/2-eR)^2, 0, 0, 2*eR*(1/2-eR), 0, 0, eR^2, 0)
  tMatrix['TKAA','TKWA', ] <- c( 0, 1/8, 1/8, 0, (1/2-eR)/2, (1/2-eR)/2, 0, eR/2, eR/2,
                                 0, (1/2-eR)^2/2, (1/2-eR)^2/2, 0, eR*(1/2-eR), eR*(1/2-eR), 0, eR^2/2, eR^2/2)
  tMatrix['TKAA','TKAA', ] <- c( 0, 0, 1/4, 0, 0, (1/2-eR), 0, 0, eR,
                                 0, 0, (1/2-eR)^2, 0, 0, 2*eR*(1/2-eR), 0, 0, eR^2)

  tMatrix['TRWW','TTWW', ] <- c( 1/2, 0, 0, 0, 0, 0, 1/2, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['TRWW','TTWA', ] <- c( 1/4, 1/4, 0, 0, 0, 0, 1/4, 1/4, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['TRWW','TTAA', ] <- c( 0, 1/2, 0, 0, 0, 0, 0, 1/2, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['TRWW','TKWW', ] <- c( 1/4, 0, 0, (1/2-eR)/2, 0, 0, (1/2+eR)/2, 0, 0,
                                 0, 0, 0, (1/2-eR)/2, 0, 0, eR/2, 0, 0)
  tMatrix['TRWW','TKWA', ] <- c( 1/8, 1/8, 0, (1/2-eR)/4, (1/2-eR)/4, 0, (1/2+eR)/4, (1/2+eR)/4, 0,
                                 0, 0, 0, (1/2-eR)/4, (1/2-eR)/4, 0, eR/4, eR/4, 0)
  tMatrix['TRWW','TKAA', ] <- c( 0, 1/4, 0, 0, (1/2-eR)/2, 0, 0, (1/2+eR)/2, 0,
                                 0, 0, 0, 0, (1/2-eR)/2, 0, 0, eR/2, 0)
  tMatrix['TRWW','TRWW', ] <- c( 1/4, 0, 0, 0, 0, 0, 1/2, 0, 0,
                                 0, 0, 0, 0, 0, 0, 1/4, 0, 0)

  tMatrix['TRWA','TTWW', ] <- c( 1/4, 1/4, 0, 0, 0, 0, 1/4, 1/4, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['TRWA','TTWA', ] <- c( 1/8, 1/4, 1/8, 0, 0, 0, 1/8, 1/4, 1/8,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['TRWA','TTAA', ] <- c( 0, 1/4, 1/4, 0, 0, 0, 0, 1/4, 1/4,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['TRWA','TKWW', ] <- c( 1/8, 1/8, 0, (1/2-eR)/4, (1/2-eR)/4, 0, (1/2+eR)/4, (1/2+eR)/4, 0,
                                 0, 0, 0, (1/2-eR)/4, (1/2-eR)/4, 0, eR/4, eR/4, 0)
  tMatrix['TRWA','TKWA', ] <- c( 1/16, 1/8, 1/16, (1/2-eR)/8, (1/2-eR)/4, (1/2-eR)/8, (1/2+eR)/8, (1/2+eR)/4, (1/2+eR)/8,
                                 0, 0, 0, (1/2-eR)/8, (1/2-eR)/4, (1/2-eR)/8, eR/8, eR/4, eR/8)
  tMatrix['TRWA','TKAA', ] <- c( 0, 1/8, 1/8, 0, (1/2-eR)/4, (1/2-eR)/4, 0, (1/2+eR)/4, (1/2+eR)/4,
                                 0, 0, 0, 0, (1/2-eR)/4, (1/2-eR)/4, 0, eR/4, eR/4)
  tMatrix['TRWA','TRWW', ] <- c( 1/8, 1/8, 0, 0, 0, 0, 1/4, 1/4, 0,
                                 0, 0, 0, 0, 0, 0, 1/8, 1/8,0)
  tMatrix['TRWA','TRWA', ] <- c( 1/16, 1/8, 1/16, 0, 0, 0, 1/8, 1/4, 1/8,
                                 0, 0, 0, 0, 0, 0, 1/16, 1/8, 1/16)

  tMatrix['TRAA','TTWW', ] <- c( 1/2, 0, 0, 0, 0, 0, 1/2, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['TRAA','TTWA', ] <- c( 0, 1/4, 1/4, 0, 0, 0, 0, 1/4, 1/4,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['TRAA','TTAA', ] <- c( 0, 0, 1/2, 0, 0, 0, 0, 0, 1/2,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['TRAA','TKWW', ] <- c( 0, 1/4, 0, 0, (1/2-eR)/2, 0, 0, (1/2+eR)/2, 0,
                                 0, 0, 0, 0, (1/2-eR)/2, 0, 0, eR/2, 0)
  tMatrix['TRAA','TKWA', ] <- c( 0, 1/8, 1/8, 0, (1/2-eR)/4, (1/2-eR)/4, 0, (1/2+eR)/4, (1/2+eR)/4,
                                 0, 0, 0, 0, (1/2-eR)/4, (1/2-eR)/4, 0, eR/4, eR/4)
  tMatrix['TRAA','TKAA', ] <- c( 0, 0, 1/4, 0, 0, (1/2-eR)/2, 0, 0, (1/2+eR)/2,
                                 0, 0, 0, 0, 0, (1/2-eR)/2, 0, 0, eR/2)
  tMatrix['TRAA','TRWW', ] <- c( 0, 1/4, 0, 0, 0, 0, 0, 1/2, 0,
                                 0, 0, 0, 0, 0, 0, 0, 1/4, 0)
  tMatrix['TRAA','TRWA', ] <- c( 0, 1/8, 1/8, 0, 0, 0, 0, 1/4, 1/4,
                                 0, 0, 0, 0, 0, 0, 0, 1/8, 1/8)
  tMatrix['TRAA','TRAA', ] <- c( 0, 0, 1/4, 0, 0, 0, 0, 0, 1/2,
                                 0, 0, 0, 0, 0, 0, 0, 0, 1/4)

  tMatrix['KKWW','TTWW', ] <- c( 0, 0, 0, 2*(1/2-eR), 0, 0, 2*eR, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['KKWW','TTWA', ] <- c( 0, 0, 0, (1/2-eR), (1/2-eR), 0, eR, eR, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['KKWW','TTAA', ] <- c( 0, 0, 0, 0, 2*(1/2-eR), 0, 0, 2*eR, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['KKWW','TKWW', ] <- c( 0, 0, 0, (1/2-eR), 0, 0, eR, 0, 0,
                                 2*(1/2-eR)^2, 0, 0, 4*eR*(1/2-eR), 0, 0, 2*eR^2, 0, 0)
  tMatrix['KKWW','TKWA', ] <- c( 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, eR/2, eR/2, 0,
                                 (1/2-eR)^2, (1/2-eR)^2, 0, 2*eR*(1/2-eR), 2*eR*(1/2-eR), 0, eR^2, eR^2, 0)
  tMatrix['KKWW','TKAA', ] <- c( 0, 0, 0, 0, (1/2-eR), 0, 0, eR, 0,
                                 0, 2*(1/2-eR)^2, 0, 0, 4*eR*(1/2-eR), 0, 0, 2*eR^2, 0)
  tMatrix['KKWW','TRWW', ] <- c( 0, 0, 0, (1/2-eR), 0, 0, eR, 0, 0,
                                 0, 0, 0, (1/2-eR), 0, 0, eR, 0, 0)
  tMatrix['KKWW','TRWA', ] <- c( 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, eR/2, eR/2, 0,
                                 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, eR/2, eR/2, 0)
  tMatrix['KKWW','TRAA', ] <- c( 0, 0, 0, 0, (1/2-eR), 0, 0, eR, 0,
                                 0, 0, 0, 0, (1/2-eR), 0, 0, eR, 0)
  tMatrix['KKWW','KKWW', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 4*(1/2-eR)^2, 0, 0, 8*eR*(1/2-eR), 0, 0, 4*eR^2, 0, 0)

  tMatrix['KKWA','TTWW', ] <- c( 0, 0, 0, (1/2-eR), (1/2-eR), 0, eR, eR, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['KKWA','TTWA', ] <- c( 0, 0, 0, (1/2-eR)/2, (1/2-eR), (1/2-eR)/2, eR/2, eR, eR/2,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['KKWA','TTAA', ] <- c( 0, 0, 0, 0, (1/2-eR), (1/2-eR), 0, eR, eR,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['KKWA','TKWW', ] <- c( 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, eR/2, eR/2, 0,
                                 (1/2-eR)^2, (1/2-eR)^2, 0, 2*eR*(1/2-eR), 2*eR*(1/2-eR), 0, eR^2, eR^2, 0)
  tMatrix['KKWA','TKWA', ] <- c( 0, 0, 0, (1/2-eR)/4, (1/2-eR)/2, (1/2-eR)/4, eR/4, eR/2, eR/4,
                                 (1/2-eR)^2/2, (1/2-eR)^2, (1/2-eR)^2/2, eR*(1/2-eR), 2*eR*(1/2-eR), eR*(1/2-eR), eR^2/2, eR^2, eR^2/2)
  tMatrix['KKWA','TKAA', ] <- c( 0, 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, eR/2, eR/2,
                                 0, (1/2-eR)^2, (1/2-eR)^2, 0, 2*eR*(1/2-eR), 2*eR*(1/2-eR), 0, eR^2, eR^2)
  tMatrix['KKWA','TRWW', ] <- c( 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, eR/2, eR/2, 0,
                                 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, eR/2, eR/2, 0)
  tMatrix['KKWA','TRWA', ] <- c( 0, 0, 0, (1/2-eR)/4, (1/2-eR)/2, (1/2-eR)/4, eR/4, eR/2, eR/4,
                                 0, 0, 0, (1/2-eR)/4, (1/2-eR)/2, (1/2-eR)/4, eR/4, eR/2, eR/4)
  tMatrix['KKWA','TRAA', ] <- c( 0, 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, eR/2, eR/2,
                                 0, 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, eR/2, eR/2)
  tMatrix['KKWA','KKWW', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 2*(1/2-eR)^2, 2*(1/2-eR)^2, 0, 4*eR*(1/2-eR), 4*eR*(1/2-eR), 0, 2*eR^2, 2*eR^2, 0)
  tMatrix['KKWA','KKWA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 (1/2-eR)^2, 2*(1/2-eR)^2, (1/2-eR)^2, 2*eR*(1/2-eR), 4*eR*(1/2-eR), 2*eR*(1/2-eR), eR^2, 2*eR^2, eR^2)

  tMatrix['KKAA','TTWW', ] <- c( 0, 0, 0, 0, 2*(1/2-eR), 0, 0, 2*eR, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['KKAA','TTWA', ] <- c( 0, 0, 0, 0, (1/2-eR), (1/2-eR), 0, eR, eR,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['KKAA','TTAA', ] <- c( 0, 0, 0, 0, 0, 2*(1/2-eR), 0, 0, 2*eR,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['KKAA','TKWW', ] <- c( 0, 0, 0, 0, (1/2-eR), 0, 0, eR, 0,
                                 0, 2*(1/2-eR)^2, 0, 0, 4*eR*(1/2-eR), 0, 0, 2*eR^2, 0)
  tMatrix['KKAA','TKWA', ] <- c( 0, 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, eR/2, eR/2,
                                 0, (1/2-eR)^2, (1/2-eR)^2, 0, 2*eR*(1/2-eR), 2*eR*(1/2-eR), 0, eR^2, eR^2)
  tMatrix['KKAA','TKAA', ] <- c( 0, 0, 0, 0, 0, (1/2-eR), 0, 0, eR,
                                 0, 0, 2*(1/2-eR)^2, 0, 0, 4*eR*(1/2-eR), 0, 0, 2*eR^2)
  tMatrix['KKAA','TRWW', ] <- c( 0, 0, 0, 0, (1/2-eR), 0, 0, eR, 0,
                                 0, 0, 0, 0, (1/2-eR), 0, 0, eR, 0)
  tMatrix['KKAA','TRWA', ] <- c( 0, 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, eR/2, eR/2,
                                 0, 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, eR/2, eR/2)
  tMatrix['KKAA','TRAA', ] <- c( 0, 0, 0, 0, 0, (1/2-eR), 0, 0, eR,
                                 0, 0, 0, 0, 0, (1/2-eR), 0, 0, eR)
  tMatrix['KKAA','KKWW', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 4*(1/2-eR)^2, 0, 0, 8*eR*(1/2-eR), 0, 0, 4*eR^2, 0)
  tMatrix['KKAA','KKWA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 2*(1/2-eR)^2, 2*(1/2-eR)^2, 0, 4*eR*(1/2-eR), 4*eR*(1/2-eR), 0, 2*eR^2, 2*eR^2)
  tMatrix['KKAA','KKAA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 4*(1/2-eR)^2, 0, 0, 8*eR*(1/2-eR), 0, 0, 4*eR^2)

  tMatrix['KRWW','TTWW', ] <- c( 0, 0, 0, 1/2-eR, 0, 0, 1/2+eR, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['KRWW','TTWA', ] <- c( 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, (1/2+eR)/2, (1/2+eR)/2, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['KRWW','TTAA', ] <- c( 0, 0, 0, 0, 1/2-eR, 0, 0, 1/2+eR, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['KRWW','TKWW', ] <- c( 0, 0, 0, (1/2-eR)/2, 0, 0, (1/2+eR)/2, 0, 0,
                                 (1/2-eR)^2, 0, 0, (1/2-eR)*(1/2+2*eR), 0, 0, eR*(1/2+eR), 0, 0)
  tMatrix['KRWW','TKWA', ] <- c( 0, 0, 0, (1/2-eR)/4, (1/2-eR)/4, 0, (1/2+eR)/4, (1/2+eR)/4, 0,
                                 (1/2-eR)^2/2, (1/2-eR)^2/2, 0, (1/2-eR)*(1/2+2*eR)/2, (1/2-eR)*(1/2+2*eR)/2, 0, eR*(1/2+eR)/2, eR*(1/2+eR)/2, 0)
  tMatrix['KRWW','TKAA', ] <- c( 0, 0, 0, 0, (1/2-eR)/2, 0, 0, (1/2+eR)/2, 0,
                                 0, (1/2-eR)^2, 0, 0, (1/2-eR)*(1/2+2*eR), 0, 0, eR*(1/2+eR), 0)
  tMatrix['KRWW','TRWW', ] <- c( 0, 0, 0, (1/2-eR)/2, 0, 0, (1/2+eR)/2, 0, 0,
                                 0, 0, 0, (1/2-eR)/2, 0, 0, (1/2+eR)/2, 0, 0)
  tMatrix['KRWW','TRWA', ] <- c( 0, 0, 0, (1/2-eR)/4, (1/2-eR)/4, 0, (1/2+eR)/4, (1/2+eR)/4, 0,
                                 0, 0, 0, (1/2-eR)/4, (1/2-eR)/4, 0, (1/2+eR)/4, (1/2+eR)/4, 0)
  tMatrix['KRWW','TRAA', ] <- c( 0, 0, 0, 0, (1/2-eR)/2, 0, 0, (1/2+eR)/2, 0,
                                 0, 0, 0, 0, (1/2-eR)/2, 0, 0, (1/2+eR)/2, 0)
  tMatrix['KRWW','KKWW', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 2*(1/2-eR)^2, 0, 0, 2*(1/2-eR)*(1/2+2*eR), 0, 0, 2*eR*(1/2+eR), 0, 0)
  tMatrix['KRWW','KKWA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 (1/2-eR)^2, (1/2-eR)^2, 0, (1/2-eR)*(1/2+2*eR), (1/2-eR)*(1/2+2*eR), 0, eR*(1/2+eR), eR*(1/2+eR), 0)
  tMatrix['KRWW','KKAA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 2*(1/2-eR)^2, 0, 0, 2*(1/2-eR)*(1/2+2*eR), 0, 0, 2*eR*(1/2+eR), 0)
  tMatrix['KRWW','KRWW', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 (1/2-eR)^2, 0, 0, 2*(1/2-eR)*(1/2+eR), 0, 0, (1/2+eR)^2, 0, 0)

  tMatrix['KRWA','TTWW', ] <- c( 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, (1/2+eR)/2, (1/2+eR)/2, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['KRWA','TTWA', ] <- c( 0, 0, 0, (1/2-eR)/4, (1/2-eR)/2, (1/2-eR)/4, (1/2+eR)/4, (1/2+eR)/2, (1/2+eR)/4,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['KRWA','TTAA', ] <- c( 0, 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, (1/2+eR)/2, (1/2+eR)/2,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['KRWA','TKWW', ] <- c( 0, 0, 0, (1/2-eR)/4, (1/2-eR)/4, 0, (1/2+eR)/4, (1/2+eR)/4, 0,
                                 (1/2-eR)^2/2, (1/2-eR)^2/2, 0, (1/2-eR)*(1/2+2*eR)/2, (1/2-eR)*(1/2+2*eR)/2, 0, eR*(1/2+eR)/2, eR*(1/2+eR)/2, 0)
  tMatrix['KRWA','TKWA', ] <- c( 0, 0, 0, (1/2-eR)/8, (1/2-eR)/4, (1/2-eR)/8, (1/2+eR)/8, (1/2+eR)/4, (1/2+eR)/8,
                                 (1/2-eR)^2/4, (1/2-eR)^2/2, (1/2-eR)^2/4, (1/2-eR)*(1/2+2*eR)/4, (1/2-eR)*(1/2+2*eR)/2, (1/2-eR)*(1/2+2*eR)/4, eR*(1/2+eR)/4, eR*(1/2+eR)/2, eR*(1/2+eR)/4)
  tMatrix['KRWA','TKAA', ] <- c( 0, 0, 0, 0, (1/2-eR)/4, (1/2-eR)/4, 0, (1/2+eR)/4, (1/2+eR)/4,
                                 0, (1/2-eR)^2/2, (1/2-eR)^2/2, 0, (1/2-eR)*(1/2+2*eR)/2, (1/2-eR)*(1/2+2*eR)/2, 0, eR*(1/2+eR)/2, eR*(1/2+eR)/2)
  tMatrix['KRWA','TRWW', ] <- c( 0, 0, 0, (1/2-eR)/4, (1/2-eR)/4, 0, (1/2+eR)/4, (1/2+eR)/4, 0,
                                 0, 0, 0, (1/2-eR)/4, (1/2-eR)/4, 0, (1/2+eR)/4, (1/2+eR)/4, 0)
  tMatrix['KRWA','TRWA', ] <- c( 0, 0, 0, (1/2-eR)/8, (1/2-eR)/4, (1/2-eR)/8, (1/2+eR)/8, (1/2+eR)/4, (1/2+eR)/8,
                                 0, 0, 0, (1/2-eR)/8, (1/2-eR)/4, (1/2-eR)/8, (1/2+eR)/8, (1/2+eR)/4, (1/2+eR)/8)
  tMatrix['KRWA','TRAA', ] <- c( 0, 0, 0, 0, (1/2-eR)/4, (1/2-eR)/4, 0, (1/2+eR)/4, (1/2+eR)/4,
                                 0, 0, 0, 0, (1/2-eR)/4, (1/2-eR)/4, 0, (1/2+eR)/4, (1/2+eR)/4)
  tMatrix['KRWA','KKWW', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 (1/2-eR)^2, (1/2-eR)^2, 0, (1/2-eR)*(1/2+2*eR), (1/2-eR)*(1/2+2*eR), 0, eR*(1/2+eR), eR*(1/2+eR), 0)
  tMatrix['KRWA','KKWA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 (1/2-eR)^2/2, (1/2-eR)^2, (1/2-eR)^2/2, (1/2-eR)*(1/2+2*eR)/2, (1/2-eR)*(1/2+2*eR), (1/2-eR)*(1/2+2*eR)/2, eR*(1/2+eR)/2, eR*(1/2+eR), eR*(1/2+eR)/2)
  tMatrix['KRWA','KKAA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, (1/2-eR)^2, (1/2-eR)^2, 0, (1/2-eR)*(1/2+2*eR), (1/2-eR)*(1/2+2*eR), 0, eR*(1/2+eR), eR*(1/2+eR))
  tMatrix['KRWA','KRWW', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 (1/2-eR)^2/2, (1/2-eR)^2/2, 0, (1/2-eR)*(1/2+eR), (1/2-eR)*(1/2+eR), 0, (1/2+eR)^2/2, (1/2+eR)^2/2, 0)
  tMatrix['KRWA','KRWA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 (1/2-eR)^2/4, (1/2-eR)^2/2, (1/2-eR)^2/4, (1/2-eR)*(1/2+eR)/2, (1/2-eR)*(1/2+eR), (1/2-eR)*(1/2+eR)/2, (1/2+eR)^2/4, (1/2+eR)^2/2, (1/2+eR)^2/4)

  tMatrix['KRAA','TTWW', ] <- c( 0, 0, 0, 0, 1/2-eR, 0, 0, 1/2+eR, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['KRAA','TTWA', ] <- c( 0, 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, (1/2+eR)/2, (1/2+eR)/2,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['KRAA','TTAA', ] <- c( 0, 0, 0, 0, 0, 1/2-eR, 0, 0, 1/2+eR,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['KRAA','TKWW', ] <- c( 0, 0, 0, 0, (1/2-eR)/2, 0, 0, (1/2+eR)/2, 0,
                                 0, (1/2-eR)^2, 0, 0, (1/2-eR)*(1/2+2*eR), 0, 0, eR*(1/2+eR), 0)
  tMatrix['KRAA','TKWA', ] <- c( 0, 0, 0, 0, (1/2-eR)/4, (1/2-eR)/4, 0, (1/2+eR)/4, (1/2+eR)/4,
                                 0, (1/2-eR)^2/2, (1/2-eR)^2/2, 0, (1/2-eR)*(1/2+2*eR)/2, (1/2-eR)*(1/2+2*eR)/2, 0, eR*(1/2+eR)/2, eR*(1/2+eR)/2)
  tMatrix['KRAA','TKAA', ] <- c( 0, 0, 0, 0, 0, (1/2-eR)/2, 0, 0, (1/2+eR)/2,
                                 0, 0, (1/2-eR)^2, 0, 0, (1/2-eR)*(1/2+2*eR), 0, 0, eR*(1/2+eR))
  tMatrix['KRAA','TRWW', ] <- c( 0, 0, 0, 0, (1/2-eR)/2, 0, 0, (1/2+eR)/2, 0,
                                 0, 0, 0, 0, (1/2-eR)/2, 0, 0, (1/2+eR)/2, 0)
  tMatrix['KRAA','TRWA', ] <- c( 0, 0, 0, 0, (1/2-eR)/4, (1/2-eR)/4, 0, (1/2+eR)/4, (1/2+eR)/4,
                                 0, 0, 0, 0, (1/2-eR)/4, (1/2-eR)/4, 0, (1/2+eR)/4, (1/2+eR)/4)
  tMatrix['KRAA','TRAA', ] <- c( 0, 0, 0, 0, 0, (1/2-eR)/2, 0, 0, (1/2+eR)/2,
                                 0, 0, 0, 0, 0, (1/2-eR)/2, 0, 0, (1/2+eR)/2)
  tMatrix['KRAA','KKWW', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 2*(1/2-eR)^2, 0, 0, 2*(1/2-eR)*(1/2+2*eR), 0, 0, 2*eR*(1/2+eR), 0)
  tMatrix['KRAA','KKWA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, (1/2-eR)^2, (1/2-eR)^2, 0, (1/2-eR)*(1/2+2*eR), (1/2-eR)*(1/2+2*eR), 0, eR*(1/2+eR), eR*(1/2+eR))
  tMatrix['KRAA','KKAA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 2*(1/2-eR)^2, 0, 0, 2*(1/2-eR)*(1/2+2*eR), 0, 0, 2*eR*(1/2+eR))
  tMatrix['KRAA','KRWW', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, (1/2-eR)^2, 0, 0, 2*(1/2-eR)*(1/2+eR), 0, 0, (1/2+eR)^2, 0)
  tMatrix['KRAA','KRWA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, (1/2-eR)^2/2, (1/2-eR)^2/2, 0, (1/2-eR)*(1/2+eR), (1/2-eR)*(1/2+eR), 0, (1/2+eR)^2/2, (1/2+eR)^2/2)
  tMatrix['KRAA','KRAA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, (1/2-eR)^2, 0, 0, 2*(1/2-eR)*(1/2+eR), 0, 0, (1/2+eR)^2)

  tMatrix['RRWW','TTWW', ] <- c( 0, 0, 0, 0, 0, 0, 1, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['RRWW','TTWA', ] <- c( 0, 0, 0, 0, 0, 0, 1/2, 1/2, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['RRWW','TTAA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 1, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['RRWW','TKWW', ] <- c( 0, 0, 0, 0, 0, 0, 1/2, 0, 0,
                                 0, 0, 0, 1/2-eR, 0, 0, eR, 0, 0)
  tMatrix['RRWW','TKWA', ] <- c( 0, 0, 0, 0, 0, 0, 1/4, 1/4, 0,
                                 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, eR/2, eR/2, 0)
  tMatrix['RRWW','TKAA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 1/2, 0,
                                 0, 0, 0, 0, 1/2-eR, 0, 0, eR, 0)
  tMatrix['RRWW','TRWW', ] <- c( 0, 0, 0, 0, 0, 0, 1/2, 0, 0,
                                 0, 0, 0, 0, 0, 0, 1/2, 0, 0)
  tMatrix['RRWW','TRWA', ] <- c( 0, 0, 0, 0, 0, 0, 1/4, 1/4, 0,
                                 0, 0, 0, 0, 0, 0, 1/4, 1/4, 0)
  tMatrix['RRWW','TRAA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 1/2, 0,
                                 0, 0, 0, 0, 0, 0, 0, 1/2, 0)
  tMatrix['RRWW','KKWW', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 1-2*eR, 0, 0, 2*eR, 0, 0)
  tMatrix['RRWW','KKWA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, (1/2-eR), (1/2-eR), 0, eR, eR, 0)
  tMatrix['RRWW','KKAA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 1-2*eR, 0, 0, 2*eR, 0)
  tMatrix['RRWW','KRWW', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 1/2-eR, 0, 0, 1/2+eR, 0, 0)
  tMatrix['RRWW','KRWA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, (1/2+eR)/2, (1/2+eR)/2, 0)
  tMatrix['RRWW','KRAA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 1/2-eR, 0, 0, 1/2+eR, 0)
  tMatrix['RRWW','RRWW', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 1, 0, 0)

  tMatrix['RRWA','TTWW', ] <- c( 0, 0, 0, 0, 0, 0, 1/2, 1/2, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['RRWA','TTWA', ] <- c( 0, 0, 0, 0, 0, 0, 1/4, 1/2, 1/4,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['RRWA','TTAA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 1/2, 1/2,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['RRWA','TKWW', ] <- c( 0, 0, 0, 0, 0, 0, 1/4, 1/4, 0,
                                 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, eR/2, eR/2, 0)
  tMatrix['RRWA','TKWA', ] <- c( 0, 0, 0, 0, 0, 0, 1/8, 1/4, 1/8,
                                 0, 0, 0, (1/2-eR)/4, (1/2-eR)/2, (1/2-eR)/4, eR/4, eR/2, eR/4)
  tMatrix['RRWA','TKAA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 1/4, 1/4,
                                 0, 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, eR/2, eR/2)
  tMatrix['RRWA','TRWW', ] <- c( 0, 0, 0, 0, 0, 0, 1/4, 1/4, 0,
                                 0, 0, 0, 0, 0, 0, 1/4, 1/4, 0)
  tMatrix['RRWA','TRWA', ] <- c( 0, 0, 0, 0, 0, 0, 1/8, 1/4, 1/8,
                                 0, 0, 0, 0, 0, 0, 1/8, 1/4, 1/8)
  tMatrix['RRWA','TRAA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 1/4, 1/4,
                                 0, 0, 0, 0, 0, 0, 0, 1/4, 1/4)
  tMatrix['RRWA','KKWW', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 1/2-eR, 1/2-eR, 0, eR, eR, 0)
  tMatrix['RRWA','KKWA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, (1/2-eR)/2, 1/2-eR, (1/2-eR)/2, eR/2, eR, eR/2)
  tMatrix['RRWA','KKAA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 1/2-eR, 1/2-eR, 0, eR, eR)
  tMatrix['RRWA','KRWW', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, (1/2+eR)/2, (1/2+eR)/2 ,0)
  tMatrix['RRWA','KRWA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, (1/2-eR)/4, (1/2-eR)/2, (1/2-eR)/4, (1/2+eR)/4, (1/2+eR)/2, (1/2+eR)/4)
  tMatrix['RRWA','KRAA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, (1/2+eR)/2, (1/2+eR)/2)
  tMatrix['RRWA','RRWW', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 1/2, 1/2, 0)
  tMatrix['RRWA','RRWA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 1/4, 1/2, 1/4)

  tMatrix['RRAA','TTWW', ] <- c( 0, 0, 0, 0, 0, 0, 0, 1, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['RRAA','TTWA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 1/2, 1/2,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['RRAA','TTAA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 1,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0)
  tMatrix['RRAA','TKWW', ] <- c( 0, 0, 0, 0, 0, 0, 0, 1/2, 0,
                                 0, 0, 0, 0, 1/2-eR, 0, 0, eR, 0)
  tMatrix['RRAA','TKWA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 1/4, 1/4,
                                 0, 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, eR/2, eR/2)
  tMatrix['RRAA','TKAA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 1/2,
                                 0, 0, 0, 0, 0, 1/2-eR, 0, 0, eR)
  tMatrix['RRAA','TRWW', ] <- c( 0, 0, 0, 0, 0, 0, 0, 1/2, 0,
                                 0, 0, 0, 0, 0, 0, 0, 1/2, 0)
  tMatrix['RRAA','TRWA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 1/4, 1/4,
                                 0, 0, 0, 0, 0, 0, 0, 1/4, 1/4)
  tMatrix['RRAA','TRAA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 1/2,
                                 0, 0, 0, 0, 0, 0, 0, 0, 1/2)
  tMatrix['RRAA','KKWW', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 1-2*eR, 0, 0, 2*eR, 0)
  tMatrix['RRAA','KKWA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, (1-2*eR)/2, (1-2*eR)/2, 0, eR, eR)
  tMatrix['RRAA','KKAA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 1-2*eR, 0, 0, 2*eR)
  tMatrix['RRAA','KRWW', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 1/2-eR, 0, 0, 1/2+eR, 0)
  tMatrix['RRAA','KRWA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, (1/2-eR)/2, (1/2-eR)/2, 0, (1/2+eR)/2, (1/2+eR)/2)
  tMatrix['RRAA','KRAA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 1/2-eR, 0, 0, 1/2+eR)
  tMatrix['RRAA','RRWW', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 1, 0)
  tMatrix['RRAA','RRWA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 1/2, 1/2)
  tMatrix['RRAA','RRAA', ] <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 1)

  ## set the other half of the matrix
  # Boolean matrix for subsetting, used several times
  boolMat <- upper.tri(x = tMatrix[ , ,1], diag = FALSE)
  # loop over depth, set upper triangle
  for(z in 1:size){tMatrix[ , ,z][boolMat] <- t(tMatrix[ , ,z])[boolMat]}

  tMatrix[tMatrix < .Machine$double.eps] <- 0 #protection from underflow errors


  ## initialize basic viability mask.
  viabilityMask <- array(data = 1, dim = c(size,size,size), dimnames = list(gtype, gtype, gtype))

  ## fill in offspring genotype specific death. Toxin and Antidote both have
  ##  predicted efficacies. The unions below represent the offspring that live.
  viabilityMask[ , ,c('TKWW', 'KRWW')] <- 1-Keff
  viabilityMask[ , ,c('TKWA', 'KRWA')] <- 1-Keff*(1-Aeff)
  viabilityMask[ , ,c('TKAA', 'KRAA')] <- 1-Keff*(1-Aeff)^2
  viabilityMask[ , ,'KKWW'] <- 2*(1-Keff)-(1-Keff)^2
  viabilityMask[ , ,'KKWA'] <- 1-Keff^2*(1-Aeff)
  viabilityMask[ , ,'KKAA'] <- 1-Keff^2*(1-Aeff)^2


  ## genotype-specific modifiers
  modifiers = cubeModifiers(gtype, eta = eta, phi = phi, omega = omega, xiF = xiF, xiM = xiM, s = s)


  ## put everything into a labeled list to return
  return(list(
    ih = tMatrix,
    tau = viabilityMask,
    genotypesID = gtype,
    genotypesN = size,
    wildType = "TTWW",
    eta = modifiers$eta,
    phi = modifiers$phi,
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = "KKAA"
  ))
}
