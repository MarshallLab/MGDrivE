###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   Split Drive with 2 Resistance Alleles and Sex-Specific homing
#   Héctor Sanchez, Jared Bennett, Sean Wu, John Marshall
#   August 2018
#   jared_bennett@berkeley.edu
#   January 2019
#    Update to reflect cutting, homing, resistance generation rates
#   June 2020
#    Updated with copy-number dependent homing
#    Updateed with copy-number dependent maternal deposition
#
###############################################################################

#' Inheritance Cube: Split CRISPR drive with 2 Resistance Alleles and male/female specific homing
#'
#' This is a sex-specific version of a split CRISPR drive. At one locus is the Cas9, inherited
#' in a Mendelian fashion. At a second, unlinked, locus are the gRNAs. When the two loci occur
#' together, the gRNAs drive, with potential damaged alleles, but the Cas9 remains
#' Mendelian. It is assumed that this is an autosomal drive.
#' This drive corresponds to the [confinable gene drive system](https://elifesciences.org/articles/51701)
#' developed by the Akbari lab.
#'
#' @param cM Cutting efficiency in males, one Cas9 allele
#' @param chM Homing efficiency in males, one Cas9 allele
#' @param crM Resistance efficiency in males, one Cas9 allele
#'
#' @param ccM Cutting efficiency in males, two Cas9 alleles
#' @param cchM Homing efficiency in males, two Cas9 alleles
#' @param ccrM Resistance efficiency in males, two Cas9 alleles
#'
#' @param cF Cutting efficiency in females, one Cas9 allele
#' @param chF Homing efficiency in females, one Cas9 allele
#' @param crF Resistance efficiency in females, one Cas9 allele
#'
#' @param ccF Cutting efficiency in females, two Cas9 alleles
#' @param cchF Homing efficiency in females, two Cas9 alleles
#' @param ccrF Resistance efficiency in females, two Cas9 alleles
#'
#' @param dW Maternal deposition cutting, one Cas9 allele
#' @param dhW Maternal deposition homing, one Cas9 allele
#' @param drW Maternal deposition resistance, one Cas9 allele
#'
#' @param ddW Maternal deposition cutting, two Cas9 alleles
#' @param ddhW Maternal deposition homing, two Cas9 alleles
#' @param ddrW Maternal deposition resistance, two Cas9 alleles
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
cubeSplitDrive <- function(cM = 1.0, chM = 0, crM = 0, ccM = cM, cchM = chM, ccrM = crM,
                           cF = 1.0, chF = 0, crF = 0, ccF = cF, cchF = chF, ccrF = crF,
                           dW = 0, dhW = 0, drW = 0, ddW = dW, ddhW = dhW, ddrW = drW,
                           eta = NULL, phi = NULL, omega = NULL, xiF = NULL, xiM = NULL, s = NULL){

  # # for testing
  # testVec <- runif(n = 18)
  #
  # cF <- testVec[1]; chF <- testVec[2]; crF <- testVec[3];
  # ccF <- testVec[4]; cchF <- testVec[5]; ccrF <- testVec[6];
  #
  # cM <- testVec[7]; chM <- testVec[8]; crM <- testVec[9];
  # ccM <- testVec[10]; cchM <- testVec[11]; ccrM <- testVec[12];
  #
  # dW <- testVec[13]; dhW <- testVec[14]; drW <- testVec[15];
  # ddW <- testVec[16]; ddhW <- testVec[17]; ddrW <- testVec[18];
  #
  # # this would run at the bottom, to check that all cells sum to 1
  # #  size is defined below
  # for(mID in 1:size){
  #   for(fID in 1:size){
  #     if((abs(sum(tMatrix[fID,mID,])-1)) > 1e-6){
  #       print(paste0("Fail: mID= ",gtype[mID],", fID= ",gtype[fID]))
  #     }
  #   }
  # }


  ## safety checks
  params <- c(cM, chM, crM, ccM, cchM, ccrM, cF, chF, crF, ccF, cchF, ccrF)
  if(any(params>1) || any(params<0)){
    stop("Parameters are rates, they must be between 0 and 1.")
  }


  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  gtype <- c('WWWW', 'WWWH', 'WWWR', 'WWWB', 'WWHH', 'WWHR', 'WWHB', 'WWRR', 'WWRB',
             'WWBB', 'WCWW', 'WCWH', 'WCWR', 'WCWB', 'WCHH', 'WCHR', 'WCHB', 'WCRR',
             'WCRB', 'WCBB', 'CCWW', 'CCWH', 'CCWR', 'CCWB', 'CCHH', 'CCHR', 'CCHB',
             'CCRR', 'CCRB', 'CCBB')
  size <- length(gtype) #because I use it several times
  tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix


  ## fill tMatrix with probabilities
  tMatrix['WWWW', 'WWWW', 'WWWW'] <- 1

  tMatrix['WWWH', 'WWWW', c('WWWW','WWWH')] <- c(0.5,0.5)
  tMatrix['WWWH', 'WWWH', c('WWWW','WWWH','WWHH')] <- c(0.25,0.5,0.25)

  tMatrix['WWWR', 'WWWW', c('WWWW','WWWR')] <- c(0.5,0.5)
  tMatrix['WWWR', 'WWWH', c('WWWW','WWWH','WWWR','WWHR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WWWR', 'WWWR', c('WWWW','WWWR','WWRR')] <- c(0.25,0.5,0.25)

  tMatrix['WWWB', 'WWWW', c('WWWW','WWWB')] <- c(0.5,0.5)
  tMatrix['WWWB', 'WWWH', c('WWWW','WWWH','WWWB','WWHB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WWWB', 'WWWR', c('WWWW','WWWR','WWWB','WWRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WWWB', 'WWWB', c('WWWW','WWWB','WWBB')] <- c(0.25,0.5,0.25)

  tMatrix['WWHH', 'WWWW', 'WWWH'] <- 1
  tMatrix['WWHH', 'WWWH', c('WWWH','WWHH')] <- c(0.5,0.5)
  tMatrix['WWHH', 'WWWR', c('WWWH','WWHR')] <- c(0.5,0.5)
  tMatrix['WWHH', 'WWWB', c('WWWH','WWHB')] <- c(0.5,0.5)
  tMatrix['WWHH', 'WWHH', 'WWHH'] <- 1

  tMatrix['WWHR', 'WWWW', c('WWWH','WWWR')] <- c(0.5,0.5)
  tMatrix['WWHR', 'WWWH', c('WWWH','WWWR','WWHH','WWHR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WWHR', 'WWWR', c('WWWH','WWWR','WWHR','WWRR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WWHR', 'WWWB', c('WWWH','WWWR','WWHB','WWRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WWHR', 'WWHH', c('WWHH','WWHR')] <- c(0.5,0.5)
  tMatrix['WWHR', 'WWHR', c('WWHH','WWHR','WWRR')] <- c(0.25,0.5,0.25)

  tMatrix['WWHB', 'WWWW', c('WWWH','WWWB')] <- c(0.5,0.5)
  tMatrix['WWHB', 'WWWH', c('WWWH','WWWB','WWHH','WWHB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WWHB', 'WWWR', c('WWWH','WWWB','WWHR','WWRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WWHB', 'WWWB', c('WWWH','WWWB','WWHB','WWBB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WWHB', 'WWHH', c('WWHH','WWHB')] <- c(0.5,0.5)
  tMatrix['WWHB', 'WWHR', c('WWHH','WWHR','WWHB','WWRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WWHB', 'WWHB', c('WWHH','WWHB','WWBB')] <- c(0.25,0.5,0.25)

  tMatrix['WWRR', 'WWWW', 'WWWR'] <- 1
  tMatrix['WWRR', 'WWWH', c('WWWR','WWHR')] <- c(0.5,0.5)
  tMatrix['WWRR', 'WWWR', c('WWWR','WWRR')] <- c(0.5,0.5)
  tMatrix['WWRR', 'WWWB', c('WWWR','WWRB')] <- c(0.5,0.5)
  tMatrix['WWRR', 'WWHH', 'WWHR'] <- 1
  tMatrix['WWRR', 'WWHR', c('WWHR','WWRR')] <- c(0.5,0.5)
  tMatrix['WWRR', 'WWHB', c('WWHR','WWRB')] <- c(0.5,0.5)
  tMatrix['WWRR', 'WWRR', 'WWRR'] <- 1

  tMatrix['WWRB', 'WWWW', c('WWWR','WWWB')] <- c(0.5,0.5)
  tMatrix['WWRB', 'WWWH', c('WWWR','WWWB','WWHR','WWHB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WWRB', 'WWWR', c('WWWR','WWWB','WWRR','WWRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WWRB', 'WWWB', c('WWWR','WWWB','WWRB','WWBB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WWRB', 'WWHH', c('WWHR','WWHB')] <- c(0.5,0.5)
  tMatrix['WWRB', 'WWHR', c('WWHR','WWHB','WWRR','WWRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WWRB', 'WWHB', c('WWHR','WWHB','WWRB','WWBB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WWRB', 'WWRR', c('WWRR','WWRB')] <- c(0.5,0.5)
  tMatrix['WWRB', 'WWRB', c('WWRR','WWRB','WWBB')] <- c(0.25,0.5,0.25)

  tMatrix['WWBB', 'WWWW', 'WWWB'] <- 1
  tMatrix['WWBB', 'WWWH', c('WWWB','WWHB')] <- c(0.5,0.5)
  tMatrix['WWBB', 'WWWR', c('WWWB','WWRB')] <- c(0.5,0.5)
  tMatrix['WWBB', 'WWWB', c('WWWB','WWBB')] <- c(0.5,0.5)
  tMatrix['WWBB', 'WWHH', 'WWHB'] <- 1
  tMatrix['WWBB', 'WWHR', c('WWHB','WWRB')] <- c(0.5,0.5)
  tMatrix['WWBB', 'WWHB', c('WWHB','WWBB')] <- c(0.5,0.5)
  tMatrix['WWBB', 'WWRR', 'WWRB'] <- 1
  tMatrix['WWBB', 'WWRB', c('WWRB','WWBB')] <- c(0.5,0.5)
  tMatrix['WWBB', 'WWBB', 'WWBB'] <- 1

  tMatrix['WCWW', 'WWWW', c('WWWW','WCWW')] <- c(0.5,0.5)
  tMatrix['WCWW', 'WWWH', c('WWWW','WWWH','WCWW','WCWH')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCWW', 'WWWR', c('WWWW','WWWR','WCWW','WCWR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCWW', 'WWWB', c('WWWW','WWWB','WCWW','WCWB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCWW', 'WWHH', c('WWWH','WCWH')] <- c(0.5,0.5)
  tMatrix['WCWW', 'WWHR', c('WWWH','WWWR','WCWH','WCWR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCWW', 'WWHB', c('WWWH','WWWB','WCWH','WCWB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCWW', 'WWRR', c('WWWR','WCWR')] <- c(0.5,0.5)
  tMatrix['WCWW', 'WWRB', c('WWWR','WWWB','WCWR','WCWB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCWW', 'WWBB', c('WWWB','WCWB')] <- c(0.5,0.5)
  tMatrix['WCWW', 'WCWW', c('WWWW','WCWW','CCWW')] <- c(0.25,0.5,0.25)

  tMatrix['WCWR', 'WWWW', c('WWWW','WWWR','WCWW','WCWR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCWR', 'WWWH', c('WWWW','WWWH','WWWR','WWHR',
                            'WCWW','WCWH','WCWR','WCHR')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WCWR', 'WWWR', c('WWWW','WWWR','WWRR',
                            'WCWW','WCWR','WCRR')] <- c(0.125,0.25,0.125,0.125,0.25,0.125)
  tMatrix['WCWR', 'WWWB', c('WWWW','WWWR','WWWB','WWRB',
                            'WCWW','WCWR','WCWB','WCRB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WCWR', 'WWHH', c('WWWH','WWHR','WCWH','WCHR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCWR', 'WWHR', c('WWWH','WWWR','WWHR','WWRR',
                            'WCWH','WCWR','WCHR','WCRR')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WCWR', 'WWHB', c('WWWH','WWWB','WWHR','WWRB',
                            'WCWH','WCWB','WCHR','WCRB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WCWR', 'WWRR', c('WWWR','WWRR','WCWR','WCRR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCWR', 'WWRB', c('WWWR','WWWB','WWRR','WWRB',
                            'WCWR','WCWB','WCRR','WCRB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WCWR', 'WWBB', c('WWWB','WWRB','WCWB','WCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCWR', 'WCWW', c('WWWW','WWWR','WCWW',
                            'WCWR','CCWW','CCWR')] <- c(0.125,0.125,0.25,0.25,0.125,0.125)
  tMatrix['WCWR', 'WCWR', c('WWWW','WWWR','WWRR',
                            'WCWW','WCWR','WCRR',
                            'CCWW','CCWR','CCRR')] <- c(1/16,0.125,1/16,0.125,0.25,0.125,1/16,0.125,1/16)

  tMatrix['WCWB', 'WWWW', c('WWWW','WWWB','WCWW','WCWB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCWB', 'WWWH', c('WWWW','WWWH','WWWB','WWHB',
                            'WCWW','WCWH','WCWB','WCHB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WCWB', 'WWWR', c('WWWW','WWWR','WWWB','WWRB',
                            'WCWW','WCWR','WCWB','WCRB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WCWB', 'WWWB', c('WWWW','WWWB','WWBB',
                            'WCWW','WCWB','WCBB')] <- c(0.125,0.25,0.125,0.125,0.25,0.125)
  tMatrix['WCWB', 'WWHH', c('WWWH','WWHB','WCWH','WCHB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCWB', 'WWHR', c('WWWH','WWWR','WWHB','WWRB',
                            'WCWH','WCWR','WCHB','WCRB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WCWB', 'WWHB', c('WWWH','WWWB','WWHB','WWBB',
                            'WCWH','WCWB','WCHB','WCBB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WCWB', 'WWRR', c('WWWR','WWRB','WCWR','WCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCWB', 'WWRB', c('WWWR','WWWB','WWRB','WWBB',
                            'WCWR','WCWB','WCRB','WCBB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WCWB', 'WWBB', c('WWWB','WWBB','WCWB','WCBB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCWB', 'WCWW', c('WWWW','WWWB','WCWW',
                            'WCWB','CCWW','CCWB')] <- c(0.125,0.125,0.25,0.25,0.125,0.125)
  tMatrix['WCWB', 'WCWR', c('WWWW','WWWR','WWWB','WWRB',
                            'WCWW','WCWR','WCWB','WCRB',
                            'CCWW','CCWR','CCWB','CCRB')] <- c(1/16,1/16,1/16,1/16,
                                                               0.125,0.125,0.125,0.125,
                                                               1/16,1/16,1/16,1/16)
  tMatrix['WCWB', 'WCWB', c('WWWW','WWWB','WWBB',
                            'WCWW','WCWB','WCBB',
                            'CCWW','CCWB','CCBB')] <- c(1/16,0.125,1/16,0.125,0.25,0.125,1/16,0.125,1/16)

  tMatrix['WCHH', 'WWHH', c('WWHH','WCHH')] <- c(0.5,0.5)
  tMatrix['WCHH', 'WWHR', c('WWHH','WWHR','WCHH','WCHR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCHH', 'WWHB', c('WWHH','WWHB','WCHH','WCHB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCHH', 'WWRR', c('WWHR','WCHR')] <- c(0.5,0.5)
  tMatrix['WCHH', 'WWRB', c('WWHR','WWHB','WCHR','WCHB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCHH', 'WWBB', c('WWHB','WCHB')] <- c(0.5,0.5)
  tMatrix['WCHH', 'WCHH', c('WWHH','WCHH','CCHH')] <- c(0.25,0.5,0.25)

  tMatrix['WCHR', 'WWHH', c('WWHH','WWHR','WCHH','WCHR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCHR', 'WWHR', c('WWHH','WWHR','WWRR',
                            'WCHH','WCHR','WCRR')] <- c(0.125,0.25,0.125,0.125,0.25,0.125)
  tMatrix['WCHR', 'WWHB', c('WWHH','WWHR','WWHB','WWRB',
                            'WCHH','WCHR','WCHB','WCRB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WCHR', 'WWRR', c('WWHR','WWRR','WCHR','WCRR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCHR', 'WWRB', c('WWHR','WWHB','WWRR','WWRB',
                            'WCHR','WCHB','WCRR','WCRB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WCHR', 'WWBB', c('WWHB','WWRB','WCHB','WCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCHR', 'WCHH', c('WWHH','WWHR','WCHH',
                            'WCHR','CCHH','CCHR')] <- c(0.125,0.25,0.125,0.125,0.25,0.125)
  tMatrix['WCHR', 'WCHR', c('WWHH','WWHR','WWRR',
                            'WCHH','WCHR','WCRR',
                            'CCHH','CCHR','CCRR')] <- c(1/16,0.125,1/16,0.125,0.25,0.125,1/16,0.125,1/16)

  tMatrix['WCHB', 'WWHH', c('WWHH','WWHB','WCHH','WCHB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCHB', 'WWHR', c('WWHH','WWHR','WWHB','WWRB',
                            'WCHH','WCHR','WCHB','WCRB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WCHB', 'WWHB', c('WWHH','WWHB','WWBB',
                            'WCHH','WCHB','WCBB')] <- c(0.125,0.25,0.125,0.125,0.25,0.125)
  tMatrix['WCHB', 'WWRR', c('WWHR','WWRB','WCHR','WCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCHB', 'WWRB', c('WWHR','WWHB','WWRB','WWBB',
                            'WCHR','WCHB','WCRB','WCBB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WCHB', 'WWBB', c('WWHB','WWBB','WCHB','WCBB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCHB', 'WCHH', c('WWHH','WWHB','WCHH',
                            'WCHB','CCHH','CCHB')] <- c(0.125,0.125,0.25,0.25,0.125,0.125)
  tMatrix['WCHB', 'WCHR', c('WWHH','WWHR','WWHB','WWRB',
                            'WCHH','WCHR','WCHB','WCRB',
                            'CCHH','CCHR','CCHB','CCRB')] <- c(1/16,1/16,1/16,1/16,
                                                               0.125,0.125,0.125,0.125,
                                                               1/16,1/16,1/16,1/16)
  tMatrix['WCHB', 'WCHB', c('WWHH','WWHB','WWBB',
                            'WCHH','WCHB','WCBB',
                            'CCHH','CCHB','CCBB')] <- c(1/16,0.125,1/16,0.125,0.25,0.125,1/16,0.125,1/16)

  tMatrix['WCRR', 'WWWW', c('WWWR','WCWR')] <- c(0.5,0.5)
  tMatrix['WCRR', 'WWWH', c('WWWR','WWHR','WCWR','WCHR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCRR', 'WWWR', c('WWWR','WWRR','WCWR','WCRR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCRR', 'WWWB', c('WWWR','WWRB','WCWR','WCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCRR', 'WWHH', c('WWHR','WCHR')] <- c(0.5,0.5)
  tMatrix['WCRR', 'WWHR', c('WWHR','WWRR','WCHR','WCRR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCRR', 'WWHB', c('WWHR','WWRB','WCHR','WCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCRR', 'WWRR', c('WWRR','WCRR')] <- c(0.5,0.5)
  tMatrix['WCRR', 'WWRB', c('WWRR','WWRB','WCRR','WCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCRR', 'WWBB', c('WWRB','WCRB')] <- c(0.5,0.5)
  tMatrix['WCRR', 'WCWW', c('WWWR','WCWR','CCWR')] <- c(0.25,0.5,0.25)
  tMatrix['WCRR', 'WCWR', c('WWWR','WWRR','WCWR',
                            'WCRR','CCWR','CCRR')] <- c(0.125,0.125,0.25,0.25,0.125,0.125)
  tMatrix['WCRR', 'WCWB', c('WWWR','WWRB','WCWR',
                            'WCRB','CCWR','CCRB')] <- c(0.125,0.125,0.25,0.25,0.125,0.125)
  tMatrix['WCRR', 'WCHH', c('WWHR','WCHR','CCHR')] <- c(0.25,0.5,0.25)
  tMatrix['WCRR', 'WCHR', c('WWHR','WWRR','WCHR',
                            'WCRR','CCHR','CCRR')] <- c(0.125,0.125,0.25,0.25,0.125,0.125)
  tMatrix['WCRR', 'WCHB', c('WWHR','WWRB','WCHR',
                            'WCRB','CCHR','CCRB')] <- c(0.125,0.125,0.25,0.25,0.125,0.125)
  tMatrix['WCRR', 'WCRR', c('WWRR','WCRR','CCRR')] <- c(0.25,0.5,0.25)

  tMatrix['WCRB', 'WWWW', c('WWWR','WWWB','WCWR','WCWB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCRB', 'WWWH', c('WWWR','WWWB','WWHR','WWHB',
                            'WCWR','WCWB','WCHR','WCHB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WCRB', 'WWWR', c('WWWR','WWWB','WWRR','WWRB',
                            'WCWR','WCWB','WCRR','WCRB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WCRB', 'WWWB', c('WWWR','WWWB','WWRB','WWBB',
                            'WCWR','WCWB','WCRB','WCBB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WCRB', 'WWHH', c('WWHR','WWHB','WCHR','WCHB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCRB', 'WWHR', c('WWHR','WWHB','WWRR','WWRB',
                            'WCHR','WCHB','WCRR','WCRB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WCRB', 'WWHB', c('WWHR','WWHB','WWRB','WWBB',
                           'WCHR','WCHB','WCRB','WCBB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WCRB', 'WWRR', c('WWRR','WWRB','WCRR','WCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCRB', 'WWRB', c('WWRR','WWRB','WWBB',
                            'WCRR','WCRB','WCBB')] <- c(0.125,0.25,0.125,0.125,0.25,0.125)
  tMatrix['WCRB', 'WWBB', c('WWRB','WWBB','WCRB','WCBB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCRB', 'WCWW', c('WWWR','WWWB','WCWR',
                            'WCWB','CCWR','CCWB')] <- c(0.125,0.125,0.25,0.25,0.125,0.125)
  tMatrix['WCRB', 'WCWR', c('WWWR','WWWB','WWRR','WWRB',
                            'WCWR','WCWB','WCRR','WCRB',
                            'CCWR','CCWB','CCRR','CCRB')] <- c(1/16,1/16,1/16,1/16,
                                                               0.125,0.125,0.125,0.125,
                                                               1/16,1/16,1/16,1/16)
  tMatrix['WCRB', 'WCWB', c('WWWR','WWWB','WWRB','WWBB',
                            'WCWR','WCWB','WCRB','WCBB',
                            'CCWR','CCWB','CCRB','CCBB')] <- c(1/16,1/16,1/16,1/16,
                                                               0.125,0.125,0.125,0.125,
                                                               1/16,1/16,1/16,1/16)
  tMatrix['WCRB', 'WCHH', c('WWHR','WWHB','WCHR',
                            'WCHB','CCHR','CCHB')] <- c(0.125,0.125,0.25,0.25,0.125,0.125)
  tMatrix['WCRB', 'WCHR', c('WWHR','WWHB','WWRR','WWRB',
                            'WCHR','WCHB','WCRR','WCRB',
                            'CCHR','CCHB','CCRR','CCRB')] <- c(1/16,1/16,1/16,1/16,
                                                               0.125,0.125,0.125,0.125,
                                                               1/16,1/16,1/16,1/16)
  tMatrix['WCRB', 'WCHB', c('WWHR','WWHB','WWRB','WWBB',
                            'WCHR','WCHB','WCRB','WCBB',
                            'CCHR','CCHB','CCRB','CCBB')] <- c(1/16,1/16,1/16,1/16,
                                                               0.125,0.125,0.125,0.125,
                                                               1/16,1/16,1/16,1/16)
  tMatrix['WCRB', 'WCRR', c('WWRR','WWRB','WCRR',
                            'WCRB','CCRR','CCRB')] <- c(0.125,0.125,0.25,0.25,0.125,0.125)
  tMatrix['WCRB', 'WCRB', c('WWRR','WWRB','WWBB',
                            'WCRR','WCRB','WCBB',
                            'CCRR','CCRB','CCBB')] <- c(1/16,0.125,1/16,0.125,0.25,0.125,1/16,0.125,1/16)

  tMatrix['WCBB', 'WWWW', c('WWWB','WCWB')] <- c(0.5,0.5)
  tMatrix['WCBB', 'WWWH', c('WWWB','WWHB','WCWB','WCHB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCBB', 'WWWR', c('WWWB','WWRB','WCWB','WCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCBB', 'WWWB', c('WWWB','WWBB','WCWB','WCBB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCBB', 'WWHH', c('WWHB','WCHB')] <- c(0.5,0.5)
  tMatrix['WCBB', 'WWHR', c('WWHB','WWRB','WCHB','WCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCBB', 'WWHB', c('WWHB','WWBB','WCHB','WCBB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCBB', 'WWRR', c('WWRB','WCRB')] <- c(0.5,0.5)
  tMatrix['WCBB', 'WWRB', c('WWRB','WWBB','WCRB','WCBB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCBB', 'WWBB', c('WWBB','WCBB')] <- c(0.5,0.5)
  tMatrix['WCBB', 'WCWW', c('WWWB','WCWB','CCWB')] <- c(0.25,0.5,0.25)
  tMatrix['WCBB', 'WCWR', c('WWWB','WWRB','WCWB',
                            'WCRB','CCWB','CCRB')] <- c(0.125,0.125,0.25,0.25,0.125,0.125)
  tMatrix['WCBB', 'WCWB', c('WWWB','WWBB','WCWB',
                            'WCBB','CCWB','CCBB')] <- c(0.125,0.125,0.25,0.25,0.125,0.125)
  tMatrix['WCBB', 'WCHH', c('WWHB','WCHB','CCHB')] <- c(0.25,0.5,0.25)
  tMatrix['WCBB', 'WCHR', c('WWHB','WWRB','WCHB',
                            'WCRB','CCHB','CCRB')] <- c(0.125,0.125,0.25,0.25,0.125,0.125)
  tMatrix['WCBB', 'WCHB', c('WWHB','WWBB','WCHB',
                            'WCBB','CCHB','CCBB')] <- c(0.125,0.125,0.25,0.25,0.125,0.125)
  tMatrix['WCBB', 'WCRR', c('WWRB','WCRB','CCRB')] <- c(0.25,0.5,0.25)
  tMatrix['WCBB', 'WCRB', c('WWRB','WWBB','WCRB',
                            'WCBB','CCRB','CCBB')] <- c(0.125,0.125,0.25,0.25,0.125,0.125)
  tMatrix['WCBB', 'WCBB', c('WWBB','WCBB','CCBB')] <- c(0.25,0.5,0.25)

  tMatrix['CCWW', 'WWWW', 'WCWW'] <- 1
  tMatrix['CCWW', 'WWWH', c('WCWW','WCWH')] <- c(0.5,0.5)
  tMatrix['CCWW', 'WWWR', c('WCWW','WCWR')] <- c(0.5,0.5)
  tMatrix['CCWW', 'WWWB', c('WCWW','WCWB')] <- c(0.5,0.5)
  tMatrix['CCWW', 'WWHH', 'WCWH'] <- 1
  tMatrix['CCWW', 'WWHR', c('WCWH','WCWR')] <- c(0.5,0.5)
  tMatrix['CCWW', 'WWHB', c('WCWH','WCWB')] <- c(0.5,0.5)
  tMatrix['CCWW', 'WWRR', 'WCWR'] <- 1
  tMatrix['CCWW', 'WWRB', c('WCWR','WCWB')] <- c(0.5,0.5)
  tMatrix['CCWW', 'WWBB', 'WCWB'] <- 1
  tMatrix['CCWW', 'WCWW', c('WCWW','CCWW')] <- c(0.5,0.5)
  tMatrix['CCWW', 'WCWR', c('WCWW','WCWR','CCWW','CCWR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWW', 'WCWB', c('WCWW','WCWB','CCWW','CCWB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWW', 'WCRR', c('WCWR','CCWR')] <- c(0.5,0.5)
  tMatrix['CCWW', 'WCRB', c('WCWR','WCWB','CCWR','CCWB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWW', 'WCBB', c('WCWB','CCWB')] <- c(0.5,0.5)
  tMatrix['CCWW', 'CCWW', 'CCWW'] <- 1

  tMatrix['CCWR', 'WWWW', c('WCWW','WCWR')] <- c(0.5,0.5)
  tMatrix['CCWR', 'WWWH', c('WCWW','WCWH','WCWR','WCHR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWR', 'WWWR', c('WCWW','WCWR','WCRR')] <- c(0.25,0.5,0.25)
  tMatrix['CCWR', 'WWWB', c('WCWW','WCWR','WCWB','WCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWR', 'WWHH', c('WCWH','WCHR')] <- c(0.5,0.5)
  tMatrix['CCWR', 'WWHR', c('WCWH','WCWR','WCHR','WCRR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWR', 'WWHB', c('WCWH','WCWB','WCHR','WCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWR', 'WWRR', c('WCWR','WCRR')] <- c(0.5,0.5)
  tMatrix['CCWR', 'WWRB', c('WCWR','WCWB','WCRR','WCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWR', 'WWBB', c('WCWB','WCRB')] <- c(0.5,0.5)
  tMatrix['CCWR', 'WCWW', c('WCWW','WCWR','CCWW','CCWR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWR', 'WCWR', c('WCWW','WCWR','WCRR',
                            'CCWW','CCWR','CCRR')] <- c(0.125,0.25,0.125,0.125,0.25,0.125)
  tMatrix['CCWR', 'WCWB', c('WCWW','WCWR','WCWB','WCRB',
                            'CCWW','CCWR','CCWB','CCRB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['CCWR', 'WCRR', c('WCWR','WCRR','CCWR','CCRR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWR', 'WCRB', c('WCWR','WCWB','WCRR','WCRB',
                            'CCWR','CCWB','CCRR','CCRB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['CCWR', 'WCBB', c('WCWB','WCRB','CCWB','CCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWR', 'CCWW', c('CCWW','CCWR')] <- c(0.5,0.5)
  tMatrix['CCWR', 'CCWR', c('CCWW','CCWR','CCRR')] <- c(0.25,0.5,0.25)

  tMatrix['CCWB', 'WWWW', c('WCWW','WCWB')] <- c(0.5,0.5)
  tMatrix['CCWB', 'WWWH', c('WCWW','WCWH','WCWB','WCHB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWB', 'WWWR', c('WCWW','WCWR','WCWB','WCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWB', 'WWWB', c('WCWW','WCWB','WCBB')] <- c(0.25,0.5,0.25)
  tMatrix['CCWB', 'WWHH', c('WCWH','WCHB')] <- c(0.5,0.5)
  tMatrix['CCWB', 'WWHR', c('WCWH','WCWR','WCHB','WCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWB', 'WWHB', c('WCWH','WCWB','WCHB','WCBB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWB', 'WWRR', c('WCWR','WCRB')] <- c(0.5,0.5)
  tMatrix['CCWB', 'WWRB', c('WCWR','WCWB','WCRB','WCBB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWB', 'WWBB', c('WCWB','WCBB')] <- c(0.5,0.5)
  tMatrix['CCWB', 'WCWW', c('WCWW','WCWB','CCWW','CCWB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWB', 'WCWR', c('WCWW','WCWR','WCWB','WCRB',
                            'CCWW','CCWR','CCWB','CCRB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['CCWB', 'WCWB', c('WCWW','WCWB','WCBB',
                            'CCWW','CCWB','CCBB')] <- c(0.125,0.25,0.125,0.125,0.25,0.125)
  tMatrix['CCWB', 'WCRR', c('WCWR','WCRB','CCWR','CCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWB', 'WCRB', c('WCWR','WCWB','WCRB','WCBB',
                            'CCWR','CCWB','CCRB','CCBB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['CCWB', 'WCBB', c('WCWB','WCBB','CCWB','CCBB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWB', 'CCWW', c('CCWW','CCWB')] <- c(0.5,0.5)
  tMatrix['CCWB', 'CCWR', c('CCWW','CCWR','CCWB','CCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWB', 'CCWB', c('CCWW','CCWB','CCBB')] <- c(0.25,0.5,0.25)

  tMatrix['CCHH', 'WWHH', 'WCHH'] <- 1
  tMatrix['CCHH', 'WWHR', c('WCHH','WCHR')] <- c(0.5,0.5)
  tMatrix['CCHH', 'WWHB', c('WCHH','WCHB')] <- c(0.5,0.5)
  tMatrix['CCHH', 'WWRR', 'WCHR'] <- 1
  tMatrix['CCHH', 'WWRB', c('WCHR','WCHB')] <- c(0.5,0.5)
  tMatrix['CCHH', 'WWBB', 'WCHB'] <- 1
  tMatrix['CCHH', 'WCHH', c('WCHH','CCHH')] <- c(0.5,0.5)
  tMatrix['CCHH', 'WCHR', c('WCHH','WCHR','CCHH','CCHR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCHH', 'WCHB', c('WCHH','WCHB','CCHH','CCHB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCHH', 'WCRR', c('WCHR','CCHR')] <- c(0.5,0.5)
  tMatrix['CCHH', 'WCRB', c('WCHR','WCHB','CCHR','CCHB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCHH', 'WCBB', c('WCHB','CCHB')] <- c(0.5,0.5)
  tMatrix['CCHH', 'CCHH', 'CCHH'] <- 1

  tMatrix['CCHR', 'WWHH', c('WCHH','WCHR')] <- c(0.5,0.5)
  tMatrix['CCHR', 'WWHR', c('WCHH','WCHR','WCRR')] <- c(0.25,0.5,0.25)
  tMatrix['CCHR', 'WWHB', c('WCHH','WCHR','WCHB','WCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCHR', 'WWRR', c('WCHR','WCRR')] <- c(0.5,0.5)
  tMatrix['CCHR', 'WWRB', c('WCHR','WCHB','WCRR','WCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCHR', 'WWBB', c('WCHB','WCRB')] <- c(0.5,0.5)
  tMatrix['CCHR', 'WCHH', c('WCHH','WCHR','CCHH','CCHR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCHR', 'WCHR', c('WCHH','WCHR','WCRR',
                            'CCHH','CCHR','CCRR')] <- c(0.125,0.25,0.125,0.125,0.25,0.125)
  tMatrix['CCHR', 'WCHB', c('WCHH','WCHR','WCHB','WCRB',
                            'CCHH','CCHR','CCHB','CCRB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['CCHR', 'WCRR', c('WCHR','WCRR','CCHR','CCRR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCHR', 'WCRB', c('WCHR','WCHB','WCRR','WCRB',
                            'CCHR','CCHB','CCRR','CCRB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['CCHR', 'WCBB', c('WCHB','WCRB','CCHB','CCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCHR', 'CCHH', c('CCHH','CCHR')] <- c(0.5,0.5)
  tMatrix['CCHR', 'CCHR', c('CCHH','CCHR','CCRR')] <- c(0.25,0.5,0.25)

  tMatrix['CCHB', 'WWHH', c('WCHH','WCHB')]<- c(0.5,0.5)
  tMatrix['CCHB', 'WWHR', c('WCHH','WCHR','WCHB','WCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCHB', 'WWHB', c('WCHH','WCHB','WCBB')] <- c(0.25,0.5,0.25)
  tMatrix['CCHB', 'WWRR', c('WCHR','WCRB')] <- c(0.5,0.5)
  tMatrix['CCHB', 'WWRB', c('WCHR','WCHB','WCRB','WCBB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCHB', 'WWBB', c('WCHB','WCBB')] <- c(0.5,0.5)
  tMatrix['CCHB', 'WCHH', c('WCHH','WCHB','CCHH','CCHB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCHB', 'WCHR', c('WCHH','WCHR','WCHB','WCRB',
                            'CCHH','CCHR','CCHB','CCRB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['CCHB', 'WCHB', c('WCHH','WCHB','WCBB',
                            'CCHH','CCHB','CCBB')] <- c(0.125,0.25,0.125,0.125,0.25,0.125)
  tMatrix['CCHB', 'WCRR', c('WCHR','WCRB','CCHR','CCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCHB', 'WCRB', c('WCHR','WCHB','WCRB','WCBB',
                            'CCHR','CCHB','CCRB','CCBB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['CCHB', 'WCBB', c('WCHB','WCBB','CCHB','CCBB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCHB', 'CCHH', c('CCHH','CCHB')] <- c(0.5,0.5)
  tMatrix['CCHB', 'CCHR', c('CCHH','CCHR','CCHB','CCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCHB', 'CCHB', c('CCHH','CCHB','CCBB')] <- c(0.25,0.5,0.25)

  tMatrix['CCRR', 'WWWW', 'WCWR'] <- 1
  tMatrix['CCRR', 'WWWH', c('WCWR','WCHR')] <- c(0.5,0.5)
  tMatrix['CCRR', 'WWWR', c('WCWR','WCRR')] <- c(0.5,0.5)
  tMatrix['CCRR', 'WWWB', c('WCWR','WCRB')] <- c(0.5,0.5)
  tMatrix['CCRR', 'WWHH', 'WCHR'] <- 1
  tMatrix['CCRR', 'WWHR', c('WCHR','WCRR')] <- c(0.5,0.5)
  tMatrix['CCRR', 'WWHB', c('WCHR','WCRB')] <- c(0.5,0.5)
  tMatrix['CCRR', 'WWRR', 'WCRR'] <- 1
  tMatrix['CCRR', 'WWRB', c('WCRR','WCRB')] <- c(0.5,0.5)
  tMatrix['CCRR', 'WWBB', 'WCRB'] <- 1
  tMatrix['CCRR', 'WCWW', c('WCWR','CCWR')] <- c(0.5,0.5)
  tMatrix['CCRR', 'WCWR', c('WCWR','WCRR','CCWR','CCRR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCRR', 'WCWB', c('WCWR','WCRB','CCWR','CCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCRR', 'WCHH', c('WCHR','CCHR')] <- c(0.5,0.5)
  tMatrix['CCRR', 'WCHR', c('WCHR','WCRR','CCHR','CCRR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCRR', 'WCHB', c('WCHR','WCRB','CCHR','CCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCRR', 'WCRR', c('WCRR','CCRR')] <- c(0.5,0.5)
  tMatrix['CCRR', 'WCRB', c('WCRR','WCRB','CCRR','CCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCRR', 'WCBB', c('WCRB','CCRB')] <- c(0.5,0.5)
  tMatrix['CCRR', 'CCWW', 'CCWR'] <- 1
  tMatrix['CCRR', 'CCWR', c('CCWR','CCRR')] <- c(0.5,0.5)
  tMatrix['CCRR', 'CCWB', c('CCWR','CCRB')] <- c(0.5,0.5)
  tMatrix['CCRR', 'CCHH', 'CCHR'] <- 1
  tMatrix['CCRR', 'CCHR', c('CCHR','CCRR')] <- c(0.5,0.5)
  tMatrix['CCRR', 'CCHB', c('CCHR','CCRB')] <- c(0.5,0.5)
  tMatrix['CCRR', 'CCRR', 'CCRR'] <- 1

  tMatrix['CCRB', 'WWWW', c('WCWR','WCWB')] <- c(0.5,0.5)
  tMatrix['CCRB', 'WWWH', c('WCWR','WCWB','WCHR','WCHB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCRB', 'WWWR', c('WCWR','WCWB','WCRR','WCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCRB', 'WWWB', c('WCWR','WCWB','WCRB','WCBB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCRB', 'WWHH', c('WCHR','WCHB')] <- c(0.5,0.5)
  tMatrix['CCRB', 'WWHR', c('WCHR','WCHB','WCRR','WCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCRB', 'WWHB', c('WCHR','WCHB','WCRB','WCBB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCRB', 'WWRR', c('WCRR','WCRB')] <- c(0.5,0.5)
  tMatrix['CCRB', 'WWRB', c('WCRR','WCRB','WCBB')] <- c(0.25,0.5,0.25)
  tMatrix['CCRB', 'WWBB', c('WCRB','WCBB')] <- c(0.5,0.5)
  tMatrix['CCRB', 'WCWW', c('WCWR','WCWB','CCWR','CCWB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCRB', 'WCWR', c('WCWR','WCWB','WCRR','WCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCRB', 'WCWB', c('WCWR','WCWB','WCRB','WCBB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCRB', 'WCHH', c('WCHR','WCHB','CCHR','CCHB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCRB', 'WCHR', c('WCHR','WCHB','WCRR','WCRB',
                            'CCHR','CCHB','CCRR','CCRB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['CCRB', 'WCHB', c('WCHR','WCHB','WCRB','WCBB',
                            'CCHR','CCHB','CCRB','CCBB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['CCRB', 'WCRR', c('WCRR','WCRB','CCRR','CCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCRB', 'WCRB', c('WCRR','WCRB','WCBB',
                            'CCRR','CCRB','CCBB')] <- c(0.125,0.25,0.125,0.125,0.25,0.125)
  tMatrix['CCRB', 'WCBB', c('WCRB','WCBB','CCRB','CCBB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCRB', 'CCWW', c('CCWR','CCWB')] <- c(0.5,0.5)
  tMatrix['CCRB', 'CCWR', c('CCWR','CCWB','CCRR','CCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCRB', 'CCWB', c('CCWR','CCWB','CCRB','CCBB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCRB', 'CCHH', c('CCHR','CCHB')] <- c(0.5,0.5)
  tMatrix['CCRB', 'CCHR', c('CCHR','CCHB','CCRR','CCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCRB', 'CCHB', c('CCHR','CCHB','CCRB','CCBB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCRB', 'CCRR', c('CCRR','CCRB')] <- c(0.5,0.5)
  tMatrix['CCRB', 'CCRB', c('CCRR','CCRB','CCBB')] <- c(0.25,0.5,0.25)

  tMatrix['CCBB', 'WWWW', 'WCWB'] <- 1
  tMatrix['CCBB', 'WWWH', c('WCWB','WCHB')] <- c(0.5,0.5)
  tMatrix['CCBB', 'WWWR', c('WCWB','WCRB')] <- c(0.5,0.5)
  tMatrix['CCBB', 'WWWB', c('WCWB','WCBB')] <- c(0.5,0.5)
  tMatrix['CCBB', 'WWHH', 'WCHB'] <- 1
  tMatrix['CCBB', 'WWHR', c('WCWB','WCRB')] <- c(0.5,0.5)
  tMatrix['CCBB', 'WWHB', c('WCHB','WCBB')] <- c(0.5,0.5)
  tMatrix['CCBB', 'WWRR', 'WCRB'] <- 1
  tMatrix['CCBB', 'WWRB', c('WCRB','WCBB')] <- c(0.5,0.5)
  tMatrix['CCBB', 'WWBB', 'WCBB'] <- 1
  tMatrix['CCBB', 'WCWW', c('WCWB','CCWB')] <- c(0.5,0.5)
  tMatrix['CCBB', 'WCWR', c('WCWB','WCRB','CCWB','CCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCBB', 'WCWB', c('WCWB','WCBB','CCWB','CCBB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCBB', 'WCHH', 'WCHB'] <- 1
  tMatrix['CCBB', 'WCHR', c('WCHB','WCRB')] <- c(0.5,0.5)
  tMatrix['CCBB', 'WCHB', c('WCHB','WCBB')] <- c(0.5,0.5)
  tMatrix['CCBB', 'WCRR', c('WCRB','CCRB')] <- c(0.5,0.5)
  tMatrix['CCBB', 'WCRB', c('WCRB','WCBB','CCRB','CCBB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCBB', 'WCBB', c('WCBB','CCBB')] <- c(0.5,0.5)
  tMatrix['CCBB', 'CCWW', 'CCWB'] <- 1
  tMatrix['CCBB', 'CCWR', c('CCWB','CCRB')] <- c(0.5,0.5)
  tMatrix['CCBB', 'CCWB', c('CCWB','CCBB')] <- c(0.5,0.5)
  tMatrix['CCBB', 'CCHH', 'CCHB'] <- 1
  tMatrix['CCBB', 'CCHR', c('CCHB','CCRB')] <- c(0.5,0.5)
  tMatrix['CCBB', 'CCHB', c('CCHB','CCBB')] <- c(0.5,0.5)
  tMatrix['CCBB', 'CCRR', 'CCRB'] <- 1
  tMatrix['CCBB', 'CCRB', c('CCRB','CCBB')] <- c(0.5,0.5)
  tMatrix['CCBB', 'CCBB', 'CCBB'] <- 1

  ## set the other half of the matrix that is symmetric
  # Boolean matrix for subsetting, used several times
  boolMat <- upper.tri(x = tMatrix[ , ,1], diag = FALSE)
  # loop over depth, set upper triangle
  for(z in 1:size){tMatrix[ , ,z][boolMat] <- t(tMatrix[ , ,z])[boolMat]}


  ## fill asymmetric parts of tMatrix
  ## non-deposition stuff
  tMatrix['WWWW', 'WCHH', c('WWWH','WCWH')] <- c(0.5,0.5)
  tMatrix['WWWH', 'WCHH', c('WWWH','WWHH','WCWH','WCHH')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WWWR', 'WCHH', c('WWWH','WWHR','WCWH','WCHR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WWWB', 'WCHH', c('WWWH','WWHB','WCWH','WCHB')] <- c(0.25,0.25,0.25,0.25)

  tMatrix['WWWW', 'WCHR', c('WWWH','WWWR','WCWH','WCWR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WWWH', 'WCHR', c('WWWH','WWWR','WWHH','WWHR',
                            'WCWH','WCWR','WCHH','WCHR')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WWWR', 'WCHR', c('WWWH','WWWR','WWHR','WWRR',
                            'WCWH','WCWR','WCHR','WCRR')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WWWB', 'WCHR', c('WWWH','WWWR','WWHB','WWRB',
                            'WCWH','WCWR','WCHB','WCRB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)

  tMatrix['WWWW', 'WCHB', c('WWWH','WWWB','WCWH','WCWB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WWWH', 'WCHB', c('WWWH','WWWB','WWHH','WWHB',
                            'WCWH','WCWB','WCHH','WCHB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WWWR', 'WCHB', c('WWWH','WWWB','WWHR','WWRB',
                            'WCWH','WCWB','WCHR','WCRB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WWWB', 'WCHB', c('WWWH','WWWB','WWHB','WWBB',
                            'WCWH','WCWB','WCHB','WCBB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)


  tMatrix['WCWW', 'WCHH', c('WWWH','WCWH','CCWH')] <- c(0.25,0.5,0.25)
  tMatrix['WCWR', 'WCHH', c('WWWH','WWHR','WCWH',
                            'WCHR','CCWH','CCHR')] <- c(0.125,0.125,0.25,0.25,0.125,0.125)
  tMatrix['WCWB', 'WCHH', c('WWWH','WWHB','WCWH',
                            'WCHB','CCWH','CCHB')] <- c(0.125,0.125,0.25,0.25,0.125,0.125)


  tMatrix['WCWW', 'WCHR', c('WWWH','WWWR','WCWH',
                            'WCWR','CCWH','CCWR')] <- c(0.125,0.125,0.25,0.25,0.125,0.125)
  tMatrix['WCWR', 'WCHR', c('WWWH','WWWR','WWHR','WWRR',
                            'WCWH','WCWR','WCHR','WCRR',
                            'CCWH','CCWR','CCHR','CCRR')] <- c(1/16,1/16,1/16,1/16,
                                                               0.125,0.125,0.125,0.125,
                                                               1/16,1/16,1/16,1/16)
  tMatrix['WCWB', 'WCHR', c('WWWH','WWWR','WWHB','WWRB',
                            'WCWH','WCWR','WCHB','WCRB',
                            'CCWH','CCWR','CCHB','CCRB')] <- c(1/16,1/16,1/16,1/16,
                                                               0.125,0.125,0.125,0.125,
                                                               1/16,1/16,1/16,1/16)


  tMatrix['WCWW', 'WCHB', c('WWWH','WWWB','WCWH',
                            'WCWB','CCWH','CCWB')] <- c(0.125,0.125,0.25,0.25,0.125,0.125)
  tMatrix['WCWR', 'WCHB', c('WWWH','WWWB','WWHR','WWRB',
                            'WCWH','WCWB','WCHR','WCRB',
                            'CCWH','CCWB','CCHR','CCRB')] <- c(1/16,1/16,1/16,1/16,
                                                               0.125,0.125,0.125,0.125,
                                                               1/16,1/16,1/16,1/16)
  tMatrix['WCWB', 'WCHB', c('WWWH','WWWB','WWHB','WWBB',
                            'WCWH','WCWB','WCHB','WCBB',
                            'CCWH','CCWB','CCHB','CCBB')] <- c(1/16,1/16,1/16,1/16,
                                                               0.125,0.125,0.125,0.125,
                                                               1/16,1/16,1/16,1/16)


  tMatrix['CCWW', 'WCHH', c('WCWH','CCWH')] <- c(0.5,0.5)
  tMatrix['CCWW', 'WCHR', c('WCWH','WCWR','CCWH','CCWR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWW', 'WCHB', c('WCWH','WCWB','CCWH','CCWB')] <- c(0.25,0.25,0.25,0.25)

  tMatrix['CCWR', 'WCHH', c('WCWH','WCHR','CCWH','CCHR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWR', 'WCHR', c('WCWH','WCWR','WCHR','WCRR',
                            'CCWH','CCWR','CCHR','CCRR')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['CCWR', 'WCHB', c('WCWH','WCWB','WCHR','WCRB',
                            'CCWH','CCWB','CCHR','CCRB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)

  tMatrix['CCWB', 'WCHH', c('WCWH','WCHB','CCWH','CCHB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWB', 'WCHR', c('WCWH','WCWR','WCHB','WCRB',
                            'CCWH','CCWR','CCHB','CCRB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['CCWB', 'WCHB', c('WCWH','WCWB','WCHB','WCBB',
                            'CCWH','CCWB','CCHB','CCBB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)


  tMatrix['WWWW', 'CCHH', 'WCWH'] <- 1
  tMatrix['WWWH', 'CCHH', c('WCWH','WCHH')] <- c(0.5,0.5)
  tMatrix['WWWR', 'CCHH', c('WCWH','WCHR')] <- c(0.5,0.5)
  tMatrix['WWWB', 'CCHH', c('WCWH','WCHB')] <- c(0.5,0.5)

  tMatrix['WWWW', 'CCHR', c('WCWH','WCWR')] <- c(0.5,0.5)
  tMatrix['WWWH', 'CCHR', c('WCWH','WCWR','WCHH','WCHR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WWWR', 'CCHR', c('WCWH','WCWR','WCHR','WCRR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WWWB', 'CCHR', c('WCWH','WCWR','WCHB','WCRB')] <- c(0.25,0.25,0.25,0.25)

  tMatrix['WWWW', 'CCHB', c('WCWH','WCWB')] <- c(0.5,0.5)
  tMatrix['WWWH', 'CCHB', c('WCWH','WCWB','WCHH','WCHB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WWWR', 'CCHB', c('WCWH','WCWB','WCHR','WCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WWWB', 'CCHB', c('WCWH','WCWB','WCHB','WCBB')] <- c(0.25,0.25,0.25,0.25)


  tMatrix['WCWW', 'CCHH', c('WCWH','CCWH')] <- c(0.5,0.5)
  tMatrix['WCWR', 'CCHH', c('WCWH','WCHR','CCWH','CCHR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCWB', 'CCHH', c('WCWH','WCHB','CCWH','CCHB')] <- c(0.25,0.25,0.25,0.25)

  tMatrix['WCWW', 'CCHR', c('WCWH','WCWR','CCWH','CCWR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCWR', 'CCHR', c('WCWH','WCWR','WCHR','WCRR',
                            'CCWH','CCWR','CCHR','CCRR')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WCWB', 'CCHR', c('WCWH','WCWR','WCHB','WCRB',
                            'CCWH','CCWR','CCHB','CCRB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)

  tMatrix['WCWW', 'CCHB', c('WCWH','WCWB','CCWH','CCWB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['WCWR', 'CCHB', c('WCWH','WCWB','WCHR','WCRB',
                            'CCWH','CCWB','CCHR','CCRB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)
  tMatrix['WCWB', 'CCHB', c('WCWH','WCWB','WCHB','WCBB',
                            'CCWH','CCWB','CCHB','CCBB')] <- c(0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125)


  tMatrix['CCWW', 'CCHH', 'CCWH'] <- 1
  tMatrix['CCWR', 'CCHH', c('CCWH','CCHR')] <- c(0.5,0.5)
  tMatrix['CCWB', 'CCHH', c('CCWH','CCHB')] <- c(0.5,0.5)

  tMatrix['CCWW', 'CCHR', c('CCWH','CCWR')] <- c(0.5,0.5)
  tMatrix['CCWR', 'CCHR', c('CCWH','CCWR','CCHR','CCRR')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWB', 'CCHR', c('CCWH','CCWR','CCHB','CCRB')] <- c(0.25,0.25,0.25,0.25)

  tMatrix['CCWW', 'CCHB', c('CCWH','CCWB')] <- c(0.5,0.5)
  tMatrix['CCWR', 'CCHB', c('CCWH','CCWB','CCHR','CCRB')] <- c(0.25,0.25,0.25,0.25)
  tMatrix['CCWB', 'CCHB', c('CCWH','CCWB','CCHB','CCBB')] <- c(0.25,0.25,0.25,0.25)


  ## deposition stuff
  eTen <- rep.int(x = 0, times = 10)

  wrbTen <- c(0, 1-dW, 0, 0, 0, dW*drW, dW*(1-drW), 0, 0, 0)
  tMatrix['WCHH', 'WWWW', ] <- c(wrbTen, wrbTen, eTen)/2
  tMatrix['WCHH', 'WCWW', ] <- c(wrbTen, 2*wrbTen, wrbTen)/4
  tMatrix['WCHH', 'CCWW', ] <- c(eTen, wrbTen, wrbTen)/2

  wrbTen <- c(0, 1-dW, 0, 0, 1+dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 0, 0, 0)
  tMatrix['WCHH', 'WWWH', ] <- c(wrbTen, wrbTen, eTen)/4

  wrbTen <- c(0, 1-dW, 0, 0, 0, 1+dW*drW, dW*(1-drW), 0, 0, 0)
  tMatrix['WCHH', 'WWWR', ] <- c(wrbTen, wrbTen, eTen)/4
  tMatrix['WCHH', 'WCWR', ] <- c(wrbTen, 2*wrbTen, wrbTen)/8
  tMatrix['WCHH', 'CCWR', ] <- c(eTen, wrbTen, wrbTen)/4

  wrbTen <- c(0, 1-dW, 0, 0, 0, dW*drW, 1+dW*(1-drW), 0, 0, 0)
  tMatrix['WCHH', 'WWWB', ] <- c(wrbTen, wrbTen, eTen)/4
  tMatrix['WCHH', 'WCWB', ] <- c(wrbTen, 2*wrbTen, wrbTen)/8
  tMatrix['WCHH', 'CCWB', ] <- c(eTen, wrbTen, wrbTen)/4


  wrbTen <- c(0, 1-dW, 1-dW, 0, 0, dW*drW, dW*(1-drW), dW*drW, dW*(1-drW), 0)
  tMatrix['WCHR', 'WWWW', ] <- c(wrbTen, wrbTen, eTen)/4
  tMatrix['WCHR', 'WCWW', ] <- c(wrbTen, 2*wrbTen, wrbTen)/8
  tMatrix['WCHR', 'CCWW', ] <- c(eTen, wrbTen, wrbTen)/4

  wrbTen <- c(0, 1-dW, 1-dW, 0, 1+dW*dhW, 1+dW*dhW + dW*(1-dhW)*drW,
              dW*(1-dhW)*(1-drW), dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 0)
  tMatrix['WCHR', 'WWWH', ] <- c(wrbTen, wrbTen, eTen)/8

  wrbTen <- c(0, 1-dW, 1-dW, 0, 0, 1+dW*drW, dW*(1-drW), 1+dW*drW, dW*(1-drW), 0)
  tMatrix['WCHR', 'WWWR', ] <- c(wrbTen, wrbTen, eTen)/8
  tMatrix['WCHR', 'WCWR', ] <- c(wrbTen, 2*wrbTen, wrbTen)/16
  tMatrix['WCHR', 'CCWR', ] <- c(eTen, wrbTen, wrbTen)/8

  wrbTen <- c(0, 1-dW, 1-dW, 0, 0, dW*drW, 1+dW*(1-drW), dW*drW, 1+dW*(1-drW), 0)
  tMatrix['WCHR', 'WWWB', ] <- c(wrbTen, wrbTen, eTen)/8
  tMatrix['WCHR', 'WCWB', ] <- c(wrbTen, 2*wrbTen, wrbTen)/16
  tMatrix['WCHR', 'CCWB', ] <- c(eTen, wrbTen, wrbTen)/8


  wrbTen <- c(0, 1-dW, 0, 1-dW, 0, dW*drW, dW*(1-drW), 0, dW*drW, dW*(1-drW))
  tMatrix['WCHB', 'WWWW', ] <- c(wrbTen, wrbTen, eTen)/4
  tMatrix['WCHB', 'WCWW', ] <- c(wrbTen, 2*wrbTen, wrbTen)/8
  tMatrix['WCHB', 'CCWW', ] <- c(eTen, wrbTen, wrbTen)/4

  wrbTen <- c(0, 1-dW, 0, 1-dW, 1+dW*dhW, dW*(1-dhW)*drW,
              1+dW*dhW + dW*(1-dhW)*(1-drW), 0, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW))
  tMatrix['WCHB', 'WWWH', ] <- c(wrbTen, wrbTen, eTen)/8

  wrbTen <- c(0, 1-dW, 0, 1-dW, 0, 1+dW*drW, dW*(1-drW), 0, 1+dW*drW, dW*(1-drW))
  tMatrix['WCHB', 'WWWR', ] <- c(wrbTen, wrbTen, eTen)/8
  tMatrix['WCHB', 'WCWR', ] <- c(wrbTen, 2*wrbTen, wrbTen)/16
  tMatrix['WCHB', 'CCWR', ] <- c(eTen, wrbTen, wrbTen)/8

  wrbTen <- c(0, 1-dW, 0, 1-dW, 0, dW*drW, 1+dW*(1-drW), 0, dW*drW, 1+dW*(1-drW))
  tMatrix['WCHB', 'WWWB', ] <- c(wrbTen, wrbTen, eTen)/8
  tMatrix['WCHB', 'WCWB', ] <- c(wrbTen, 2*wrbTen, wrbTen)/16
  tMatrix['WCHB', 'CCWB', ] <- c(eTen, wrbTen, wrbTen)/8


  # double deposition stuff
  wrbTen <- c(0, 1-ddW, 0, 0, 0, ddW*ddrW, ddW*(1-ddrW), 0, 0, 0)
  tMatrix['CCHH', 'WWWW', ] <- c(eTen, wrbTen, eTen)
  tMatrix['CCHH', 'WCWW', ] <- c(eTen, wrbTen, wrbTen)/2
  tMatrix['CCHH', 'CCWW', ] <- c(eTen, eTen, wrbTen)

  wrbTen <- c(0, 1-ddW, 0, 0, 1+ddW*ddhW, ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW), 0, 0, 0)
  tMatrix['CCHH', 'WWWH', ] <- c(eTen, wrbTen, eTen)/2

  wrbTen <- c(0, 1-ddW, 0, 0, 0, 1+ddW*ddrW, ddW*(1-ddrW), 0, 0, 0)
  tMatrix['CCHH', 'WWWR', ] <- c(eTen, wrbTen, eTen)/2
  tMatrix['CCHH', 'WCWR', ] <- c(eTen, wrbTen, wrbTen)/4
  tMatrix['CCHH', 'CCWR', ] <- c(eTen, eTen, wrbTen)/2

  wrbTen <- c(0, 1-ddW, 0, 0, 0, ddW*ddrW, 1+ddW*(1-ddrW), 0, 0, 0)
  tMatrix['CCHH', 'WWWB', ] <- c(eTen, wrbTen, eTen)/2
  tMatrix['CCHH', 'WCWB', ] <- c(eTen, wrbTen, wrbTen)/4
  tMatrix['CCHH', 'CCWB', ] <- c(eTen, eTen, wrbTen)/2


  wrbTen <- c(0, 1-ddW, 1-ddW, 0, 0, ddW*ddrW, ddW*(1-ddrW), ddW*ddrW, ddW*(1-ddrW), 0)
  tMatrix['CCHR', 'WWWW', ] <- c(eTen, wrbTen, eTen)/2
  tMatrix['CCHR', 'WCWW', ] <- c(eTen, wrbTen, wrbTen)/4
  tMatrix['CCHR', 'CCWW', ] <- c(eTen, eTen, wrbTen)/2

  wrbTen <- c(0, 1-ddW, 1-ddW, 0, 1+ddW*ddhW, 1+ddW*ddhW + ddW*(1-ddhW)*ddrW,
              ddW*(1-ddhW)*(1-ddrW), ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW), 0)
  tMatrix['CCHR', 'WWWH', ] <- c(eTen, wrbTen, eTen)/4

  wrbTen <- c(0, 1-ddW, 1-ddW, 0, 0, 1+ddW*ddrW, ddW*(1-ddrW), 1+ddW*ddrW, ddW*(1-ddrW), 0)
  tMatrix['CCHR', 'WWWR', ] <- c(eTen, wrbTen, eTen)/4
  tMatrix['CCHR', 'WCWR', ] <- c(eTen, wrbTen, wrbTen)/8
  tMatrix['CCHR', 'CCWR', ] <- c(eTen, eTen, wrbTen)/4

  wrbTen <- c(0, 1-ddW, 1-ddW, 0, 0, ddW*ddrW, 1+ddW*(1-ddrW), ddW*ddrW, 1+ddW*(1-ddrW), 0)
  tMatrix['CCHR', 'WWWB', ] <- c(eTen, wrbTen, eTen)/4
  tMatrix['CCHR', 'WCWB', ] <- c(eTen, wrbTen, wrbTen)/8
  tMatrix['CCHR', 'CCWB', ] <- c(eTen, eTen, wrbTen)/4


  wrbTen <- c(0, 1-ddW, 0, 1-ddW, 0, ddW*ddrW, ddW*(1-ddrW), 0, ddW*ddrW, ddW*(1-ddrW))
  tMatrix['CCHB', 'WWWW', ] <- c(eTen, wrbTen, eTen)/2
  tMatrix['CCHB', 'WCWW', ] <- c(eTen, wrbTen, wrbTen)/4
  tMatrix['CCHB', 'CCWW', ] <- c(eTen, eTen, wrbTen)/2

  wrbTen <- c(0, 1-ddW, 0, 1-ddW, 1+ddW*ddhW, ddW*(1-ddhW)*ddrW,
              1+ddW*ddhW + ddW*(1-ddhW)*(1-ddrW), 0, ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW))
  tMatrix['CCHB', 'WWWH', ] <- c(eTen, wrbTen, eTen)/4

  wrbTen <- c(0, 1-ddW, 0, 1-ddW, 0, 1+ddW*ddrW, ddW*(1-ddrW), 0, 1+ddW*ddrW, ddW*(1-ddrW))
  tMatrix['CCHB', 'WWWR', ] <- c(eTen, wrbTen, eTen)/4
  tMatrix['CCHB', 'WCWR', ] <- c(eTen, wrbTen, wrbTen)/8
  tMatrix['CCHB', 'CCWR', ] <- c(eTen, eTen, wrbTen)/4

  wrbTen <- c(0, 1-ddW, 0, 1-ddW, 0, ddW*ddrW, 1+ddW*(1-ddrW), 0, ddW*ddrW, 1+ddW*(1-ddrW))
  tMatrix['CCHB', 'WWWB', ] <- c(eTen, wrbTen, eTen)/4
  tMatrix['CCHB', 'WCWB', ] <- c(eTen, wrbTen, wrbTen)/8
  tMatrix['CCHB', 'CCWB', ] <- c(eTen, eTen, wrbTen)/4


  ## female specific homing (with any deposition stuff)
  wrbTen <- c((1-cF)*(1-dW), (1+cF*chF)*(1-dW), (1-cF)*dW*drW + cF*(1-chF)*crF*(1-dW),
              (1-cF)*dW*(1-drW) + cF*(1-chF)*(1-crF)*(1-dW), 0, (1+cF*chF)*dW*drW,
              (1+cF*chF)*dW*(1-drW), cF*(1-chF)*crF*dW*drW,
              cF*(1-chF)*crF*dW*(1-drW) + cF*(1-chF)*(1-crF)*dW*drW, cF*(1-chF)*(1-crF)*dW*(1-drW) )
  tMatrix['WCWH', 'WWWW', ] <- c(wrbTen, wrbTen, eTen)/4
  tMatrix['WCWH', 'WCWW', ] <- c(wrbTen, 2*wrbTen, wrbTen)/8
  tMatrix['WCWH', 'CCWW', ] <- c(eTen, wrbTen, wrbTen)/4

  wrbTen <- c((1-cF)*(1-dW), (1-cF)*(1+dW*dhW) + (1+cF*chF)*(1-dW), (1-cF)*dW*(1-dhW)*drW + cF*(1-chF)*crF*(1-dW),
              (1-cF)*dW*(1-dhW)*(1-drW) + cF*(1-chF)*(1-crF)*(1-dW), (1+cF*chF)*(1+dW*dhW),
              (1+cF*chF)*dW*(1-dhW)*drW + cF*(1-chF)*crF*(1+dW*dhW),
              (1+cF*chF)*dW*(1-dhW)*(1-drW) + cF*(1-chF)*(1-crF)*(1+dW*dhW), cF*(1-chF)*crF*dW*(1-dhW)*drW,
              cF*(1-chF)*crF*dW*(1-dhW)*(1-drW) + cF*(1-chF)*(1-crF)*dW*(1-dhW)*drW,
              cF*(1-chF)*(1-crF)*dW*(1-dhW)*(1-drW) )
  tMatrix['WCWH', 'WWWH', ] <- c(wrbTen, wrbTen, eTen)/8

  wrbTen <- c((1-cF)*(1-dW), (1+cF*chF)*(1-dW), (1-cF)*(1+dW*drW) + cF*(1-chF)*crF*(1-dW),
              (1-cF)*dW*(1-drW) + cF*(1-chF)*(1-crF)*(1-dW), 0, (1+cF*chF)*(1+dW*drW),
              (1+cF*chF)*dW*(1-drW), cF*(1-chF)*crF*(1+dW*drW),
              cF*(1-chF)*crF*dW*(1-drW) + cF*(1-chF)*(1-crF)*(1+dW*drW), cF*(1-chF)*(1-crF)*dW*(1-drW) )
  tMatrix['WCWH', 'WWWR', ] <- c(wrbTen, wrbTen, eTen)/8
  tMatrix['WCWH', 'WCWR', ] <- c(wrbTen, 2*wrbTen, wrbTen)/16
  tMatrix['WCWH', 'CCWR', ] <- c(eTen, wrbTen, wrbTen)/8

  wrbTen <- c((1-cF)*(1-dW), (1+cF*chF)*(1-dW), (1-cF)*dW*drW + cF*(1-chF)*crF*(1-dW),
              (1-cF)*(1+dW*(1-drW)) + cF*(1-chF)*(1-crF)*(1-dW), 0, (1+cF*chF)*dW*drW,
              (1+cF*chF)*(1+dW*(1-drW)), cF*(1-chF)*crF*dW*drW,
              cF*(1-chF)*crF*(1+dW*(1-drW)) + cF*(1-chF)*(1-crF)*dW*drW, cF*(1-chF)*(1-crF)*(1+dW*(1-drW)) )
  tMatrix['WCWH', 'WWWB', ] <- c(wrbTen, wrbTen, eTen)/8
  tMatrix['WCWH', 'WCWB', ] <- c(wrbTen, 2*wrbTen, wrbTen)/16
  tMatrix['WCWH', 'CCWB', ] <- c(eTen, wrbTen, wrbTen)/8

  tMatrix['WCWH', 'WWHH', c('WWWH','WWHH','WWHR','WWHB',
                            'WCWH','WCHH','WCHR','WCHB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF),
                                                               1-cF, 1+cF*chF, cF*(1-chF)*crF,cF*(1-chF)*(1-crF))/4
  tMatrix['WCWH', 'WWHR', c('WWWH','WWHH','WWHR','WWHB','WWWR','WWRR','WWRB',
                            'WCWH','WCHH','WCHR','WCHB','WCWR','WCRR','WCRB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF + 1+cF*chF,  cF*(1-chF)*(1-crF), 1-cF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF),
                                                                                    1-cF, 1+cF*chF, cF*(1-chF)*crF + 1+cF*chF,  cF*(1-chF)*(1-crF), 1-cF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/8
  tMatrix['WCWH', 'WWHB', c('WWWH','WWHH','WWHR','WWHB','WWWB','WWRB','WWBB',
                            'WCWH','WCHH','WCHR','WCHB','WCWB','WCRB','WCBB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF) + 1+cF*chF, 1-cF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF),
                                                                                    1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF) + 1+cF*chF, 1-cF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/8
  tMatrix['WCWH', 'WWRR', c('WWWR','WWHR','WWRR','WWRB',
                            'WCWR','WCHR','WCRR','WCRB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF),
                                                               1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/4
  tMatrix['WCWH', 'WWRB', c('WWWR','WWHR','WWRR','WWRB','WWWB','WWHB','WWBB',
                            'WCWR','WCHR','WCRR','WCRB','WCWB','WCHB','WCBB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF) + cF*(1-chF)*crF, 1-cF, 1+cF*chF, cF*(1-chF)*(1-crF),
                                                                                    1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF) + cF*(1-chF)*crF, 1-cF, 1+cF*chF, cF*(1-chF)*(1-crF))/8
  tMatrix['WCWH', 'WWBB', c('WWWB','WWHB','WWRB','WWBB',
                            'WCWB','WCHB','WCRB','WCBB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF),
                                                               1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/4

  tMatrix['WCWH', 'WCHH', c('WWWH','WWHH','WWHR','WWHB',
                            'WCWH','WCHH','WCHR','WCHB',
                            'CCWH','CCHH','CCHR','CCHB')] <- c(c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/2,
                                                               1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF),
                                                               c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/2)/4
  tMatrix['WCWH', 'WCHR', c('WWWH','WWHH','WWHR','WWHB','WWWR','WWRR','WWRB',
                            'WCWH','WCHH','WCHR','WCHB','WCWR','WCRR','WCRB',
                            'CCWH','CCHH','CCHR','CCHB','CCWR','CCRR','CCRB')] <- c(c(1-cF, 1+cF*chF, cF*(1-chF)*crF + 1+cF*chF, cF*(1-chF)*(1-crF), 1-cF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/2,
                                                                                    1-cF, 1+cF*chF, cF*(1-chF)*crF + 1+cF*chF, cF*(1-chF)*(1-crF), 1-cF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF),
                                                                                    c(1-cF, 1+cF*chF, cF*(1-chF)*crF + 1+cF*chF, cF*(1-chF)*(1-crF), 1-cF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/2)/8
  tMatrix['WCWH', 'WCHB', c('WWWH','WWHH','WWHR','WWHB','WWWB','WWRB','WWBB',
                            'WCWH','WCHH','WCHR','WCHB','WCWB','WCRB','WCBB',
                            'CCWH','CCHH','CCHR','CCHB','CCWB','CCRB','CCBB')] <- c(c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF) + 1+cF*chF, 1-cF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/2,
                                                                                    1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF) + 1+cF*chF, 1-cF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF),
                                                                                    c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF) + 1+cF*chF, 1-cF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/2)/8
  tMatrix['WCWH', 'WCRR', c('WWWR','WWHR','WWRR','WWRB',
                            'WCWR','WCHR','WCRR','WCRB',
                            'CCWR','CCHR','CCRR','CCRB')] <- c(c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/2,
                                                               1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF),
                                                               c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/2)/4
  tMatrix['WCWH', 'WCRB', c('WWWR','WWHR','WWRR','WWRB','WWWB','WWHB','WWBB',
                            'WCWR','WCHR','WCRR','WCRB','WCWB','WCHB','WCBB',
                            'CCWR','CCHR','CCRR','CCRB','CCWB','CCHB','CCBB')] <- c(c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF)+ cF*(1-chF)*crF, 1-cF, 1+cF*chF, cF*(1-chF)*(1-crF))/2,
                                                                                    1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF)+ cF*(1-chF)*crF, 1-cF, 1+cF*chF, cF*(1-chF)*(1-crF),
                                                                                    c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF)+ cF*(1-chF)*crF, 1-cF, 1+cF*chF, cF*(1-chF)*(1-crF))/2)/8
  tMatrix['WCWH', 'WCBB', c('WWWB','WWHB','WWRB','WWBB',
                            'WCWB','WCHB','WCRB','WCBB',
                            'CCWB','CCHB','CCRB','CCBB')] <- c(c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/2,
                                                               1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF),
                                                               c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/2)/4

  tMatrix['WCWH', 'CCHH', c('WCWH','WCHH','WCHR','WCHB',
                            'CCWH','CCHH','CCHR','CCHB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF),
                                                               1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/4
  tMatrix['WCWH', 'CCHR', c('WCWH','WCHH','WCHR','WCHB','WCWR','WCRR','WCRB',
                            'CCWH','CCHH','CCHR','CCHB','CCWR','CCRR','CCRB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF + 1+cF*chF, cF*(1-chF)*(1-crF), 1-cF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF),
                                                                                    1-cF, 1+cF*chF, cF*(1-chF)*crF + 1+cF*chF, cF*(1-chF)*(1-crF), 1-cF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/8
  tMatrix['WCWH', 'CCHB', c('WCWH','WCHH','WCHR','WCHB','WCWB','WCRB','WCBB',
                            'CCWH','CCHH','CCHR','CCHB','CCWB','CCRB','CCBB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF) + 1+cF*chF, 1-cF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF),
                                                                                    1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF) + 1+cF*chF, 1-cF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/8
  tMatrix['WCWH', 'CCRR', c('WCWR','WCHR','WCRR','WCRB',
                            'CCWR','CCHR','CCRR','CCRB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF),
                                                               1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/4
  tMatrix['WCWH', 'CCRB', c('WCWR','WCHR','WCRR','WCRB','WCWB','WCHB','WCBB',
                            'CCWR','CCHR','CCRR','CCRB','CCWB','CCHB','CCBB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF) + cF*(1-chF)*crF, 1-cF, 1+cF*chF, cF*(1-chF)*(1-crF),
                                                                                    1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF) + cF*(1-chF)*crF, 1-cF, 1+cF*chF, cF*(1-chF)*(1-crF))/8
  tMatrix['WCWH', 'CCBB', c('WCWB','WCHB','WCRB','WCBB',
                            'CCWB','CCHB','CCRB','CCBB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF),
                                                               1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/4


  wrbTen <- c((1-ccF)*(1-ddW), (1+ccF*cchF)*(1-ddW), (1-ccF)*ddW*ddrW + ccF*(1-cchF)*ccrF*(1-ddW),
              (1-ccF)*ddW*(1-ddrW) + ccF*(1-cchF)*(1-ccrF)*(1-ddW), 0, (1+ccF*cchF)*ddW*ddrW,
              (1+ccF*cchF)*ddW*(1-ddrW), ccF*(1-cchF)*ccrF*ddW*ddrW,
              ccF*(1-cchF)*ccrF*ddW*(1-ddrW) + ccF*(1-cchF)*(1-ccrF)*ddW*ddrW, ccF*(1-cchF)*(1-ccrF)*ddW*(1-ddrW) )
  tMatrix['CCWH', 'WWWW', ] <- c(eTen, wrbTen, eTen)/2
  tMatrix['CCWH', 'WCWW', ] <- c(eTen, wrbTen, wrbTen)/4
  tMatrix['CCWH', 'CCWW', ] <- c(eTen, eTen, wrbTen)/2

  wrbTen <- c((1-ccF)*(1-ddW), (1-ccF)*(1+ddW*ddhW) + (1+ccF*cchF)*(1-ddW), (1-ccF)*ddW*(1-ddhW)*ddrW + ccF*(1-cchF)*ccrF*(1-ddW),
              (1-ccF)*ddW*(1-ddhW)*(1-ddrW) + ccF*(1-cchF)*(1-ccrF)*(1-ddW), (1+ccF*cchF)*(1+ddW*ddhW),
              (1+ccF*cchF)*ddW*(1-ddhW)*ddrW + ccF*(1-cchF)*ccrF*(1+ddW*ddhW),
              (1+ccF*cchF)*ddW*(1-ddhW)*(1-ddrW) + ccF*(1-cchF)*(1-ccrF)*(1+ddW*ddhW), ccF*(1-cchF)*ccrF*ddW*(1-ddhW)*ddrW,
              ccF*(1-cchF)*ccrF*ddW*(1-ddhW)*(1-ddrW) + ccF*(1-cchF)*(1-ccrF)*ddW*(1-ddhW)*ddrW,
              ccF*(1-cchF)*(1-ccrF)*ddW*(1-ddhW)*(1-ddrW) )
  tMatrix['CCWH', 'WWWH', ] <- c(eTen, wrbTen, eTen)/4

  wrbTen <- c((1-ccF)*(1-ddW), (1+ccF*cchF)*(1-ddW), (1-ccF)*(1+ddW*ddrW) + ccF*(1-cchF)*ccrF*(1-ddW),
              (1-ccF)*ddW*(1-ddrW) + ccF*(1-cchF)*(1-ccrF)*(1-ddW), 0, (1+ccF*cchF)*(1+ddW*ddrW),
              (1+ccF*cchF)*ddW*(1-ddrW), ccF*(1-cchF)*ccrF*(1+ddW*ddrW),
              ccF*(1-cchF)*ccrF*ddW*(1-ddrW) + ccF*(1-cchF)*(1-ccrF)*(1+ddW*ddrW), ccF*(1-cchF)*(1-ccrF)*ddW*(1-ddrW) )
  tMatrix['CCWH', 'WWWR', ] <- c(eTen, wrbTen, eTen)/4
  tMatrix['CCWH', 'WCWR', ] <- c(eTen, wrbTen, wrbTen)/8
  tMatrix['CCWH', 'CCWR', ] <- c(eTen, eTen, wrbTen)/4

  wrbTen <- c((1-ccF)*(1-ddW), (1+ccF*cchF)*(1-ddW), (1-ccF)*ddW*ddrW + ccF*(1-cchF)*ccrF*(1-ddW),
              (1-ccF)*(1+ddW*(1-ddrW)) + ccF*(1-cchF)*(1-ccrF)*(1-ddW), 0, (1+ccF*cchF)*ddW*ddrW,
              (1+ccF*cchF)*(1+ddW*(1-ddrW)), ccF*(1-cchF)*ccrF*ddW*ddrW,
              ccF*(1-cchF)*ccrF*(1+ddW*(1-ddrW)) + ccF*(1-cchF)*(1-ccrF)*ddW*ddrW, ccF*(1-cchF)*(1-ccrF)*(1+ddW*(1-ddrW)) )
  tMatrix['CCWH', 'WWWB', ] <- c(eTen, wrbTen, eTen)/4
  tMatrix['CCWH', 'WCWB', ] <- c(eTen, wrbTen, wrbTen)/8
  tMatrix['CCWH', 'CCWB', ] <- c(eTen, eTen, wrbTen)/4


  tMatrix['CCWH', 'WWHH', c('WCWH','WCHH','WCHR','WCHB')] <- c(1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF))/2
  tMatrix['CCWH', 'WWHR', c('WCWH','WCHH','WCHR','WCHB','WCWR','WCRR','WCRB')] <- c(1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF + 1+ccF*cchF, ccF*(1-cchF)*(1-ccrF), 1-ccF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF))/4
  tMatrix['CCWH', 'WWHB', c('WCWH','WCHH','WCHR','WCHB','WCWB','WCRB','WCBB')] <- c(1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF) + 1+ccF*cchF, 1-ccF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF))/4
  tMatrix['CCWH', 'WWRR', c('WCWR','WCHR','WCRR','WCRB')] <- c(1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF))/2
  tMatrix['CCWH', 'WWRB', c('WCWR','WCHR','WCRR','WCRB','WCWB','WCHB','WCBB')] <- c(1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF) + ccF*(1-cchF)*ccrF, 1-ccF, 1+ccF*cchF, ccF*(1-cchF)*(1-ccrF))/4
  tMatrix['CCWH', 'WWBB', c('WCWB','WCHB','WCRB','WCBB')] <- c(1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF))/2


  tMatrix['CCWH', 'WCHH', c('WCWH','WCHH','WCHR','WCHB',
                            'CCWH','CCHH','CCHR','CCHB')] <- c(1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF),
                                                               1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF))/4
  tMatrix['CCWH', 'WCHR', c('WCWH','WCHH','WCHR','WCHB','WCWR','WCRR','WCRB',
                            'CCWH','CCHH','CCHR','CCHB','CCWR','CCRR','CCRB')] <- c(1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF + 1+ccF*cchF, ccF*(1-cchF)*(1-ccrF), 1-ccF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF),
                                                                                    1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF + 1+ccF*cchF, ccF*(1-cchF)*(1-ccrF), 1-ccF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF))/8
  tMatrix['CCWH', 'WCHB', c('WCWH','WCHH','WCHR','WCHB','WCWB','WCRB','WCBB',
                            'CCWH','CCHH','CCHR','CCHB','CCWB','CCRB','CCBB')] <- c(1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF) + 1+ccF*cchF, 1-ccF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF),
                                                                                    1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF) + 1+ccF*cchF, 1-ccF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF))/8
  tMatrix['CCWH', 'WCRR', c('WCWR','WCHR','WCRR','WCRB',
                            'CCWR','CCHR','CCRR','CCRB')] <- c(1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF),
                                                               1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF))/4
  tMatrix['CCWH', 'WCRB', c('WCWR','WCHR','WCRR','WCRB','WCWB','WCHB','WCBB',
                            'CCWR','CCHR','CCRR','CCRB','CCWB','CCHB','CCBB')] <- c(1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF) + ccF*(1-cchF)*ccrF, 1-ccF, 1+ccF*cchF, ccF*(1-cchF)*(1-ccrF),
                                                                                    1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF) + ccF*(1-cchF)*ccrF, 1-ccF, 1+ccF*cchF, ccF*(1-cchF)*(1-ccrF))/8
  tMatrix['CCWH', 'WCBB', c('WCWB','WCHB','WCRB','WCBB',
                            'CCWB','CCHB','CCRB','CCBB')] <- c(1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF),
                                                               1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF))/4


  tMatrix['CCWH', 'CCHH', c('CCWH','CCHH','CCHR','CCHB')] <- c(1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF))/2
  tMatrix['CCWH', 'CCHR', c('CCWH','CCHH','CCHR','CCHB','CCWR','CCRR','CCRB')] <- c(1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF + 1+ccF*cchF, ccF*(1-cchF)*(1-ccrF), 1-ccF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF))/4
  tMatrix['CCWH', 'CCHB', c('CCWH','CCHH','CCHR','CCHB','CCWB','CCRB','CCBB')] <- c(1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF) + 1+ccF*cchF, 1-ccF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF))/4
  tMatrix['CCWH', 'CCRR', c('CCWR','CCHR','CCRR','CCRB')] <- c(1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF))/2
  tMatrix['CCWH', 'CCRB', c('CCWR','CCHR','CCRR','CCRB','CCWB','CCHB','CCBB')] <- c(1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF) + ccF*(1-cchF)*ccrF, 1-ccF, 1+ccF*cchF, ccF*(1-cchF)*(1-ccrF))/4
  tMatrix['CCWH', 'CCBB', c('CCWB','CCHB','CCRB','CCBB')] <- c(1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF))/2


  # mixed sex homing/deposition stuff
  wrbTen <- c((1-cF)*(1-cM)*(1-dW), (1+cF*chF)*(1-cM)*(1-dW) + (1-cF)*(1+cM*chM),
              (1-cF)*(1-cM)*dW*drW + (cF*(1-chF)*crF)*(1-cM)*(1-dW) + (1-cF)*(cM*(1-chM)*crM),
              (1-cF)*(1-cM)*dW*(1-drW) + (cF*(1-chF)*(1-crF))*(1-cM)*(1-dW) + (1-cF)*(cM*(1-chM)*(1-crM)),
              (1+cF*chF)*(1-cM)*dW*dhW + (1+cF*chF)*(1+cM*chM),
              (1+cF*chF)*(1-cM)*dW*(1-dhW)*drW + (cF*(1-chF)*crF)*(1+cM*chM) + (1+cF*chF)*(cM*(1-chM)*crM),
              (1+cF*chF)*(1-cM)*dW*(1-dhW)*(1-drW) + (cF*(1-chF)*(1-crF))*(1+cM*chM) + (1+cF*chF)*(cM*(1-chM)*(1-crM)),
              (cF*(1-chF)*crF)*(1-cM)*dW*drW + (cF*(1-chF)*crF)*(cM*(1-chM)*crM),
              (cF*(1-chF)*crF)*(1-cM)*dW*(1-drW) + (cF*(1-chF)*(1-crF))*(1-cM)*dW*drW + (cF*(1-chF)*(1-crF))*(cM*(1-chM)*crM) + (cF*(1-chF)*crF)*(cM*(1-chM)*(1-crM)),
              (cF*(1-chF)*(1-crF))*(1-cM)*dW*(1-drW) + (cF*(1-chF)*(1-crF))*(cM*(1-chM)*(1-crM)) )
  tMatrix['WCWH', 'WCWH', ] <- c(wrbTen, 2*wrbTen, wrbTen)/16

  wrbTen <- c((1-cF)*(1-ccM)*(1-dW), (1+cF*chF)*(1-ccM)*(1-dW) + (1-cF)*(1+ccM*cchM),
              (1-cF)*(1-ccM)*dW*drW + (cF*(1-chF)*crF)*(1-ccM)*(1-dW) + (1-cF)*(ccM*(1-cchM)*ccrM),
              (1-cF)*(1-ccM)*dW*(1-drW) + (cF*(1-chF)*(1-crF))*(1-ccM)*(1-dW) + (1-cF)*(ccM*(1-cchM)*(1-ccrM)),
              (1+cF*chF)*(1-ccM)*dW*dhW + (1+cF*chF)*(1+ccM*cchM),
              (1+cF*chF)*(1-ccM)*dW*(1-dhW)*drW + (cF*(1-chF)*crF)*(1+ccM*cchM) + (1+cF*chF)*(ccM*(1-cchM)*ccrM),
              (1+cF*chF)*(1-ccM)*dW*(1-dhW)*(1-drW) + (cF*(1-chF)*(1-crF))*(1+ccM*cchM) + (1+cF*chF)*(ccM*(1-cchM)*(1-ccrM)),
              (cF*(1-chF)*crF)*(1-ccM)*dW*drW + (cF*(1-chF)*crF)*(ccM*(1-cchM)*ccrM),
              (cF*(1-chF)*crF)*(1-ccM)*dW*(1-drW) + (cF*(1-chF)*(1-crF))*(1-ccM)*dW*drW + (cF*(1-chF)*(1-crF))*(ccM*(1-cchM)*ccrM) + (cF*(1-chF)*crF)*(ccM*(1-cchM)*(1-ccrM)),
              (cF*(1-chF)*(1-crF))*(1-ccM)*dW*(1-drW) + (cF*(1-chF)*(1-crF))*(ccM*(1-cchM)*(1-ccrM)) )
  tMatrix['WCWH', 'CCWH', ] <- c(eTen, wrbTen, wrbTen)/8

  wrbTen <- c((1-ccF)*(1-cM)*(1-ddW), (1+ccF*cchF)*(1-cM)*(1-ddW) + (1-ccF)*(1+cM*chM),
              (1-ccF)*(1-cM)*ddW*ddrW + (ccF*(1-cchF)*ccrF)*(1-cM)*(1-ddW) + (1-ccF)*(cM*(1-chM)*crM),
              (1-ccF)*(1-cM)*ddW*(1-ddrW) + (ccF*(1-cchF)*(1-ccrF))*(1-cM)*(1-ddW) + (1-ccF)*(cM*(1-chM)*(1-crM)),
              (1+ccF*cchF)*(1-cM)*ddW*ddhW + (1+ccF*cchF)*(1+cM*chM),
              (1+ccF*cchF)*(1-cM)*ddW*(1-ddhW)*ddrW + (ccF*(1-cchF)*ccrF)*(1+cM*chM) + (1+ccF*cchF)*(cM*(1-chM)*crM),
              (1+ccF*cchF)*(1-cM)*ddW*(1-ddhW)*(1-ddrW) + (ccF*(1-cchF)*(1-ccrF))*(1+cM*chM) + (1+ccF*cchF)*(cM*(1-chM)*(1-crM)),
              (ccF*(1-cchF)*ccrF)*(1-cM)*ddW*ddrW + (ccF*(1-cchF)*ccrF)*(cM*(1-chM)*crM),
              (ccF*(1-cchF)*ccrF)*(1-cM)*ddW*(1-ddrW) + (ccF*(1-cchF)*(1-ccrF))*(1-cM)*ddW*ddrW + (ccF*(1-cchF)*(1-ccrF))*(cM*(1-chM)*crM) + (ccF*(1-cchF)*ccrF)*(cM*(1-chM)*(1-crM)),
              (ccF*(1-cchF)*(1-ccrF))*(1-cM)*ddW*(1-ddrW) + (ccF*(1-cchF)*(1-ccrF))*(cM*(1-chM)*(1-crM)) )
  tMatrix['CCWH', 'WCWH', ] <- c(eTen, wrbTen, wrbTen)/8

  wrbTen <- c((1-ccF)*(1-ccM)*(1-ddW), (1+ccF*cchF)*(1-ccM)*(1-ddW) + (1-ccF)*(1+ccM*cchM),
              (1-ccF)*(1-ccM)*ddW*ddrW + (ccF*(1-cchF)*ccrF)*(1-ccM)*(1-ddW) + (1-ccF)*(ccM*(1-cchM)*ccrM),
              (1-ccF)*(1-ccM)*ddW*(1-ddrW) + (ccF*(1-cchF)*(1-ccrF))*(1-ccM)*(1-ddW) + (1-ccF)*(ccM*(1-cchM)*(1-ccrM)),
              (1+ccF*cchF)*(1-ccM)*ddW*ddhW + (1+ccF*cchF)*(1+ccM*cchM),
              (1+ccF*cchF)*(1-ccM)*ddW*(1-ddhW)*ddrW + (ccF*(1-cchF)*ccrF)*(1+ccM*cchM) + (1+ccF*cchF)*(ccM*(1-cchM)*ccrM),
              (1+ccF*cchF)*(1-ccM)*ddW*(1-ddhW)*(1-ddrW) + (ccF*(1-cchF)*(1-ccrF))*(1+ccM*cchM) + (1+ccF*cchF)*(ccM*(1-cchM)*(1-ccrM)),
              (ccF*(1-cchF)*ccrF)*(1-ccM)*ddW*ddrW + (ccF*(1-cchF)*ccrF)*(ccM*(1-cchM)*ccrM),
              (ccF*(1-cchF)*ccrF)*(1-ccM)*ddW*(1-ddrW) + (ccF*(1-cchF)*(1-ccrF))*(1-ccM)*ddW*ddrW + (ccF*(1-cchF)*(1-ccrF))*(ccM*(1-cchM)*ccrM) + (ccF*(1-cchF)*ccrF)*(ccM*(1-cchM)*(1-ccrM)),
              (ccF*(1-cchF)*(1-ccrF))*(1-ccM)*ddW*(1-ddrW) + (ccF*(1-cchF)*(1-ccrF))*(ccM*(1-cchM)*(1-ccrM)) )
  tMatrix['CCWH', 'CCWH', ] <- c(eTen, eTen, wrbTen/4)


  ## male specific homing (with any deposition stuff)
  wrbTen <- c(0, (1-cM)*(1-dW), 0, 0, (1-cM)*dW*dhW + 1+cM*chM,
              (1-cM)*dW*(1-dhW)*drW + cM*(1-chM)*crM,
              (1-cM)*dW*(1-dhW)*(1-drW) + cM*(1-chM)*(1-crM), 0, 0, 0)
  tMatrix['WCHH', 'WCWH', ] <- c(wrbTen, 2*wrbTen, wrbTen)/8

  wrbTen <- c(0, (1-ccM)*(1-dW), 0, 0, (1-ccM)*dW*dhW + 1+ccM*cchM,
              (1-ccM)*dW*(1-dhW)*drW + ccM*(1-cchM)*ccrM,
              (1-ccM)*dW*(1-dhW)*(1-drW) + ccM*(1-cchM)*(1-ccrM), 0, 0, 0)
  tMatrix['WCHH', 'CCWH', ] <- c(eTen, wrbTen, wrbTen)/4

  wrbTen <- c(0, (1-cM)*(1-dW), (1-cM)*(1-dW), 0, (1-cM)*dW*dhW + 1+cM*chM,
              (1-cM)*dW*(1-dhW)*drW + cM*(1-chM)*crM  + 1+cM*chM,
              (1-cM)*dW*(1-dhW)*(1-drW) + cM*(1-chM)*(1-crM),
              (1-cM)*dW*drW + cM*(1-chM)*crM, (1-cM)*dW*(1-drW) + cM*(1-chM)*(1-crM), 0)
  tMatrix['WCHR', 'WCWH', ] <- c(wrbTen, 2*wrbTen, wrbTen)/16

  wrbTen <- c(0, (1-ccM)*(1-dW), (1-ccM)*(1-dW), 0, (1-ccM)*dW*dhW + 1+ccM*cchM,
              (1-ccM)*dW*(1-dhW)*drW + ccM*(1-cchM)*ccrM  + 1+ccM*cchM,
              (1-ccM)*dW*(1-dhW)*(1-drW) + ccM*(1-cchM)*(1-ccrM),
              (1-ccM)*dW*drW + ccM*(1-cchM)*ccrM, (1-ccM)*dW*(1-drW) + ccM*(1-cchM)*(1-ccrM), 0)
  tMatrix['WCHR', 'CCWH', ] <- c(eTen, wrbTen, wrbTen)/8

  wrbTen <- c(0, (1-cM)*(1-dW), 0, (1-cM)*(1-dW), (1-cM)*dW*dhW + 1+cM*chM,
              (1-cM)*dW*(1-dhW)*drW + cM*(1-chM)*crM,
              (1-cM)*dW*(1-dhW)*(1-drW) + cM*(1-chM)*(1-crM) + 1+cM*chM,
              0, cM*(1-chM)*crM + (1-cM)*dW*drW, cM*(1-chM)*(1-crM) + (1-cM)*dW*(1-drW) )
  tMatrix['WCHB', 'WCWH', ] <- c(wrbTen, 2*wrbTen, wrbTen)/16

  wrbTen <- c(0, (1-ccM)*(1-dW), 0, (1-ccM)*(1-dW), (1-ccM)*dW*dhW + 1+ccM*cchM,
              (1-ccM)*dW*(1-dhW)*drW + ccM*(1-cchM)*ccrM,
              (1-ccM)*dW*(1-dhW)*(1-drW) + ccM*(1-cchM)*(1-ccrM) + 1+ccM*cchM,
              0, ccM*(1-cchM)*ccrM + (1-ccM)*dW*drW, ccM*(1-cchM)*(1-ccrM) + (1-ccM)*dW*(1-drW) )
  tMatrix['WCHB', 'CCWH', ] <- c(eTen, wrbTen, wrbTen)/8


  wrbTen <- c(0, (1-cM)*(1-ddW), 0, 0, (1-cM)*ddW*ddhW + 1+cM*chM,
              (1-cM)*ddW*(1-ddhW)*ddrW + cM*(1-chM)*crM,
              (1-cM)*ddW*(1-ddhW)*(1-ddrW) + cM*(1-chM)*(1-crM), 0, 0, 0)
  tMatrix['CCHH', 'WCWH', ] <- c(eTen, wrbTen, wrbTen)/4

  wrbTen <- c(0, (1-ccM)*(1-ddW), 0, 0, (1-ccM)*ddW*ddhW + 1+ccM*cchM,
              (1-ccM)*ddW*(1-ddhW)*ddrW + ccM*(1-cchM)*ccrM,
              (1-ccM)*ddW*(1-ddhW)*(1-ddrW) + ccM*(1-cchM)*(1-ccrM), 0, 0, 0)
  tMatrix['CCHH', 'CCWH', ] <- c(eTen, eTen, wrbTen)/2

  wrbTen <- c(0, (1-cM)*(1-ddW), (1-cM)*(1-ddW), 0, (1-cM)*ddW*ddhW + 1+cM*chM,
              (1-cM)*ddW*(1-ddhW)*ddrW + cM*(1-chM)*crM  + 1+cM*chM,
              (1-cM)*ddW*(1-ddhW)*(1-ddrW) + cM*(1-chM)*(1-crM),
              (1-cM)*ddW*ddrW + cM*(1-chM)*crM, (1-cM)*ddW*(1-ddrW) + cM*(1-chM)*(1-crM), 0)
  tMatrix['CCHR', 'WCWH', ] <- c(eTen, wrbTen, wrbTen)/8

  wrbTen <- c(0, (1-ccM)*(1-ddW), (1-ccM)*(1-ddW), 0, (1-ccM)*ddW*ddhW + 1+ccM*cchM,
              (1-ccM)*ddW*(1-ddhW)*ddrW + ccM*(1-cchM)*ccrM  + 1+ccM*cchM,
              (1-ccM)*ddW*(1-ddhW)*(1-ddrW) + ccM*(1-cchM)*(1-ccrM),
              (1-ccM)*ddW*ddrW + ccM*(1-cchM)*ccrM, (1-ccM)*ddW*(1-ddrW) + ccM*(1-cchM)*(1-ccrM), 0)
  tMatrix['CCHR', 'CCWH', ] <- c(eTen, eTen, wrbTen)/4

  wrbTen <- c(0, (1-cM)*(1-ddW), 0, (1-cM)*(1-ddW), (1-cM)*ddW*ddhW + 1+cM*chM,
              (1-cM)*ddW*(1-ddhW)*ddrW + cM*(1-chM)*crM,
              (1-cM)*ddW*(1-ddhW)*(1-ddrW) + cM*(1-chM)*(1-crM) + 1+cM*chM,
              0, cM*(1-chM)*crM + (1-cM)*ddW*ddrW, cM*(1-chM)*(1-crM) + (1-cM)*ddW*(1-ddrW) )
  tMatrix['CCHB', 'WCWH', ] <- c(eTen, wrbTen, wrbTen)/8

  wrbTen <- c(0, (1-ccM)*(1-ddW), 0, (1-ccM)*(1-ddW), (1-ccM)*ddW*ddhW + 1+ccM*cchM,
              (1-ccM)*ddW*(1-ddhW)*ddrW + ccM*(1-cchM)*ccrM,
              (1-ccM)*ddW*(1-ddhW)*(1-ddrW) + ccM*(1-cchM)*(1-ccrM) + 1+ccM*cchM,
              0, ccM*(1-cchM)*ccrM + (1-ccM)*ddW*ddrW, ccM*(1-cchM)*(1-ccrM) + (1-ccM)*ddW*(1-ddrW) )
  tMatrix['CCHB', 'CCWH', ] <- c(eTen, eTen, wrbTen)/4


  tMatrix['WWWW', 'WCWH', c('WWWW','WWWH','WWWR','WWWB',
                            'WCWW','WCWH','WCWR','WCWB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/4
  tMatrix['WWWH', 'WCWH', c('WWWW','WWWH','WWWR','WWWB','WWHH','WWHR','WWHB',
                            'WCWW','WCWH','WCWR','WCWB','WCHH','WCHR','WCHB')] <- c(1-cM, 1+cM*chM + 1-cM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM), 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                                                    1-cM, 1+cM*chM + 1-cM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM), 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['WWWR', 'WCWH', c('WWWW','WWWH','WWWR','WWWB','WWHR','WWRR','WWRB',
                            'WCWW','WCWH','WCWR','WCWB','WCHR','WCRR','WCRB')]<- c(1-cM, 1+cM*chM, cM*(1-chM)*crM + 1-cM, cM*(1-chM)*(1-crM), 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                                                   1-cM, 1+cM*chM, cM*(1-chM)*crM + 1-cM, cM*(1-chM)*(1-crM), 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['WWWB', 'WCWH', c('WWWW','WWWH','WWWR','WWWB','WWHB','WWRB','WWBB',
                            'WCWW','WCWH','WCWR','WCWB','WCHB','WCRB','WCBB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + 1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                                                    1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + 1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['WWHH', 'WCWH', c('WWWH','WWHH','WWHR','WWHB',
                            'WCWH','WCHH','WCHR','WCHB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM,cM*(1-chM)*(1-crM))/4
  tMatrix['WWHR', 'WCWH',c('WWWH','WWHH','WWHR','WWHB','WWWR','WWRR','WWRB',
                            'WCWH','WCHH','WCHR','WCHB','WCWR','WCRR','WCRB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM + 1+cM*chM,  cM*(1-chM)*(1-crM), 1-cM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                                                    1-cM, 1+cM*chM, cM*(1-chM)*crM + 1+cM*chM,  cM*(1-chM)*(1-crM), 1-cM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['WWHB', 'WCWH', c('WWWH','WWHH','WWHR','WWHB','WWWB','WWRB','WWBB',
                            'WCWH','WCHH','WCHR','WCHB','WCWB','WCRB','WCBB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + 1+cM*chM, 1-cM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                                                    1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + 1+cM*chM, 1-cM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['WWRR', 'WCWH', c('WWWR','WWHR','WWRR','WWRB',
                            'WCWR','WCHR','WCRR','WCRB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/4
  tMatrix['WWRB', 'WCWH', c('WWWR','WWHR','WWRR','WWRB','WWWB','WWHB','WWBB',
                            'WCWR','WCHR','WCRR','WCRB','WCWB','WCHB','WCBB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + cM*(1-chM)*crM, 1-cM, 1+cM*chM, cM*(1-chM)*(1-crM),
                                                                                    1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + cM*(1-chM)*crM, 1-cM, 1+cM*chM, cM*(1-chM)*(1-crM))/8
  tMatrix['WWBB', 'WCWH', c('WWWB','WWHB','WWRB','WWBB',
                            'WCWB','WCHB','WCRB','WCBB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/4
  tMatrix['WCWW', 'WCWH', c('WWWW','WWWH','WWWR','WWWB',
                            'WCWW','WCWH','WCWR','WCWB',
                            'CCWW','CCWH','CCWR','CCWB')] <- c(c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/2,
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/2)/4
  tMatrix['WCWR', 'WCWH', c('WWWW','WWWH','WWWR','WWWB','WWHR','WWRR','WWRB',
                            'WCWW','WCWH','WCWR','WCWB','WCHR','WCRR','WCRB',
                            'CCWW','CCWH','CCWR','CCWB','CCHR','CCRR','CCRB')] <- c(c(1-cM, 1+cM*chM, cM*(1-chM)*crM + 1-cM, cM*(1-chM)*(1-crM), 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/2,
                                                                                    1-cM, 1+cM*chM, cM*(1-chM)*crM + 1-cM, cM*(1-chM)*(1-crM), 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                                                    c(1-cM, 1+cM*chM, cM*(1-chM)*crM + 1-cM, cM*(1-chM)*(1-crM), 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/2)/8
  tMatrix['WCWB', 'WCWH', c('WWWW','WWWH','WWWR','WWWB','WWHB','WWRB','WWBB',
                            'WCWW','WCWH','WCWR','WCWB','WCHB','WCRB','WCBB',
                            'CCWW','CCWH','CCWR','CCWB','CCHB','CCRB','CCBB')] <- c(c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + 1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/2,
                                                                                    1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + 1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                                                    c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + 1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/2)/8
  tMatrix['WCRR', 'WCWH', c('WWWR','WWHR','WWRR','WWRB',
                            'WCWR','WCHR','WCRR','WCRB',
                            'CCWR','CCHR','CCRR','CCRB')] <- c(c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/2,
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/2)/4
  tMatrix['WCRB', 'WCWH', c('WWWR','WWHR','WWRR','WWRB','WWWB','WWHB','WWBB',
                            'WCWR','WCHR','WCRR','WCRB','WCWB','WCHB','WCBB',
                            'CCWR','CCHR','CCRR','CCRB','CCWB','CCHB','CCBB')] <- c(c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM)+ cM*(1-chM)*crM, 1-cM, 1+cM*chM, cM*(1-chM)*(1-crM))/2,
                                                                                    1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM)+ cM*(1-chM)*crM, 1-cM, 1+cM*chM, cM*(1-chM)*(1-crM),
                                                                                    c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM)+ cM*(1-chM)*crM, 1-cM, 1+cM*chM, cM*(1-chM)*(1-crM))/2)/8
  tMatrix['WCBB', 'WCWH', c('WWWB','WWHB','WWRB','WWBB',
                            'WCWB','WCHB','WCRB','WCBB',
                            'CCWB','CCHB','CCRB','CCBB')] <- c(c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/2,
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/2)/4


  tMatrix['CCWW', 'WCWH', c('WCWW','WCWH','WCWR','WCWB',
                            'CCWW','CCWH','CCWR','CCWB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/4
  tMatrix['CCWR', 'WCWH', c('WCWW','WCWH','WCWR','WCWB','WCHR','WCRR','WCRB',
                            'CCWW','CCWH','CCWR','CCWB','CCHR','CCRR','CCRB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM + 1-cM, cM*(1-chM)*(1-crM), 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                                                    1-cM, 1+cM*chM, cM*(1-chM)*crM + 1-cM, cM*(1-chM)*(1-crM), 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['CCWB', 'WCWH', c('WCWW','WCWH','WCWR','WCWB','WCHB','WCRB','WCBB',
                            'CCWW','CCWH','CCWR','CCWB','CCHB','CCRB','CCBB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + 1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                                                    1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + 1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['CCRR', 'WCWH', c('WCWR','WCHR','WCRR','WCRB',
                            'CCWR','CCHR','CCRR','CCRB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/4
  tMatrix['CCRB', 'WCWH', c('WCWR','WCHR','WCRR','WCRB','WCWB','WCHB','WCBB',
                            'CCWR','CCHR','CCRR','CCRB','CCWB','CCHB','CCBB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + cM*(1-chM)*crM, 1-cM, 1+cM*chM, cM*(1-chM)*(1-crM),
                                                                                    1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + cM*(1-chM)*crM, 1-cM, 1+cM*chM, cM*(1-chM)*(1-crM))/8
  tMatrix['CCBB', 'WCWH', c('WCWB','WCHB','WCRB','WCBB',
                            'CCWB','CCHB','CCRB','CCBB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/4


  tMatrix['WWWW', 'CCWH', c('WCWW','WCWH','WCWR','WCWB')] <- c(1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM))/2
  tMatrix['WWWH', 'CCWH', c('WCWW','WCWH','WCWR','WCWB','WCHH','WCHR','WCHB')] <- c(1-ccM, 1+ccM*cchM + 1-ccM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM), 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM))/4
  tMatrix['WWWR', 'CCWH', c('WCWW','WCWH','WCWR','WCWB','WCHR','WCRR','WCRB')] <- c(1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM + 1-ccM, ccM*(1-cchM)*(1-ccrM), 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM))/4
  tMatrix['WWWB', 'CCWH', c('WCWW','WCWH','WCWR','WCWB','WCHB','WCRB','WCBB')] <- c(1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM) + 1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM))/4
  tMatrix['WWHH', 'CCWH', c('WCWH','WCHH','WCHR','WCHB')] <- c(1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM))/2
  tMatrix['WWHR', 'CCWH', c('WCWH','WCHH','WCHR','WCHB','WCWR','WCRR','WCRB')] <- c(1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM + 1+ccM*cchM, ccM*(1-cchM)*(1-ccrM), 1-ccM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM))/4
  tMatrix['WWHB', 'CCWH', c('WCWH','WCHH','WCHR','WCHB','WCWB','WCRB','WCBB')] <- c(1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM) + 1+ccM*cchM, 1-ccM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM))/4
  tMatrix['WWRR', 'CCWH', c('WCWR','WCHR','WCRR','WCRB')] <- c(1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM))/2
  tMatrix['WWRB', 'CCWH', c('WCWR','WCHR','WCRR','WCRB','WCWB','WCHB','WCBB')] <- c(1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM) + ccM*(1-cchM)*ccrM, 1-ccM, 1+ccM*cchM, ccM*(1-cchM)*(1-ccrM))/4
  tMatrix['WWBB', 'CCWH', c('WCWB','WCHB','WCRB','WCBB')] <- c(1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM))/2


  tMatrix['WCWW', 'CCWH', c('WCWW','WCWH','WCWR','WCWB',
                            'CCWW','CCWH','CCWR','CCWB')] <- c(1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM),
                                                               1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM))/4
  tMatrix['WCWR', 'CCWH', c('WCWW','WCWH','WCWR','WCWB','WCHR','WCRR','WCRB',
                            'CCWW','CCWH','CCWR','CCWB','CCHR','CCRR','CCRB')] <- c(1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM + 1-ccM, ccM*(1-cchM)*(1-ccrM), 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM),
                                                                                    1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM + 1-ccM, ccM*(1-cchM)*(1-ccrM), 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM))/8
  tMatrix['WCWB', 'CCWH', c('WCWW','WCWH','WCWR','WCWB','WCHB','WCRB','WCBB',
                            'CCWW','CCWH','CCWR','CCWB','CCHB','CCRB','CCBB')] <- c(1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM) + 1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM),
                                                                                    1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM) + 1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM))/8
  tMatrix['WCRR', 'CCWH', c('WCWR','WCHR','WCRR','WCRB',
                            'CCWR','CCHR','CCRR','CCRB')] <- c(1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM),
                                                               1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM))/4
  tMatrix['WCRB', 'CCWH', c('WCWR','WCHR','WCRR','WCRB','WCWB','WCHB','WCBB',
                            'CCWR','CCHR','CCRR','CCRB','CCWB','CCHB','CCBB')] <- c(1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM) + ccM*(1-cchM)*ccrM, 1-ccM, 1+ccM*cchM, ccM*(1-cchM)*(1-ccrM),
                                                                                    1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM) + ccM*(1-cchM)*ccrM, 1-ccM, 1+ccM*cchM, ccM*(1-cchM)*(1-ccrM))/8
  tMatrix['WCBB', 'CCWH', c('WCWB','WCHB','WCRB','WCBB',
                            'CCWB','CCHB','CCRB','CCBB')] <- c(1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM),
                                                               1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM))/4


  tMatrix['CCWW', 'CCWH', c('CCWW','CCWH','CCWR','CCWB')] <- c(1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM))/2
  tMatrix['CCWR', 'CCWH', c('CCWW','CCWH','CCWR','CCWB','CCHR','CCRR','CCRB')] <- c(1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM + 1-ccM, ccM*(1-cchM)*(1-ccrM), 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM))/4
  tMatrix['CCWB', 'CCWH', c('CCWW','CCWH','CCWR','CCWB','CCHB','CCRB','CCBB')] <- c(1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM) + 1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM))/4
  tMatrix['CCRR', 'CCWH', c('CCWR','CCHR','CCRR','CCRB')] <- c(1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM))/2
  tMatrix['CCRB', 'CCWH', c('CCWR','CCHR','CCRR','CCRB','CCWB','CCHB','CCBB')] <- c(1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM) + ccM*(1-cchM)*ccrM, 1-ccM, 1+ccM*cchM, ccM*(1-cchM)*(1-ccrM))/4
  tMatrix['CCBB', 'CCWH', c('CCWB','CCHB','CCRB','CCBB')] <- c(1-ccM, 1+ccM*cchM, ccM*(1-cchM)*ccrM, ccM*(1-cchM)*(1-ccrM))/2


  #protection from underflow errors
  tMatrix[tMatrix < .Machine$double.eps] <- 0


  ## initialize viability mask. No mother/father-specific death, so use basic mask
  viabilityMask <- array(data = 1, dim = c(size,size,size), dimnames = list(gtype, gtype, gtype))

  ## genotype-specific modifiers
  modifiers = cubeModifiers(gtype, eta = eta, phi = phi, omega = omega, xiF = xiF, xiM = xiM, s = s)

  ## put everytWing into a labeled list to return
  return(list(
    ih = tMatrix,
    tau = viabilityMask,
    genotypesID = gtype,
    genotypesN = size,
    wildType = "WWWW",
    eta = modifiers$eta,
    phi = modifiers$phi,
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = "CCHH"
  ))

}
