###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   Y-linked Split Drive with 2 Resistance Alleles
#   Héctor Sanchez, Jared Bennett, Sean Wu, John Marshall
#   October 2020
#   jared_bennett@berkeley.edu
#
###############################################################################

#' Inheritance Cube: Y-linked Split CRISPR Drive with 2 Resistance Alleles
#'
#' This is a Y-linked version of a split CRISPR drive. At the Y-locus is the Cas9,
#' inherited in a Mendelian fashion. At a second, unlinked, autosomal locus are
#' the gRNAs. When the two loci occur together (i.e. in males), the gRNAs drive,
#' with potential damaged alleles, but the Cas9 remains Mendelian.
#' This drive has 2 loci:
#'  * "Locus" 1, sex chromosomes, has 3 alleles:
#'    * X: Wild-type X chromosome
#'    * Y: Wild-type Y chromosome
#'    * C: Y chromosome with Cas9
#'  * Locus 2, autosomal locus, has 4 alleles:
#'    * W: Wild-type allele
#'    * G: gRNA allele
#'    * R: Functional or low-cost resistance allele
#'    * B: Non-functional or high-cost resistance allele
#'
#' @param cM Cutting efficiency in males
#' @param chM Homing efficiency in males
#' @param crM Resistance efficiency in males
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
cubeSplitDriveY <- function(cM = 1.0, chM = 0, crM = 0,
                            eta = NULL, phi = NULL, omega = NULL,
                            xiF = NULL, xiM = NULL, s = NULL){

  # # for testing
  # testVec <- runif(n = 3)
  #
  # cM <- testVec[1]; chM <- testVec[2]; crM <- testVec[3];
  #
  # # # this would run at the bottom, to check that all cells sum to 1
  # # #  size is defined below
  # # for(mID in 1:size){
  # #   for(fID in 1:size){
  # #     if((abs(sum(tMatrix[fID,mID,])-1)) > 1e-6){
  # #       print(paste0("Fail: mID= ",gtype[mID],", fID= ",gtype[fID]))
  # #     }
  # #   }
  # # }


  ## safety checks
  params <- c(cM, chM, crM)
  if(any(params>1) || any(params<0)){
    stop("Parameters are rates, they must be between 0 and 1.")
  }


  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  gtype <- c('XXWW','XXWG','XXWR','XXWB','XXGG','XXGR','XXGB','XXRR','XXRB','XXBB',
             'XYWW','XYWG','XYWR','XYWB','XYGG','XYGR','XYGB','XYRR','XYRB','XYBB',
             'XCWW','XCWG','XCWR','XCWB','XCGG','XCGR','XCGB','XCRR','XCRB','XCBB')
  size <- length(gtype) #because I use it several times
  tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix


  ## fill tMatrix with probabilities
  tMatrix['XXWW', 'XYWW', c('XXWW','XYWW')] <- c(1,1)/2
  tMatrix['XXWW', 'XYWG', c('XXWW','XXWG',
                            'XYWW','XYWG')] <- c(1,1,1,1)/4
  tMatrix['XXWW', 'XYWR', c('XXWW','XXWR',
                            'XYWW','XYWR')] <- c(1,1,1,1)/4
  tMatrix['XXWW', 'XYWB', c('XXWW','XXWB',
                            'XYWW','XYWB')] <- c(1,1,1,1)/4
  tMatrix['XXWW', 'XYGG', c('XXWG','XYWG')] <- c(1,1)/2
  tMatrix['XXWW', 'XYGR', c('XXWG','XXWR',
                            'XYWG','XYWR')] <- c(1,1,1,1)/4
  tMatrix['XXWW', 'XYGB', c('XXWG','XXWB',
                            'XYWG','XYWB')] <- c(1,1,1,1)/4
  tMatrix['XXWW', 'XYRR', c('XXWR','XYWR')] <- c(1,1)/2
  tMatrix['XXWW', 'XYRB', c('XXWR','XXWB',
                            'XYWR','XYWB')] <- c(1,1,1,1)/4
  tMatrix['XXWW', 'XYBB', c('XXWB','XYWB')] <- c(1,1)/2

  tMatrix['XXWW', 'XCWW', c('XXWW','XCWW')] <- c(1,1)/2
  tMatrix['XXWW', 'XCWG', c('XXWW','XXWG','XXWR','XXWB',
                            'XCWW','XCWG','XCWR','XCWB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/4
  tMatrix['XXWW', 'XCWR', c('XXWW','XXWR',
                            'XCWW','XCWR')] <- c(1,1,1,1)/4
  tMatrix['XXWW', 'XCWB', c('XXWW','XXWB',
                            'XCWW','XCWB')] <- c(1,1,1,1)/4
  tMatrix['XXWW', 'XCGG', c('XXWG','XCWG')] <- c(1,1)/2
  tMatrix['XXWW', 'XCGR', c('XXWG','XXWR',
                            'XCWG','XCWR')] <- c(1,1,1,1)/4
  tMatrix['XXWW', 'XCGB', c('XXWG','XXWB',
                            'XCWG','XCWB')] <- c(1,1,1,1)/4
  tMatrix['XXWW', 'XCRR', c('XXWR','XCWR')] <- c(1,1)/2
  tMatrix['XXWW', 'XCRB', c('XXWR','XXWB',
                            'XCWR','XCWB')] <- c(1,1,1,1)/4
  tMatrix['XXWW', 'XCBB', c('XXWB','XCWB')] <- c(1,1)/2


  tMatrix['XXWG', 'XYWW', c('XXWW','XXWG',
                            'XYWW','XYWG')] <- c(1,1,1,1)/4
  tMatrix['XXWG', 'XYWG', c('XXWW','XXWG','XXGG',
                            'XYWW','XYWG','XYGG')] <- c(1,2,1,1,2,1)/8
  tMatrix['XXWG', 'XYWR', c('XXWW','XXWG','XXWR','XXGR',
                            'XYWW','XYWG','XYWR','XYGR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWG', 'XYWB', c('XXWW','XXWG','XXWB','XXGB',
                            'XYWW','XYWG','XYWB','XYGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWG', 'XYGG', c('XXWG','XXGG',
                            'XYWG','XYGG')] <- c(1,1,1,1)/4
  tMatrix['XXWG', 'XYGR', c('XXWG','XXGG','XXWR','XXGR',
                            'XYWG','XYGG','XYWR','XYGR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWG', 'XYGB', c('XXWG','XXGG','XXWB','XXGB',
                            'XYWG','XYGG','XYWB','XYGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWG', 'XYRR', c('XXWR','XXGR',
                            'XYWR','XYGR')] <- c(1,1,1,1)/4
  tMatrix['XXWG', 'XYRB', c('XXWR','XXGR','XXWB','XXGB',
                            'XYWR','XYGR','XYWB','XYGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWG', 'XYBB', c('XXWB','XXGB',
                            'XYWB','XYGB')] <- c(1,1,1,1)/4

  tMatrix['XXWG', 'XCWW', c('XXWW','XXWG',
                            'XCWW','XCWG')] <- c(1,1,1,1)/4
  tMatrix['XXWG', 'XCWG', c('XXWW','XXWG','XXWR','XXWB','XXGG','XXGR','XXGB',
                            'XCWW','XCWG','XCWR','XCWB','XCGG','XCGR','XCGB')] <- c(1-cM, 1+cM*chM + 1-cM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM), 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                                                    1-cM, 1+cM*chM + 1-cM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM), 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['XXWG', 'XCWR', c('XXWW','XXWG','XXWR','XXGR',
                            'XCWW','XCWG','XCWR','XCGR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWG', 'XCWB', c('XXWW','XXWG','XXWB','XXGB',
                            'XCWW','XCWG','XCWB','XCGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWG', 'XCGG', c('XXWG','XXGG',
                            'XCWG','XCGG')] <- c(1,1,1,1)/4
  tMatrix['XXWG', 'XCGR', c('XXWG','XXGG','XXWR','XXGR',
                            'XCWG','XCGG','XCWR','XCGR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWG', 'XCGB', c('XXWG','XXGG','XXWB','XXGB',
                            'XCWG','XCGG','XCWB','XCGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWG', 'XCRR', c('XXWR','XXGR',
                            'XCWR','XCGR')] <- c(1,1,1,1)/4
  tMatrix['XXWG', 'XCRB', c('XXWR','XXGR','XXWB','XXGB',
                            'XCWR','XCGR','XCWB','XCGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWG', 'XCBB', c('XXWB','XXGB',
                            'XCWB','XCGB')] <- c(1,1,1,1)/4


  tMatrix['XXWR', 'XYWW', c('XXWW','XXWR',
                            'XYWW','XYWR')] <- c(1,1,1,1)/4
  tMatrix['XXWR', 'XYWG', c('XXWW','XXWR','XXWG','XXGR',
                            'XYWW','XYWR','XYWG','XYGR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWR', 'XYWR', c('XXWW','XXWR','XXRR',
                            'XYWW','XYWR','XYRR')] <- c(1,2,1,1,2,1)/8
  tMatrix['XXWR', 'XYWB', c('XXWW','XXWR','XXWB','XXRB',
                            'XYWW','XYWR','XYWB','XYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWR', 'XYGG', c('XXWG','XXGR',
                            'XYWG','XYGR')] <- c(1,1,1,1)/4
  tMatrix['XXWR', 'XYGR', c('XXWG','XXGR','XXWR','XXRR',
                            'XYWG','XYGR','XYWR','XYRR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWR', 'XYGB', c('XXWG','XXGR','XXWB','XXRB',
                            'XYWG','XYGR','XYWB','XYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWR', 'XYRR', c('XXWR','XXRR',
                            'XYWR','XYRR')] <- c(1,1,1,1)/4
  tMatrix['XXWR', 'XYRB', c('XXWR','XXRR','XXWB','XXRB',
                            'XYWR','XYRR','XYWB','XYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWR', 'XYBB', c('XXWB','XXRB',
                            'XYWB','XYRB')] <- c(1,1,1,1)/4

  tMatrix['XXWR', 'XCWW', c('XXWW','XXWR',
                            'XCWW','XCWR')] <- c(1,1,1,1)/4
  tMatrix['XXWR', 'XCWG', c('XXWW','XXWG','XXWR','XXWB','XXGR','XXRR','XXRB',
                            'XCWW','XCWG','XCWR','XCWB','XCGR','XCRR','XCRB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM + 1-cM, cM*(1-chM)*(1-crM), 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                                                    1-cM, 1+cM*chM, cM*(1-chM)*crM + 1-cM, cM*(1-chM)*(1-crM), 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['XXWR', 'XCWR', c('XXWW','XXWR','XXRR',
                            'XCWW','XCWR','XCRR')] <- c(1,2,1,1,2,1)/8
  tMatrix['XXWR', 'XCWB', c('XXWW','XXWR','XXWB','XXRB',
                            'XCWW','XCWR','XCWB','XCRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWR', 'XCGG', c('XXWG','XXGR',
                            'XCWG','XCGR')] <- c(1,1,1,1)/4
  tMatrix['XXWR', 'XCGR', c('XXWG','XXGR','XXWR','XXRR',
                            'XCWG','XCGR','XCWR','XCRR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWR', 'XCGB', c('XXWG','XXGR','XXWB','XXRB',
                            'XCWG','XCGR','XCWB','XCRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWR', 'XCRR', c('XXWR','XXRR',
                            'XCWR','XCRR')] <- c(1,1,1,1)/4
  tMatrix['XXWR', 'XCRB', c('XXWR','XXRR','XXWB','XXRB',
                            'XCWR','XCRR','XCWB','XCRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWR', 'XCBB', c('XXWB','XXRB',
                            'XCWB','XCRB')] <- c(1,1,1,1)/4


  tMatrix['XXWB', 'XYWW', c('XXWW','XXWB',
                            'XYWW','XYWB')] <- c(1,1,1,1)/4
  tMatrix['XXWB', 'XYWG', c('XXWW','XXWB','XXWG','XXGB',
                            'XYWW','XYWB','XYWG','XYGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWB', 'XYWR', c('XXWW','XXWB','XXWR','XXRB',
                            'XYWW','XYWB','XYWR','XYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWB', 'XYWB', c('XXWW','XXWB','XXBB',
                            'XYWW','XYWB','XYBB')] <- c(1,2,1,1,2,1)/8
  tMatrix['XXWB', 'XYGG', c('XXWG','XXGB',
                            'XYWG','XYGB')] <- c(1,1,1,1)/4
  tMatrix['XXWB', 'XYGR', c('XXWG','XXGB','XXWR','XXGR',
                            'XYWG','XYGB','XYWR','XYGR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWB', 'XYGB', c('XXWG','XXGB','XXWB','XXBB',
                            'XYWG','XYGB','XYWB','XYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWB', 'XYRR', c('XXWR','XXRB',
                            'XYWR','XYRB')] <- c(1,1,1,1)/4
  tMatrix['XXWB', 'XYRB', c('XXWR','XXRB','XXWB','XXBB',
                            'XYWR','XYRB','XYWB','XYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWB', 'XYBB', c('XXWB','XXBB',
                            'XYWB','XYBB')] <- c(1,1,1,1)/4

  tMatrix['XXWB', 'XCWW', c('XXWW','XXWB',
                            'XCWW','XCWB')] <- c(1,1,1,1)/4
  tMatrix['XXWB', 'XCWG', c('XXWW','XXWG','XXWR','XXWB','XXGB','XXRB','XXBB',
                            'XCWW','XCWG','XCWR','XCWB','XCGB','XCRB','XCBB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + 1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                                                    1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + 1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['XXWB', 'XCWR', c('XXWW','XXWB','XXWR','XXRB',
                            'XCWW','XCWB','XCWR','XCRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWB', 'XCWB', c('XXWW','XXWB','XXBB',
                            'XCWW','XCWB','XCBB')] <- c(1,2,1,1,2,1)/8
  tMatrix['XXWB', 'XCGG', c('XXWG','XXGB',
                            'XCWG','XCGB')] <- c(1,1,1,1)/4
  tMatrix['XXWB', 'XCGR', c('XXWG','XXGB','XXWR','XXGR',
                            'XCWG','XCGB','XCWR','XCGR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWB', 'XCGB', c('XXWG','XXGB','XXWB','XXBB',
                            'XCWG','XCGB','XCWB','XCBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWB', 'XCRR', c('XXWR','XXRB',
                            'XCWR','XCRB')] <- c(1,1,1,1)/4
  tMatrix['XXWB', 'XCRB', c('XXWR','XXRB','XXWB','XXBB',
                            'XCWR','XCRB','XCWB','XCBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWB', 'XCBB', c('XXWB','XXBB',
                            'XCWB','XCBB')] <- c(1,1,1,1)/4


  tMatrix['XXGG', 'XYWW', c('XXWG','XYWG')] <- c(1,1)/2
  tMatrix['XXGG', 'XYWG', c('XXWG','XXGG',
                            'XYWG','XYGG')] <- c(1,1,1,1)/4
  tMatrix['XXGG', 'XYWR', c('XXWG','XXGR',
                            'XYWG','XYGR')] <- c(1,1,1,1)/4
  tMatrix['XXGG', 'XYWB', c('XXWG','XXGB',
                            'XYWG','XYGB')] <- c(1,1,1,1)/4
  tMatrix['XXGG', 'XYGG', c('XXGG','XYGG')] <- c(1,1)/2
  tMatrix['XXGG', 'XYGR', c('XXGG','XXGR',
                            'XYGG','XYGR')] <- c(1,1,1,1)/4
  tMatrix['XXGG', 'XYGB', c('XXGG','XXGB',
                            'XYGG','XYGB')] <- c(1,1,1,1)/4
  tMatrix['XXGG', 'XYRR', c('XXGR','XYGR')] <- c(1,1)/2
  tMatrix['XXGG', 'XYRB', c('XXGR','XXGB',
                            'XYGR','XYGB')] <- c(1,1,1,1)/4
  tMatrix['XXGG', 'XYBB', c('XXGB','XYGB')] <- c(1,1)/2

  tMatrix['XXGG', 'XCWW', c('XXWG','XCWG')] <- c(1,1)/2
  tMatrix['XXGG', 'XCWG', c('XXWG','XXGG','XXGR','XXGB',
                            'XCWG','XCGG','XCGR','XCGB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/4
  tMatrix['XXGG', 'XCWR', c('XXWG','XXGR',
                            'XCWG','XCGR')] <- c(1,1,1,1)/4
  tMatrix['XXGG', 'XCWB', c('XXWG','XXGB',
                            'XCWG','XCGB')] <- c(1,1,1,1)/4
  tMatrix['XXGG', 'XCGG', c('XXGG','XCGG')] <- c(1,1)/2
  tMatrix['XXGG', 'XCGR', c('XXGG','XXGR',
                            'XCGG','XCGR')] <- c(1,1,1,1)/4
  tMatrix['XXGG', 'XCGB', c('XXGG','XXGB',
                            'XCGG','XCGB')] <- c(1,1,1,1)/4
  tMatrix['XXGG', 'XCRR', c('XXGR','XCGR')] <- c(1,1)/2
  tMatrix['XXGG', 'XCRB', c('XXGR','XXGB',
                            'XCGR','XCGB')] <- c(1,1,1,1)/4
  tMatrix['XXGG', 'XCBB', c('XXGB','XCGB')] <- c(1,1)/2


  tMatrix['XXGR', 'XYWW', c('XXWG','XXWR',
                            'XYWG','XYWR')] <- c(1,1,1,1)/4
  tMatrix['XXGR', 'XYWG', c('XXWG','XXWR','XXGG','XXGR',
                            'XYWG','XYWR','XYGG','XYGR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGR', 'XYWR', c('XXWG','XXWR','XXGR','XXRR',
                            'XYWG','XYWR','XYGR','XYRR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGR', 'XYWB', c('XXWG','XXWR','XXGB','XXRB',
                            'XYWG','XYWR','XYGB','XYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGR', 'XYGG', c('XXGG','XXGR',
                            'XYGG','XYGR')] <- c(1,1,1,1)/4
  tMatrix['XXGR', 'XYGR', c('XXGG','XXGR','XXRR',
                            'XYGG','XYGR','XYRR')] <- c(1,2,1,1,2,1)/8
  tMatrix['XXGR', 'XYGB', c('XXGG','XXGR','XXGB','XXRB',
                            'XYGG','XYGR','XYGB','XYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGR', 'XYRR', c('XXGR','XXRR',
                            'XYGR','XYRR')] <- c(1,1,1,1)/4
  tMatrix['XXGR', 'XYRB', c('XXGR','XXRR','XXGB','XXRB',
                            'XYGR','XYRR','XYGB','XYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGR', 'XYBB', c('XXGB','XXRB',
                            'XYGB','XYRB')] <- c(1,1,1,1)/4

  tMatrix['XXGR', 'XCWW', c('XXWG','XXWR',
                            'XCWG','XCWR')] <- c(1,1,1,1)/4
  tMatrix['XXGR', 'XCWG', c('XXWG','XXGG','XXGR','XXGB','XXWR','XXRR','XXRB',
                            'XCWG','XCGG','XCGR','XCGB','XCWR','XCRR','XCRB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM + 1+cM*chM, cM*(1-chM)*(1-crM), 1-cM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                                                    1-cM, 1+cM*chM, cM*(1-chM)*crM + 1+cM*chM, cM*(1-chM)*(1-crM), 1-cM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['XXGR', 'XCWR', c('XXWG','XXWR','XXGR','XXRR',
                            'XCWG','XCWR','XCGR','XCRR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGR', 'XCWB', c('XXWG','XXWR','XXGB','XXRB',
                            'XCWG','XCWR','XCGB','XCRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGR', 'XCGG', c('XXGG','XXGR',
                            'XCGG','XCGR')] <- c(1,1,1,1)/4
  tMatrix['XXGR', 'XCGR', c('XXGG','XXGR','XXRR',
                            'XCGG','XCGR','XCRR')] <- c(1,2,1,1,2,1)/8
  tMatrix['XXGR', 'XCGB', c('XXGG','XXGR','XXGB','XXRB',
                            'XCGG','XCGR','XCGB','XCRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGR', 'XCRR', c('XXGR','XXRR',
                            'XCGR','XCRR')] <- c(1,1,1,1)/4
  tMatrix['XXGR', 'XCRB', c('XXGR','XXRR','XXGB','XXRB',
                            'XCGR','XCRR','XCGB','XCRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGR', 'XCBB', c('XXGB','XXRB',
                            'XCGB','XCRB')] <- c(1,1,1,1)/4


  tMatrix['XXGB', 'XYWW', c('XXWG','XXWB',
                            'XYWG','XYWB')] <- c(1,1,1,1)/4
  tMatrix['XXGB', 'XYWG', c('XXWG','XXWB','XXGG','XXGB',
                            'XYWG','XYWB','XYGG','XYGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGB', 'XYWR', c('XXWG','XXWB','XXGR','XXRB',
                            'XYWG','XYWB','XYGR','XYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGB', 'XYWB', c('XXWG','XXWB','XXGB','XXBB',
                            'XYWG','XYWB','XYGB','XYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGB', 'XYGG', c('XXGG','XXGB',
                            'XYGG','XYGB')] <- c(1,1,1,1)/4
  tMatrix['XXGB', 'XYGR', c('XXGG','XXGB','XXGR','XXRB',
                            'XYGG','XYGB','XYGR','XYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGB', 'XYGB', c('XXGG','XXGB','XXBB',
                            'XYGG','XYGB','XYBB')] <- c(1,2,1,1,2,1)/8
  tMatrix['XXGB', 'XYRR', c('XXGR','XXRB',
                            'XYGR','XYRB')] <- c(1,1,1,1)/4
  tMatrix['XXGB', 'XYRB', c('XXGR','XXRB','XXGB','XXBB',
                            'XYGR','XYRB','XYGB','XYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGB', 'XYBB', c('XXGB','XXBB',
                            'XYGB','XYBB')] <- c(1,1,1,1)/4

  tMatrix['XXGB', 'XCWW', c('XXWG','XXWB',
                            'XCWG','XCWB')] <- c(1,1,1,1)/4
  tMatrix['XXGB', 'XCWG', c('XXWG','XXGG','XXGR','XXGB','XXWB','XXRB','XXBB',
                            'XCWG','XCGG','XCGR','XCGB','XCWB','XCRB','XCBB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + 1+cM*chM, 1-cM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                                                    1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + 1+cM*chM, 1-cM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['XXGB', 'XCWR', c('XXWG','XXWB','XXGR','XXRB',
                            'XCWG','XCWB','XCGR','XCRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGB', 'XCWB', c('XXWG','XXWB','XXGB','XXBB',
                            'XCWG','XCWB','XCGB','XCBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGB', 'XCGG', c('XXGG','XXGB',
                            'XCGG','XCGB')] <- c(1,1,1,1)/4
  tMatrix['XXGB', 'XCGR', c('XXGG','XXGB','XXGR','XXRB',
                            'XCGG','XCGB','XCGR','XCRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGB', 'XCGB', c('XXGG','XXGB','XXBB',
                            'XCGG','XCGB','XCBB')] <- c(1,2,1,1,2,1)/8
  tMatrix['XXGB', 'XCRR', c('XXGR','XXRB',
                            'XCGR','XCRB')] <- c(1,1,1,1)/4
  tMatrix['XXGB', 'XCRB', c('XXGR','XXRB','XXGB','XXBB',
                            'XCGR','XCRB','XCGB','XCBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGB', 'XCBB', c('XXGB','XXBB',
                            'XCGB','XCBB')] <- c(1,1,1,1)/4


  tMatrix['XXRR', 'XYWW', c('XXWR','XYWR')] <- c(1,1)/2
  tMatrix['XXRR', 'XYWG', c('XXWR','XXGR',
                            'XYWR','XYGR')] <- c(1,1,1,1)/4
  tMatrix['XXRR', 'XYWR', c('XXWR','XXRR',
                            'XYWR','XYRR')] <- c(1,1,1,1)/4
  tMatrix['XXRR', 'XYWB', c('XXWR','XXRB',
                            'XYWR','XYRB')] <- c(1,1,1,1)/4
  tMatrix['XXRR', 'XYGG', c('XXGR','XYGR')] <- c(1,1)/2
  tMatrix['XXRR', 'XYGR', c('XXGR','XXRR',
                            'XYGR','XYRR')] <- c(1,1,1,1)/4
  tMatrix['XXRR', 'XYGB', c('XXGR','XXRB',
                            'XYGR','XYRB')] <- c(1,1,1,1)/4
  tMatrix['XXRR', 'XYRR', c('XXRR','XYRR')] <- c(1,1)/2
  tMatrix['XXRR', 'XYRB', c('XXRR','XXRB',
                            'XYRR','XYRB')] <- c(1,1,1,1)/4
  tMatrix['XXRR', 'XYBB', c('XXRB','XYRB')] <- c(1,1)/2

  tMatrix['XXRR', 'XCWW', c('XXWR','XCWR')] <- c(1,1)/2
  tMatrix['XXRR', 'XCWG', c('XXWR','XXGR','XXRR','XXRB',
                            'XCWR','XCGR','XCRR','XCRB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/4
  tMatrix['XXRR', 'XCWR', c('XXWR','XXRR',
                            'XCWR','XCRR')] <- c(1,1,1,1)/4
  tMatrix['XXRR', 'XCWB', c('XXWR','XXRB',
                            'XCWR','XCRB')] <- c(1,1,1,1)/4
  tMatrix['XXRR', 'XCGG', c('XXGR','XCGR')] <- c(1,1)/2
  tMatrix['XXRR', 'XCGR', c('XXGR','XXRR',
                            'XCGR','XCRR')] <- c(1,1,1,1)/4
  tMatrix['XXRR', 'XCGB', c('XXGR','XXRB',
                            'XCGR','XCRB')] <- c(1,1,1,1)/4
  tMatrix['XXRR', 'XCRR', c('XXRR','XCRR')] <- c(1,1)/2
  tMatrix['XXRR', 'XCRB', c('XXRR','XXRB',
                            'XCRR','XCRB')] <- c(1,1,1,1)/4
  tMatrix['XXRR', 'XCBB', c('XXRB','XCRB')] <- c(1,1)/2


  tMatrix['XXRB', 'XYWW', c('XXWR','XXWB',
                            'XYWR','XYWB')] <- c(1,1,1,1)/4
  tMatrix['XXRB', 'XYWG', c('XXWR','XXWB','XXGR','XXGB',
                            'XYWR','XYWB','XYGR','XYGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXRB', 'XYWR', c('XXWR','XXWB','XXRR','XXRB',
                            'XYWR','XYWB','XYRR','XYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXRB', 'XYWB', c('XXWR','XXWB','XXRB','XXBB',
                            'XYWR','XYWB','XYRB','XYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXRB', 'XYGG', c('XXGR','XXGB',
                            'XYGR','XYGB')] <- c(1,1,1,1)/4
  tMatrix['XXRB', 'XYGR', c('XXGR','XXGB','XXRR','XXRB',
                            'XYGR','XYGB','XYRR','XYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXRB', 'XYGB', c('XXGR','XXGB','XXRB','XXBB',
                            'XYGR','XYGB','XYRB','XYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXRB', 'XYRR', c('XXRR','XXRB',
                            'XYRR','XYRB')] <- c(1,1,1,1)/4
  tMatrix['XXRB', 'XYRB', c('XXRR','XXRB','XXBB',
                            'XYRR','XYRB','XYBB')] <- c(1,2,1,1,2,1)/8
  tMatrix['XXRB', 'XYBB', c('XXRB','XXBB',
                            'XYRB','XYBB')] <- c(1,1,1,1)/4

  tMatrix['XXRB', 'XCWW', c('XXWR','XXWB',
                            'XCWR','XCWB')] <- c(1,1,1,1)/4
  tMatrix['XXRB', 'XCWG', c('XXWR','XXGR','XXRR','XXRB','XXWB','XXGB','XXBB',
                            'XCWR','XCGR','XCRR','XCRB','XCWB','XCGB','XCBB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + cM*(1-chM)*crM, 1-cM, 1+cM*chM, cM*(1-chM)*(1-crM),
                                                                                    1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + cM*(1-chM)*crM, 1-cM, 1+cM*chM, cM*(1-chM)*(1-crM))/8
  tMatrix['XXRB', 'XCWR', c('XXWR','XXWB','XXRR','XXRB',
                            'XCWR','XCWB','XCRR','XCRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXRB', 'XCWB', c('XXWR','XXWB','XXRB','XXBB',
                            'XCWR','XCWB','XCRB','XCBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXRB', 'XCGG', c('XXGR','XXGB',
                            'XCGR','XCGB')] <- c(1,1,1,1)/4
  tMatrix['XXRB', 'XCGR', c('XXGR','XXGB','XXRR','XXRB',
                            'XCGR','XCGB','XCRR','XCRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXRB', 'XCGB', c('XXGR','XXGB','XXRB','XXBB',
                            'XCGR','XCGB','XCRB','XCBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXRB', 'XCRR', c('XXRR','XXRB',
                            'XCRR','XCRB')] <- c(1,1,1,1)/4
  tMatrix['XXRB', 'XCRB', c('XXRR','XXRB','XXBB',
                            'XCRR','XCRB','XCBB')] <- c(1,2,1,1,2,1)/8
  tMatrix['XXRB', 'XCBB', c('XXRB','XXBB',
                            'XCRB','XCBB')] <- c(1,1,1,1)/4


  tMatrix['XXBB', 'XYWW', c('XXWB','XYWB')] <- c(1,1)/2
  tMatrix['XXBB', 'XYWG', c('XXWB','XXGB',
                            'XYWB','XYGB')] <- c(1,1,1,1)/4
  tMatrix['XXBB', 'XYWR', c('XXWB','XXRB',
                            'XYWB','XYRB')] <- c(1,1,1,1)/4
  tMatrix['XXBB', 'XYWB', c('XXWB','XXBB',
                            'XYWB','XYBB')] <- c(1,1,1,1)/4
  tMatrix['XXBB', 'XYGG', c('XXGB','XYGB')] <- c(1,1)/2
  tMatrix['XXBB', 'XYGR', c('XXGB','XXRB',
                            'XYGB','XYRB')] <- c(1,1,1,1)/4
  tMatrix['XXBB', 'XYGB', c('XXGB','XXBB',
                            'XYGB','XYBB')] <- c(1,1,1,1)/4
  tMatrix['XXBB', 'XYRR', c('XXRB','XYRB')] <- c(1,1)/2
  tMatrix['XXBB', 'XYRB', c('XXRB','XXBB',
                            'XYRB','XYBB')] <- c(1,1,1,1)/4
  tMatrix['XXBB', 'XYBB', c('XXBB','XYBB')] <- c(1,1)/2

  tMatrix['XXBB', 'XCWW', c('XXWB','XCWB')] <- c(1,1)/2
  tMatrix['XXBB', 'XCWG', c('XXWB','XXGB','XXRB','XXBB',
                            'XCWB','XCGB','XCRB','XCBB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/4
  tMatrix['XXBB', 'XCWR', c('XXWB','XXRB',
                            'XCWB','XCRB')] <- c(1,1,1,1)/4
  tMatrix['XXBB', 'XCWB', c('XXWB','XXBB',
                            'XCWB','XCBB')] <- c(1,1,1,1)/4
  tMatrix['XXBB', 'XCGG', c('XXGB','XCGB')] <- c(1,1)/2
  tMatrix['XXBB', 'XCGR', c('XXGB','XXRB',
                            'XCGB','XCRB')] <- c(1,1,1,1)/4
  tMatrix['XXBB', 'XCGB', c('XXGB','XXBB',
                            'XCGB','XCBB')] <- c(1,1,1,1)/4
  tMatrix['XXBB', 'XCRR', c('XXRB','XCRB')] <- c(1,1)/2
  tMatrix['XXBB', 'XCRB', c('XXRB','XXBB',
                            'XCRB','XCBB')] <- c(1,1,1,1)/4
  tMatrix['XXBB', 'XCBB', c('XXBB','XCBB')] <- c(1,1)/2


  #protection from underflow errors
  tMatrix[tMatrix < .Machine$double.eps] <- 0


  ## initialize viability mask. No mother/father-specific death, so use basic mask
  viabilityMask <- array(data = 1, dim = c(size,size,size), dimnames = list(gtype, gtype, gtype))

  ## genotype-specific modifiers
  if(!is.null(phi)){
    stop("This cube has a special phi, due to being male/female specific.
         Please edit it after the cube is built, if you are absolutely sure it needs to change.")
  }
  phi = setNames(object = c(rep.int(x = 1, times = 10),rep.int(x = 0, times = 20)), nm = gtype)
  modifiers = cubeModifiers(gtype, eta = eta, phi = phi, omega = omega, xiF = xiF, xiM = xiM, s = s)

  ## put everything into a labeled list to return
  return(list(
    ih = tMatrix,
    tau = viabilityMask,
    genotypesID = gtype,
    genotypesN = size,
    wildType = c('XXWW','XYWW'),
    eta = modifiers$eta,
    phi = modifiers$phi,
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = "XCGG"
  ))

}
