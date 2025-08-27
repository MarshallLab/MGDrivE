###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   X-linked Split Drive with 2 Resistance Alleles and Sex-Specific homing
#   Héctor Sanchez, Jared Bennett, Sean Wu, John Marshall
#   October 2020
#   jared_bennett@berkeley.edu
#
###############################################################################

#' Inheritance Cube: X-linked Split CRISPR Drive with 2 Resistance Alleles and male/female specific homing
#'
#' This is a X-linked, sex-specific version of a split CRISPR drive. At the X locus
#' is the Cas9, inherited in a Mendelian fashion. At a second, unlinked, autosomal
#' locus are the gRNAs. When the two loci occur together, the gRNAs drive, with
#' potential damaged alleles, but the Cas9 remains Mendelian.
#' Deposition in this cube is performed when both pieces come together in females.
#' This drive has 2 loci:
#'  * "Locus" 1, sex chromosomes, has 3 alleles
#'    * X: Wild-type X chromosome
#'    * C: X-chromosome carrying a Cas9 construct
#'    * Y: Wild-type Y chromosome
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
cubeSplitDriveX <- function(cM = 1.0, chM = 0, crM = 0,
                            cF = 1.0, chF = 0, crF = 0, ccF = cF, cchF = chF, ccrF = crF,
                            dW = 0, dhW = 0, drW = 0, ddW = dW, ddhW = dhW, ddrW = drW,
                            eta = NULL, phi = NULL, omega = NULL, xiF = NULL, xiM = NULL, s = NULL){

  # # for testing
  # testVec <- runif(n = 15)
  #
  # cF <- testVec[1]; chF <- testVec[2]; crF <- testVec[3];
  # ccF <- testVec[4]; cchF <- testVec[5]; ccrF <- testVec[6];
  #
  # cM <- testVec[7]; chM <- testVec[8]; crM <- testVec[9];
  #
  # dW <- testVec[13]; dhW <- testVec[14]; drW <- testVec[15];
  # ddW <- testVec[10]; ddhW <- testVec[11]; ddrW <- testVec[12];
  #
  # # this would run at the bottom, to check that all cells sum to 1
  # #  size is defined below
  # for(mID in c('XYWW','XYWG','XYWR','XYWB','XYGG','XYGR','XYGB','XYRR','XYRB','XYBB',
  #              'CYWW','CYWG','CYWR','CYWB','CYGG','CYGR','CYGB','CYRR','CYRB','CYBB')){
  #   for(fID in c('XXWW','XXWG','XXWR','XXWB','XXGG','XXGR','XXGB','XXRR','XXRB','XXBB',
  #                'XCWW','XCWG','XCWR','XCWB','XCGG','XCGR','XCGB','XCRR','XCRB','XCBB',
  #                'CCWW','CCWG','CCWR','CCWB','CCGG','CCGR','CCGB','CCRR','CCRB','CCBB')){
  #     if((abs(sum(tMatrix[fID,mID,])-1)) > 1e-6){
  #       print(paste0("Fail: mID= ",mID,", fID= ",fID))
  #     }
  #   }
  # }


  ## safety checks
  params <- c(cM, chM, crM, cF, chF, crF, ccF, cchF, ccrF, dW, dhW, drW, ddW, ddhW, ddrW)
  if(any(params>1) || any(params<0)){
    stop("Parameters are rates, they must be between 0 and 1.")
  }


  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  gtype <- c('XXWW','XXWG','XXWR','XXWB','XXGG','XXGR','XXGB','XXRR','XXRB','XXBB',
             'XCWW','XCWG','XCWR','XCWB','XCGG','XCGR','XCGB','XCRR','XCRB','XCBB',
             'CCWW','CCWG','CCWR','CCWB','CCGG','CCGR','CCGB','CCRR','CCRB','CCBB',
             'XYWW','XYWG','XYWR','XYWB','XYGG','XYGR','XYGB','XYRR','XYRB','XYBB',
             'CYWW','CYWG','CYWR','CYWB','CYGG','CYGR','CYGB','CYRR','CYRB','CYBB')
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

  tMatrix['XXWW', 'CYWW', c('XCWW','XYWW')] <- c(1,1)/2
  tMatrix['XXWW', 'CYWG', c('XCWW','XCWG','XCWR','XCWB',
                            'XYWW','XYWG','XYWR','XYWB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/4
  tMatrix['XXWW', 'CYWR', c('XCWW','XCWR',
                            'XYWW','XYWR')] <- c(1,1,1,1)/4
  tMatrix['XXWW', 'CYWB', c('XCWW','XCWB',
                            'XYWW','XYWB')] <- c(1,1,1,1)/4
  tMatrix['XXWW', 'CYGG', c('XCWG','XYWG')] <- c(1,1)/2
  tMatrix['XXWW', 'CYGR', c('XCWG','XCWR',
                            'XYWG','XYWR')] <- c(1,1,1,1)/4
  tMatrix['XXWW', 'CYGB', c('XCWG','XCWB',
                            'XYWG','XYWB')] <- c(1,1,1,1)/4
  tMatrix['XXWW', 'CYRR', c('XCWR','XYWR')] <- c(1,1)/2
  tMatrix['XXWW', 'CYRB', c('XCWR','XCWB',
                            'XYWR','XYWB')] <- c(1,1,1,1)/4
  tMatrix['XXWW', 'CYBB', c('XCWB','XYWB')] <- c(1,1)/2


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

  tMatrix['XXWG', 'CYWW', c('XCWW','XCWG',
                            'XYWW','XYWG')] <- c(1,1,1,1)/4
  tMatrix['XXWG', 'CYWG', c('XCWW','XCWG','XCWR','XCWB','XCGG','XCGR','XCGB',
                            'XYWW','XYWG','XYWR','XYWB','XYGG','XYGR','XYGB')] <- c(1-cM, 1+cM*chM + 1-cM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),  1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                                                    1-cM, 1+cM*chM + 1-cM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),  1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['XXWG', 'CYWR', c('XCWW','XCWG','XCWR','XCGR',
                            'XYWW','XYWG','XYWR','XYGR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWG', 'CYWB', c('XCWW','XCWG','XCWB','XCGB',
                            'XYWW','XYWG','XYWB','XYGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWG', 'CYGG', c('XCWG','XCGG',
                            'XYWG','XYGG')] <- c(1,1,1,1)/4
  tMatrix['XXWG', 'CYGR', c('XCWG','XCGG','XCWR','XCGR',
                            'XYWG','XYGG','XYWR','XYGR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWG', 'CYGB', c('XCWG','XCGG','XCWB','XCGB',
                            'XYWG','XYGG','XYWB','XYGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWG', 'CYRR', c('XCWR','XCGR',
                            'XYWR','XYGR')] <- c(1,1,1,1)/4
  tMatrix['XXWG', 'CYRB', c('XCWR','XCGR','XCWB','XCGB',
                            'XYWR','XYGR','XYWB','XYGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWG', 'CYBB', c('XCWB','XCGB',
                            'XYWB','XYGB')] <- c(1,1,1,1)/4


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

  tMatrix['XXWR', 'CYWW', c('XCWW','XCWR',
                            'XYWW','XYWR')] <- c(1,1,1,1)/4
  tMatrix['XXWR', 'CYWG', c('XCWW','XCWG','XCWR','XCWB','XCGR','XCRR','XCRB',
                            'XYWW','XYWG','XYWR','XYWB','XYGR','XYRR','XYRB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM + 1-cM, cM*(1-chM)*(1-crM), 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                                                    1-cM, 1+cM*chM, cM*(1-chM)*crM + 1-cM, cM*(1-chM)*(1-crM), 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['XXWR', 'CYWR', c('XCWW','XCWR','XCRR',
                            'XYWW','XYWR','XYRR')] <- c(1,2,1,1,2,1)/8
  tMatrix['XXWR', 'CYWB', c('XCWW','XCWR','XCWB','XCRB',
                            'XYWW','XYWR','XYWB','XYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWR', 'CYGG', c('XCWG','XCGR',
                            'XYWG','XYGR')] <- c(1,1,1,1)/4
  tMatrix['XXWR', 'CYGR', c('XCWG','XCGR','XCWR','XCRR',
                            'XYWG','XYGR','XYWR','XYRR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWR', 'CYGB', c('XCWG','XCGR','XCWB','XCRB',
                            'XYWG','XYGR','XYWB','XYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWR', 'CYRR', c('XCWR','XCRR',
                            'XYWR','XYRR')] <- c(1,1,1,1)/4
  tMatrix['XXWR', 'CYRB', c('XCWR','XCRR','XCWB','XCRB',
                            'XYWR','XYRR','XYWB','XYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWR', 'CYBB', c('XCWB','XCRB',
                            'XYWB','XYRB')] <- c(1,1,1,1)/4


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

  tMatrix['XXWB', 'CYWW', c('XCWW','XCWB',
                            'XYWW','XYWB')] <- c(1,1,1,1)/4
  tMatrix['XXWB', 'CYWG', c('XCWW','XCWG','XCWR','XCWB','XCGB','XCRB','XCBB',
                            'XYWW','XYWG','XYWR','XYWB','XYGB','XYRB','XYBB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + 1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                                                    1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + 1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['XXWB', 'CYWR', c('XCWW','XCWB','XCWR','XCRB',
                            'XYWW','XYWB','XYWR','XYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWB', 'CYWB', c('XCWW','XCWB','XCBB',
                            'XYWW','XYWB','XYBB')] <- c(1,2,1,1,2,1)/8
  tMatrix['XXWB', 'CYGG', c('XCWG','XCGB',
                            'XYWG','XYGB')] <- c(1,1,1,1)/4
  tMatrix['XXWB', 'CYGR', c('XCWG','XCGB','XCWR','XCGR',
                            'XYWG','XYGB','XYWR','XYGR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWB', 'CYGB', c('XCWG','XCGB','XCWB','XCBB',
                            'XYWG','XYGB','XYWB','XYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWB', 'CYRR', c('XCWR','XCRB',
                            'XYWR','XYRB')] <- c(1,1,1,1)/4
  tMatrix['XXWB', 'CYRB', c('XCWR','XCRB','XCWB','XCBB',
                            'XYWR','XYRB','XYWB','XYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXWB', 'CYBB', c('XCWB','XCBB',
                            'XYWB','XYBB')] <- c(1,1,1,1)/4


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

  tMatrix['XXGG', 'CYWW', c('XCWG','XYWG')] <- c(1,1)/2
  tMatrix['XXGG', 'CYWG', c('XCWG','XCGG','XCGR','XCGB',
                            'XYWG','XYGG','XYGR','XYGB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/4
  tMatrix['XXGG', 'CYWR', c('XCWG','XCGR',
                            'XYWG','XYGR')] <- c(1,1,1,1)/4
  tMatrix['XXGG', 'CYWB', c('XCWG','XCGB',
                            'XYWG','XYGB')] <- c(1,1,1,1)/4
  tMatrix['XXGG', 'CYGG', c('XCGG','XYGG')] <- c(1,1)/2
  tMatrix['XXGG', 'CYGR', c('XCGG','XCGR',
                            'XYGG','XYGR')] <- c(1,1,1,1)/4
  tMatrix['XXGG', 'CYGB', c('XCGG','XCGB',
                            'XYGG','XYGB')] <- c(1,1,1,1)/4
  tMatrix['XXGG', 'CYRR', c('XCGR','XYGR')] <- c(1,1)/2
  tMatrix['XXGG', 'CYRB', c('XCGR','XCGB',
                            'XYGR','XYGB')] <- c(1,1,1,1)/4
  tMatrix['XXGG', 'CYBB', c('XCGB','XYGB')] <- c(1,1)/2


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

  tMatrix['XXGR', 'CYWW', c('XCWG','XCWR',
                            'XYWG','XYWR')] <- c(1,1,1,1)/4
  tMatrix['XXGR', 'CYWG', c('XCWG','XCGG','XCGR','XCGB','XCWR','XCRR','XCRB',
                            'XYWG','XYGG','XYGR','XYGB','XYWR','XYRR','XYRB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM + 1+cM*chM, cM*(1-chM)*(1-crM), 1-cM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                                                    1-cM, 1+cM*chM, cM*(1-chM)*crM + 1+cM*chM, cM*(1-chM)*(1-crM), 1-cM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['XXGR', 'CYWR', c('XCWG','XCWR','XCGR','XCRR',
                            'XYWG','XYWR','XYGR','XYRR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGR', 'CYWB', c('XCWG','XCWR','XCGB','XCRB',
                            'XYWG','XYWR','XYGB','XYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGR', 'CYGG', c('XCGG','XCGR',
                            'XYGG','XYGR')] <- c(1,1,1,1)/4
  tMatrix['XXGR', 'CYGR', c('XCGG','XCGR','XCRR',
                            'XYGG','XYGR','XYRR')] <- c(1,2,1,1,2,1)/8
  tMatrix['XXGR', 'CYGB', c('XCGG','XCGR','XCGB','XCRB',
                            'XYGG','XYGR','XYGB','XYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGR', 'CYRR', c('XCGR','XCRR',
                            'XYGR','XYRR')] <- c(1,1,1,1)/4
  tMatrix['XXGR', 'CYRB', c('XCGR','XCRR','XCGB','XCRB',
                            'XYGR','XYRR','XYGB','XYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGR', 'CYBB', c('XCGB','XCRB',
                            'XYGB','XYRB')] <- c(1,1,1,1)/4


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

  tMatrix['XXGB', 'CYWW', c('XCWG','XCWB',
                            'XYWG','XYWB')] <- c(1,1,1,1)/4
  tMatrix['XXGB', 'CYWG', c('XCWG','XCGG','XCGR','XCGB','XCWB','XCRB','XCBB',
                            'XYWG','XYGG','XYGR','XYGB','XYWB','XYRB','XYBB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + 1+cM*chM, 1-cM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                                                    1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + 1+cM*chM, 1-cM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['XXGB', 'CYWR', c('XCWG','XCWB','XCGR','XCRB',
                            'XYWG','XYWB','XYGR','XYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGB', 'CYWB', c('XCWG','XCWB','XCGB','XCBB',
                            'XYWG','XYWB','XYGB','XYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGB', 'CYGG', c('XCGG','XCGB',
                            'XYGG','XYGB')] <- c(1,1,1,1)/4
  tMatrix['XXGB', 'CYGR', c('XCGG','XCGB','XCGR','XCRB',
                            'XYGG','XYGB','XYGR','XYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGB', 'CYGB', c('XCGG','XCGB','XCBB',
                            'XYGG','XYGB','XYBB')] <- c(1,2,1,1,2,1)/8
  tMatrix['XXGB', 'CYRR', c('XCGR','XCRB',
                            'XYGR','XYRB')] <- c(1,1,1,1)/4
  tMatrix['XXGB', 'CYRB', c('XCGR','XCRB','XCGB','XCBB',
                            'XYGR','XYRB','XYGB','XYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXGB', 'CYBB', c('XCGB','XCBB',
                            'XYGB','XYBB')] <- c(1,1,1,1)/4


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

  tMatrix['XXRR', 'CYWW', c('XCWR','XYWR')] <- c(1,1)/2
  tMatrix['XXRR', 'CYWG', c('XCWR','XCGR','XCRR','XCRB',
                            'XYWR','XYGR','XYRR','XYRB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/4
  tMatrix['XXRR', 'CYWR', c('XCWR','XCRR',
                            'XYWR','XYRR')] <- c(1,1,1,1)/4
  tMatrix['XXRR', 'CYWB', c('XCWR','XCRB',
                            'XYWR','XYRB')] <- c(1,1,1,1)/4
  tMatrix['XXRR', 'CYGG', c('XCGR','XYGR')] <- c(1,1)/2
  tMatrix['XXRR', 'CYGR', c('XCGR','XCRR',
                            'XYGR','XYRR')] <- c(1,1,1,1)/4
  tMatrix['XXRR', 'CYGB', c('XCGR','XCRB',
                            'XYGR','XYRB')] <- c(1,1,1,1)/4
  tMatrix['XXRR', 'CYRR', c('XCRR','XYRR')] <- c(1,1)/2
  tMatrix['XXRR', 'CYRB', c('XCRR','XCRB',
                            'XYRR','XYRB')] <- c(1,1,1,1)/4
  tMatrix['XXRR', 'CYBB', c('XCRB','XYRB')] <- c(1,1)/2


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

  tMatrix['XXRB', 'CYWW', c('XCWR','XCWB',
                            'XYWR','XYWB')] <- c(1,1,1,1)/4
  tMatrix['XXRB', 'CYWG', c('XCWR','XCGR','XCRR','XCRB','XCWB','XCGB','XCBB',
                            'XYWR','XYGR','XYRR','XYRB','XYWB','XYGB','XYBB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + cM*(1-chM)*crM, 1-cM, 1+cM*chM, cM*(1-chM)*(1-crM),
                                                                                    1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + cM*(1-chM)*crM, 1-cM, 1+cM*chM, cM*(1-chM)*(1-crM))/8
  tMatrix['XXRB', 'CYWR', c('XCWR','XCWB','XCRR','XCRB',
                            'XYWR','XYWB','XYRR','XYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXRB', 'CYWB', c('XCWR','XCWB','XCRB','XCBB',
                            'XYWR','XYWB','XYRB','XYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXRB', 'CYGG', c('XCGR','XCGB',
                            'XYGR','XYGB')] <- c(1,1,1,1)/4
  tMatrix['XXRB', 'CYGR', c('XCGR','XCGB','XCRR','XCRB',
                            'XYGR','XYGB','XYRR','XYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXRB', 'CYGB', c('XCGR','XCGB','XCRB','XCBB',
                            'XYGR','XYGB','XYRB','XYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XXRB', 'CYRR', c('XCRR','XCRB',
                            'XYRR','XYRB')] <- c(1,1,1,1)/4
  tMatrix['XXRB', 'CYRB', c('XCRR','XCRB','XCBB',
                            'XYRR','XYRB','XYBB')] <- c(1,2,1,1,2,1)/8
  tMatrix['XXRB', 'CYBB', c('XCRB','XCBB',
                            'XYRB','XYBB')] <- c(1,1,1,1)/4


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

  tMatrix['XXBB', 'CYWW', c('XCWB','XYWB')] <- c(1,1)/2
  tMatrix['XXBB', 'CYWG', c('XCWB','XCGB','XCRB','XCBB',
                            'XYWB','XYGB','XYRB','XYBB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/4
  tMatrix['XXBB', 'CYWR', c('XCWB','XCRB',
                            'XYWB','XYRB')] <- c(1,1,1,1)/4
  tMatrix['XXBB', 'CYWB', c('XCWB','XCBB',
                            'XYWB','XYBB')] <- c(1,1,1,1)/4
  tMatrix['XXBB', 'CYGG', c('XCGB','XYGB')] <- c(1,1)/2
  tMatrix['XXBB', 'CYGR', c('XCGB','XCRB',
                            'XYGB','XYRB')] <- c(1,1,1,1)/4
  tMatrix['XXBB', 'CYGB', c('XCGB','XCBB',
                            'XYGB','XYBB')] <- c(1,1,1,1)/4
  tMatrix['XXBB', 'CYRR', c('XCRB','XYRB')] <- c(1,1)/2
  tMatrix['XXBB', 'CYRB', c('XCRB','XCBB',
                            'XYRB','XYBB')] <- c(1,1,1,1)/4
  tMatrix['XXBB', 'CYBB', c('XCBB','XYBB')] <- c(1,1)/2


  tMatrix['XCWW', 'XYWW', c('XXWW','XCWW',
                            'XYWW','CYWW')] <- c(1,1,1,1)/4
  tMatrix['XCWW', 'XYWG', c('XXWW','XXWG','XCWW','XCWG',
                            'XYWW','XYWG','CYWW','CYWG')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCWW', 'XYWR', c('XXWW','XXWR','XCWW','XCWR',
                            'XYWW','XYWR','CYWW','CYWR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCWW', 'XYWB', c('XXWW','XXWB','XCWW','XCWB',
                            'XYWW','XYWB','CYWW','CYWB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCWW', 'XYGG', c('XXWG','XCWG',
                            'XYWG','CYWG')] <- c(1,1,1,1)/4
  tMatrix['XCWW', 'XYGR', c('XXWG','XXWR','XCWG','XCWR',
                            'XYWG','XYWR','CYWG','CYWR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCWW', 'XYGB', c('XXWG','XXWB','XCWG','XCWB',
                            'XYWG','XYWB','CYWG','CYWB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCWW', 'XYRR', c('XXWR','XCWR',
                            'XYWR','CYWR')] <- c(1,1,1,1)/4
  tMatrix['XCWW', 'XYRB', c('XXWR','XXWB','XCWR','XCWB',
                            'XYWR','XYWB','CYWR','CYWB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCWW', 'XYBB', c('XXWB','XCWB',
                            'XYWB','CYWB')] <- c(1,1,1,1)/4

  tMatrix['XCWW', 'CYWW', c('XCWW','CCWW',
                            'XYWW','CYWW')] <- c(1,1,1,1)/4
  tMatrix['XCWW', 'CYWG', c('XCWW','XCWG','XCWR','XCWB',
                            'CCWW','CCWG','CCWR','CCWB',
                            'XYWW','XYWG','XYWR','XYWB',
                            'CYWW','CYWG','CYWR','CYWB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['XCWW', 'CYWR', c('XCWW','XCWR','CCWW','CCWR',
                            'XYWW','XYWR','CYWW','CYWR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCWW', 'CYWB', c('XCWW','XCWB','CCWW','CCWB',
                            'XYWW','XYWB','CYWW','CYWB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCWW', 'CYGG', c('XCWG','CCWG',
                            'XYWG','CYWG')] <- c(1,1,1,1)/4
  tMatrix['XCWW', 'CYGR', c('XCWG','XCWR','CCWG','CCWR',
                            'XYWG','XYWR','CYWG','CYWR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCWW', 'CYGB', c('XCWG','XCWB','CCWG','CCWB',
                            'XYWG','XYWB','CYWG','CYWB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCWW', 'CYRR', c('XCWR','CCWR',
                            'XYWR','CYWR')] <- c(1,1,1,1)/4
  tMatrix['XCWW', 'CYRB', c('XCWR','XCWB','CCWR','CCWB',
                            'XYWR','XYWB','CYWR','CYWB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCWW', 'CYBB', c('XCWB','CCWB',
                            'XYWB','CYWB')] <- c(1,1,1,1)/4


  eTen <- rep.int(x = 0, times = 10)
  wrbTen <- c((1-cF)*(1-dW), (1+cF*chF)*(1-dW), (cF*(1-chF)*crF)*(1-dW) + (1-cF)*dW*drW,
              (cF*(1-chF)*(1-crF))*(1-dW) + (1-cF)*dW*(1-drW), (1+cF*chF)*dW*dhW,
              (1+cF*chF)*dW*(1-dhW)*drW, (1+cF*chF)*dW*(1-dhW)*(1-drW),
              (cF*(1-chF)*crF)*dW*drW, (cF*(1-chF)*crF)*dW*(1-drW) + (cF*(1-chF)*(1-crF))*dW*drW,
              (cF*(1-chF)*(1-crF))*dW*(1-drW))/8
  tMatrix['XCWG', 'XYWW', ] <- c(wrbTen, wrbTen, eTen, wrbTen, wrbTen)
  tMatrix['XCWG', 'CYWW', ] <- c(eTen, wrbTen, wrbTen, wrbTen, wrbTen)

  wrbTen <- c((1-cF)*(1-dW), 1-cF + (1+cF*chF)*(1-dW), (cF*(1-chF)*crF)*(1-dW) + (1-cF)*dW*drW,
              (cF*(1-chF)*(1-crF))*(1-dW) + (1-cF)*dW*(1-drW), 1+cF*chF + (1+cF*chF)*dW*dhW,
              cF*(1-chF)*crF + (1+cF*chF)*dW*(1-dhW)*drW, cF*(1-chF)*(1-crF) + (1+cF*chF)*dW*(1-dhW)*(1-drW),
              (cF*(1-chF)*crF)*dW*drW, (cF*(1-chF)*crF)*dW*(1-drW) + (cF*(1-chF)*(1-crF))*dW*drW,
              (cF*(1-chF)*(1-crF))*dW*(1-drW))/16
  tMatrix['XCWG', 'XYWG', ] <- c(wrbTen, wrbTen, eTen, wrbTen, wrbTen)

  wrbTen <- c((1-cF)*(1-dW), (1+cF*chF)*(1-dW), 1-cF + (cF*(1-chF)*crF)*(1-dW) + (1-cF)*dW*drW,
              (cF*(1-chF)*(1-crF))*(1-dW) + (1-cF)*dW*(1-drW), (1+cF*chF)*dW*dhW,
              1+cF*chF + (1+cF*chF)*dW*(1-dhW)*drW, (1+cF*chF)*dW*(1-dhW)*(1-drW), cF*(1-chF)*crF + (cF*(1-chF)*crF)*dW*drW,
              cF*(1-chF)*(1-crF) + (cF*(1-chF)*crF)*dW*(1-drW) + (cF*(1-chF)*(1-crF))*dW*drW,
              (cF*(1-chF)*(1-crF))*dW*(1-drW))/16
  tMatrix['XCWG', 'XYWR', ] <- c(wrbTen, wrbTen, eTen, wrbTen, wrbTen)
  tMatrix['XCWG', 'CYWR', ] <- c(eTen, wrbTen, wrbTen, wrbTen, wrbTen)

  wrbTen <- c((1-cF)*(1-dW), (1+cF*chF)*(1-dW), (cF*(1-chF)*crF)*(1-dW) + (1-cF)*dW*drW,
              1-cF + (cF*(1-chF)*(1-crF))*(1-dW) + (1-cF)*dW*(1-drW),
              (1+cF*chF)*dW*dhW, (1+cF*chF)*dW*(1-dhW)*drW, 1+cF*chF + (1+cF*chF)*dW*(1-dhW)*(1-drW),
              (cF*(1-chF)*crF)*dW*drW, cF*(1-chF)*crF + (cF*(1-chF)*crF)*dW*(1-drW) + (cF*(1-chF)*(1-crF))*dW*drW,
              cF*(1-chF)*(1-crF) + (cF*(1-chF)*(1-crF))*dW*(1-drW))/16
  tMatrix['XCWG', 'XYWB', ] <- c(wrbTen, wrbTen, eTen, wrbTen, wrbTen)
  tMatrix['XCWG', 'CYWB', ] <- c(eTen, wrbTen, wrbTen, wrbTen, wrbTen)

  wrbTen <- c(0, 1-cF, 0, 0, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF), 0, 0, 0)/8
  tMatrix['XCWG', 'XYGG', ] <- c(wrbTen, wrbTen, eTen, wrbTen, wrbTen)
  tMatrix['XCWG', 'CYGG', ] <- c(eTen, wrbTen, wrbTen, wrbTen, wrbTen)

  wrbTen <- c(0, 1-cF, 1-cF, 0, 1+cF*chF, cF*(1-chF)*crF + 1+cF*chF,
              cF*(1-chF)*(1-crF),cF*(1-chF)*crF, cF*(1-chF)*(1-crF), 0)/16
  tMatrix['XCWG', 'XYGR', ] <- c(wrbTen, wrbTen, eTen, wrbTen, wrbTen)
  tMatrix['XCWG', 'CYGR', ] <- c(eTen, wrbTen, wrbTen, wrbTen, wrbTen)

  wrbTen <- c(0, 1-cF, 0, 1-cF, 1+cF*chF, cF*(1-chF)*crF,
              cF*(1-chF)*(1-crF) + 1+cF*chF, 0, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/16
  tMatrix['XCWG', 'XYGB', ] <- c(wrbTen, wrbTen, eTen, wrbTen, wrbTen)
  tMatrix['XCWG', 'CYGB', ] <- c(eTen, wrbTen, wrbTen, wrbTen, wrbTen)

  wrbTen <- c(0, 0, 1-cF, 0, 0, 1+cF*chF, 0, cF*(1-chF)*crF, cF*(1-chF)*(1-crF), 0)/8
  tMatrix['XCWG', 'XYRR', ] <- c(wrbTen, wrbTen, eTen, wrbTen, wrbTen)
  tMatrix['XCWG', 'CYRR', ] <- c(eTen, wrbTen, wrbTen, wrbTen, wrbTen)

  wrbTen <- c(0, 0, 1-cF, 1-cF, 0, 1+cF*chF, 1+cF*chF, cF*(1-chF)*crF,
              cF*(1-chF)*(1-crF) + cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/16
  tMatrix['XCWG', 'XYRB', ] <- c(wrbTen, wrbTen, eTen, wrbTen, wrbTen)
  tMatrix['XCWG', 'CYRB', ] <- c(eTen, wrbTen, wrbTen, wrbTen, wrbTen)

  wrbTen <- c(0, 0, 0, 1-cF, 0, 0, 1+cF*chF, 0, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/8
  tMatrix['XCWG', 'XYBB', ] <- c(wrbTen, wrbTen, eTen, wrbTen, wrbTen)
  tMatrix['XCWG', 'CYBB', ] <- c(eTen, wrbTen, wrbTen, wrbTen, wrbTen)

  wrbTen <- c((1-cF)*(1-cM)*(1-dW), (1+cF*chF)*(1-cM)*(1-dW) + (1-cF)*(1+cM*chM),
              (1-cF)*(1-cM)*dW*drW + (cF*(1-chF)*crF)*(1-cM)*(1-dW) + (1-cF)*(cM*(1-chM)*crM),
              (1-cF)*(1-cM)*dW*(1-drW) + (cF*(1-chF)*(1-crF))*(1-cM)*(1-dW) + (1-cF)*(cM*(1-chM)*(1-crM)),
              (1+cF*chF)*(1-cM)*dW*dhW + (1+cF*chF)*(1+cM*chM),
              (1+cF*chF)*(1-cM)*dW*(1-dhW)*drW + (cF*(1-chF)*crF)*(1+cM*chM) + (1+cF*chF)*(cM*(1-chM)*crM),
              (1+cF*chF)*(1-cM)*dW*(1-dhW)*(1-drW) + (cF*(1-chF)*(1-crF))*(1+cM*chM) + (1+cF*chF)*(cM*(1-chM)*(1-crM)),
              (cF*(1-chF)*crF)*(1-cM)*dW*drW + (cF*(1-chF)*crF)*(cM*(1-chM)*crM),
              (cF*(1-chF)*crF)*(1-cM)*dW*(1-drW) + (cF*(1-chF)*(1-crF))*(1-cM)*dW*drW + (cF*(1-chF)*(1-crF))*(cM*(1-chM)*crM) + (cF*(1-chF)*crF)*(cM*(1-chM)*(1-crM)),
              (cF*(1-chF)*(1-crF))*(1-cM)*dW*(1-drW) + (cF*(1-chF)*(1-crF))*(cM*(1-chM)*(1-crM)))/16
  tMatrix['XCWG', 'CYWG', ] <- c(eTen, wrbTen, wrbTen, wrbTen, wrbTen)


  tMatrix['XCWR', 'XYWW', c('XXWW','XXWR','XCWW','XCWR',
                            'XYWW','XYWR','CYWW','CYWR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCWR', 'XYWG', c('XXWW','XXWG','XXWR','XXGR',
                            'XCWW','XCWG','XCWR','XCGR',
                            'XYWW','XYWG','XYWR','XYGR',
                            'CYWW','CYWG','CYWR','CYGR')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCWR', 'XYWR', c('XXWW','XXWR','XXRR',
                            'XCWW','XCWR','XCRR',
                            'XYWW','XYWR','XYRR',
                            'CYWW','CYWR','CYRR')] <- c(1,2,1,1,2,1,1,2,1,1,2,1)/16
  tMatrix['XCWR', 'XYWB', c('XXWW','XXWB','XXWR','XXRB',
                            'XCWW','XCWB','XCWR','XCRB',
                            'XYWW','XYWB','XYWR','XYRB',
                            'CYWW','CYWB','CYWR','CYRB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCWR', 'XYGG', c('XXWG','XXGR','XCWG','XCGR',
                            'XYWG','XYGR','CYWG','CYGR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCWR', 'XYGR', c('XXWG','XXWR','XXGR','XXRR',
                            'XCWG','XCWR','XCGR','XCRR',
                            'XYWG','XYWR','XYGR','XYRR',
                            'CYWG','CYWR','CYGR','CYRR')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCWR', 'XYGB', c('XXWG','XXWB','XXGR','XXRB',
                            'XCWG','XCWB','XCGR','XCRB',
                            'XYWG','XYWB','XYGR','XYRB',
                            'CYWG','CYWB','CYGR','CYRB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCWR', 'XYRR', c('XXWR','XXRR','XCWR','XCRR',
                            'XYWR','XYRR','CYWR','CYRR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCWR', 'XYRB', c('XXWR','XXWB','XXRR','XXRB',
                            'XCWR','XCWB','XCRR','XCRB',
                            'XYWR','XYWB','XYRR','XYRB',
                            'CYWR','CYWB','CYRR','CYRB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCWR', 'XYBB', c('XXWB','XXRB','XCWB','XCRB',
                            'XYWB','XYRB','CYWB','CYRB')] <- c(1,1,1,1,1,1,1,1)/8

  tMatrix['XCWR', 'CYWW', c('XCWW','XCWR','CCWW','CCWR',
                            'XYWW','XYWR','CYWW','CYWR')] <- c(1,1,1,1,1,1,1,1)/8
  wrb7 <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM + 1-cM, cM*(1-chM)*(1-crM), 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/16
  tMatrix['XCWR', 'CYWG', c('XCWW','XCWG','XCWR','XCWB','XCGR','XCRR','XCRB',
                            'CCWW','CCWG','CCWR','CCWB','CCGR','CCRR','CCRB',
                            'XYWW','XYWG','XYWR','XYWB','XYGR','XYRR','XYRB',
                            'CYWW','CYWG','CYWR','CYWB','CYGR','CYRR','CYRB')] <- c(wrb7,wrb7,wrb7,wrb7)
  tMatrix['XCWR', 'CYWR', c('XCWW','XCWR','XCRR',
                            'CCWW','CCWR','CCRR',
                            'XYWW','XYWR','XYRR',
                            'CYWW','CYWR','CYRR')] <- c(1,2,1,1,2,1,1,2,1,1,2,1)/16
  tMatrix['XCWR', 'CYWB', c('XCWW','XCWB','XCWR','XCRB',
                            'CCWW','CCWB','CCWR','CCRB',
                            'XYWW','XYWB','XYWR','XYRB',
                            'CYWW','CYWB','CYWR','CYRB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCWR', 'CYGG', c('XCWG','XCGR','CCWG','CCGR',
                            'XYWG','XYGR','CYWG','CYGR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCWR', 'CYGR', c('XCWG','XCWR','XCGR','XCRR',
                            'CCWG','CCWR','CCGR','CCRR',
                            'XYWG','XYWR','XYGR','XYRR',
                            'CYWG','CYWR','CYGR','CYRR')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCWR', 'CYGB', c('XCWG','XCWB','XCGR','XCRB',
                            'CCWG','CCWB','CCGR','CCRB',
                            'XYWG','XYWB','XYGR','XYRB',
                            'CYWG','CYWB','CYGR','CYRB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCWR', 'CYRR', c('XCWR','XCRR','CCWR','CCRR',
                            'XYWR','XYRR','CYWR','CYRR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCWR', 'CYRB', c('XCWR','XCWB','XCRR','XCRB',
                            'CCWR','CCWB','CCRR','CCRB',
                            'XYWR','XYWB','XYRR','XYRB',
                            'CYWR','CYWB','CYRR','CYRB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCWR', 'CYBB', c('XCWB','XCRB','CCWB','CCRB',
                            'XYWB','XYRB','CYWB','CYRB')] <- c(1,1,1,1,1,1,1,1)/8


  tMatrix['XCWB', 'XYWW', c('XXWW','XXWB','XCWW','XCWB',
                            'XYWW','XYWB','CYWW','CYWB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCWB', 'XYWG', c('XXWW','XXWG','XXWB','XXGB',
                            'XCWW','XCWG','XCWB','XCGB',
                            'XYWW','XYWG','XYWB','XYGB',
                            'CYWW','CYWG','CYWB','CYGB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCWB', 'XYWR', c('XXWW','XXWB','XXWR','XXRB',
                            'XCWW','XCWB','XCWR','XCRB',
                            'XYWW','XYWB','XYWR','XYRB',
                            'CYWW','CYWB','CYWR','CYRB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCWB', 'XYWB', c('XXWW','XXWB','XXBB',
                            'XCWW','XCWB','XCBB',
                            'XYWW','XYWB','XYBB',
                            'CYWW','CYWB','CYBB')] <- c(1,2,1,1,2,1,1,2,1,1,2,1)/16
  tMatrix['XCWB', 'XYGG', c('XXWG','XXGB','XCWG','XCGB',
                            'XYWG','XYGB','CYWG','CYGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCWB', 'XYGR', c('XXWG','XXWR','XXGB','XXRB',
                            'XCWG','XCWR','XCGB','XCRB',
                            'XYWG','XYWR','XYGB','XYRB',
                            'CYWG','CYWR','CYGB','CYRB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCWB', 'XYGB', c('XXWG','XXWB','XXGB','XXBB',
                            'XCWG','XCWB','XCGB','XCBB',
                            'XYWG','XYWB','XYGB','XYBB',
                            'CYWG','CYWB','CYGB','CYBB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCWB', 'XYRR', c('XXWR','XXRB','XCWR','XCRB',
                            'XYWR','XYRB','CYWR','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCWB', 'XYRB', c('XXWR','XXWB','XXRB','XXBB',
                            'XCWR','XCWB','XCRB','XCBB',
                            'XYWR','XYWB','XYRB','XYBB',
                            'CYWR','CYWB','CYRB','CYBB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCWB', 'XYBB', c('XXWB','XXBB','XCWB','XCBB',
                            'XYWB','XYBB','CYWB','CYBB')] <- c(1,1,1,1,1,1,1,1)/8

  tMatrix['XCWB', 'CYWW', c('XCWW','XCWB','CCWW','CCWB',
                            'XYWW','XYWB','CYWW','CYWB')] <- c(1,1,1,1,1,1,1,1)/8
  wrb7 <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + 1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/16
  tMatrix['XCWB', 'CYWG', c('XCWW','XCWG','XCWR','XCWB','XCGB','XCRB','XCBB',
                            'CCWW','CCWG','CCWR','CCWB','CCGB','CCRB','CCBB',
                            'XYWW','XYWG','XYWR','XYWB','XYGB','XYRB','XYBB',
                            'CYWW','CYWG','CYWR','CYWB','CYGB','CYRB','CYBB')] <- c(wrb7,wrb7,wrb7,wrb7)
  tMatrix['XCWB', 'CYWR', c('XCWW','XCWB','XCWR','XCRB',
                            'CCWW','CCWB','CCWR','CCRB',
                            'XYWW','XYWB','XYWR','XYRB',
                            'CYWW','CYWB','CYWR','CYRB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCWB', 'CYWB', c('XCWW','XCWB','XCBB',
                            'CCWW','CCWB','CCBB',
                            'XYWW','XYWB','XYBB',
                            'CYWW','CYWB','CYBB')] <- c(1,2,1,1,2,1,1,2,1,1,2,1)/16
  tMatrix['XCWB', 'CYGG', c('XCWG','XCGB','CCWG','CCGB',
                            'XYWG','XYGB','CYWG','CYGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCWB', 'CYGR', c('XCWG','XCWR','XCGB','XCRB',
                            'CCWG','CCWR','CCGB','CCRB',
                            'XYWG','XYWR','XYGB','XYRB',
                            'CYWG','CYWR','CYGB','CYRB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCWB', 'CYGB', c('XCWG','XCWB','XCGB','XCBB',
                            'CCWG','CCWB','CCGB','CCBB',
                            'XYWG','XYWB','XYGB','XYBB',
                            'CYWG','CYWB','CYGB','CYBB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCWB', 'CYRR', c('XCWR','XCRB','CCWR','CCRB',
                            'XYWR','XYRB','CYWR','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCWB', 'CYRB', c('XCWR','XCWB','XCRB','XCBB',
                            'CCWR','CCWB','CCRB','CCBB',
                            'XYWR','XYWB','XYRB','XYBB',
                            'CYWR','CYWB','CYRB','CYBB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCWB', 'CYBB', c('XCWB','XCBB','CCWB','CCBB',
                            'XYWB','XYBB','CYWB','CYBB')] <- c(1,1,1,1,1,1,1,1)/8


  tMatrix['XCGG', 'XYWW', c('XXWG','XXGG','XXGR','XXGB',
                            'XCWG','XCGG','XCGR','XCGB',
                            'XYWG','XYGG','XYGR','XYGB',
                            'CYWG','CYGG','CYGR','CYGB')] <- c(1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW),
                                                               1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW),
                                                               1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW),
                                                               1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW))/4
  tMatrix['XCGG', 'XYWG', c('XXWG','XXGG','XXGR','XXGB',
                            'XCWG','XCGG','XCGR','XCGB',
                            'XYWG','XYGG','XYGR','XYGB',
                            'CYWG','CYGG','CYGR','CYGB')] <- c(1-dW, 1 + dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW),
                                                               1-dW, 1 + dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW),
                                                               1-dW, 1 + dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW),
                                                               1-dW, 1 + dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW))/8
  tMatrix['XCGG', 'XYWR', c('XXWG','XXGG','XXGR','XXGB',
                            'XCWG','XCGG','XCGR','XCGB',
                            'XYWG','XYGG','XYGR','XYGB',
                            'CYWG','CYGG','CYGR','CYGB')] <- c(1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW),
                                                               1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW),
                                                               1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW),
                                                               1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW))/8
  tMatrix['XCGG', 'XYWB', c('XXWG','XXGG','XXGR','XXGB',
                            'XCWG','XCGG','XCGR','XCGB',
                            'XYWG','XYGG','XYGR','XYGB',
                            'CYWG','CYGG','CYGR','CYGB')] <- c(1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW),
                                                               1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW),
                                                               1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW),
                                                               1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW))/8
  tMatrix['XCGG', 'XYGG', c('XXGG','XCGG',
                            'XYGG','CYGG')] <- c(1,1,1,1)/4
  tMatrix['XCGG', 'XYGR', c('XXGG','XXGR','XCGG','XCGR',
                            'XYGG','XYGR','CYGG','CYGR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCGG', 'XYGB', c('XXGG','XXGB','XCGG','XCGB',
                            'XYGG','XYGB','CYGG','CYGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCGG', 'XYRR', c('XXGR','XCGR',
                            'XYGR','CYGR')] <- c(1,1,1,1)/4
  tMatrix['XCGG', 'XYRB', c('XXGR','XXGB','XCGR','XCGB',
                            'XYGR','XYGB','CYGR','CYGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCGG', 'XYBB', c('XXGB','XCGB',
                            'XYGB','CYGB')] <- c(1,1,1,1)/4

  tMatrix['XCGG', 'CYWW', c('XCWG','XCGG','XCGR','XCGB',
                            'CCWG','CCGG','CCGR','CCGB',
                            'XYWG','XYGG','XYGR','XYGB',
                            'CYWG','CYGG','CYGR','CYGB')] <- c(1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW),
                                                               1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW),
                                                               1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW),
                                                               1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW))/4
  tMatrix['XCGG', 'CYWG', c('XCWG','XCGG','XCGR','XCGB',
                            'CCWG','CCGG','CCGR','CCGB',
                            'XYWG','XYGG','XYGR','XYGB',
                            'CYWG','CYGG','CYGR','CYGB')] <- c((1-cM)*(1-dW), 1+cM*chM + (1-cM)*dW*dhW, cM*(1-chM)*crM + (1-cM)*dW*(1-dhW)*drW, cM*(1-chM)*(1-crM) + (1-cM)*dW*(1-dhW)*(1-drW),
                                                               (1-cM)*(1-dW), 1+cM*chM + (1-cM)*dW*dhW, cM*(1-chM)*crM + (1-cM)*dW*(1-dhW)*drW, cM*(1-chM)*(1-crM) + (1-cM)*dW*(1-dhW)*(1-drW),
                                                               (1-cM)*(1-dW), 1+cM*chM + (1-cM)*dW*dhW, cM*(1-chM)*crM + (1-cM)*dW*(1-dhW)*drW, cM*(1-chM)*(1-crM) + (1-cM)*dW*(1-dhW)*(1-drW),
                                                               (1-cM)*(1-dW), 1+cM*chM + (1-cM)*dW*dhW, cM*(1-chM)*crM + (1-cM)*dW*(1-dhW)*drW, cM*(1-chM)*(1-crM) + (1-cM)*dW*(1-dhW)*(1-drW))/8
  tMatrix['XCGG', 'CYWR', c('XCWG','XCGG','XCGR','XCGB',
                            'CCWG','CCGG','CCGR','CCGB',
                            'XYWG','XYGG','XYGR','XYGB',
                            'CYWG','CYGG','CYGR','CYGB')] <- c(1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW),
                                                               1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW),
                                                               1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW),
                                                               1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW))/8
  tMatrix['XCGG', 'CYWB', c('XCWG','XCGG','XCGR','XCGB',
                            'CCWG','CCGG','CCGR','CCGB',
                            'XYWG','XYGG','XYGR','XYGB',
                            'CYWG','CYGG','CYGR','CYGB')] <- c(1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW),
                                                               1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW),
                                                               1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW),
                                                               1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW))/8
  tMatrix['XCGG', 'CYGG', c('XCGG','CCGG',
                            'XYGG','CYGG')] <- c(1,1,1,1)/4
  tMatrix['XCGG', 'CYGR', c('XCGG','XCGR','CCGG','CCGR',
                            'XYGG','XYGR','CYGG','CYGR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCGG', 'CYGB', c('XCGG','XCGB','CCGG','CCGB',
                            'XYGG','XYGB','CYGG','CYGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCGG', 'CYRR', c('XCGR','CCGR',
                            'XYGR','CYGR')] <- c(1,1,1,1)/4
  tMatrix['XCGG', 'CYRB', c('XCGR','XCGB','CCGR','CCGB',
                            'XYGR','XYGB','CYGR','CYGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCGG', 'CYBB', c('XCGB','CCGB',
                            'XYGB','CYGB')] <- c(1,1,1,1)/4


  tMatrix['XCGR', 'XYWW', c('XXWG','XXGG','XXGR','XXGB','XXWR','XXRR','XXRB',
                            'XCWG','XCGG','XCGR','XCGB','XCWR','XCRR','XCRB',
                            'XYWG','XYGG','XYGR','XYGB','XYWR','XYRR','XYRB',
                            'CYWG','CYGG','CYGR','CYGB','CYWR','CYRR','CYRB')] <- c(1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW))/8
  tMatrix['XCGR', 'XYWG', c('XXWG','XXGG','XXGR','XXGB','XXWR','XXRR','XXRB',
                            'XCWG','XCGG','XCGR','XCGB','XCWR','XCRR','XCRB',
                            'XYWG','XYGG','XYGR','XYGB','XYWR','XYRR','XYRB',
                            'CYWG','CYGG','CYGR','CYGB','CYWR','CYRR','CYRB')] <- c(1-dW, 1 + dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW),
                                                                                    1-dW, 1 + dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW),
                                                                                    1-dW, 1 + dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW),
                                                                                    1-dW, 1 + dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW))/16
  tMatrix['XCGR', 'XYWR', c('XXWG','XXGG','XXGR','XXGB','XXWR','XXRR','XXRB',
                            'XCWG','XCGG','XCGR','XCGB','XCWR','XCRR','XCRB',
                            'XYWG','XYGG','XYGR','XYGB','XYWR','XYRR','XYRB',
                            'CYWG','CYGG','CYGR','CYGB','CYWR','CYRR','CYRB')] <- c(1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, 1 + dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, 1 + dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, 1 + dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, 1 + dW*drW, dW*(1-drW))/16
  tMatrix['XCGR', 'XYWB', c('XXWG','XXGG','XXGR','XXGB','XXWR','XXRR','XXRB',
                            'XCWG','XCGG','XCGR','XCGB','XCWR','XCRR','XCRB',
                            'XYWG','XYGG','XYGR','XYGB','XYWR','XYRR','XYRB',
                            'CYWG','CYGG','CYGR','CYGB','CYWR','CYRR','CYRB')] <- c(1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW), 1-dW, dW*drW, 1 + dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW), 1-dW, dW*drW, 1 + dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW), 1-dW, dW*drW, 1 + dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW), 1-dW, dW*drW, 1 + dW*(1-drW))/16
  tMatrix['XCGR', 'XYGG', c('XXGG','XXGR','XCGG','XCGR',
                            'XYGG','XYGR','CYGG','CYGR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCGR', 'XYGR', c('XXGG','XXGR','XXRR',
                            'XCGG','XCGR','XCRR',
                            'XYGG','XYGR','XYRR',
                            'CYGG','CYGR','CYRR')] <- c(1,2,1,1,2,1,1,2,1,1,2,1)/16
  tMatrix['XCGR', 'XYGB', c('XXGG','XXGR','XXGB','XXRB',
                            'XCGG','XCGR','XCGB','XCRB',
                            'XYGG','XYGR','XYGB','XYRB',
                            'CYGG','CYGR','CYGB','CYRB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCGR', 'XYRR', c('XXGR','XXRR','XCGR','XCRR',
                            'XYGR','XYRR','CYGR','CYRR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCGR', 'XYRB', c('XXGR','XXGB','XXRR','XXRB',
                            'XCGR','XCGB','XCRR','XCRB',
                            'XYGR','XYGB','XYRR','XYRB',
                            'CYGR','CYGB','CYRR','CYRB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCGR', 'XYBB', c('XXGB','XXRB','XCGB','XCRB',
                            'XYGB','XYRB','CYGB','CYRB')] <- c(1,1,1,1,1,1,1,1)/8

  tMatrix['XCGR', 'CYWW', c('XCWG','XCGG','XCGR','XCGB','XCWR','XCRR','XCRB',
                            'CCWG','CCGG','CCGR','CCGB','CCWR','CCRR','CCRB',
                            'XYWG','XYGG','XYGR','XYGB','XYWR','XYRR','XYRB',
                            'CYWG','CYGG','CYGR','CYGB','CYWR','CYRR','CYRB')] <- c(1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW))/8
  wrb7 <- c((1-cM)*(1-dW), 1+cM*chM + (1-cM)*dW*dhW, cM*(1-chM)*crM + 1+cM*chM + (1-cM)*dW*(1-dhW)*drW, cM*(1-chM)*(1-crM) + (1-cM)*dW*(1-dhW)*(1-drW), (1-cM)*(1-dW), cM*(1-chM)*crM + (1-cM)*dW*drW, cM*(1-chM)*(1-crM) + (1-cM)*dW*(1-drW))/16
  tMatrix['XCGR', 'CYWG', c('XCWG','XCGG','XCGR','XCGB','XCWR','XCRR','XCRB',
                            'CCWG','CCGG','CCGR','CCGB','CCWR','CCRR','CCRB',
                            'XYWG','XYGG','XYGR','XYGB','XYWR','XYRR','XYRB',
                            'CYWG','CYGG','CYGR','CYGB','CYWR','CYRR','CYRB')] <- c(wrb7,wrb7,wrb7,wrb7)
  tMatrix['XCGR', 'CYWR', c('XCWG','XCGG','XCGR','XCGB','XCWR','XCRR','XCRB',
                            'CCWG','CCGG','CCGR','CCGB','CCWR','CCRR','CCRB',
                            'XYWG','XYGG','XYGR','XYGB','XYWR','XYRR','XYRB',
                            'CYWG','CYGG','CYGR','CYGB','CYWR','CYRR','CYRB')] <- c(1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, 1 + dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, 1 + dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, 1 + dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, 1 + dW*drW, dW*(1-drW))/16
  tMatrix['XCGR', 'CYWB', c('XCWG','XCGG','XCGR','XCGB','XCWR','XCRR','XCRB',
                            'CCWG','CCGG','CCGR','CCGB','CCWR','CCRR','CCRB',
                            'XYWG','XYGG','XYGR','XYGB','XYWR','XYRR','XYRB',
                            'CYWG','CYGG','CYGR','CYGB','CYWR','CYRR','CYRB')] <- c(1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW), 1-dW, dW*drW, 1 + dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW), 1-dW, dW*drW, 1 + dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW), 1-dW, dW*drW, 1 + dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW), 1-dW, dW*drW, 1 + dW*(1-drW))/16
  tMatrix['XCGR', 'CYGG', c('XCGG','XCGR','CCGG','CCGR',
                            'XYGG','XYGR','CYGG','CYGR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCGR', 'CYGR', c('XCGG','XCGR','XCRR',
                            'CCGG','CCGR','CCRR',
                            'XYGG','XYGR','XYRR',
                            'CYGG','CYGR','CYRR')] <- c(1,2,1,1,2,1,1,2,1,1,2,1)/16
  tMatrix['XCGR', 'CYGB', c('XCGG','XCGR','XCGB','XCRB',
                            'CCGG','CCGR','CCGB','CCRB',
                            'XYGG','XYGR','XYGB','XYRB',
                            'CYGG','CYGR','CYGB','CYRB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCGR', 'CYRR', c('XCGR','XCRR','CCGR','CCRR',
                            'XYGR','XYRR','CYGR','CYRR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCGR', 'CYRB', c('XCGR','XCGB','XCRR','XCRB',
                            'CCGR','CCGB','CCRR','CCRB',
                            'XYGR','XYGB','XYRR','XYRB',
                            'CYGR','CYGB','CYRR','CYRB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCGR', 'CYBB', c('XCGB','XCRB','CCGB','CCRB',
                            'XYGB','XYRB','CYGB','CYRB')] <- c(1,1,1,1,1,1,1,1)/8

  tMatrix['XCGB', 'XYWW', c('XXWG','XXGG','XXGR','XXGB','XXWB','XXRB','XXBB',
                            'XCWG','XCGG','XCGR','XCGB','XCWB','XCRB','XCBB',
                            'XYWG','XYGG','XYGR','XYGB','XYWB','XYRB','XYBB',
                            'CYWG','CYGG','CYGR','CYGB','CYWB','CYRB','CYBB')] <- c(1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW))/8
  tMatrix['XCGB', 'XYWG', c('XXWG','XXGG','XXGR','XXGB','XXWB','XXRB','XXBB',
                            'XCWG','XCGG','XCGR','XCGB','XCWB','XCRB','XCBB',
                            'XYWG','XYGG','XYGR','XYGB','XYWB','XYRB','XYBB',
                            'CYWG','CYGG','CYGR','CYGB','CYWB','CYRB','CYBB')] <- c(1-dW, 1 + dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW),
                                                                                    1-dW, 1 + dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW),
                                                                                    1-dW, 1 + dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW),
                                                                                    1-dW, 1 + dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW))/16
  tMatrix['XCGB', 'XYWR', c('XXWG','XXGG','XXGR','XXGB','XXWB','XXRB','XXBB',
                            'XCWG','XCGG','XCGR','XCGB','XCWB','XCRB','XCBB',
                            'XYWG','XYGG','XYGR','XYGB','XYWB','XYRB','XYBB',
                            'CYWG','CYGG','CYGR','CYGB','CYWB','CYRB','CYBB')] <- c(1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, 1 + dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, 1 + dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, 1 + dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, 1 + dW*drW, dW*(1-drW))/16
  tMatrix['XCGB', 'XYWB', c('XXWG','XXGG','XXGR','XXGB','XXWB','XXRB','XXBB',
                            'XCWG','XCGG','XCGR','XCGB','XCWB','XCRB','XCBB',
                            'XYWG','XYGG','XYGR','XYGB','XYWB','XYRB','XYBB',
                            'CYWG','CYGG','CYGR','CYGB','CYWB','CYRB','CYBB')] <- c(1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW), 1-dW, dW*drW, 1 + dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW), 1-dW, dW*drW, 1 + dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW), 1-dW, dW*drW, 1 + dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW), 1-dW, dW*drW, 1 + dW*(1-drW))/16
  tMatrix['XCGB', 'XYGG', c('XXGG','XXGB','XCGG','XCGB',
                            'XYGG','XYGB','CYGG','CYGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCGB', 'XYGR', c('XXGG','XXGR','XXGB','XXRB',
                            'XCGG','XCGR','XCGB','XCRB',
                            'XYGG','XYGR','XYGB','XYRB',
                            'CYGG','CYGR','CYGB','CYRB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCGB', 'XYGB', c('XXGG','XXGB','XXBB',
                            'XCGG','XCGB','XCBB',
                            'XYGG','XYGB','XYBB',
                            'CYGG','CYGB','CYBB')] <- c(1,2,1,1,2,1,1,2,1,1,2,1)/16
  tMatrix['XCGB', 'XYRR', c('XXGR','XXRB','XCGR','XCRB',
                            'XYGR','XYRB','CYGR','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCGB', 'XYRB', c('XXGR','XXGB','XXRB','XXBB',
                            'XCGR','XCGB','XCRB','XCBB',
                            'XYGR','XYGB','XYRB','XYBB',
                            'CYGR','CYGB','CYRB','CYBB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCGB', 'XYBB', c('XXGB','XXBB','XCGB','XCBB',
                            'XYGB','XYBB','CYGB','CYBB')] <- c(1,1,1,1,1,1,1,1)/8

  tMatrix['XCGB', 'CYWW', c('XCWG','XCGG','XCGR','XCGB','XCWB','XCRB','XCBB',
                            'CCWG','CCGG','CCGR','CCGB','CCWB','CCRB','CCBB',
                            'XYWG','XYGG','XYGR','XYGB','XYWB','XYRB','XYBB',
                            'CYWG','CYGG','CYGR','CYGB','CYWB','CYRB','CYBB')] <- c(1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, dW*drW, dW*(1-drW))/8
  wrb7 <- c((1-cM)*(1-dW), 1+cM*chM + (1-cM)*dW*dhW, cM*(1-chM)*crM + (1-cM)*dW*(1-dhW)*drW, cM*(1-chM)*(1-crM) + 1+cM*chM + (1-cM)*dW*(1-dhW)*(1-drW), (1-cM)*(1-dW), cM*(1-chM)*crM + (1-cM)*dW*drW, cM*(1-chM)*(1-crM) + (1-cM)*dW*(1-drW))/16
  tMatrix['XCGB', 'CYWG', c('XCWG','XCGG','XCGR','XCGB','XCWB','XCRB','XCBB',
                            'CCWG','CCGG','CCGR','CCGB','CCWB','CCRB','CCBB',
                            'XYWG','XYGG','XYGR','XYGB','XYWB','XYRB','XYBB',
                            'CYWG','CYGG','CYGR','CYGB','CYWB','CYRB','CYBB')] <- c(wrb7,wrb7,wrb7,wrb7)
  tMatrix['XCGB', 'CYWR', c('XCWG','XCGG','XCGR','XCGB','XCWB','XCRB','XCBB',
                            'CCWG','CCGG','CCGR','CCGB','CCWB','CCRB','CCBB',
                            'XYWG','XYGG','XYGR','XYGB','XYWB','XYRB','XYBB',
                            'CYWG','CYGG','CYGR','CYGB','CYWB','CYRB','CYBB')] <- c(1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, 1 + dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, 1 + dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, 1 + dW*drW, dW*(1-drW),
                                                                                    1-dW, dW*dhW, 1 + dW*(1-dhW)*drW, dW*(1-dhW)*(1-drW), 1-dW, 1 + dW*drW, dW*(1-drW))/16
  tMatrix['XCGB', 'CYWB', c('XCWG','XCGG','XCGR','XCGB','XCWB','XCRB','XCBB',
                            'CCWG','CCGG','CCGR','CCGB','CCWB','CCRB','CCBB',
                            'XYWG','XYGG','XYGR','XYGB','XYWB','XYRB','XYBB',
                            'CYWG','CYGG','CYGR','CYGB','CYWB','CYRB','CYBB')] <- c(1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW), 1-dW, dW*drW, 1 + dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW), 1-dW, dW*drW, 1 + dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW), 1-dW, dW*drW, 1 + dW*(1-drW),
                                                                                    1-dW, dW*dhW, dW*(1-dhW)*drW, 1 + dW*(1-dhW)*(1-drW), 1-dW, dW*drW, 1 + dW*(1-drW))/16
  tMatrix['XCGB', 'CYGG', c('XCGG','XCGB','CCGG','CCGB',
                            'XYGG','XYGB','CYGG','CYGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCGB', 'CYGR', c('XCGG','XCGR','XCGB','XCRB',
                            'CCGG','CCGR','CCGB','CCRB',
                            'XYGG','XYGR','XYGB','XYRB',
                            'CYGG','CYGR','CYGB','CYRB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCGB', 'CYGB', c('XCGG','XCGB','XCBB',
                            'CCGG','CCGB','CCBB',
                            'XYGG','XYGB','XYBB',
                            'CYGG','CYGB','CYBB')] <- c(1,2,1,1,2,1,1,2,1,1,2,1)/16
  tMatrix['XCGB', 'CYRR', c('XCGR','XCRB','CCGR','CCRB',
                            'XYGR','XYRB','CYGR','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCGB', 'CYRB', c('XCGR','XCGB','XCRB','XCBB',
                            'CCGR','CCGB','CCRB','CCBB',
                            'XYGR','XYGB','XYRB','XYBB',
                            'CYGR','CYGB','CYRB','CYBB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCGB', 'CYBB', c('XCGB','XCBB','CCGB','CCBB',
                            'XYGB','XYBB','CYGB','CYBB')] <- c(1,1,1,1,1,1,1,1)/8


  tMatrix['XCRR', 'XYWW', c('XXWR','XCWR',
                            'XYWR','CYWR')] <- c(1,1,1,1)/4
  tMatrix['XCRR', 'XYWG', c('XXWR','XXGR','XCWR','XCGR',
                            'XYWR','XYGR','CYWR','CYGR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCRR', 'XYWR', c('XXWR','XXRR','XCWR','XCRR',
                            'XYWR','XYRR','CYWR','CYRR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCRR', 'XYWB', c('XXWR','XXRB','XCWR','XCRB',
                            'XYWR','XYRB','CYWR','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCRR', 'XYGG', c('XXGR','XCGR',
                            'XYGR','CYGR')] <- c(1,1,1,1)/4
  tMatrix['XCRR', 'XYGR', c('XXGR','XXRR','XCGR','XCRR',
                            'XYGR','XYRR','CYGR','CYRR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCRR', 'XYGB', c('XXGR','XXRB','XCGR','XCRB',
                            'XYGR','XYRB','CYGR','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCRR', 'XYRR', c('XXRR','XCRR',
                            'XYRR','CYRR')] <- c(1,1,1,1)/4
  tMatrix['XCRR', 'XYRB', c('XXRR','XXRB','XCRR','XCRB',
                            'XYRR','XYRB','CYRR','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCRR', 'XYBB', c('XXRB','XCRB',
                            'XYRB','CYRB')] <- c(1,1,1,1)/4

  tMatrix['XCRR', 'CYWW', c('XCWR','CCWR',
                            'XYWR','CYWR')] <- c(1,1,1,1)/4
  tMatrix['XCRR', 'CYWG', c('XCWR','XCGR','XCRR','XCRB',
                            'CCWR','CCGR','CCRR','CCRB',
                            'XYWR','XYGR','XYRR','XYRB',
                            'CYWR','CYGR','CYRR','CYRB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['XCRR', 'CYWR', c('XCWR','XCRR','CCWR','CCRR',
                            'XYWR','XYRR','CYWR','CYRR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCRR', 'CYWB', c('XCWR','XCRB','CCWR','CCRB',
                            'XYWR','XYRB','CYWR','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCRR', 'CYGG', c('XCGR','CCGR',
                            'XYGR','CYGR')] <- c(1,1,1,1)/4
  tMatrix['XCRR', 'CYGR', c('XCGR','XCRR','CCGR','CCRR',
                            'XYGR','XYRR','CYGR','CYRR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCRR', 'CYGB', c('XCGR','XCRB','CCGR','CCRB',
                            'XYGR','XYRB','CYGR','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCRR', 'CYRR', c('XCRR','CCRR',
                            'XYRR','CYRR')] <- c(1,1,1,1)/4
  tMatrix['XCRR', 'CYRB', c('XCRR','XCRB','CCRR','CCRB',
                            'XYRR','XYRB','CYRR','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCRR', 'CYBB', c('XCRB','CCRB',
                            'XYRB','CYRB')] <- c(1,1,1,1)/4


  tMatrix['XCRB', 'XYWW', c('XXWR','XXWB','XCWR','XCWB',
                            'XYWR','XYWB','CYWR','CYWB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCRB', 'XYWG', c('XXWR','XXWB','XXGR','XXGB',
                            'XCWR','XCWB','XCGR','XCGB',
                            'XYWR','XYWB','XYGR','XYGB',
                            'CYWR','CYWB','CYGR','CYGB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCRB', 'XYWR', c('XXWR','XXWB','XXRR','XXRB',
                            'XCWR','XCWB','XCRR','XCRB',
                            'XYWR','XYWB','XYRR','XYRB',
                            'CYWR','CYWB','CYRR','CYRB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCRB', 'XYWB', c('XXWR','XXWB','XXRB','XXBB',
                            'XCWR','XCWB','XCRB','XCBB',
                            'XYWR','XYWB','XYRB','XYBB',
                            'CYWR','CYWB','CYRB','CYBB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCRB', 'XYGG', c('XXGR','XXGB','XCGR','XCGB',
                            'XYGR','XYGB','CYGR','CYGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCRB', 'XYGR', c('XXGR','XXGB','XXRR','XXRB',
                            'XCGR','XCGB','XCRR','XCRB',
                            'XYGR','XYGB','XYRR','XYRB',
                            'CYGR','CYGB','CYRR','CYRB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCRB', 'XYGB', c('XXGR','XXGB','XXRB','XXBB',
                            'XCGR','XCGB','XCRB','XCBB',
                            'XYGR','XYGB','XYRB','XYBB',
                            'CYGR','CYGB','CYRB','CYBB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCRB', 'XYRR', c('XXRR','XXRB','XCRR','XCRB',
                            'XYRR','XYRB','CYRR','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCRB', 'XYRB', c('XXRR','XXRB','XXBB',
                            'XCRR','XCRB','XCBB',
                            'XYRR','XYRB','XYBB',
                            'CYRR','CYRB','CYBB')] <- c(1,2,1,1,2,1,1,2,1,1,2,1)/16
  tMatrix['XCRB', 'XYBB', c('XXRB','XXBB','XCRB','XCBB',
                            'XYRB','XYBB','CYRB','CYBB')] <- c(1,1,1,1,1,1,1,1)/8

  tMatrix['XCRB', 'CYWW', c('XCWR','XCWB','CCWR','CCWB',
                            'XYWR','XYWB','CYWR','CYWB')] <- c(1,1,1,1,1,1,1,1)/8
  wrb7 <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + cM*(1-chM)*crM, 1-cM, 1+cM*chM, cM*(1-chM)*(1-crM))/16
  tMatrix['XCRB', 'CYWG', c('XCWR','XCGR','XCRR','XCRB','XCWB','XCGB','XCBB',
                            'CCWR','CCGR','CCRR','CCRB','CCWB','CCGB','CCBB',
                            'XYWR','XYGR','XYRR','XYRB','XYWB','XYGB','XYBB',
                            'CYWR','CYGR','CYRR','CYRB','CYWB','CYGB','CYBB')] <- c(wrb7,wrb7,wrb7,wrb7)
  tMatrix['XCRB', 'CYWR', c('XCWR','XCWB','XCRR','XCRB',
                            'CCWR','CCWB','CCRR','CCRB',
                            'XYWR','XYWB','XYRR','XYRB',
                            'CYWR','CYWB','CYRR','CYRB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCRB', 'CYWB', c('XCWR','XCWB','XCRB','XCBB',
                            'CCWR','CCWB','CCRB','CCBB',
                            'XYWR','XYWB','XYRB','XYBB',
                            'CYWR','CYWB','CYRB','CYBB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCRB', 'CYGG', c('XCGR','XCGB','CCGR','CCGB',
                            'XYGR','XYGB','CYGR','CYGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCRB', 'CYGR', c('XCGR','XCGB','XCRR','XCRB',
                            'CCGR','CCGB','CCRR','CCRB',
                            'XYGR','XYGB','XYRR','XYRB',
                            'CYGR','CYGB','CYRR','CYRB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCRB', 'CYGB', c('XCGR','XCGB','XCRB','XCBB',
                            'CCGR','CCGB','CCRB','CCBB',
                            'XYGR','XYGB','XYRB','XYBB',
                            'CYGR','CYGB','CYRB','CYBB')] <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/16
  tMatrix['XCRB', 'CYRR', c('XCRR','XCRB','CCRR','CCRB',
                            'XYRR','XYRB','CYRR','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCRB', 'CYRB', c('XCRR','XCRB','XCBB',
                            'CCRR','CCRB','CCBB',
                            'XYRR','XYRB','XYBB',
                            'CYRR','CYRB','CYBB')] <- c(1,2,1,1,2,1,1,2,1,1,2,1)/16
  tMatrix['XCRB', 'CYBB', c('XCRB','XCBB','CCRB','CCBB',
                            'XYRB','XYBB','CYRB','CYBB')] <- c(1,1,1,1,1,1,1,1)/8


  tMatrix['XCBB', 'XYWW', c('XXWB','XCWB',
                            'XYWB','CYWB')] <- c(1,1,1,1)/4
  tMatrix['XCBB', 'XYWG', c('XXWB','XXGB','XCWB','XCGB',
                            'XYWB','XYGB','CYWB','CYGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCBB', 'XYWR', c('XXWB','XXRB','XCWB','XCRB',
                            'XYWB','XYRB','CYWB','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCBB', 'XYWB', c('XXWB','XXBB','XCWB','XCBB',
                            'XYWB','XYBB','CYWB','CYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCBB', 'XYGG', c('XXGB','XCGB',
                            'XYGB','CYGB')] <- c(1,1,1,1)/4
  tMatrix['XCBB', 'XYGR', c('XXGB','XXRB','XCGB','XCRB',
                            'XYGB','XYRB','CYGB','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCBB', 'XYGB', c('XXGB','XXBB','XCGB','XCBB',
                            'XYGB','XYBB','CYGB','CYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCBB', 'XYRR', c('XXRB','XCRB',
                            'XYRB','CYRB')] <- c(1,1,1,1)/4
  tMatrix['XCBB', 'XYRB', c('XXRB','XXBB','XCRB','XCBB',
                            'XYRB','XYBB','CYRB','CYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCBB', 'XYBB', c('XXBB','XCBB',
                            'XYBB','CYBB')] <- c(1,1,1,1)/4

  tMatrix['XCBB', 'CYWW', c('XCWB','CCWB',
                            'XYWB','CYWB')] <- c(1,1,1,1)/4
  tMatrix['XCBB', 'CYWG', c('XCWB','XCGB','XCRB','XCBB',
                            'CCWB','CCGB','CCRB','CCBB',
                            'XYWB','XYGB','XYRB','XYBB',
                            'CYWB','CYGB','CYRB','CYBB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['XCBB', 'CYWR', c('XCWB','XCRB','CCWB','CCRB',
                            'XYWB','XYRB','CYWB','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCBB', 'CYWB', c('XCWB','XCBB','CCWB','CCBB',
                            'XYWB','XYBB','CYWB','CYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCBB', 'CYGG', c('XCGB','CCGB',
                            'XYGB','CYGB')] <- c(1,1,1,1)/4
  tMatrix['XCBB', 'CYGR', c('XCGB','XCRB','CCGB','CCRB',
                            'XYGB','XYRB','CYGB','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCBB', 'CYGB', c('XCGB','XCBB','CCGB','CCBB',
                            'XYGB','XYBB','CYGB','CYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCBB', 'CYRR', c('XCRB','CCRB',
                            'XYRB','CYRB')] <- c(1,1,1,1)/4
  tMatrix['XCBB', 'CYRB', c('XCRB','XCBB','CCRB','CCBB',
                            'XYRB','XYBB','CYRB','CYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['XCBB', 'CYBB', c('XCBB','CCBB',
                            'XYBB','CYBB')] <- c(1,1,1,1)/4


  tMatrix['CCWW', 'XYWW', c('XCWW','CYWW')] <- c(1,1)/2
  tMatrix['CCWW', 'XYWG', c('XCWW','XCWG',
                            'CYWW','CYWG')] <- c(1,1,1,1)/4
  tMatrix['CCWW', 'XYWR', c('XCWW','XCWR',
                            'CYWW','CYWR')] <- c(1,1,1,1)/4
  tMatrix['CCWW', 'XYWB', c('XCWW','XCWB',
                            'CYWW','CYWB')] <- c(1,1,1,1)/4
  tMatrix['CCWW', 'XYGG', c('XCWG','CYWG')] <- c(1,1)/2
  tMatrix['CCWW', 'XYGR', c('XCWG','XCWR',
                            'CYWG','CYWR')] <- c(1,1,1,1)/4
  tMatrix['CCWW', 'XYGB', c('XCWG','XCWB',
                            'CYWG','CYWB')] <- c(1,1,1,1)/4
  tMatrix['CCWW', 'XYRR', c('XCWR','CYWR')] <- c(1,1)/2
  tMatrix['CCWW', 'XYRB', c('XCWR','XCWB',
                            'CYWR','CYWB')] <- c(1,1,1,1)/4
  tMatrix['CCWW', 'XYBB', c('XCWB','CYWB')] <- c(1,1)/2

  tMatrix['CCWW', 'CYWW', c('CCWW','CYWW')] <- c(1,1)/2
  tMatrix['CCWW', 'CYWG', c('CCWW','CCWG','CCWR','CCWB',
                            'CYWW','CYWG','CYWR','CYWB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/4
  tMatrix['CCWW', 'CYWR', c('CCWW','CCWR',
                            'CYWW','CYWR')] <- c(1,1,1,1)/4
  tMatrix['CCWW', 'CYWB', c('CCWW','CCWB',
                            'CYWW','CYWB')] <- c(1,1,1,1)/4
  tMatrix['CCWW', 'CYGG', c('CCWG','CYWG')] <- c(1,1)/2
  tMatrix['CCWW', 'CYGR', c('CCWG','CCWR',
                            'CYWG','CYWR')] <- c(1,1,1,1)/4
  tMatrix['CCWW', 'CYGB', c('CCWG','CCWB',
                            'CYWG','CYWB')] <- c(1,1,1,1)/4
  tMatrix['CCWW', 'CYRR', c('CCWR','CYWR')] <- c(1,1)/2
  tMatrix['CCWW', 'CYRB', c('CCWR','CCWB',
                            'CYWR','CYWB')] <- c(1,1,1,1)/4
  tMatrix['CCWW', 'CYBB', c('CCWB','CYWB')] <- c(1,1)/2


  wrbTen <- c((1-ccF)*(1-ddW), (1+ccF*cchF)*(1-ddW), (ccF*(1-cchF)*ccrF)*(1-ddW) + (1-ccF)*ddW*ddrW,
              (ccF*(1-cchF)*(1-ccrF))*(1-ddW) + (1-ccF)*ddW*(1-ddrW), (1+ccF*cchF)*ddW*ddhW,
              (1+ccF*cchF)*ddW*(1-ddhW)*ddrW, (1+ccF*cchF)*ddW*(1-ddhW)*(1-ddrW),
              (ccF*(1-cchF)*ccrF)*ddW*ddrW, (ccF*(1-cchF)*ccrF)*ddW*(1-ddrW) + (ccF*(1-cchF)*(1-ccrF))*ddW*ddrW,
              (ccF*(1-cchF)*(1-ccrF))*ddW*(1-ddrW))/4
  tMatrix['CCWG', 'XYWW', ] <- c(eTen, wrbTen, eTen, eTen, wrbTen)
  tMatrix['CCWG', 'CYWW', ] <- c(eTen, eTen, wrbTen, eTen, wrbTen)

  wrbTen <- c((1-ccF)*(1-ddW), 1-ccF + (1+ccF*cchF)*(1-ddW), (ccF*(1-cchF)*ccrF)*(1-ddW) + (1-ccF)*ddW*ddrW,
              (ccF*(1-cchF)*(1-ccrF))*(1-ddW) + (1-ccF)*ddW*(1-ddrW), 1+ccF*cchF + (1+ccF*cchF)*ddW*ddhW,
              ccF*(1-cchF)*ccrF + (1+ccF*cchF)*ddW*(1-ddhW)*ddrW, ccF*(1-cchF)*(1-ccrF) + (1+ccF*cchF)*ddW*(1-ddhW)*(1-ddrW), (ccF*(1-cchF)*ccrF)*ddW*ddrW,
              (ccF*(1-cchF)*ccrF)*ddW*(1-ddrW) + (ccF*(1-cchF)*(1-ccrF))*ddW*ddrW, (ccF*(1-cchF)*(1-ccrF))*ddW*(1-ddrW))/8
  tMatrix['CCWG', 'XYWG', ] <- c(eTen, wrbTen, eTen, eTen, wrbTen)

  wrbTen <- c((1-ccF)*(1-ddW), (1+ccF*cchF)*(1-ddW), 1-ccF + (ccF*(1-cchF)*ccrF)*(1-ddW) + (1-ccF)*ddW*ddrW,
              (ccF*(1-cchF)*(1-ccrF))*(1-ddW) + (1-ccF)*ddW*(1-ddrW), (1+ccF*cchF)*ddW*ddhW,
              1+ccF*cchF + (1+ccF*cchF)*ddW*(1-ddhW)*ddrW, (1+ccF*cchF)*ddW*(1-ddhW)*(1-ddrW), ccF*(1-cchF)*ccrF + (ccF*(1-cchF)*ccrF)*ddW*ddrW,
              ccF*(1-cchF)*(1-ccrF) + (ccF*(1-cchF)*ccrF)*ddW*(1-ddrW) + (ccF*(1-cchF)*(1-ccrF))*ddW*ddrW, (ccF*(1-cchF)*(1-ccrF))*ddW*(1-ddrW))/8
  tMatrix['CCWG', 'XYWR', ] <- c(eTen, wrbTen, eTen, eTen, wrbTen)
  tMatrix['CCWG', 'CYWR', ] <- c(eTen, eTen, wrbTen, eTen, wrbTen)

  wrbTen <- c((1-ccF)*(1-ddW), (1+ccF*cchF)*(1-ddW), (ccF*(1-cchF)*ccrF)*(1-ddW) + (1-ccF)*ddW*ddrW,
              1-ccF + (ccF*(1-cchF)*(1-ccrF))*(1-ddW) + (1-ccF)*ddW*(1-ddrW), (1+ccF*cchF)*ddW*ddhW,
              (1+ccF*cchF)*ddW*(1-ddhW)*ddrW, 1+ccF*cchF + (1+ccF*cchF)*ddW*(1-ddhW)*(1-ddrW), (ccF*(1-cchF)*ccrF)*ddW*ddrW,
              ccF*(1-cchF)*ccrF + (ccF*(1-cchF)*ccrF)*ddW*(1-ddrW) + (ccF*(1-cchF)*(1-ccrF))*ddW*ddrW, ccF*(1-cchF)*(1-ccrF) + (ccF*(1-cchF)*(1-ccrF))*ddW*(1-ddrW))/8
  tMatrix['CCWG', 'XYWB', ] <- c(eTen, wrbTen, eTen, eTen, wrbTen)
  tMatrix['CCWG', 'CYWB', ] <- c(eTen, eTen, wrbTen, eTen, wrbTen)

  wrbTen <- c(0, 1-ccF, 0, 0, 1+ccF*cchF, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF), 0, 0, 0)/4
  tMatrix['CCWG', 'XYGG', ] <- c(eTen, wrbTen, eTen, eTen, wrbTen)
  tMatrix['CCWG', 'CYGG', ] <- c(eTen, eTen, wrbTen, eTen, wrbTen)

  wrbTen <- c(0, 1-ccF, 1-ccF, 0, 1+ccF*cchF, ccF*(1-cchF)*ccrF + 1+ccF*cchF,
              ccF*(1-cchF)*(1-ccrF),ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF), 0)/8
  tMatrix['CCWG', 'XYGR', ] <- c(eTen, wrbTen, eTen, eTen, wrbTen)
  tMatrix['CCWG', 'CYGR', ] <- c(eTen, eTen, wrbTen, eTen, wrbTen)

  wrbTen <- c(0, 1-ccF, 0, 1-ccF, 1+ccF*cchF, ccF*(1-cchF)*ccrF,
              ccF*(1-cchF)*(1-ccrF) + 1+ccF*cchF, 0, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF))/8
  tMatrix['CCWG', 'XYGB', ] <- c(eTen, wrbTen, eTen, eTen, wrbTen)
  tMatrix['CCWG', 'CYGB', ] <- c(eTen, eTen, wrbTen, eTen, wrbTen)

  wrbTen <- c(0, 0, 1-ccF, 0, 0, 1+ccF*cchF, 0, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF), 0)/4
  tMatrix['CCWG', 'XYRR', ] <- c(eTen, wrbTen, eTen, eTen, wrbTen)
  tMatrix['CCWG', 'CYRR', ] <- c(eTen, eTen, wrbTen, eTen, wrbTen)

  wrbTen <- c(0, 0, 1-ccF, 1-ccF, 0, 1+ccF*cchF, 1+ccF*cchF, ccF*(1-cchF)*ccrF,
              ccF*(1-cchF)*(1-ccrF) + ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF))/8
  tMatrix['CCWG', 'XYRB', ] <- c(eTen, wrbTen, eTen, eTen, wrbTen)
  tMatrix['CCWG', 'CYRB', ] <- c(eTen, eTen, wrbTen, eTen, wrbTen)

  wrbTen <- c(0, 0, 0, 1-ccF, 0, 0, 1+ccF*cchF, 0, ccF*(1-cchF)*ccrF, ccF*(1-cchF)*(1-ccrF))/4
  tMatrix['CCWG', 'XYBB', ] <- c(eTen, wrbTen, eTen, eTen, wrbTen)
  tMatrix['CCWG', 'CYBB', ] <- c(eTen, eTen, wrbTen, eTen, wrbTen)

  wrbTen <- c((1-ccF)*(1-cM)*(1-ddW), (1+ccF*cchF)*(1-cM)*(1-ddW) + (1-ccF)*(1+cM*chM),
              (1-ccF)*(1-cM)*ddW*ddrW + (ccF*(1-cchF)*ccrF)*(1-cM)*(1-ddW) + (1-ccF)*(cM*(1-chM)*crM),
              (1-ccF)*(1-cM)*ddW*(1-ddrW) + (ccF*(1-cchF)*(1-ccrF))*(1-cM)*(1-ddW) + (1-ccF)*(cM*(1-chM)*(1-crM)),
              (1+ccF*cchF)*(1-cM)*ddW*ddhW + (1+ccF*cchF)*(1+cM*chM),
              (1+ccF*cchF)*(1-cM)*ddW*(1-ddhW)*ddrW + (ccF*(1-cchF)*ccrF)*(1+cM*chM) + (1+ccF*cchF)*(cM*(1-chM)*crM),
              (1+ccF*cchF)*(1-cM)*ddW*(1-ddhW)*(1-ddrW) + (ccF*(1-cchF)*(1-ccrF))*(1+cM*chM) + (1+ccF*cchF)*(cM*(1-chM)*(1-crM)),
              (ccF*(1-cchF)*ccrF)*(1-cM)*ddW*ddrW + (ccF*(1-cchF)*ccrF)*(cM*(1-chM)*crM),
              (ccF*(1-cchF)*ccrF)*(1-cM)*ddW*(1-ddrW) + (ccF*(1-cchF)*(1-ccrF))*(1-cM)*ddW*ddrW + (ccF*(1-cchF)*(1-ccrF))*(cM*(1-chM)*crM) + (ccF*(1-cchF)*ccrF)*(cM*(1-chM)*(1-crM)),
              (ccF*(1-cchF)*(1-ccrF))*(1-cM)*ddW*(1-ddrW) + (ccF*(1-cchF)*(1-ccrF))*(cM*(1-chM)*(1-crM)))/8
  tMatrix['CCWG', 'CYWG', ] <- c(eTen, eTen, wrbTen, eTen, wrbTen)


  tMatrix['CCWR', 'XYWW', c('XCWW','XCWR',
                            'CYWW','CYWR')] <- c(1,1,1,1)/4
  tMatrix['CCWR', 'XYWG', c('XCWW','XCWR','XCWG','XCGR',
                            'CYWW','CYWR','CYWG','CYGR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCWR', 'XYWR', c('XCWW','XCWR','XCRR',
                            'CYWW','CYWR','CYRR')] <- c(1,2,1,1,2,1)/8
  tMatrix['CCWR', 'XYWB', c('XCWW','XCWR','XCWB','XCRB',
                            'CYWW','CYWR','CYWB','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCWR', 'XYGG', c('XCWG','XCGR',
                            'CYWG','CYGR')] <- c(1,1,1,1)/4
  tMatrix['CCWR', 'XYGR', c('XCWG','XCGR','XCWR','XCRR',
                            'CYWG','CYGR','CYWR','CYRR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCWR', 'XYGB', c('XCWG','XCGR','XCWB','XCRB',
                            'CYWG','CYGR','CYWB','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCWR', 'XYRR', c('XCWR','XCRR',
                            'CYWR','CYRR')] <- c(1,1,1,1)/4
  tMatrix['CCWR', 'XYRB', c('XCWR','XCRR','XCWB','XCRB',
                            'CYWR','CYRR','CYWB','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCWR', 'XYBB', c('XCWB','XCRB',
                            'CYWB','CYRB')] <- c(1,1,1,1)/4

  tMatrix['CCWR', 'CYWW', c('CCWW','CCWR',
                            'CYWW','CYWR')] <- c(1,1,1,1)/4
  tMatrix['CCWR', 'CYWG', c('CCWW','CCWG','CCWR','CCWB','CCGR','CCRR','CCRB',
                            'CYWW','CYWG','CYWR','CYWB','CYGR','CYRR','CYRB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM + 1-cM, cM*(1-chM)*(1-crM), 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                                                    1-cM, 1+cM*chM, cM*(1-chM)*crM + 1-cM, cM*(1-chM)*(1-crM), 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['CCWR', 'CYWR', c('CCWW','CCWR','CCRR',
                            'CYWW','CYWR','CYRR')] <- c(1,2,1,1,2,1)/8
  tMatrix['CCWR', 'CYWB', c('CCWW','CCWR','CCWB','CCRB',
                            'CYWW','CYWR','CYWB','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCWR', 'CYGG', c('CCWG','CCGR',
                            'CYWG','CYGR')] <- c(1,1,1,1)/4
  tMatrix['CCWR', 'CYGR', c('CCWG','CCGR','CCWR','CCRR',
                            'CYWG','CYGR','CYWR','CYRR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCWR', 'CYGB', c('CCWG','CCGR','CCWB','CCRB',
                            'CYWG','CYGR','CYWB','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCWR', 'CYRR', c('CCWR','CCRR',
                            'CYWR','CYRR')] <- c(1,1,1,1)/4
  tMatrix['CCWR', 'CYRB', c('CCWR','CCRR','CCWB','CCRB',
                            'CYWR','CYRR','CYWB','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCWR', 'CYBB', c('CCWB','CCRB',
                            'CYWB','CYRB')] <- c(1,1,1,1)/4


  tMatrix['CCWB', 'XYWW', c('XCWW','XCWB',
                            'CYWW','CYWB')] <- c(1,1,1,1)/4
  tMatrix['CCWB', 'XYWG', c('XCWW','XCWB','XCWG','XCGB',
                            'CYWW','CYWB','CYWG','CYGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCWB', 'XYWR', c('XCWW','XCWB','XCWR','XCRB',
                            'CYWW','CYWB','CYWR','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCWB', 'XYWB', c('XCWW','XCWB','XCBB',
                            'CYWW','CYWB','CYBB')] <- c(1,2,1,1,2,1)/8
  tMatrix['CCWB', 'XYGG', c('XCWG','XCGB',
                            'CYWG','CYGB')] <- c(1,1,1,1)/4
  tMatrix['CCWB', 'XYGR', c('XCWG','XCGB','XCWR','XCGR',
                            'CYWG','CYGB','CYWR','CYGR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCWB', 'XYGB', c('XCWG','XCGB','XCWB','XCBB',
                            'CYWG','CYGB','CYWB','CYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCWB', 'XYRR', c('XCWR','XCRB',
                            'CYWR','CYRB')] <- c(1,1,1,1)/4
  tMatrix['CCWB', 'XYRB', c('XCWR','XCRB','XCWB','XCBB',
                            'CYWR','CYRB','CYWB','CYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCWB', 'XYBB', c('XCWB','XCBB',
                            'CYWB','CYBB')] <- c(1,1,1,1)/4

  tMatrix['CCWB', 'CYWW', c('CCWW','CCWB',
                            'CYWW','CYWB')] <- c(1,1,1,1)/4
  tMatrix['CCWB', 'CYWG', c('CCWW','CCWG','CCWR','CCWB','CCGB','CCRB','CCBB',
                            'CYWW','CYWG','CYWR','CYWB','CYGB','CYRB','CYBB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + 1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                                                    1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + 1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['CCWB', 'CYWR', c('CCWW','CCWB','CCWR','CCRB',
                            'CYWW','CYWB','CYWR','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCWB', 'CYWB', c('CCWW','CCWB','CCBB',
                            'CYWW','CYWB','CYBB')] <- c(1,2,1,1,2,1)/8
  tMatrix['CCWB', 'CYGG', c('CCWG','CCGB',
                            'CYWG','CYGB')] <- c(1,1,1,1)/4
  tMatrix['CCWB', 'CYGR', c('CCWG','CCGB','CCWR','CCGR',
                            'CYWG','CYGB','CYWR','CYGR')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCWB', 'CYGB', c('CCWG','CCGB','CCWB','CCBB',
                            'CYWG','CYGB','CYWB','CYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCWB', 'CYRR', c('CCWR','CCRB',
                            'CYWR','CYRB')] <- c(1,1,1,1)/4
  tMatrix['CCWB', 'CYRB', c('CCWR','CCRB','CCWB','CCBB',
                            'CYWR','CYRB','CYWB','CYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCWB', 'CYBB', c('CCWB','CCBB',
                            'CYWB','CYBB')] <- c(1,1,1,1)/4


  tMatrix['CCGG', 'XYWW', c('XCWG','XCGG','XCGR','XCGB',
                            'CYWG','CYGG','CYGR','CYGB')] <- c(1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW),
                                                               1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW))/2
  tMatrix['CCGG', 'XYWG', c('XCWG','XCGG','XCGR','XCGB',
                            'CYWG','CYGG','CYGR','CYGB')] <- c(1-ddW, 1 + ddW*ddhW, ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW),
                                                               1-ddW, 1 + ddW*ddhW, ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW))/4
  tMatrix['CCGG', 'XYWR', c('XCWG','XCGG','XCGR','XCGB',
                            'CYWG','CYGG','CYGR','CYGB')] <- c(1-ddW, ddW*ddhW, 1 + ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW),
                                                               1-ddW, ddW*ddhW, 1 + ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW))/4
  tMatrix['CCGG', 'XYWB', c('XCWG','XCGG','XCGR','XCGB',
                            'CYWG','CYGG','CYGR','CYGB')] <- c(1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, 1 + ddW*(1-ddhW)*(1-ddrW),
                                                               1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, 1 + ddW*(1-ddhW)*(1-ddrW))/4
  tMatrix['CCGG', 'XYGG', c('XCGG','CYGG')] <- c(1,1)/2
  tMatrix['CCGG', 'XYGR', c('XCGG','XCGR',
                            'CYGG','CYGR')] <- c(1,1,1,1)/4
  tMatrix['CCGG', 'XYGB', c('XCGG','XCGB',
                            'CYGG','CYGB')] <- c(1,1,1,1)/4
  tMatrix['CCGG', 'XYRR', c('XCGR','CYGR')] <- c(1,1)/2
  tMatrix['CCGG', 'XYRB', c('XCGR','XCGB',
                            'CYGR','CYGB')] <- c(1,1,1,1)/4
  tMatrix['CCGG', 'XYBB', c('XCGB','CYGB')] <- c(1,1)/2

  tMatrix['CCGG', 'CYWW', c('CCWG','CCGG','CCGR','CCGB',
                            'CYWG','CYGG','CYGR','CYGB')] <- c(1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW),
                                                               1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW))/2
  tMatrix['CCGG', 'CYWG', c('CCWG','CCGG','CCGR','CCGB',
                            'CYWG','CYGG','CYGR','CYGB')] <- c((1-cM)*(1-ddW), 1+cM*chM + (1-cM)*ddW*ddhW, cM*(1-chM)*crM + (1-cM)*ddW*(1-ddhW)*ddrW, cM*(1-chM)*(1-crM) + (1-cM)*ddW*(1-ddhW)*(1-ddrW),
                                                               (1-cM)*(1-ddW), 1+cM*chM + (1-cM)*ddW*ddhW, cM*(1-chM)*crM + (1-cM)*ddW*(1-ddhW)*ddrW, cM*(1-chM)*(1-crM) + (1-cM)*ddW*(1-ddhW)*(1-ddrW))/4
  tMatrix['CCGG', 'CYWR', c('CCWG','CCGG','CCGR','CCGB',
                            'CYWG','CYGG','CYGR','CYGB')] <- c(1-ddW, ddW*ddhW, 1 + ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW),
                                                               1-ddW, ddW*ddhW, 1 + ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW))/4
  tMatrix['CCGG', 'CYWB', c('CCWG','CCGG','CCGR','CCGB',
                            'CYWG','CYGG','CYGR','CYGB')] <- c(1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, 1 + ddW*(1-ddhW)*(1-ddrW),
                                                               1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, 1 + ddW*(1-ddhW)*(1-ddrW))/4
  tMatrix['CCGG', 'CYGG', c('CCGG','CYGG')] <- c(1,1)/2
  tMatrix['CCGG', 'CYGR', c('CCGG','CCGR',
                            'CYGG','CYGR')] <- c(1,1,1,1)/4
  tMatrix['CCGG', 'CYGB', c('CCGG','CCGB',
                            'CYGG','CYGB')] <- c(1,1,1,1)/4
  tMatrix['CCGG', 'CYRR', c('CCGR','CYGR')] <- c(1,1)/2
  tMatrix['CCGG', 'CYRB', c('CCGR','CCGB',
                            'CYGR','CYGB')] <- c(1,1,1,1)/4
  tMatrix['CCGG', 'CYBB', c('CCGB','CYGB')] <- c(1,1)/2


  tMatrix['CCGR', 'XYWW', c('XCWG','XCGG','XCGR','XCGB','XCWR','XCRR','XCRB',
                            'CYWG','CYGG','CYGR','CYGB','CYWR','CYRR','CYRB')] <- c(1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW), 1-ddW, ddW*ddrW, ddW*(1-ddrW),
                                                                                    1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW), 1-ddW, ddW*ddrW, ddW*(1-ddrW))/4
  tMatrix['CCGR', 'XYWG', c('XCWG','XCGG','XCGR','XCGB','XCWR','XCRR','XCRB',
                            'CYWG','CYGG','CYGR','CYGB','CYWR','CYRR','CYRB')] <- c(1-ddW, 1 + ddW*ddhW, 1 + ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW), 1-ddW, ddW*ddrW, ddW*(1-ddrW),
                                                                                    1-ddW, 1 + ddW*ddhW, 1 + ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW), 1-ddW, ddW*ddrW, ddW*(1-ddrW))/8
  tMatrix['CCGR', 'XYWR', c('XCWG','XCGG','XCGR','XCGB','XCWR','XCRR','XCRB',
                            'CYWG','CYGG','CYGR','CYGB','CYWR','CYRR','CYRB')] <- c(1-ddW, ddW*ddhW, 1 + ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW), 1-ddW, 1 + ddW*ddrW, ddW*(1-ddrW),
                                                                                    1-ddW, ddW*ddhW, 1 + ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW), 1-ddW, 1 + ddW*ddrW, ddW*(1-ddrW))/8
  tMatrix['CCGR', 'XYWB', c('XCWG','XCGG','XCGR','XCGB','XCWR','XCRR','XCRB',
                            'CYWG','CYGG','CYGR','CYGB','CYWR','CYRR','CYRB')] <- c(1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, 1 + ddW*(1-ddhW)*(1-ddrW), 1-ddW, ddW*ddrW, 1 + ddW*(1-ddrW),
                                                                                    1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, 1 + ddW*(1-ddhW)*(1-ddrW), 1-ddW, ddW*ddrW, 1 + ddW*(1-ddrW))/8
  tMatrix['CCGR', 'XYGG', c('XCGG','XCGR',
                            'CYGG','CYGR')] <- c(1,1,1,1)/4
  tMatrix['CCGR', 'XYGR', c('XCGG','XCGR','XCRR',
                            'CYGG','CYGR','CYRR')] <- c(1,2,1,1,2,1)/8
  tMatrix['CCGR', 'XYGB', c('XCGG','XCGR','XCGB','XCRB',
                            'CYGG','CYGR','CYGB','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCGR', 'XYRR', c('XCGR','XCRR',
                            'CYGR','CYRR')] <- c(1,1,1,1)/4
  tMatrix['CCGR', 'XYRB', c('XCGR','XCGB','XCRR','XCRB',
                            'CYGR','CYGB','CYRR','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCGR', 'XYBB', c('XCGB','XCRB',
                            'CYGB','CYRB')] <- c(1,1,1,1)/4

  tMatrix['CCGR', 'CYWW', c('CCWG','CCGG','CCGR','CCGB','CCWR','CCRR','CCRB',
                            'CYWG','CYGG','CYGR','CYGB','CYWR','CYRR','CYRB')] <- c(1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW), 1-ddW, ddW*ddrW, ddW*(1-ddrW),
                                                                                    1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW), 1-ddW, ddW*ddrW, ddW*(1-ddrW))/4
  wrb7 <- c((1-cM)*(1-ddW), 1+cM*chM + (1-cM)*ddW*ddhW, cM*(1-chM)*crM + 1+cM*chM + (1-cM)*ddW*(1-ddhW)*ddrW, cM*(1-chM)*(1-crM) + (1-cM)*ddW*(1-ddhW)*(1-ddrW), (1-cM)*(1-ddW), cM*(1-chM)*crM + (1-cM)*ddW*ddrW, cM*(1-chM)*(1-crM) + (1-cM)*ddW*(1-ddrW))/8
  tMatrix['CCGR', 'CYWG', c('CCWG','CCGG','CCGR','CCGB','CCWR','CCRR','CCRB',
                            'CYWG','CYGG','CYGR','CYGB','CYWR','CYRR','CYRB')] <- c(wrb7,wrb7)
  tMatrix['CCGR', 'CYWR', c('CCWG','CCGG','CCGR','CCGB','CCWR','CCRR','CCRB',
                            'CYWG','CYGG','CYGR','CYGB','CYWR','CYRR','CYRB')] <- c(1-ddW, ddW*ddhW, 1 + ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW), 1-ddW, 1 + ddW*ddrW, ddW*(1-ddrW),
                                                                                    1-ddW, ddW*ddhW, 1 + ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW), 1-ddW, 1 + ddW*ddrW, ddW*(1-ddrW))/8
  tMatrix['CCGR', 'CYWB', c('CCWG','CCGG','CCGR','CCGB','CCWR','CCRR','CCRB',
                            'CYWG','CYGG','CYGR','CYGB','CYWR','CYRR','CYRB')] <- c(1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, 1 + ddW*(1-ddhW)*(1-ddrW), 1-ddW, ddW*ddrW, 1 + ddW*(1-ddrW),
                                                                                    1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, 1 + ddW*(1-ddhW)*(1-ddrW), 1-ddW, ddW*ddrW, 1 + ddW*(1-ddrW))/8
  tMatrix['CCGR', 'CYGG', c('CCGG','CCGR',
                            'CYGG','CYGR')] <- c(1,1,1,1)/4
  tMatrix['CCGR', 'CYGR', c('CCGG','CCGR','CCRR',
                            'CYGG','CYGR','CYRR')] <- c(1,2,1,1,2,1)/8
  tMatrix['CCGR', 'CYGB', c('CCGG','CCGR','CCGB','CCRB',
                            'CYGG','CYGR','CYGB','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCGR', 'CYRR', c('CCGR','CCRR',
                            'CYGR','CYRR')] <- c(1,1,1,1)/4
  tMatrix['CCGR', 'CYRB', c('CCGR','CCGB','CCRR','CCRB',
                            'CYGR','CYGB','CYRR','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCGR', 'CYBB', c('CCGB','CCRB',
                            'CYGB','CYRB')] <- c(1,1,1,1)/4


  tMatrix['CCGB', 'XYWW', c('XCWG','XCGG','XCGR','XCGB','XCWB','XCRB','XCBB',
                            'CYWG','CYGG','CYGR','CYGB','CYWB','CYRB','CYBB')] <- c(1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW), 1-ddW, ddW*ddrW, ddW*(1-ddrW),
                                                                                    1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW), 1-ddW, ddW*ddrW, ddW*(1-ddrW))/4
  tMatrix['CCGB', 'XYWG', c('XCWG','XCGG','XCGR','XCGB','XCWB','XCRB','XCBB',
                            'CYWG','CYGG','CYGR','CYGB','CYWB','CYRB','CYBB')] <- c(1-ddW, 1 + ddW*ddhW, 1 + ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW), 1-ddW, ddW*ddrW, ddW*(1-ddrW),
                                                                                    1-ddW, 1 + ddW*ddhW, 1 + ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW), 1-ddW, ddW*ddrW, ddW*(1-ddrW))/8
  tMatrix['CCGB', 'XYWR', c('XCWG','XCGG','XCGR','XCGB','XCWB','XCRB','XCBB',
                            'CYWG','CYGG','CYGR','CYGB','CYWB','CYRB','CYBB')] <- c(1-ddW, ddW*ddhW, 1 + ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW), 1-ddW, 1 + ddW*ddrW, ddW*(1-ddrW),
                                                                                    1-ddW, ddW*ddhW, 1 + ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW), 1-ddW, 1 + ddW*ddrW, ddW*(1-ddrW))/8
  tMatrix['CCGB', 'XYWB', c('XCWG','XCGG','XCGR','XCGB','XCWB','XCRB','XCBB',
                            'CYWG','CYGG','CYGR','CYGB','CYWB','CYRB','CYBB')] <- c(1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, 1 + ddW*(1-ddhW)*(1-ddrW), 1-ddW, ddW*ddrW, 1 + ddW*(1-ddrW),
                                                                                    1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, 1 + ddW*(1-ddhW)*(1-ddrW), 1-ddW, ddW*ddrW, 1 + ddW*(1-ddrW))/8
  tMatrix['CCGB', 'XYGG', c('XCGG','XCGB',
                            'CYGG','CYGB')] <- c(1,1,1,1)/4
  tMatrix['CCGB', 'XYGR', c('XCGG','XCGR','XCGB','XCRB',
                            'CYGG','CYGR','CYGB','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCGB', 'XYGB', c('XCGG','XCGB','XCBB',
                            'CYGG','CYGB','CYBB')] <- c(1,2,1,1,2,1)/8
  tMatrix['CCGB', 'XYRR', c('XCGR','XCRB',
                            'CYGR','CYRB')] <- c(1,1,1,1)/4
  tMatrix['CCGB', 'XYRB', c('XCGR','XCGB','XCRB','XCBB',
                            'CYGR','CYGB','CYRB','CYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCGB', 'XYBB', c('XCGB','XCBB',
                            'CYGB','CYBB')] <- c(1,1,1,1)/4

  tMatrix['CCGB', 'CYWW', c('CCWG','CCGG','CCGR','CCGB','CCWB','CCRB','CCBB',
                            'CYWG','CYGG','CYGR','CYGB','CYWB','CYRB','CYBB')] <- c(1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW), 1-ddW, ddW*ddrW, ddW*(1-ddrW),
                                                                                    1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW), 1-ddW, ddW*ddrW, ddW*(1-ddrW))/4
  wrb7 <- c((1-cM)*(1-ddW), 1+cM*chM + (1-cM)*ddW*ddhW, cM*(1-chM)*crM + (1-cM)*ddW*(1-ddhW)*ddrW, cM*(1-chM)*(1-crM) + 1+cM*chM + (1-cM)*ddW*(1-ddhW)*(1-ddrW), (1-cM)*(1-ddW), cM*(1-chM)*crM + (1-cM)*ddW*ddrW, cM*(1-chM)*(1-crM) + (1-cM)*ddW*(1-ddrW))/8
  tMatrix['CCGB', 'CYWG', c('CCWG','CCGG','CCGR','CCGB','CCWB','CCRB','CCBB',
                            'CYWG','CYGG','CYGR','CYGB','CYWB','CYRB','CYBB')] <- c(wrb7,wrb7)
  tMatrix['CCGB', 'CYWR', c('CCWG','CCGG','CCGR','CCGB','CCWB','CCRB','CCBB',
                            'CYWG','CYGG','CYGR','CYGB','CYWB','CYRB','CYBB')] <- c(1-ddW, ddW*ddhW, 1 + ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW), 1-ddW, 1 + ddW*ddrW, ddW*(1-ddrW),
                                                                                    1-ddW, ddW*ddhW, 1 + ddW*(1-ddhW)*ddrW, ddW*(1-ddhW)*(1-ddrW), 1-ddW, 1 + ddW*ddrW, ddW*(1-ddrW))/8
  tMatrix['CCGB', 'CYWB', c('CCWG','CCGG','CCGR','CCGB','CCWB','CCRB','CCBB',
                            'CYWG','CYGG','CYGR','CYGB','CYWB','CYRB','CYBB')] <- c(1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, 1 + ddW*(1-ddhW)*(1-ddrW), 1-ddW, ddW*ddrW, 1 + ddW*(1-ddrW),
                                                                                    1-ddW, ddW*ddhW, ddW*(1-ddhW)*ddrW, 1 + ddW*(1-ddhW)*(1-ddrW), 1-ddW, ddW*ddrW, 1 + ddW*(1-ddrW))/8
  tMatrix['CCGB', 'CYGG', c('CCGG','CCGB',
                            'CYGG','CYGB')] <- c(1,1,1,1)/4
  tMatrix['CCGB', 'CYGR', c('CCGG','CCGR','CCGB','CCRB',
                            'CYGG','CYGR','CYGB','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCGB', 'CYGB', c('CCGG','CCGB','CCBB',
                            'CYGG','CYGB','CYBB')] <- c(1,2,1,1,2,1)/8
  tMatrix['CCGB', 'CYRR', c('CCGR','CCRB',
                            'CYGR','CYRB')] <- c(1,1,1,1)/4
  tMatrix['CCGB', 'CYRB', c('CCGR','CCGB','CCRB','CCBB',
                            'CYGR','CYGB','CYRB','CYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCGB', 'CYBB', c('CCGB','CCBB',
                            'CYGB','CYBB')] <- c(1,1,1,1)/4


  tMatrix['CCRR', 'XYWW', c('XCWR','CYWR')] <- c(1,1)/2
  tMatrix['CCRR', 'XYWG', c('XCWR','XCGR',
                            'CYWR','CYGR')] <- c(1,1,1,1)/4
  tMatrix['CCRR', 'XYWR', c('XCWR','XCRR',
                            'CYWR','CYRR')] <- c(1,1,1,1)/4
  tMatrix['CCRR', 'XYWB', c('XCWR','XCRB',
                            'CYWR','CYRB')] <- c(1,1,1,1)/4
  tMatrix['CCRR', 'XYGG', c('XCGR','CYGR')] <- c(1,1)/2
  tMatrix['CCRR', 'XYGR', c('XCGR','XCRR',
                            'CYGR','CYRR')] <- c(1,1,1,1)/4
  tMatrix['CCRR', 'XYGB', c('XCGR','XCRB',
                            'CYGR','CYRB')] <- c(1,1,1,1)/4
  tMatrix['CCRR', 'XYRR', c('XCRR','CYRR')] <- c(1,1)/2
  tMatrix['CCRR', 'XYRB', c('XCRR','XCRB',
                            'CYRR','CYRB')] <- c(1,1,1,1)/4
  tMatrix['CCRR', 'XYBB', c('XCRB','CYRB')] <- c(1,1)/2

  tMatrix['CCRR', 'CYWW', c('CCWR','CYWR')] <- c(1,1)/2
  tMatrix['CCRR', 'CYWG', c('CCWR','CCGR','CCRR','CCRB',
                            'CYWR','CYGR','CYRR','CYRB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/4
  tMatrix['CCRR', 'CYWR', c('CCWR','CCRR',
                            'CYWR','CYRR')] <- c(1,1,1,1)/4
  tMatrix['CCRR', 'CYWB', c('CCWR','CCRB',
                            'CYWR','CYRB')] <- c(1,1,1,1)/4
  tMatrix['CCRR', 'CYGG', c('CCGR','CYGR')] <- c(1,1)/2
  tMatrix['CCRR', 'CYGR', c('CCGR','CCRR',
                            'CYGR','CYRR')] <- c(1,1,1,1)/4
  tMatrix['CCRR', 'CYGB', c('CCGR','CCRB',
                            'CYGR','CYRB')] <- c(1,1,1,1)/4
  tMatrix['CCRR', 'CYRR', c('CCRR','CYRR')] <- c(1,1)/2
  tMatrix['CCRR', 'CYRB', c('CCRR','CCRB',
                            'CYRR','CYRB')] <- c(1,1,1,1)/4
  tMatrix['CCRR', 'CYBB', c('CCRB','CYRB')] <- c(1,1)/2


  tMatrix['CCRB', 'XYWW', c('XCWR','XCWB',
                            'CYWR','CYWB')] <- c(1,1,1,1)/4
  tMatrix['CCRB', 'XYWG', c('XCWR','XCWB','XCGR','XCGB',
                            'CYWR','CYWB','CYGR','CYGB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCRB', 'XYWR', c('XCWR','XCWB','XCRR','XCRB',
                            'CYWR','CYWB','CYRR','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCRB', 'XYWB', c('XCWR','XCWB','XCRB','XCBB',
                            'CYWR','CYWB','CYRB','CYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCRB', 'XYGG', c('XCGR','XCGB',
                            'CYGR','CYGB')] <- c(1,1,1,1)/4
  tMatrix['CCRB', 'XYGR', c('XCGR','XCGB','XCRR','XCRB',
                            'CYGR','CYGB','CYRR','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCRB', 'XYGB', c('XCGR','XCGB','XCRB','XCBB',
                            'CYGR','CYGB','CYRB','CYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCRB', 'XYRR', c('XCRR','XCRB',
                          'CYRR','CYRB')] <- c(1,1,1,1)/4
  tMatrix['CCRB', 'XYRB', c('XCRR','XCRB','XCBB',
                            'CYRR','CYRB','CYBB')] <- c(1,2,1,1,2,1)/8
  tMatrix['CCRB', 'XYBB', c('XCRB','XCBB',
                            'CYRB','CYBB')] <- c(1,1,1,1)/4

  tMatrix['CCRB', 'CYWW', c('CCWR','CCWB',
                            'CYWR','CYWB')] <- c(1,1,1,1)/4
  tMatrix['CCRB', 'CYWG', c('CCWR','CCGR','CCRR','CCRB','CCWB','CCGB','CCBB',
                            'CYWR','CYGR','CYRR','CYRB','CYWB','CYGB','CYBB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + cM*(1-chM)*crM, 1-cM, 1+cM*chM, cM*(1-chM)*(1-crM),
                                                                                    1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + cM*(1-chM)*crM, 1-cM, 1+cM*chM, cM*(1-chM)*(1-crM))/8
  tMatrix['CCRB', 'CYWR', c('CCWR','CCWB','CCRR','CCRB',
                            'CYWR','CYWB','CYRR','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCRB', 'CYWB', c('CCWR','CCWB','CCRB','CCBB',
                            'CYWR','CYWB','CYRB','CYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCRB', 'CYGG', c('CCGR','CCGB',
                            'CYGR','CYGB')] <- c(1,1,1,1)/4
  tMatrix['CCRB', 'CYGR', c('CCGR','CCGB','CCRR','CCRB',
                            'CYGR','CYGB','CYRR','CYRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCRB', 'CYGB', c('CCGR','CCGB','CCRB','CCBB',
                            'CYGR','CYGB','CYRB','CYBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['CCRB', 'CYRR', c('CCRR','CCRB',
                            'CYRR','CYRB')] <- c(1,1,1,1)/4
  tMatrix['CCRB', 'CYRB', c('CCRR','CCRB','CCBB',
                            'CYRR','CYRB','CYBB')] <- c(1,2,1,1,2,1)/8
  tMatrix['CCRB', 'CYBB', c('CCRB','CCBB',
                            'CYRB','CYBB')] <- c(1,1,1,1)/4


  tMatrix['CCBB', 'XYWW', c('XCWB','CYWB')] <- c(1,1)/2
  tMatrix['CCBB', 'XYWG', c('XCWB','XCGB',
                            'CYWB','CYGB')] <- c(1,1,1,1)/4
  tMatrix['CCBB', 'XYWR', c('XCWB','XCRB',
                            'CYWB','CYRB')] <- c(1,1,1,1)/4
  tMatrix['CCBB', 'XYWB', c('XCWB','XCBB',
                            'CYWB','CYBB')] <- c(1,1,1,1)/4
  tMatrix['CCBB', 'XYGG', c('XCGB','CYGB')] <- c(1,1)/2
  tMatrix['CCBB', 'XYGR', c('XCGB','XCRB',
                            'CYGB','CYRB')] <- c(1,1,1,1)/4
  tMatrix['CCBB', 'XYGB', c('XCGB','XCBB',
                            'CYGB','CYBB')] <- c(1,1,1,1)/4
  tMatrix['CCBB', 'XYRR', c('XCRB','CYRB')] <- c(1,1)/2
  tMatrix['CCBB', 'XYRB', c('XCRB','XCBB',
                            'CYRB','CYBB')] <- c(1,1,1,1)/4
  tMatrix['CCBB', 'XYBB', c('XCBB','CYBB')] <- c(1,1)/2

  tMatrix['CCBB', 'CYWW', c('CCWB','CYWB')] <- c(1,1)/2
  tMatrix['CCBB', 'CYWG', c('CCWB','CCGB','CCRB','CCBB',
                            'CYWB','CYGB','CYRB','CYBB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM),
                                                               1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/4
  tMatrix['CCBB', 'CYWR', c('CCWB','CCRB',
                            'CYWB','CYRB')] <- c(1,1,1,1)/4
  tMatrix['CCBB', 'CYWB', c('CCWB','CCBB',
                            'CYWB','CYBB')] <- c(1,1,1,1)/4
  tMatrix['CCBB', 'CYGG', c('CCGB','CYGB')] <- c(1,1)/2
  tMatrix['CCBB', 'CYGR', c('CCGB','CCRB',
                            'CYGB','CYRB')] <- c(1,1,1,1)/4
  tMatrix['CCBB', 'CYGB', c('CCGB','CCBB',
                            'CYGB','CYBB')] <- c(1,1,1,1)/4
  tMatrix['CCBB', 'CYRR', c('CCRB','CYRB')] <- c(1,1)/2
  tMatrix['CCBB', 'CYRB', c('CCRB','CCBB',
                            'CYRB','CYBB')] <- c(1,1,1,1)/4
  tMatrix['CCBB', 'CYBB', c('CCBB','CYBB')] <- c(1,1)/2


  #protection from underflow errors
  tMatrix[tMatrix < .Machine$double.eps] <- 0


  ## initialize viability mask. No mother/father-specific death, so use basic mask
  viabilityMask <- array(data = 1, dim = c(size,size,size), dimnames = list(gtype, gtype, gtype))

  ## genotype-specific modifiers
  if(!is.null(phi)){
    stop("This cube has a special phi, due to being male/female specific.
         Please edit it after the cube is built, if you are absolutely sure it needs to change.")
  }
  phi = setNames(object = c(rep.int(x = 1, times = 30),rep.int(x = 0, times = 20)), nm = gtype)
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
    releaseType = 'CYGG'
  ))

}
