###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   Immunizing Reversal - X-Linked
#   Héctor Sanchez, Jared Bennett, Sean Wu, John Marshall
#   jared_bennett@berkeley.edu
#   August 2017
#   December 2018
#    Update to reflect cutting, homing, resistance generation rates
#    Removed male deposition stuff, it doesn't happen.
#    Added female deposition against W allele, not sure why it wasn't there before.
#
###############################################################################

#' Inheritance Cube: Immunizing Reversal
#'
#' This function creates an Immunizing Reversal construct, it has 5 alleles at 1 locus
#'  * W: Wild-type
#'  * H: Homing allele
#'  * E: Eraser allele
#'  * R: No-cost resistance allele
#'  * B: Detrimental resistance allele
#'
#'
#'
#'  This is the general form for an immunizing reversal drive. If the c_EW and d_EW,
#'  parameters are all 0, then this simplifies to a basic reversal drive.
#'  This drive represents an X-linked IR drive.
#'
#' @param cHW Cutting efficiency of H into W
#' @param cEW Cutting efficiency of E into W
#' @param cEH Cutting efficiency of E into H
#' @param chHW Homing efficiency of H at W
#' @param crHW Resistance generation efficiency of H at W
#' @param ceEW Homing efficiency of E at W
#' @param crEW Resistance generation efficiency of E at W
#' @param ceEH Homing efficiency of E at H
#' @param crEH Resistance efficiency of E at H
#' @param dHW H deposition efficiency against W
#' @param dEW E deposition efficiency against W
#' @param dEH E deposition efficiency against H
#' @param dhHW H deposition homing efficiency against W
#' @param drHW H deposition resistance efficiency against W
#' @param deEW E deposition homing efficiency against W
#' @param drEW E deposition resistance efficiency against W
#' @param deEH E deposition homing efficiency against H
#' @param drEH E deposition resistance efficiency against H
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
cubeImmunizingReversalX <- function(cHW=1.0, cEW=1.0, cEH=1.0, chHW=0, crHW=0,
                                     ceEW=0, crEW=0, ceEH=0, crEH=0, dHW=0, dEW=0,
                                     dEH=0, dhHW=0, drHW=0, deEW=0, drEW=0,
                                     deEH=0, drEH=0, eta=NULL, phi=NULL,
                                     omega=NULL, xiF=NULL, xiM=NULL, s=NULL){

  ## safety checks
  if(any(c(cHW,cEW,cEH,chHW,crHW,ceEW,crEW,ceEH,crEH,dHW,dEW,dEH,dhHW,drHW,deEW,drEW,deEH,drEH)>1) ||
     any(c(cHW,cEW,cEH,chHW,crHW,ceEW,crEW,ceEH,crEH,dHW,dEW,dEH,dhHW,drHW,deEW,drEW,deEH,drEH)<0)){
    stop("Parameters are rates.\n0 <= x <= 1")
  }


  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  gtype <- c('WW', 'WH', 'WE', 'WY', 'WR', 'WB', 'HH', 'HE', 'HY', 'HR', 'HB',
             'EE', 'EY', 'ER', 'EB', 'YR', 'YB', 'RR', 'RB', 'BB')
  size <- length(gtype) #because I use it several times
  tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix


  ## fill tMatrix with probabilities
  #( 'WW', 'WH', 'WE', 'WY', 'WR', 'WB',
  #  'HH', 'HE', 'HY', 'HR', 'HB', 'EE',
  #  'EY', 'ER', 'EB', 'YR', 'YB',
  #  'RR', 'RB', 'BB')
  tMatrix['WW','WY', c('WW', 'WY')] <- c( 1, 1)/2
  tMatrix['WW','HY', c('WH', 'WY')] <- c( 1, 1)/2
  tMatrix['WW','EY', c('WE', 'WY')] <- c( 1, 1)/2
  tMatrix['WW','YR', c('WY', 'WR')] <- c( 1, 1)/2
  tMatrix['WW','YB', c('WY', 'WB')] <- c( 1, 1)/2

  tMatrix['WH','WY', ] <- c( (1-cHW)*(1-dHW), (1+cHW*chHW)*(1-dHW), 0, 1-cHW, (1-cHW)*dHW*drHW + (cHW*(1-chHW)*crHW)*(1-dHW), (1-cHW)*dHW*(1-drHW) + (cHW*(1-chHW)*(1-crHW))*(1-dHW),
                             (1+cHW*chHW)*dHW*dhHW, 0, 1+cHW*chHW, (1+cHW*chHW)*dHW*(1-dhHW)*drHW, (1+cHW*chHW)*dHW*(1-dhHW)*(1-drHW), 0,
                             0, 0, 0, cHW*(1-chHW)*crHW, cHW*(1-chHW)*(1-crHW),
                             (cHW*(1-chHW)*crHW)*dHW*drHW, (cHW*(1-chHW)*crHW)*dHW*(1-drHW) + (cHW*(1-chHW)*(1-crHW))*dHW*drHW, (cHW*(1-chHW)*(1-crHW))*dHW*(1-drHW))/4
  tMatrix['WH','HY', c('WH', 'HH', 'HR', 'HB',
                       'WY', 'HY', 'YR', 'YB')] <- c( (1-cHW), 1+cHW*chHW, cHW*(1-chHW)*crHW, cHW*(1-chHW)*(1-crHW),
                                                      (1-cHW), 1+cHW*chHW, cHW*(1-chHW)*crHW, cHW*(1-chHW)*(1-crHW))/4
  tMatrix['WH','EY', c('WE', 'HE', 'ER', 'EB',
                       'WY', 'HY', 'YR', 'YB')] <- c( (1-cHW), 1+cHW*chHW, cHW*(1-chHW)*crHW, cHW*(1-chHW)*(1-crHW),
                                                      (1-cHW), 1+cHW*chHW, cHW*(1-chHW)*crHW, cHW*(1-chHW)*(1-crHW))/4
  tMatrix['WH','YR', c('WY', 'HY', 'YR', 'YB',
                       'WR', 'HR', 'RR', 'RB')] <- c( (1-cHW), 1+cHW*chHW, cHW*(1-chHW)*crHW, cHW*(1-chHW)*(1-crHW),
                                                      (1-cHW), 1+cHW*chHW, cHW*(1-chHW)*crHW, cHW*(1-chHW)*(1-crHW))/4
  tMatrix['WH','YB', c('WY', 'HY', 'YR', 'YB',
                       'WB', 'HB', 'RB', 'BB')] <- c( (1-cHW), 1+cHW*chHW, cHW*(1-chHW)*crHW, cHW*(1-chHW)*(1-crHW),
                                                      (1-cHW), 1+cHW*chHW, cHW*(1-chHW)*crHW, cHW*(1-chHW)*(1-crHW))/4

  tMatrix['WE','WY', ] <- c( (1-cEW)*(1-dEW), 0, (1+cEW*ceEW)*(1-dEW), 1-cEW, (1-cEW)*dEW*drEW + (cEW*(1-ceEW)*crEW)*(1-dEW), (1-cEW)*dEW*(1-drEW) + (cEW*(1-ceEW)*(1-crEW))*(1-dEW),
                             0, 0, 0, 0, 0, (1+cEW*ceEW)*dEW*deEW,
                             1+cEW*ceEW, (1+cEW*ceEW)*dEW*(1-deEW)*drEW, (1+cEW*ceEW)*dEW*(1-deEW)*(1-drEW), cEW*(1-ceEW)*crEW, cEW*(1-ceEW)*(1-crEW),
                             (cEW*(1-ceEW)*crEW)*dEW*drEW, (cEW*(1-ceEW)*crEW)*dEW*(1-drEW) + (cEW*(1-ceEW)*(1-crEW))*dEW*drEW, (cEW*(1-ceEW)*(1-crEW))*dEW*(1-drEW))/4
  tMatrix['WE','HY', ] <- c( 0, (1-cEW)*(1-dEH), 0, 1-cEW, (1-cEW)*dEH*drEH, (1-cEW)*dEH*(1-drEH),
                             0, (1+cEW*ceEW)*(1-dEH), 0, (cEW*(1-ceEW)*crEW)*(1-dEH), (cEW*(1-ceEW)*(1-crEW))*(1-dEH), (1+cEW*ceEW)*dEH*deEH,
                             1+cEW*ceEW, (1+cEW*ceEW)*dEH*(1-deEH)*drEH, (1+cEW*ceEW)*dEH*(1-deEH)*(1-drEH), cEW*(1-ceEW)*crEW, cEW*(1-ceEW)*(1-crEW),
                             (cEW*(1-ceEW)*crEW)*dEH*drEH, (cEW*(1-ceEW)*crEW)*dEH*(1-drEH) + (cEW*(1-ceEW)*(1-crEW))*dEH*drEH,(cEW*(1-ceEW)*(1-crEW))*dEH*(1-drEH))/4
  tMatrix['WE','EY', c('WE', 'EE', 'ER', 'EB',
                       'WY', 'EY', 'YR', 'YB')] <- c( (1-cEW), 1+cEW*ceEW, cEW*(1-ceEW)*crEW, cEW*(1-ceEW)*(1-crEW),
                                                      (1-cEW), 1+cEW*ceEW, cEW*(1-ceEW)*crEW, cEW*(1-ceEW)*(1-crEW))/4
  tMatrix['WE','YR', c('WY', 'EY', 'YR', 'YB',
                       'WR', 'ER', 'RR', 'RB')] <- c( (1-cEW), 1+cEW*ceEW, cEW*(1-ceEW)*crEW, cEW*(1-ceEW)*(1-crEW),
                                                      (1-cEW), 1+cEW*ceEW, cEW*(1-ceEW)*crEW, cEW*(1-ceEW)*(1-crEW))/4
  tMatrix['WE','YB', c('WY', 'EY', 'YR', 'YB',
                       'WB', 'EB', 'RB', 'BB')] <- c( (1-cEW), 1+cEW*ceEW, cEW*(1-ceEW)*crEW, cEW*(1-ceEW)*(1-crEW),
                                                      (1-cEW), 1+cEW*ceEW, cEW*(1-ceEW)*crEW, cEW*(1-ceEW)*(1-crEW))/4

  tMatrix['WR','WY', c('WW', 'WY', 'WR', 'YR')] <- c( 1, 1, 1, 1)/4
  tMatrix['WR','HY', c('WH', 'WY', 'HR', 'YR')] <- c( 1, 1, 1, 1)/4
  tMatrix['WR','EY', c('WE', 'WY', 'ER', 'YR')] <- c( 1, 1, 1, 1)/4
  tMatrix['WR','YR', c('WY', 'WR', 'YR', 'RR')] <- c( 1, 1, 1, 1)/4
  tMatrix['WR','YB', c('WY', 'WB', 'YR', 'RB')] <- c( 1, 1, 1, 1)/4

  tMatrix['WB','WY', c('WW', 'WY', 'WB', 'YB')] <- c( 1, 1, 1, 1)/4
  tMatrix['WB','HY', c('WH', 'WY', 'HB', 'YB')] <- c( 1, 1, 1, 1)/4
  tMatrix['WB','EY', c('WE', 'WY', 'EB', 'YB')] <- c( 1, 1, 1, 1)/4
  tMatrix['WB','YR', c('WY', 'WR', 'YB', 'RB')] <- c( 1, 1, 1, 1)/4
  tMatrix['WB','YB', c('WY', 'WB', 'YB', 'BB')] <- c( 1, 1, 1, 1)/4

  tMatrix['HH','WY', c('WH', 'HH', 'HR', 'HB', 'HY')] <- c( 1-dHW, dHW*dhHW, dHW*(1-dhHW)*drHW, dHW*(1-dhHW)*(1-drHW), 1)/2
  tMatrix['HH','HY', c('HH', 'HY')] <- c( 1, 1)/2
  tMatrix['HH','EY', c('HE', 'HY')] <- c( 1, 1)/2
  tMatrix['HH','YR', c('HY', 'HR')] <- c( 1, 1)/2
  tMatrix['HH','YB', c('HY', 'HB')] <- c( 1, 1)/2

  # There is an assumption of no H deposition into the W allele here.
  #  This could be wrong if the E allele is terrible at cutting/homing/breaking/everytyihng
  tMatrix['HE','WY', ] <- c(0, (1-cEH)*(1-dEW), (1+cEH*ceEH)*(1-dEW), 0, (cEH*(1-ceEH)*crEH)*(1-dEW), (cEH*(1-ceEH)*(1-crEH))*(1-dEW),
                            0, 0, 1-cEH, (1-cEH)*dEW*drEW, (1-cEH)*dEW*(1-drEW), (1+cEH*ceEH)*dEW*deEW,
                            1+cEH*ceEH, (1+cEH*ceEH)*dEW*(1-deEW)*drEW, (1+cEH*ceEH)*dEW*(1-deEW)*(1-drEW), cEH*(1-ceEH)*crEH, cEH*(1-ceEH)*(1-crEH),
                            (cEH*(1-ceEH)*crEH)*dEW*drEW, (cEH*(1-ceEH)*crEH)*dEW*(1-drEW) + (cEH*(1-ceEH)*(1-crEH))*dEW*drEW, (cEH*(1-ceEH)*(1-crEH))*dEW*(1-drEW))/4
  tMatrix['HE','HY',] <- c( 0, 0, 0, 0, 0, 0,
                           (1-cEH)*(1-dEH), (1+cEH*ceEH)*(1-dEH), 1-cEH, (1-cEH)*dEH*drEH + (cEH*(1-ceEH)*crEH)*(1-dEH), (1-cEH)*dEH*(1-drEH) + (cEH*(1-ceEH)*(1-crEH))*(1-dEH), (1+cEH*ceEH)*dEH*deEH,
                           1+cEH*ceEH, (1+cEH*ceEH)*dEH*(1-deEH)*drEH, (1+cEH*ceEH)*dEH*(1-deEH)*(1-drEH), cEH*(1-ceEH)*crEH, cEH*(1-ceEH)*(1-crEH),
                           (cEH*(1-ceEH)*crEH)*dEH*drEH, (cEH*(1-ceEH)*crEH)*dEH*(1-drEH) + (cEH*(1-ceEH)*(1-crEH))*dEH*drEH, (cEH*(1-ceEH)*(1-crEH))*dEH*(1-drEH))/4
  tMatrix['HE','EY', c('HE', 'EE', 'ER', 'EB',
                       'HY', 'EY', 'YR', 'YB')] <- c( 1-cEH, 1+cEH*ceEH, cEH*(1-ceEH)*crEH, cEH*(1-ceEH)*(1-crEH),
                                                      1-cEH, 1+cEH*ceEH, cEH*(1-ceEH)*crEH, cEH*(1-ceEH)*(1-crEH))/4
  tMatrix['HE','YR', c('HY', 'EY', 'YR', 'YB',
                       'HR', 'ER', 'RR', 'RB')] <- c( 1-cEH, 1+cEH*ceEH, cEH*(1-ceEH)*crEH, cEH*(1-ceEH)*(1-crEH),
                                                      1-cEH, 1+cEH*ceEH, cEH*(1-ceEH)*crEH, cEH*(1-ceEH)*(1-crEH))/4
  tMatrix['HE','YB', c('HY', 'EY', 'YR', 'YB',
                       'HB', 'EB', 'RB', 'BB')] <- c( 1-cEH, 1+cEH*ceEH, cEH*(1-ceEH)*crEH, cEH*(1-ceEH)*(1-crEH),
                                                      1-cEH, 1+cEH*ceEH, cEH*(1-ceEH)*crEH, cEH*(1-ceEH)*(1-crEH))/4

  tMatrix['HR','WY', c('WH', 'HH', 'HR', 'HB', 'HY',
                       'WR', 'RR', 'RB', 'YR')] <- c( 1-dHW, dHW*dhHW, dHW*(1-dhHW)*drHW, dHW*(1-dhHW)*(1-drHW), 1,
                                                      1-dHW, dHW*drHW, dHW*(1-drHW), 1)/4
  tMatrix['HR','HY', c('HH', 'HY', 'HR', 'YR')] <- c( 1, 1, 1, 1)/4
  tMatrix['HR','EY', c('HE', 'HY', 'ER', 'YR')] <- c( 1, 1, 1, 1)/4
  tMatrix['HR','YR', c('HY', 'HR', 'YR', 'RR')] <- c( 1, 1, 1, 1)/4
  tMatrix['HR','YB', c('HY', 'HB', 'YR', 'RB')] <- c( 1, 1, 1, 1)/4

  tMatrix['HB','WY', c('WH', 'HH', 'HR', 'HB', 'HY',
                       'WB', 'RB', 'BB', 'YB')] <- c( 1-dHW, dHW*dhHW, dHW*(1-dhHW)*drHW, dHW*(1-dhHW)*(1-drHW), 1,
                                                      1-dHW, dHW*drHW, dHW*(1-drHW), 1)/4
  tMatrix['HB','HY', c('HH', 'HY', 'HB', 'YB')] <- c( 1, 1, 1, 1)/4
  tMatrix['HB','EY', c('HE', 'HY', 'EB', 'YB')] <- c( 1, 1, 1, 1)/4
  tMatrix['HB','YR', c('HY', 'HR', 'YB', 'RB')] <- c( 1, 1, 1, 1)/4
  tMatrix['HB','YB', c('HY', 'HB', 'YB', 'BB')] <- c( 1, 1, 1, 1)/4

  tMatrix['EE','WY', c('WE', 'EE', 'ER', 'EB', 'EY')] <- c( 1-dEW, dEW*deEW, dEW*(1-deEW)*drEW, dEW*(1-deEW)*(1-drEW), 1)/2
  tMatrix['EE','HY', c('HE', 'EE', 'ER', 'EB', 'EY')] <- c( 1-dEH, dEH*deEH, dEH*(1-deEH)*drEH, dEH*(1-deEH)*(1-drEH), 1)/2
  tMatrix['EE','EY', c('EE', 'EY')] <- c( 1, 1)/2
  tMatrix['EE','YR', c('EY', 'ER')] <- c( 1, 1)/2
  tMatrix['EE','YB', c('EY', 'EB')] <- c( 1, 1)/2

  tMatrix['ER','WY', c('WE', 'EE', 'ER', 'EB', 'EY',
                       'WR', 'RR', 'RB', 'YR')] <- c( 1-dEW, dEW*deEW, dEW*(1-deEW)*drEW, dEW*(1-deEW)*(1-drEW), 1,
                                                      1-dEW, dEW*drEW, dEW*(1-drEW), 1)/4
  tMatrix['ER','HY', c('HE', 'EE', 'ER', 'EB', 'EY',
                       'HR', 'RR', 'RB', 'YR')] <- c( 1-dEH, dEH*deEH, dEH*(1-deEH)*drEH, dEH*(1-deEH)*(1-drEH), 1,
                                                      1-dEH, dEH*drEH, dEH*(1-drEH), 1)/4
  tMatrix['ER','EY', c('EE', 'EY', 'ER', 'YR')] <- c( 1, 1, 1, 1)/4
  tMatrix['ER','YR', c('EY', 'ER', 'YR', 'RR')] <- c( 1, 1, 1, 1)/4
  tMatrix['ER','YB', c('EY', 'EB', 'YR', 'RB')] <- c( 1, 1, 1, 1)/4

  tMatrix['EB','WY', c('WE', 'EE', 'ER', 'EB', 'EY',
                       'WB', 'RB', 'BB', 'YB')] <- c( 1-dEW, dEW*deEW, dEW*(1-deEW)*drEW, dEW*(1-deEW)*(1-drEW), 1,
                                                      1-dEW, dEW*drEW, dEW*(1-drEW), 1)/4
  tMatrix['EB','HY', c('HE', 'EE', 'ER', 'EB', 'EY',
                       'HB', 'RB', 'BB', 'YB')] <- c( 1-dEH, dEH*deEH, dEH*(1-deEH)*drEH, dEH*(1-deEH)*(1-drEH), 1,
                                                      1-dEH, dEH*drEH, dEH*(1-drEH), 1)/4
  tMatrix['EB','EY', c('EE', 'EY', 'EB', 'YB')] <- c( 1, 1, 1, 1)/4
  tMatrix['EB','YR', c('EY', 'ER', 'YB', 'RB')] <- c( 1, 1, 1, 1)/4
  tMatrix['EB','YB', c('EY', 'EB', 'YB', 'BB')] <- c( 1, 1, 1, 1)/4

  tMatrix['RR','WY', c('WR', 'YR')] <- c( 1, 1)/2
  tMatrix['RR','HY', c('HR', 'YR')] <- c( 1, 1)/2
  tMatrix['RR','EY', c('ER', 'YR')] <- c( 1, 1)/2
  tMatrix['RR','YR', c('YR', 'RR')] <- c( 1, 1)/2
  tMatrix['RR','YB', c('YR', 'RB')] <- c( 1, 1)/2

  tMatrix['RB','WY', c('WR', 'WB', 'YR', 'YB')] <- c( 1, 1, 1, 1)/4
  tMatrix['RB','HY', c('HR', 'HB', 'YR', 'YB')] <- c( 1, 1, 1, 1)/4
  tMatrix['RB','EY', c('ER', 'EB', 'YR', 'YB')] <- c( 1, 1, 1, 1)/4
  tMatrix['RB','YR', c('YR', 'YB', 'RR', 'RB')] <- c( 1, 1, 1, 1)/4
  tMatrix['RB','YB', c('YR', 'YB', 'RB', 'BB')] <- c( 1, 1, 1, 1)/4

  tMatrix['BB','WY', c('WB', 'YB')] <- c( 1, 1)/2
  tMatrix['BB','HY', c('HB', 'YB')] <- c( 1, 1)/2
  tMatrix['BB','EY', c('EB', 'YB')] <- c( 1, 1)/2
  tMatrix['BB','YR', c('YB', 'RB')] <- c( 1, 1)/2
  tMatrix['BB','YB', c('YB', 'BB')] <- c( 1, 1)/2


  ## protection from underflow errors
  tMatrix[tMatrix < .Machine$double.eps] <- 0


  ## initialize viability mask. No mother-specific death.
  viabilityMask <- array(data = 1, dim = c(size,size,size), dimnames = list(gtype, gtype, gtype))

  ## genotype-specific modifiers
  modifiers = cubeModifiers(gtype, eta = eta, phi = phi, omega = omega, xiF = xiF, xiM = xiM, s = s)
  modifiers$phi[] <- 1
  modifiers$phi[c('WY', 'HY', 'EY', 'YR', 'YB')] <- 0 #These genotypes can only be male. The rest only female.


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
    releaseType = "HH"
  ))

}
