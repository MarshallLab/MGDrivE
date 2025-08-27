###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   Immunizing Reversal
#   Héctor Sanchez, Jared Bennett, Sean Wu, John Marshall
#   jared_bennett@berkeley.edu
#   August 2017
#   December 2018
#    Update to reflect cutting, homing, resistance generation rates
#    Added female deposition. Made the assumption that if a cross is similar to
#    HExWW, there is only deposition from the E construct.
#    This breaks down if there is heavy deposition from the H construct and the
#    E construct sucks at the same time.
#
###############################################################################

#' Inheritance Cube: Immunizing Reversal/Basic Reversal
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
#'  This is the general form for an immunizing reversal drive. If the *EW* terms are 0,
#'  then this simplifies to a basic reversal drive.This drive handles different
#'  male and female homing rates, and female deposition from each allele,
#'  signifying differential expression from an autosome.
#'
#' @param cHWM Cutting efficiency of H into W in males
#' @param cHWF Cutting efficiency of H into W in females
#' @param cEWM Cutting efficiency of E into W in males
#' @param cEWF Cutting efficiency of E into W in females
#' @param cEHM Cutting efficiency of E into H in males
#' @param cEHF Cutting efficiency of E into H in females
#' @param chHWM Homing efficiency of H into W in males
#' @param chHWF Homing efficiency of H into W in females
#' @param crHWM Resistance efficiency of H into W in males
#' @param crHWF Resistance efficiency of H into W in females
#' @param ceEWM Homing efficiency of E into W in males
#' @param ceEWF Homing efficiency of E into W in females
#' @param crEWM Resistance efficiency of E into W in males
#' @param crEWF Resistance efficiency of E into W in females
#' @param ceEHM Homing efficiency of E into H in males
#' @param ceEHF Homing efficiency of E into H in females
#' @param crEHM Resistance efficiency of E into H in males
#' @param crEHF Resistance efficiency of E into H in females
#' @param dHW Deposition cutting efficiency of H into W
#' @param dEW Deposition cutting efficiency of E in to W
#' @param dEH Deposition cutting efficiency of E into H
#' @param dhHW Deposition homing efficiency of H into W
#' @param drHW Deposition resistance efficiency of H into W
#' @param deEW Deposition homing efficiency of E into W
#' @param drEW Deposition resistance efficiency of E into W
#' @param deEH Deposition homing efficiency of E into H
#' @param drEH Deposition resistance efficiency of E into H
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
cubeImmunizingReversalMF <- function(cHWM=1.0, cHWF=1.0, cEWM=1.0, cEWF=1.0,
                                      cEHM=1.0, cEHF=1.0, chHWM=0, chHWF=0,
                                      crHWM=0, crHWF=0, ceEWM=0, ceEWF=0,
                                      crEWM=0, crEWF=0, ceEHM=0, ceEHF=0,
                                      crEHM=0, crEHF=0, dHW=0, dEW=0, dEH=0,
                                      dhHW=0, drHW=0, deEW=0, drEW=0, deEH=0,
                                      drEH=0, eta=NULL, phi=NULL,omega=NULL,
                                      xiF=NULL, xiM=NULL, s=NULL){

  #This cube has female and male specific homing.
  # the general form is a full immunizing reversal drive
  # setting eEWM=eEWF=0, and all other *EW*=0, this reduces to the basic reversal drive.


  ## safety checks
  if(any(c(cHWM,cHWF,cEWM,cEWF,cEHM,cEHF,chHWM,chHWF,crHWM,crHWF,ceEWM,ceEWF,crEWM,crEWF,ceEHM,ceEHF,crEHM,crEHF,dHW,dEW,dEH,dhHW,drHW,deEW,drEW,deEH,drEH)>1) ||
     any(c(cHWM,cHWF,cEWM,cEWF,cEHM,cEHF,chHWM,chHWF,crHWM,crHWF,ceEWM,ceEWF,crEWM,crEWF,ceEHM,ceEHF,crEHM,crEHF,dHW,dEW,dEH,dhHW,drHW,deEW,drEW,deEH,drEH)<0)){
    stop("Parameters are rates.\n0 <= X <= 1")
  }

  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  gtype <- c('WW', 'WH', 'WE', 'WR', 'WB', 'HH', 'HE', 'HR', 'HB', 'EE', 'ER', 'EB', 'RR', 'RB', 'BB')
  size <- length(gtype) #because I use it several times
  tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix


  ## fill tMatrix with symmetric probabilities
  #( 'WW', 'WH', 'WE', 'WR', 'WB',
  #  'HH', 'HE', 'HR', 'HB', 'EE',
  #  'ER', 'EB', 'RR', 'RB', 'BB')
  tMatrix['WW','WW', 'WW'] <- 1

  tMatrix['WR','WW', c('WW', 'WR')] <- c( 1, 1)/2
  tMatrix['WR','WR', c('WW', 'WR', 'RR')] <- c( 1/2, 1, 1/2)/2

  tMatrix['WB','WW', c('WW', 'WB')] <- c( 1, 1)/2
  tMatrix['WB','WR', c('WW', 'WR', 'WB', 'RB')] <- c( 1, 1, 1, 1)/4
  tMatrix['WB','WB', c('WW', 'WB', 'BB')] <- c( 1/2, 1, 1/2)/2

  tMatrix['HH','HH', 'HH'] <- 1

  tMatrix['HR','HH', c('HH', 'HR')] <- c( 1, 1)/2
  tMatrix['HR','HR', c('HH', 'HR', 'RR')] <- c( 1/4, 1/2, 1/4)

  tMatrix['HB','HH', c('HH', 'HB')] <- c( 1, 1)/2
  tMatrix['HB','HR', c('HH', 'HR', 'HB', 'RB')] <- c( 1, 1, 1, 1)/4
  tMatrix['HB','HB', c('HH', 'HB', 'BB')] <- c( 1/4, 1/2, 1/4)

  tMatrix['EE','EE', 'EE'] <- 1

  tMatrix['ER','EE', c('EE', 'ER')] <- c( 1, 1)/2
  tMatrix['ER','ER', c('EE', 'ER', 'RR')] <- c( 1/4, 1/2, 1/4)

  tMatrix['EB','EE', c('EE', 'EB')] <- c( 1, 1)/2
  tMatrix['EB','ER', c('EE', 'ER', 'EB', 'RB')] <- c( 1, 1, 1, 1)/4
  tMatrix['EB','EB', c('EE', 'EB', 'BB')] <- c( 1/4, 1/2, 1/4)

  tMatrix['RR','WW', 'WR'] <- 1
  tMatrix['RR','WR', c('WR', 'RR')] <- c( 1, 1)/2
  tMatrix['RR','WB', c('WR', 'RB')] <- c( 1, 1)/2
  tMatrix['RR','HH', 'HR'] <- 1
  tMatrix['RR','HR', c('HR', 'RR')] <- c( 1, 1)/2
  tMatrix['RR','HB', c('HR', 'RB')] <- c( 1, 1)/2
  tMatrix['RR','EE', c('ER')] <- 1
  tMatrix['RR','ER', c('ER', 'RR')] <- c( 1, 1)/2
  tMatrix['RR','EB', c('ER', 'RB')] <- c( 1, 1)/2
  tMatrix['RR','RR', 'RR'] <- 1

  tMatrix['RB','WW', c('WR', 'WB')] <- c( 1, 1)/2
  tMatrix['RB','WR', c('WR', 'WB', 'RR', 'RB')] <- c( 1, 1, 1, 1)/4
  tMatrix['RB','WB', c('WR', 'WB', 'RB', 'BB')] <- c( 1, 1, 1, 1)/4
  tMatrix['RB','HH', c('HR', 'HB')] <- c( 1, 1)/2
  tMatrix['RB','HR', c('HR', 'HB', 'RR', 'RB')] <- c( 1, 1, 1, 1)/4
  tMatrix['RB','HB', c('HR', 'HB', 'RB', 'BB')] <- c( 1, 1, 1, 1)/4
  tMatrix['RB','EE', c('ER', 'EB')] <- c( 1, 1)/2
  tMatrix['RB','ER', c('ER', 'EB', 'RR', 'RB')] <- c( 1, 1, 1, 1)/4
  tMatrix['RB','EB', c('ER', 'EB', 'RB', 'BB')] <- c( 1, 1, 1, 1)/4
  tMatrix['RB','RR', c('RR', 'RB')] <- c( 1, 1)/2
  tMatrix['RB','RB', c('RR', 'RB', 'BB')] <- c( 1/4, 1/2, 1/4)

  tMatrix['BB','WW', 'WB'] <- 1
  tMatrix['BB','WR', c('WB', 'RB')] <- c( 1, 1)/2
  tMatrix['BB','WB', c('WB', 'BB')] <- c( 1, 1)/2
  tMatrix['BB','HH', 'HB'] <- 1
  tMatrix['BB','HR', c('HB', 'RB')] <- c( 1, 1)/2
  tMatrix['BB','HB', c('HB', 'BB')] <- c( 1, 1)/2
  tMatrix['BB','EE', 'EB'] <- 1
  tMatrix['BB','ER', c('EB', 'RB')] <- c( 1, 1)/2
  tMatrix['BB','EB', c('EB', 'BB')] <- c( 1, 1)/2
  tMatrix['BB','RR', 'RB'] <- 1
  tMatrix['BB','RB', c('RB', 'BB')] <- c( 1, 1)/2
  tMatrix['BB','BB', 'BB'] <- 1

  ## set the other half of the symmetric part of tMatrix
  # Boolean matrix for subsetting, used several times
  boolMat <- upper.tri(x = tMatrix[ , ,1], diag = FALSE)
  # loop over depth, set upper triangle
  for(z in 1:size){tMatrix[ , ,z][boolMat] <- t(tMatrix[ , ,z])[boolMat]}


  ## fill tMatrix with asymmetric probabilities
  ## This section corresponding to deposition only, no homing

  # female
  tMatrix['HH','WW', c('WH', 'HH', 'HR', 'HB')] <- c( 1-dHW, dHW*dhHW, dHW*(1-dhHW)*drHW, dHW*(1-dhHW)*(1-drHW))
  tMatrix['HH','WR', c('WH', 'HH', 'HR', 'HB')] <- c( 1-dHW, dHW*dhHW, dHW*(1-dhHW)*drHW + 1, dHW*(1-dhHW)*(1-drHW))/2
  tMatrix['HH','WB', c('WH', 'HH', 'HR', 'HB')] <- c( 1-dHW, dHW*dhHW, dHW*(1-dhHW)*drHW, dHW*(1-dhHW)*(1-drHW) + 1)/2

  # male
  tMatrix['WW','HH', 'WH'] <- 1
  tMatrix['WR','HH', c('WH', 'HR')] <- c( 1, 1)/2
  tMatrix['WB','HH', c('WH', 'HB')] <- c( 1, 1)/2


  # female
  tMatrix['HR','WW', c('WH', 'HH', 'HR', 'HB',
                       'WR', 'RR', 'RB')] <- c( 1-dHW, dHW*dhHW, dHW*(1-dhHW)*drHW, dHW*(1-dhHW)*(1-drHW),
                                                1-dHW, dHW*drHW, dHW*(1-drHW))/2
  tMatrix['HR','WR', c('WH', 'HH', 'HR', 'HB',
                       'WR', 'RR', 'RB')] <- c( 1-dHW, dHW*dhHW, dHW*(1-dhHW)*drHW + 1, dHW*(1-dhHW)*(1-drHW),
                                                1-dHW, dHW*drHW + 1, dHW*(1-drHW))/4
  tMatrix['HR','WB', c('WH', 'HH', 'HR', 'HB',
                       'WR', 'RR', 'RB')] <- c( 1-dHW, dHW*dhHW, dHW*(1-dhHW)*drHW, dHW*(1-dhHW)*(1-drHW) + 1,
                                                1-dHW, dHW*drHW, dHW*(1-drHW) + 1)/4

  # male
  tMatrix['WW','HR', c('WH', 'WR')] <- c( 1/2, 1/2)
  tMatrix['WR','HR', c('WH', 'WR', 'HR', 'RR')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['WB','HR', c('WH', 'WR', 'HB', 'RB')] <- c( 1/4, 1/4, 1/4, 1/4)


  # Female
  tMatrix['HB','WW', c('WH', 'HH', 'HR', 'HB',
                       'WB', 'RB', 'BB')] <- c( 1-dHW, dHW*dhHW, dHW*(1-dhHW)*drHW, dHW*(1-dhHW)*(1-drHW),
                                                1-dHW, dHW*drHW, dHW*(1-drHW))/2
  tMatrix['HB','WR', c('WH', 'HH', 'HR', 'HB',
                       'WB', 'RB', 'BB')] <- c( 1-dHW, dHW*dhHW, dHW*(1-dhHW)*drHW + 1, dHW*(1-dhHW)*(1-drHW),
                                                1-dHW, dHW*drHW + 1, dHW*(1-drHW))/4
  tMatrix['HB','WB', c('WH', 'HH', 'HR', 'HB',
                       'WB', 'RB', 'BB')] <- c( 1-dHW, dHW*dhHW, dHW*(1-dhHW)*drHW, dHW*(1-dhHW)*(1-drHW) + 1,
                                                1-dHW, dHW*drHW, dHW*(1-drHW) + 1)/4

  # male
  tMatrix['WW','HB', c('WH', 'WB')] <- c( 1/2, 1/2)
  tMatrix['WR','HB', c('WH', 'WB', 'HR', 'RB')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['WB','HB', c('WH', 'WB', 'HB', 'BB')] <- c( 1/4, 1/4, 1/4, 1/4)


  # female
  tMatrix['EE','WW', c('WE', 'EE', 'ER', 'EB')] <- c( 1-dEW, dEW*deEW, dEW*(1-deEW)*drEW, dEW*(1-deEW)*(1-drEW))
  tMatrix['EE','WR', c('WE', 'EE', 'ER', 'EB')] <- c( 1-dEW, dEW*deEW, dEW*(1-deEW)*drEW + 1, dEW*(1-deEW)*(1-drEW))/2
  tMatrix['EE','WB', c('WE', 'EE', 'ER', 'EB')] <- c( 1-dEW, dEW*deEW, dEW*(1-deEW)*drEW, dEW*(1-deEW)*(1-drEW) + 1)/2
  tMatrix['EE','HH', c('HE', 'EE', 'ER', 'EB')] <- c( 1-dEH, dEH*deEH, dEH*(1-deEH)*drEH, dEH*(1-deEH)*(1-drEH))
  tMatrix['EE','HR', c('HE', 'EE', 'ER', 'EB')] <- c( 1-dEH, dEH*deEH, dEH*(1-deEH)*drEH + 1, dEH*(1-deEH)*(1-drEH))/2
  tMatrix['EE','HB', c('HE', 'EE', 'ER', 'EB')] <- c( 1-dEH, dEH*deEH, dEH*(1-deEH)*drEH, dEH*(1-deEH)*(1-drEH) + 1)/2

  # male
  tMatrix['WW','EE', 'WE'] <- 1
  tMatrix['WR','EE', c('WE', 'ER')] <- c( 1, 1)/2
  tMatrix['WB','EE', c('WE', 'EB')] <- c( 1, 1)/2
  tMatrix['HH','EE', 'HE'] <- 1
  tMatrix['HR','EE', c('HE', 'ER')] <- c( 1, 1)/2
  tMatrix['HB','EE', c('HE', 'EB')] <- c( 1, 1)/2


  # female
  tMatrix['ER','WW', c('WE', 'EE', 'ER', 'EB',
                       'WR', 'RR', 'RB')] <- c( 1-dEW, dEW*deEW, dEW*(1-deEW)*drEW, dEW*(1-deEW)*(1-drEW),
                                                1-dEW, dEW*drEW, dEW*(1-drEW))/2
  tMatrix['ER','WR', c('WE', 'EE', 'ER', 'EB',
                       'WR', 'RR', 'RB')] <- c( 1-dEW, dEW*deEW, dEW*(1-deEW)*drEW + 1, dEW*(1-deEW)*(1-drEW),
                                                1-dEW, dEW*drEW + 1, dEW*(1-drEW))/4
  tMatrix['ER','WB', c('WE', 'EE', 'ER', 'EB',
                       'WR', 'RR', 'RB')] <- c( 1-dEW, dEW*deEW, dEW*(1-deEW)*drEW, dEW*(1-deEW)*(1-drEW) + 1,
                                                1-dEW, dEW*drEW, dEW*(1-drEW) + 1)/4
  tMatrix['ER','HH', c('HE', 'EE', 'ER', 'EB',
                       'HR', 'RR', 'RB')] <- c( 1-dEH, dEH*deEH, dEH*(1-deEH)*drEH, dEH*(1-deEH)*(1-drEH),
                                                1-dEH, dEH*drEH, dEH*(1-drEH))/2
  tMatrix['ER','HR', c('HE', 'EE', 'ER', 'EB',
                       'HR', 'RR', 'RB')] <- c( 1-dEH, dEH*deEH, dEH*(1-deEH)*drEH + 1, dEH*(1-deEH)*(1-drEH),
                                                1-dEH, dEH*drEH + 1, dEH*(1-drEH))/4
  tMatrix['ER','HB', c('HE', 'EE', 'ER', 'EB',
                       'HR', 'RR', 'RB')] <- c( 1-dEH, dEH*deEH, dEH*(1-deEH)*drEH, dEH*(1-deEH)*(1-drEH) + 1,
                                                1-dEH, dEH*drEH, dEH*(1-drEH) + 1)/4

  # male
  tMatrix['WW','ER', c('WE', 'WR')] <- c( 1/2, 1/2)
  tMatrix['WR','ER', c('WE', 'WR', 'ER', 'RR')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['WB','ER', c('WE', 'WR', 'EB', 'RB')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['HH','ER', c('HE', 'HR')] <- c( 1/2, 1/2)
  tMatrix['HR','ER', c('HE', 'HR', 'ER', 'RR')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['HB','ER', c('HE', 'HR', 'EB', 'RB')] <- c( 1/4, 1/4, 1/4, 1/4)


  # female
  tMatrix['EB','WW', c('WE', 'EE', 'ER', 'EB',
                       'WB', 'RB', 'BB')] <- c( 1-dEW, dEW*deEW, dEW*(1-deEW)*drEW, dEW*(1-deEW)*(1-drEW),
                                                1-dEW, dEW*drEW, dEW*(1-drEW))/2
  tMatrix['EB','WR', c('WE', 'EE', 'ER', 'EB',
                       'WB', 'RB', 'BB')] <- c( 1-dEW, dEW*deEW, dEW*(1-deEW)*drEW + 1, dEW*(1-deEW)*(1-drEW),
                                                1-dEW, dEW*drEW + 1, dEW*(1-drEW))/4
  tMatrix['EB','WB', c('WE', 'EE', 'ER', 'EB',
                       'WB', 'RB', 'BB')] <- c( 1-dEW, dEW*deEW, dEW*(1-deEW)*drEW, dEW*(1-deEW)*(1-drEW) + 1,
                                                1-dEW, dEW*drEW, dEW*(1-drEW) + 1)/4
  tMatrix['EB','HH', c('HE', 'EE', 'ER', 'EB',
                       'HB', 'RB', 'BB')] <- c( 1-dEH, dEH*deEH, dEH*(1-deEH)*drEH, dEH*(1-deEH)*(1-drEH),
                                                1-dEH, dEH*drEH, dEH*(1-drEH))/2
  tMatrix['EB','HR', c('HE', 'EE', 'ER', 'EB',
                       'HB', 'RB', 'BB')] <- c( 1-dEH, dEH*deEH, dEH*(1-deEH)*drEH + 1, dEH*(1-deEH)*(1-drEH),
                                                1-dEH, dEH*drEH + 1, dEH*(1-drEH))/4
  tMatrix['EB','HB', c('HE', 'EE', 'ER', 'EB',
                       'HB', 'RB', 'BB')] <- c( 1-dEH, dEH*deEH, dEH*(1-deEH)*drEH, dEH*(1-deEH)*(1-drEH) + 1,
                                                1-dEH, dEH*drEH, dEH*(1-drEH) + 1)/4

  # male
  tMatrix['WW','EB', c('WE', 'WB')] <- c( 1/2, 1/2)
  tMatrix['WR','EB', c('WE', 'WB', 'ER', 'RB')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['WB','EB', c('WE', 'WB', 'EB', 'BB')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['HH','EB', c('HE', 'HB')] <- c( 1/2, 1/2)
  tMatrix['HR','EB', c('HE', 'HB', 'ER', 'RB')] <- c( 1/4, 1/4, 1/4, 1/4)
  tMatrix['HB','EB', c('HE', 'HB', 'EB', 'BB')] <- c( 1/4, 1/4, 1/4, 1/4)


  ## These correspond to homing EFFICIENCIES THAT DIFFER IN MALES AND FEMALES

  # females
  tMatrix['WH','WW', c('WW', 'WH', 'WR', 'WB',
                       'HH', 'HR', 'HB',
                       'RR', 'RB', 'BB')] <- c( (1-cHWF)*(1-dHW), (1+cHWF*chHWF)*(1-dHW), (1-cHWF)*dHW*drHW + (cHWF*(1-chHWF)*crHWF)*(1-dHW), (1-cHWF)*dHW*(1-drHW) + (cHWF*(1-chHWF)*(1-crHWF))*(1-dHW),
                                                (1+cHWF*chHWF)*dHW*dhHW, (1+cHWF*chHWF)*dHW*(1-dhHW)*drHW, (1+cHWF*chHWF)*dHW*(1-dhHW)*(1-drHW),
                                                (cHWF*(1-chHWF)*crHWF)*dHW*drHW, (cHWF*(1-chHWF)*crHWF)*dHW*(1-drHW) + (cHWF*(1-chHWF)*(1-crHWF))*dHW*drHW, (cHWF*(1-chHWF)*(1-crHWF))*dHW*(1-drHW))/2
  tMatrix['WE','WW', c('WW', 'WE', 'WR', 'WB',
                       'EE', 'ER', 'EB',
                       'RR', 'RB', 'BB')] <- c( (1-cEWF)*(1-dEW), (1+cEWF*ceEWF)*(1-dEW), (1-cEWF)*dEW*drEW + (cEWF*(1-ceEWF)*crEWF)*(1-dEW), (1-cEWF)*dEW*(1-drEW) + (cEWF*(1-ceEWF)*(1-crEWF))*(1-dEW),
                                                (1+cEWF*ceEWF)*dEW*deEW, (1+cEWF*ceEWF)*dEW*(1-deEW)*drEW, (1+cEWF*ceEWF)*dEW*(1-deEW)*(1-drEW),
                                                (cEWF*(1-ceEWF)*crEWF)*dEW*drEW, (cEWF*(1-ceEWF)*crEWF)*dEW*(1-drEW) + (cEWF*(1-ceEWF)*(1-crEWF))*dEW*drEW, (cEWF*(1-ceEWF)*(1-crEWF))*dEW*(1-drEW))/2
  tMatrix['HE','WW', c('WH', 'WE', 'WR', 'WB',
                       'EE', 'ER', 'EB', 'HR', 'HB',
                       'RR', 'RB', 'BB')] <- c( (1-cEHF)*(1-dEW), (1+cEHF*ceEHF)*(1-dEW), (cEHF*(1-ceEHF)*crEHF)*(1-dEW), (cEHF*(1-ceEHF)*(1-crEHF))*(1-dEW),
                                                (1+cEHF*ceEHF)*dEW*deEW, (1+cEHF*ceEHF)*dEW*(1-deEW)*drEW, (1+cEHF*ceEHF)*dEW*(1-deEW)*(1-drEW), (1-cEHF)*dEW*drEW, (1-cEHF)*dEW*(1-drEW),
                                                (cEHF*(1-ceEHF)*crEHF)*dEW*drEW, (cEHF*(1-ceEHF)*crEHF)*dEW*(1-drEW) + (cEHF*(1-ceEHF)*(1-crEHF))*dEW*drEW, (cEHF*(1-ceEHF)*(1-crEHF))*dEW*(1-drEW))/2

  # males
  tMatrix['WW','WH', c('WW', 'WH', 'WR', 'WB')] <- c( 1-cHWM, 1+cHWM*chHWM, cHWM*(1-chHWM)*crHWM, cHWM*(1-chHWM)*(1-crHWM))/2
  tMatrix['WW','WE', c('WW', 'WE', 'WR', 'WB')] <- c( 1-cEWM, 1+cEWM*ceEWM, cEWM*(1-ceEWM)*crEWM, cEWM*(1-ceEWM)*(1-crEWM))/2
  tMatrix['WW','HE', c('WH', 'WE', 'WR', 'WB')] <- c( 1-cEHM, 1+cEHM*ceEHM, cEHM*(1-ceEHM)*crEHM, cEHM*(1-ceEHM)*(1-crEHM))/2


  # both
  tMatrix['WH','WH', ] <- c( (1-cHWF)*(1-cHWM)*(1-dHW), (1+cHWF*chHWF)*(1-cHWM)*(1-dHW) + (1-cHWF)*(1+cHWM*chHWM), 0, (1-cHWF)*(1-cHWM)*dHW*drHW + (cHWF*(1-chHWF)*crHWF)*(1-cHWM)*(1-dHW) + (1-cHWF)*(cHWM*(1-chHWM)*crHWM),
                             (1-cHWF)*(1-cHWM)*dHW*(1-drHW) + (cHWF*(1-chHWF)*(1-crHWF))*(1-cHWM)*(1-dHW) + (1-cHWF)*(cHWM*(1-chHWM)*(1-crHWM)), (1+cHWF*chHWF)*(1-cHWM)*dHW*dhHW + (1+cHWF*chHWF)*(1+cHWM*chHWM),
                             0, (1+cHWF*chHWF)*(1-cHWM)*dHW*(1-dhHW)*drHW + (cHWF*(1-chHWF)*crHWF)*(1+cHWM*chHWM) + (1+cHWF*chHWF)*(cHWM*(1-chHWM)*crHWM),
                             (1+cHWF*chHWF)*(1-cHWM)*dHW*(1-dhHW)*(1-drHW) + (cHWF*(1-chHWF)*(1-crHWF))*(1+cHWM*chHWM) + (1+cHWF*chHWF)*(cHWM*(1-chHWM)*(1-crHWM)),
                             0, 0, 0, (cHWF*(1-chHWF)*crHWF)*(1-cHWM)*dHW*drHW + (cHWF*(1-chHWF)*crHWF)*(cHWM*(1-chHWM)*crHWM),
                             (cHWF*(1-chHWF)*crHWF)*(1-cHWM)*dHW*(1-drHW) + (cHWF*(1-chHWF)*(1-crHWF))*(1-cHWM)*dHW*drHW + (cHWF*(1-chHWF)*(1-crHWF))*(cHWM*(1-chHWM)*crHWM) + (cHWF*(1-chHWF)*crHWF)*(cHWM*(1-chHWM)*(1-crHWM)),
                             (cHWF*(1-chHWF)*(1-crHWF))*(1-cHWM)*dHW*(1-drHW) + (cHWF*(1-chHWF)*(1-crHWF))*(cHWM*(1-chHWM)*(1-crHWM)))/4

  # females
  tMatrix['WH','WE', ] <- c( (1-cHWF)*(1-cEWM)*(1-dHW), (1+cHWF*chHWF)*(1-cEWM)*(1-dHW), (1-cHWF)*(1+cEWM*ceEWM),
                             (1-cHWF)*(1-cEWM)*dHW*drHW + (1-cHWF)*cEWM*(1-ceEWM)*crEWM + (cHWF*(1-chHWF)*crHWF)*(1-cEWM)*(1-dHW), (1-cHWF)*(1-cEWM)*dHW*(1-drHW) + (1-cHWF)*cEWM*(1-ceEWM)*(1-crEWM) + cHWF*(1-chHWF)*(1-crHWF)*(1-cEWM)*(1-dHW),
                             (1+cHWF*chHWF)*(1-cEWM)*dHW*dhHW, (1+cHWF*chHWF)*(1+cEWM*ceEWM), (1+cHWF*chHWF)*(1-cEWM)*dHW*(1-dhHW)*drHW + (1+cHWF*chHWF)*cEWM*(1-ceEWM)*crEWM,
                             (1+cHWF*chHWF)*(1-cEWM)*dHW*(1-dhHW)*(1-drHW) + (1+cHWF*chHWF)*cEWM*(1-ceEWM)*(1-crEWM), 0,
                             (cHWF*(1-chHWF)*crHWF)*(1+cEWM*ceEWM), cHWF*(1-chHWF)*(1-crHWF)*(1+cEWM*ceEWM), (cHWF*(1-chHWF)*crHWF)*(1-cEWM)*dHW*drHW + (cHWF*(1-chHWF)*crHWF)*cEWM*(1-ceEWM)*crEWM,
                             (cHWF*(1-chHWF)*crHWF)*(1-cEWM)*dHW*(1-drHW) + (cHWF*(1-chHWF)*crHWF)*cEWM*(1-ceEWM)*(1-crEWM) + cHWF*(1-chHWF)*(1-crHWF)*(1-cEWM)*dHW*drHW + cHWF*(1-chHWF)*(1-crHWF)*cEWM*(1-ceEWM)*crEWM,
                             cHWF*(1-chHWF)*(1-crHWF)*(1-cEWM)*dHW*(1-drHW) + cHWF*(1-chHWF)*(1-crHWF)*cEWM*(1-ceEWM)*(1-crEWM))/4
  tMatrix['WH','WR', c('WW', 'WH', 'WR', 'WB',
                       'HH', 'HR', 'HB',
                       'RR', 'RB', 'BB')] <- c( (1-cHWF)*(1-dHW), (1+cHWF*chHWF)*(1-dHW), (1-cHWF)*dHW*drHW + (cHWF*(1-chHWF)*crHWF)*(1-dHW) + (1-cHWF), (1-cHWF)*dHW*(1-drHW) + (cHWF*(1-chHWF)*(1-crHWF))*(1-dHW),
                                                (1+cHWF*chHWF)*dHW*dhHW, (1+cHWF*chHWF)*dHW*(1-dhHW)*drHW + 1+cHWF*chHWF, (1+cHWF*chHWF)*dHW*(1-dhHW)*(1-drHW),
                                                (cHWF*(1-chHWF)*crHWF)*dHW*drHW + cHWF*(1-chHWF)*crHWF, (cHWF*(1-chHWF)*crHWF)*dHW*(1-drHW) + (cHWF*(1-chHWF)*(1-crHWF))*dHW*drHW + cHWF*(1-chHWF)*(1-crHWF),
                                                (cHWF*(1-chHWF)*(1-crHWF))*dHW*(1-drHW))/4
  tMatrix['WH','WB', c('WW', 'WH', 'WR', 'WB',
                       'HH', 'HR', 'HB',
                       'RR', 'RB', 'BB')] <- c( (1-cHWF)*(1-dHW), (1+cHWF*chHWF)*(1-dHW), (1-cHWF)*dHW*drHW + (cHWF*(1-chHWF)*crHWF)*(1-dHW), (1-cHWF)*dHW*(1-drHW) + (cHWF*(1-chHWF)*(1-crHWF))*(1-dHW) + (1-cHWF),
                                                (1+cHWF*chHWF)*dHW*dhHW, (1+cHWF*chHWF)*dHW*(1-dhHW)*drHW, (1+cHWF*chHWF)*dHW*(1-dhHW)*(1-drHW) + 1+cHWF*chHWF,
                                                (cHWF*(1-chHWF)*crHWF)*dHW*drHW, (cHWF*(1-chHWF)*crHWF)*dHW*(1-drHW) + (cHWF*(1-chHWF)*(1-crHWF))*dHW*drHW + cHWF*(1-chHWF)*crHWF,
                                                (cHWF*(1-chHWF)*(1-crHWF))*dHW*(1-drHW) + cHWF*(1-chHWF)*(1-crHWF))/4
  tMatrix['WH','HH', c('WH', 'HH', 'HR', 'HB')] <- c( (1-cHWF), (1+cHWF*chHWF), cHWF*(1-chHWF)*crHWF, cHWF*(1-chHWF)*(1-crHWF))/2
  tMatrix['WH','HE', ] <- c( 0, (1-cHWF)*(1-cEHM), (1-cHWF)*(1+cEHM*ceEHM), (1-cHWF)*cEHM*(1-ceEHM)*crEHM, (1-cHWF)*cEHM*(1-ceEHM)*(1-crEHM),
                             (1+cHWF*chHWF)*(1-cEHM), (1+cHWF*chHWF)*(1+cEHM*ceEHM), (1+cHWF*chHWF)*cEHM*(1-ceEHM)*crEHM + (cHWF*(1-chHWF)*crHWF)*(1-cEHM), (1+cHWF*chHWF)*cEHM*(1-ceEHM)*(1-crEHM) + (cHWF*(1-chHWF)*(1-crHWF))*(1-cEHM), 0,
                             (cHWF*(1-chHWF)*crHWF)*(1+cEHM*ceEHM), (cHWF*(1-chHWF)*(1-crHWF))*(1+cEHM*ceEHM), (cHWF*(1-chHWF)*crHWF)*cEHM*(1-ceEHM)*crEHM + (cHWF*(1-chHWF)*(1-crHWF))*cEHM*(1-ceEHM)*crEHM,
                             (cHWF*(1-chHWF)*crHWF)*cEHM*(1-ceEHM)*(1-crEHM), (cHWF*(1-chHWF)*(1-crHWF))*cEHM*(1-ceEHM)*(1-crEHM))/4
  tMatrix['WH','HR', c('WH', 'HH', 'HR', 'HB',
                       'WR', 'RR', 'RB')] <- c( 1-cHWF, 1+cHWF*chHWF, cHWF*(1-chHWF)*crHWF + 1+cHWF*chHWF, cHWF*(1-chHWF)*(1-crHWF),
                                                1-cHWF, cHWF*(1-chHWF)*crHWF, cHWF*(1-chHWF)*(1-crHWF))/4
  tMatrix['WH','HB', c('WH', 'HH', 'HR', 'HB',
                       'WB', 'RB', 'BB')] <- c( 1-cHWF, 1+cHWF*chHWF, cHWF*(1-chHWF)*crHWF, cHWF*(1-chHWF)*(1-crHWF) + 1+cHWF*chHWF,
                                                1-cHWF, cHWF*(1-chHWF)*crHWF, cHWF*(1-chHWF)*(1-crHWF))/4
  tMatrix['WH','EE', c('WE', 'HE', 'ER', 'EB')] <- c( 1-cHWF, 1+cHWF*chHWF, cHWF*(1-chHWF)*crHWF, cHWF*(1-chHWF)*(1-crHWF))/2
  tMatrix['WH','ER', c('WE', 'HE', 'ER', 'EB',
                       'WR', 'HR', 'RR', 'RB')] <- c( 1-cHWF, 1+cHWF*chHWF, cHWF*(1-chHWF)*crHWF, cHWF*(1-chHWF)*(1-crHWF),
                                                      1-cHWF, 1+cHWF*chHWF, cHWF*(1-chHWF)*crHWF, cHWF*(1-chHWF)*(1-crHWF))/4
  tMatrix['WH','EB', c('WE', 'HE', 'ER', 'EB',
                       'WB', 'HB', 'RB', 'BB')] <- c( 1-cHWF, 1+cHWF*chHWF, cHWF*(1-chHWF)*crHWF, cHWF*(1-chHWF)*(1-crHWF),
                                                      1-cHWF, 1+cHWF*chHWF, cHWF*(1-chHWF)*crHWF, cHWF*(1-chHWF)*(1-crHWF))/4
  tMatrix['WH','RR', c('WR', 'HR', 'RR', 'RB')] <- c( 1-cHWF, 1+cHWF*chHWF, cHWF*(1-chHWF)*crHWF, cHWF*(1-chHWF)*(1-crHWF))/2
  tMatrix['WH','RB', c('WR', 'HR', 'RR', 'RB',
                       'WB', 'HB', 'BB')] <- c( 1-cHWF, 1+cHWF*chHWF, cHWF*(1-chHWF)*crHWF, cHWF*(1-chHWF)*(1-crHWF) + cHWF*(1-chHWF)*crHWF,
                                                1-cHWF, 1+cHWF*chHWF, cHWF*(1-chHWF)*(1-crHWF))/4
  tMatrix['WH','BB', c('WB', 'HB', 'RB', 'BB')] <- c( 1-cHWF, 1+cHWF*chHWF, cHWF*(1-chHWF)*crHWF, cHWF*(1-chHWF)*(1-crHWF))/2

  # males
  tMatrix['WE','WH', ] <- c( (1-cEWF)*(1-cHWM)*(1-dEW), (1-cEWF)*(1+cHWM*chHWM)*(1-dEH), (1+cEWF*ceEWF)*(1-cHWM)*(1-dEW),
                             (1-cEWF)*(1-cHWM)*dEW*drEW + (cEWF*(1-ceEWF)*crEWF)*(1-cHWM)*(1-dEW) + (1-cEWF)*(1+cHWM*chHWM)*dEH*drEH + (1-cEWF)*(cHWM*(1-chHWM)*crHWM),
                             (1-cEWF)*(1-cHWM)*dEW*(1-drEW) + (cEWF*(1-ceEWF)*(1-crEWF))*(1-cHWM)*(1-dEW) + (1-cEWF)*(1+cHWM*chHWM)*dEH*(1-drEH) + (1-cEWF)*(cHWM*(1-chHWM)*(1-crHWM)),
                             0, (1+cEWF*ceEWF)*(1+cHWM*chHWM)*(1-dEH), (cEWF*(1-ceEWF)*crEWF)*(1+cHWM*chHWM)*(1-dEH), (cEWF*(1-ceEWF)*(1-crEWF))*(1+cHWM*chHWM)*(1-dEH),
                             (1+cEWF*ceEWF)*(1-cHWM)*dEW*deEW + (1+cEWF*ceEWF)*(1+cHWM*chHWM)*dEH*deEH,
                             (1+cEWF*ceEWF)*(1-cHWM)*dEW*(1-deEW)*drEW + (1+cEWF*ceEWF)*(1+cHWM*chHWM)*dEH*(1-deEH)*drEH + (1+cEWF*ceEWF)*(cHWM*(1-chHWM)*crHWM),
                             (1+cEWF*ceEWF)*(1-cHWM)*dEW*(1-deEW)*(1-drEW) + (1+cEWF*ceEWF)*(1+cHWM*chHWM)*dEH*(1-deEH)*(1-drEH) + (1+cEWF*ceEWF)*(cHWM*(1-chHWM)*(1-crHWM)),
                             (cEWF*(1-ceEWF)*crEWF)*(1-cHWM)*dEW*drEW + (cEWF*(1-ceEWF)*crEWF)*(1+cHWM*chHWM)*dEH*drEH + (cEWF*(1-ceEWF)*crEWF)*(cHWM*(1-chHWM)*crHWM),
                             (cEWF*(1-ceEWF)*crEWF)*(1-cHWM)*dEW*(1-drEW) + (cEWF*(1-ceEWF)*(1-crEWF))*(1-cHWM)*dEW*drEW + (cEWF*(1-ceEWF)*crEWF)*(1+cHWM*chHWM)*dEH*(1-drEH) + (cEWF*(1-ceEWF)*(1-crEWF))*(1+cHWM*chHWM)*dEH*drEH + (cEWF*(1-ceEWF)*(1-crEWF))*(cHWM*(1-chHWM)*crHWM) + (cEWF*(1-ceEWF)*crEWF)*(cHWM*(1-chHWM)*(1-crHWM)),
                             (cEWF*(1-ceEWF)*(1-crEWF))*(1-cHWM)*dEW*(1-drEW) + (cEWF*(1-ceEWF)*(1-crEWF))*(1+cHWM*chHWM)*dEH*(1-drEH) + (cEWF*(1-ceEWF)*(1-crEWF))*(cHWM*(1-chHWM)*(1-crHWM)))/4
  tMatrix['WR','WH', c('WW', 'WH', 'WR', 'WB',
                       'HR', 'RR', 'RB')] <- c( 1-cHWM, 1+cHWM*chHWM, cHWM*(1-chHWM)*crHWM + (1-cHWM), cHWM*(1-chHWM)*(1-crHWM),
                                                1+cHWM*chHWM, cHWM*(1-chHWM)*crHWM, cHWM*(1-chHWM)*(1-crHWM))/4
  tMatrix['WB','WH', c('WW', 'WH', 'WR', 'WB',
                       'HB', 'RB', 'BB')] <- c( 1-cHWM, 1+cHWM*chHWM, cHWM*(1-chHWM)*crHWM, cHWM*(1-chHWM)*(1-crHWM) + (1-cHWM),
                                                1+cHWM*chHWM, cHWM*(1-chHWM)*crHWM, cHWM*(1-chHWM)*(1-crHWM))/4
  tMatrix['HH','WH', c('WH', 'HH', 'HR', 'HB')] <- c( (1-cHWM)*(1-dHW), (1-cHWM)*dHW*dhHW + (1+cHWM*chHWM), (1-cHWM)*dHW*(1-dhHW)*drHW + cHWM*(1-chHWM)*crHWM, (1-cHWM)*dHW*(1-dhHW)*(1-drHW) + cHWM*(1-chHWM)*(1-crHWM))/2
  tMatrix['HE','WH', ] <- c( 0, (1-cEHF)*(1-cHWM)*(1-dEW), (1+cEHF*ceEHF)*(1-cHWM)*(1-dEW), (cEHF*(1-ceEHF)*crEHF)*(1-cHWM)*(1-dEW), (cEHF*(1-ceEHF)*(1-crEHF))*(1-cHWM)*(1-dEW), (1-cEHF)*(1+cHWM*chHWM)*(1-dEH),
                             (1+cEHF*ceEHF)*(1+cHWM*chHWM)*(1-dEH), (1-cEHF)*(1-cHWM)*dEW*drEW + (1-cEHF)*(1+cHWM*chHWM)*dEH*drEH + (cEHF*(1-ceEHF)*crEHF)*(1+cHWM*chHWM)*(1-dEH) + (1-cEHF)*(cHWM*(1-chHWM)*crHWM),
                             (1-cEHF)*(1-cHWM)*dEW*(1-drEW) + (1-cEHF)*(1+cHWM*chHWM)*dEH*(1-drEH) + (cEHF*(1-ceEHF)*(1-crEHF))*(1+cHWM*chHWM)*(1-dEH) + (1-cEHF)*(cHWM*(1-chHWM)*(1-crHWM)),
                             (1+cEHF*ceEHF)*(1-cHWM)*dEW*deEW + (1+cEHF*ceEHF)*(1+cHWM*chHWM)*dEH*deEH, (1+cEHF*ceEHF)*(1-cHWM)*dEW*(1-deEW)*drEW + (1+cEHF*ceEHF)*(1+cHWM*chHWM)*dEH*(1-deEH)*drEH + (1+cEHF*ceEHF)*(cHWM*(1-chHWM)*crHWM),
                             (1+cEHF*ceEHF)*(1-cHWM)*dEW*(1-deEW)*(1-drEW) + (1+cEHF*ceEHF)*(1+cHWM*chHWM)*dEH*(1-deEH)*(1-drEH) + (1+cEHF*ceEHF)*(cHWM*(1-chHWM)*(1-crHWM)),
                             (cEHF*(1-ceEHF)*crEHF)*(1-cHWM)*dEW*drEW + (cEHF*(1-ceEHF)*crEHF)*(1+cHWM*chHWM)*dEH*drEH + (cEHF*(1-ceEHF)*crEHF)*(cHWM*(1-chHWM)*crHWM),
                             (cEHF*(1-ceEHF)*crEHF)*(1-cHWM)*dEW*(1-drEW) + (cEHF*(1-ceEHF)*(1-crEHF))*(1-cHWM)*dEW*drEW + (cEHF*(1-ceEHF)*crEHF)*(1+cHWM*chHWM)*dEH*(1-drEH) + (cEHF*(1-ceEHF)*(1-crEHF))*(1+cHWM*chHWM)*dEH*drEH + (cEHF*(1-ceEHF)*(1-crEHF))*(cHWM*(1-chHWM)*crHWM) + (cEHF*(1-ceEHF)*crEHF)*(cHWM*(1-chHWM)*(1-crHWM)),
                             (cEHF*(1-ceEHF)*(1-crEHF))*(1-cHWM)*dEW*(1-drEW) + (cEHF*(1-ceEHF)*(1-crEHF))*(1+cHWM*chHWM)*dEH*(1-drEH) + (cEHF*(1-ceEHF)*(1-crEHF))*(cHWM*(1-chHWM)*(1-crHWM)))/4
  tMatrix['HR','WH', c('WH', 'HH', 'HR', 'HB',
                       'WR', 'RR', 'RB')] <- c( (1-cHWM)*(1-dHW), (1-cHWM)*dHW*dhHW + (1+cHWM*chHWM), (1-cHWM)*dHW*(1-dhHW)*drHW + cHWM*(1-chHWM)*crHWM + (1+cHWM*chHWM), (1-cHWM)*dHW*(1-dhHW)*(1-drHW) + cHWM*(1-chHWM)*(1-crHWM),
                                                (1-cHWM)*(1-dHW), (1-cHWM)*dHW*drHW + cHWM*(1-chHWM)*crHWM, (1-cHWM)*dHW*(1-drHW) + cHWM*(1-chHWM)*(1-crHWM))/4
  tMatrix['HB','WH', c('WH', 'HH', 'HR', 'HB',
                       'WB', 'RB', 'BB')] <- c( (1-cHWM)*(1-dHW), (1-cHWM)*dHW*dhHW + (1+cHWM*chHWM), (1-cHWM)*dHW*(1-dhHW)*drHW + cHWM*(1-chHWM)*crHWM, (1-cHWM)*dHW*(1-dhHW)*(1-drHW) + cHWM*(1-chHWM)*(1-crHWM) + (1+cHWM*chHWM),
                                                (1-cHWM)*(1-dHW), (1-cHWM)*dHW*drHW + cHWM*(1-chHWM)*crHWM, (1-cHWM)*dHW*(1-drHW) + cHWM*(1-chHWM)*(1-crHWM))/4
  tMatrix['EE','WH', c('WE', 'HE', 'ER',
                       'EB', 'EE')] <- c( (1-cHWM)*(1-dEW), (1+cHWM*chHWM)*(1-dEH), (1-cHWM)*dEW*(1-deEW)*drEW + (1+cHWM*chHWM)*dEH*(1-deEH)*drEH + cHWM*(1-chHWM)*crHWM,
                                          (1-cHWM)*dEW*(1-deEW)*(1-drEW) + (1+cHWM*chHWM)*dEH*(1-deEH)*(1-drEH) + cHWM*(1-chHWM)*(1-crHWM), (1-cHWM)*dEW*deEW + (1+cHWM*chHWM)*dEH*deEH)/2
  tMatrix['ER','WH', c('WE', 'HE', 'ER',
                       'EB', 'EE', 'WR',
                       'HR', 'RR', 'RB')] <- c( (1-cHWM)*(1-dEW), (1+cHWM*chHWM)*(1-dEH), (1-cHWM)*dEW*(1-deEW)*drEW + (1+cHWM*chHWM)*dEH*(1-deEH)*drEH + cHWM*(1-chHWM)*crHWM,
                                                (1-cHWM)*dEW*(1-deEW)*(1-drEW) + (1+cHWM*chHWM)*dEH*(1-deEH)*(1-drEH) + cHWM*(1-chHWM)*(1-crHWM), (1-cHWM)*dEW*deEW + (1+cHWM*chHWM)*dEH*deEH, (1-cHWM)*(1-dEW),
                                                (1+cHWM*chHWM)*(1-dEH), (1-cHWM)*dEW*drEW + (1+cHWM*chHWM)*dEH*drEH + cHWM*(1-chHWM)*crHWM, (1-cHWM)*dEW*(1-drEW) + (1+cHWM*chHWM)*dEH*(1-drEH) + cHWM*(1-chHWM)*(1-crHWM))/4
  tMatrix['EB','WH', c('WE', 'HE', 'ER',
                       'EB', 'EE', 'WB',
                       'HB', 'RB', 'BB')] <- c( (1-cHWM)*(1-dEW), (1+cHWM*chHWM)*(1-dEH), (1-cHWM)*dEW*(1-deEW)*drEW + (1+cHWM*chHWM)*dEH*(1-deEH)*drEH + cHWM*(1-chHWM)*crHWM,
                                                (1-cHWM)*dEW*(1-deEW)*(1-drEW) + (1+cHWM*chHWM)*dEH*(1-deEH)*(1-drEH) + cHWM*(1-chHWM)*(1-crHWM), (1-cHWM)*dEW*deEW + (1+cHWM*chHWM)*dEH*deEH, (1-cHWM)*(1-dEW),
                                                (1+cHWM*chHWM)*(1-dEH), (1-cHWM)*dEW*drEW + (1+cHWM*chHWM)*dEH*drEH + cHWM*(1-chHWM)*crHWM, (1-cHWM)*dEW*(1-drEW) + (1+cHWM*chHWM)*dEH*(1-drEH) + cHWM*(1-chHWM)*(1-crHWM))/4
  tMatrix['RR','WH', c('WR', 'HR', 'RR', 'RB')] <- c( 1-cHWM, 1+cHWM*chHWM, cHWM*(1-chHWM)*crHWM, cHWM*(1-chHWM)*(1-crHWM))/2
  tMatrix['RB','WH', c('WR', 'HR', 'RR', 'RB',
                       'WB', 'HB', 'BB')] <- c( 1-cHWM, 1+cHWM*chHWM, cHWM*(1-chHWM)*crHWM, cHWM*(1-chHWM)*(1-crHWM) + cHWM*(1-chHWM)*crHWM,
                                                1-cHWM, 1+cHWM*chHWM, cHWM*(1-chHWM)*(1-crHWM))/4
  tMatrix['BB','WH', c('WB', 'HB', 'RB', 'BB')] <- c( 1-cHWM, 1+cHWM*chHWM, cHWM*(1-chHWM)*crHWM, cHWM*(1-chHWM)*(1-crHWM))/2


  # both
  tMatrix['WE','WE', ] <- c( (1-cEWF)*(1-cEWM)*(1-dEW), 0, (1+cEWF*ceEWF)*(1-cEWM)*(1-dEW) + (1-cEWF)*(1+cEWM*ceEWM), (1-cEWF)*(1-cEWM)*dEW*drEW + (cEWF*(1-ceEWF)*crEWF)*(1-cEWM)*(1-dEW) + (1-cEWF)*(cEWM*(1-ceEWM)*crEWM),
                             (1-cEWF)*(1-cEWM)*dEW*(1-drEW) + (cEWF*(1-ceEWF)*(1-crEWF))*(1-cEWM)*(1-dEW) + (1-cEWF)*(cEWM*(1-ceEWM)*(1-crEWM)), 0, 0, 0, 0,
                             (1+cEWF*ceEWF)*(1-cEWM)*dEW*deEW + (1+cEWF*ceEWF)*(1+cEWM*ceEWM), (1+cEWF*ceEWF)*(1-cEWM)*dEW*(1-deEW)*drEW + (cEWF*(1-ceEWF)*crEWF)*(1+cEWM*ceEWM) + (1+cEWF*ceEWF)*(cEWM*(1-ceEWM)*crEWM),
                             (1+cEWF*ceEWF)*(1-cEWM)*dEW*(1-deEW)*(1-drEW) + (cEWF*(1-ceEWF)*(1-crEWF))*(1+cEWM*ceEWM) + (1+cEWF*ceEWF)*(cEWM*(1-ceEWM)*(1-crEWM)),
                             (cEWF*(1-ceEWF)*crEWF)*(1-cEWM)*dEW*drEW + (cEWF*(1-ceEWF)*crEWF)*(cEWM*(1-ceEWM)*crEWM),
                             (cEWF*(1-ceEWF)*crEWF)*(1-cEWM)*dEW*(1-drEW) + (cEWF*(1-ceEWF)*(1-crEWF))*(1-cEWM)*dEW*drEW + (cEWF*(1-ceEWF)*(1-crEWF))*(cEWM*(1-ceEWM)*crEWM) + (cEWF*(1-ceEWF)*crEWF)*(cEWM*(1-ceEWM)*(1-crEWM)),
                             (cEWF*(1-ceEWF)*(1-crEWF))*(1-cEWM)*dEW*(1-drEW) + (cEWF*(1-ceEWF)*(1-crEWF))*(cEWM*(1-ceEWM)*(1-crEWM)))/4

  # females
  tMatrix['WE','WR', c('WW', 'WE', 'WR', 'WB', 'EE',
                       'ER', 'EB', 'RR', 'RB', 'BB')] <- c( (1-cEWF)*(1-dEW), (1+cEWF*ceEWF)*(1-dEW), (1-cEWF)*dEW*drEW + (cEWF*(1-ceEWF)*crEWF)*(1-dEW) + (1-cEWF), (1-cEWF)*dEW*(1-drEW) + (cEWF*(1-ceEWF)*(1-crEWF))*(1-dEW),
                                                            (1+cEWF*ceEWF)*dEW*deEW, (1+cEWF*ceEWF)*dEW*(1-deEW)*drEW + (1+cEWF*ceEWF), (1+cEWF*ceEWF)*dEW*(1-deEW)*(1-drEW), (cEWF*(1-ceEWF)*crEWF)*dEW*drEW + (cEWF*(1-ceEWF)*crEWF),
                                                            (cEWF*(1-ceEWF)*crEWF)*dEW*(1-drEW) + (cEWF*(1-ceEWF)*(1-crEWF))*dEW*drEW + (cEWF*(1-ceEWF)*(1-crEWF)), (cEWF*(1-ceEWF)*(1-crEWF))*dEW*(1-drEW))/4
  tMatrix['WE','WB', c('WW', 'WE', 'WR', 'WB', 'EE',
                       'ER', 'EB', 'RR', 'RB', 'BB')] <- c( (1-cEWF)*(1-dEW), (1+cEWF*ceEWF)*(1-cEWF), (1-cEWF)*dEW*drEW + (cEWF*(1-ceEWF)*crEWF)*(1-dEW), (1-cEWF)*dEW*(1-drEW) + (cEWF*(1-ceEWF)*(1-crEWF))*(1-dEW) + (1-cEWF),
                                                            (1+cEWF*ceEWF)*dEW*deEW, (1+cEWF*ceEWF)*dEW*(1-deEW)*drEW, (1+cEWF*ceEWF)*dEW*(1-deEW)*(1-drEW) + (1+cEWF*ceEWF), (cEWF*(1-ceEWF)*crEWF)*dEW*drEW,
                                                            (cEWF*(1-ceEWF)*crEWF)*dEW*(1-drEW) + (cEWF*(1-ceEWF)*(1-crEWF))*dEW*drEW + (cEWF*(1-ceEWF)*crEWF), (cEWF*(1-ceEWF)*(1-crEWF))*dEW*(1-drEW) + (cEWF*(1-ceEWF)*(1-crEWF)))/4
  tMatrix['WE','HH', c('WH', 'HE', 'HR', 'HB', 'WR', 'WB',
                       'EE', 'ER', 'EB', 'RR', 'RB', 'BB')] <- c( (1-cEWF)*(1-dEW), (1+cEWF*ceEWF)*(1-dEW), (cEWF*(1-ceEWF)*crEWF)*(1-dEW), (cEWF*(1-ceEWF)*(1-crEWF))*(1-dEW), (1-cEWF)*dEW*drEW,
                                                                  (1-cEWF)*dEW*(1-drEW), (1+cEWF*ceEWF)*dEW*deEW, (1+cEWF*ceEWF)*dEW*(1-deEW)*drEW, (1+cEWF*ceEWF)*dEW*(1-deEW)*(1-drEW),
                                                                  (cEWF*(1-ceEWF)*crEWF)*dEW*drEW, (cEWF*(1-ceEWF)*crEWF)*dEW*(1-drEW) + (cEWF*(1-ceEWF)*(1-crEWF))*dEW*drEW, (cEWF*(1-ceEWF)*(1-crEWF))*dEW*(1-drEW))/2
  tMatrix['WE','HE', ] <- c( 0, (1-cEWF)*(1-cEHM)*(1-dEH), (1-cEWF)*(1+cEHM*ceEHM), (1-cEWF)*(1-cEHM)*dEH*drEH + (1-cEWF)*(cEHM*(1-ceEHM)*crEHM), (1-cEWF)*(1-cEHM)*dEH*(1-drEH) + (1-cEWF)*(cEHM*(1-ceEHM)*(1-crEHM)),
                             0, (1+cEWF*ceEWF)*(1-cEHM)*(1-dEH), (cEWF*(1-ceEWF)*crEWF)*(1-cEHM)*(1-dEH), (cEWF*(1-ceEWF)*(1-crEWF))*(1-cEHM)*(1-dEH), (1+cEWF*ceEWF)*(1-cEHM)*dEH*deEH + (1+cEWF*ceEWF)*(1+cEHM*ceEHM),
                             (1+cEWF*ceEWF)*(1-cEHM)*dEH*(1-deEH)*drEH + (cEWF*(1-ceEWF)*crEWF)*(1+cEHM*ceEHM) + (1+cEWF*ceEWF)*(cEHM*(1-ceEHM)*crEHM),
                             (1+cEWF*ceEWF)*(1-cEHM)*dEH*(1-deEH)*(1-drEH) + (cEWF*(1-ceEWF)*(1-crEWF))*(1+cEHM*ceEHM) + (1+cEWF*ceEWF)*(cEHM*(1-ceEHM)*(1-crEHM)),
                             (cEWF*(1-ceEWF)*crEWF)*(1-cEHM)*dEH*drEH + (cEWF*(1-ceEWF)*crEWF)*(cEHM*(1-ceEHM)*crEHM),
                             (cEWF*(1-ceEWF)*crEWF)*(1-cEHM)*dEH*(1-drEH) + (cEWF*(1-ceEWF)*(1-crEWF))*(1-cEHM)*dEH*drEH + (cEWF*(1-ceEWF)*(1-crEWF))*(cEHM*(1-ceEHM)*crEHM) + (cEWF*(1-ceEWF)*crEWF)*(cEHM*(1-ceEHM)*(1-crEHM)),
                             (cEWF*(1-ceEWF)*(1-crEWF))*(1-cEHM)*dEH*(1-drEH) + (cEWF*(1-ceEWF)*(1-crEWF))*(cEHM*(1-ceEHM)*(1-crEHM)))/4
  tMatrix['WE','HR', c('WH', 'HE', 'HR', 'HB', 'WR', 'WB',
                       'EE', 'ER', 'EB', 'RR', 'RB', 'BB')] <- c( (1-cEWF)*(1-dEW), (1+cEWF*ceEWF)*(1-dEW), (cEWF*(1-ceEWF)*crEWF)*(1-dEW), (cEWF*(1-ceEWF)*(1-crEWF))*(1-dEW), (1-cEWF)*dEW*drEW + (1-cEWF), (1-cEWF)*dEW*(1-drEW),
                                                                  (1+cEWF*ceEWF)*dEW*deEW, (1+cEWF*ceEWF)*dEW*(1-deEW)*drEW + (1+cEWF*ceEWF), (1+cEWF*ceEWF)*dEW*(1-deEW)*(1-drEW), (cEWF*(1-ceEWF)*crEWF)*dEW*drEW + (cEWF*(1-ceEWF)*crEWF),
                                                                  (cEWF*(1-ceEWF)*crEWF)*dEW*(1-drEW) + (cEWF*(1-ceEWF)*(1-crEWF))*dEW*drEW + (cEWF*(1-ceEWF)*(1-crEWF)), (cEWF*(1-ceEWF)*(1-crEWF))*dEW*(1-drEW))/4
  tMatrix['WE','HB', c('WH', 'HE', 'HR', 'HB', 'WR', 'WB',
                       'EE', 'ER', 'EB', 'RR', 'RB', 'BB')] <- c( (1-cEWF)*(1-dEW), (1+cEWF*ceEWF)*(1-dEW), (cEWF*(1-ceEWF)*crEWF)*(1-dEW), (cEWF*(1-ceEWF)*(1-crEWF))*(1-dEW), (1-cEWF)*dEW*drEW, (1-cEWF)*dEW*(1-drEW) + (1-cEWF),
                                                                  (1+cEWF*ceEWF)*dEW*deEW, (1+cEWF*ceEWF)*dEW*(1-deEW)*drEW, (1+cEWF*ceEWF)*dEW*(1-deEW)*(1-drEW) + (1+cEWF*ceEWF), (cEWF*(1-ceEWF)*crEWF)*dEW*drEW,
                                                                  (cEWF*(1-ceEWF)*crEWF)*dEW*(1-drEW) + (cEWF*(1-ceEWF)*(1-crEWF))*dEW*drEW + (cEWF*(1-ceEWF)*crEWF), (cEWF*(1-ceEWF)*(1-crEWF))*dEW*(1-drEW) + (cEWF*(1-ceEWF)*(1-crEWF)))/4
  tMatrix['WE','EE', c('WE', 'EE', 'ER', 'EB')] <- c( 1-cEWF, 1+cEWF*ceEWF, cEWF*(1-ceEWF)*crEWF, cEWF*(1-ceEWF)*(1-crEWF))/2
  tMatrix['WE','ER', c('WE', 'EE', 'ER', 'EB',
                       'WR', 'RR', 'RB')] <- c( 1-cEWF, 1+cEWF*ceEWF, cEWF*(1-ceEWF)*crEWF + 1+cEWF*ceEWF, cEWF*(1-ceEWF)*(1-crEWF),
                                                1-cEWF, cEWF*(1-ceEWF)*crEWF, cEWF*(1-ceEWF)*(1-crEWF))/4
  tMatrix['WE','EB', c('WE', 'EE', 'ER', 'EB',
                       'WB', 'RB', 'BB')] <- c( 1-cEWF, 1+cEWF*ceEWF, cEWF*(1-ceEWF)*crEWF, cEWF*(1-ceEWF)*(1-crEWF) + 1+cEWF*ceEWF,
                                                1-cEWF, cEWF*(1-ceEWF)*crEWF, cEWF*(1-ceEWF)*(1-crEWF))/4
  tMatrix['WE','RR', c('WR', 'ER', 'RR', 'RB')] <- c( 1-cEWF, 1+cEWF*ceEWF, cEWF*(1-ceEWF)*crEWF, cEWF*(1-ceEWF)*(1-crEWF))/2
  tMatrix['WE','RB', c('WR', 'ER', 'RR', 'RB',
                       'WB', 'EB', 'BB')] <- c( 1-cEWF, 1+cEWF*ceEWF, cEWF*(1-ceEWF)*crEWF, cEWF*(1-ceEWF)*(1-crEWF) + cEWF*(1-ceEWF)*crEWF,
                                                1-cEWF, 1+cEWF*ceEWF, cEWF*(1-ceEWF)*(1-crEWF))/4
  tMatrix['WE','BB', c('WB', 'EB', 'RB', 'BB')] <- c( 1-cEWF, 1+cEWF*ceEWF, cEWF*(1-ceEWF)*crEWF, cEWF*(1-ceEWF)*(1-crEWF))/2

  # males
  tMatrix['WR','WE', c('WW', 'WE', 'WR', 'WB',
                       'ER', 'RR', 'RB')] <- c( 1-cEWM, 1+cEWM*ceEWM, cEWM*(1-ceEWM)*crEWM + 1-cEWM, cEWM*(1-ceEWM)*(1-crEWM),
                                                1+cEWM*ceEWM, cEWM*(1-ceEWM)*crEWM, cEWM*(1-ceEWM)*(1-crEWM))/4
  tMatrix['WB','WE', c('WW', 'WE', 'WR', 'WB',
                       'EB', 'RB', 'BB')] <- c( 1-cEWM, 1+cEWM*ceEWM, cEWM*(1-ceEWM)*crEWM, cEWM*(1-ceEWM)*(1-crEWM) + 1-cEWM,
                                                1+cEWM*ceEWM, cEWM*(1-ceEWM)*crEWM, cEWM*(1-ceEWM)*(1-crEWM))/4
  tMatrix['HH','WE', c('WH', 'HE', 'HR', 'HB')] <- c( (1-cEWM)*(1-dHW), 1+cEWM*ceEWM, (1-cEWM)*dHW*drHW + cEWM*(1-ceEWM)*crEWM, (1-cEWM)*dHW*(1-drHW) + cEWM*(1-ceEWM)*(1-crEWM))/2
  tMatrix['HE','WE', ] <- c( 0, (1-cEHF)*(1-cEWM)*(1-dEW), (1+cEHF*ceEHF)*(1-cEWM)*(1-dEW), (cEHF*(1-ceEHF)*crEHF)*(1-cEWM)*(1-dEW), (cEHF*(1-ceEHF)*(1-crEHF))*(1-cEWM)*(1-dEW),
                             0, (1-cEHF)*(1+cEWM*ceEWM), (1-cEHF)*(1-cEWM)*dEW*drEW + (1-cEHF)*(cEWM*(1-ceEWM)*crEWM), (1-cEHF)*(1-cEWM)*dEW*(1-drEW) + (1-cEHF)*(cEWM*(1-ceEWM)*(1-crEWM)),
                             (1+cEHF*ceEHF)*(1-cEWM)*dEW*deEW + (1+cEHF*ceEHF)*(1+cEWM*ceEWM), (1+cEHF*ceEHF)*(1-cEWM)*dEW*(1-deEW)*drEW + (cEHF*(1-ceEHF)*crEHF)*(1+cEWM*ceEWM) + (1+cEHF*ceEHF)*(cEWM*(1-ceEWM)*crEWM),
                             (1+cEHF*ceEHF)*(1-cEWM)*dEW*(1-deEW)*(1-drEW) + (cEHF*(1-ceEHF)*(1-crEHF))*(1+cEWM*ceEWM) + (1+cEHF*ceEHF)*(cEWM*(1-ceEWM)*(1-crEWM)),
                             (cEHF*(1-ceEHF)*crEHF)*(1-cEWM)*dEW*drEW + (cEHF*(1-ceEHF)*crEHF)*(cEWM*(1-ceEWM)*crEWM),
                             (cEHF*(1-ceEHF)*crEHF)*(1-cEWM)*dEW*(1-drEW) + (cEHF*(1-ceEHF)*(1-crEHF))*(1-cEWM)*dEW*drEW + (cEHF*(1-ceEHF)*(1-crEHF))*(cEWM*(1-ceEWM)*crEWM) + (cEHF*(1-ceEHF)*crEHF)*(cEWM*(1-ceEWM)*(1-crEWM)),
                             (cEHF*(1-ceEHF)*(1-crEHF))*(1-cEWM)*dEW*(1-drEW) + (cEHF*(1-ceEHF)*(1-crEHF))*(cEWM*(1-ceEWM)*(1-crEWM)))/4
  tMatrix['HR','WE', c('WH', 'HE', 'HR', 'HB',
                       'WR', 'ER', 'RR', 'RB')] <- c( (1-cEWM)*(1-dHW), 1+cEWM*ceEWM, (1-cEWM)*dHW*drHW + cEWM*(1-ceEWM)*crEWM, (1-cEWM)*dHW*(1-drHW) + cEWM*(1-ceEWM)*(1-crEWM),
                                                      (1-cEWM)*(1-dHW), 1+cEWM*ceEWM, (1-cEWM)*dHW*drHW + cEWM*(1-ceEWM)*crEWM, (1-cEWM)*dHW*(1-drHW) + cEWM*(1-ceEWM)*(1-crEWM))/4
  tMatrix['HB','WE', c('WH', 'HE', 'HR', 'HB',
                       'WB', 'EB', 'RB', 'BB')] <- c( (1-cEWM)*(1-dHW), 1+cEWM*ceEWM, (1-cEWM)*dHW*drHW + cEWM*(1-ceEWM)*crEWM, (1-cEWM)*dHW*(1-drHW) + cEWM*(1-ceEWM)*(1-crEWM),
                                                      (1-cEWM)*(1-dHW), 1+cEWM*ceEWM, (1-cEWM)*dHW*drHW + cEWM*(1-ceEWM)*crEWM, (1-cEWM)*dHW*(1-drHW) + cEWM*(1-ceEWM)*(1-crEWM))/4
  tMatrix['EE','WE', c('WE', 'EE', 'ER', 'EB')] <- c( (1-cEWM)*(1-dEW), (1-cEWM)*dEW*deEW + (1+cEWM*ceEWM), (1-cEWM)*dEW*(1-deEW)*drEW + cEWM*(1-ceEWM)*crEWM, (1-cEWM)*dEW*(1-deEW)*(1-drEW) + cEWM*(1-ceEWM)*(1-crEWM))/2
  tMatrix['ER','WE', c('WE', 'EE', 'ER', 'EB',
                       'WR', 'RR', 'RB')] <- c( (1-cEWM)*(1-dEW), (1-cEWM)*dEW*deEW + (1+cEWM*ceEWM), (1-cEWM)*dEW*(1-deEW)*drEW + cEWM*(1-ceEWM)*crEWM + (1+cEWM*ceEWM), (1-cEWM)*dEW*(1-deEW)*(1-drEW) + cEWM*(1-ceEWM)*(1-crEWM),
                                                (1-cEWM)*(1-dEW), (1-cEWM)*dEW*drEW + cEWM*(1-ceEWM)*crEWM, (1-cEWM)*dEW*(1-drEW) + cEWM*(1-ceEWM)*(1-crEWM))/4
  tMatrix['EB','WE', c('WE', 'EE', 'ER', 'EB',
                       'WB', 'RB', 'BB')] <- c( (1-cEWM)*(1-dEW), (1-cEWM)*dEW*deEW + (1+cEWM*ceEWM), (1-cEWM)*dEW*(1-deEW)*drEW + cEWM*(1-ceEWM)*crEWM, (1-cEWM)*dEW*(1-deEW)*(1-drEW) + cEWM*(1-ceEWM)*(1-crEWM) + (1+cEWM*ceEWM),
                                                (1-cEWM)*(1-dEW), (1-cEWM)*dEW*drEW + cEWM*(1-ceEWM)*crEWM, (1-cEWM)*dEW*(1-drEW) + cEWM*(1-ceEWM)*(1-crEWM))/4
  tMatrix['RR','WE', c('WR', 'ER', 'RR', 'RB')] <- c( 1-cEWM, 1+cEWM*ceEWM, cEWM*(1-ceEWM)*crEWM, cEWM*(1-ceEWM)*(1-crEWM))/2
  tMatrix['RB','WE', c('WR', 'ER', 'RR', 'RB',
                       'WB', 'EB', 'BB')] <- c( 1-cEWM, 1+cEWM*ceEWM, cEWM*(1-ceEWM)*crEWM, cEWM*(1-ceEWM)*(1-crEWM) + cEWM*(1-ceEWM)*crEWM,
                                                1-cEWM, 1+cEWM*ceEWM, cEWM*(1-ceEWM)*(1-crEWM))/4
  tMatrix['BB','WE', c('WB', 'EB', 'RB', 'BB')] <- c( 1-cEWM, 1+cEWM*ceEWM, cEWM*(1-ceEWM)*crEWM, cEWM*(1-ceEWM)*(1-crEWM))/2


  # both
  tMatrix['HE','HE', ] <- c( 0, 0, 0, 0, 0, (1-cEHF)*(1-cEHM)*(1-dEH), (1+cEHF*ceEHF)*(1-cEHM)*(1-dEH) + (1-cEHF)*(1+cEHM*ceEHM), (1-cEHF)*(1-cEHM)*dEH*drEH + (cEHF*(1-ceEHF)*crEHF)*(1-cEHM)*(1-dEH) + (1-cEHF)*(cEHM*(1-ceEHM)*crEHM),
                             (1-cEHF)*(1-cEHM)*dEH*(1-drEH) + (cEHF*(1-ceEHF)*(1-crEHF))*(1-cEHM)*(1-dEH) + (1-cEHF)*(cEHM*(1-ceEHM)*(1-crEHM)), (1+cEHF*ceEHF)*(1-cEHM)*dEH*deEH + (1+cEHF*ceEHF)*(1+cEHM*ceEHM),
                             (1+cEHF*ceEHF)*(1-cEHM)*dEH*(1-deEH)*drEH + (cEHF*(1-ceEHF)*crEHF)*(1+cEHM*ceEHM) + (1+cEHF*ceEHF)*(cEHM*(1-ceEHM)*crEHM),
                             (1+cEHF*ceEHF)*(1-cEHM)*dEH*(1-deEH)*(1-drEH) + (cEHF*(1-ceEHF)*(1-crEHF))*(1+cEHM*ceEHM) + (1+cEHF*ceEHF)*(cEHM*(1-ceEHM)*(1-crEHM)),
                             (cEHF*(1-ceEHF)*crEHF)*(1-cEHM)*dEH*drEH + (cEHF*(1-ceEHF)*crEHF)*(cEHM*(1-ceEHM)*crEHM),
                             (cEHF*(1-ceEHF)*crEHF)*(1-cEHM)*dEH*(1-drEH) + (cEHF*(1-ceEHF)*(1-crEHF))*(1-cEHM)*dEH*drEH + (cEHF*(1-ceEHF)*(1-crEHF))*(cEHM*(1-ceEHM)*crEHM) + (cEHF*(1-ceEHF)*crEHF)*(cEHM*(1-ceEHM)*(1-crEHM)),
                             (cEHF*(1-ceEHF)*(1-crEHF))*(1-cEHM)*dEH*(1-drEH) + (cEHF*(1-ceEHF)*(1-crEHF))*(cEHM*(1-ceEHM)*(1-crEHM)))/4

  # females
  tMatrix['HE','WR', ] <- c( 0, (1-cEHF)*(1-dEH), (1+cEHF*ceEHF)*(1-dEH), (cEHF*(1-ceEHF)*crEHF)*(1-dEH), (cEHF*(1-ceEHF)*(1-crEHF))*(1-dEH), 0, 0, (1-cEHF)*dEH*drEH + (1-cEHF), (1-cEHF)*dEH*(1-drEH), (1+cEHF*ceEHF)*dEH*deEH,
                             (1+cEHF*ceEHF)*dEH*(1-deEH)*drEH + (1+cEHF*ceEHF), (1+cEHF*ceEHF)*dEH*(1-deEH)*(1-drEH), (cEHF*(1-ceEHF)*crEHF)*dEH*drEH + (cEHF*(1-ceEHF)*crEHF),
                             (cEHF*(1-ceEHF)*crEHF)*dEH*(1-drEH) + (cEHF*(1-ceEHF)*(1-crEHF))*dEH*drEH + (cEHF*(1-ceEHF)*(1-crEHF)), (cEHF*(1-ceEHF)*(1-crEHF))*dEH*(1-drEH))/4
  tMatrix['HE','WB', ] <- c( 0, (1-cEHF)*(1-dEH), (1+cEHF*ceEHF)*(1-dEH), (cEHF*(1-ceEHF)*crEHF)*(1-dEH), (cEHF*(1-ceEHF)*(1-crEHF))*(1-dEH), 0, 0, (1-cEHF)*dEH*drEH, (1-cEHF)*dEH*(1-drEH) + (1-cEHF), (1+cEHF*ceEHF)*dEH*deEH,
                             (1+cEHF*ceEHF)*dEH*(1-deEH)*drEH, (1+cEHF*ceEHF)*dEH*(1-deEH)*(1-drEH) + (1+cEHF*ceEHF), (cEHF*(1-ceEHF)*crEHF)*dEH*drEH,
                             (cEHF*(1-ceEHF)*crEHF)*dEH*(1-drEH) + (cEHF*(1-ceEHF)*(1-crEHF))*dEH*drEH + (cEHF*(1-ceEHF)*crEHF), (cEHF*(1-ceEHF)*(1-crEHF))*dEH*(1-drEH) + (cEHF*(1-ceEHF)*(1-crEHF)))/4
  tMatrix['HE','HH', c('HH', 'HE', 'HR', 'HB', 'EE',
                       'ER', 'EB', 'RR', 'RB', 'BB')] <- c( (1-cEHF)*(1-dEH), (1+cEHF*ceEHF)*(1-dEH), (1-cEHF)*dEH*drEH + (cEHF*(1-ceEHF)*crEHF)*(1-dEH), (1-cEHF)*dEH*(1-drEH) + (cEHF*(1-ceEHF)*(1-crEHF))*(1-dEH),
                                                            (1+cEHF*ceEHF)*dEH*deEH, (1+cEHF*ceEHF)*dEH*(1-deEH)*drEH, (1+cEHF*ceEHF)*dEH*(1-deEH)*(1-drEH), (cEHF*(1-ceEHF)*crEHF)*dEH*drEH,
                                                            (cEHF*(1-ceEHF)*crEHF)*dEH*(1-drEH) + (cEHF*(1-ceEHF)*(1-crEHF))*dEH*drEH, (cEHF*(1-ceEHF)*(1-crEHF))*dEH*(1-drEH))/2
  tMatrix['HE','HR', c('HH', 'HE', 'HR', 'HB', 'EE',
                       'ER', 'EB', 'RR', 'RB', 'BB')] <- c( (1-cEHF)*(1-dEH), (1+cEHF*ceEHF)*(1-dEH), (1-cEHF)*dEH*drEH + (cEHF*(1-ceEHF)*crEHF)*(1-dEH) + (1-cEHF), (1-cEHF)*dEH*(1-drEH) + (cEHF*(1-ceEHF)*(1-crEHF))*(1-dEH),
                                                            (1+cEHF*ceEHF)*dEH*deEH, (1+cEHF*ceEHF)*dEH*(1-deEH)*drEH + (1+cEHF*ceEHF), (1+cEHF*ceEHF)*dEH*(1-deEH)*(1-drEH), (cEHF*(1-ceEHF)*crEHF)*dEH*drEH + (cEHF*(1-ceEHF)*crEHF),
                                                            (cEHF*(1-ceEHF)*crEHF)*dEH*(1-drEH) + (cEHF*(1-ceEHF)*(1-crEHF))*dEH*drEH + (cEHF*(1-ceEHF)*(1-crEHF)), (cEHF*(1-ceEHF)*(1-crEHF))*dEH*(1-drEH))/4
  tMatrix['HE','HB', c('HH', 'HE', 'HR', 'HB', 'EE',
                       'ER', 'EB', 'RR', 'RB', 'BB')] <- c( (1-cEHF)*(1-dEH), (1+cEHF*ceEHF)*(1-dEH), (1-cEHF)*dEH*drEH + (cEHF*(1-ceEHF)*crEHF)*(1-dEH), (1-cEHF)*dEH*(1-drEH) + (cEHF*(1-ceEHF)*(1-crEHF))*(1-dEH) + (1-cEHF),
                                                            (1+cEHF*ceEHF)*dEH*deEH, (1+cEHF*ceEHF)*dEH*(1-deEH)*drEH, (1+cEHF*ceEHF)*dEH*(1-deEH)*(1-drEH) + (1+cEHF*ceEHF), (cEHF*(1-ceEHF)*crEHF)*dEH*drEH,
                                                            (cEHF*(1-ceEHF)*crEHF)*dEH*(1-drEH) + (cEHF*(1-ceEHF)*(1-crEHF))*dEH*drEH + (cEHF*(1-ceEHF)*crEHF), (cEHF*(1-ceEHF)*(1-crEHF))*dEH*(1-drEH) + (cEHF*(1-ceEHF)*(1-crEHF)))/4
  tMatrix['HE','EE', c('HE', 'EE', 'ER', 'EB')] <- c( 1-cEHF, 1+cEHF*ceEHF, cEHF*(1-ceEHF)*crEHF, cEHF*(1-ceEHF)*(1-crEHF))/2
  tMatrix['HE','ER', c('HE', 'EE', 'ER', 'EB',
                       'HR', 'RR', 'RB')] <- c( 1-cEHF, 1+cEHF*ceEHF, cEHF*(1-ceEHF)*crEHF + 1+cEHF*ceEHF, cEHF*(1-ceEHF)*(1-crEHF),
                                                1-cEHF, cEHF*(1-ceEHF)*crEHF, cEHF*(1-ceEHF)*(1-crEHF))/4
  tMatrix['HE','EB', c('HE', 'EE', 'ER', 'EB',
                       'HB', 'RB', 'BB')] <- c( 1-cEHF, 1+cEHF*ceEHF, cEHF*(1-ceEHF)*crEHF, cEHF*(1-ceEHF)*(1-crEHF) + 1+cEHF*ceEHF,
                                                1-cEHF, cEHF*(1-ceEHF)*crEHF, cEHF*(1-ceEHF)*(1-crEHF))/4
  tMatrix['HE','RR', c('HR', 'ER', 'RR', 'RB')] <- c( 1-cEHF, 1+cEHF*ceEHF, cEHF*(1-ceEHF)*crEHF, cEHF*(1-ceEHF)*(1-crEHF))/2
  tMatrix['HE','RB', c('HR', 'ER', 'RR', 'RB',
                       'HB', 'EB', 'BB')] <- c( 1-cEHF, 1+cEHF*ceEHF, cEHF*(1-ceEHF)*crEHF, cEHF*(1-ceEHF)*(1-crEHF) + cEHF*(1-ceEHF)*crEHF,
                                                1-cEHF, 1+cEHF*ceEHF, cEHF*(1-ceEHF)*(1-crEHF))/4
  tMatrix['HE','BB', c('HB', 'EB', 'RB', 'BB')] <- c( 1-cEHF, 1+cEHF*ceEHF, cEHF*(1-ceEHF)*crEHF, cEHF*(1-ceEHF)*(1-crEHF))/2

  # males
  tMatrix['WR','HE', c('WH', 'WE', 'WR', 'WB',
                       'HR', 'ER', 'RR', 'RB')] <- c( 1-cEHM, 1+cEHM*ceEHM, cEHM*(1-ceEHM)*crEHM, cEHM*(1-ceEHM)*(1-crEHM),
                                                      1-cEHM, 1+cEHM*ceEHM, cEHM*(1-ceEHM)*crEHM, cEHM*(1-ceEHM)*(1-crEHM))/4
  tMatrix['WB','HE', c('WH', 'WE', 'WR', 'WB',
                       'HB', 'EB', 'RB', 'BB')] <- c( 1-cEHM, 1+cEHM*ceEHM, cEHM*(1-ceEHM)*crEHM, cEHM*(1-ceEHM)*(1-crEHM),
                                                      1-cEHM, 1+cEHM*ceEHM, cEHM*(1-ceEHM)*crEHM, cEHM*(1-ceEHM)*(1-crEHM))/4
  tMatrix['HH','HE', c('HH', 'HE', 'HR', 'HB')] <- c( 1-cEHM, 1+cEHM*ceEHM, cEHM*(1-ceEHM)*crEHM, cEHM*(1-ceEHM)*(1-crEHM))/2
  tMatrix['HR','HE', c('HH', 'HE', 'HR', 'HB',
                       'ER', 'RR', 'RB')] <- c( 1-cEHM, 1+cEHM*ceEHM, cEHM*(1-ceEHM)*crEHM + 1-cEHM, cEHM*(1-ceEHM)*(1-crEHM),
                                                1+cEHM*ceEHM, cEHM*(1-ceEHM)*crEHM, cEHM*(1-ceEHM)*(1-crEHM))/4
  tMatrix['HB','HE', c('HH', 'HE', 'HR', 'HB',
                       'EB', 'RB', 'BB')] <- c( 1-cEHM, 1+cEHM*ceEHM, cEHM*(1-ceEHM)*crEHM, cEHM*(1-ceEHM)*(1-crEHM) + 1-cEHM,
                                                1+cEHM*ceEHM, cEHM*(1-ceEHM)*crEHM, cEHM*(1-ceEHM)*(1-crEHM))/4
  tMatrix['EE','HE', c('HE', 'EE', 'ER', 'EB')] <- c( (1-cEHM)*(1-dEH), (1-cEHM)*dEH*deEH + (1+cEHM*ceEHM), (1-cEHM)*dEH*(1-deEH)*drEH + (cEHM*(1-ceEHM)*crEHM), (1-cEHM)*dEH*(1-deEH)*(1-drEH) + (cEHM*(1-ceEHM)*(1-crEHM)))/2
  tMatrix['ER','HE', c('HE', 'EE', 'ER', 'EB',
                       'HR', 'RR', 'RB')] <- c( (1-cEHM)*(1-dEH), (1-cEHM)*dEH*deEH + (1+cEHM*ceEHM), (1-cEHM)*dEH*(1-deEH)*drEH + (cEHM*(1-ceEHM)*crEHM) + (1+cEHM*ceEHM), (1-cEHM)*dEH*(1-deEH)*(1-drEH) + (cEHM*(1-ceEHM)*(1-crEHM)),
                                                (1-cEHM)*(1-dEH), (1-cEHM)*dEH*drEH + (cEHM*(1-ceEHM)*crEHM), (1-cEHM)*dEH*(1-drEH) + (cEHM*(1-ceEHM)*(1-crEHM)))/4
  tMatrix['EB','HE', c('HE', 'EE', 'ER', 'EB',
                       'HB', 'RB', 'BB')] <- c( (1-cEHM)*(1-dEH), (1-cEHM)*dEH*deEH + (1+cEHM*ceEHM), (1-cEHM)*dEH*(1-deEH)*drEH + (cEHM*(1-ceEHM)*crEHM), (1-cEHM)*dEH*(1-deEH)*(1-drEH) + (cEHM*(1-ceEHM)*(1-crEHM)) + (1+cEHM*ceEHM),
                                                (1-cEHM)*(1-dEH), (1-cEHM)*dEH*drEH + (cEHM*(1-ceEHM)*crEHM), (1-cEHM)*dEH*(1-drEH) + (cEHM*(1-ceEHM)*(1-crEHM)))/4
  tMatrix['RR','HE', c('HR', 'ER', 'RR', 'RB')] <- c( 1-cEHM, 1+cEHM*ceEHM, cEHM*(1-ceEHM)*crEHM, cEHM*(1-ceEHM)*(1-crEHM))/2
  tMatrix['RB','HE', c('HR', 'ER', 'RR', 'RB',
                       'HB', 'EB', 'BB')] <- c( 1-cEHM, 1+cEHM*ceEHM, cEHM*(1-ceEHM)*crEHM, cEHM*(1-ceEHM)*(1-crEHM) + cEHM*(1-ceEHM)*crEHM,
                                                1-cEHM, 1+cEHM*ceEHM, cEHM*(1-ceEHM)*(1-crEHM))/4
  tMatrix['BB','HE', c('HB', 'EB', 'RB', 'BB')] <- c( 1-cEHM, 1+cEHM*ceEHM, cEHM*(1-ceEHM)*crEHM, cEHM*(1-ceEHM)*(1-crEHM))/2


  ## protection from underflow errors
  tMatrix[tMatrix < .Machine$double.eps] <- 0


  ## initialize viability mask. No mother-specific death.
  viabilityMask <- array(data = 1, dim = c(size,size,size), dimnames = list(gtype, gtype, gtype))

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
    releaseType = "HH"
  ))

}
