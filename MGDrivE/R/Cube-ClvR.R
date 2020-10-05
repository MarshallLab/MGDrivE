###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDD: Mosquito Gene Drive Explorer - Discrete
#   ClvR (Cleave and Rescue)
#    - sex specific homing
#    - copy-number dependent female deposition
#   June 2020
#   jared_bennett@berkeley.edu
#
###############################################################################

#' Inheritance Cube: ClvR (Cleave and Rescue)
#'
#' Based on the Cleave-and-Rescue system of [Oberhofer](https://doi.org/10.1073/pnas.1921698117),
#' this is a 2-locus Cas9-based toxin-antidote system. The first locus carries the
#' Cas9, gRNAs, and a recoded copy of an essential gene. The second locus is the
#' targeted essential gene. This gene can be completely haplosufficient (\code{hSuf} = 1)
#' or completely haploinsufficient (\code{hSuf} = 0). It is assumed that having
#' 2 copies of the gene (be it wild-type at the second locus or recoded at the first)
#' confers complete viability.
#'
#'
#' @param cF Female cutting rate, one ClvR allele
#' @param crF Female functional resistance rate, one ClvR allele
#' @param ccF Female cutting rate, two ClvR alleles
#' @param ccrF Female functional resistance rate, two ClvR alleles
#'
#' @param cM Male cutting rate, one ClvR allele
#' @param crM Male functional resistance rate, one ClvR allele
#' @param ccM Male cutting rate, two ClvR alleles
#' @param ccrM Male functional resistance rate, two ClvR alleles
#'
#' @param dW Female deposition cutting rate
#' @param drW Female deposition functional resistance rate
#' @param ddW Female deposition (HH) cutting rate
#' @param ddrW Female deposition (HH) functional resistance rate
#' @param hSuf Haplosufficiency level, default is completely sufficient
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
cubeClvR <- function(cF = 1, crF = 0, ccF = cF, ccrF = crF,
                     cM = 1, crM = 0, ccM = cM, ccrM = crM,
                     dW = 0, drW = 0, ddW = dW, ddrW = drW, hSuf = 1,
                     eta = NULL, phi = NULL,
                     omega = NULL, xiF = NULL, xiM = NULL, s = NULL){

  # # for testing
  # testVec <- runif(n = 12)
  #
  # cF <- testVec[1]; crF <- testVec[2];
  # ccF <- testVec[9]; ccrF <- testVec[10];
  #
  # cM <- testVec[3]; crM <- testVec[4];
  # ccM <- testVec[11]; ccrM <- testVec[12];
  #
  # dW <- testVec[5]; drW <- testVec[6];
  # ddW <- testVec[7]; ddrW <- testVec[8];
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
  params <- c(cF,crF,ccF,ccrF,cM,crM,ccM,ccrM,dW,drW,ddW,ddrW,hSuf)
  if(any(params>1) || any(params<0)){
    stop("Parameters are rates, they must be between 0 and 1.")
  }


  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  gtype <- c('WWWW','WWWR','WWWB','WWRR','WWRB','WWBB',
             'WHWW','WHWR','WHWB','WHRR','WHRB','WHBB',
             'HHWW','HHWR','HHWB','HHRR','HHRB','HHBB')
  size <- length(gtype) #because I use it several times
  tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix


  ## fill tMatrix with symmetric probs
  tMatrix['WWWW','WWWW','WWWW'] <- 1

  tMatrix['WWWR','WWWW',c('WWWW','WWWR')] <- c(1,1)/2
  tMatrix['WWWR','WWWR',c('WWWW','WWWR','WWRR')] <- c(1,2,1)/4

  tMatrix['WWWB','WWWW',c('WWWW','WWWB')] <- c(1,1)/2
  tMatrix['WWWB','WWWR',c('WWWW','WWWR','WWWB','WWRB')] <- c(1,1,1,1)/4
  tMatrix['WWWB','WWWB',c('WWWW','WWWB','WWBB')] <- c(1,2,1)/4

  tMatrix['WWRR','WWWW',c('WWWR')] <- 1
  tMatrix['WWRR','WWWR',c('WWWR','WWRR')] <- c(1,1)/2
  tMatrix['WWRR','WWWB',c('WWWR','WWRB')] <- c(1,1)/2
  tMatrix['WWRR','WWRR',c('WWRR')] <- 1

  tMatrix['WWRB','WWWW',c('WWWR','WWWB')] <- c(1,1)/2
  tMatrix['WWRB','WWWR',c('WWWR','WWWB','WWRR','WWRB')] <- c(1,1,1,1)/4
  tMatrix['WWRB','WWWB',c('WWWR','WWWB','WWRB','WWBB')] <- c(1,1,1,1)/4
  tMatrix['WWRB','WWRR',c('WWRR','WWRB')] <- c(1,1)/2
  tMatrix['WWRB','WWRB',c('WWRR','WWRB','WWBB')] <- c(1,2,1)/4

  tMatrix['WWBB','WWWW','WWWB'] <- 1
  tMatrix['WWBB','WWWR',c('WWWB','WWRB')] <- c(1,1)/2
  tMatrix['WWBB','WWWB',c('WWWB','WWBB')] <- c(1,1)/2
  tMatrix['WWBB','WWRR','WWRB'] <- 1
  tMatrix['WWBB','WWRB',c('WWRB','WWBB')] <- c(1,1)/2
  tMatrix['WWBB','WWBB','WWBB'] <- 1

  tMatrix['WHRR','WWRR',c('WWRR','WHRR')] <- c(1,1)/2
  tMatrix['WHRR','WWRB',c('WWRR','WHRR','WWRB','WHRB')] <- c(1,1,1,1)/4
  tMatrix['WHRR','WWBB',c('WWRB','WHRB')] <- c(1,1)/2
  tMatrix['WHRR','WHRR',c('WWRR','WHRR','HHRR')] <- c(1,2,1)/4

  tMatrix['WHRB','WWRR',c('WWRR','WWRB','WHRR','WHRB')] <- c(1,1,1,1)/4
  tMatrix['WHRB','WWRB',c('WWRR','WWRB','WWBB',
                          'WHRR','WHRB','WHBB')] <- c(1,2,1,1,2,1)/8
  tMatrix['WHRB','WWBB',c('WWRB','WWBB','WHRB','WHBB')] <- c(1,1,1,1)/4
  tMatrix['WHRB','WHRR',c('WWRR','WWRB','WHRR',
                          'WHRB','HHRR','HHRB')] <- c(1,1,2,2,1,1)/8
  tMatrix['WHRB','WHRB',c('WWRR','WWRB','WWBB',
                          'WHRR','WHRB','WHBB',
                          'HHRR','HHRB','HHBB')] <- c(1,2,1, 2,4,2, 1,2,1)/16

  tMatrix['WHBB','WWRR',c('WWRB','WHRB')] <- c(1,1)/2
  tMatrix['WHBB','WWRB',c('WWRB','WWBB','WHRB','WHBB')] <- c(1,1,1,1)/4
  tMatrix['WHBB','WWBB',c('WWBB','WHBB')] <- c(1,1)/2
  tMatrix['WHBB','WHRR',c('WWRB','WHRB','HHRB')] <- c(1,2,1)/4
  tMatrix['WHBB','WHRB',c('WWRB','WWBB',
                          'WHRB','WHBB',
                          'HHRB','HHBB')] <- c(1,1,2,2,1,1)/8
  tMatrix['WHBB','WHBB',c('WWBB','WHBB','HHBB')] <- c(1,2,1)/4

  tMatrix['HHRR','WWRR','WHRR'] <- 1
  tMatrix['HHRR','WWRB',c('WHRR','WHRB')] <- c(1,1)/2
  tMatrix['HHRR','WWBB','WHRB'] <- 1
  tMatrix['HHRR','WHRR',c('WHRR','HHRR')] <- c(1,1)/2
  tMatrix['HHRR','WHRB',c('WHRR','WHRB','HHRR','HHRB')] <- c(1,1,1,1)/4
  tMatrix['HHRR','WHBB',c('WHRB','HHRB')] <- c(1,1)/2
  tMatrix['HHRR','HHRR','HHRR'] <- 1

  tMatrix['HHRB','WWRR',c('WHRR','WHRB')] <- c(1,1)/2
  tMatrix['HHRB','WWRB',c('WHRR','WHRB','WHBB')] <- c(1,2,1)/4
  tMatrix['HHRB','WWBB',c('WHRB','WHBB')] <- c(1,1)/2
  tMatrix['HHRB','WHRR',c('WHRR','WHRB','HHRR','HHRB')] <- c(1,1,1,1)/4
  tMatrix['HHRB','WHRB',c('WHRR','WHRB','WHBB',
                          'HHRR','HHRB','HHBB')] <- c(1,2,1, 1,2,1)/8
  tMatrix['HHRB','WHBB',c('WHRB','WHBB','HHRB','HHBB')] <- c(1,1,1,1)/4
  tMatrix['HHRB','HHRR',c('HHRR','HHRB')] <- c(1,1)/2
  tMatrix['HHRB','HHRB',c('HHRR','HHRB','HHBB')] <- c(1,2,1)/4


  tMatrix['HHBB','WWRR','WHRB'] <- 1
  tMatrix['HHBB','WWRB',c('WHRB','WHBB')] <- c(1,1)/2
  tMatrix['HHBB','WWBB','WHBB'] <- 1
  tMatrix['HHBB','WHRR',c('WHRB','HHRB')] <- c(1,1)/2
  tMatrix['HHBB','WHRB',c('WHRB','WHBB','HHRB','HHBB')] <- c(1,1,1,1)/4
  tMatrix['HHBB','WHBB',c('WHBB','HHBB')] <- c(1,1)/2
  tMatrix['HHBB','HHRR','HHRB'] <- 1
  tMatrix['HHBB','HHRB',c('HHRB','HHBB')] <- c(1,1)/2
  tMatrix['HHBB','HHBB','HHBB'] <- 1


  ## set the other half of the matrix that is symmetric
  # Boolean matrix for subsetting, used several times
  boolMat <- upper.tri(x = tMatrix[ , ,1], diag = FALSE)
  # loop over depth, set upper triangle
  for(z in 1:size){tMatrix[ , ,z][boolMat] <- t(tMatrix[ , ,z])[boolMat]}


  ## non-deposition stuff
  tMatrix['WWWW','WHRR',c('WWWR','WHWR')] <- c(1,1)/2
  tMatrix['WWWW','WHRB',c('WWWR','WWWB','WHWR','WHWB')] <- c(1,1,1,1)/4
  tMatrix['WWWW','WHBB',c('WWWB','WHWB')] <- c(1,1)/2
  tMatrix['WWWW','HHRR','WHWR'] <- 1
  tMatrix['WWWW','HHRB',c('WHWR','WHWB')] <- c(1,1)/2
  tMatrix['WWWW','HHBB','WHWB'] <- 1

  tMatrix['WWWR','WHRR',c('WWWR','WWRR','WHWR','WHRR')] <- c(1,1,1,1)/4
  tMatrix['WWWR','WHRB',c('WWWR','WWWB','WWRR','WWRB',
                          'WHWR','WHWB','WHRR','WHRB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['WWWR','WHBB',c('WWWB','WWRB','WHWB','WHRB')] <- c(1,1,1,1)/4
  tMatrix['WWWR','HHRR',c('WHWR','WHRR')] <- c(1,1)/2
  tMatrix['WWWR','HHRB',c('WHWR','WHWB','WHRR','WHRB')] <- c(1,1,1,1)/4
  tMatrix['WWWR','HHBB',c('WHWB','WHRB')] <- c(1,1)/2

  tMatrix['WWWB','WHRR',c('WWWR','WWRB','WHWR','WHRB')] <- c(1,1,1,1)/4
  tMatrix['WWWB','WHRB',c('WWWR','WWWB','WWRB','WWBB',
                          'WHWR','WHWB','WHRB','WHBB')] <- c(1,1,1,1,1,1,1,1)/8
  tMatrix['WWWB','WHBB',c('WWWB','WWBB','WHWB','WHBB')] <- c(1,1,1,1)/4
  tMatrix['WWWB','HHRR',c('WHWR','WHRB')] <- c(1,1)/2
  tMatrix['WWWB','HHRB',c('WHWR','WHWB','WHRB','WHBB')] <- c(1,1,1,1)/4
  tMatrix['WWWB','HHBB',c('WHWB','WHBB')] <- c(1,1)/2


  ## deposition stuff
  tMatrix['WHRR','WWWW',c('WWWR','WWRR','WWRB',
                          'WHWR','WHRR','WHRB')] <- c(1-dW,dW*drW,dW*(1-drW),
                                                      1-dW,dW*drW,dW*(1-drW))/2
  tMatrix['WHRR','WWWR',c('WWWR','WWRR','WWRB',
                          'WHWR','WHRR','WHRB')] <- c(1-dW,1+dW*drW,dW*(1-drW),
                                                      1-dW,1+dW*drW,dW*(1-drW))/4
  tMatrix['WHRR','WWWB',c('WWWR','WWRR','WWRB',
                          'WHWR','WHRR','WHRB')] <- c(1-dW,dW*drW,1+dW*(1-drW),
                                                      1-dW,dW*drW,1+dW*(1-drW))/4

  tMatrix['WHRB','WWWW',c('WWWR','WWRR','WWRB','WWWB','WWBB',
                          'WHWR','WHRR','WHRB','WHWB','WHBB')] <- c(1-dW,dW*drW,dW*(1-drW) + dW*drW,1-dW,dW*(1-drW),
                                                                    1-dW,dW*drW,dW*(1-drW) + dW*drW,1-dW,dW*(1-drW))/4
  tMatrix['WHRB','WWWR',c('WWWR','WWRR','WWRB','WWWB','WWBB',
                          'WHWR','WHRR','WHRB','WHWB','WHBB')] <- c(1-dW,1+dW*drW,dW*(1-drW) + 1+dW*drW,1-dW,dW*(1-drW),
                                                                    1-dW,1+dW*drW,dW*(1-drW) + 1+dW*drW,1-dW,dW*(1-drW))/8
  tMatrix['WHRB','WWWB',c('WWWR','WWRR','WWRB','WWWB','WWBB',
                          'WHWR','WHRR','WHRB','WHWB','WHBB')] <- c(1-dW,dW*drW,1+dW*(1-drW) + dW*drW,1-dW,1+dW*(1-drW),
                                                                    1-dW,dW*drW,1+dW*(1-drW) + dW*drW,1-dW,1+dW*(1-drW))/8

  tMatrix['WHBB','WWWW',c('WWWB','WWRB','WWBB',
                          'WHWB','WHRB','WHBB')] <- c(1-dW,dW*drW,dW*(1-drW),
                                                      1-dW,dW*drW,dW*(1-drW))/2
  tMatrix['WHBB','WWWR',c('WWWB','WWRB','WWBB',
                          'WHWB','WHRB','WHBB')] <- c(1-dW,1+dW*drW,dW*(1-drW),
                                                      1-dW,1+dW*drW,dW*(1-drW))/4
  tMatrix['WHBB','WWWB',c('WWWB','WWRB','WWBB',
                          'WHWB','WHRB','WHBB')] <- c(1-dW,dW*drW,1+dW*(1-drW),
                                                      1-dW,dW*drW,1+dW*(1-drW))/4

  tMatrix['HHRR','WWWW',c('WHWR','WHRR','WHRB')] <- c(1-ddW,ddW*ddrW,ddW*(1-ddrW))
  tMatrix['HHRR','WWWR',c('WHWR','WHRR','WHRB')] <- c(1-ddW,1+ddW*ddrW,ddW*(1-ddrW))/2
  tMatrix['HHRR','WWWB',c('WHWR','WHRR','WHRB')] <- c(1-ddW,ddW*ddrW,1+ddW*(1-ddrW))/2

  tMatrix['HHRB','WWWW',c('WHWR','WHRR','WHRB','WHWB','WHBB')] <- c(1-ddW,ddW*ddrW,ddW*(1-ddrW) + ddW*ddrW,1-ddW,ddW*(1-ddrW))/2
  tMatrix['HHRB','WWWR',c('WHWR','WHRR','WHRB','WHWB','WHBB')] <- c(1-ddW,1+ddW*ddrW,ddW*(1-ddrW) + 1+ddW*ddrW,1-ddW,ddW*(1-ddrW))/4
  tMatrix['HHRB','WWWB',c('WHWR','WHRR','WHRB','WHWB','WHBB')] <- c(1-ddW,ddW*ddrW,1+ddW*(1-ddrW) + ddW*ddrW,1-ddW,1+ddW*(1-ddrW))/4

  tMatrix['HHBB','WWWW',c('WHWB','WHRB','WHBB')] <- c(1-ddW,ddW*ddrW,ddW*(1-ddrW))
  tMatrix['HHBB','WWWR',c('WHWB','WHRB','WHBB')] <- c(1-ddW,1+ddW*ddrW,ddW*(1-ddrW))/2
  tMatrix['HHBB','WWWB',c('WHWB','WHRB','WHBB')] <- c(1-ddW,ddW*ddrW,1+ddW*(1-ddrW))/2


  ## Female specific homing (with any deposition stuff)
  eSix <- rep.int(x = 0, times = 6)
  wrbSix <- c((1-cF)*(1-dW), (1-cF)*dW*drW + cF*crF*(1-dW), (1-cF)*dW*(1-drW) + cF*(1-crF)*(1-dW),
              cF*crF*dW*drW, cF*crF*dW*(1-drW) + cF*(1-crF)*dW*drW, cF*(1-crF)*dW*(1-drW))
  tMatrix['WHWW','WWWW', ] <- c(wrbSix,wrbSix,eSix)/2

  wrbSix <- c((1-cF)*(1-dW), (1-cF)*(1+dW*drW) + cF*crF*(1-dW), (1-cF)*dW*(1-drW) + cF*(1-crF)*(1-dW),
              cF*crF*(1+dW*drW), cF*crF*dW*(1-drW) + cF*(1-crF)*(1+dW*drW), cF*(1-crF)*dW*(1-drW))
  tMatrix['WHWW','WWWR', ] <- c(wrbSix,wrbSix,eSix)/4

  wrbSix <- c((1-cF)*(1-dW), (1-cF)*dW*drW + cF*crF*(1-dW), (1-cF)*(1+dW*(1-drW)) + cF*(1-crF)*(1-dW),
              cF*crF*dW*drW, cF*crF*(1+dW*(1-drW)) + cF*(1-crF)*dW*drW, cF*(1-crF)*(1+dW*(1-drW)))
  tMatrix['WHWW','WWWB', ] <- c(wrbSix,wrbSix,eSix)/4

  wrbSix <- c(0, 1-cF, 0, cF*crF, cF*(1-crF), 0)
  tMatrix['WHWW','WWRR', ] <- c(wrbSix,wrbSix,eSix)/2
  tMatrix['WHWW','WHRR', ] <- c(wrbSix,2*wrbSix,wrbSix)/4
  tMatrix['WHWW','HHRR', ] <- c(eSix,wrbSix,wrbSix)/2

  wrbSix <- c(0, 1-cF, 1-cF, cF*crF, cF*crF + cF*(1-crF), cF*(1-crF))
  tMatrix['WHWW','WWRB', ] <- c(wrbSix,wrbSix,eSix)/4
  tMatrix['WHWW','WHRB', ] <- c(wrbSix,2*wrbSix,wrbSix)/8
  tMatrix['WHWW','HHRB', ] <- c(eSix,wrbSix,wrbSix)/4

  wrbSix <- c(0, 0, 1-cF, 0, cF*crF, cF*(1-crF))
  tMatrix['WHWW','WWBB', ] <- c(wrbSix,wrbSix,eSix)/2
  tMatrix['WHWW','WHBB', ] <- c(wrbSix,2*wrbSix,wrbSix)/4
  tMatrix['WHWW','HHBB', ] <- c(eSix,wrbSix,wrbSix)/2

  wrbSix <- c((1-cF)*(1-dW), (1-cF)*dW*drW + (1+cF*crF)*(1-dW), (1-cF)*dW*(1-drW) + cF*(1-crF)*(1-dW),
              (1+cF*crF)*dW*drW, (1+cF*crF)*dW*(1-drW) + cF*(1-crF)*dW*drW, cF*(1-crF)*dW*(1-drW))
  tMatrix['WHWR','WWWW', ] <- c(wrbSix,wrbSix,eSix)/4

  wrbSix <- c((1-cF)*(1-dW), (1-cF)*(1+dW*drW) + (1+cF*crF)*(1-dW), (1-cF)*dW*(1-drW) + cF*(1-crF)*(1-dW),
              (1+cF*crF)*(1+dW*drW), (1+cF*crF)*dW*(1-drW) + cF*(1-crF)*(1+dW*drW), cF*(1-crF)*dW*(1-drW))
  tMatrix['WHWR','WWWR', ] <- c(wrbSix,wrbSix,eSix)/8

  wrbSix <- c((1-cF)*(1-dW), (1-cF)*dW*drW + (1+cF*crF)*(1-dW), (1-cF)*(1+dW*(1-drW)) + cF*(1-crF)*(1-dW),
              (1+cF*crF)*dW*drW, (1+cF*crF)*(1+dW*(1-drW)) + cF*(1-crF)*dW*drW, cF*(1-crF)*(1+dW*(1-drW)))
  tMatrix['WHWR','WWWB', ] <- c(wrbSix,wrbSix,eSix)/8

  wrbSix <- c(0, 1-cF, 0, 1+cF*crF, cF*(1-crF), 0)
  tMatrix['WHWR','WWRR', ] <- c(wrbSix,wrbSix,eSix)/4
  tMatrix['WHWR','WHRR', ] <- c(wrbSix,2*wrbSix,wrbSix)/8
  tMatrix['WHWR','HHRR', ] <- c(eSix,wrbSix,wrbSix)/4

  wrbSix <- c(0, 1-cF, 1-cF, 1+cF*crF, 1+cF*crF + cF*(1-crF), cF*(1-crF))
  tMatrix['WHWR','WWRB', ] <- c(wrbSix,wrbSix,eSix)/8
  tMatrix['WHWR','WHRB', ] <- c(wrbSix,2*wrbSix,wrbSix)/16
  tMatrix['WHWR','HHRB', ] <- c(eSix,wrbSix,wrbSix)/8

  wrbSix <- c(0, 0, 1-cF, 0, 1+cF*crF, cF*(1-crF))
  tMatrix['WHWR','WWBB', ] <- c(wrbSix,wrbSix,eSix)/4
  tMatrix['WHWR','WHBB', ] <- c(wrbSix,2*wrbSix,wrbSix)/8
  tMatrix['WHWR','HHBB', ] <- c(eSix,wrbSix,wrbSix)/4

  wrbSix <- c((1-cF)*(1-dW), (1-cF)*dW*drW + cF*crF*(1-dW), (1-cF)*dW*(1-drW) + (1+cF*(1-crF))*(1-dW),
              cF*crF*dW*drW, cF*crF*dW*(1-drW) + (1+cF*(1-crF))*dW*drW, (1+cF*(1-crF))*dW*(1-drW))
  tMatrix['WHWB','WWWW', ] <- c(wrbSix,wrbSix,eSix)/4

  wrbSix <- c((1-cF)*(1-dW), (1-cF)*(1+dW*drW) + cF*crF*(1-dW), (1-cF)*dW*(1-drW) + (1+cF*(1-crF))*(1-dW),
              cF*crF*(1+dW*drW), cF*crF*dW*(1-drW) + (1+cF*(1-crF))*(1+dW*drW), (1+cF*(1-crF))*dW*(1-drW))
  tMatrix['WHWB','WWWR', ] <- c(wrbSix,wrbSix,eSix)/8

  wrbSix <- c((1-cF)*(1-dW), (1-cF)*dW*drW + cF*crF*(1-dW), (1-cF)*(1+dW*(1-drW)) + (1+cF*(1-crF))*(1-dW),
              cF*crF*dW*drW, cF*crF*(1+dW*(1-drW)) + (1+cF*(1-crF))*dW*drW, (1+cF*(1-crF))*(1+dW*(1-drW)))
  tMatrix['WHWB','WWWB', ] <- c(wrbSix,wrbSix,eSix)/8

  wrbSix <- c(0, 1-cF, 0, cF*crF, 1+cF*(1-crF), 0)
  tMatrix['WHWB','WWRR', ] <- c(wrbSix,wrbSix,eSix)/4
  tMatrix['WHWB','WHRR', ] <- c(wrbSix,2*wrbSix,wrbSix)/8
  tMatrix['WHWB','HHRR', ] <- c(eSix,wrbSix,wrbSix)/4

  wrbSix <- c(0, 1-cF, 1-cF, cF*crF, cF*crF + 1+cF*(1-crF), 1+cF*(1-crF))
  tMatrix['WHWB','WWRB', ] <- c(wrbSix,wrbSix,eSix)/8
  tMatrix['WHWB','WHRB', ] <- c(wrbSix,2*wrbSix,wrbSix)/16
  tMatrix['WHWB','HHRB', ] <- c(eSix,wrbSix,wrbSix)/8

  wrbSix <- c(0, 0, 1-cF, 0, cF*crF, 1+cF*(1-crF))
  tMatrix['WHWB','WWBB', ] <- c(wrbSix,wrbSix,eSix)/4
  tMatrix['WHWB','WHBB', ] <- c(wrbSix,2*wrbSix,wrbSix)/8
  tMatrix['WHWB','HHBB', ] <- c(eSix,wrbSix,wrbSix)/4

  wrbSix <- c((1-ccF)*(1-ddW), (1-ccF)*ddW*ddrW + ccF*ccrF*(1-ddW), (1-ccF)*ddW*(1-ddrW) + ccF*(1-ccrF)*(1-ddW),
              ccF*ccrF*ddW*ddrW, ccF*ccrF*ddW*(1-ddrW) + ccF*(1-ccrF)*ddW*ddrW, ccF*(1-ccrF)*ddW*(1-ddrW))
  tMatrix['HHWW','WWWW', ] <- c(eSix,wrbSix,eSix)

  wrbSix <- c((1-ccF)*(1-ddW), (1-ccF)*(1+ddW*ddrW) + ccF*ccrF*(1-ddW), (1-ccF)*ddW*(1-ddrW) + ccF*(1-ccrF)*(1-ddW),
              ccF*ccrF*(1+ddW*ddrW), ccF*ccrF*ddW*(1-ddrW) + ccF*(1-ccrF)*(1+ddW*ddrW), ccF*(1-ccrF)*ddW*(1-ddrW))
  tMatrix['HHWW','WWWR', ] <- c(eSix,wrbSix,eSix)/2

  wrbSix <- c((1-ccF)*(1-ddW), (1-ccF)*ddW*ddrW + ccF*ccrF*(1-ddW), (1-ccF)*(1+ddW*(1-ddrW)) + ccF*(1-ccrF)*(1-ddW),
              ccF*ccrF*ddW*ddrW, ccF*ccrF*(1+ddW*(1-ddrW)) + ccF*(1-ccrF)*ddW*ddrW, ccF*(1-ccrF)*(1+ddW*(1-ddrW)))
  tMatrix['HHWW','WWWB', ] <- c(eSix,wrbSix,eSix)/2

  wrbSix <- c(0, 1-ccF, 0, ccF*ccrF, ccF*(1-ccrF), 0)
  tMatrix['HHWW','WWRR', ] <- c(eSix,wrbSix,eSix)
  tMatrix['HHWW','WHRR', ] <- c(eSix,wrbSix,wrbSix)/2
  tMatrix['HHWW','HHRR', ] <- c(eSix,eSix,wrbSix)

  wrbSix <- c(0, 1-ccF, 1-ccF, ccF*ccrF, ccF*ccrF + ccF*(1-ccrF), ccF*(1-ccrF))
  tMatrix['HHWW','WWRB', ] <- c(eSix,wrbSix,eSix)/2
  tMatrix['HHWW','WHRB', ] <- c(eSix,wrbSix,wrbSix)/4
  tMatrix['HHWW','HHRB', ] <- c(eSix,eSix,wrbSix)/2

  wrbSix <- c(0, 0, 1-ccF, 0, ccF*ccrF, ccF*(1-ccrF))
  tMatrix['HHWW','WWBB', ] <- c(eSix,wrbSix,eSix)
  tMatrix['HHWW','WHBB', ] <- c(eSix,wrbSix,wrbSix)/2
  tMatrix['HHWW','HHBB', ] <- c(eSix,eSix,wrbSix)

  wrbSix <- c((1-ccF)*(1-ddW), (1-ccF)*ddW*ddrW + (1+ccF*ccrF)*(1-ddW), (1-ccF)*ddW*(1-ddrW) + ccF*(1-ccrF)*(1-ddW),
              (1+ccF*ccrF)*ddW*ddrW, (1+ccF*ccrF)*ddW*(1-ddrW) + ccF*(1-ccrF)*ddW*ddrW, ccF*(1-ccrF)*ddW*(1-ddrW))
  tMatrix['HHWR','WWWW', ] <- c(eSix,wrbSix,eSix)/2

  wrbSix <- c((1-ccF)*(1-ddW), (1-ccF)*(1+ddW*ddrW) + (1+ccF*ccrF)*(1-ddW), (1-ccF)*ddW*(1-ddrW) + ccF*(1-ccrF)*(1-ddW),
              (1+ccF*ccrF)*(1+ddW*ddrW), (1+ccF*ccrF)*ddW*(1-ddrW) + ccF*(1-ccrF)*(1+ddW*ddrW), ccF*(1-ccrF)*ddW*(1-ddrW))
  tMatrix['HHWR','WWWR', ] <- c(eSix,wrbSix,eSix)/4

  wrbSix <- c((1-ccF)*(1-ddW), (1-ccF)*ddW*ddrW + (1+ccF*ccrF)*(1-ddW), (1-ccF)*(1+ddW*(1-ddrW)) + ccF*(1-ccrF)*(1-ddW),
              (1+ccF*ccrF)*ddW*ddrW, (1+ccF*ccrF)*(1+ddW*(1-ddrW)) + ccF*(1-ccrF)*ddW*ddrW, ccF*(1-ccrF)*(1+ddW*(1-ddrW)))
  tMatrix['HHWR','WWWB', ] <- c(eSix,wrbSix,eSix)/4

  wrbSix <- c(0, 1-ccF, 0, (1+ccF*ccrF), ccF*(1-ccrF), 0)
  tMatrix['HHWR','WWRR', ] <- c(eSix,wrbSix,eSix)/2
  tMatrix['HHWR','WHRR', ] <- c(eSix,wrbSix,wrbSix)/4
  tMatrix['HHWR','HHRR', ] <- c(eSix,eSix,wrbSix)/2

  wrbSix <- c(0, 1-ccF, 1-ccF, (1+ccF*ccrF), (1+ccF*ccrF) + ccF*(1-ccrF), ccF*(1-ccrF))
  tMatrix['HHWR','WWRB', ] <- c(eSix,wrbSix,eSix)/4
  tMatrix['HHWR','WHRB', ] <- c(eSix,wrbSix,wrbSix)/8
  tMatrix['HHWR','HHRB', ] <- c(eSix,eSix,wrbSix)/4

  wrbSix <- c(0, 0, 1-ccF, 0, (1+ccF*ccrF), ccF*(1-ccrF))
  tMatrix['HHWR','WWBB', ] <- c(eSix,wrbSix,eSix)/2
  tMatrix['HHWR','WHBB', ] <- c(eSix,wrbSix,wrbSix)/4
  tMatrix['HHWR','HHBB', ] <- c(eSix,eSix,wrbSix)/2

  wrbSix <- c((1-ccF)*(1-ddW), (1-ccF)*ddW*ddrW + ccF*ccrF*(1-ddW), (1-ccF)*ddW*(1-ddrW) + (1+ccF*(1-ccrF))*(1-ddW),
              ccF*ccrF*ddW*ddrW, ccF*ccrF*ddW*(1-ddrW) + (1+ccF*(1-ccrF))*ddW*ddrW, (1+ccF*(1-ccrF))*ddW*(1-ddrW))
  tMatrix['HHWB','WWWW', ] <- c(eSix,wrbSix,eSix)/2

  wrbSix <- c((1-ccF)*(1-ddW), (1-ccF)*(1+ddW*ddrW) + ccF*ccrF*(1-ddW), (1-ccF)*ddW*(1-ddrW) + (1+ccF*(1-ccrF))*(1-ddW),
              ccF*ccrF*(1+ddW*ddrW), ccF*ccrF*ddW*(1-ddrW) + (1+ccF*(1-ccrF))*(1+ddW*ddrW), (1+ccF*(1-ccrF))*ddW*(1-ddrW))
  tMatrix['HHWB','WWWR', ] <- c(eSix,wrbSix,eSix)/4

  wrbSix <- c((1-ccF)*(1-ddW), (1-ccF)*ddW*ddrW + ccF*ccrF*(1-ddW), (1-ccF)*(1+ddW*(1-ddrW)) + (1+ccF*(1-ccrF))*(1-ddW),
              ccF*ccrF*ddW*ddrW, ccF*ccrF*(1+ddW*(1-ddrW)) + (1+ccF*(1-ccrF))*ddW*ddrW, (1+ccF*(1-ccrF))*(1+ddW*(1-ddrW)))
  tMatrix['HHWB','WWWB', ] <- c(eSix,wrbSix,eSix)/4

  wrbSix <- c(0, 1-ccF, 0, ccF*ccrF, (1+ccF*(1-ccrF)), 0)
  tMatrix['HHWB','WWRR', ] <- c(eSix,wrbSix,eSix)/2
  tMatrix['HHWB','WHRR', ] <- c(eSix,wrbSix,wrbSix)/4
  tMatrix['HHWB','HHRR', ] <- c(eSix,eSix,wrbSix)/2

  wrbSix <- c(0, 1-ccF, 1-ccF, ccF*ccrF, ccF*ccrF + (1+ccF*(1-ccrF)), (1+ccF*(1-ccrF)))
  tMatrix['HHWB','WWRB', ] <- c(eSix,wrbSix,eSix)/4
  tMatrix['HHWB','WHRB', ] <- c(eSix,wrbSix,wrbSix)/8
  tMatrix['HHWB','HHRB', ] <- c(eSix,eSix,wrbSix)/4

  wrbSix <- c(0, 0, 1-ccF, 0, ccF*ccrF, (1+ccF*(1-ccrF)))
  tMatrix['HHWB','WWBB', ] <- c(eSix,wrbSix,eSix)/2
  tMatrix['HHWB','WHBB', ] <- c(eSix,wrbSix,wrbSix)/4
  tMatrix['HHWB','HHBB', ] <- c(eSix,eSix,wrbSix)/2


  ## male specific homing (with any deposition stuff)
  wrbSix <- c(1-cM, cM*crM, cM*(1-crM), 0, 0, 0)
  tMatrix['WWWW','WHWW', ] <- c(wrbSix,wrbSix,eSix)/2
  wrbSix <- c(1-ccM, ccM*ccrM, ccM*(1-ccrM), 0, 0, 0)
  tMatrix['WWWW','HHWW', ] <- c(eSix,wrbSix,eSix)

  wrbSix <- c(1-cM, 1+cM*crM, cM*(1-crM), 0, 0, 0)
  tMatrix['WWWW','WHWR', ] <- c(wrbSix,wrbSix,eSix)/4
  wrbSix <- c(1-ccM, 1+ccM*ccrM, ccM*(1-ccrM), 0, 0, 0)
  tMatrix['WWWW','HHWR', ] <- c(eSix,wrbSix,eSix)/2

  wrbSix <- c(1-cM, cM*crM, 1+cM*(1-crM), 0, 0, 0)
  tMatrix['WWWW','WHWB', ] <- c(wrbSix,wrbSix,eSix)/4
  wrbSix <- c(1-ccM, ccM*ccrM, 1+ccM*(1-ccrM), 0, 0, 0)
  tMatrix['WWWW','HHWB', ] <- c(eSix,wrbSix,eSix)/2

  wrbSix <- c(1-cM, cM*crM + 1-cM, cM*(1-crM), cM*crM, cM*(1-crM), 0)
  tMatrix['WWWR','WHWW', ] <- c(wrbSix,wrbSix,eSix)/4
  wrbSix <- c(1-ccM, ccM*ccrM + 1-ccM, ccM*(1-ccrM), ccM*ccrM, ccM*(1-ccrM), 0)
  tMatrix['WWWR','HHWW', ] <- c(eSix,wrbSix,eSix)/2

  wrbSix <- c(1-cM, 1+cM*crM + 1-cM, cM*(1-crM), 1+cM*crM, cM*(1-crM), 0)
  tMatrix['WWWR','WHWR', ] <- c(wrbSix,wrbSix,eSix)/8
  wrbSix <- c(1-ccM, 1+ccM*ccrM + 1-ccM, ccM*(1-ccrM), 1+ccM*ccrM, ccM*(1-ccrM), 0)
  tMatrix['WWWR','HHWR', ] <- c(eSix,wrbSix,eSix)/4

  wrbSix <- c(1-cM, cM*crM + 1-cM, 1+cM*(1-crM), cM*crM, 1+cM*(1-crM), 0)
  tMatrix['WWWR','WHWB', ] <- c(wrbSix,wrbSix,eSix)/8
  wrbSix <- c(1-ccM, ccM*ccrM + 1-ccM, 1+ccM*(1-ccrM), ccM*ccrM, 1+ccM*(1-ccrM), 0)
  tMatrix['WWWR','HHWB', ] <- c(eSix,wrbSix,eSix)/4

  wrbSix <- c(1-cM, cM*crM, cM*(1-crM) + 1-cM, 0, cM*crM, cM*(1-crM))
  tMatrix['WWWB','WHWW', ] <- c(wrbSix,wrbSix,eSix)/4
  wrbSix <- c(1-ccM, ccM*ccrM, ccM*(1-ccrM) + 1-ccM, 0, ccM*ccrM, ccM*(1-ccrM))
  tMatrix['WWWB','HHWW', ] <- c(eSix,wrbSix,eSix)/2

  wrbSix <- c(1-cM, 1+cM*crM, cM*(1-crM) + 1-cM, 0, 1+cM*crM, cM*(1-crM))
  tMatrix['WWWB','WHWR', ] <- c(wrbSix,wrbSix,eSix)/8
  wrbSix <- c(1-ccM, 1+ccM*ccrM, ccM*(1-ccrM) + 1-ccM, 0, 1+ccM*ccrM, ccM*(1-ccrM))
  tMatrix['WWWB','HHWR', ] <- c(eSix,wrbSix,eSix)/4

  wrbSix <- c(1-cM, cM*crM, 1+cM*(1-crM) + 1-cM, 0, cM*crM, 1+cM*(1-crM))
  tMatrix['WWWB','WHWB', ] <- c(wrbSix,wrbSix,eSix)/8
  wrbSix <- c(1-ccM, ccM*ccrM, 1+ccM*(1-ccrM) + 1-ccM, 0, ccM*ccrM, 1+ccM*(1-ccrM))
  tMatrix['WWWB','HHWB', ] <- c(eSix,wrbSix,eSix)/4

  wrbSix <- c(0, 1-cM, 0, cM*crM, cM*(1-crM), 0)
  tMatrix['WWRR','WHWW', ] <- c(wrbSix,wrbSix,eSix)/2
  wrbSix <- c(0, 1-ccM, 0, ccM*ccrM, ccM*(1-ccrM), 0)
  tMatrix['WWRR','HHWW', ] <- c(eSix,wrbSix,eSix)

  wrbSix <- c(0, 1-cM, 0, 1+cM*crM, cM*(1-crM), 0)
  tMatrix['WWRR','WHWR', ] <- c(wrbSix,wrbSix,eSix)/4
  wrbSix <- c(0, 1-ccM, 0, 1+ccM*ccrM, ccM*(1-ccrM), 0)
  tMatrix['WWRR','HHWR', ] <- c(eSix,wrbSix,eSix)/2

  wrbSix <- c(0, 1-cM, 0, cM*crM, 1+cM*(1-crM), 0)
  tMatrix['WWRR','WHWB', ] <- c(wrbSix,wrbSix,eSix)/4
  wrbSix <- c(0, 1-ccM, 0, ccM*ccrM, 1+ccM*(1-ccrM), 0)
  tMatrix['WWRR','HHWB', ] <- c(eSix,wrbSix,eSix)/2

  wrbSix <- c(0, 1-cM, 1-cM, cM*crM, cM*(1-crM) + cM*crM, cM*(1-crM))
  tMatrix['WWRB','WHWW', ] <- c(wrbSix,wrbSix,eSix)/4
  wrbSix <- c(0, 1-ccM, 1-ccM, ccM*ccrM, ccM*(1-ccrM) + ccM*ccrM, ccM*(1-ccrM))
  tMatrix['WWRB','HHWW', ] <- c(eSix,wrbSix,eSix)/2

  wrbSix <- c(0, 1-cM, 1-cM, 1+cM*crM, cM*(1-crM) + 1+cM*crM, cM*(1-crM))
  tMatrix['WWRB','WHWR', ] <- c(wrbSix,wrbSix,eSix)/8
  wrbSix <- c(0, 1-ccM, 1-ccM, 1+ccM*ccrM, ccM*(1-ccrM) + 1+ccM*ccrM, ccM*(1-ccrM))
  tMatrix['WWRB','HHWR', ] <- c(eSix,wrbSix,eSix)/4

  wrbSix <- c(0, 1-cM, 1-cM, cM*crM, 1+cM*(1-crM) + cM*crM, 1+cM*(1-crM))
  tMatrix['WWRB','WHWB', ] <- c(wrbSix,wrbSix,eSix)/8
  wrbSix <- c(0, 1-ccM, 1-ccM, ccM*ccrM, 1+ccM*(1-ccrM) + ccM*ccrM, 1+ccM*(1-ccrM))
  tMatrix['WWRB','HHWB', ] <- c(eSix,wrbSix,eSix)/4

  wrbSix <- c(0, 0, 1-cM, 0, cM*crM, cM*(1-crM))
  tMatrix['WWBB','WHWW', ] <- c(wrbSix,wrbSix,eSix)/2
  wrbSix <- c(0, 0, 1-ccM, 0, ccM*ccrM, ccM*(1-ccrM))
  tMatrix['WWBB','HHWW', ] <- c(eSix,wrbSix,eSix)

  wrbSix <- c(0, 0, 1-cM, 0, 1+cM*crM, cM*(1-crM))
  tMatrix['WWBB','WHWR', ] <- c(wrbSix,wrbSix,eSix)/4
  wrbSix <- c(0, 0, 1-ccM, 0, 1+ccM*ccrM, ccM*(1-ccrM))
  tMatrix['WWBB','HHWR', ] <- c(eSix,wrbSix,eSix)/2

  wrbSix <- c(0, 0, 1-cM, 0, cM*crM, 1+cM*(1-crM))
  tMatrix['WWBB','WHWB', ] <- c(wrbSix,wrbSix,eSix)/4
  wrbSix <- c(0, 0, 1-ccM, 0, ccM*ccrM, 1+ccM*(1-ccrM))
  tMatrix['WWBB','HHWB', ] <- c(eSix,wrbSix,eSix)/2

  wrbSix <- c(0, (1-cM)*(1-dW), 0, cM*crM + (1-cM)*dW*drW, cM*(1-crM) + (1-cM)*dW*(1-drW), 0)
  tMatrix['WHRR','WHWW', ] <- c(wrbSix,2*wrbSix,wrbSix)/4
  wrbSix <- c(0, (1-ccM)*(1-dW), 0, ccM*ccrM + (1-ccM)*dW*drW, ccM*(1-ccrM) + (1-ccM)*dW*(1-drW), 0)
  tMatrix['WHRR','HHWW', ] <- c(eSix,wrbSix,wrbSix)/2

  wrbSix <- c(0, (1-cM)*(1-dW), 0, 1+cM*crM + (1-cM)*dW*drW, cM*(1-crM) + (1-cM)*dW*(1-drW), 0)
  tMatrix['WHRR','WHWR', ] <- c(wrbSix,2*wrbSix,wrbSix)/8
  wrbSix <- c(0, (1-ccM)*(1-dW), 0, 1+ccM*ccrM + (1-ccM)*dW*drW, ccM*(1-ccrM) + (1-ccM)*dW*(1-drW), 0)
  tMatrix['WHRR','HHWR', ] <- c(eSix,wrbSix,wrbSix)/4

  wrbSix <- c(0, (1-cM)*(1-dW), 0, cM*crM + (1-cM)*dW*drW, 1+cM*(1-crM) + (1-cM)*dW*(1-drW), 0)
  tMatrix['WHRR','WHWB', ] <- c(wrbSix,2*wrbSix,wrbSix)/8
  wrbSix <- c(0, (1-ccM)*(1-dW), 0, ccM*ccrM + (1-ccM)*dW*drW, 1+ccM*(1-ccrM) + (1-ccM)*dW*(1-drW), 0)
  tMatrix['WHRR','HHWB', ] <- c(eSix,wrbSix,wrbSix)/4

  wrbSix <- c(0, (1-cM)*(1-dW), (1-cM)*(1-dW), cM*crM + (1-cM)*dW*drW,
              cM*(1-crM) + (1-cM)*dW*(1-drW) + cM*crM + (1-cM)*dW*drW, cM*(1-crM) + (1-cM)*dW*(1-drW))
  tMatrix['WHRB','WHWW', ] <- c(wrbSix,2*wrbSix,wrbSix)/8
  wrbSix <- c(0, (1-ccM)*(1-dW), (1-ccM)*(1-dW), ccM*ccrM + (1-ccM)*dW*drW,
              ccM*(1-ccrM) + (1-ccM)*dW*(1-drW) + ccM*ccrM + (1-ccM)*dW*drW, ccM*(1-ccrM) + (1-ccM)*dW*(1-drW))
  tMatrix['WHRB','HHWW', ] <- c(eSix,wrbSix,wrbSix)/4

  wrbSix <- c(0, (1-cM)*(1-dW), (1-cM)*(1-dW), 1+cM*crM + (1-cM)*dW*drW,
              cM*(1-crM) + (1-cM)*dW*(1-drW) + 1+cM*crM + (1-cM)*dW*drW, cM*(1-crM) + (1-cM)*dW*(1-drW))
  tMatrix['WHRB','WHWR', ] <- c(wrbSix,2*wrbSix,wrbSix)/16
  wrbSix <- c(0, (1-ccM)*(1-dW), (1-ccM)*(1-dW), 1+ccM*ccrM + (1-ccM)*dW*drW,
              ccM*(1-ccrM) + (1-ccM)*dW*(1-drW) + 1+ccM*ccrM + (1-ccM)*dW*drW, ccM*(1-ccrM) + (1-ccM)*dW*(1-drW))
  tMatrix['WHRB','HHWR', ] <- c(eSix,wrbSix,wrbSix)/8

  wrbSix <- c(0, (1-cM)*(1-dW), (1-cM)*(1-dW), cM*crM + (1-cM)*dW*drW,
              1+cM*(1-crM) + (1-cM)*dW*(1-drW) + cM*crM + (1-cM)*dW*drW, 1+cM*(1-crM) + (1-cM)*dW*(1-drW))
  tMatrix['WHRB','WHWB', ] <- c(wrbSix,2*wrbSix,wrbSix)/16
  wrbSix <- c(0, (1-ccM)*(1-dW), (1-ccM)*(1-dW), ccM*ccrM + (1-ccM)*dW*drW,
              1+ccM*(1-ccrM) + (1-ccM)*dW*(1-drW) + ccM*ccrM + (1-ccM)*dW*drW, 1+ccM*(1-ccrM) + (1-ccM)*dW*(1-drW))
  tMatrix['WHRB','HHWB', ] <- c(eSix,wrbSix,wrbSix)/8

  wrbSix <- c(0, 0, (1-cM)*(1-dW), 0, cM*crM + (1-cM)*dW*drW, cM*(1-crM) + (1-cM)*dW*(1-drW))
  tMatrix['WHBB','WHWW', ] <- c(wrbSix,2*wrbSix,wrbSix)/4
  wrbSix <- c(0, 0, (1-ccM)*(1-dW), 0, ccM*ccrM + (1-ccM)*dW*drW, ccM*(1-ccrM) + (1-ccM)*dW*(1-drW))
  tMatrix['WHBB','HHWW', ] <- c(eSix,wrbSix,wrbSix)/2

  wrbSix <- c(0, 0, (1-cM)*(1-dW), 0, 1+cM*crM + (1-cM)*dW*drW, cM*(1-crM) + (1-cM)*dW*(1-drW))
  tMatrix['WHBB','WHWR', ] <- c(wrbSix,2*wrbSix,wrbSix)/8
  wrbSix <- c(0, 0, (1-ccM)*(1-dW), 0, 1+ccM*ccrM + (1-ccM)*dW*drW, ccM*(1-ccrM) + (1-ccM)*dW*(1-drW))
  tMatrix['WHBB','HHWR', ] <- c(eSix,wrbSix,wrbSix)/4

  wrbSix <- c(0, 0, (1-cM)*(1-dW), 0, cM*crM + (1-cM)*dW*drW, 1+cM*(1-crM) + (1-cM)*dW*(1-drW))
  tMatrix['WHBB','WHWB', ] <- c(wrbSix,2*wrbSix,wrbSix)/8
  wrbSix <- c(0, 0, (1-ccM)*(1-dW), 0, ccM*ccrM + (1-ccM)*dW*drW, 1+ccM*(1-ccrM) + (1-ccM)*dW*(1-drW))
  tMatrix['WHBB','HHWB', ] <- c(eSix,wrbSix,wrbSix)/4

  wrbSix <- c(0, (1-cM)*(1-ddW), 0, cM*crM + (1-cM)*ddW*ddrW, cM*(1-crM) + (1-cM)*ddW*(1-ddrW), 0)
  tMatrix['HHRR','WHWW', ] <- c(eSix,wrbSix,wrbSix)/2
  wrbSix <- c(0, (1-ccM)*(1-ddW), 0, ccM*ccrM + (1-ccM)*ddW*ddrW, ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW), 0)
  tMatrix['HHRR','HHWW', ] <- c(eSix,eSix,wrbSix)

  wrbSix <- c(0, (1-cM)*(1-ddW), 0, 1+cM*crM + (1-cM)*ddW*ddrW, cM*(1-crM) + (1-cM)*ddW*(1-ddrW), 0)
  tMatrix['HHRR','WHWR', ] <- c(eSix,wrbSix,wrbSix)/4
  wrbSix <- c(0, (1-ccM)*(1-ddW), 0, 1+ccM*ccrM + (1-ccM)*ddW*ddrW, ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW), 0)
  tMatrix['HHRR','HHWR', ] <- c(eSix,eSix,wrbSix)/2

  wrbSix <- c(0, (1-cM)*(1-ddW), 0, cM*crM + (1-cM)*ddW*ddrW, 1+cM*(1-crM) + (1-cM)*ddW*(1-ddrW), 0)
  tMatrix['HHRR','WHWB', ] <- c(eSix,wrbSix,wrbSix)/4
  wrbSix <- c(0, (1-ccM)*(1-ddW), 0, ccM*ccrM + (1-ccM)*ddW*ddrW, 1+ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW), 0)
  tMatrix['HHRR','HHWB', ] <- c(eSix,eSix,wrbSix)/2

  wrbSix <- c(0, (1-cM)*(1-ddW), (1-cM)*(1-ddW), cM*crM + (1-cM)*ddW*ddrW,
              cM*(1-crM) + (1-cM)*ddW*(1-ddrW) + cM*crM + (1-cM)*ddW*ddrW, cM*(1-crM) + (1-cM)*ddW*(1-ddrW))
  tMatrix['HHRB','WHWW', ] <- c(eSix,wrbSix,wrbSix)/4
  wrbSix <- c(0, (1-ccM)*(1-ddW), (1-ccM)*(1-ddW), ccM*ccrM + (1-ccM)*ddW*ddrW,
              ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW) + ccM*ccrM + (1-ccM)*ddW*ddrW, ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW))
  tMatrix['HHRB','HHWW', ] <- c(eSix,eSix,wrbSix)/2

  wrbSix <- c(0, (1-cM)*(1-ddW), (1-cM)*(1-ddW), 1+cM*crM + (1-cM)*ddW*ddrW,
              cM*(1-crM) + (1-cM)*ddW*(1-ddrW) + 1+cM*crM + (1-cM)*ddW*ddrW, cM*(1-crM) + (1-cM)*ddW*(1-ddrW))
  tMatrix['HHRB','WHWR', ] <- c(eSix,wrbSix,wrbSix)/8
  wrbSix <- c(0, (1-ccM)*(1-ddW), (1-ccM)*(1-ddW), 1+ccM*ccrM + (1-ccM)*ddW*ddrW,
              ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW) + 1+ccM*ccrM + (1-ccM)*ddW*ddrW, ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW))
  tMatrix['HHRB','HHWR', ] <- c(eSix,eSix,wrbSix)/4

  wrbSix <- c(0, (1-cM)*(1-ddW), (1-cM)*(1-ddW), cM*crM + (1-cM)*ddW*ddrW,
              1+cM*(1-crM) + (1-cM)*ddW*(1-ddrW) + cM*crM + (1-cM)*ddW*ddrW, 1+cM*(1-crM) + (1-cM)*ddW*(1-ddrW))
  tMatrix['HHRB','WHWB', ] <- c(eSix,wrbSix,wrbSix)/8
  wrbSix <- c(0, (1-ccM)*(1-ddW), (1-ccM)*(1-ddW), ccM*ccrM + (1-ccM)*ddW*ddrW,
              1+ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW) + ccM*ccrM + (1-ccM)*ddW*ddrW, 1+ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW))
  tMatrix['HHRB','HHWB', ] <- c(eSix,eSix,wrbSix)/4

  wrbSix <- c(0, 0, (1-cM)*(1-ddW), 0, cM*crM + (1-cM)*ddW*ddrW, cM*(1-crM) + (1-cM)*ddW*(1-ddrW))
  tMatrix['HHBB','WHWW', ] <- c(eSix,wrbSix,wrbSix)/2
  wrbSix <- c(0, 0, (1-ccM)*(1-ddW), 0, ccM*ccrM + (1-ccM)*ddW*ddrW, ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW))
  tMatrix['HHBB','HHWW', ] <- c(eSix,eSix,wrbSix)

  wrbSix <- c(0, 0, (1-cM)*(1-ddW), 0, 1+cM*crM + (1-cM)*ddW*ddrW, cM*(1-crM) + (1-cM)*ddW*(1-ddrW))
  tMatrix['HHBB','WHWR', ] <- c(eSix,wrbSix,wrbSix)/4
  wrbSix <- c(0, 0, (1-ccM)*(1-ddW), 0, 1+ccM*ccrM + (1-ccM)*ddW*ddrW, ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW))
  tMatrix['HHBB','HHWR', ] <- c(eSix,eSix,wrbSix)/2

  wrbSix <- c(0, 0, (1-cM)*(1-ddW), 0, cM*crM + (1-cM)*ddW*ddrW, 1+cM*(1-crM) + (1-cM)*ddW*(1-ddrW))
  tMatrix['HHBB','WHWB', ] <- c(eSix,wrbSix,wrbSix)/4
  wrbSix <- c(0, 0, (1-ccM)*(1-ddW), 0, ccM*ccrM + (1-ccM)*ddW*ddrW, 1+ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW))
  tMatrix['HHBB','HHWB', ] <- c(eSix,eSix,wrbSix)/2


  ## mixed homing/deposition stuff
  wrbSix <- c((1-cF)*((1-cM)*(1-dW)), (1-cF)*(cM*crM + (1-cM)*dW*drW) + cF*crF*((1-cM)*(1-dW)),
              (1-cF)*(cM*(1-crM) + (1-cM)*dW*(1-drW)) + cF*(1-crF)*((1-cM)*(1-dW)),
              cF*crF*(cM*crM + (1-cM)*dW*drW),
              cF*crF*(cM*(1-crM) + (1-cM)*dW*(1-drW)) + cF*(1-crF)*(cM*crM + (1-cM)*dW*drW),
              cF*(1-crF)*(cM*(1-crM) + (1-cM)*dW*(1-drW)))
  tMatrix['WHWW','WHWW', ] <- c(wrbSix,2*wrbSix,wrbSix)/4
  wrbSix <- c((1-cF)*((1-ccM)*(1-dW)), (1-cF)*(ccM*ccrM + (1-ccM)*dW*drW) + cF*crF*((1-ccM)*(1-dW)),
              (1-cF)*(ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)) + cF*(1-crF)*((1-ccM)*(1-dW)),
              cF*crF*(ccM*ccrM + (1-ccM)*dW*drW),
              cF*crF*(ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)) + cF*(1-crF)*(ccM*ccrM + (1-ccM)*dW*drW),
              cF*(1-crF)*(ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)))
  tMatrix['WHWW','HHWW', ] <- c(eSix,wrbSix,wrbSix)/2

  wrbSix <- c((1-cF)*((1-cM)*(1-dW)), (1-cF)*(1+cM*crM + (1-cM)*dW*drW) + cF*crF*((1-cM)*(1-dW)),
              (1-cF)*(cM*(1-crM) + (1-cM)*dW*(1-drW)) + cF*(1-crF)*((1-cM)*(1-dW)),
              cF*crF*(1+cM*crM + (1-cM)*dW*drW),
              cF*crF*(cM*(1-crM) + (1-cM)*dW*(1-drW)) + cF*(1-crF)*(1+cM*crM + (1-cM)*dW*drW),
              cF*(1-crF)*(cM*(1-crM) + (1-cM)*dW*(1-drW)))
  tMatrix['WHWW','WHWR', ] <- c(wrbSix,2*wrbSix,wrbSix)/8
  wrbSix <- c((1-cF)*((1-ccM)*(1-dW)), (1-cF)*(1+ccM*ccrM + (1-ccM)*dW*drW) + cF*crF*((1-ccM)*(1-dW)),
              (1-cF)*(ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)) + cF*(1-crF)*((1-ccM)*(1-dW)),
              cF*crF*(1+ccM*ccrM + (1-ccM)*dW*drW),
              cF*crF*(ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)) + cF*(1-crF)*(1+ccM*ccrM + (1-ccM)*dW*drW),
              cF*(1-crF)*(ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)))
  tMatrix['WHWW','HHWR', ] <- c(eSix,wrbSix,wrbSix)/4

  wrbSix <- c((1-cF)*((1-cM)*(1-dW)), (1-cF)*(cM*crM + (1-cM)*dW*drW) + cF*crF*((1-cM)*(1-dW)),
              (1-cF)*(1+cM*(1-crM) + (1-cM)*dW*(1-drW)) + cF*(1-crF)*((1-cM)*(1-dW)),
              cF*crF*(cM*crM + (1-cM)*dW*drW),
              cF*crF*(1+cM*(1-crM) + (1-cM)*dW*(1-drW)) + cF*(1-crF)*(cM*crM + (1-cM)*dW*drW),
              cF*(1-crF)*(1+cM*(1-crM) + (1-cM)*dW*(1-drW)))
  tMatrix['WHWW','WHWB', ] <- c(wrbSix,2*wrbSix,wrbSix)/8
  wrbSix <- c((1-cF)*((1-ccM)*(1-dW)), (1-cF)*(ccM*ccrM + (1-ccM)*dW*drW) + cF*crF*((1-ccM)*(1-dW)),
              (1-cF)*(1+ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)) + cF*(1-crF)*((1-ccM)*(1-dW)),
              cF*crF*(ccM*ccrM + (1-ccM)*dW*drW),
              cF*crF*(1+ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)) + cF*(1-crF)*(ccM*ccrM + (1-ccM)*dW*drW),
              cF*(1-crF)*(1+ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)))
  tMatrix['WHWW','HHWB', ] <- c(eSix,wrbSix,wrbSix)/4

  wrbSix <- c((1-cF)*((1-cM)*(1-dW)), (1-cF)*(cM*crM + (1-cM)*dW*drW) + (1+cF*crF)*((1-cM)*(1-dW)),
              (1-cF)*(cM*(1-crM) + (1-cM)*dW*(1-drW)) + cF*(1-crF)*((1-cM)*(1-dW)),
              (1+cF*crF)*(cM*crM + (1-cM)*dW*drW),
              (1+cF*crF)*(cM*(1-crM) + (1-cM)*dW*(1-drW)) + cF*(1-crF)*(cM*crM + (1-cM)*dW*drW),
              cF*(1-crF)*(cM*(1-crM) + (1-cM)*dW*(1-drW)))
  tMatrix['WHWR','WHWW', ] <- c(wrbSix,2*wrbSix,wrbSix)/8
  wrbSix <- c((1-cF)*((1-ccM)*(1-dW)), (1-cF)*(ccM*ccrM + (1-ccM)*dW*drW) + (1+cF*crF)*((1-ccM)*(1-dW)),
              (1-cF)*(ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)) + cF*(1-crF)*((1-ccM)*(1-dW)),
              (1+cF*crF)*(ccM*ccrM + (1-ccM)*dW*drW),
              (1+cF*crF)*(ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)) + cF*(1-crF)*(ccM*ccrM + (1-ccM)*dW*drW),
              cF*(1-crF)*(ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)))
  tMatrix['WHWR','HHWW', ] <- c(eSix,wrbSix,wrbSix)/4

  wrbSix <- c((1-cF)*((1-cM)*(1-dW)), (1-cF)*(1+cM*crM + (1-cM)*dW*drW) + (1+cF*crF)*((1-cM)*(1-dW)),
              (1-cF)*(cM*(1-crM) + (1-cM)*dW*(1-drW)) + cF*(1-crF)*((1-cM)*(1-dW)),
              (1+cF*crF)*(1+cM*crM + (1-cM)*dW*drW),
              (1+cF*crF)*(cM*(1-crM) + (1-cM)*dW*(1-drW)) + cF*(1-crF)*(1+cM*crM + (1-cM)*dW*drW),
              cF*(1-crF)*(cM*(1-crM) + (1-cM)*dW*(1-drW)))
  tMatrix['WHWR','WHWR', ] <- c(wrbSix,2*wrbSix,wrbSix)/16
  wrbSix <- c((1-cF)*((1-ccM)*(1-dW)), (1-cF)*(1+ccM*ccrM + (1-ccM)*dW*drW) + (1+cF*crF)*((1-ccM)*(1-dW)),
              (1-cF)*(ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)) + cF*(1-crF)*((1-ccM)*(1-dW)),
              (1+cF*crF)*(1+ccM*ccrM + (1-ccM)*dW*drW),
              (1+cF*crF)*(ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)) + cF*(1-crF)*(1+ccM*ccrM + (1-ccM)*dW*drW),
              cF*(1-crF)*(ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)))
  tMatrix['WHWR','HHWR', ] <- c(eSix,wrbSix,wrbSix)/8

  wrbSix <- c((1-cF)*((1-cM)*(1-dW)), (1-cF)*(cM*crM + (1-cM)*dW*drW) + (1+cF*crF)*((1-cM)*(1-dW)),
              (1-cF)*(1+cM*(1-crM) + (1-cM)*dW*(1-drW)) + cF*(1-crF)*((1-cM)*(1-dW)),
              (1+cF*crF)*(cM*crM + (1-cM)*dW*drW),
              (1+cF*crF)*(1+cM*(1-crM) + (1-cM)*dW*(1-drW)) + cF*(1-crF)*(cM*crM + (1-cM)*dW*drW),
              cF*(1-crF)*(1+cM*(1-crM) + (1-cM)*dW*(1-drW)))
  tMatrix['WHWR','WHWB', ] <- c(wrbSix,2*wrbSix,wrbSix)/16
  wrbSix <- c((1-cF)*((1-ccM)*(1-dW)), (1-cF)*(ccM*ccrM + (1-ccM)*dW*drW) + (1+cF*crF)*((1-ccM)*(1-dW)),
              (1-cF)*(1+ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)) + cF*(1-crF)*((1-ccM)*(1-dW)),
              (1+cF*crF)*(ccM*ccrM + (1-ccM)*dW*drW),
              (1+cF*crF)*(1+ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)) + cF*(1-crF)*(ccM*ccrM + (1-ccM)*dW*drW),
              cF*(1-crF)*(1+ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)))
  tMatrix['WHWR','HHWB', ] <- c(eSix,wrbSix,wrbSix)/8

  wrbSix <- c((1-cF)*((1-cM)*(1-dW)), (1-cF)*(cM*crM + (1-cM)*dW*drW) + cF*crF*((1-cM)*(1-dW)),
              (1-cF)*(cM*(1-crM) + (1-cM)*dW*(1-drW)) + (1+cF*(1-crF))*((1-cM)*(1-dW)),
              cF*crF*(cM*crM + (1-cM)*dW*drW),
              cF*crF*(cM*(1-crM) + (1-cM)*dW*(1-drW)) + (1+cF*(1-crF))*(cM*crM + (1-cM)*dW*drW),
              (1+cF*(1-crF))*(cM*(1-crM) + (1-cM)*dW*(1-drW)))
  tMatrix['WHWB','WHWW', ] <- c(wrbSix,2*wrbSix,wrbSix)/8
  wrbSix <- c((1-cF)*((1-ccM)*(1-dW)), (1-cF)*(ccM*ccrM + (1-ccM)*dW*drW) + cF*crF*((1-ccM)*(1-dW)),
              (1-cF)*(ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)) + (1+cF*(1-crF))*((1-ccM)*(1-dW)),
              cF*crF*(ccM*ccrM + (1-ccM)*dW*drW),
              cF*crF*(ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)) + (1+cF*(1-crF))*(ccM*ccrM + (1-ccM)*dW*drW),
              (1+cF*(1-crF))*(ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)))
  tMatrix['WHWB','HHWW', ] <- c(eSix,wrbSix,wrbSix)/4

  wrbSix <- c((1-cF)*((1-cM)*(1-dW)), (1-cF)*(1+cM*crM + (1-cM)*dW*drW) + cF*crF*((1-cM)*(1-dW)),
              (1-cF)*(cM*(1-crM) + (1-cM)*dW*(1-drW)) + (1+cF*(1-crF))*((1-cM)*(1-dW)),
              cF*crF*(1+cM*crM + (1-cM)*dW*drW),
              cF*crF*(cM*(1-crM) + (1-cM)*dW*(1-drW)) + (1+cF*(1-crF))*(1+cM*crM + (1-cM)*dW*drW),
              (1+cF*(1-crF))*(cM*(1-crM) + (1-cM)*dW*(1-drW)))
  tMatrix['WHWB','WHWR', ] <- c(wrbSix,2*wrbSix,wrbSix)/16
  wrbSix <- c((1-cF)*((1-ccM)*(1-dW)), (1-cF)*(1+ccM*ccrM + (1-ccM)*dW*drW) + cF*crF*((1-ccM)*(1-dW)),
              (1-cF)*(ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)) + (1+cF*(1-crF))*((1-ccM)*(1-dW)),
              cF*crF*(1+ccM*ccrM + (1-ccM)*dW*drW),
              cF*crF*(ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)) + (1+cF*(1-crF))*(1+ccM*ccrM + (1-ccM)*dW*drW),
              (1+cF*(1-crF))*(ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)))
  tMatrix['WHWB','HHWR', ] <- c(eSix,wrbSix,wrbSix)/8

  wrbSix <- c((1-cF)*((1-cM)*(1-dW)), (1-cF)*(cM*crM + (1-cM)*dW*drW) + cF*crF*((1-cM)*(1-dW)),
              (1-cF)*(1+cM*(1-crM) + (1-cM)*dW*(1-drW)) + (1+cF*(1-crF))*((1-cM)*(1-dW)),
              cF*crF*(cM*crM + (1-cM)*dW*drW),
              cF*crF*(1+cM*(1-crM) + (1-cM)*dW*(1-drW)) + (1+cF*(1-crF))*(cM*crM + (1-cM)*dW*drW),
              (1+cF*(1-crF))*(1+cM*(1-crM) + (1-cM)*dW*(1-drW)))
  tMatrix['WHWB','WHWB', ] <- c(wrbSix,2*wrbSix,wrbSix)/16
  wrbSix <- c((1-cF)*((1-ccM)*(1-dW)), (1-cF)*(ccM*ccrM + (1-ccM)*dW*drW) + cF*crF*((1-ccM)*(1-dW)),
              (1-cF)*(1+ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)) + (1+cF*(1-crF))*((1-ccM)*(1-dW)),
              cF*crF*(ccM*ccrM + (1-ccM)*dW*drW),
              cF*crF*(1+ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)) + (1+cF*(1-crF))*(ccM*ccrM + (1-ccM)*dW*drW),
              (1+cF*(1-crF))*(1+ccM*(1-ccrM) + (1-ccM)*dW*(1-drW)))
  tMatrix['WHWB','HHWB', ] <- c(eSix,wrbSix,wrbSix)/8

  wrbSix <- c((1-ccF)*((1-cM)*(1-ddW)), (1-ccF)*(cM*crM + (1-cM)*ddW*ddrW) + ccF*ccrF*((1-cM)*(1-ddW)),
              (1-ccF)*(cM*(1-crM) + (1-cM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*((1-cM)*(1-ddW)),
              ccF*ccrF*(cM*crM + (1-cM)*ddW*ddrW),
              ccF*ccrF*(cM*(1-crM) + (1-cM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*(cM*crM + (1-cM)*ddW*ddrW),
              ccF*(1-ccrF)*(cM*(1-crM) + (1-cM)*ddW*(1-ddrW)))
  tMatrix['HHWW','WHWW', ] <- c(eSix,wrbSix,wrbSix)/2
  wrbSix <- c((1-ccF)*((1-ccM)*(1-ddW)), (1-ccF)*(ccM*ccrM + (1-ccM)*ddW*ddrW) + ccF*ccrF*((1-ccM)*(1-ddW)),
              (1-ccF)*(ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*((1-ccM)*(1-ddW)),
              ccF*ccrF*(ccM*ccrM + (1-ccM)*ddW*ddrW),
              ccF*ccrF*(ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*(ccM*ccrM + (1-ccM)*ddW*ddrW),
              ccF*(1-ccrF)*(ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)))
  tMatrix['HHWW','HHWW', ] <- c(eSix,eSix,wrbSix)

  wrbSix <- c((1-ccF)*((1-cM)*(1-ddW)), (1-ccF)*(1+cM*crM + (1-cM)*ddW*ddrW) + ccF*ccrF*((1-cM)*(1-ddW)),
              (1-ccF)*(cM*(1-crM) + (1-cM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*((1-cM)*(1-ddW)),
              ccF*ccrF*(1+cM*crM + (1-cM)*ddW*ddrW),
              ccF*ccrF*(cM*(1-crM) + (1-cM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*(1+cM*crM + (1-cM)*ddW*ddrW),
              ccF*(1-ccrF)*(cM*(1-crM) + (1-cM)*ddW*(1-ddrW)))
  tMatrix['HHWW','WHWR', ] <- c(eSix,wrbSix,wrbSix)/4
  wrbSix <- c((1-ccF)*((1-ccM)*(1-ddW)), (1-ccF)*(1+ccM*ccrM + (1-ccM)*ddW*ddrW) + ccF*ccrF*((1-ccM)*(1-ddW)),
              (1-ccF)*(ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*((1-ccM)*(1-ddW)),
              ccF*ccrF*(1+ccM*ccrM + (1-ccM)*ddW*ddrW),
              ccF*ccrF*(ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*(1+ccM*ccrM + (1-ccM)*ddW*ddrW),
              ccF*(1-ccrF)*(ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)))
  tMatrix['HHWW','HHWR', ] <- c(eSix,eSix,wrbSix)/2

  wrbSix <- c((1-ccF)*((1-cM)*(1-ddW)), (1-ccF)*(cM*crM + (1-cM)*ddW*ddrW) + ccF*ccrF*((1-cM)*(1-ddW)),
              (1-ccF)*(1+cM*(1-crM) + (1-cM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*((1-cM)*(1-ddW)),
              ccF*ccrF*(cM*crM + (1-cM)*ddW*ddrW),
              ccF*ccrF*(1+cM*(1-crM) + (1-cM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*(cM*crM + (1-cM)*ddW*ddrW),
              ccF*(1-ccrF)*(1+cM*(1-crM) + (1-cM)*ddW*(1-ddrW)))
  tMatrix['HHWW','WHWB', ] <- c(eSix,wrbSix,wrbSix)/4
  wrbSix <- c((1-ccF)*((1-ccM)*(1-ddW)), (1-ccF)*(ccM*ccrM + (1-ccM)*ddW*ddrW) + ccF*ccrF*((1-ccM)*(1-ddW)),
              (1-ccF)*(1+ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*((1-ccM)*(1-ddW)),
              ccF*ccrF*(ccM*ccrM + (1-ccM)*ddW*ddrW),
              ccF*ccrF*(1+ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*(ccM*ccrM + (1-ccM)*ddW*ddrW),
              ccF*(1-ccrF)*(1+ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)))
  tMatrix['HHWW','HHWB', ] <- c(eSix,eSix,wrbSix)/2

  wrbSix <- c((1-ccF)*((1-cM)*(1-ddW)), (1-ccF)*(cM*crM + (1-cM)*ddW*ddrW) + (1+ccF*ccrF)*((1-cM)*(1-ddW)),
              (1-ccF)*(cM*(1-crM) + (1-cM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*((1-cM)*(1-ddW)),
              (1+ccF*ccrF)*(cM*crM + (1-cM)*ddW*ddrW),
              (1+ccF*ccrF)*(cM*(1-crM) + (1-cM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*(cM*crM + (1-cM)*ddW*ddrW),
              ccF*(1-ccrF)*(cM*(1-crM) + (1-cM)*ddW*(1-ddrW)))
  tMatrix['HHWR','WHWW', ] <- c(eSix,wrbSix,wrbSix)/4
  wrbSix <- c((1-ccF)*((1-ccM)*(1-ddW)), (1-ccF)*(ccM*ccrM + (1-ccM)*ddW*ddrW) + (1+ccF*ccrF)*((1-ccM)*(1-ddW)),
              (1-ccF)*(ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*((1-ccM)*(1-ddW)),
              (1+ccF*ccrF)*(ccM*ccrM + (1-ccM)*ddW*ddrW),
              (1+ccF*ccrF)*(ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*(ccM*ccrM + (1-ccM)*ddW*ddrW),
              ccF*(1-ccrF)*(ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)))
  tMatrix['HHWR','HHWW', ] <- c(eSix,eSix,wrbSix)/2

  wrbSix <- c((1-ccF)*((1-cM)*(1-ddW)), (1-ccF)*(1+cM*crM + (1-cM)*ddW*ddrW) + (1+ccF*ccrF)*((1-cM)*(1-ddW)),
              (1-ccF)*(cM*(1-crM) + (1-cM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*((1-cM)*(1-ddW)),
              (1+ccF*ccrF)*(1+cM*crM + (1-cM)*ddW*ddrW),
              (1+ccF*ccrF)*(cM*(1-crM) + (1-cM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*(1+cM*crM + (1-cM)*ddW*ddrW),
              ccF*(1-ccrF)*(cM*(1-crM) + (1-cM)*ddW*(1-ddrW)))
  tMatrix['HHWR','WHWR', ] <- c(eSix,wrbSix,wrbSix)/8
  wrbSix <- c((1-ccF)*((1-ccM)*(1-ddW)), (1-ccF)*(1+ccM*ccrM + (1-ccM)*ddW*ddrW) + (1+ccF*ccrF)*((1-ccM)*(1-ddW)),
              (1-ccF)*(ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*((1-ccM)*(1-ddW)),
              (1+ccF*ccrF)*(1+ccM*ccrM + (1-ccM)*ddW*ddrW),
              (1+ccF*ccrF)*(ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*(1+ccM*ccrM + (1-ccM)*ddW*ddrW),
              ccF*(1-ccrF)*(ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)))
  tMatrix['HHWR','HHWR', ] <- c(eSix,eSix,wrbSix)/4

  wrbSix <- c((1-ccF)*((1-cM)*(1-ddW)), (1-ccF)*(cM*crM + (1-cM)*ddW*ddrW) + (1+ccF*ccrF)*((1-cM)*(1-ddW)),
              (1-ccF)*(1+cM*(1-crM) + (1-cM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*((1-cM)*(1-ddW)),
              (1+ccF*ccrF)*(cM*crM + (1-cM)*ddW*ddrW),
              (1+ccF*ccrF)*(1+cM*(1-crM) + (1-cM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*(cM*crM + (1-cM)*ddW*ddrW),
              ccF*(1-ccrF)*(1+cM*(1-crM) + (1-cM)*ddW*(1-ddrW)))
  tMatrix['HHWR','WHWB', ] <- c(eSix,wrbSix,wrbSix)/8
  wrbSix <- c((1-ccF)*((1-ccM)*(1-ddW)), (1-ccF)*(ccM*ccrM + (1-ccM)*ddW*ddrW) + (1+ccF*ccrF)*((1-ccM)*(1-ddW)),
              (1-ccF)*(1+ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*((1-ccM)*(1-ddW)),
              (1+ccF*ccrF)*(ccM*ccrM + (1-ccM)*ddW*ddrW),
              (1+ccF*ccrF)*(1+ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)) + ccF*(1-ccrF)*(ccM*ccrM + (1-ccM)*ddW*ddrW),
              ccF*(1-ccrF)*(1+ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)))
  tMatrix['HHWR','HHWB', ] <- c(eSix,eSix,wrbSix)/4

  wrbSix <- c((1-ccF)*((1-cM)*(1-ddW)), (1-ccF)*(cM*crM + (1-cM)*ddW*ddrW) + ccF*ccrF*((1-cM)*(1-ddW)),
              (1-ccF)*(cM*(1-crM) + (1-cM)*ddW*(1-ddrW)) + (1+ccF*(1-ccrF))*((1-cM)*(1-ddW)),
              ccF*ccrF*(cM*crM + (1-cM)*ddW*ddrW),
              ccF*ccrF*(cM*(1-crM) + (1-cM)*ddW*(1-ddrW)) + (1+ccF*(1-ccrF))*(cM*crM + (1-cM)*ddW*ddrW),
              (1+ccF*(1-ccrF))*(cM*(1-crM) + (1-cM)*ddW*(1-ddrW)))
  tMatrix['HHWB','WHWW', ] <- c(eSix,wrbSix,wrbSix)/4
  wrbSix <- c((1-ccF)*((1-ccM)*(1-ddW)), (1-ccF)*(ccM*ccrM + (1-ccM)*ddW*ddrW) + ccF*ccrF*((1-ccM)*(1-ddW)),
              (1-ccF)*(ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)) + (1+ccF*(1-ccrF))*((1-ccM)*(1-ddW)),
              ccF*ccrF*(ccM*ccrM + (1-ccM)*ddW*ddrW),
              ccF*ccrF*(ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)) + (1+ccF*(1-ccrF))*(ccM*ccrM + (1-ccM)*ddW*ddrW),
              (1+ccF*(1-ccrF))*(ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)))
  tMatrix['HHWB','HHWW', ] <- c(eSix,eSix,wrbSix)/2

  wrbSix <- c((1-ccF)*((1-cM)*(1-ddW)), (1-ccF)*(1+cM*crM + (1-cM)*ddW*ddrW) + ccF*ccrF*((1-cM)*(1-ddW)),
              (1-ccF)*(cM*(1-crM) + (1-cM)*ddW*(1-ddrW)) + (1+ccF*(1-ccrF))*((1-cM)*(1-ddW)),
              ccF*ccrF*(1+cM*crM + (1-cM)*ddW*ddrW),
              ccF*ccrF*(cM*(1-crM) + (1-cM)*ddW*(1-ddrW)) + (1+ccF*(1-ccrF))*(1+cM*crM + (1-cM)*ddW*ddrW),
              (1+ccF*(1-ccrF))*(cM*(1-crM) + (1-cM)*ddW*(1-ddrW)))
  tMatrix['HHWB','WHWR', ] <- c(eSix,wrbSix,wrbSix)/8
  wrbSix <- c((1-ccF)*((1-ccM)*(1-ddW)), (1-ccF)*(1+ccM*ccrM + (1-ccM)*ddW*ddrW) + ccF*ccrF*((1-ccM)*(1-ddW)),
              (1-ccF)*(ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)) + (1+ccF*(1-ccrF))*((1-ccM)*(1-ddW)),
              ccF*ccrF*(1+ccM*ccrM + (1-ccM)*ddW*ddrW),
              ccF*ccrF*(ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)) + (1+ccF*(1-ccrF))*(1+ccM*ccrM + (1-ccM)*ddW*ddrW),
              (1+ccF*(1-ccrF))*(ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)))
  tMatrix['HHWB','HHWR', ] <- c(eSix,eSix,wrbSix)/4

  wrbSix <- c((1-ccF)*((1-cM)*(1-ddW)), (1-ccF)*(cM*crM + (1-cM)*ddW*ddrW) + ccF*ccrF*((1-cM)*(1-ddW)),
              (1-ccF)*(1+cM*(1-crM) + (1-cM)*ddW*(1-ddrW)) + (1+ccF*(1-ccrF))*((1-cM)*(1-ddW)),
              ccF*ccrF*(cM*crM + (1-cM)*ddW*ddrW),
              ccF*ccrF*(1+cM*(1-crM) + (1-cM)*ddW*(1-ddrW)) + (1+ccF*(1-ccrF))*(cM*crM + (1-cM)*ddW*ddrW),
              (1+ccF*(1-ccrF))*(1+cM*(1-crM) + (1-cM)*ddW*(1-ddrW)))
  tMatrix['HHWB','WHWB', ] <- c(eSix,wrbSix,wrbSix)/8
  wrbSix <- c((1-ccF)*((1-ccM)*(1-ddW)), (1-ccF)*(ccM*ccrM + (1-ccM)*ddW*ddrW) + ccF*ccrF*((1-ccM)*(1-ddW)),
              (1-ccF)*(1+ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)) + (1+ccF*(1-ccrF))*((1-ccM)*(1-ddW)),
              ccF*ccrF*(ccM*ccrM + (1-ccM)*ddW*ddrW),
              ccF*ccrF*(1+ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)) + (1+ccF*(1-ccrF))*(ccM*ccrM + (1-ccM)*ddW*ddrW),
              (1+ccF*(1-ccrF))*(1+ccM*(1-ccrM) + (1-ccM)*ddW*(1-ddrW)))
  tMatrix['HHWB','HHWB', ] <- c(eSix,eSix,wrbSix)/4


  ## protection from underflow errors
  tMatrix[tMatrix < .Machine$double.eps] <- 0


  ## initialize viability mask.
  # Anything with 2 copies of the gene is fine, be it wild-type or recoded
  #   HH at first locus all survive, WW/WR/RR at second locus all survive
  # Anything with 0 copies never survives
  #   WWBB is always dead
  # hSuf defines heterozygous sufficiency, ie is it dominant
  # aperm allows us to fill the array by depth, which isn't normally possible
  #   resize keeps dim names
  viabilityVec <- c(1,1,hSuf,1,hSuf,0, 1,1,1,1,1,hSuf, 1,1,1,1,1,1)
  viabilityMask <- aperm(a = array(data = viabilityVec, dim = c(size,size,size),
                                   dimnames = list(gtype, gtype, gtype)),
                         perm = c(3,2,1),resize = TRUE)


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
    releaseType = "HHWW"
  ))

}
