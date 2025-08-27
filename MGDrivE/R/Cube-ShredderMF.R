###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   Autosomal X-Shredder
#   Héctor Sanchez, Jared Bennett, Sean Wu, John Marshall
#   October 2018
#   jared_bennett@berkeley.edu
#
###############################################################################

#' Inheritance Cube: Autosomal X-Shredder
#'
#' This function creates an inheritance cube to model an autosomal X-Shredder construct.
#' This construct resides on an autosomal chromosome, and chops the X chromosome into
#' many pieces during gametogenesis, destroying the X chromosome. Thus, males
#' may only produce Y gametes and females can become sterile. \cr
#' This drive has 2 loci:
#'  * Locus 1, the autosomal locus, has 3 alleles:
#'    * W: Wild-type allele
#'    * A: Attacking allele, contains the shredder construct
#'    * B: Broken attacking allele, shredder construct is defunct
#'  * Locus 2, the sex locus, has 3 alleles:
#'    * X: Wild-type X allele
#'    * R: X-allele resistant to cleavage
#'    * Y: Wild-type Y allele
#'
#' @param cM Rate of X shredding in males (default is 1, complete shredding)
#' @param cF Rate of X shredding in females (default is 1, complete shredding)
#' @param crM Rate of resistance chromosome generation in males (default is 0)
#' @param crF Rate of resistance chromosome generation in females (default is 0)
#' @param cbM Rate of shredder construct breakdown in males (default is 0)
#' @param cbF Rate of shredder construct breakdown in females (default is 0)
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
cubeXShredderMF <- function(cM= 1, cF = 1, crM = 0, crF = 0,
                            cbM = 0, cbF = 0,
                            eta = NULL, phi = NULL, omega = NULL,
                            xiF = NULL, xiM = NULL, s = NULL){

  ## safety checks
  if(any(c(cM, cF, crM, crF, cbM, cbF)>1) || any(c(cM, cF, crM, crF, cbM, cbF)<0)){
    stop("Parameters are rates, they must be between 0 and 1.")
  }

  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  gtype <- c('WWXX','WAXX','WBXX','AAXX','ABXX','BBXX','WWXR','WAXR','WBXR',
             'AAXR','ABXR','BBXR','WWRR','WARR','WBRR','AARR','ABRR','BBRR',
             'WWXY','WAXY','WBXY','AAXY','ABXY','BBXY','WWRY','WARY','WBRY',
             'AARY','ABRY','BBRY')
  size <- length(gtype) #because I use it several times
  tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix


  ## fill tMatrix with probabilities
  tMatrix['WWXX','WWXY',c('WWXX','WWXY')] <- c( 1, 1)/2
  tMatrix['WWXX','WAXY',c('WWXX','WWXR','WWXY',
                          'WAXX','WAXR','WAXY',
                          'WBXX','WBXR','WBXY')] <- c(1-cM, cM*crM, 1,
                                                      (1-cbM)*(1-cM), (1-cbM)*cM*crM, (1-cbM),
                                                      cbM*(1-cM), cbM*cM*crM, cbM)/(4+2*cM*(crM-1))
  tMatrix['WWXX','WBXY',c('WWXX','WWXY','WBXX','WBXY')] <- c( 1, 1, 1, 1)/4
  tMatrix['WWXX','AAXY',c('WAXX','WAXR','WAXY',
                          'WBXX','WBXR','WBXY')] <- c((1-cbM)*(1-cM), (1-cbM)*cM*crM, (1-cbM),
                                                      cbM*(1-cM), cbM*cM*crM, cbM)/(2+cM*(crM-1))
  tMatrix['WWXX','ABXY',c('WAXX','WAXR','WAXY',
                          'WBXX','WBXR','WBXY')] <-c((1-cbM)*(1-cM), (1-cbM)*cM*crM, (1-cbM),
                                                     (1+cbM)*(1-cM), (1+cbM)*cM*crM, 1+cbM)/(4+2*cM*(crM-1))
  tMatrix['WWXX','BBXY',c('WBXX','WBXY')] <- c( 1, 1)/2
  tMatrix['WWXX','WWRY',c('WWXR','WWXY')] <- c( 1, 1)/2
  tMatrix['WWXX','WARY',c('WWXR','WWXY','WAXR','WAXY',
                          'WBXR','WBXY')] <- c( 1, 1, 1-cbM, 1-cbM, cbM, cbM)/4
  tMatrix['WWXX','WBRY',c('WWXR','WWXY','WBXR','WBXY')] <- c( 1, 1, 1, 1)/4
  tMatrix['WWXX','AARY',c('WAXR','WAXY','WBXR','WBXY')] <- c( 1-cbM, 1-cbM, cbM, cbM)/2
  tMatrix['WWXX','ABRY',c('WAXR','WAXY','WBXR','WBXY')] <- c( 1-cbM, 1-cbM, 1+cbM, 1+cbM)/4
  tMatrix['WWXX','BBRY',c('WBXR','WBXY')] <- c( 1, 1)/2


  # these stay the same for XXxXY
  normalizer <- (1+cF*(crF-1))*(2+cM*(crM-1))
  XRY <- c(((1-cF)*(1-cM)), (cF*crF*(1-cM) + (1-cF)*cM*crM), (1-cF), (cF*crF*cM*crM), (cF*crF))

  tMatrix['WAXX','WWXY',c('WWXX','WWXY','WWXR','WWRY',
                          'WAXX','WAXY','WAXR','WARY',
                          'WBXX','WBXY','WBXR','WBRY')] <- c(1-cF, 1-cF, cF*crF, cF*crF,
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*cF*crF, (1-cbF)*cF*crF,
                                                             cbF*(1-cF), cbF*(1-cF), cbF*cF*crF, cbF*cF*crF)/(4 + 4*cF*(crF-1))
  WAB <- c(1, (1-cbF) + (1-cbM), cbF + cbM, (1-cbF)*(1-cbM), cbF*(1-cbM) + (1-cbF)*cbM, cbF*cbM)
  tMatrix['WAXX','WAXY',c('WWXX','WWXR','WWXY','WWRR','WWRY',
                          'WAXX','WAXR','WAXY','WARR','WARY',
                          'WBXX','WBXR','WBXY','WBRR','WBRY',
                          'AAXX','AAXR','AAXY','AARR','AARY',
                          'ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- t(outer(X = WAB, Y = XRY, FUN = "*")/(4*normalizer))
  tMatrix['WAXX','WBXY',c('WWXX','WWXY','WWXR','WWRY',
                          'WAXX','WAXY','WAXR','WARY',
                          'WBXX','WBXY','WBXR','WBRY',
                          'ABXX','ABXY','ABXR','ABRY',
                          'BBXX','BBXY','BBXR','BBRY')] <- c(1-cF, 1-cF, cF*crF, cF*crF,
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*cF*crF, (1-cbF)*cF*crF,
                                                             cbF*(1-cF) + 1-cF, cbF*(1-cF) + 1-cF, cbF*cF*crF + cF*crF, cbF*cF*crF + cF*crF,
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*cF*crF, (1-cbF)*cF*crF,
                                                             cbF*(1-cF), cbF*(1-cF), cbF*cF*crF, cbF*cF*crF)/(8 + 8*cF*(crF-1))
  WAB <- c((1-cbM), cbM, (1-cbF)*(1-cbM), cbF*(1-cbM) + (1-cbF)*cbM, cbF*cbM)
  tMatrix['WAXX','AAXY',c('WAXX','WAXR','WAXY','WARR','WARY',
                          'WBXX','WBXR','WBXY','WBRR','WBRY',
                          'AAXX','AAXR','AAXY','AARR','AARY',
                          'ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- t(outer(X = WAB, Y = XRY, FUN = "*")/(2*normalizer))
  WAB <- c((1-cbM), 1+cbM, (1-cbF)*(1-cbM), cbF*(1-cbM) + (1-cbF)*(1+cbM), cbF*(1+cbM))
  tMatrix['WAXX','ABXY',c('WAXX','WAXR','WAXY','WARR','WARY',
                          'WBXX','WBXR','WBXY','WBRR','WBRY',
                          'AAXX','AAXR','AAXY','AARR','AARY',
                          'ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- t(outer(X = WAB, Y = XRY, FUN = "*")/(4*normalizer))
  tMatrix['WAXX','BBXY',c('WBXX','WBXY','WBXR','WBRY',
                          'ABXX','ABXY','ABXR','ABRY',
                          'BBXX','BBXY','BBXR','BBRY')] <- c(1-cF, 1-cF, cF*crF, cF*crF,
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*cF*crF, (1-cbF)*cF*crF,
                                                             cbF*(1-cF), cbF*(1-cF), cbF*cF*crF, cbF*cF*crF)/(4 + 4*cF*(crF-1))
  tMatrix['WAXX','WWRY',c('WWXR','WWXY','WWRR','WWRY',
                          'WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY')] <- c(1-cF, 1-cF, cF*crF, cF*crF,
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*cF*crF, (1-cbF)*cF*crF,
                                                             cbF*(1-cF), cbF*(1-cF), cbF*cF*crF, cbF*cF*crF)/(4 + 4*cF*(crF-1))
  tMatrix['WAXX','WARY',c('WWXR','WWXY','WWRR','WWRY',
                          'WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY',
                          'AAXR','AAXY','AARR','AARY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c(1-cF, 1-cF, cF*crF, cF*crF,
                                                             ((1-cbF) + (1-cbM))*(1-cF), ((1-cbF) + (1-cbM))*(1-cF), ((1-cbF) + (1-cbM))*cF*crF, ((1-cbF) + (1-cbM))*cF*crF,
                                                             (cbF + cbM)*(1-cF), (cbF + cbM)*(1-cF), (cbF + cbM)*cF*crF, (cbF + cbM)*cF*crF,
                                                             ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*cF*crF, ((1-cbF)*(1-cbM))*cF*crF,
                                                             (cbF*(1-cbM) + (1-cbF)*cbM)*(1-cF), (cbF*(1-cbM) + (1-cbF)*cbM)*(1-cF), (cbF*(1-cbM) + (1-cbF)*cbM)*cF*crF, (cbF*(1-cbM) + (1-cbF)*cbM)*cF*crF,
                                                             (cbF*cbM)*(1-cF), (cbF*cbM)*(1-cF), (cbF*cbM)*cF*crF, (cbF*cbM)*cF*crF)/(8 + 8*cF*(crF-1))
  tMatrix['WAXX','WBRY',c('WWXR','WWXY','WWRR','WWRY',
                          'WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c(1-cF, 1-cF, cF*crF, cF*crF,
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*cF*crF, (1-cbF)*cF*crF,
                                                             cbF*(1-cF) + 1-cF, cbF*(1-cF) + 1-cF, cbF*cF*crF + cF*crF, cbF*cF*crF + cF*crF,
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*cF*crF, (1-cbF)*cF*crF,
                                                             cbF*(1-cF), cbF*(1-cF), cbF*cF*crF, cbF*cF*crF)/(8 + 8*cF*(crF-1))
  tMatrix['WAXX','AARY',c('WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY',
                          'AAXR','AAXY','AARR','AARY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c(((1-cbM))*(1-cF), ((1-cbM))*(1-cF), ((1-cbM))*cF*crF, ((1-cbM))*cF*crF,
                                                             (cbM)*(1-cF), (cbM)*(1-cF), (cbM)*cF*crF, (cbM)*cF*crF,
                                                             ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*cF*crF, ((1-cbF)*(1-cbM))*cF*crF,
                                                             (cbF*(1-cbM) + (1-cbF)*cbM)*(1-cF), (cbF*(1-cbM) + (1-cbF)*cbM)*(1-cF), (cbF*(1-cbM) + (1-cbF)*cbM)*cF*crF, (cbF*(1-cbM) + (1-cbF)*cbM)*cF*crF,
                                                             (cbF*cbM)*(1-cF), (cbF*cbM)*(1-cF), (cbF*cbM)*cF*crF, (cbF*cbM)*cF*crF)/(4 + 4*cF*(crF-1))
  tMatrix['WAXX','ABRY',c('WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY',
                          'AAXR','AAXY','AARR','AARY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c(((1-cbM))*(1-cF), ((1-cbM))*(1-cF), ((1-cbM))*cF*crF, ((1-cbM))*cF*crF,
                                                             (1+cbM)*(1-cF), (1+cbM)*(1-cF), (1+cbM)*cF*crF, (1+cbM)*cF*crF,
                                                             ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*cF*crF, ((1-cbF)*(1-cbM))*cF*crF,
                                                             (cbF*(1-cbM) + (1-cbF)*(1+cbM))*(1-cF), (cbF*(1-cbM) + (1-cbF)*(1+cbM))*(1-cF), (cbF*(1-cbM) + (1-cbF)*(1+cbM))*cF*crF, (cbF*(1-cbM) + (1-cbF)*(1+cbM))*cF*crF,
                                                             (cbF*(1+cbM))*(1-cF), (cbF*(1+cbM))*(1-cF), (cbF*(1+cbM))*cF*crF, (cbF*(1+cbM))*cF*crF)/(8 + 8*cF*(crF-1))
  tMatrix['WAXX','BBRY',c('WBXR','WBXY','WBRR','WBRY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c(1-cF, 1-cF, cF*crF, cF*crF,
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*cF*crF, (1-cbF)*cF*crF,
                                                             cbF*(1-cF), cbF*(1-cF), cbF*cF*crF, cbF*cF*crF)/(4 + 4*cF*(crF-1))


  tMatrix['WBXX','WWXY',c('WWXX','WWXY','WBXX','WBXY')] <- c( 1, 1, 1, 1)/4
  tMatrix['WBXX','WAXY',c('WWXX','WWXR','WWXY',
                          'WAXX','WAXR','WAXY',
                          'WBXX','WBXR','WBXY',
                          'ABXX','ABXR','ABXY',
                          'BBXX','BBXR','BBXY')] <- c(1-cM, cM*crM, 1,
                                                      (1-cbM)*(1-cM), (1-cbM)*cM*crM, (1-cbM),
                                                      (1+cbM)*(1-cM), (1+cbM)*cM*crM, (1+cbM),
                                                      (1-cbM)*(1-cM), (1-cbM)*cM*crM, (1-cbM),
                                                      cbM*(1-cM), cbM*cM*crM, cbM)/(8+4*cM*(crM-1))
  tMatrix['WBXX','WBXY',c('WWXX','WWXY','WBXX','WBXY',
                          'BBXX','BBXY')] <- c( 1, 1, 2, 2, 1, 1)/8
  tMatrix['WBXX','AAXY',c('WAXX','WAXR','WAXY',
                          'WBXX','WBXR','WBXY',
                          'ABXX','ABXR','ABXY',
                          'BBXX','BBXR','BBXY')] <- c((1-cbM)*(1-cM), (1-cbM)*cM*crM, (1-cbM),
                                                      cbM*(1-cM), cbM*cM*crM, cbM,
                                                      (1-cbM)*(1-cM), (1-cbM)*cM*crM, (1-cbM),
                                                      cbM*(1-cM), cbM*cM*crM, cbM)/(4+2*cM*(crM-1))
  tMatrix['WBXX','ABXY',c('WAXX','WAXR','WAXY',
                          'WBXX','WBXR','WBXY',
                          'ABXX','ABXR','ABXY',
                          'BBXX','BBXR','BBXY')] <- c((1-cbM)*(1-cM), (1-cbM)*cM*crM, (1-cbM),
                                                      (1+cbM)*(1-cM), (1+cbM)*cM*crM, (1+cbM),
                                                      (1-cbM)*(1-cM), (1-cbM)*cM*crM, (1-cbM),
                                                      (1+cbM)*(1-cM), (1+cbM)*cM*crM, (1+cbM))/(8+4*cM*(crM-1))
  tMatrix['WBXX','BBXY',c('WBXX','WBXY','BBXX','BBXY')] <- c( 1, 1, 1, 1)/4
  tMatrix['WBXX','WWRY',c('WWXR','WWXY','WBXR','WBXY')] <- c( 1, 1, 1, 1)/4
  tMatrix['WBXX','WARY',c('WWXR','WWXY','WAXR','WAXY',
                          'WBXR','WBXY','ABXR','ABXY',
                          'BBXR','BBXY')] <- c( 1, 1, 1-cbM, 1-cbM,
                                                1+cbM, 1+cbM, 1-cbM, 1-cbM, cbM, cbM)/8
  tMatrix['WBXX','WBRY',c('WWXR','WWXY','WBXR','WBXY',
                          'BBXR','BBXY')] <- c( 1, 1, 2, 2, 1, 1)/8
  tMatrix['WBXX','AARY',c('WAXR','WAXY','WBXR','WBXY',
                          'ABXR','ABXY','BBXR','BBXY')] <- c( 1-cbM, 1-cbM, cbM, cbM,
                                                              1-cbM, 1-cbM, cbM, cbM)/4
  tMatrix['WBXX','ABRY',c('WAXR','WAXY','WBXR','WBXY',
                          'ABXR','ABXY','BBXR','BBXY')] <- c( 1-cbM, 1-cbM, 1+cbM, 1+cbM,
                                                              1-cbM, 1-cbM, 1+cbM, 1+cbM)/8
  tMatrix['WBXX','BBRY',c('WBXR','WBXY','BBXR','BBXY')] <- c( 1, 1, 1, 1)/4


  tMatrix['AAXX','WWXY',c('WAXX','WAXY','WAXR','WARY',
                          'WBXX','WBXY','WBXR','WBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(cF*crF), (1-cbF)*(cF*crF),
                                                             cbF*(1-cF), cbF*(1-cF), cbF*(cF*crF), cbF*(cF*crF))/(2 + 2*cF*(crF-1))
  WAB <- c(1-cbF, cbF, (1-cbF)*(1-cbM), cbF*(1-cbM) + (1-cbF)*cbM, cbF*cbM)
  tMatrix['AAXX','WAXY',c('WAXX','WAXR','WAXY','WARR','WARY',
                          'WBXX','WBXR','WBXY','WBRR','WBRY',
                          'AAXX','AAXR','AAXY','AARR','AARY',
                          'ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- t(outer(X = WAB, Y = XRY, FUN = "*")/(2*normalizer))
  tMatrix['AAXX','WBXY',c('WAXX','WAXY','WAXR','WARY',
                          'WBXX','WBXY','WBXR','WBRY',
                          'ABXX','ABXY','ABXR','ABRY',
                          'BBXX','BBXY','BBXR','BBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(cF*crF), (1-cbF)*(cF*crF),
                                                             cbF*(1-cF), cbF*(1-cF), cbF*(cF*crF), cbF*(cF*crF),
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(cF*crF), (1-cbF)*(cF*crF),
                                                             cbF*(1-cF), cbF*(1-cF), cbF*(cF*crF), cbF*(cF*crF))/(4 + 4*cF*(crF-1))
  WAB <- c((1-cbF)*(1-cbM), cbF*(1-cbM) + (1-cbF)*cbM, cbF*cbM)
  tMatrix['AAXX','AAXY',c('AAXX','AAXR','AAXY','AARR','AARY',
                          'ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- t(outer(X = WAB, Y = XRY, FUN = "*")/(normalizer))
  WAB <- c((1-cbF)*(1-cbM), cbF*(1-cbM) + (1-cbF)*(1+cbM), cbF*(1+cbM))
  tMatrix['AAXX','ABXY',c('AAXX','AAXR','AAXY','AARR','AARY',
                          'ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- t(outer(X = WAB, Y = XRY, FUN = "*")/(2*normalizer))
  tMatrix['AAXX','BBXY',c('ABXX','ABXY','ABXR','ABRY',
                          'BBXX','BBXY','BBXR','BBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(cF*crF), (1-cbF)*(cF*crF),
                                                             cbF*(1-cF), cbF*(1-cF), cbF*(cF*crF), cbF*(cF*crF))/(2 + 2*cF*(crF-1))
  tMatrix['AAXX','WWRY',c('WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(cF*crF), (1-cbF)*(cF*crF),
                                                             cbF*(1-cF), cbF*(1-cF), cbF*(cF*crF), cbF*(cF*crF))/(2 + 2*cF*(crF-1))
  tMatrix['AAXX','WARY',c('WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY',
                          'AAXR','AAXY','AARR','AARY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(cF*crF), (1-cbF)*(cF*crF),
                                                             (cbF)*(1-cF), (cbF)*(1-cF), (cbF)*(cF*crF), (cbF)*(cF*crF),
                                                             ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(cF*crF), ((1-cbF)*(1-cbM))*(cF*crF),
                                                             (cbF*(1-cbM) + (1-cbF)*cbM)*(1-cF), (cbF*(1-cbM) + (1-cbF)*cbM)*(1-cF), (cbF*(1-cbM) + (1-cbF)*cbM)*(cF*crF), (cbF*(1-cbM) + (1-cbF)*cbM)*(cF*crF),
                                                             (cbF*cbM)*(1-cF), (cbF*cbM)*(1-cF), (cbF*cbM)*(cF*crF), (cbF*cbM)*(cF*crF))/(4 + 4*cF*(crF-1))
  tMatrix['AAXX','WBRY',c('WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(cF*crF), (1-cbF)*(cF*crF),
                                                             cbF*(1-cF), cbF*(1-cF), cbF*(cF*crF), cbF*(cF*crF),
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(cF*crF), (1-cbF)*(cF*crF),
                                                             cbF*(1-cF), cbF*(1-cF), cbF*(cF*crF), cbF*(cF*crF))/(4 + 4*cF*(crF-1))
  tMatrix['AAXX','AARY',c('AAXR','AAXY','AARR','AARY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c(((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(cF*crF), ((1-cbF)*(1-cbM))*(cF*crF),
                                                             (cbF*(1-cbM) + (1-cbF)*cbM)*(1-cF), (cbF*(1-cbM) + (1-cbF)*cbM)*(1-cF), (cbF*(1-cbM) + (1-cbF)*cbM)*(cF*crF), (cbF*(1-cbM) + (1-cbF)*cbM)*(cF*crF),
                                                             (cbF*cbM)*(1-cF), (cbF*cbM)*(1-cF), (cbF*cbM)*(cF*crF), (cbF*cbM)*(cF*crF))/(2 + 2*cF*(crF-1))
  tMatrix['AAXX','ABRY',c('AAXR','AAXY','AARR','AARY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c(((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(cF*crF), ((1-cbF)*(1-cbM))*(cF*crF),
                                                             (cbF*(1-cbM) + (1-cbF)*(1+cbM))*(1-cF), (cbF*(1-cbM) + (1-cbF)*(1+cbM))*(1-cF), (cbF*(1-cbM) + (1-cbF)*(1+cbM))*(cF*crF), (cbF*(1-cbM) + (1-cbF)*(1+cbM))*(cF*crF),
                                                             (cbF*(1+cbM))*(1-cF), (cbF*(1+cbM))*(1-cF), (cbF*(1+cbM))*(cF*crF), (cbF*(1+cbM))*(cF*crF))/(4 + 4*cF*(crF-1))
  tMatrix['AAXX','BBRY',c('ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(cF*crF), (1-cbF)*(cF*crF),
                                                             cbF*(1-cF), cbF*(1-cF), cbF*(cF*crF), cbF*(cF*crF))/(2 + 2*cF*(crF-1))


  tMatrix['ABXX','WWXY',c('WAXX','WAXY','WAXR','WARY',
                          'WBXX','WBXY','WBXR','WBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(cF*crF), (1-cbF)*(cF*crF),
                                                             (1+cbF)*(1-cF), (1+cbF)*(1-cF), (1+cbF)*(cF*crF), (1+cbF)*(cF*crF))/(4 + 4*cF*(crF-1))
  WAB <- c(1-cbF, 1+cbF, (1-cbF)*(1-cbM), (1+cbF)*(1-cbM) + (1-cbF)*cbM, (1+cbF)*cbM)
  tMatrix['ABXX','WAXY',c('WAXX','WAXR','WAXY','WARR','WARY',
                          'WBXX','WBXR','WBXY','WBRR','WBRY',
                          'AAXX','AAXR','AAXY','AARR','AARY',
                          'ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- t(outer(X = WAB, Y = XRY, FUN = "*")/(4*normalizer))
  tMatrix['ABXX','WBXY',c('WAXX','WAXY','WAXR','WARY',
                          'WBXX','WBXY','WBXR','WBRY',
                          'ABXX','ABXY','ABXR','ABRY',
                          'BBXX','BBXY','BBXR','BBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(cF*crF), (1-cbF)*(cF*crF),
                                                             (1+cbF)*(1-cF), (1+cbF)*(1-cF), (1+cbF)*(cF*crF), (1+cbF)*(cF*crF),
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(cF*crF), (1-cbF)*(cF*crF),
                                                             (1+cbF)*(1-cF), (1+cbF)*(1-cF), (1+cbF)*(cF*crF), (1+cbF)*(cF*crF))/(8 + 8*cF*(crF-1))
  WAB <- c((1-cbF)*(1-cbM), (1+cbF)*(1-cbM) + (1-cbF)*cbM, (1+cbF)*cbM)
  tMatrix['ABXX','AAXY',c('AAXX','AAXR','AAXY','AARR','AARY',
                          'ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- t(outer(X = WAB, Y = XRY, FUN = "*")/(2*normalizer))
  WAB <- c((1-cbF)*(1-cbM), (1+cbF)*(1-cbM) + (1-cbF)*(1+cbM), (1+cbF)*(1+cbM))
  tMatrix['ABXX','ABXY',c('AAXX','AAXR','AAXY','AARR','AARY',
                          'ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- t(outer(X = WAB, Y = XRY, FUN = "*")/(4*normalizer))
  tMatrix['ABXX','BBXY',c('ABXX','ABXY','ABXR','ABRY',
                          'BBXX','BBXY','BBXR','BBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(cF*crF), (1-cbF)*(cF*crF),
                                                             (1+cbF)*(1-cF), (1+cbF)*(1-cF), (1+cbF)*(cF*crF), (1+cbF)*(cF*crF))/(4 + 4*cF*(crF-1))
  tMatrix['ABXX','WWRY',c('WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(cF*crF), (1-cbF)*(cF*crF),
                                                             (1+cbF)*(1-cF), (1+cbF)*(1-cF), (1+cbF)*(cF*crF), (1+cbF)*(cF*crF))/(4 + 4*cF*(crF-1))
  tMatrix['ABXX','WARY',c('WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY',
                          'AAXR','AAXY','AARR','AARY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(cF*crF), (1-cbF)*(cF*crF),
                                                             (1+cbF)*(1-cF), (1+cbF)*(1-cF), (1+cbF)*(cF*crF), (1+cbF)*(cF*crF),
                                                             (1-cbF)*(1-cbM)*(1-cF), (1-cbF)*(1-cbM)*(1-cF), (1-cbF)*(1-cbM)*(cF*crF), (1-cbF)*(1-cbM)*(cF*crF),
                                                             ((1+cbF)*(1-cbM) + (1-cbF)*cbM)*(1-cF), ((1+cbF)*(1-cbM) + (1-cbF)*cbM)*(1-cF), ((1+cbF)*(1-cbM) + (1-cbF)*cbM)*(cF*crF), ((1+cbF)*(1-cbM) + (1-cbF)*cbM)*(cF*crF),
                                                             ((1+cbF)*cbM)*(1-cF), ((1+cbF)*cbM)*(1-cF), ((1+cbF)*cbM)*(cF*crF), ((1+cbF)*cbM)*(cF*crF))/(8 + 8*cF*(crF-1))
  tMatrix['ABXX','WBRY',c('WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(cF*crF), (1-cbF)*(cF*crF),
                                                             (1+cbF)*(1-cF), (1+cbF)*(1-cF), (1+cbF)*(cF*crF), (1+cbF)*(cF*crF),
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(cF*crF), (1-cbF)*(cF*crF),
                                                             (1+cbF)*(1-cF), (1+cbF)*(1-cF), (1+cbF)*(cF*crF), (1+cbF)*(cF*crF))/(8 + 8*cF*(crF-1))
  tMatrix['ABXX','AARY',c('AAXR','AAXY','AARR','AARY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c((1-cbF)*(1-cbM)*(1-cF), (1-cbF)*(1-cbM)*(1-cF), (1-cbF)*(1-cbM)*(cF*crF), (1-cbF)*(1-cbM)*(cF*crF),
                                                             ((1+cbF)*(1-cbM) + (1-cbF)*cbM)*(1-cF), ((1+cbF)*(1-cbM) + (1-cbF)*cbM)*(1-cF), ((1+cbF)*(1-cbM) + (1-cbF)*cbM)*(cF*crF), ((1+cbF)*(1-cbM) + (1-cbF)*cbM)*(cF*crF),
                                                             ((1+cbF)*cbM)*(1-cF), ((1+cbF)*cbM)*(1-cF), ((1+cbF)*cbM)*(cF*crF), ((1+cbF)*cbM)*(cF*crF))/(4 + 4*cF*(crF-1))
  tMatrix['ABXX','ABRY',c('AAXR','AAXY','AARR','AARY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c((1-cbF)*(1-cbM)*(1-cF), (1-cbF)*(1-cbM)*(1-cF), (1-cbF)*(1-cbM)*(cF*crF), (1-cbF)*(1-cbM)*(cF*crF),
                                                             ((1+cbF)*(1-cbM) + (1-cbF)*(1+cbM))*(1-cF), ((1+cbF)*(1-cbM) + (1-cbF)*(1+cbM))*(1-cF), ((1+cbF)*(1-cbM) + (1-cbF)*(1+cbM))*(cF*crF), ((1+cbF)*(1-cbM) + (1-cbF)*(1+cbM))*(cF*crF),
                                                             ((1+cbF)*(1+cbM))*(1-cF), ((1+cbF)*(1+cbM))*(1-cF), ((1+cbF)*(1+cbM))*(cF*crF), ((1+cbF)*(1+cbM))*(cF*crF))/(8 + 8*cF*(crF-1))
  tMatrix['ABXX','BBRY',c('ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(cF*crF), (1-cbF)*(cF*crF),
                                                             (1+cbF)*(1-cF), (1+cbF)*(1-cF), (1+cbF)*(cF*crF), (1+cbF)*(cF*crF))/(4 + 4*cF*(crF-1))


  tMatrix['BBXX','WWXY',c('WBXX','WBXY')] <- c( 1, 1)/2
  tMatrix['BBXX','WAXY',c('WBXX','WBXR','WBXY',
                          'ABXX','ABXR','ABXY',
                          'BBXX','BBXR','BBXY')] <- c(1-cM, cM*crM, 1,
                                                      (1-cbM)*(1-cM), (1-cbM)*cM*crM, (1-cbM),
                                                      cbM*(1-cM), cbM*cM*crM, cbM)/(4+2*cM*(crM-1))
  tMatrix['BBXX','WBXY',c('WBXX','WBXY','BBXX','BBXY')] <- c( 1, 1, 1, 1)/4
  tMatrix['BBXX','AAXY',c('ABXX','ABXR','ABXY',
                          'BBXX','BBXR','BBXY')] <-c((1-cbM)*(1-cM), (1-cbM)*cM*crM, (1-cbM),
                                                     cbM*(1-cM), cbM*cM*crM, cbM)/(2+cM*(crM-1))
  tMatrix['BBXX','ABXY',c('ABXX','ABXR','ABXY',
                          'BBXX','BBXR','BBXY')] <-c((1-cbM)*(1-cM), (1-cbM)*cM*crM, (1-cbM),
                                                     (1+cbM)*(1-cM), (1+cbM)*cM*crM, 1+cbM)/(4+2*cM*(crM-1))
  tMatrix['BBXX','BBXY',c('BBXX','BBXY')] <- c( 1, 1)/2
  tMatrix['BBXX','WWRY',c('WBXR','WBXY')] <- c( 1, 1)/2
  tMatrix['BBXX','WARY',c('WBXR','WBXY','ABXR','ABXY',
                          'BBXR','BBXY')] <- c( 1, 1, 1-cbM, 1-cbM, cbM, cbM)/4
  tMatrix['BBXX','WBRY',c('WBXR','WBXY','BBXR','BBXY')] <- c( 1, 1, 1, 1)/4
  tMatrix['BBXX','AARY',c('ABXR','ABXY','BBXR','BBXY')] <- c( 1-cbM, 1-cbM, cbM, cbM)/2
  tMatrix['BBXX','ABRY',c('ABXR','ABXY','BBXR','BBXY')] <- c( 1-cbM, 1-cbM, 1+cbM, 1+cbM)/4
  tMatrix['BBXX','BBRY',c('BBXR','BBXY')] <- c( 1, 1)/2


  tMatrix['WWXR','WWXY',c('WWXX','WWXY','WWXR','WWRY')] <- c( 1, 1, 1, 1)/4
  tMatrix['WWXR','WAXY',c('WWXX','WWXR','WWXY','WWRR','WWRY',
                          'WAXX','WAXR','WAXY','WARR','WARY',
                          'WBXX','WBXR','WBXY','WBRR','WBRY')] <- c(1-cM, 1+cM*(crM-1), 1, cM*crM, 1,
                                                                    (1-cbM)*(1-cM), (1-cbM)*(1+cM*(crM-1)), (1-cbM), (1-cbM)*cM*crM, (1-cbM),
                                                                    cbM*(1-cM), cbM*(1+cM*(crM-1)), cbM, cbM*cM*crM, cbM)/(8+4*cM*(crM-1))
  tMatrix['WWXR','WBXY',c('WWXX','WWXY','WWXR','WWRY',
                          'WBXX','WBXY','WBXR','WBRY')] <- c( 1, 1, 1, 1, 1, 1, 1, 1)/8
  tMatrix['WWXR','AAXY',c('WAXX','WAXR','WAXY','WARR','WARY',
                          'WBXX','WBXR','WBXY','WBRR','WBRY')] <- c((1-cbM)*(1-cM), (1-cbM)*(1+cM*(crM-1)), (1-cbM), (1-cbM)*cM*crM, (1-cbM),
                                                                   cbM*(1-cM), cbM*(1+cM*(crM-1)), cbM, cbM*cM*crM, cbM)/(4+2*cM*(crM-1))
  tMatrix['WWXR','ABXY',c('WAXX','WAXR','WAXY','WARR','WARY',
                          'WBXX','WBXR','WBXY','WBRR','WBRY')] <- c((1-cbM)*(1-cM), (1-cbM)*(1+cM*(crM-1)), (1-cbM), (1-cbM)*cM*crM, (1-cbM),
                                                                    (1+cbM)*(1-cM), (1+cbM)*(1+cM*(crM-1)), 1+cbM, (1+cbM)*cM*crM, 1+cbM)/(8+4*cM*(crM-1))
  tMatrix['WWXR','BBXY',c('WBXX','WBXY','WBXR','WBRY')] <- c( 1, 1, 1, 1)/4
  tMatrix['WWXR','WWRY',c('WWXR','WWXY','WWRR','WWRY')] <- c( 1, 1, 1, 1)/4
  tMatrix['WWXR','WARY',c('WWXR','WWXY','WWRR','WWRY',
                          'WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY')] <- c( 1, 1, 1, 1,
                                                              1-cbM, 1-cbM, 1-cbM, 1-cbM,
                                                              cbM, cbM, cbM, cbM)/8
  tMatrix['WWXR','WBRY',c('WWXR','WWXY','WWRR','WWRY',
                          'WBXR','WBXY','WBRR','WBRY')] <- c( 1, 1, 1, 1,
                                                              1, 1, 1, 1)/8
  tMatrix['WWXR','AARY',c('WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY')] <- c( 1-cbM, 1-cbM, 1-cbM, 1-cbM,
                                                              cbM, cbM, cbM, cbM)/4
  tMatrix['WWXR','ABRY',c('WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY')] <- c( 1-cbM, 1-cbM, 1-cbM, 1-cbM,
                                                              1+cbM, 1+cbM, 1+cbM, 1+cbM)/8
  tMatrix['WWXR','BBRY',c('WBXR','WBXY','WBRR','WBRY')] <- c( 1, 1, 1, 1)/4


  # change normalizer and XRY vectors for rest of cube
  normalizer <- (2+cF*(crF-1))*(2+cM*(crM-1))
  XRY <- c((1-cF)*(1-cM), (1+cF*crF)*(1-cM) + (1-cF)*(cM*crM), 1-cF, (1+cF*crF)*(cM*crM), 1+cF*crF)

  tMatrix['WAXR','WWXY',c('WWXX','WWXY','WWXR','WWRY',
                          'WAXX','WAXY','WAXR','WARY',
                          'WBXX','WBXY','WBXR','WBRY')] <- c(1-cF, 1-cF, 1+cF*crF, 1+cF*crF,
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             cbF*(1-cF), cbF*(1-cF), cbF*(1+cF*crF), cbF*(1+cF*crF))/(8 + 4*cF*(crF-1))
  WAB <- c(1, (1-cbF) + (1-cbM), cbF + cbM, (1-cbF)*(1-cbM), cbF*(1-cbM) + (1-cbF)*cbM, cbF*cbM)
  tMatrix['WAXR','WAXY',c('WWXX','WWXR','WWXY','WWRR','WWRY',
                          'WAXX','WAXR','WAXY','WARR','WARY',
                          'WBXX','WBXR','WBXY','WBRR','WBRY',
                          'AAXX','AAXR','AAXY','AARR','AARY',
                          'ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- t(outer(X = WAB, Y = XRY, FUN = "*")/(4*normalizer))
  tMatrix['WAXR','WBXY',c('WWXX','WWXY','WWXR','WWRY',
                          'WAXX','WAXY','WAXR','WARY',
                          'WBXX','WBXY','WBXR','WBRY',
                          'ABXX','ABXY','ABXR','ABRY',
                          'BBXX','BBXY','BBXR','BBRY')] <- c(1-cF, 1-cF, 1+cF*crF, 1+cF*crF,
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             (1+cbF)*(1-cF), (1+cbF)*(1-cF), (1+cbF)*(1+cF*crF), (1+cbF)*(1+cF*crF),
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             cbF*(1-cF), cbF*(1-cF), cbF*(1+cF*crF), cbF*(1+cF*crF))/(16 + 8*cF*(crF-1))
  WAB <- c((1-cbM), cbM, (1-cbF)*(1-cbM), cbF*(1-cbM) + (1-cbF)*cbM, cbF*cbM)
  tMatrix['WAXR','AAXY',c('WAXX','WAXR','WAXY','WARR','WARY',
                          'WBXX','WBXR','WBXY','WBRR','WBRY',
                          'AAXX','AAXR','AAXY','AARR','AARY',
                          'ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- t(outer(X = WAB, Y = XRY, FUN = "*")/(2*normalizer))
  WAB <- c((1-cbM), 1+cbM, (1-cbF)*(1-cbM), cbF*(1-cbM) + (1-cbF)*(1+cbM), cbF*(1+cbM))
  tMatrix['WAXR','ABXY',c('WAXX','WAXR','WAXY','WARR','WARY',
                          'WBXX','WBXR','WBXY','WBRR','WBRY',
                          'AAXX','AAXR','AAXY','AARR','AARY',
                          'ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- t(outer(X = WAB, Y = XRY, FUN = "*")/(4*normalizer))
  tMatrix['WAXR','BBXY',c('WBXX','WBXY','WBXR','WBRY',
                          'ABXX','ABXY','ABXR','ABRY',
                          'BBXX','BBXY','BBXR','BBRY')] <- c(1-cF, 1-cF, 1+cF*crF, 1+cF*crF,
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             cbF*(1-cF), cbF*(1-cF), cbF*(1+cF*crF), cbF*(1+cF*crF))/(8 + 4*cF*(crF-1))
  tMatrix['WAXR','WWRY',c('WWXR','WWXY','WWRR','WWRY',
                          'WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY')] <- c(1-cF, 1-cF, 1+cF*crF, 1+cF*crF,
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             cbF*(1-cF), cbF*(1-cF), cbF*(1+cF*crF), cbF*(1+cF*crF))/(8 + 4*cF*(crF-1))
  tMatrix['WAXR','WARY',c('WWXR','WWXY','WWRR','WWRY',
                          'WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY',
                          'AAXR','AAXY','AARR','AARY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c(1-cF, 1-cF, 1+cF*crF, 1+cF*crF,
                                                             (1-cbF + 1-cbM)*(1-cF), (1-cbF + 1-cbM)*(1-cF), (1-cbF + 1-cbM)*(1+cF*crF), (1-cbF + 1-cbM)*(1+cF*crF),
                                                             (cbF + cbM)*(1-cF), (cbF + cbM)*(1-cF), (cbF + cbM)*(1+cF*crF), (cbF + cbM)*(1+cF*crF),
                                                             ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1+cF*crF), ((1-cbF)*(1-cbM))*(1+cF*crF),
                                                             (cbF*(1-cbM) + (1-cbF)*cbM)*(1-cF), (cbF*(1-cbM) + (1-cbF)*cbM)*(1-cF), (cbF*(1-cbM) + (1-cbF)*cbM)*(1+cF*crF), (cbF*(1-cbM) + (1-cbF)*cbM)*(1+cF*crF),
                                                             (cbF*cbM)*(1-cF), (cbF*cbM)*(1-cF), (cbF*cbM)*(1+cF*crF), (cbF*cbM)*(1+cF*crF))/(16 + 8*cF*(crF-1))
  tMatrix['WAXR','WBRY',c('WWXR','WWXY','WWRR','WWRY',
                          'WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c(1-cF, 1-cF, 1+cF*crF, 1+cF*crF,
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             (1+cbF)*(1-cF), (1+cbF)*(1-cF), (1+cbF)*(1+cF*crF), (1+cbF)*(1+cF*crF),
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             cbF*(1-cF), cbF*(1-cF), cbF*(1+cF*crF), cbF*(1+cF*crF))/(16 + 8*cF*(crF-1))
  tMatrix['WAXR','AARY',c('WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY',
                          'AAXR','AAXY','AARR','AARY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c((1-cbM)*(1-cF), (1-cbM)*(1-cF), (1-cbM)*(1+cF*crF), (1-cbM)*(1+cF*crF),
                                                             (cbM)*(1-cF), (cbM)*(1-cF), (cbM)*(1+cF*crF), (cbM)*(1+cF*crF),
                                                             ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1+cF*crF), ((1-cbF)*(1-cbM))*(1+cF*crF),
                                                             (cbF*(1-cbM) + (1-cbF)*cbM)*(1-cF), (cbF*(1-cbM) + (1-cbF)*cbM)*(1-cF), (cbF*(1-cbM) + (1-cbF)*cbM)*(1+cF*crF), (cbF*(1-cbM) + (1-cbF)*cbM)*(1+cF*crF),
                                                             (cbF*cbM)*(1-cF), (cbF*cbM)*(1-cF), (cbF*cbM)*(1+cF*crF), (cbF*cbM)*(1+cF*crF))/(8 + 4*cF*(crF-1))
  tMatrix['WAXR','ABRY',c('WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY',
                          'AAXR','AAXY','AARR','AARY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c((1-cbM)*(1-cF), (1-cbM)*(1-cF), (1-cbM)*(1+cF*crF), (1-cbM)*(1+cF*crF),
                                                             (1+cbM)*(1-cF), (1+cbM)*(1-cF), (1+cbM)*(1+cF*crF), (1+cbM)*(1+cF*crF),
                                                             ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1+cF*crF), ((1-cbF)*(1-cbM))*(1+cF*crF),
                                                             (cbF*(1-cbM) + (1-cbF)*(1+cbM))*(1-cF), (cbF*(1-cbM) + (1-cbF)*(1+cbM))*(1-cF), (cbF*(1-cbM) + (1-cbF)*(1+cbM))*(1+cF*crF), (cbF*(1-cbM) + (1-cbF)*(1+cbM))*(1+cF*crF),
                                                             (cbF*(1+cbM))*(1-cF), (cbF*(1+cbM))*(1-cF), (cbF*(1+cbM))*(1+cF*crF), (cbF*(1+cbM))*(1+cF*crF))/(16 + 8*cF*(crF-1))
  tMatrix['WAXR','BBRY',c('WBXR','WBXY','WBRR','WBRY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c(1-cF, 1-cF, 1+cF*crF, 1+cF*crF,
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             cbF*(1-cF), cbF*(1-cF), cbF*(1+cF*crF), cbF*(1+cF*crF))/(8 + 4*cF*(crF-1))


  tMatrix['WBXR','WWXY',c('WWXX','WWXY','WWXR','WWRY',
                          'WBXX','WBXY','WBXR','WBRY')] <- c( 1, 1, 1, 1, 1, 1, 1, 1)/8
  tMatrix['WBXR','WAXY',c('WWXX','WWXR','WWXY','WWRR','WWRY',
                          'WAXX','WAXR','WAXY','WARR','WARY',
                          'WBXX','WBXR','WBXY','WBRR','WBRY',
                          'ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- c((1-cM), (1+cM*(crM-1)), 1, cM*crM, 1,
                                                                    (1-cbM)*(1-cM), (1-cbM)*(1+cM*(crM-1)), (1-cbM), (1-cbM)*cM*crM, (1-cbM),
                                                                    (1+cbM)*(1-cM), (1+cbM)*(1+cM*(crM-1)), (1+cbM), (1+cbM)*cM*crM, (1+cbM),
                                                                    (1-cbM)*(1-cM), (1-cbM)*(1+cM*(crM-1)), (1-cbM), (1-cbM)*cM*crM, (1-cbM),
                                                                    cbM*(1-cM), cbM*(1+cM*(crM-1)), cbM, cbM*cM*crM, cbM)/(16 + 8*cM*(crM-1))
  tMatrix['WBXR','WBXY',c('WWXX','WWXY','WWXR','WWRY',
                          'WBXX','WBXY','WBXR','WBRY',
                          'BBXX','BBXY','BBXR','BBRY')] <- c( 1, 1, 1, 1,
                                                              2, 2, 2, 2,
                                                              1, 1, 1, 1)/16
  tMatrix['WBXR','AAXY',c('WAXX','WAXR','WAXY','WARR','WARY',
                          'WBXX','WBXR','WBXY','WBRR','WBRY',
                          'ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- c((1-cbM)*(1-cM), (1-cbM)*(1+cM*(crM-1)), (1-cbM), (1-cbM)*cM*crM, (1-cbM),
                                                                    cbM*(1-cM), cbM*(1+cM*(crM-1)), cbM, cbM*cM*crM, cbM,
                                                                    (1-cbM)*(1-cM), (1-cbM)*(1+cM*(crM-1)), (1-cbM), (1-cbM)*cM*crM, (1-cbM),
                                                                    cbM*(1-cM), cbM*(1+cM*(crM-1)), cbM, cbM*cM*crM, cbM)/(8 + 4*cM*(crM-1))
  tMatrix['WBXR','ABXY',c('WAXX','WAXR','WAXY','WARR','WARY',
                          'WBXX','WBXR','WBXY','WBRR','WBRY',
                          'ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- c((1-cbM)*(1-cM), (1-cbM)*(1+cM*(crM-1)), (1-cbM), (1-cbM)*cM*crM, (1-cbM),
                                                                    (1+cbM)*(1-cM), (1+cbM)*(1+cM*(crM-1)), (1+cbM), (1+cbM)*cM*crM, (1+cbM),
                                                                    (1-cbM)*(1-cM), (1-cbM)*(1+cM*(crM-1)), (1-cbM), (1-cbM)*cM*crM, (1-cbM),
                                                                    (1+cbM)*(1-cM), (1+cbM)*(1+cM*(crM-1)), (1+cbM), (1+cbM)*cM*crM, (1+cbM))/(16 + 8*cM*(crM-1))
  tMatrix['WBXR','BBXY',c('WBXX','WBXY','WBXR','WBRY',
                          'BBXX','BBXY','BBXR','BBRY')] <- c( 1, 1, 1, 1,
                                                              1, 1, 1, 1)/8
  tMatrix['WBXR','WWRY',c('WWXR','WWXY','WWRR','WWRY',
                          'WBXR','WBXY','WBRR','WBRY')] <- c( 1, 1, 1, 1,
                                                              1, 1, 1, 1)/8
  tMatrix['WBXR','WARY',c('WWXR','WWXY','WWRR','WWRY',
                          'WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c( 1, 1, 1, 1,
                                                              1-cbM, 1-cbM, 1-cbM, 1-cbM,
                                                              1+cbM, 1+cbM, 1+cbM, 1+cbM,
                                                              1-cbM, 1-cbM, 1-cbM, 1-cbM,
                                                              cbM, cbM, cbM, cbM)/16
  tMatrix['WBXR','WBRY',c('WWXR','WWXY','WWRR','WWRY',
                          'WBXR','WBXY','WBRR','WBRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c( 1, 1, 1, 1,
                                                              2, 2, 2, 2,
                                                              1, 1, 1, 1)/16
  tMatrix['WBXR','AARY',c('WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c( 1-cbM, 1-cbM, 1-cbM, 1-cbM,
                                                              cbM, cbM, cbM, cbM,
                                                              1-cbM, 1-cbM, 1-cbM, 1-cbM,
                                                              cbM, cbM, cbM, cbM)/8
  tMatrix['WBXR','ABRY',c('WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c( 1-cbM, 1-cbM, 1-cbM, 1-cbM,
                                                              1+cbM, 1+cbM, 1+cbM, 1+cbM,
                                                              1-cbM, 1-cbM, 1-cbM, 1-cbM,
                                                              1+cbM, 1+cbM, 1+cbM, 1+cbM)/16
  tMatrix['WBXR','BBRY',c('WBXR','WBXY','WBRR','WBRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c( 1, 1, 1, 1,
                                                              1, 1, 1, 1)/8


  tMatrix['AAXR','WWXY',c('WAXX','WAXY','WAXR','WARY',
                          'WBXX','WBXY','WBXR','WBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             cbF*(1-cF), cbF*(1-cF), cbF*(1+cF*crF), cbF*(1+cF*crF))/(4 + 2*cF*(crF-1))
  WAB <- c(1-cbF, cbF, (1-cbF)*(1-cbM), cbF*(1-cbM) + (1-cbF)*cbM, cbF*cbM)
  tMatrix['AAXR','WAXY',c('WAXX','WAXR','WAXY','WARR','WARY',
                          'WBXX','WBXR','WBXY','WBRR','WBRY',
                          'AAXX','AAXR','AAXY','AARR','AARY',
                          'ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- t(outer(X = WAB, Y = XRY, FUN = "*")/(2*normalizer))
  tMatrix['AAXR','WBXY',c('WAXX','WAXY','WAXR','WARY',
                          'WBXX','WBXY','WBXR','WBRY',
                          'ABXX','ABXY','ABXR','ABRY',
                          'BBXX','BBXY','BBXR','BBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             cbF*(1-cF), cbF*(1-cF), cbF*(1+cF*crF), cbF*(1+cF*crF),
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             cbF*(1-cF), cbF*(1-cF), cbF*(1+cF*crF), cbF*(1+cF*crF))/(8 + 4*cF*(crF-1))
  WAB <- c((1-cbF)*(1-cbM), cbF*(1-cbM) + (1-cbF)*cbM, cbF*cbM)
  tMatrix['AAXR','AAXY',c('AAXX','AAXR','AAXY','AARR','AARY',
                          'ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- t(outer(X = WAB, Y = XRY, FUN = "*")/(normalizer))
  WAB <- c((1-cbF)*(1-cbM), cbF*(1-cbM) + (1-cbF)*(1+cbM), cbF*(1+cbM))
  tMatrix['AAXR','ABXY',c('AAXX','AAXR','AAXY','AARR','AARY',
                          'ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- t(outer(X = WAB, Y = XRY, FUN = "*")/(2*normalizer))
  tMatrix['AAXR','BBXY',c('ABXX','ABXY','ABXR','ABRY',
                          'BBXX','BBXY','BBXR','BBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             cbF*(1-cF), cbF*(1-cF), cbF*(1+cF*crF), cbF*(1+cF*crF))/(4 + 2*cF*(crF-1))
  tMatrix['AAXR','WWRY',c('WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             cbF*(1-cF), cbF*(1-cF), cbF*(1+cF*crF), cbF*(1+cF*crF))/(4 + 2*cF*(crF-1))
  tMatrix['AAXR','WARY',c('WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY',
                          'AAXR','AAXY','AARR','AARY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             cbF*(1-cF), cbF*(1-cF), cbF*(1+cF*crF), cbF*(1+cF*crF),
                                                             ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1+cF*crF), ((1-cbF)*(1-cbM))*(1+cF*crF),
                                                             (cbF*(1-cbM) + (1-cbF)*cbM)*(1-cF), (cbF*(1-cbM) + (1-cbF)*cbM)*(1-cF), (cbF*(1-cbM) + (1-cbF)*cbM)*(1+cF*crF), (cbF*(1-cbM) + (1-cbF)*cbM)*(1+cF*crF),
                                                             (cbF*cbM)*(1-cF), (cbF*cbM)*(1-cF), (cbF*cbM)*(1+cF*crF), (cbF*cbM)*(1+cF*crF))/(8 + 4*cF*(crF-1))
  tMatrix['AAXR','WBRY',c('WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             cbF*(1-cF), cbF*(1-cF), cbF*(1+cF*crF), cbF*(1+cF*crF),
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             cbF*(1-cF), cbF*(1-cF), cbF*(1+cF*crF), cbF*(1+cF*crF))/(8 + 4*cF*(crF-1))
  tMatrix['AAXR','AARY',c('AAXR','AAXY','AARR','AARY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c(((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1+cF*crF), ((1-cbF)*(1-cbM))*(1+cF*crF),
                                                             (cbF*(1-cbM) + (1-cbF)*cbM)*(1-cF), (cbF*(1-cbM) + (1-cbF)*cbM)*(1-cF), (cbF*(1-cbM) + (1-cbF)*cbM)*(1+cF*crF), (cbF*(1-cbM) + (1-cbF)*cbM)*(1+cF*crF),
                                                             (cbF*cbM)*(1-cF), (cbF*cbM)*(1-cF), (cbF*cbM)*(1+cF*crF), (cbF*cbM)*(1+cF*crF))/(4 + 2*cF*(crF-1))
  tMatrix['AAXR','ABRY',c('AAXR','AAXY','AARR','AARY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c(((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1+cF*crF), ((1-cbF)*(1-cbM))*(1+cF*crF),
                                                             (cbF*(1-cbM) + (1-cbF)*(1+cbM))*(1-cF), (cbF*(1-cbM) + (1-cbF)*(1+cbM))*(1-cF), (cbF*(1-cbM) + (1-cbF)*(1+cbM))*(1+cF*crF), (cbF*(1-cbM) + (1-cbF)*(1+cbM))*(1+cF*crF),
                                                             (cbF*(1+cbM))*(1-cF), (cbF*(1+cbM))*(1-cF), (cbF*(1+cbM))*(1+cF*crF), (cbF*(1+cbM))*(1+cF*crF))/(8 + 4*cF*(crF-1))
  tMatrix['AAXR','BBRY',c('ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             cbF*(1-cF), cbF*(1-cF), cbF*(1+cF*crF), cbF*(1+cF*crF))/(4 + 2*cF*(crF-1))


  tMatrix['ABXR','WWXY',c('WAXX','WAXY','WAXR','WARY',
                          'WBXX','WBXY','WBXR','WBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             (1+cbF)*(1-cF), (1+cbF)*(1-cF), (1+cbF)*(1+cF*crF), (1+cbF)*(1+cF*crF))/(8 + 4*cF*(crF-1))
  WAB <- c(1-cbF, 1+cbF, (1-cbF)*(1-cbM), (1+cbF)*(1-cbM) + (1-cbF)*cbM, (1+cbF)*cbM)
  tMatrix['ABXR','WAXY',c('WAXX','WAXR','WAXY','WARR','WARY',
                          'WBXX','WBXR','WBXY','WBRR','WBRY',
                          'AAXX','AAXR','AAXY','AARR','AARY',
                          'ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- t(outer(X = WAB, Y = XRY, FUN = "*")/(4*normalizer))
  tMatrix['ABXR','WBXY',c('WAXX','WAXY','WAXR','WARY',
                          'WBXX','WBXY','WBXR','WBRY',
                          'ABXX','ABXY','ABXR','ABRY',
                          'BBXX','BBXY','BBXR','BBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             (1+cbF)*(1-cF), (1+cbF)*(1-cF), (1+cbF)*(1+cF*crF), (1+cbF)*(1+cF*crF),
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             (1+cbF)*(1-cF), (1+cbF)*(1-cF), (1+cbF)*(1+cF*crF), (1+cbF)*(1+cF*crF))/(16 + 8*cF*(crF-1))
  WAB <- c((1-cbF)*(1-cbM), (1+cbF)*(1-cbM) + (1-cbF)*cbM, (1+cbF)*cbM)
  tMatrix['ABXR','AAXY',c('AAXX','AAXR','AAXY','AARR','AARY',
                          'ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- t(outer(X = WAB, Y = XRY, FUN = "*")/(2*normalizer))
  WAB <- c((1-cbF)*(1-cbM), (1+cbF)*(1-cbM) + (1-cbF)*(1+cbM), (1+cbF)*(1+cbM))
  tMatrix['ABXR','ABXY',c('AAXX','AAXR','AAXY','AARR','AARY',
                          'ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- t(outer(X = WAB, Y = XRY, FUN = "*")/(4*normalizer))
  tMatrix['ABXR','BBXY',c('ABXX','ABXY','ABXR','ABRY',
                          'BBXX','BBXY','BBXR','BBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             (1+cbF)*(1-cF), (1+cbF)*(1-cF), (1+cbF)*(1+cF*crF), (1+cbF)*(1+cF*crF))/(8 + 4*cF*(crF-1))
  tMatrix['ABXR','WWRY',c('WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             (1+cbF)*(1-cF), (1+cbF)*(1-cF), (1+cbF)*(1+cF*crF), (1+cbF)*(1+cF*crF))/(8 + 4*cF*(crF-1))
  tMatrix['ABXR','WARY',c('WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY',
                          'AAXR','AAXY','AARR','AARY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             (1+cbF)*(1-cF), (1+cbF)*(1-cF), (1+cbF)*(1+cF*crF), (1+cbF)*(1+cF*crF),
                                                             ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1+cF*crF), ((1-cbF)*(1-cbM))*(1+cF*crF),
                                                             ((1+cbF)*(1-cbM) + (1-cbF)*cbM)*(1-cF), ((1+cbF)*(1-cbM) + (1-cbF)*cbM)*(1-cF), ((1+cbF)*(1-cbM) + (1-cbF)*cbM)*(1+cF*crF), ((1+cbF)*(1-cbM) + (1-cbF)*cbM)*(1+cF*crF),
                                                             ((1+cbF)*cbM)*(1-cF), ((1+cbF)*cbM)*(1-cF), ((1+cbF)*cbM)*(1+cF*crF), ((1+cbF)*cbM)*(1+cF*crF))/(16 + 8*cF*(crF-1))
  tMatrix['ABXR','WBRY',c('WAXR','WAXY','WARR','WARY',
                          'WBXR','WBXY','WBRR','WBRY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             (1+cbF)*(1-cF), (1+cbF)*(1-cF), (1+cbF)*(1+cF*crF), (1+cbF)*(1+cF*crF),
                                                             (1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             (1+cbF)*(1-cF), (1+cbF)*(1-cF), (1+cbF)*(1+cF*crF), (1+cbF)*(1+cF*crF))/(16 + 8*cF*(crF-1))
  tMatrix['ABXR','AARY',c('AAXR','AAXY','AARR','AARY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c(((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1+cF*crF), ((1-cbF)*(1-cbM))*(1+cF*crF),
                                                             ((1+cbF)*(1-cbM) + (1-cbF)*cbM)*(1-cF), ((1+cbF)*(1-cbM) + (1-cbF)*cbM)*(1-cF), ((1+cbF)*(1-cbM) + (1-cbF)*cbM)*(1+cF*crF), ((1+cbF)*(1-cbM) + (1-cbF)*cbM)*(1+cF*crF),
                                                             ((1+cbF)*cbM)*(1-cF), ((1+cbF)*cbM)*(1-cF), ((1+cbF)*cbM)*(1+cF*crF), ((1+cbF)*cbM)*(1+cF*crF))/(8 + 4*cF*(crF-1))
  tMatrix['ABXR','ABRY',c('AAXR','AAXY','AARR','AARY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c(((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1-cF), ((1-cbF)*(1-cbM))*(1+cF*crF), ((1-cbF)*(1-cbM))*(1+cF*crF),
                                                             ((1+cbF)*(1-cbM) + (1-cbF)*(1+cbM))*(1-cF), ((1+cbF)*(1-cbM) + (1-cbF)*(1+cbM))*(1-cF), ((1+cbF)*(1-cbM) + (1-cbF)*(1+cbM))*(1+cF*crF), ((1+cbF)*(1-cbM) + (1-cbF)*(1+cbM))*(1+cF*crF),
                                                             ((1+cbF)*(1+cbM))*(1-cF), ((1+cbF)*(1+cbM))*(1-cF), ((1+cbF)*(1+cbM))*(1+cF*crF), ((1+cbF)*(1+cbM))*(1+cF*crF))/(16 + 8*cF*(crF-1))
  tMatrix['ABXR','BBRY',c('ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c((1-cbF)*(1-cF), (1-cbF)*(1-cF), (1-cbF)*(1+cF*crF), (1-cbF)*(1+cF*crF),
                                                             (1+cbF)*(1-cF), (1+cbF)*(1-cF), (1+cbF)*(1+cF*crF), (1+cbF)*(1+cF*crF))/(8 + 4*cF*(crF-1))


  tMatrix['BBXR','WWXY',c('WBXX','WBXY','WBXR','WBRY')] <- c( 1, 1, 1, 1)/4
  tMatrix['BBXR','WAXY',c('WBXX','WBXR','WBXY','WBRR','WBRY',
                          'ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- c(1-cM, 1+cM*(crM-1), 1, cM*crM, 1,
                                                                   (1-cbM)*(1-cM), (1-cbM)*(1+cM*(crM-1)), (1-cbM), (1-cbM)*cM*crM, (1-cbM),
                                                                   cbM*(1-cM), cbM*(1+cM*(crM-1)), cbM, cbM*cM*crM, cbM)/(8+4*cM*(crM-1))
  tMatrix['BBXR','WBXY',c('WBXX','WBXY','WBXR','WBRY',
                          'BBXX','BBXY','BBXR','BBRY')] <- c( 1, 1, 1, 1, 1, 1, 1, 1)/8
  tMatrix['BBXR','AAXY',c('ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- c((1-cbM)*(1-cM), (1-cbM)*(1+cM*(crM-1)), (1-cbM), (1-cbM)*cM*crM, (1-cbM),
                                                                    cbM*(1-cM), cbM*(1+cM*(crM-1)), cbM, cbM*cM*crM, cbM)/(4+2*cM*(crM-1))
  tMatrix['BBXR','ABXY',c('ABXX','ABXR','ABXY','ABRR','ABRY',
                          'BBXX','BBXR','BBXY','BBRR','BBRY')] <- c((1-cbM)*(1-cM), (1-cbM)*(1+cM*(crM-1)), (1-cbM), (1-cbM)*cM*crM, (1-cbM),
                                                                    (1+cbM)*(1-cM), (1+cbM)*(1+cM*(crM-1)), 1+cbM, (1+cbM)*cM*crM, 1+cbM)/(8+4*cM*(crM-1))
  tMatrix['BBXR','BBXY',c('BBXX','BBXY','BBXR','BBRY')] <- c( 1, 1, 1, 1)/4
  tMatrix['BBXR','WWRY',c('WBXR','WBXY','WBRR','WBRY')] <- c( 1, 1, 1, 1)/4
  tMatrix['BBXR','WARY',c('WBXR','WBXY','WBRR','WBRY',
                          'ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c( 1, 1, 1, 1,
                                                              1-cbM, 1-cbM, 1-cbM, 1-cbM,
                                                              cbM, cbM, cbM, cbM)/8
  tMatrix['BBXR','WBRY',c('WBXR','WBXY','WBRR','WBRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c( 1, 1, 1, 1,
                                                              1, 1, 1, 1)/8
  tMatrix['BBXR','AARY',c('ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c( 1-cbM, 1-cbM, 1-cbM, 1-cbM,
                                                              cbM, cbM, cbM, cbM)/4
  tMatrix['BBXR','ABRY',c('ABXR','ABXY','ABRR','ABRY',
                          'BBXR','BBXY','BBRR','BBRY')] <- c( 1-cbM, 1-cbM, 1-cbM, 1-cbM,
                                                              1+cbM, 1+cbM, 1+cbM, 1+cbM)/8
  tMatrix['BBXR','BBRY',c('BBXR','BBXY','BBRR','BBRY')] <- c( 1, 1, 1, 1)/4


  tMatrix['WWRR','WWXY',c('WWXR','WWRY')] <- c( 1, 1)/2
  tMatrix['WWRR','WAXY',c('WWXR','WWRY',
                          'WAXR','WARY',
                          'WBXR','WBRY')] <- c((1 + cM*(crM-1)), 1,
                                               (1-cbM)*(1 + cM*(crM-1)), 1-cbM,
                                               cbM*(1 + cM*(crM-1)), cbM)/(4 + 2*cM*(crM-1))
  tMatrix['WWRR','WBXY',c('WWXR','WWRY','WBXR','WBRY')] <- c( 1, 1, 1, 1)/4
  tMatrix['WWRR','AAXY',c('WAXR','WARY','WBXR','WBRY')] <- c((1-cbM)*(1 + cM*(crM-1)), 1-cbM, cbM*(1 + cM*(crM-1)), cbM)/(2 + cM*(crM-1))
  tMatrix['WWRR','ABXY',c('WAXR','WARY','WBXR','WBRY')] <- c((1-cbM)*(1 + cM*(crM-1)), 1-cbM, (1+cbM)*(1 + cM*(crM-1)), 1+cbM)/(4 + 2*cM*(crM-1))
  tMatrix['WWRR','BBXY',c('WBXR','WBRY')] <- c( 1, 1)/2
  tMatrix['WWRR','WWRY',c('WWRR','WWRY')] <- c( 1, 1)/2
  tMatrix['WWRR','WARY',c('WWRR','WWRY','WARR','WARY',
                          'WBRR','WBRY')] <- c( 1, 1, 1-cbM, 1-cbM, cbM, cbM)/4
  tMatrix['WWRR','WBRY',c('WWRR','WWRY','WBRR','WBRY')] <- c( 1, 1, 1, 1)/4
  tMatrix['WWRR','AARY',c('WARR','WARY','WBRR','WBRY')] <- c( 1-cbM, 1-cbM, cbM, cbM)/2
  tMatrix['WWRR','ABRY',c('WARR','WARY','WBRR','WBRY')] <- c( 1-cbM, 1-cbM, 1+cbM, 1+cbM)/4
  tMatrix['WWRR','BBRY',c('WBRR','WBRY')] <- c( 1, 1)/2


  tMatrix['WARR','WWXY',c('WWXR','WWRY','WAXR','WARY',
                          'WBXR','WBRY')] <- c( 1, 1, 1-cbF, 1-cbF, cbF, cbF)/4
  tMatrix['WARR','WAXY',c('WWXR','WWRY',
                          'WAXR','WARY',
                          'WBXR','WBRY',
                          'AAXR','AARY',
                          'ABXR','ABRY',
                          'BBXR','BBRY')] <- c((1+cM*(crM-1)), 1,
                                               (2-cbM-cbF)*(1+cM*(crM-1)), 2-cbM-cbF,
                                               (cbM+cbF)*(1+cM*(crM-1)), cbM+cbF,
                                               (1-cbF)*(1-cbM)*(1+cM*(crM-1)), (1-cbF)*(1-cbM),
                                               (cbF*(1-cbM) + cbM*(1-cbF))* (1+cM*(crM-1)), cbF*(1-cbM) + cbM*(1-cbF),
                                               cbF*cbM*(1+cM*(crM-1)), cbF*cbM)/(8 + 4*cM*(crM-1))
  tMatrix['WARR','WBXY',c('WWXR','WWRY','WAXR','WARY',
                          'WBXR','WBRY','ABXR','ABRY',
                          'BBXR','BBRY')] <- c( 1, 1, 1-cbF, 1-cbF,
                                                1+cbF, 1+cbF, 1-cbF, 1-cbF, cbF, cbF)/8
  tMatrix['WARR','AAXY',c('WAXR','WARY',
                          'WBXR','WBRY',
                          'AAXR','AARY',
                          'ABXR','ABRY',
                          'BBXR','BBRY')] <- c((1-cbM)*(1+cM*(crM-1)), 1-cbM,
                                               cbM*(1+cM*(crM-1)), cbM,
                                               (1-cbF)*(1-cbM)*(1+cM*(crM-1)), (1-cbF)*(1-cbM),
                                               (cbF*(1-cbM) + cbM*(1-cbF))* (1+cM*(crM-1)), cbF*(1-cbM) + cbM*(1-cbF),
                                               cbF*cbM*(1+cM*(crM-1)), cbF*cbM)/(4 + 2*cM*(crM-1))
  tMatrix['WARR','ABXY',c('WAXR','WARY',
                          'WBXR','WBRY',
                          'AAXR','AARY',
                          'ABXR','ABRY',
                          'BBXR','BBRY')] <- c((1-cbM)*(1+cM*(crM-1)), 1-cbM,
                                               (1+cbM)*(1+cM*(crM-1)), 1+cbM,
                                               (1-cbF)*(1-cbM)*(1+cM*(crM-1)), (1-cbF)*(1-cbM),
                                               (cbF*(1-cbM) + (1+cbM)*(1-cbF))* (1+cM*(crM-1)), cbF*(1-cbM) + (1+cbM)*(1-cbF),
                                               cbF*(1+cbM)*(1+cM*(crM-1)), cbF*(1+cbM))/(8 + 4*cM*(crM-1))
  tMatrix['WARR','BBXY',c('WBXR','WBRY','ABXR','ABRY',
                          'BBXR','BBRY')] <- c( 1, 1, 1-cbF, 1-cbF, cbF, cbF)/4
  tMatrix['WARR','WWRY',c('WWRR','WWRY','WARR','WARY',
                          'WBRR','WBRY')] <- c( 1, 1, 1-cbF, 1-cbF, cbF, cbF)/4
  tMatrix['WARR','WARY',c('WWRR','WWRY','WARR','WARY',
                          'WBRR','WBRY','AARR','AARY',
                          'ABRR','ABRY','BBRR','BBRY')] <- c( 1, 1, (1-cbF)+(1-cbM), (1-cbF)+(1-cbM),
                                                              cbF+cbM, cbF+cbM, (1-cbF)*(1-cbM), (1-cbF)*(1-cbM),
                                                              cbF*(1-cbM) + (1-cbF)*cbM, cbF*(1-cbM) + (1-cbF)*cbM, cbF*cbM, cbF*cbM)/8
  tMatrix['WARR','WBRY',c('WWRR','WWRY','WARR','WARY',
                          'WBRR','WBRY','ABRR','ABRY',
                          'BBRR','BBRY')] <- c( 1, 1, 1-cbF, 1-cbF,
                                                1+cbF, 1+cbF, 1-cbF, 1-cbF,
                                                cbF, cbF)/8
  tMatrix['WARR','AARY',c('WARR','WARY','AARR','AARY',
                          'ABRR','ABRY','WBRR','WBRY',
                          'BBRR','BBRY')] <- c( 1-cbM, 1-cbM, (1-cbF)*(1-cbM), (1-cbF)*(1-cbM),
                                                cbF*(1-cbM)+(1-cbF)*cbM, cbF*(1-cbM)+(1-cbF)*cbM, cbM, cbM,
                                                cbF*cbM, cbF*cbM)/4
  tMatrix['WARR','ABRY',c('WARR','WARY','AARR','AARY',
                          'ABRR','ABRY','WBRR','WBRY',
                          'BBRR','BBRY')] <- c( 1-cbM, 1-cbM, (1-cbF)*(1-cbM), (1-cbF)*(1-cbM),
                                                cbF*(1-cbM)+(1-cbF)*(1+cbM), cbF*(1-cbM)+(1-cbF)*(1+cbM), 1+cbM, 1+cbM,
                                                cbF*(1+cbM), cbF*(1+cbM))/8
  tMatrix['WARR','BBRY',c('WBRR','WBRY','ABRR','ABRY',
                          'BBRR','BBRY')] <- c( 1, 1, 1-cbF, 1-cbF, cbF, cbF)/4


  tMatrix['WBRR','WWXY',c('WWXR','WWRY','WBXR','WBRY')] <- c( 1, 1, 1, 1)/4
  tMatrix['WBRR','WAXY',c('WWXR','WWRY','WAXR','WARY',
                          'WBXR','WBRY','ABXR','ABRY',
                          'BBXR','BBRY')] <- c((1+cM*(crM-1)), 1,(1-cbM)*(1+cM*(crM-1)), 1-cbM,
                                               (1+cbM)*(1+cM*(crM-1)), 1+cbM,(1-cbM)*(1+cM*(crM-1)), 1-cbM,
                                               cbM*(1+cM*(crM-1)), cbM)/(8 + 4*cM*(crM-1))
  tMatrix['WBRR','WBXY',c('WWXR','WWRY','WBXR','WBRY',
                          'BBXR','BBRY')] <- c( 1, 1, 2, 2, 1, 1)/8
  tMatrix['WBRR','AAXY',c('WAXR','WARY','WBXR','WBRY',
                          'ABXR','ABRY','BBXR','BBRY')] <- c((1-cbM)*(1+cM*(crM-1)), 1-cbM,cbM*(1+cM*(crM-1)), cbM,
                                                             (1-cbM)*(1+cM*(crM-1)), 1-cbM,cbM*(1+cM*(crM-1)), cbM)/(4 + 2*cM*(crM-1))
  tMatrix['WBRR','ABXY',c('WAXR','WARY','WBXR','WBRY',
                          'ABXR','ABRY','BBXR','BBRY')] <- c((1-cbM)*(1+cM*(crM-1)), 1-cbM, (1+cbM)*(1+cM*(crM-1)), 1+cbM,
                                                             (1-cbM)*(1+cM*(crM-1)), 1-cbM, (1+cbM)*(1+cM*(crM-1)), 1+cbM)/(8 + 4*cM*(crM-1))
  tMatrix['WBRR','BBXY',c('WBXR','WBRY','BBXR','BBRY')] <- c( 1, 1, 1, 1)/4
  tMatrix['WBRR','WWRY',c('WWRR','WWRY','WBRR','WBRY')] <- c( 1, 1, 1, 1)/4
  tMatrix['WBRR','WARY',c('WWRR','WWRY','WBRR','WBRY',
                          'WARR','WARY','ABRR','ABRY',
                          'BBRR','BBRY')] <- c( 1, 1, 1+cbM, 1+cbM,
                                                1-cbM, 1-cbM, 1-cbM, 1-cbM,
                                                cbM, cbM)/8
  tMatrix['WBRR','WBRY',c('WWRR','WWRY','WBRR','WBRY',
                          'BBRR','BBRY')] <- c( 1, 1, 2, 2, 1, 1)/8
  tMatrix['WBRR','AARY',c('WARR','WARY','ABRR','ABRY',
                          'WBRR','WBRY','BBRR','BBRY')] <- c( 1-cbM, 1-cbM, 1-cbM, 1-cbM,
                                                              cbM, cbM, cbM, cbM)/4
  tMatrix['WBRR','ABRY',c('WARR','WARY','ABRR','ABRY',
                          'WBRR','WBRY','BBRR','BBRY')] <- c( 1-cbM, 1-cbM, 1-cbM, 1-cbM,
                                                              1+cbM, 1+cbM, 1+cbM, 1+cbM)/8
  tMatrix['WBRR','BBRY',c('WBRR','WBRY','BBRR','BBRY')] <- c( 1, 1, 1, 1)/4


  tMatrix['AARR','WWXY',c('WAXR','WARY','WBXR','WBRY')] <- c( 1-cbF, 1-cbF, cbF, cbF)/2
  tMatrix['AARR','WAXY',c('WAXR','WARY',
                          'WBXR','WBRY',
                          'AAXR','AARY',
                          'ABXR','ABRY',
                          'BBXR','BBRY')] <- c((1-cbF)*(1+cM*(crM-1)), 1-cbF,
                                               cbF*(1+cM*(crM-1)), cbF,
                                               (1-cbF)*(1-cbM)*(1+cM*(crM-1)), (1-cbF)*(1-cbM),
                                               (cbF*(1-cbM) + cbM*(1-cbF))* (1+cM*(crM-1)), cbF*(1-cbM) + cbM*(1-cbF),
                                               cbF*cbM*(1+cM*(crM-1)), cbF*cbM)/(4 + 2*cM*(crM-1))
  tMatrix['AARR','WBXY',c('WAXR','WARY','WBXR','WBRY',
                          'ABXR','ABRY','BBXR','BBRY')] <- c( 1-cbF, 1-cbF, cbF, cbF,
                                                               1-cbF, 1-cbF, cbF, cbF)/4
  tMatrix['AARR','AAXY',c('AAXR','AARY',
                          'ABXR','ABRY',
                          'BBXR','BBRY')] <- c((1-cbF)*(1-cbM)*(1+cM*(crM-1)), (1-cbF)*(1-cbM),
                                               (cbF*(1-cbM) + cbM*(1-cbF))* (1+cM*(crM-1)), cbF*(1-cbM) + cbM*(1-cbF),
                                               cbF*cbM*(1+cM*(crM-1)), cbF*cbM)/(2 + cM*(crM-1))
  tMatrix['AARR','ABXY',c('AAXR','AARY',
                          'ABXR','ABRY',
                          'BBXR','BBRY')] <- c((1-cbF)*(1-cbM)*(1+cM*(crM-1)), (1-cbF)*(1-cbM),
                                               (cbF*(1-cbM) + (1+cbM)*(1-cbF))* (1+cM*(crM-1)), cbF*(1-cbM) + (1+cbM)*(1-cbF),
                                               cbF*(1+cbM)*(1+cM*(crM-1)), cbF*(1+cbM))/(4 + 2*cM*(crM-1))
  tMatrix['AARR','BBXY',c('ABXR','ABRY','BBXR','BBRY')] <- c( 1-cbF, 1-cbF, cbF, cbF)/2
  tMatrix['AARR','WWRY',c('WARR','WARY','WBRR','WBRY')] <- c( 1-cbF, 1-cbF, cbF, cbF)/2
  tMatrix['AARR','WARY',c('WARR','WARY','WBRR','WBRY',
                          'AARR','AARY','ABRR','ABRY',
                          'BBRR','BBRY')] <- c( 1-cbF, 1-cbF, cbF, cbF,
                                                (1-cbF)*(1-cbM), (1-cbF)*(1-cbM), cbF*(1-cbM) + (1-cbF)*cbM, cbF*(1-cbM) + (1-cbF)*cbM,
                                                cbF*cbM, cbF*cbM)/4
  tMatrix['AARR','WBRY',c('WARR','WARY','WBRR','WBRY',
                          'ABRR','ABRY','BBRR','BBRY')] <- c( 1-cbF, 1-cbF, cbF, cbF,
                                                              1-cbF, 1-cbF, cbF, cbF)/4
  tMatrix['AARR','AARY',c('AARR','AARY','ABRR','ABRY',
                          'BBRR','BBRY')] <- c( (1-cbF)*(1-cbM), (1-cbF)*(1-cbM), cbF*(1-cbM) + (1-cbF)*cbM, cbF*(1-cbM) + (1-cbF)*cbM,
                                                cbF*cbM, cbF*cbM)/2
  tMatrix['AARR','ABRY',c('AARR','AARY','ABRR','ABRY',
                          'BBRR','BBRY')] <- c( (1-cbF)*(1-cbM), (1-cbF)*(1-cbM), cbF*(1-cbM) + (1-cbF)*(1+cbM), cbF*(1-cbM) + (1-cbF)*(1+cbM),
                                                cbF*(1+cbM), cbF*(1+cbM))/4
  tMatrix['AARR','BBRY',c('ABRR','ABRY','BBRR','BBRY')] <- c( 1-cbF, 1-cbF, cbF, cbF)/2


  tMatrix['ABRR','WWXY',c('WAXR','WARY','WBXR','WBRY')] <- c( 1-cbF, 1-cbF, 1+cbF, 1+cbF)/4
  tMatrix['ABRR','WAXY',c('WAXR','WARY',
                          'WBXR','WBRY',
                          'AAXR','AARY',
                          'ABXR','ABRY',
                          'BBXR','BBRY')] <- c((1-cbF)*(1+cM*(crM-1)), 1-cbF,
                                               (1+cbF)*(1+cM*(crM-1)), (1+cbF),
                                               (1-cbF)*(1-cbM)*(1+cM*(crM-1)), (1-cbF)*(1-cbM),
                                               ((1+cbF)*(1-cbM) + cbM*(1-cbF))* (1+cM*(crM-1)), (1+cbF)*(1-cbM) + cbM*(1-cbF),
                                               (1+cbF)*cbM*(1+cM*(crM-1)), (1+cbF)*cbM)/(8 + 4*cM*(crM-1))
  tMatrix['ABRR','WBXY',c('WAXR','WARY','WBXR','WBRY',
                          'ABXR','ABRY','BBXR','BBRY')] <- c( 1-cbF, 1-cbF, 1+cbF, 1+cbF,
                                                               1-cbF, 1-cbF, 1+cbF, 1+cbF)/8
  tMatrix['ABRR','AAXY',c('AAXR','AARY',
                          'ABXR','ABRY',
                          'BBXR','BBRY')] <- c((1-cbF)*(1-cbM)*(1+cM*(crM-1)), (1-cbF)*(1-cbM),
                                               ((1+cbF)*(1-cbM) + cbM*(1-cbF))* (1+cM*(crM-1)), (1+cbF)*(1-cbM) + cbM*(1-cbF),
                                               (1+cbF)*cbM*(1+cM*(crM-1)), (1+cbF)*cbM)/(4 + 2*cM*(crM-1))
  tMatrix['ABRR','ABXY',c('AAXR','AARY',
                          'ABXR','ABRY',
                          'BBXR','BBRY')] <- c((1-cbF)*(1-cbM)*(1+cM*(crM-1)), (1-cbF)*(1-cbM),
                                               ((1+cbF)*(1-cbM) + (1+cbM)*(1-cbF))* (1+cM*(crM-1)), (1+cbF)*(1-cbM) + (1+cbM)*(1-cbF),
                                               (1+cbF)*(1+cbM)*(1+cM*(crM-1)), (1+cbF)*(1+cbM))/(8 + 4*cM*(crM-1))
  tMatrix['ABRR','BBXY',c('ABXR','ABRY','BBXR','BBRY')] <- c( 1-cbF, 1-cbF, 1+cbF, 1+cbF)/4
  tMatrix['ABRR','WWRY',c('WARR','WARY','WBRR','WBRY')] <- c( 1-cbF, 1-cbF, 1+cbF, 1+cbF)/4
  tMatrix['ABRR','WARY',c('WARR','WARY','WBRR','WBRY',
                          'AARR','AARY','ABRR','ABRY',
                          'BBRR','BBRY')] <- c( 1-cbF, 1-cbF, 1+cbF, 1+cbF,
                                                (1-cbF)*(1-cbM), (1-cbF)*(1-cbM), (1+cbF)*(1-cbM) + (1-cbF)*cbM, (1+cbF)*(1-cbM) + (1-cbF)*cbM,
                                                (1+cbF)*cbM, (1+cbF)*cbM)/8
  tMatrix['ABRR','WBRY',c('WARR','WARY','WBRR','WBRY',
                          'ABRR','ABRY','BBRR','BBRY')] <- c( 1-cbF, 1-cbF, 1+cbF, 1+cbF,
                                                              1-cbF, 1-cbF, 1+cbF, 1+cbF)/8
  tMatrix['ABRR','AARY',c('AARR','AARY','ABRR','ABRY',
                          'BBRR','BBRY')] <- c( (1-cbF)*(1-cbM), (1-cbF)*(1-cbM), (1+cbF)*(1-cbM) + (1-cbF)*cbM, (1+cbF)*(1-cbM) + (1-cbF)*cbM,
                                                (1+cbF)*cbM, (1+cbF)*cbM)/4
  tMatrix['ABRR','ABRY',c('AARR','AARY','ABRR','ABRY',
                          'BBRR','BBRY')] <- c( (1-cbF)*(1-cbM), (1-cbF)*(1-cbM), (1+cbF)*(1-cbM) + (1-cbF)*(1+cbM), (1+cbF)*(1-cbM) + (1-cbF)*(1+cbM),
                                                (1+cbF)*(1+cbM), (1+cbF)*(1+cbM))/8
  tMatrix['ABRR','BBRY',c('ABRR','ABRY','BBRR','BBRY')] <- c( 1-cbF, 1-cbF, 1+cbF, 1+cbF)/4


  tMatrix['BBRR','WWXY',c('WBXR','WBRY')] <- c( 1, 1)/2
  tMatrix['BBRR','WAXY',c('WBXR','WBRY',
                          'ABXR','ABRY',
                          'BBXR','BBRY')] <- c((1 + cM*(crM-1)), 1,
                                              (1-cbM)*(1 + cM*(crM-1)), 1-cbM,
                                              cbM*(1 + cM*(crM-1)), cbM)/(4 + 2*cM*(crM-1))
  tMatrix['BBRR','WBXY',c('WBXR','WBRY','BBXR','BBRY')] <- c( 1, 1, 1, 1)/4
  tMatrix['BBRR','AAXY',c('ABXR','ABRY','BBXR','BBRY')] <- c((1-cbM)*(1 + cM*(crM-1)), 1-cbM, cbM*(1 + cM*(crM-1)), cbM)/(2 + cM*(crM-1))
  tMatrix['BBRR','ABXY',c('ABXR','ABRY','BBXR','BBRY')] <- c((1-cbM)*(1 + cM*(crM-1)), 1-cbM, (1+cbM)*(1 + cM*(crM-1)), 1+cbM)/(4 + 2*cM*(crM-1))
  tMatrix['BBRR','BBXY',c('BBXR','BBRY')] <- c( 1, 1)/2
  tMatrix['BBRR','WWRY',c('WBRR','WBRY')] <- c( 1, 1)/2
  tMatrix['BBRR','WARY',c('WBRR','WBRY','ABRR','ABRY',
                          'BBRR','BBRY')] <- c( 1, 1, 1-cbM, 1-cbM, cbM, cbM)/4
  tMatrix['BBRR','WBRY',c('WBRR','WBRY','BBRR','BBRY')] <- c( 1, 1, 1, 1)/4
  tMatrix['BBRR','AARY',c('ABRR','ABRY','BBRR','BBRY')] <- c( 1-cbM, 1-cbM, cbM, cbM)/2
  tMatrix['BBRR','ABRY',c('ABRR','ABRY','BBRR','BBRY')] <- c( 1-cbM, 1-cbM, 1+cbM, 1+cbM)/4
  tMatrix['BBRR','BBRY',c('BBRR','BBRY')] <- c( 1, 1)/2


  ## take care of NaNs
  # As the shredding goes to 100%, the fecundity will go to 0, and some of the cube
  # entries will tend towards infinity, which can't be masked by s parameter
  # So, if the fecundity tends to 0, we will just set those cube entries to zero
  #  which gets handled (skipped) in the simulation
  fecundityXX <- 1 + cF*(crF - 1)
  fecundityXR <- 1 + 0.5*cF*(crF - 1)

  if(fecundityXX < sqrt(.Machine$double.eps)){
    # if small, set to 0
    fecundityXX <- 0
    # set cube entries that suffer to 0
    tMatrix[c('WAXX','AAXX','WBXX'), , ] <- 0
  }


  ## initialize viability mask.
  viabilityMask <- array(data = 1, dim = c(size,size,size), dimnames = list(gtype, gtype, gtype))


  ## genotype-specific modifiers
  if(!is.null(phi)){
    stop("This cube has a special phi, due to being male/female specific.
         Please edit it after the cube is built, if you are absolutely sure it needs to change.")
  }
  if(!is.null(s)){
    stop("This cube has a special fecundity (s), due to the X-shredding in females.
         Please edit it after the cube is built, if you are absolutely sure it needs to change.")
  }

  s = c('WAXX'=fecundityXX,'AAXX'=fecundityXX,'ABXX'=fecundityXX,
        'WAXR'=fecundityXR,'AAXR'=fecundityXR,'ABXR'=fecundityXR)
  phi = setNames(object = c(rep.int(x = 1, times = 18),rep.int(x = 0, times = 12)), nm = gtype)

  modifiers = cubeModifiers(gtype, eta = eta, phi = phi, omega = omega, xiF = xiF, xiM = xiM, s = s)


  ## put everything into a labeled list to return
  return(list(
    ih = tMatrix,
    tau = viabilityMask,
    genotypesID = gtype,
    genotypesN = size,
    wildType = c('WWXX','WWXY'),
    eta = modifiers$eta,
    phi = modifiers$phi,
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = 'AAXY' # these are males only. Females would be AAXX, but then have 0 offspring
  ))
}
