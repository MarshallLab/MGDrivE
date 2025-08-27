###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   Split-Drive for Shih-Che Weng (s2weng@ucsd.edu)
#   Héctor Sanchez, Jared Bennett, John Marshall
#   jared_bennett@berkeley.edu
#   June 2025
#
###############################################################################

#' Inheritance Cube: 3-Piece Allele Sail
#'
#' A generalized implementation of the [Allele Sail](https://doi.org/10.1038/s41467-024-50992-9)
#' idea.
#'
#' This is an autosomal, 3-locus system. The first locus contains the Cas9
#' allele, the second locus carries the gRNA, and the third locus is the target.
#' All loci can be linked/unlinked to the locus before it (so, 1 to 2 or 2 to 3).
#' Cas9 efficacy due to provenance (mother vs father) is included.
#'
#' This construct is very similar to our [2-locus Cleave and Rescue][MGDrivE::cubeClvR2()] design for
#' [Oberhofer](https://doi.org/10.1101/2020.07.09.196253).
#'
#' This construct has 3 alleles at the first locus, 2 alleles at the second locus,
#' and 3 alleles at the third locus.
#'  * Locus 1
#'    * W: Wild-type
#'    * P: Paternal Cas9
#'    * M: Maternal Cas9
#'  * Locus 2
#'    * W: Wild-type
#'    * G: gRNAs
#'  * Locus 3
#'    * W: Wild-type
#'    * R: Resistant 1
#'    * B: Resistant 2
#'
#' Female deposition is implemented incorrectly. Right now, it is performed on
#' male alleles prior to zygote formation - it should happen post-zygote formation.
#' Since this construct doesn't have HDR, this should be fine. \cr
#'
#' @param cMM Cutting efficacy of maternally-inherited Cas9 in males
#' @param crMM Resistance rate of maternally-inherited Cas9 in males
#'
#' @param cPM Cutting efficacy of paternally-inherited Cas9 in males
#' @param crPM Resistance rate of paternally-inherited Cas9 in males
#'
#' @param cMF Cutting efficacy of maternally-inherited Cas9 in females
#' @param crMF Resistance rate of maternally-inherited Cas9 in females
#'
#' @param cPF Cutting efficacy of paternally-inherited Cas9 in females
#' @param crPF Resistance rate of paternally-inherited Cas9 in females
#'
#' @param dMW Female deposition cutting rate, maternal Cas9
#' @param dMrW Female deposition functional resistance rate, maternal Cas9
#' @param dPW Female deposition (HH) cutting rate, paternal Cas9
#' @param dPrW Female deposition (HH) functional resistance rate, paternal Cas9
#'
#' @param crF12 Female crossover rate between loci 1 and 2, 0 is completely linked, 0.5 is unlinked, 1.0 is complete divergence
#' @param crM12 Male crossover rate between loci 1 and 2, 0 is completely linked, 0.5 is unlinked, 1.0 is complete divergence
#' @param crF23 Female crossover rate between loci 2 and 3, 0 is completely linked, 0.5 is unlinked, 1.0 is complete divergence
#' @param crM23 Male crossover rate between loci 2 and 3, 0 is completely linked, 0.5 is unlinked, 1.0 is complete divergence
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
cubeAlleleSail <- function(cMM=0, crMM=0,
                         cPM=0, crPM=0,
                         cMF=0, crMF=0,
                         cPF=0, crPF=0,
                         dMW=0, dMrW=0, dPW=0, dPrW=0,
                         crF12=0.5, crM12=0.5, crF23=0.5, crM23=0.5,
                         eta=NULL, phi=NULL,omega=NULL, xiF=NULL, xiM=NULL, s=NULL){

  # This was stolen from the ASmidler cube.
  # Super similar to ClvR2 cube.
  # Same deposition warnings - not sure there's a way to fix them in this kind of
  #  cube build implementation.
  #  For ease-of-finding, the deposition issues is that deposition should occur
  #   in the zygote, when it has both copies of its alleles. This defines the
  #   repair pattern.
  #   In this implementation, deposition occurs in the male gamete alone, so the
  #   "repair" process is wrong.
  #   In this drive, we have no HDR, just cleavage, so no sister allele to repair
  #   from, therefore this is fine.
  # IDK why previous ones had copy-number-dependent actions - we've never known
  #  that info.
  #  Not gonna consider it here.


  # # for testing
  # testVec <- numeric(16)+1 #runif(n = 16)
  # cMM <- testVec[1]
  # crMM <- testVec[2]
  # cPM <- testVec[3]
  # crPM <- testVec[4]
  # cMF <- testVec[5]
  # crMF <- testVec[6]
  # cPF <- testVec[7]
  # crPF <- testVec[8]
  # dMW <- 0 #testVec[9]
  # dMrW <- testVec[10]
  # dPW <- 0 #testVec[11]
  # dPrW <- testVec[12]
  # crF12 <- 0 #testVec[13]
  # crM12 <- 0 #testVec[14]
  # crF23 <- 0 #testVec[15]
  # crM23 <- 0 #testVec[16]

  # # this would run at the bottom, to check that all cells sum to 1
  # #  numGen is defined below
  # for(mID in 1:numGen){
  #   for(fID in 1:numGen){
  #     if((abs(sum(tMatrix[fID,mID,])-1)) > 1e-6){
  #       print(paste0("Fail: mID= ",gtype[mID],", fID= ",gtype[fID]))
  #     }
  #   }
  # }


  ## safety checks
  params <- c(cMM,crMM,cPM,crPM,cMF,crMF,cPF,crPF,
              dMW,dMrW,dPW,dPrW,crF12,crM12,crF23,crM23)
  if(any(params>1) || any(params<0)){
    stop('Parameters are rates.\n0 <= x <= 1')
  }


  #############################################################################
  ## generate all genotypes, set up vectors and matrices
  #############################################################################
  #list of possible alleles at each locus
  # the first locus has the Cas9
  # the second locus has the gRNAs
  # the third locus has the target site
  gTypes <- list(c('W', 'P', 'M'), c('W', 'G'), c('W', 'R', 'B'))


  # # generate alleles
  # # this generates all combinations of Target Sites 1 and 2 for one allele,
  # #  then mixes alleles, so each half is a different allele.
  # # expand combinations of target site 1 and 2 on one allele
  # hold <- expand.grid(gTypes,KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  # # paste them all together
  # hold <- do.call(what = paste0, args = list(hold[,1], hold[,2], hold[,3]))
  # #get all combinations of both loci
  # openAlleles <- expand.grid(hold, hold, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  # # sort both alleles to remove uniques
  # # ie, each target site is unique and can't be moved
  # # but, the alleles aren't unique, and can occur in any order
  # openAlleles <- apply(X = openAlleles, MARGIN = 1, FUN = sort)
  # # paste together and get unique ones only
  # genotypes <- unique(do.call(what = paste0, args = list(openAlleles[1, ], openAlleles[2, ])))
  # # you can't get both alleles from mother or both from father, it must be mixed
  # #  this removes those doubles
  # genotypes <- grep(pattern = 'P..P..|M..M..', x = genotypes, value = TRUE, invert = TRUE)

  genotypes <- c('WWWWWW','PWWWWW','MWWWWW','WGWWWW','PGWWWW','MGWWWW','WWRWWW',
                'PWRWWW','MWRWWW','WGRWWW','PGRWWW','MGRWWW','WWBWWW','PWBWWW',
                'MWBWWW','WGBWWW','PGBWWW','MGBWWW','MWWPWW','PWWWGW','MGWPWW',
                'PWWWWR','MWRPWW','PWWWGR','MGRPWW','PWWWWB','MWBPWW','PWWWGB',
                'MGBPWW','MWWWGW','MWWPGW','MWWWWR','MWWPWR','MWWWGR','MWWPGR',
                'MWWWWB','MWWPWB','MWWWGB','MWWPGB','WGWWGW','PGWWGW','MGWWGW',
                'WGWWWR','PWRWGW','MWRWGW','WGRWGW','PGRWGW','MGRWGW','WGWWWB',
                'PWBWGW','MWBWGW','WGBWGW','PGBWGW','MGBWGW','MGWPGW','PGWWWR',
                'MWRPGW','PGWWGR','MGRPGW','PGWWWB','MWBPGW','PGWWGB','MGBPGW',
                'MGWWWR','MGWPWR','MGWWGR','MGWPGR','MGWWWB','MGWPWB','MGWWGB',
                'MGWPGB','WWRWWR','PWRWWR','MWRWWR','WGRWWR','PGRWWR','MGRWWR',
                'WWBWWR','PWBWWR','MWBWWR','WGBWWR','PGBWWR','MGBWWR','MWRPWR',
                'PWRWGR','MGRPWR','PWRWWB','MWBPWR','PWRWGB','MGBPWR','MWRWGR',
                'MWRPGR','MWRWWB','MWRPWB','MWRWGB','MWRPGB','WGRWGR','PGRWGR',
                'MGRWGR','WGRWWB','PWBWGR','MWBWGR','WGBWGR','PGBWGR','MGBWGR',
                'MGRPGR','PGRWWB','MWBPGR','PGRWGB','MGBPGR','MGRWWB','MGRPWB',
                'MGRWGB','MGRPGB','WWBWWB','PWBWWB','MWBWWB','WGBWWB','PGBWWB',
                'MGBWWB','MWBPWB','PWBWGB','MGBPWB','MWBWGB','MWBPGB','WGBWGB',
                'PGBWGB','MGBWGB','MGBPGB')


  #############################################################################
  ## setup all probability lists
  #############################################################################
  # loci 1 and 2 are both Mendelian only
  #  So, we don't need separate things for them

  # Female Mendelian
  fMendelian <- list('W'=c('W'=1),
                     'P'=c('M'=1),
                     'M'=c('M'=1),
                     'G'=c('G'=1),
                     'R'=c('R'=1),
                     'B'=c('B'=1))

  # Male Mendelian
  mMendelian <- list('W'=c('W'=1),
                     'P'=c('P'=1),
                     'M'=c('P'=1),
                     'G'=c('G'=1),
                     'R'=c('R'=1),
                     'B'=c('B'=1))

  # Female Locus 3
  fLocus3 <- list()

  # the nothing case is Mendelian
  # we don't do deposition in female alleles - assumption is that you can't
  #  differentiate it from proper cleavage, b/c everything is there within the egg.
  #  Ethan has made some noise about distribution of repair processes at each
  #  stage, but ... meh, ignoring that >.<
  fLocus3$Maternal <- list('W'=c('W'=1-cMF,
                                'R'=cMF*crMF,
                                'B'=cMF*(1-crMF)),
                          'R'=c('R'=1),
                          'B'=c('B'=1))

  fLocus3$Paternal <- list('W'=c('W'=1-cPF,
                                'R'=cPF*crPF,
                                'B'=cPF*(1-crPF)),
                          'R'=c('R'=1),
                          'B'=c('B'=1))

  # Male Locus 3
  mLocus3 <- list()

  # No male homing
  #  the nothing case is Mendelian
  mLocus3$MDep <- list('W'=c('W'=1-dMW,
                                'R'=dMW*dMrW,
                                'B'=dMW*(1-dMrW)),
                          'R'=c('R'=1),
                          'B'=c('B'=1))
  mLocus3$PDep <- list('W'=c('W'=1-dPW,
                                'R'=dPW*dPrW,
                                'B'=dPW*(1-dPrW)),
                          'R'=c('R'=1),
                          'B'=c('B'=1))

  # Maternal - 3 cases
  mLocus3$M <- list('W'=c('W'=1-cMM,
                          'R'=cMM*crMM,
                          'B'=cMM*(1-crMM)),
                    'R'=c('R'=1),
                    'B'=c('B'=1))
  mLocus3$MMDep <- list('W'=c('W'=(1-cMM)*(1-dMW),
                                'R'=cMM*crMM + (1-cMM)*dMW*dMrW,
                                'B'=cMM*(1-crMM) + (1-cMM)*dMW*(1-dMrW)),
                          'R'=c('R'=1),
                          'B'=c('B'=1))
  mLocus3$MPDep <- list('W'=c('W'=(1-cMM)*(1-dPW),
                                'R'=cMM*crMM + (1-cMM)*dPW*dPrW,
                                'B'=cMM*(1-crMM) + (1-cMM)*dPW*(1-dPrW)),
                          'R'=c('R'=1),
                          'B'=c('B'=1))

  # Paternal - 3 cases
  mLocus3$P <- list('W'=c('W'=1-cPM,
                          'R'=cPM*crPM,
                          'B'=cPM*(1-crPM)),
                    'R'=c('R'=1),
                    'B'=c('B'=1))
  mLocus3$PMDep <- list('W'=c('W'=(1-cPM)*(1-dMW),
                          'R'=cPM*crPM + (1-cPM)*dMW*dMrW,
                          'B'=cPM*(1-crPM) + (1-cPM)*dMW*(1-dMrW)),
                    'R'=c('R'=1),
                    'B'=c('B'=1))
  mLocus3$PPDep <- list('W'=c('W'=(1-cPM)*(1-dMW),
                          'R'=cPM*crPM + (1-cPM)*dPW*dPrW,
                          'B'=cPM*(1-crPM) + (1-cPM)*dPW*(1-dPrW)),
                    'R'=c('R'=1),
                    'B'=c('B'=1))


  #############################################################################
  ## fill transition matrix
  #############################################################################
  #use this many times down below
  numGen <- length(genotypes)

  #create transition matrix to fill
  tMatrix <- array(data = 0,
                   dim = c(numGen,numGen,numGen),
                   dimnames = list(genotypes,genotypes,genotypes))

  #number of alleles, set score vectors
  numLoci <- length(gTypes)
  numAlleles <- 2
  fScore <- mScore <- logical(length = numAlleles)


  #############################################################################
  ## loop over all matings, female outer loop
  #############################################################################
  for(fi in 1:numGen){
    #do female stuff here
    # This splits all characters.
    fSplit <- strsplit(x = genotypes[fi], split = "", useBytes = TRUE)[[1]]
    # Matrix of alleles and target sites
    #  Each row is an allele, there are 2 alleles
    #  Each column is the target site on that allele, there are 3 target sites
    #  this is because target sites are linked on an allele
    fAlleleMat <- matrix(data = fSplit, nrow = numAlleles, ncol = numLoci, byrow = TRUE)
    #Score them
    fMAllele <- any("M" == fAlleleMat)
    fScore[1] <- any("P" == fAlleleMat) || fMAllele
    fScore[2] <- any("G" == fAlleleMat)

    #setup offspring allele lists
    fPHold <- rep(x = list(list()), times = numAlleles)
    fAllele <- character(length = 0L)
    fProbs <- numeric(length = 0L)

    ##########
    # Base Inheritance
    ##########
    # check if there will be homing
    if(fScore[1] && fScore[2]){
      # There is a Cas9 - either maternal or paternal grandparent of origin
      # There are gRNAs
      # Thus, there is homing

      # loop over alleles, must keep target sites linked
      for(allele in 1:numAlleles){
        # loci 1/2 are always mendelian
        fPHold[[allele]][[1]] <- fMendelian[[ fAlleleMat[allele,1] ]]
        fPHold[[allele]][[2]] <- fMendelian[[ fAlleleMat[allele,2] ]]

        # check if maternal homing
        if(fMAllele){
          # at least one of the Cas9 proteins is from the grandmother
          fPHold[[allele]][[3]] <- fLocus3$Maternal[[ fAlleleMat[allele,3] ]]
        } else {
          # the Cas9 protein is from the grandfather
          fPHold[[allele]][[3]] <- fLocus3$Paternal[[ fAlleleMat[allele,3] ]]
        } # end M/P homing

      } # end loop over alleles

    } else {
      # either no Cas9 or no gRNAs
      # all inheritance is Mendelian

      # loop over alleles, must keep target sites linked
      for(allele in 1:numAlleles){
        # set target site 1, then target site 2
        fPHold[[allele]][[1]] <- fMendelian[[ fAlleleMat[allele,1] ]]
        fPHold[[allele]][[2]] <- fMendelian[[ fAlleleMat[allele,2] ]]
        fPHold[[allele]][[3]] <- fMendelian[[ fAlleleMat[allele,3] ]]
      } # end loop over alleles

    } # end female scoring

    ##########
    # perform cross-overs
    ##########
    # holder objects
    hold1 <- fPHold[[1]][[1]]
    hold3 <- fPHold[[1]][[3]]
    # crossover between loci 1 and 2
    fPHold[[1]][[1]] <- c((1-crF12)*hold1, crF12*fPHold[[2]][[1]])
    fPHold[[2]][[1]] <- c((1-crF12)*fPHold[[2]][[1]], crF12*hold1)
    # crossover between loci 2 and 3 (well, 3 and 2 is how it's written)
    fPHold[[1]][[3]] <- c((1-crF23)*hold3, crF23*fPHold[[2]][[3]])
    fPHold[[2]][[3]] <- c((1-crF23)*fPHold[[2]][[3]], crF23*hold3)

    ##########
    # all combinations of female alleles
    ##########
    for(allele in 1:numAlleles){
      # make combinations of the allele, then store those combinations to mix
      #  with the next allele
      # expand combinations
      holdProbs <- expand.grid(fPHold[[allele]],KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      holdAllele <- expand.grid(names(fPHold[[allele]][[1]]),
                                names(fPHold[[allele]][[2]]),
                                names(fPHold[[allele]][[3]]),
                                KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

      # contract (paste or multiply) combinations and store
      fProbs <- c(fProbs, holdProbs[,1]*holdProbs[,2]*holdProbs[,3])
      fAllele <- c(fAllele, file.path(holdAllele[,1], holdAllele[,2], holdAllele[,3], fsep = ""))
    }

    # remove zeros
    fAllele <- fAllele[fProbs!=0]
    fProbs <- fProbs[fProbs!=0]


    ###########################################################################
    ## loop over male mate. This is the inner loop
    ###########################################################################
    # we make the assumption that maternal homing is more important than paternal
    # ie, if we have M and P at locus one on both alleles, the M takes precedence
    #  in terms of homing.
    for(mi in 1:numGen){
      # isn't symmetric

      #do male stuff here
      # split male genotype
      # This splits all characters.
      mSplit <- strsplit(x = genotypes[mi], split = "", useBytes = TRUE)[[1]]
      # Matrix of alleles and target sites
      #  Each row is an allele, there are 2 alleles
      #  Each column is the target site on that allele, there are 3 target sites
      #  this is because target sites are linked on an allele
      mAlleleMat <- matrix(data = mSplit, nrow = numAlleles, ncol = numLoci, byrow = TRUE)
      #Score them
      mMAllele <- any("M" == mAlleleMat)
      mScore[1] <- any("P" == mAlleleMat) || mMAllele
      mScore[2] <- any("G" == mAlleleMat)

      #setup offspring allele lists
      mPHold <- rep(x = list(list()), times = numAlleles)
      mAllele <- character(length = 0L)
      mProbs <- numeric(length = 0L)

      ##########
      # Base Inheritance
      ##########
      # check if there will be homing
      if(mScore[1] && mScore[2]){
        # Father has Cas9 - either maternal or paternal grandparent of origin
        # There are gRNAs
        # Thus, there is homing
        # check for maternal deposition

        # 2 possible scenarios
        # mother has cas9 (father already has gRNA)
        # mother does not have cas9 - mendelian

        # loop over alleles, must keep target sites linked
        for(allele in 1:numAlleles){
          # loci 1/2 are always mendelian
          mPHold[[allele]][[1]] <- mMendelian[[ mAlleleMat[allele,1] ]]
          mPHold[[allele]][[2]] <- mMendelian[[ mAlleleMat[allele,2] ]]

          # check for deposition
          if(fScore[1]){
            # Deposition occurs

            # check origin of Cas9 for deposition
            if(fMAllele){
              # Maternal deposition, Cas9 from grandmother

              # check Cas9 origin
              if(mMAllele){
                # at least one of the Cas9 proteins is from the grandmother
                mPHold[[allele]][[3]] <- mLocus3$MMDep[[ mAlleleMat[allele,3] ]]
              } else {
                # the Cas9 protein is from the grandfather
                mPHold[[allele]][[3]] <- mLocus3$PMDep[[ mAlleleMat[allele,3] ]]
              } # end M/P homing

            } else {
              # Maternal deposition, Cas9 from grandfather

              # check Cas9 origin
              if(mMAllele){
                # at least one of the Cas9 proteins is from the grandmother
                mPHold[[allele]][[3]] <- mLocus3$MPDep[[ mAlleleMat[allele,3] ]]
              } else {
                # the Cas9 protein is from the grandfather
                mPHold[[allele]][[3]] <- mLocus3$PPDep[[ mAlleleMat[allele,3] ]]
              } # end M/P homing

            } # end deposition origin check

          } else {
            # No deposition
            # paternal homing only

            # check Cas9 origin
            if(mMAllele){
              # at least one of the Cas9 proteins is from the grandmother
              mPHold[[allele]][[3]] <- mLocus3$M[[ mAlleleMat[allele,3] ]]
            } else {
              # the Cas9 protein is from the grandfather
              mPHold[[allele]][[3]] <- mLocus3$P[[ mAlleleMat[allele,3] ]]
            } # end M/P homing

          } # end deposition check
        } # end loop over alleles

      } else {
        # father has either no Cas9 or no gRNAs
        # inheritance is Mendelian
        # check for maternal deposition

        # 3 possible scenarios
        # mother has cas9 and gRNA (father doesn't matter)
        # mother has cas9, father has gRNA
        # mother does not have cas9 - mendelian

        # check for deposition
        if(fScore[1] && (fScore[2] || mScore[2])){
          # Deposition occurs
          #  Mother has cas9
          #  Mother or father has gRNA

          # loop over alleles, must keep target sites linked
          for(allele in 1:numAlleles){
            # loci 1/2 are always mendelian
            mPHold[[allele]][[1]] <- mMendelian[[ mAlleleMat[allele,1] ]]
            mPHold[[allele]][[2]] <- mMendelian[[ mAlleleMat[allele,2] ]]

            # check origin of Cas9
            if(fMAllele){
              # Maternal deposition, Cas9 from grandmother
              mPHold[[allele]][[3]] <- mLocus3$MDep[[ mAlleleMat[allele,3] ]]
            } else {
              # Maternal deposition, Cas9 from grandfather
              mPHold[[allele]][[3]] <- mLocus3$PDep[[ mAlleleMat[allele,3] ]]
            } # end deposition origin check

          } # end loop over alleles

        } else {
          # No deposition

          # loop over alleles, must keep target sites linked
          for(allele in 1:numAlleles){
            mPHold[[allele]][[1]] <- mMendelian[[ mAlleleMat[allele,1] ]]
            mPHold[[allele]][[2]] <- mMendelian[[ mAlleleMat[allele,2] ]]
            mPHold[[allele]][[3]] <- mMendelian[[ mAlleleMat[allele,3] ]]
          } # end loop over alleles

        } # end deposition-only check

      } # end male scoring

      ##########
      # perform cross-overs
      ##########
      # holder objects
      hold1 <- mPHold[[1]][[1]]
      hold3 <- mPHold[[1]][[3]]
      # crossover between loci 1 and 2
      mPHold[[1]][[1]] <- c((1-crM12)*hold1, crM12*mPHold[[2]][[1]])
      mPHold[[2]][[1]] <- c((1-crM12)*mPHold[[2]][[1]], crM12*hold1)
      # crossover between loci 3 and 2
      mPHold[[1]][[3]] <- c((1-crM23)*hold3, crM23*mPHold[[2]][[3]])
      mPHold[[2]][[3]] <- c((1-crM23)*mPHold[[2]][[3]], crM23*hold3)

      ##########
      # all combinations of male alleles
      ##########
      for(allele in 1:numAlleles){
        # make combinations of the allele, then store those combinations to mix
        #  with the next allele
        # expand combinations
        holdProbs <- expand.grid(mPHold[[allele]],KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        holdAllele <- expand.grid(names(mPHold[[allele]][[1]]),
                                  names(mPHold[[allele]][[2]]),
                                  names(mPHold[[allele]][[3]]),
                                  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

        # contract (paste or multiply) combinations and store
        mProbs <- c(mProbs, holdProbs[,1]*holdProbs[,2]*holdProbs[,3])
        mAllele <- c(mAllele, file.path(holdAllele[,1], holdAllele[,2], holdAllele[,3], fsep = ""))
      }

      # remove zeros
      mAllele <- mAllele[mProbs!=0]
      mProbs <- mProbs[mProbs!=0]


      #########################################################################
      ## Get combinations and put them in the tMatrix. This must be done
      ##  inside the inner loop
      #########################################################################
      # male and female alleles/probs are already combined by target sites, and
      #  we assume alleles segregate independently, so we just have to get combinations
      #  of male vs female for the offspring
      #  NOTE: this is why deposition is technically wrong in cubes implemented like this.
      holdProbs <- expand.grid(fProbs, mProbs, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      holdAllele <- expand.grid(fAllele, mAllele, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

      # sort alleles into order
      holdAllele <- apply(X = holdAllele, MARGIN = 1, FUN = sort.int)

      # contract (paste or multiply) combinations
      holdProbs <- holdProbs[ ,1]*holdProbs[ ,2]
      holdAllele <- file.path(holdAllele[1,], holdAllele[2,], fsep = "")

      #aggregate duplicate genotypes
      reducedP <- vapply(X = unique(holdAllele), FUN = function(x){
        sum(holdProbs[holdAllele==x])},
        FUN.VALUE = numeric(length = 1L))

      #normalize
      reducedP <- reducedP/sum(reducedP)

      #set values in tMatrix
      tMatrix[fi,mi, names(reducedP) ] <- reducedP

    }# end male loop
  }# end female loop

  ## set the other half of the matrix
  tMatrix[tMatrix < .Machine$double.eps] <- 0 #protection from underflow errors

  ## initialize viability mask. No mother-specific death.
  viabilityMask <- array(data = 1, dim = c(numGen,numGen,numGen),
                         dimnames = list(genotypes, genotypes, genotypes))

  ## genotype-specific modifiers
  modifiers = cubeModifiers(genotypes, eta = eta, phi = phi,
                            omega = omega,xiF = xiF, xiM = xiM, s = s)

  ## put everything into a labeled list to return
  return(list(
    ih = tMatrix,
    tau = viabilityMask,
    genotypesID = genotypes,
    genotypesN = numGen,
    wildType = "WWWWWW",
    eta = modifiers$eta,
    phi = modifiers$phi,
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = "MGWPGW"
  ))

}
