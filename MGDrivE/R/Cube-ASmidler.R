###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   Split-Drive for Andrea Smidler (asmidler@ucsd.edu)
#   Héctor Sanchez, Jared Bennett, Sean Wu, John Marshall
#   jared_bennett@berkeley.edu
#   March 2021
#
###############################################################################

#' Inheritance Cube: Split-Drive for Andrea Smidler
#'
#' This is a modified split-drive construct for split-suppression drives.
#' This is an autosomal, 2-locus drive. The first locus contains the Cas9 allele
#' and the second locus has the gRNA construct. Cas9 efficacy in each sex is
#' dependent upon parent of inheritance.
#' This construct has 3 alleles at the first locus and 4 alleles at the second.
#'  * Locus 1
#'    * W: Wild-type
#'    * P: Paternal Cas9
#'    * M: Maternal Cas9
#'  * Locus 2
#'    * W: Wild-type
#'    * G: gRNAs
#'    * R: Resistant 1
#'    * B: Resistant 2
#'
#'
#' @param cMM Cutting efficacy of maternally-inherited Cas9 in males
#' @param chMM Homing efficacy of maternally-inherited Cas9 in males
#' @param crMM Resistance rate of maternally-inherited Cas9 in males
#'
#' @param cPM Cutting efficacy of paternally-inherited Cas9 in males
#' @param chPM Homing efficacy of paternally-inherited Cas9 in males
#' @param crPM Resistance rate of paternally-inherited Cas9 in males
#'
#' @param cMF Cutting efficacy of maternally-inherited Cas9 in females
#' @param chMF Homing efficacy of maternally-inherited Cas9 in females
#' @param crMF Resistance rate of maternally-inherited Cas9 in females
#'
#' @param cPF Cutting efficacy of paternally-inherited Cas9 in females
#' @param chPF Homing efficacy of paternally-inherited Cas9 in females
#' @param crPF Resistance rate of paternally-inherited Cas9 in females
#'
#' @param crF Female crossover rate, 0 is completely linked, 0.5 is unlinked, 1.0 is complete divergence
#' @param crM Male crossover rate, 0 is completely linked, 0.5 is unlinked, 1.0 is complete divergence
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
cubeASmidler <- function(cMM=0, chMM=0, crMM=0,
                         cPM=0, chPM=0, crPM=0,
                         cMF=0, chMF=0, crMF=0,
                         cPF=0, chPF=0, crPF=0,
                         crF=0.5, crM=0.5,
                         eta=NULL, phi=NULL,omega=NULL, xiF=NULL, xiM=NULL, s=NULL){


  # This is heavily influenced by the tGD cube.
  # It's basically stolen straight from that, with some reductions in locus 1
  #  stuff and increase in male/female based parameters.
  # Deposition was removed because females with both pieces are sterile in
  # the current experiments - this may not always be true.
  #  Also, see other cubes with notes about deposition being wrong via this design.
  # There is no known copy-number-dependent action, so it has been left out.


  # # for testing
  # testVec <- runif(n = 14)
  # cMM <- testVec[1]
  # chMM <- testVec[2]
  # crMM <- testVec[3]
  # cPM <- testVec[4]
  # chPM <- testVec[5]
  # crPM <- testVec[6]
  # cMF <- testVec[7]
  # chMF <- testVec[8]
  # crMF <- testVec[9]
  # cPF <- testVec[10]
  # chPF <- testVec[11]
  # crPF <- testVec[12]
  # crF <- testVec[13]
  # crM <- testVec[14]

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
  params <- c(cMM,chMM,crMM,cPM,chPM,crPM,cMF,chMF,crMF,cPF,chPF,crPF,crF,crM)
  if(any(params>1) || any(params<0)){
    stop("Parameters are rates.\n0 <= x <= 1")
  }


  #############################################################################
  ## generate all genotypes, set up vectors and matrices
  #############################################################################

  #list of possible alleles at each locus
  # the first locus has the Cas9
  # the second locus has the gRNAs (and assumed effector molecule)
  gTypes <- list(c('W', 'P', 'M'), c('W', 'G', 'R', 'B'))

  # # generate alleles
  # # this generates all combinations of Target Sites 1 and 2 for one allele,
  # #  then mixes alleles, so each half is a different allele.
  # # expand combinations of target site 1 and 2 on one allele
  # hold <- expand.grid(gTypes,KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  # # paste them all together
  # hold <- do.call(what = paste0, args = list(hold[,1], hold[,2]))
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
  # genotypes <- grep(pattern = "P.P.|M.M.", x = genotypes, value = TRUE, invert = TRUE)


  genotypes <- c('WWWW','PWWW','MWWW','WGWW','PGWW','MGWW','WRWW','PRWW','MRWW','WBWW','PBWW',
                  'MBWW','MWPW','PWWG','MGPW','PWWR','MRPW','PWWB','MBPW','MWWG','MWPG','MWWR',
                  'MWPR','MWWB','MWPB','WGWG','PGWG','MGWG','WGWR','PRWG','MRWG','WBWG','PBWG',
                  'MBWG','MGPG','PGWR','MRPG','PGWB','MBPG','MGWR','MGPR','MGWB','MGPB','WRWR',
                  'PRWR','MRWR','WBWR','PBWR','MBWR','MRPR','PRWB','MBPR','MRWB','MRPB','WBWB',
                  'PBWB','MBWB','MBPB')


  #############################################################################
  ## setup all probability lists
  #############################################################################

  # Female Mendelian
  fMendelian <- list("W"=c("W"=1),
                     "P"=c("M"=1),
                     "M"=c("M"=1),
                     "G"=c("G"=1),
                     "R"=c("R"=1),
                     "B"=c("B"=1))

  # Male Mendelian
  mMendelian <- list("W"=c("W"=1),
                     "P"=c("P"=1),
                     "M"=c("P"=1),
                     "G"=c("G"=1),
                     "R"=c("R"=1),
                     "B"=c("B"=1))

  # Female Locus 2
  fLocus2 <- list()
  fLocus2$Maternal <- list("W"=c("W"=1-cMF,
                                "G"=cMF*chMF,
                                "R"=cMF*(1-chMF)*crMF,
                                "B"=cMF*(1-chMF)*(1-crMF)),
                          "G"=c("G"=1),
                          "R"=c("R"=1),
                          "B"=c("B"=1))

  fLocus2$Paternal <- list("W"=c("W"=1-cPF,
                                "G"=cPF*chPF,
                                "R"=cPF*(1-chPF)*crPF,
                                "B"=cPF*(1-chPF)*(1-crPF)),
                          "G"=c("G"=1),
                          "R"=c("R"=1),
                          "B"=c("B"=1))

  # Male Locus 2
  mLocus2 <- list()
  mLocus2$Maternal <- list("W"=c("W"=1-cMM,
                                "G"=cMM*chMM,
                                "R"=cMM*(1-chMM)*crMM,
                                "B"=cMM*(1-chMM)*(1-crMM)),
                          "G"=c("G"=1),
                          "R"=c("R"=1),
                          "B"=c("B"=1))

  mLocus2$Paternal <- list("W"=c("W"=1-cPM,
                                "G"=cPM*chPM,
                                "R"=cPM*(1-chPM)*crPM,
                                "B"=cPM*(1-chPM)*(1-crPM)),
                          "G"=c("G"=1),
                          "R"=c("R"=1),
                          "B"=c("B"=1))


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
  numAlleles <- length(gTypes)
  fScore <- mScore <- logical(length = numAlleles)

  #############################################################################
  ## loop over all matings, female outer loop
  #############################################################################
  for(fi in 1:numGen){
    #do female stuff here
    #This splits all characters.
    fSplit <- strsplit(x = genotypes[fi], split = "", useBytes = TRUE)[[1]]
    #Matrix of alleles and target sites
    #  Each row is an allele, there are 2 alleles
    #  Each column is the target site on that allele, there are 2 target sites
    #  this is because target sites are linked on an allele
    fAlleleMat <- matrix(data = fSplit, nrow = numAlleles, ncol = numAlleles, byrow = TRUE)
    #Score them
    fScore[1] <- any("P" == fAlleleMat) || any("M" == fAlleleMat)
    fScore[2] <- any("G" == fAlleleMat)

    #setup offspring allele lists
    fPHold <- rep(x = list(list()), times = numAlleles) #vector(mode = "list", length = numAlleles)
    fAllele <- character(length = 0L)
    fProbs <- numeric(length = 0L)


    # check if there is a Cas9 from either parent
    if(fScore[1] && fScore[2]){
      # There is a Cas9 from one parent or the other
      # There are gRNAs
      # Thus, there is homing

      # loop over alleles, must keep target sites linked
      for(allele in 1:numAlleles){
        # we make the assumption that female homing is more important than male
        # ie, if we have M and P at locus one on both alleles, the M takes precedence
        #  in terms of homing.

        # check if maternal homing
        if(any(fAlleleMat[ ,1]=="M")){
          # at least one of the Cas9 proteins is from the mother
          fPHold[[allele]][[1]] <- fMendelian[[ fAlleleMat[allele,1] ]]
          fPHold[[allele]][[2]] <- fLocus2$Maternal[[ fAlleleMat[allele,2] ]]
        } else {
          # the Cas9 protein is from the father
          fPHold[[allele]][[1]] <- fMendelian[[ fAlleleMat[allele,1] ]]
          fPHold[[allele]][[2]] <- fLocus2$Paternal[[ fAlleleMat[allele,2] ]]
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

      } # end loop over alleles
    } # end female scoring


    # perform cross-overs
    hold1 <- fPHold[[1]][[1]] # need

    fPHold[[1]][[1]] <- c((1-crF)*hold1, crF*fPHold[[2]][[1]])
    fPHold[[2]][[1]] <- c((1-crF)*fPHold[[2]][[1]], crF*hold1)

    # all combinations of female alleles.
    for(allele in 1:numAlleles){
      # make combinations of the allele, then store those combinations to mix
      #  with the next allele
      # expand combinations
      holdProbs <- expand.grid(fPHold[[allele]],KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      holdAllele <- expand.grid(names(fPHold[[allele]][[1]]), names(fPHold[[allele]][[2]]),
                                KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

      # contract (paste or multiply) combinations and store
      fProbs <- c(fProbs, holdProbs[,1]*holdProbs[,2])
      fAllele <- c(fAllele, file.path(holdAllele[,1], holdAllele[,2], fsep = ""))
    }

    # remove zeros
    fAllele <- fAllele[fProbs!=0]
    fProbs <- fProbs[fProbs!=0]


    ###########################################################################
    ## loop over male mate. This is the inner loop
    ###########################################################################
    for(mi in 1:numGen){
      # isn't symmetric

      #do male stuff here
      #split male genotype
      #This splits all characters.
      mSplit <- strsplit(x = genotypes[mi], split = "", useBytes = TRUE)[[1]]
      #Matrix of alleles and target sites
      #  Each row is an allele, there are 2 alleles
      #  Each column is the target site on that allele, there are 2 target sites
      #  this is because target sites are linked on an allele
      mAlleleMat <- matrix(data = mSplit, nrow = numAlleles, ncol = numAlleles, byrow = TRUE)
      #Score them
      mScore[1] <- any("P" == mAlleleMat) || any("M" == mAlleleMat)
      mScore[2] <- any("G" == mAlleleMat)

      #setup offspring allele lists
      mPHold <- rep(x = list(list()), times = numAlleles) #vector(mode = "list", length = numAlleles)
      mAllele <- character(length = 0L)
      mProbs <- numeric(length = 0L)


      # check if there is a Cas9 from either parent
      if(mScore[1] && mScore[2]){
        # There is a Cas9 from one parent or the other
        # There are gRNAs
        # Thus, there is homing

        # loop over alleles, must keep target sites linked
        for(allele in 1:numAlleles){
          # we make the assumption that female homing is more important than male
          # ie, if we have M and P at locus one on both alleles, the M takes precedence
          #  in terms of homing.

          # check if maternal homing
          if(any(mAlleleMat[ ,1]=="M")){
            # at least one of the Cas9 proteins is from the mother
            mPHold[[allele]][[1]] <- mMendelian[[ mAlleleMat[allele,1] ]]
            mPHold[[allele]][[2]] <- mLocus2$Maternal[[ mAlleleMat[allele,2] ]]
          } else {
            # the Cas9 protein is from the father
            mPHold[[allele]][[1]] <- mMendelian[[ mAlleleMat[allele,1] ]]
            mPHold[[allele]][[2]] <- mLocus2$Paternal[[ mAlleleMat[allele,2] ]]
          } # end M/P homing

        } # end loop over alleles
      } else {
        # either no Cas9 or no gRNAs
        # all inheritance is Mendelian

        # loop over alleles, must keep target sites linked
        for(allele in 1:numAlleles){
          # set target site 1, then target site 2
          mPHold[[allele]][[1]] <- mMendelian[[ mAlleleMat[allele,1] ]]
          mPHold[[allele]][[2]] <- mMendelian[[ mAlleleMat[allele,2] ]]

        } # end loop over alleles
      } # end male scoring


      # perform cross-overs
      hold1 <- mPHold[[1]][[1]] # need

      mPHold[[1]][[1]] <- c((1-crM)*hold1, crM*mPHold[[2]][[1]])
      mPHold[[2]][[1]] <- c((1-crM)*mPHold[[2]][[1]], crM*hold1)

      # all combinations of female alleles.
      for(allele in 1:numAlleles){
        # make combinations of the allele, then store those combinations to mix
        #  with the next allele
        # expand combinations
        holdProbs <- expand.grid(mPHold[[allele]],KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        holdAllele <- expand.grid(names(mPHold[[allele]][[1]]), names(mPHold[[allele]][[2]]),
                                  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

        # contract (paste or multiply) combinations and store
        mProbs <- c(mProbs, holdProbs[,1]*holdProbs[,2])
        mAllele <- c(mAllele, file.path(holdAllele[,1], holdAllele[,2], fsep = ""))
      }

      # remove zeros for test
      mAllele <- mAllele[mProbs!=0]
      mProbs <- mProbs[mProbs!=0]


      #########################################################################
      ## Get combinations and put them in the tMatrix. This must be done
      ##  inside the inner loop
      #########################################################################

      # male and female alleles/probs are allready combined by target sites, and
      #  we assume alleles segregate indepenently, so we just have to get combinations
      #  of male vs female for the offspring
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
  modifiers = cubeModifiers(genotypes, eta = eta, phi = phi, omega = omega, xiF = xiF, xiM = xiM, s = s)

  ## put everything into a labeled list to return
  return(list(
    ih = tMatrix,
    tau = viabilityMask,
    genotypesID = genotypes,
    genotypesN = numGen,
    wildType = "WWWW",
    eta = modifiers$eta,
    phi = modifiers$phi,
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = "MGPG"
  ))

}
