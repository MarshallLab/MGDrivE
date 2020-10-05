###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   tGD (trans-complementing gene drive )
#   Héctor Sanchez, Jared Bennett, Sean Wu, John Marshall
#   jared_bennett@berkeley.edu
#   February 2018
#
###############################################################################

#' Inheritance Cube: tGD
#'
#' The trans-complementing Gene Drive (tGD) is a 1-locus, 2 target site drive.
#' The first target site corresponds to the Cas protein, the second to an effector
#' gene and the gRNAs. There are two sets of gRNAs, because each target site may have
#' different cutting/homing/resistance rates, and each sex can have different rates
#' for all of those things. Additionally, the parent that you receive your Cas from
#' dictates its efficiency.
#' Therefor, this construct has 5 alleles at the first locus and 4 alleles at the second.
#'  * Locus 1
#'    * W: Wild-type
#'    * P: Paternal Cas9
#'    * M: Maternal Cas9
#'    * R: Resistant allele 1
#'    * B: Resistant allele 2
#'  * Locus 2
#'    * W: Wild-type
#'    * G: gRNAs
#'    * R: Resistant 1
#'    * B: Resistant 2
#'
#' This drive corresponds to the [transcomplementing gene drive](https://www.nature.com/articles/s41467-019-13977-7)
#' developed by the Gantz and Bier lab.
#'
#' @param cM1 Maternally inherited Cas9 cutting rate at locus 1
#' @param cM2 Maternally inherited Cas9 cutting rate at locus 2
#' @param cP1 Paternally inherited Cas9 cutting rate at locus 1
#' @param cP2 Paternally inherited Cas9 cutting rate at locus 2
#' @param hM1 Maternally inherited Cas9 homing efficiency at locus 1
#' @param hM2 Maternally inherited Cas9 homing efficiency at locus 2
#' @param hP1 Paternally inherited Cas9 homing efficiency at locus 1
#' @param hP2 Paternally inherited Cas9 homing efficiency at locus 2
#' @param rM1 Maternally inherited Cas9 resistance efficiency at locus 1
#' @param rM2 Maternally inherited Cas9 resistance efficiency at locus 2
#' @param rP1 Paternally inherited Cas9 resistance efficiency at locus 1
#' @param rP2 Paternally inherited Cas9 resistance efficiency at locus 2
#' @param crM Maternal crossover rate, 0 is completely linked, 0.5 is unlinked, 1.0 is complete divergence
#' @param crP Paternal crossover rate, 0 is completely linked, 0.5 is unlinked, 1.0 is complete divergence
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
cubeTGD <- function(cM1=0, cM2=0, cP1=0, cP2=0,
                     hM1=0, hM2=0, hP1=0, hP2=0,
                     rM1=0, rM2=0, rP1=0, rP2=0, crM=0, crP=0,
                     eta=NULL, phi=NULL,omega=NULL, xiF=NULL, xiM=NULL, s=NULL){

  ## safety checks
  if(any(c(cM1,cM2,cP1,cP2,hM1,hM2,hP1,hP2,rM1,rM2,rP1,rP2,crM,crP)>1) || any(c(cM1,cM2,cP1,cP2,hM1,hM2,hP1,hP2,rM1,rM2,rP1,rP2,crM,crP)<0)){
    stop("Parameters are rates.\n0 <= x <= 1")
  }


  #############################################################################
  ## generate all genotypes, set up vectors and matrices
  #############################################################################

  #list of possible alleles at each locus
  # the first locus has the Cas9
  # the second locus has the gRNAs (and assumed effector molecule)
  gTypes <- list(c('W', 'P', 'M', 'R', 'B'), c('W', 'G', 'R', 'B'))

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

  genotypes <- c("WWWW","PWWW","MWWW","RWWW","BWWW","WGWW","PGWW","MGWW","RGWW",
                 "BGWW","WRWW","PRWW","MRWW","RRWW","BRWW","WBWW","PBWW","MBWW",
                 "RBWW","BBWW","PWPW","MWPW","PWRW","BWPW","PWWG","PGPW","MGPW",
                 "PWRG","BGPW","PWWR","PRPW","MRPW","PWRR","BRPW","PWWB","PBPW",
                 "MBPW","PWRB","BBPW","MWMW","MWRW","BWMW","MWWG","MWPG","MGMW",
                 "MWRG","BGMW","MWWR","MWPR","MRMW","MWRR","BRMW","MWWB","MWPB",
                 "MBMW","MWRB","BBMW","RWRW","BWRW","RWWG","PGRW","MGRW","RGRW",
                 "BGRW","RWWR","PRRW","MRRW","RRRW","BRRW","RWWB","PBRW","MBRW",
                 "RBRW","BBRW","BWBW","BWWG","BWPG","BWMG","BWRG","BGBW","BWWR",
                 "BWPR","BWMR","BWRR","BRBW","BWWB","BWPB","BWMB","BWRB","BBBW",
                 "WGWG","PGWG","MGWG","RGWG","BGWG","WGWR","PRWG","MRWG","RRWG",
                 "BRWG","WBWG","PBWG","MBWG","RBWG","BBWG","PGPG","MGPG","PGRG",
                 "BGPG","PGWR","PGPR","MRPG","PGRR","BRPG","PGWB","PBPG","MBPG",
                 "PGRB","BBPG","MGMG","MGRG","BGMG","MGWR","MGPR","MGMR","MGRR",
                 "BRMG","MGWB","MGPB","MBMG","MGRB","BBMG","RGRG","BGRG","RGWR",
                 "PRRG","MRRG","RGRR","BRRG","RGWB","PBRG","MBRG","RBRG","BBRG",
                 "BGBG","BGWR","BGPR","BGMR","BGRR","BGBR","BGWB","BGPB","BGMB",
                 "BGRB","BBBG","WRWR","PRWR","MRWR","RRWR","BRWR","WBWR","PBWR",
                 "MBWR","RBWR","BBWR","PRPR","MRPR","PRRR","BRPR","PRWB","PBPR",
                 "MBPR","PRRB","BBPR","MRMR","MRRR","BRMR","MRWB","MRPB","MBMR",
                 "MRRB","BBMR","RRRR","BRRR","RRWB","PBRR","MBRR","RBRR","BBRR",
                 "BRBR","BRWB","BRPB","BRMB","BRRB","BBBR","WBWB","PBWB","MBWB",
                 "RBWB","BBWB","PBPB","MBPB","PBRB","BBPB","MBMB","MBRB","BBMB",
                 "RBRB","BBRB","BBBB")


  #############################################################################
  ## setup all probability lists
  #############################################################################

  # Female Locus 1
  fLocus1 <- list()
  fLocus1$mendelian <- list("W"=c("W"=1),
                            "P"=c("M"=1),
                            "M"=c("M"=1),
                            "R"=c("R"=1),
                            "B"=c("B"=1))

  fLocus1$paternal <- list("W"=c("W"=1-cP1,
                                 "M"=cP1*hP1,
                                 "R"=cP1*(1-hP1)*rP1,
                                 "B"=cP1*(1-hP1)*(1-rP1)),
                           "P"=c("M"=1),
                           "M"=c("M"=1),
                           "R"=c("R"=1),
                           "B"=c("B"=1))

  fLocus1$maternal <- list("W"=c("W"=1-cM1,
                                 "M"=cM1*hM1,
                                 "R"=cM1*(1-hM1)*rM1,
                                 "B"=cM1*(1-hM1)*(1-rM1)),
                           "P"=c("M"=1),
                           "M"=c("M"=1),
                           "R"=c("R"=1),
                           "B"=c("B"=1))

  # Male Locus 1
  mLocus1 <- list()
  mLocus1$mendelian <- list("W"=c("W"=1),
                            "P"=c("P"=1),
                            "M"=c("P"=1),
                            "R"=c("R"=1),
                            "B"=c("B"=1))

  mLocus1$paternal <- list("W"=c("W"=1-cP1,
                                 "P"=cP1*hP1,
                                 "R"=cP1*(1-hP1)*rP1,
                                 "B"=cP1*(1-hP1)*(1-rP1)),
                           "P"=c("P"=1),
                           "M"=c("P"=1),
                           "R"=c("R"=1),
                           "B"=c("B"=1))

  mLocus1$maternal <- list("W"=c("W"=1-cM1,
                                 "P"=cM1*hM1,
                                 "R"=cM1*(1-hM1)*rM1,
                                 "B"=cM1*(1-hM1)*(1-rM1)),
                           "P"=c("P"=1),
                           "M"=c("P"=1),
                           "R"=c("R"=1),
                           "B"=c("B"=1))


  # Locus 2
  Locus2 <- list()
  Locus2$mendelian <- list("W"=c("W"=1),
                           "G"=c("G"=1),
                           "R"=c("R"=1),
                           "B"=c("B"=1))

  Locus2$paternal <- list("W"=c("W"=1-cP2,
                                "G"=cP2*hP2,
                                "R"=cP2*(1-hP2)*rP2,
                                "B"=cP2*(1-hP2)*(1-rP2)),
                          "G"=c("G"=1),
                          "R"=c("R"=1),
                          "B"=c("B"=1))

  Locus2$maternal <- list("W"=c("W"=1-cM2,
                                "G"=cM2*hM2,
                                "R"=cM2*(1-hM2)*rM2,
                                "B"=cM2*(1-hM2)*(1-rM2)),
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
  for (fi in 1:numGen){
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
          fPHold[[allele]][[1]] <- fLocus1$maternal[[ fAlleleMat[allele,1] ]]
          fPHold[[allele]][[2]] <- Locus2$maternal[[ fAlleleMat[allele,2] ]]
        } else {
          # the Cas9 protein is from the father
          fPHold[[allele]][[1]] <- fLocus1$paternal[[ fAlleleMat[allele,1] ]]
          fPHold[[allele]][[2]] <- Locus2$paternal[[ fAlleleMat[allele,2] ]]
        } # end M/P homing

      } # end loop over alleles
    } else {
      # either no Cas9 or no gRNAs
      # all inheritance is Mendelian

      # loop over alleles, must keep target sites linked
      for(allele in 1:numAlleles){
        # set target site 1, then target site 2
        fPHold[[allele]][[1]] <- fLocus1$mendelian[[ fAlleleMat[allele,1] ]]
        fPHold[[allele]][[2]] <- Locus2$mendelian[[ fAlleleMat[allele,2] ]]

      } # end loop over alleles
    } # end female scoring


    # perform cross-overs
    hold1 <- fPHold[[1]][[1]] # need

    fPHold[[1]][[1]] <- c((1-crM)*hold1, crM*fPHold[[2]][[1]])
    fPHold[[2]][[1]] <- c((1-crM)*fPHold[[2]][[1]], crM*hold1)

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
    for (mi in 1:numGen){
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
            mPHold[[allele]][[1]] <- mLocus1$maternal[[ mAlleleMat[allele,1] ]]
            mPHold[[allele]][[2]] <- Locus2$maternal[[ mAlleleMat[allele,2] ]]
          } else {
            # the Cas9 protein is from the father
            mPHold[[allele]][[1]] <- mLocus1$paternal[[ mAlleleMat[allele,1] ]]
            mPHold[[allele]][[2]] <- Locus2$paternal[[ mAlleleMat[allele,2] ]]
          } # end M/P homing

        } # end loop over alleles
      } else {
        # either no Cas9 or no gRNAs
        # all inheritance is Mendelian

        # loop over alleles, must keep target sites linked
        for(allele in 1:numAlleles){
          # set target site 1, then target site 2
          mPHold[[allele]][[1]] <- mLocus1$mendelian[[ mAlleleMat[allele,1] ]]
          mPHold[[allele]][[2]] <- Locus2$mendelian[[ mAlleleMat[allele,2] ]]

        } # end loop over alleles
      } # end male scoring


      # perform cross-overs
      hold1 <- mPHold[[1]][[1]] # need

      mPHold[[1]][[1]] <- c((1-crP)*hold1, crP*mPHold[[2]][[1]])
      mPHold[[2]][[1]] <- c((1-crP)*mPHold[[2]][[1]], crP*hold1)

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
