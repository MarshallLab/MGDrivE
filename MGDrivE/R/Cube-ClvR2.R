###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDD: Mosquito Gene Drive Explorer - Discrete
#   2-Locus ClvR (Cleave and Rescue)
#    - sex specific homing
#    - copy-number dependent female deposition
#      - only correct for broken alleles, no proper homing from deposition*
#    - crossover
#   August 2020
#   jared_bennett@berkeley.edu
#
###############################################################################

#' Inheritance Cube: 2-Locus ClvR (Cleave and Rescue)
#'
#' Based on the Cleave-and-Rescue system of [Oberhofer](https://doi.org/10.1101/2020.07.09.196253),
#' this is a 3-locus Cas9-based toxin-antidote system. The first locus carries the
#' Cas9, the second locus carries the gRNAs, and a recoded copy of an essential gene.
#' The third locus is the targeted essential gene. This gene can be completely
#' haplosufficient (\code{hSuf} = 1) or completely haploinsufficient (\code{hSuf} = 0).
#' It is assumed that having 2 copies of the gene (be it wild-type at the second
#' locus or recoded at the first) confers complete viability. Additionally, loci
#' 1 and 2 can be linked, given \code{crM} and \code{crF}, imitating the original
#' 2-locus ClvR system.
#' For this construct, the first locus will have 2 alleles, the second will have 2
#' alleles, and the third will have 3 alleles:
#'  * Locus 1
#'    * W: Wild-type
#'    * C: Cas9
#'  * Locus 2
#'    * W: Wild-type
#'    * G: gRNAs and recoded essential gene
#'  * Locus 3
#'    * W: Wild-type
#'    * R: Functional resistant
#'    * B: Non-functional resistant
#'
#' Female deposition is implemented incorrectly. Right now, it is performed on
#' male alleles prior to zygote formation - it should happen post-zygote formation.
#' Since this construct doesn't have HDR, this should be fine. \cr
#' Additionally, it is assumed that deposition requries loaded Cas9-RNP complexes
#' from the mother, having Cas9 and no maternal gRNA, even in the presence of
#' paternal gRNA, will not result in maternal deposition mediated cleavage.
#'
#' Copy-number dependent rates are based on Cas9, not gRNA. The assumption is that
#' RNA is easier to produce, and therefore won't limit cleavage by Cas9.
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
#'
#' @param crossF Crossover rate in females, 0 is completely linked, 0.5 is unlinked, 1.0 is complete divergence
#' @param crossM Crossover rate in males, 0 is completely linked, 0.5 is unlinked, 1.0 is complete divergence
#'
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
cubeClvR2 <- function(cF = 1, crF = 0, ccF = cF, ccrF = crF,
                      cM = 1, crM = 0, ccM = cM, ccrM = crM,
                      dW = 0, drW = 0, ddW = dW, ddrW = drW,
                      hSuf = 1, crossF = 0.5, crossM = 0.5,
                      eta=NULL, phi=NULL,omega=NULL, xiF=NULL, xiM=NULL, s=NULL){

  # # for testing
  # testVec <- runif(n = 15)
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
  # crossF <- testVec[13]; crossM <- testVec[14];
  #
  # hSuf <- testVec[15];
  #
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
  params <- c(cF,crF,ccF,ccrF,cM,crM,ccM,ccrM,dW,drW,ddW,ddrW,hSuf,crossF,crossM)
  if(any(params>1) || any(params<0)){
    stop("Parameters are rates, they must be between 0 and 1.")
  }


  #############################################################################
  ## generate all genotypes, set up vectors and matrices
  #############################################################################

  #list of possible alleles at each locus
  # the first locus has the Cas9
  # the second locus has the gRNAs (and assumed effector molecule)
  gTypes <- list(c('W', 'C'), c('W', 'G'))
  gTypes2 <- c('W', 'R', 'B')

  # # generate alleles
  # # this generates all combinations of Target Sites 1 and 2 for one allele,
  # #  then mixes alleles, so each half is a different allele.
  # # expand combinations of target site 1 and 2 on one allele
  # hold <- expand.grid(gTypes,KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  # hold1 <- expand.grid(gTypes2,gTypes2,KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  # # paste them all together
  # hold <- do.call(what = paste0, args = list(hold[,1], hold[,2]))
  # #get all combinations of both loci
  # openAlleles <- expand.grid(hold, hold, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  # # sort both alleles to remove uniques
  # # ie, each target site is unique and can't be moved
  # # but, the alleles aren't unique, and can occur in any order
  # openAlleles <- apply(X = openAlleles, MARGIN = 1, FUN = sort)
  # openAlleles1 <- apply(X = hold1, MARGIN = 1, FUN = sort)
  # # paste together and get unique ones only
  # genotypes <- unique(do.call(what = paste0, args = list(openAlleles[1, ], openAlleles[2, ])))
  # genotypes1 <- unique(do.call(what = paste0, args = list(openAlleles1[1, ], openAlleles1[2, ])))
  # # Finally, combine 1/2 and 3 loci
  # genTot <- expand.grid(genotypes1, genotypes, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  # genTot <- do.call(what = paste0, args = list(genTot[ ,2], genTot[ ,1]))

  genotypes <- c('WWWWWW','WWWWRW','WWWWBW','WWWWRR','WWWWBR','WWWWBB',
                 'CWWWWW','CWWWRW','CWWWBW','CWWWRR','CWWWBR','CWWWBB',
                 'WGWWWW','WGWWRW','WGWWBW','WGWWRR','WGWWBR','WGWWBB',
                 'CGWWWW','CGWWRW','CGWWBW','CGWWRR','CGWWBR','CGWWBB',
                 'CWCWWW','CWCWRW','CWCWBW','CWCWRR','CWCWBR','CWCWBB',
                 'CWWGWW','CWWGRW','CWWGBW','CWWGRR','CWWGBR','CWWGBB',
                 'CGCWWW','CGCWRW','CGCWBW','CGCWRR','CGCWBR','CGCWBB',
                 'WGWGWW','WGWGRW','WGWGBW','WGWGRR','WGWGBR','WGWGBB',
                 'CGWGWW','CGWGRW','CGWGBW','CGWGRR','CGWGBR','CGWGBB',
                 'CGCGWW','CGCGRW','CGCGBW','CGCGRR','CGCGBR','CGCGBB')


  #############################################################################
  ## setup all probability lists
  #############################################################################

  # Locus 1/2
  locus12 <- list("W"=c("W"=1),
                  "C"=c("C"=1),
                  "G"=c("G"=1))


  # Female Locus 3
  fLocus3 <- list()
  fLocus3$mendelian <- list("W"=c("W"=1),
                            "R"=c("R"=1),
                            "B"=c("B"=1))

  fLocus3$oneCas <- list("W"=c("W"=1-cF,
                               "R"=cF*crF,
                               "B"=cF*(1-crF)),
                          "R"=c("R"=1),
                          "B"=c("B"=1))

  fLocus3$twoCas <- list("W"=c("W"=1-ccF,
                               "R"=ccF*ccrF,
                               "B"=ccF*(1-ccrF)),
                            "R"=c("R"=1),
                            "B"=c("B"=1))


  # Male Locus 3
  mLocus3 <- list()

  mLocus3$zero <- list("W"=c("W"=1),
                       "R"=c("R"=1),
                       "B"=c("B"=1))
  mLocus3$zeroOne <- list("W"=c("W"=1-dW,
                                "R"=dW*drW,
                                "B"=dW*(1-drW)),
                          "R"=c("R"=1),
                          "B"=c("B"=1))
  mLocus3$zeroTwo <- list("W"=c("W"=1-ddW,
                                "R"=ddW*ddrW,
                                "B"=ddW*(1-ddrW)),
                          "R"=c("R"=1),
                          "B"=c("B"=1))

  mLocus3$one <- list("W"=c("W"=1-cM,
                            "R"=cM*crM,
                            "B"=cM*(1-crM)),
                      "R"=c("R"=1),
                      "B"=c("B"=1))
  mLocus3$oneOne <- list("W"=c("W"=(1-cM)*(1-dW),
                                "R"=cM*crM + dW*drW,
                                "B"=cM*(1-crM) + dW*(1-drW)),
                         "R"=c("R"=1),
                         "B"=c("B"=1))
  mLocus3$oneTwo <- list("W"=c("W"=(1-cM)*(1-ddW),
                               "R"=cM*crM + ddW*ddrW,
                               "B"=cM*(1-crM) + ddW*(1-ddrW)),
                         "R"=c("R"=1),
                         "B"=c("B"=1))

  mLocus3$two <- list("W"=c("W"=1-ccM,
                            "R"=ccM*ccrM,
                            "B"=ccM*(1-ccrM)),
                      "R"=c("R"=1),
                      "B"=c("B"=1))
  mLocus3$twoOne <- list("W"=c("W"=(1-ccM)*(1-dW),
                               "R"=ccM*ccrM + dW*drW,
                               "B"=ccM*(1-ccrM) + dW*(1-drW)),
                         "R"=c("R"=1),
                         "B"=c("B"=1))
  mLocus3$twoTwo <- list("W"=c("W"=(1-ccM)*(1-ddW),
                               "R"=ccM*ccrM + ddW*ddrW,
                               "B"=ccM*(1-ccrM) + ddW*(1-ddrW)),
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
  numAlleles <- 2
  numLoci <- 3
  fScore <- mScore <- integer(length = numAlleles)


  #############################################################################
  ## loop over all matings, female outer loop
  #############################################################################
  for(fi in 1:numGen){
    #do female stuff here
    #This splits all characters.
    fSplit <- strsplit(x = genotypes[fi], split = "", useBytes = TRUE)[[1]]
    #Matrix of alleles and target sites
    #  Each row is an allele, there are 2 alleles
    #  Each column is the target locus on that allele, there are 3 target sites
    #  this is because target sites are linked on an allele, for the first 2
    #  the third one is independent, and it's just weird here.
    #    Locus1, Locus1, Locus3
    #    Locus2, Locus2, Locus3
    fAlleleMat <- matrix(data = fSplit, nrow = numAlleles, ncol = numLoci)
    #Score them
    fScore[1] <- sum('C' == fAlleleMat)
    fScore[2] <- any('G' == fAlleleMat)

    #setup offspring allele lists
    fPHold <- rep(x = list(list()), times = numLoci)
    fAllele <- character(length = 0L) # these 2 things need cleared each loop, hence this line
    fProbs <- numeric(length = 0L)


    # loci 1/2 are always mendelian
    fPHold[[1]] <- locus12[ fAlleleMat[ ,1] ]
    fPHold[[2]] <- locus12[ fAlleleMat[ ,2] ]

    # check if there is Cas9 and gRNAs for locus 3
    if(all(fScore)){
      # there is cas9 and gRNA

      # how many Cas9
      if(fScore[1] == 1){
        # 1 cas9
        fPHold[[3]] <- fLocus3$oneCas[ fAlleleMat[ ,3] ]

      } else {
        # 2 cas9
        fPHold[[3]] <- fLocus3$twoCas[ fAlleleMat[ ,3] ]

      } # end Cas9 number

    } else {
      # either cas9 or gRNAs are missing
      # mendelian inheritance only
      fPHold[[3]] <- fLocus3$mendelian[ fAlleleMat[ ,3] ]

    } # end female alleles


    # perform cross-overs
    hold1 <- fPHold[[1]][[1]] # need

    fPHold[[1]][[1]] <- c((1-crossF)*hold1, crossF*fPHold[[2]][[1]])
    fPHold[[2]][[1]] <- c((1-crossF)*fPHold[[2]][[1]], crossF*hold1)

    # all combinations of female alleles at loci 1/2
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
      #do male stuff here
      #split male genotype
      #This splits all characters.
      mSplit <- strsplit(x = genotypes[mi], split = "", useBytes = TRUE)[[1]]
      #Matrix of alleles and target sites
      #  Each row is an allele, there are 2 alleles
      #  Each column is the target locus on that allele, there are 3 target sites
      #  this is because target sites are linked on an allele, for the first 2
      #  the third one is independent, and it's just weird here.
      #    Locus1, Locus1, Locus3
      #    Locus2, Locus2, Locus3
      mAlleleMat <- matrix(data = mSplit, nrow = numAlleles, ncol = numLoci)
      #Score them
      mScore[1] <- sum('C' == mAlleleMat)
      mScore[2] <- any('G' == mAlleleMat)

      #setup offspring allele lists
      mPHold <- rep(x = list(list()), times = numLoci)
      mAllele <- character(length = 0L)
      mProbs <- numeric(length = 0L)


      # loci 1/2 are always mendelian
      mPHold[[1]] <- locus12[ mAlleleMat[ ,1] ]
      mPHold[[2]] <- locus12[ mAlleleMat[ ,2] ]


      # check if there is Cas9 and gRNAs for locus 3
      if(all(mScore)){
        # there is cas9 and gRNA

        # how many Cas9
        if(mScore[1] == 1){
          # 1 cas9 for male homing

          # is there female deposition
          if(all(fScore)){
            # there is female deposition

            # check type of deposition
            if(fScore[1]== 1){
              # 1 allele deposition
              mPHold[[3]] <- mLocus3$oneOne[ mAlleleMat[ ,3] ]

            } else {
              # 2 allele deposition
              mPHold[[3]] <- mLocus3$oneTwo[ mAlleleMat[ ,3] ]

            } # end deposition check

          } else {
            # no female deposition
            # just male homing
            mPHold[[3]] <- mLocus3$one[ mAlleleMat[ ,3] ]

          } # 1 allele male homing

        } else {
          # 2 cas9 for male homing

          # is there female deposition
          if(all(fScore)){
            # there is female deposition

            # check type of deposition
            if(fScore[1]== 1){
              # 1 allele deposition
              mPHold[[3]] <- mLocus3$twoOne[ mAlleleMat[ ,3] ]

            } else {
              # 2 allele deposition
              mPHold[[3]] <- mLocus3$twoTwo[ mAlleleMat[ ,3] ]

            } # end deposition check

          } else {
            # no female deposition
            # just male homing
            mPHold[[3]] <- mLocus3$two[ mAlleleMat[ ,3] ]

          } # 2 allele male homing

        } # end Cas9 number

      } else {
        # either cas9 or gRNAs are missing
        # mendelian inheritance only

        # is there female deposition
        if(all(fScore)){
          # there is female deposition

          # check type of deposition
          if(fScore[1]== 1){
            # 1 allele deposition
            mPHold[[3]] <- mLocus3$zeroOne[ mAlleleMat[ ,3] ]

          } else {
            # 2 allele deposition
            mPHold[[3]] <- mLocus3$zeroTwo[ mAlleleMat[ ,3] ]

          } # end deposition check

        } else {
          # no female deposition
          # no male homing
          mPHold[[3]] <- mLocus3$zero[ mAlleleMat[ ,3] ]

        } # male mendelian

      } # end male alleles


      # perform cross-overs
      hold1 <- mPHold[[1]][[1]] # need

      mPHold[[1]][[1]] <- c((1-crossM)*hold1, crossM*mPHold[[2]][[1]])
      mPHold[[2]][[1]] <- c((1-crossM)*mPHold[[2]][[1]], crossM*hold1)

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
      ##########
      ## Locus 1/2
      ##########
      # male and female alleles/probs are already combined by target sites, and
      #  we assume alleles segregate independently, so we just have to get combinations
      #  of male vs female for the offspring
      holdProbs <- expand.grid(fProbs, mProbs, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      holdAllele <- expand.grid(fAllele, mAllele, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

      # sort alleles into order
      holdAllele <- apply(X = holdAllele, MARGIN = 1, FUN = sort.int)

      # contract (paste or multiply) combinations
      holdProbs <- holdProbs[ ,1]*holdProbs[ ,2]
      holdAllele <- file.path(holdAllele[1,], holdAllele[2,], fsep = "")

      #aggregate duplicate genotypes
      reducedOneTwo <- vapply(X = unique(holdAllele), FUN = function(x){
        sum(holdProbs[holdAllele==x])},
        FUN.VALUE = numeric(length = 1L))


      ##########
      ## Locus 3
      ##########
      # get combinationes of 3rd locus, which is independent of the other 1/2
      # this is hideous, I'm sorry.
      holdProbs <- expand.grid(unlist(fPHold[[3]], use.names = FALSE), unlist(mPHold[[3]], use.names = FALSE),
                               KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      holdAllele <- expand.grid(names(c(fPHold[[3]][[1]], fPHold[[3]][[2]])),
                                names(c(mPHold[[3]][[1]], mPHold[[3]][[2]])),
                                KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

      # sort alleles into order
      holdAllele <- apply(X = holdAllele, MARGIN = 1, FUN = sort.int)

      # contract (paste or multiply) combinations
      holdProbs <- holdProbs[ ,1]*holdProbs[ ,2]
      holdAllele <- file.path(holdAllele[1,], holdAllele[2,], fsep = "")

      #aggregate duplicate genotypes
      reducedThree <- vapply(X = unique(holdAllele), FUN = function(x){
        sum(holdProbs[holdAllele==x])},
        FUN.VALUE = numeric(length = 1L))


      ##########
      ## All Loci
      ##########
      # here, combine the alleles at loci 1/2 and 3
      # all combinations of locus 1/2 and 3
      alleles <- as.vector(outer(X = names(reducedOneTwo), Y = names(reducedThree), FUN = file.path, fsep = ""))
      probs <- as.vector(outer(X = reducedOneTwo, Y = reducedThree, FUN = '*'))

      # normalize probs
      probs <- probs/sum(probs)

      #set values in tMatrix
      tMatrix[fi,mi, alleles ] <- probs

    }# end male loop
  }# end female loop


  ## set the other half of the matrix
  tMatrix[tMatrix < .Machine$double.eps] <- 0 #protection from underflow errors

  ## initialize viability mask.
  # Anything with 2 copies of the gene is fine, be it wild-type or recoded
  #   GG at second locus all survive, WW/WR/RR at third locus all survive
  # Anything with 0 copies never survives
  #   WWWWBB is always dead
  # hSuf defines heterozygous sufficiency, ie is it dominant
  # aperm allows us to fill the array by depth, which isn't normally possible
  #   resize keeps dim names
  viabilityVec <- c(1,1,hSuf,1,hSuf,0, 1,1,hSuf,1,hSuf,0, 1,1,1,1,1,hSuf,
                    1,1,1,1,1,hSuf, 1,1,hSuf,1,hSuf,0, 1,1,1,1,1,hSuf,
                    1,1,1,1,1,hSuf, 1,1,1,1,1,1, 1,1,1,1,1,1, 1,1,1,1,1,1)
  viabilityMask <- aperm(a = array(data = viabilityVec, dim = c(numGen,numGen,numGen),
                                   dimnames = list(genotypes, genotypes, genotypes)),
                         perm = c(3,2,1),resize = TRUE)

  ## genotype-specific modifiers
  modifiers = cubeModifiers(genotypes, eta = eta, phi = phi, omega = omega, xiF = xiF, xiM = xiM, s = s)

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
    releaseType = "CGCGWW"
  ))

}
