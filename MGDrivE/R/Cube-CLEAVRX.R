###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   CLEAVR Cube - X-linked version
#   Héctor Sanchez, Jared Bennett, Sean Wu, John Marshall
#   February 2019
#   March 2019
#     Updated to reflect new desires from Ana
#
###############################################################################

#' Inheritance Cube: CLEAVR - Cleave and Rescue
#'
#' This is a novel cube from the Akbari lab. There are 2 loci: the first is a female
#' fertility locus (e.g. doubleSex), where the Cas, gRNAs, and a recoded essential
#' gene go. This locus is inherited in a Mendelian fashion, but is also targeted for
#' destruction by the homing allele. The second locus involves
#' an essential gene, for both males and females, and this is the target of the
#' gRNAs at the first locus. No homing is performed, it is simply destroyed. There
#' is different cutting rates in males and females, with no possibility for a rescuing
#' resistant allele. Females homozygous for the H or B alleles at locus 1 are viable but infertile,
#' while males are unaffected. All animals homozygous at locus two must contain the
#' recoded copy at locus 1 to be viable. This corresponds to the homing construct
#' being X-linked.
#'
#' @param cM1 Male cutting rate at first locus IS NOT USED IN THIS CUBE
#' @param cM2 Male cutting rate at second locus
#' @param cF1 Female cutting rate at first locus
#' @param cF2 Female cutting rate at second locus
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
cubeCLEAVRX <- function(cM1 = 1.0, cM2 = 1.0, cF1 = 1.0, cF2=1.0, eta = NULL,
                         phi = NULL, omega = NULL, xiF = NULL, xiM = NULL, s = NULL){

  ## safety checks
  if(any(c(cM1, cM2, cF1, cF2)>1) || any(c(cM1, cM2, cF1, cF2)<0)){
    stop("Parameters are rates, they must be between 0 and 1.")
  }

  # locus 1 has homing allele(H), wild-type(W), and borked (B) alleles, and Y since it is x-linked
  #  genotypes are WW, WH, WB, HH, HB, BB, WY, HY, BY
  #  inheritance is Mendelian, but the H allele promotes destruction of the W allele,
  #  ie, WH, with perfect homing, results in H and B gametes
  # Locus 2 has wild-type(W) and broken (B) alleles
  #  genotypes are WW, WB, BB
  #  inheritance is Mendelian, with cutting based on locus 1 genotype


  #############################################################################
  ## generate all genotypes, set up vectors and matrices
  #############################################################################

  #list of possible alleles at each locus
  gTypes <- list(c('W', 'H', 'B','Y'), c('W', 'B'))

  # # generate alleles
  # alleleList <- vector(mode = "list", length = 2)
  # # loop over both loci
  # for(i in 1:2){
  #   #get all combinations of alleles at locus 1
  #   hold <- expand.grid(gTypes[[i]], gTypes[[i]],KEEP.OUT.ATTRS = FALSE,
  #                       stringsAsFactors = FALSE)
  #   #sort them, paste them, keep the unique ones
  #   alleleList[[i]] <- unique(vapply(X = 1:dim(hold)[1],
  #                              FUN = function(x){
  #                                paste0(sort(x = hold[x, ]), collapse = "")},
  #                              FUN.VALUE = character(1)
  #                              )
  #                             )
  # }
  # # remove YY criter, not alive
  # alleleList[[1]] <- alleleList[[1]][-length(alleleList[[1]])]
  # #get all combinations of both loci
  # openAlleles <- expand.grid(alleleList, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  # #paste them together. This is the final list of genotypes
  # genotypes <- file.path(openAlleles[,1], openAlleles[,2], fsep = "")

  femGen <- c('WWWW','HWWW','BWWW','HHWW','BHWW','BBWW',
              'WWBW','HWBW','BWBW','HHBW','BHBW','BBBW',
              'WWBB','HWBB','BWBB','HHBB','BHBB','BBBB')

  malGen <- c('WYWW','HYWW','BYWW','WYBW','HYBW','BYBW','WYBB','HYBB','BYBB')


  #############################################################################
  ## setup all probability lists
  #############################################################################

  # Mendelian
  mend <- list()
  mend$one <- list('W'=c('W'=1),
                   'B'=c('B'=1))
  mend$two <- list('W'=c('W'=1),
                   'B'=c('B'=1))

  # Female homing
  matern <- list()
  matern$one <- list('W'=c('W'=1-cF1,
                           'B'=cF1),
                     'H'=c('H'=1),
                     'B'=c('B'=1))
  matern$two <- list('W'=c('W'=1-cF2,
                           'B'=cF2),
                     'B'=c('B'=1))

  # Male homing
  patern <- list()
  patern$two <- list('W'=c('W'=1-cM2,
                           'B'=cM2),
                     'B'=c('B'=1))


  #############################################################################
  ## fill transition matrix
  #############################################################################

  #use this many times down below
  numGen <- length(femGen) + length(malGen)

  #create transition matrix to fill
  tMatrix <- array(data = 0,
                   dim = c(numGen,numGen,numGen),
                   dimnames = list(c(femGen,malGen),c(femGen,malGen),c(femGen,malGen)))

  #number of alleles, set score lists
  numAlleles <- length(gTypes)
  fScore <- mScore <- logical(length = 1)


  #############################################################################
  ## loop over all matings, female outer loop
  #############################################################################

  # so, future, these probs and the male probs don't change.
  #  ie, we could probably calculate all female things, then all male things, and
  #  store both of those in lists. Finally, combining all elements of those lists
  #  would generate the cube, while calculating each male/female part only once.

  for(fi in femGen){
    #do female stuff here
    #This splits all characters.
    fSplit <- strsplit(x = fi, split = "")[[1]]
    #Score them: is there an H present on allele 1
    fScore[] <- any(grepl(pattern = 'H', x = fSplit[1:2], fixed = TRUE))


    #setup offspring allele lists
    fPHold <- vector(mode = "list", length = numAlleles)

    # check if there is an H allele
    if(fScore){
      # there is an H allele present at allele 1, homing occurs

      # get inheritance things
      fPHold[[1]] <- c(matern$one[[ fSplit[1] ]], matern$one[[ fSplit[2] ]])
      fPHold[[2]] <- c(matern$two[[ fSplit[3] ]], matern$two[[ fSplit[4] ]])

    } else {
      # there is not an H present, so no homing

      # get inheritance things
      fPHold[[1]] <- c(mend$one[[ fSplit[1] ]], mend$one[[ fSplit[2] ]])
      fPHold[[2]] <- c(mend$two[[ fSplit[3] ]], mend$two[[ fSplit[4] ]])

    }# end female scoring


    ###########################################################################
    ## loop over male mate. This is the inner loop
    ###########################################################################
    for(mi in malGen){
      #do male stuff here
      #This splits all characters.
      mSplit <- strsplit(x = mi, split = "")[[1]]
      #Score them: is there an H present on allele 1
      mScore[] <- grepl(pattern = 'H', x = mSplit[1], fixed = TRUE)


      #setup offspring allele lists
      mPHold <- vector(mode = "list", length = numAlleles)

      # check if there is an H allele
      if(mScore){
        # there is an H allele present at allele 1, homing occurs

        # get inheritance things
        mPHold[[1]] <- c('H'=1, 'Y'=1)
        mPHold[[2]] <- c(patern$two[[ mSplit[3] ]], patern$two[[ mSplit[4] ]])

      } else {
        # there is not an H present, so no homing

        # get inheritance things
        mPHold[[1]] <- c(mend$one[[ mSplit[1] ]], 'Y'=1)
        mPHold[[2]] <- c(mend$two[[ mSplit[3] ]], mend$two[[ mSplit[4] ]])

      }# end male scoring


      #########################################################################
      ## Get combinations and put them in the tMatrix. This must be done
      ##  inside the inner loop
      #########################################################################

      holdProbs <- vector(mode = "list", length = numAlleles)

      # loop over alleles
      # This gets all combinations of male/female at both alleles
      for(allele in 1:numAlleles){
        # get combinations at each locus
        hProbs <- expand.grid(fPHold[[allele]], mPHold[[allele]], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        hAllele <- expand.grid(names(fPHold[[allele]]), names(mPHold[[allele]]), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

        # Combine and sort
        hProbs <- hProbs[,1]*hProbs[,2]
        hAllele <- apply(X = hAllele, MARGIN = 1, FUN = sort)
        hAllele <- file.path(hAllele[1,],hAllele[2,], fsep = "")

        # reduce and store
        holdProbs[[allele]] <- vapply(X = unique(hAllele), FUN = function(x){sum(hProbs[hAllele==x])},
                                      FUN.VALUE = numeric(length = 1L))
      }# end loop over alleles

      # expand
      finalProbs <- expand.grid(holdProbs, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      finalGens <- expand.grid(names(holdProbs[[1]]), names(holdProbs[[2]]), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

      # combine
      finalProbs <- finalProbs[ ,1]*finalProbs[ ,2]
      finalGens <- file.path(finalGens[ ,1], finalGens[ ,2], fsep = "")

      # normalize
      finalProbs <- finalProbs/16

      # store
      tMatrix[fi,mi, finalGens ] <- finalProbs

    }# end male loop
  }# end female loop


  #protection from underflow errors
  tMatrix[tMatrix < .Machine$double.eps] <- 0


  ## initialize viability mask. No mother/father-specific death, so use basic mask
  viabilityMask <- array(data = 1, dim = c(numGen,numGen,numGen),
                         dimnames = list(c(femGen,malGen), c(femGen,malGen), c(femGen,malGen)))
  viabilityMask[ , ,c('WWBB','BWBB','BBBB','WYBB','BYBB')] <- 0


  ## genotype-specific modifiers
  if(!is.null(phi)){
    stop("This cube has a special phi, due to being male/female specific.
         Please edit it after the cube is built, if you are absolutely sure it needs to change.")
  }

  phi = setNames(object = c(rep.int(x = 1, times = length(femGen)),
                            rep.int(x = 0, times = length(malGen))),
                 nm = c(femGen,malGen))

  modifiers = cubeModifiers(c(femGen,malGen), eta = eta, phi = phi, omega = omega, xiF = xiF, xiM = xiM, s = s)


  ## put everytWing into a labeled list to return
  return(list(
    ih = tMatrix,
    tau = viabilityMask,
    genotypesID = c(femGen,malGen),
    genotypesN = numGen,
    wildType = c("WWWW","WYWW"),
    eta = modifiers$eta,
    phi = modifiers$phi,
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = "HYWW"
  ))

}
