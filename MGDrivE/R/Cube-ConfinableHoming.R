###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   Confinable Homing
#   Héctor Sanchez, Jared Bennett, Sean Wu, John Marshall
#   jared_bennett@berkeley.edu
#   August 2017
#   December 2018
#    Update to reflect cutting, homing, resistance generation rates
#    Also added male/female homing and crossover differences
#    Did not add male/female deposition.
#
###############################################################################

#' Inheritance Cube: Confinable Homing
#'
#' This function creates a confinable homing construct, it has 4 alleles at the first locus
#' and 3 alleles at the second.
#'  * W: Wild-type
#'  * H: Homing allele
#'  * A: Antidote allele
#'  * R: No-cost resistance allele
#'  * B: Detrimental resistance allele
#'
#' @param cF Cutting efficiency of drive allele at locus 1 in females
#' @param cM Cutting efficiency of drive allele at locus 1 in males
#' @param chF Homing efficiency of drive allele at locus 1 in females
#' @param crF Resistance allele generation rate at locus 1 in females
#' @param chM Homing efficiency of drive allele at locus 1 in males
#' @param crM Resistance allele generation rate at locus 1 in males
#' @param dR Background mutation rate from W and H into R allele in males and females
#' @param dB Background mutation rate from A into B allele in males and females
#' @param crossF Female crossover rate
#' @param crossM Male crossover rate
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
cubeConfinableHoming <- function(cF=1.0, cM=1.0, chF=0, crF=0, chM=0, crM=0,
                                   dR=0, dB=0, crossF=0, crossM=0,
                                   eta=NULL, phi=NULL,omega=NULL, xiF=NULL, xiM=NULL, s=NULL){

  ## safety checks
  if(any(c(cF,cM,chF,crF,chM,crM,dR,dB,crossF,crossM)>1) || any(c(cF,cM,chF,crF,chM,crM,dR,dB,crossF,crossM)<0)){
    stop("Parameters are rates.
         0 <= x <= 1")
  }


  #############################################################################
  ## generate all genotypes, set up vectors and matrices
  #############################################################################

  #list of possible alleles at each locus
  #the first locus has the drive, and can be erased
  # the second locus is the "CHACR" element. No drive, but eracing piece
  gTypes <- list(c("W", "H", "R", "B"), c("W", "A", "B"))

  # this generates genotypes from gTypes. Only needed done once.
  # #get all combinations of locus 1 and locus 2
  # alleles <- expand.grid(gTypes,KEEP.OUT.ATTRS = FALSE,
  #                        stringsAsFactors = FALSE)
  # #paste them together.
  # alleles <- do.call(what = paste0, args = list(alleles[,1], alleles[,2]))
  # #expand all combinations of alleles with locus 1 and locus 2
  # hold <- expand.grid(alleles,alleles,KEEP.OUT.ATTRS = FALSE,
  #                     stringsAsFactors = FALSE)
  # #sort, because order of alleles doesn't matter, paste, and keep unique
  # genotypes <- unique(vapply(X = 1:dim(hold)[1],
  #                            FUN = function(x){
  #                              paste0(sort(x = hold[x, ]), collapse = "")},
  #                            FUN.VALUE = character(1)
  # ))

  genotypes <- c("WWWW","HWWW","RWWW","BWWW","WAWW","HAWW","RAWW","BAWW","WBWW",
                 "HBWW","RBWW","BBWW","HWHW","HWRW","BWHW","HWWA","HAHW","HWRA",
                 "BAHW","HWWB","HBHW","HWRB","BBHW","RWRW","BWRW","RWWA","HARW",
                 "RARW","BARW","RWWB","HBRW","RBRW","BBRW","BWBW","BWWA","BWHA",
                 "BWRA","BABW","BWWB","BWHB","BWRB","BBBW","WAWA","HAWA","RAWA",
                 "BAWA","WAWB","HBWA","RBWA","BBWA","HAHA","HARA","BAHA","HAWB",
                 "HAHB","HARB","BBHA","RARA","BARA","RAWB","HBRA","RARB","BBRA",
                 "BABA","BAWB","BAHB","BARB","BABB","WBWB","HBWB","RBWB","BBWB",
                 "HBHB","HBRB","BBHB","RBRB","BBRB","BBBB")

  #############################################################################
  ## setup all probability lists
  #############################################################################

  #set probabilities for the first locus, the Homing element
  locus1Probs <- list()

  locus1Probs$mendelian <- list("W"=c("W"=1-dR,"R"=dR),
                                "H"=c("H"=1-dR,"R"=dR),
                                "R"=c("R"=1),
                                "B"=c("B"=1))

  locus1Probs$homingM <- list("W"=c("W"=(1-dR)*(1-cM),
                                    "H"=(1-dR)*cM*chM,
                                    "R"=(1-dR)*cM*(1-chM)*crM + dR,
                                    "B"=(1-dR)*cM*(1-chM)*(1-crM)),
                              "H"=c("H"=1-dR,"R"=dR),
                              "R"=c("R"=1),
                              "B"=c("B"=1))

  locus1Probs$homingF <- list("W"=c("W"=(1-dR)*(1-cF),
                                    "H"=(1-dR)*cF*chF,
                                    "R"=(1-dR)*cF*(1-chF)*crF + dR,
                                    "B"=(1-dR)*cF*(1-chF)*(1-crF)),
                              "H"=c("H"=1-dR,"R"=dR),
                              "R"=c("R"=1),
                              "B"=c("B"=1))

  #set probabilities for the second locus, the antidote element
  locus2Probs <- list("W"=c("W"=1),
                      "A"=c("A"=1-dB,"B"=dB),
                      "B"=c("B"=1))



  #############################################################################
  ## fill transition matrix
  #############################################################################

  #use this many times down below
  numGen <- length(genotypes)
  #number of alleles, set score lists
  numAlleles <- 2

  #create transition matrix to fill
  tMatrix <- array(data = 0,
                   dim = c(numGen,numGen,numGen),
                   dimnames = list(genotypes,genotypes,genotypes))


  #############################################################################
  ## loop over all matings, female outer loop
  #############################################################################
  for (fi in 1:numGen){
    #do female stuff here
    #This splits all characters.
    fSplit <- strsplit(x = genotypes[fi], split = "")[[1]]
    #make a list of each allele at every locus. This list is length nmPlex, and each
    # sublist has length 2
    momAlleles <- list(fSplit[1:2], fSplit[3:4])
    #Score for H in either allele, always at locus 1 though
    fScore <- grepl(pattern = "H", x = genotypes[fi], fixed = TRUE)
    #setup offspring allele lists
    fAllele <- rep(x = list(vector(mode = "list", 2)), numAlleles)
    fProbs <- rep(x = list(vector(mode = "list", 2)), numAlleles)


    #set mother alleles
    if(fScore){
      #homing in females
      #lists of allele letter and probabilities
      for (i in 1:numAlleles){

        #female locus 1
        fProbs[[i]][[1]] <- locus1Probs$homingF[[ momAlleles[[i]][[1]] ]]
        fAllele[[i]][[1]] <- names(fProbs[[i]][[1]])

        #female locus 2
        fProbs[[i]][[2]] <- locus2Probs[[ momAlleles[[i]][[2]] ]]
        fAllele[[i]][[2]] <- names(fProbs[[i]][[2]])

      }#end loop over alleles
    } else {
      #no homing in females
      #lists of allele letter and probabilities
      for (i in 1:numAlleles){

        #female locus 1
        fProbs[[i]][[1]] <- locus1Probs$mendelian[[momAlleles[[i]][[1]] ]]
        fAllele[[i]][[1]] <- names(fProbs[[i]][[1]])

        #female locus 2
        fProbs[[i]][[2]] <- locus2Probs[[ momAlleles[[i]][[2]] ]]
        fAllele[[i]][[2]] <- names(fProbs[[i]][[2]])

      }#end loop over alleles
    }#end female check



    #perform female crossover
    #because this is done recursively, I need to hold the value out
    holdA <- fAllele[[1]][[2]]
    holdP <- fProbs[[1]][[2]]

    #swap allele 2's in both
    fAllele[[1]][[2]] <- c(holdA, fAllele[[2]][[2]])
    fAllele[[2]][[2]] <- c(fAllele[[2]][[2]], holdA)

    #set the probs for the new allele 2's in both
    fProbs[[1]][[2]] <- c(holdP*(1-crossF), fProbs[[2]][[2]]*crossF)
    fProbs[[2]][[2]] <- c(fProbs[[2]][[2]]*(1-crossF), holdP*crossF)



    #combine loci for each female allele.
    # This requires looping through each locus, getting all combinations of
    fLociA <- fLociP <- vector(mode = "list", length = numAlleles)

    for( i in 1:numAlleles){

      #get all combinations of loci for each allele, keep male/female separate still
      holdAllF <- expand.grid(fAllele[[i]], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      holdProbF <- expand.grid(fProbs[[i]], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

      #paste alleles together and store in list
      fLociA[[i]] <- do.call(what = "paste0", list(holdAllF[ ,1], holdAllF[ ,2]))
      fLociP[[i]] <- holdProbF[ ,1]*holdProbF[ ,2]
    }





    ###########################################################################
    ## loop over male mate. This is the inner loop
    ###########################################################################
    for (mi in 1:numGen){
      #do male stuff here
      #split male genotype
      #This splits all characters.
      mSplit <- strsplit(x = genotypes[mi], split = "")[[1]]

      #make a list of each allele at every locus. This list is length nmPlex, and each
      # sublist has length 2
      dadAlleles <- list(mSplit[1:2], mSplit[3:4])

      #Score them
      mScore <- grepl(pattern = "H", x = genotypes[mi], fixed = TRUE)

      #setup offspring allele lists
      mAllele <- rep(x = list(vector(mode = "list", 2)), numAlleles)
      mProbs <- rep(x = list(vector(mode = "list", 2)), numAlleles)


      #set father alleles
      if(mScore){
        #homing in males
        #lists of allele letter and probabilities
        for (i in 1:numAlleles){

          #male locus 1
          mProbs[[i]][[1]] <- locus1Probs$homingM[[ dadAlleles[[i]][[1]] ]]
          mAllele[[i]][[1]] <- names(mProbs[[i]][[1]])

          #male locus 2
          mProbs[[i]][[2]] <- locus2Probs[[ dadAlleles[[i]][[2]] ]]
          mAllele[[i]][[2]] <- names(mProbs[[i]][[2]])

        }#end loop over alleles
      } else {
        #no homing in males
        #lists of allele letter and probabilities
        for (i in 1:numAlleles){

          #male locus 1
          mProbs[[i]][[1]] <- locus1Probs$mendelian[[ dadAlleles[[i]][[1]] ]]
          mAllele[[i]][[1]] <- names(mProbs[[i]][[1]])

          #male locus 2
          mProbs[[i]][[2]] <- locus2Probs[[ dadAlleles[[i]][[2]] ]]
          mAllele[[i]][[2]] <- names(mProbs[[i]][[2]])

        }#end loop over alleles
      }#end male checks

      #perform male crossover
      #because this is done recursively, I need to hold the value out
      holdA <- mAllele[[1]][[2]]
      holdP <- mProbs[[1]][[2]]

      #swap allele 2's in both
      mAllele[[1]][[2]] <- c(holdA, mAllele[[2]][[2]])
      mAllele[[2]][[2]] <- c(mAllele[[2]][[2]], holdA)

      #set the probs for the new allele 2's in both
      mProbs[[1]][[2]] <- c(holdP*(1-crossM), mProbs[[2]][[2]]*crossM)
      mProbs[[2]][[2]] <- c(mProbs[[2]][[2]]*(1-crossM), holdP*crossM)



      #########################################################################
      ## Get combinations and put them in the tMatrix. This must be done
      ##  inside the inner loop
      #########################################################################

      #combine loci for each male allele.
      # This requires looping through each locus, getting all combinations of
      mLociA <- mLociP <- vector(mode = "list", length = numAlleles)

      for( i in 1:numAlleles){

        #get all combinations of loci for each allele, keep male/female separate still
        holdAllM <- expand.grid(mAllele[[i]], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        holdProbM <- expand.grid(mProbs[[i]], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

        #paste alleles together and store in list
        mLociA[[i]] <- do.call(what = "paste0", list(holdAllM[ ,1], holdAllM[ ,2]))
        mLociP[[i]] <- holdProbM[ ,1]*holdProbM[ ,2]
      }


      #combine each Allele into single lists, so that alleles within one parent can't
      # be combined with each other, but do get combined with alleles for the
      # other parent
      # ie, unlist the parent alleles, then combine with the other parent
      holdAllOne <- expand.grid(unlist(fLociA), unlist(mLociA), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      holdProbOne <- expand.grid(unlist(fLociP), unlist(mLociP), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

      #sort each combination so they are the same.
      holdAllOne <- apply(X = holdAllOne, MARGIN = 1, FUN = sort, method = 'radix')

      #paste alleles together
      holdAllTwo <- do.call(what = "paste0", list(holdAllOne[1, ], holdAllOne[2, ]))
      holdProbTwo <- holdProbOne[ ,1]*holdProbOne[ ,2]

      #aggregate and return
      aggregateHold <- vapply(X = unique(holdAllTwo), FUN = function(x){
        sum(holdProbTwo[holdAllTwo==x])},
        FUN.VALUE = numeric(length = 1L))

      #normalize
      if(!sum(aggregateHold)==0){
        aggregateHold <- aggregateHold/sum(aggregateHold)
      }

      #set values in tMatrix
      tMatrix[fi,mi, names(aggregateHold) ] <- aggregateHold

    }# end male loop
  }# end female loop



  ## set the other half of the matrix
  tMatrix[tMatrix < .Machine$double.eps] <- 0 #protection from underflow errors

  ## initialize viability mask. No mother-specific death.
  viabilityMask <- array(data = 1, dim = c(numGen,numGen,numGen), dimnames = list(genotypes, genotypes, genotypes))

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
    releaseType = "HAHA"
  ))
}










# older version
# may be helpful if/when translation into Cpp

# cubeConfineableHoming <- function(eH=1.0, pRH=0, pBH=0, dR=0, dB=0, cross=0,
#                                    eta=NULL, phi=NULL,omega=NULL, xiF=NULL, xiM=NULL, s=NULL){
#
#   ## safety checks in case someone is dumb
#   if(any(c(eH,pRH,pBH,dR,dB)>1) || any(c(eH,pRH,pBH,dR,dB)<0)){
#     stop("Parameters are rates.
#          0 <= x <= 1")
#   }
#   if(eH > 1){
#     stop("The homing parameter must be less than 1.")
#   }
#   if(pRH+pBH>1){
#     stop("NHEJ rates (R and B alleles) must sum to less than 1.")
#   }
#
#   #############################################################################
#   ## generate all genotypes, set up vectors and matrices
#   #############################################################################
#
#   #list of possible alleles at each locus
#   #the first locus has the drive, and can be erased
#   # the second locus is the "CHACR" element. No drive, but eracing piece
#   gTypes <- list(c("W", "H", "R", "B"), c("W", "A", "B"))
#
#   #get all combinations of locus 1 and locus 2
#   alleles <- expand.grid(gTypes,KEEP.OUT.ATTRS = FALSE,
#                          stringsAsFactors = FALSE)
#   #paste them together.
#   alleles <- do.call(what = paste0, args = list(alleles[,1], alleles[,2]))
#   #expand all combinations of alleles with locus 1 and locus 2
#   hold <- expand.grid(alleles,alleles,KEEP.OUT.ATTRS = FALSE,
#                       stringsAsFactors = FALSE)
#   #sort, because order of alleles doesn't matter, paste, and keep unique
#   genotypes <- unique(vapply(X = 1:dim(hold)[1],
#                              FUN = function(x){
#                                paste0(sort(x = hold[x, ]), collapse = "")},
#                              FUN.VALUE = character(1)
#   ))
#
#
#   #############################################################################
#   ## setup all probability lists
#   #############################################################################
#
#   #set probabilities for the first locus, the Homing and Erasing element
#   locus1Probs <- list()
#
#   locus1Probs$mendelian$W <- setNames(object = c(1-dR,dR), nm = c("W", "R"))
#   locus1Probs$mendelian$H <- setNames(object = c(1-dR,dR), nm = c("H", "R"))
#   locus1Probs$mendelian$R <- setNames(object = 1, nm = "R")
#   locus1Probs$mendelian$B <- setNames(object = 1, nm = "B")
#
#   locus1Probs$homing$W <- setNames(object = c(1-dR-eH, eH*(1-pRH-pBH),dR+eH*pRH, eH*pBH),
#                                    nm = c("W", "H", "R", "B"))
#   locus1Probs$homing$H <- setNames(object = c(1-dR,dR), nm = c("H", "R"))
#   locus1Probs$homing$R <- setNames(object = 1, nm = "R")
#   locus1Probs$homing$B <- setNames(object = 1, nm = "B")
#
#   #set probabilities for the second locus, the CHACR element
#   locus2Probs <- list()
#
#   locus2Probs$W <- setNames(object = 1, nm = "W")
#   locus2Probs$A <- setNames(object = c(1-dB, dB), nm = c("A", "B"))
#   locus2Probs$B <- setNames(object = 1, nm = "B")
#
#
#   #############################################################################
#   ## fill transition matrix
#   #############################################################################
#
#   #use this many times down below
#   numGen <- length(genotypes)
#   #number of alleles, set score lists
#   numAlleles <- length(gTypes)
#
#   #create transition matrix to fill
#   tMatrix <- array(data = 0,
#                    dim = c(numGen,numGen,numGen),
#                    dimnames = list(genotypes,genotypes,genotypes))
#
#
#   #############################################################################
#   ## loop over all matings, female outer loop
#   #############################################################################
#   for (fi in 1:numGen){
#     #do female stuff here
#     #This splits all characters.
#     fSplit <- strsplit(x = genotypes[fi], split = "")[[1]]
#     #make a list of each allele at every locus. This list is length nmPlex, and each
#     # sublist has length 2
#     momAlleles <- list(fSplit[1:2], fSplit[3:4])
#     #Score for H in either allele, always at locus 1 though
#     fScore <- grepl(pattern = "H", x = genotypes[fi], fixed = TRUE)
#     #setup offspring allele lists
#     fAllele <- rep(x = list(vector(mode = "list", 2)), numAlleles)
#     fProbs <- rep(x = list(vector(mode = "list", 2)), numAlleles)
#
#
#     #set mother alleles
#     if(fScore){
#       #homing in females
#       #lists of allele letter and probabilities
#       for (i in 1:numAlleles){
#
#         #female locus 1
#         if(momAlleles[[i]][1]=="W"){
#           fAllele[[i]][[1]] <- c("W", "H", "R", "B")
#           fProbs[[i]][[1]] <- locus1Probs$homing$W
#         } else if(momAlleles[[i]][1]=="H"){
#           fAllele[[i]][[1]] <- c("H", "R")
#           fProbs[[i]][[1]] <- locus1Probs$homing$H
#         } else if(momAlleles[[i]][1]=="R"){
#           fAllele[[i]][[1]] <- "R"
#           fProbs[[i]][[1]] <- locus1Probs$homing$R
#         } else if(momAlleles[[i]][1]=="B"){
#           fAllele[[i]][[1]] <- "B"
#           fProbs[[i]][[1]] <- locus1Probs$homing$B
#         }#end locus 1
#
#         #female locus 2
#         if(momAlleles[[i]][2]=="W"){
#           fAllele[[i]][[2]] <- "W"
#           fProbs[[i]][[2]] <- locus2Probs$W
#         } else if(momAlleles[[i]][2]=="A"){
#           fAllele[[i]][[2]] <- c("A", "B")
#           fProbs[[i]][[2]] <- locus2Probs$A
#         } else if(momAlleles[[i]][2]=="B"){
#           fAllele[[i]][[2]] <- "B"
#           fProbs[[i]][[2]] <- locus2Probs$B
#         }#end locus 2
#
#       }#end loop over alleles
#     } else {
#       #no homing in females
#       #lists of allele letter and probabilities
#       for (i in 1:numAlleles){
#
#         #female locus 1
#         if(momAlleles[[i]][[1]]=="W"){
#           fAllele[[i]][[1]] <- c("W", "R")
#           fProbs[[i]][[1]] <- locus1Probs$mendelian$W
#         } else if(momAlleles[[i]][1]=="R"){
#           fAllele[[i]][[1]] <- "R"
#           fProbs[[i]][[1]] <- locus1Probs$mendelian$R
#         } else if(momAlleles[[i]][1]=="B"){
#           fAllele[[i]][[1]] <- "B"
#           fProbs[[i]][[1]] <- locus1Probs$mendelian$B
#         }#end locus 1
#
#         #female locus 2
#         if(momAlleles[[i]][2]=="W"){
#           fAllele[[i]][[2]] <- "W"
#           fProbs[[i]][[2]] <- locus2Probs$W
#         } else if(momAlleles[[i]][2]=="A"){
#           fAllele[[i]][[2]] <- c("A", "B")
#           fProbs[[i]][[2]] <- locus2Probs$A
#         } else if(momAlleles[[i]][2]=="B"){
#           fAllele[[i]][[2]] <- "B"
#           fProbs[[i]][[2]] <- locus2Probs$B
#         }#end locus 2
#
#       }#end loop over alleles
#     }#end female check
#
#     #perform female crossover
#     #because this is done recursively, I need to hold the value out
#     holdA <- fAllele[[1]][[2]]
#     holdP <- fProbs[[1]][[2]]
#
#     #swap allele 2's in both
#     fAllele[[1]][[2]] <- c(holdA, fAllele[[2]][[2]])
#     fAllele[[2]][[2]] <- c(fAllele[[2]][[2]], holdA)
#
#     #set the probs for the new allele 2's in both
#     fProbs[[1]][[2]] <- c(holdP*(1-cross), fProbs[[2]][[2]]*cross)
#     fProbs[[2]][[2]] <- c(fProbs[[2]][[2]]*(1-cross), holdP*cross)
#
#
#
#
#
#     ###########################################################################
#     ## loop over male mate. This is the inner loop
#     ###########################################################################
#     for (mi in 1:fi){ #make this mi in 1 to fi
#       #do male stuff here
#       #split male genotype
#       #This splits all characters.
#       mSplit <- strsplit(x = genotypes[mi], split = "")[[1]]
#
#       #make a list of each allele at every locus. This list is length nmPlex, and each
#       # sublist has length 2
#       dadAlleles <- list(mSplit[1:2], mSplit[3:4])
#
#       #Score them
#       mScore <- grepl(pattern = "H", x = genotypes[mi], fixed = TRUE)
#
#       #setup offspring allele lists
#       mAllele <- rep(x = list(vector(mode = "list", 2)), numAlleles)
#       mProbs <- rep(x = list(vector(mode = "list", 2)), numAlleles)
#
#
#       #set father alleles
#       if(mScore){
#         #homing in females
#         #lists of allele letter and probabilities
#         for (i in 1:numAlleles){
#
#           #male locus 1
#           if(dadAlleles[[i]][1]=="W"){
#             mAllele[[i]][[1]] <- c("W", "H", "R", "B")
#             mProbs[[i]][[1]] <- locus1Probs$homing$W
#           } else if(dadAlleles[[i]][1]=="H"){
#             mAllele[[i]][[1]] <- c("H", "R")
#             mProbs[[i]][[1]] <- locus1Probs$homing$H
#           } else if(dadAlleles[[i]][1]=="R"){
#             mAllele[[i]][[1]] <- "R"
#             mProbs[[i]][[1]] <- locus1Probs$homing$R
#           } else if(dadAlleles[[i]][1]=="B"){
#             mAllele[[i]][[1]] <- "B"
#             mProbs[[i]][[1]] <- locus1Probs$homing$B
#           }#end locus 1
#
#           #male locus 2
#           if(dadAlleles[[i]][2]=="W"){
#             mAllele[[i]][[2]] <- "W"
#             mProbs[[i]][[2]] <- locus2Probs$W
#           } else if(dadAlleles[[i]][2]=="A"){
#             mAllele[[i]][[2]] <- c("A", "B")
#             mProbs[[i]][[2]] <- locus2Probs$A
#           } else if(dadAlleles[[i]][2]=="B"){
#             mAllele[[i]][[2]] <- "B"
#             mProbs[[i]][[2]] <- locus2Probs$B
#           }#end locus 2
#
#         }#end loop over alleles
#       } else {
#         #no homing in females
#         #lists of allele letter and probabilities
#         for (i in 1:numAlleles){
#
#           #male locus 1
#           if(dadAlleles[[i]][1]=="W"){
#             mAllele[[i]][[1]] <- c("W", "R")
#             mProbs[[i]][[1]] <- locus1Probs$mendelian$W
#           } else if(dadAlleles[[i]][1]=="R"){
#             mAllele[[i]][[1]] <- "R"
#             mProbs[[i]][[1]] <- locus1Probs$mendelian$R
#           } else if(dadAlleles[[i]][1]=="B"){
#             mAllele[[i]][[1]] <- "B"
#             mProbs[[i]][[1]] <- locus1Probs$mendelian$B
#           }#end locus 1
#
#           #male locus 2
#           if(dadAlleles[[i]][2]=="W"){
#             mAllele[[i]][[2]] <- "W"
#             mProbs[[i]][[2]] <- locus2Probs$W
#           } else if(dadAlleles[[i]][2]=="A"){
#             mAllele[[i]][[2]] <- c("A", "B")
#             mProbs[[i]][[2]] <- locus2Probs$A
#           } else if(dadAlleles[[i]][2]=="B"){
#             mAllele[[i]][[2]] <- "B"
#             mProbs[[i]][[2]] <- locus2Probs$B
#           }#end locus 2
#
#         }#end loop over alleles
#       }#end male checks
#
#       #perform male crossover
#       #because this is done recursively, I need to hold the value out
#       holdA <- mAllele[[1]][[2]]
#       holdP <- mProbs[[1]][[2]]
#
#       #swap allele 2's in both
#       mAllele[[1]][[2]] <- c(holdA, mAllele[[2]][[2]])
#       mAllele[[2]][[2]] <- c(mAllele[[2]][[2]], holdA)
#
#       #set the probs for the new allele 2's in both
#       mProbs[[1]][[2]] <- c(holdP*(1-cross), mProbs[[2]][[2]]*cross)
#       mProbs[[2]][[2]] <- c(mProbs[[2]][[2]]*(1-cross), holdP*cross)
#
#
#
#       #########################################################################
#       ## Get combinations and put them in the tMatrix. This must be done
#       ##  inside the inner loop
#       #########################################################################
#
#       #combine loci for each allele.
#       # This requires looping through each locus, getting all combinations of
#       fLociA <- fLociP <- mLociA <- mLociP <- vector(mode = "list", length = numAlleles)
#
#       for( i in 1:numAlleles){
#
#         #get all combinations of loci for each allele, keep male/female separate still
#         holdAllF <- expand.grid(fAllele[[i]], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
#         holdProbF <- expand.grid(fProbs[[i]], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
#
#         holdAllM <- expand.grid(mAllele[[i]], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
#         holdProbM <- expand.grid(mProbs[[i]], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
#
#         #paste alleles together and store in list
#         fLociA[[i]] <- do.call(what = "paste0", list(holdAllF[ ,1], holdAllF[ ,2]))
#         fLociP[[i]] <- holdProbF[ ,1]*holdProbF[ ,2]
#
#         mLociA[[i]] <- do.call(what = "paste0", list(holdAllM[ ,1], holdAllM[ ,2]))
#         mLociP[[i]] <- holdProbM[ ,1]*holdProbM[ ,2]
#       }
#
#
#       #combine each Allele into single lists, so that alleles within one parent can't
#       # be combined with each other, but do get combined with alleles for the
#       # other parent
#       # ie, unlist the parent alleles, then combine with the other parent
#       holdAllOne <- expand.grid(unlist(fLociA), unlist(mLociA), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
#       holdProbOne <- expand.grid(unlist(fLociP), unlist(mLociP), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
#
#       #sort each combination so they are the same.
#       holdAllOne <- apply(X = holdAllOne, MARGIN = 1, FUN = sort)
#
#       #paste alleles togheter
#       holdAllTwo <- do.call(what = "paste0", list(holdAllOne[1, ], holdAllOne[2, ]))
#       holdProbTwo <- holdProbOne[ ,1]*holdProbOne[ ,2]
#
#       #aggregate and return
#       aggregateHold <- vapply(X = unique(holdAllTwo), FUN = function(x){
#         sum(holdProbTwo[holdAllTwo==x])},
#         FUN.VALUE = numeric(length = 1L))
#
#       #normalize
#       if(!sum(aggregateHold)==0){
#         aggregateHold <- aggregateHold/sum(aggregateHold)
#       }
#
#       #set values in tMatrix
#       tMatrix[fi,mi, names(aggregateHold) ] <- aggregateHold
#
#     }# end male loop
#   }# end female loop
#
#
#
#   ## set the other half of the matrix
#   symCube(lowerMat = tMatrix)
#   tMatrix[tMatrix < .Machine$double.eps] <- 0 #protection from underflow errors
#
#   ## initialize viability mask. No mother-specific death.
#   viabilityMask <- array(data = 1L, dim = c(numGen,numGen,numGen), dimnames = list(genotypes, genotypes, genotypes))
#
#   ## genotype-specific modifiers
#   modifiers = cubeModifiers(genotypes, eta = eta, phi = phi, omega = omega, xiF = xiF, xiM = xiM, s = s)
#
#   ## put everything into a labeled list to return
#   return(list(
#     ih = tMatrix,
#     tau = viabilityMask,
#     genotypesID = genotypes,
#     genotypesN = numGen,
#     wildType = "WWWW",
#     eta = modifiers$eta,
#     phi = modifiers$phi,
#     omega = modifiers$omega,
#     xiF = modifiers$xiF,
#     xiM = modifiers$xiM,
#     s = modifiers$s,
#     releaseType = "HAHA"
#   ))
# }
