###############################################################################
#      ______      __
#     / ____/_  __/ /_  ___
#    / /   / / / / __ \/ _ \
#   / /___/ /_/ / /_/ /  __/
#   \____/\__,_/_.___/\___/
#
#   MGDrivE: Mosquito Gene Drive Explorer
#   ECHACR
#   Héctor Sanchez, Jared Bennett, Sean Wu, John Marshall
#   jared_bennett@berkeley.edu
#   August 2017
#   December 2018
#    Update to reflect cutting, homing, resistance generation rates
#    Removed many of the locus specific parameters, not sure we can estimate them
#
###############################################################################

#' Inheritance Cube: ECHACR
#'
#' This function creates an ECHACR construct, it has 5 alleles at the first locus
#' and 4 alleles at the second.
#'  * W: Wild-type
#'  * H: Homing allele
#'  * E: Eraser allele
#'  * R: No-cost resistance allele
#'  * B: Detrimental resistance allele
#'  * cHW: Rate of homing from H, W -> H transition
#'  * cEH: Rate of homing from E, H -> E transition
#'  * cEW: Rate of homing from E, W -> E transition
#'
#' This inheritance pattern corresponds to the [Active Genetic Neutralizing Elements for Halting or Deleting Gene Drives](https://doi.org/10.1016/j.molcel.2020.09.003) publication.
#'
#' @param cHW Cutting efficiency of drive allele at locus 1
#' @param cEW Cutting efficiency of ECHACR element into W
#' @param cEH Cutting efficiency of ECHACR element into H
#' @param chHW Homing efficiency of drive allele at locus 1
#' @param crHW Resistance allele efficiency of drive allele at locus 1
#' @param ceEW Homing efficiency of ECHACR element into W
#' @param crEW Resistance allele efficiency of ECHACR element into W
#' @param ceEH Homing efficiency of ECHACR element into H
#' @param crEH Resistance allele efficiency of ECHACR element into H
#' @param d1 Background mutation rate from W into R allele
#' @param d2 Background mutation rate from H into R allele
#' @param d3 Background mutation rate from E into R allele
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
cubeECHACR <- function(cHW=1.0, cEW=1.0, cEH=1.0, chHW=0, crHW=0, ceEW=0, crEW=0,
                        ceEH=0, crEH=0, d1=0, d2=0, d3=0,
                        eta=NULL, phi=NULL,omega=NULL, xiF=NULL, xiM=NULL, s=NULL){


  ## safety checks
  if(any(c(cHW,cEW,cEH,chHW,crHW,ceEW,crEW,ceEH,crEH,d1,d2,d3)>1) || any(c(cHW,cEW,cEH,chHW,crHW,ceEW,crEW,ceEH,crEH,d1,d2,d3)<0)){
    stop("Parameters are rates.\n0 <= x <= 1")
  }


  #############################################################################
  ## generate all genotypes, set up vectors and matrices
  #############################################################################
  # 20200116
  # Notes about this construct
  #  The drive that Ethan et al built for an X-linked ECHACR is substantially different
  #  than the construct that was previously described. Below is an autosomal version
  #  of the first x-linked version. This version assumes independence between the loci,
  #  ie, they are on different chromosomes or more than 50cM apart on the same chromosome.
  # Genotypes are 4 characters long: the first 2 characters are the first locus, the second
  #  2 characters are the second locus.
  # Homing occurs if the H is present at the first locus, and it only targets the
  #  first locus.
  # ERACING/CHACING occurs if E is present at the second locus. It copies itself, whole,
  #  to the second locus, and disrupts the first locus - assumption is it cuts out
  #  the construct.
  #
  # This is very different from how the new X-linked version internally. The logic
  #  is the same, but to allow linkage and crossovers, the implementation is different
  #  and the genotypes are organized differently.



  #list of possible alleles at each locus
  #the first locus has the drive, and can be erased
  # the second locus is the "CHACR" element. No drive, but eracing piece
  gTypes <- list(c("W", "H", "E", "R", "B"), c("W", "E", "R", "B"))

  # this generates genotypes from gTypes. Only needed done once.
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
  #
  # #get all combinations of both loci
  # openAlleles <- expand.grid(alleleList, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  # #paste them together. This is the final list of genotypes
  # genotypes <- do.call(what = paste0, args = list(openAlleles[,1], openAlleles[,2]))

  genotypes <- c("WWWW","HWWW","EWWW","RWWW","BWWW","HHWW","EHWW","HRWW","BHWW",
                 "EEWW","ERWW","BEWW","RRWW","BRWW","BBWW","WWEW","HWEW","EWEW",
                 "RWEW","BWEW","HHEW","EHEW","HREW","BHEW","EEEW","EREW","BEEW",
                 "RREW","BREW","BBEW","WWRW","HWRW","EWRW","RWRW","BWRW","HHRW",
                 "EHRW","HRRW","BHRW","EERW","ERRW","BERW","RRRW","BRRW","BBRW",
                 "WWBW","HWBW","EWBW","RWBW","BWBW","HHBW","EHBW","HRBW","BHBW",
                 "EEBW","ERBW","BEBW","RRBW","BRBW","BBBW","WWEE","HWEE","EWEE",
                 "RWEE","BWEE","HHEE","EHEE","HREE","BHEE","EEEE","EREE","BEEE",
                 "RREE","BREE","BBEE","WWER","HWER","EWER","RWER","BWER","HHER",
                 "EHER","HRER","BHER","EEER","ERER","BEER","RRER","BRER","BBER",
                 "WWBE","HWBE","EWBE","RWBE","BWBE","HHBE","EHBE","HRBE","BHBE",
                 "EEBE","ERBE","BEBE","RRBE","BRBE","BBBE","WWRR","HWRR","EWRR",
                 "RWRR","BWRR","HHRR","EHRR","HRRR","BHRR","EERR","ERRR","BERR",
                 "RRRR","BRRR","BBRR","WWBR","HWBR","EWBR","RWBR","BWBR","HHBR",
                 "EHBR","HRBR","BHBR","EEBR","ERBR","BEBR","RRBR","BRBR","BBBR",
                 "WWBB","HWBB","EWBB","RWBB","BWBB","HHBB","EHBB","HRBB","BHBB",
                 "EEBB","ERBB","BEBB","RRBB","BRBB","BBBB")


  #############################################################################
  ## setup all probability lists
  #############################################################################

  #set probabilities for the first locus, the Homing and Erasing element
  locus1Probs <- list()

  locus1Probs$mendelian <- list("W"=c("W"=1-d1,"R"=d1),
                                "H"=c("H"=1-d2,"R"=d2),
                                "E"=c("E"=1-d3,"R"=d3),
                                "R"=c("R"=1),
                                "B"=c("B"=1))

  locus1Probs$homingOne <- list("W"=c("W"=(1-d1)*(1-cHW),
                                      "H"=(1-d1)*cHW*chHW,
                                      "R"=d1 + (1-d1)*cHW*(1-chHW)*crHW,
                                      "B"=(1-d1)*(1-chHW)*cHW*(1-chHW)*(1-crHW)),
                                "H"=c("H"=1-d2,"R"=d2),
                                "E"=c("E"=1-d3,"R"=d3),
                                "R"=c("R"=1),
                                "B"=c("B"=1))

  locus1Probs$homingMixed <- list("W"=c("W"=(1-d1)*(1-cHW)*(1-cEW),
                                        "H"=(1-d1)*cHW*chHW*(1-cEH),
                                        "E"=(1-d1)*(1-cHW)*cEW*ceEW + (1-d1)*cHW*chHW*cEH*ceEH,
                                        "R"=d1 + (1-d1)*cHW*(1-chHW)*crHW + (1-d1)*(1-cHW)*cEW*(1-ceEW)*crEW + (1-d1)*cHW*chHW*cEH*(1-ceEH)*crEH,
                                        "B"=(1-d1)*cHW*(1-chHW)*(1-crHW) + (1-d1)*(1-cHW)*cEW*(1-ceEW)*(1-crEW) + (1-d1)*cHW*chHW*cEH*(1-ceEH)*(1-crEH)),
                                  "H"=c("H"=(1-d2)*(1-cEH),
                                        "E"=(1-d2)*cEH*ceEH,
                                        "R"=d2 + (1-d2)*cEH*(1-ceEH)*crEH,
                                        "B"=(1-d2)*cEH*(1-ceEH)*(1-crEH)),
                                  "E"=c("E"=1-d3,"R"=d3),
                                  "R"=c("R"=1),
                                  "B"=c("B"=1))

  locus1Probs$homingTwo <- list("W"=c("W"=1-d1,"R"=d1),
                                "H"=c("H"=(1-d2)*(1-cEH),
                                      "E"=(1-d2)*cEH*ceEH,
                                      "R"=d2 + (1-d2)*cEH*(1-ceEH)*crEH,
                                      "B"=(1-d2)*cEH*(1-ceEH)*(1-crEH)),
                                "E"=c("E"=1-d3,"R"=d3),
                                "R"=c("R"=1),
                                "B"=c("B"=1))


  #set probabilities for the second locus, the CHACR element
  locus2Probs <- list()

  locus2Probs$mendelian <- list("W"=c("W"=1-d1,"R"=d1),
                                "E"=c("E"=1-d3,"R"=d3),
                                "R"=c("R"=1),
                                "B"=c("B"=1))

  locus2Probs$homing <- list("W"=c("W"=(1-d1)*(1-cEW),
                                   "E"=(1-d1)*cEW*ceEW,
                                   "R"=d1 + (1-d1)*cEW*(1-ceEW)*crEW,
                                   "B"=(1-d1)*cEW*(1-ceEW)*(1-crEW)),
                             "E"=c("E"=1-d3,"R"=d3),
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

  #number of alleles, set score lists
  numAlleles <- length(gTypes)
  fScore <- mScore <- vector(mode = "list", length = numAlleles)

  #############################################################################
  ## loop over all matings, female outer loop
  #############################################################################
  for (fi in 1:length(genotypes)){
    #do female stuff here
    #This splits all characters.
    fSplit <- strsplit(x = genotypes[fi], split = "")[[1]]
    #make a list of each allele at every locus. This list is length nmPlex, and each
    # sublist has length 2
    momAlleles <- list(fSplit[1:2], fSplit[3:4])
    #Score them
    fScore[[1]] <- grepl(pattern = "H", x = fSplit[1:2], fixed = TRUE)
    fScore[[2]] <- grepl(pattern = "E", x = fSplit[3:4], fixed = TRUE)
    #setup offspring allele lists
    fAllele <- rep(x = list(vector(mode = "list", 2)), numAlleles)
    fProbs <- rep(x = list(vector(mode = "list", 2)), numAlleles)



    #set mother alleles
    if(all(!fScore[[1]])){
      #FF and FF||FT||TF||TT
      #loop over 2 alleles at each locus
      for(j in 1:2){
        #at locus one, what is there, fill it

        if(momAlleles[[1]][[j]]=="W"){
          fAllele[[1]][[j]] <- c("W", "R")
          fProbs[[1]][[j]] <- locus1Probs$mendelian$W
        } else if(momAlleles[[1]][[j]]=="E"){
          fAllele[[1]][[j]] <- c("E", "R")
          fProbs[[1]][[j]] <- locus1Probs$mendelian$E
        } else if(momAlleles[[1]][[j]]=="R"){
          fAllele[[1]][[j]] <- "R"
          fProbs[[1]][[j]] <- locus1Probs$mendelian$R
        }else if(momAlleles[[1]][[j]]=="B"){
          fAllele[[1]][[j]] <- "B"
          fProbs[[1]][[j]] <- locus1Probs$mendelian$B
        }#end locus 1 if statements

        if(momAlleles[[2]][[j]]=="W"){
          fAllele[[2]][[j]] <- c("W", "R")
          fProbs[[2]][[j]] <- locus2Probs$mendelian$W
        } else if(momAlleles[[2]][[j]]=="E"){
          fAllele[[2]][[j]] <- c("E", "R")
          fProbs[[2]][[j]] <- locus2Probs$mendelian$E
        } else if(momAlleles[[2]][[j]]=="R"){
          fAllele[[2]][[j]] <- "R"
          fProbs[[2]][[j]] <- locus2Probs$mendelian$R
        }else if(momAlleles[[2]][[j]]=="B"){
          fAllele[[2]][[j]] <- "B"
          fProbs[[2]][[j]] <- locus2Probs$mendelian$B
        }#end locus 2 if statements
      }#end loop over alleles
    } else if(xor(fScore[[1]][1],fScore[[1]][2]) && all(!fScore[[2]])){
      #TF||FT and FF
      #loop over 2 alleles at each locus
      for(j in 1:2){
        #at locus one, what is there, fill it

        if(momAlleles[[1]][[j]]=="W"){
          fAllele[[1]][[j]] <- c("W", "H", "R", "B")
          fProbs[[1]][[j]] <- locus1Probs$homingOne$W
        } else if(momAlleles[[1]][[j]]=="H"){
          fAllele[[1]][[j]] <- c("H","R")
          fProbs[[1]][[j]] <- locus1Probs$homingOne$H
        } else if(momAlleles[[1]][[j]]=="E"){
          fAllele[[1]][[j]] <- c("E", "R")
          fProbs[[1]][[j]] <- locus1Probs$homingOne$E
        } else if(momAlleles[[1]][[j]]=="R"){
          fAllele[[1]][[j]] <- "R"
          fProbs[[1]][[j]] <- locus1Probs$homingOne$R
        }else if(momAlleles[[1]][[j]]=="B"){
          fAllele[[1]][[j]] <- "B"
          fProbs[[1]][[j]] <- locus1Probs$homingOne$B
        }#end locus 1 if statements

        if(momAlleles[[2]][[j]]=="W"){
          fAllele[[2]][[j]] <- c("W", "R")
          fProbs[[2]][[j]] <- locus2Probs$mendelian$W
        } else if(momAlleles[[2]][[j]]=="R"){
          fAllele[[2]][[j]] <- "R"
          fProbs[[2]][[j]] <- locus2Probs$mendelian$R
        }else if(momAlleles[[2]][[j]]=="B"){
          fAllele[[2]][[j]] <- "B"
          fProbs[[2]][[j]] <- locus2Probs$mendelian$B
        }#end locus 2 if statements
      }#end loop over alleles
    } else if(xor(fScore[[1]][1],fScore[[1]][2]) && xor(fScore[[2]][1],fScore[[2]][2])){
      #TF||FT and TF||FT
      #loop over 2 alleles at each locus
      for(j in 1:2){
        #at locus one, what is there, fill it

        if(momAlleles[[1]][[j]]=="W"){
          fAllele[[1]][[j]] <- c("W", "H", "E", "R", "B")
          fProbs[[1]][[j]] <- locus1Probs$homingMixed$W
        } else if(momAlleles[[1]][[j]]=="H"){
          fAllele[[1]][[j]] <- c("H", "E", "R", "B")
          fProbs[[1]][[j]] <- locus1Probs$homingMixed$H
        } else if(momAlleles[[1]][[j]]=="E"){
          fAllele[[1]][[j]] <- c("E", "R")
          fProbs[[1]][[j]] <- locus1Probs$homingMixed$E
        } else if(momAlleles[[1]][[j]]=="R"){
          fAllele[[1]][[j]] <- "R"
          fProbs[[1]][[j]] <- locus1Probs$homingMixed$R
        } else if(momAlleles[[1]][[j]]=="B"){
          fAllele[[1]][[j]] <- "B"
          fProbs[[1]][[j]] <- locus1Probs$homingMixed$B
        }#end locus 1 if statements

        if(momAlleles[[2]][[j]]=="W"){
          fAllele[[2]][[j]] <- c("W", "E", "R", "B")
          fProbs[[2]][[j]] <- locus2Probs$homing$W
        } else if(momAlleles[[2]][[j]]=="E"){
          fAllele[[2]][[j]] <- c("E", "R")
          fProbs[[2]][[j]] <- locus2Probs$homing$E
        } else if(momAlleles[[2]][[j]]=="R"){
          fAllele[[2]][[j]] <- "R"
          fProbs[[2]][[j]] <- locus2Probs$homing$R
        }else if(momAlleles[[2]][[j]]=="B"){
          fAllele[[2]][[j]] <- "B"
          fProbs[[2]][[j]] <- locus2Probs$homing$B
        }#end locus 2 if statements
      }#end loop over alleles
    } else if(xor(fScore[[1]][1],fScore[[1]][2]) && all(fScore[[2]])){
      #TF||FT and TT
      #loop over 2 alleles at each locus
      for(j in 1:2){
        #at locus one, what is there, fill it

        if(momAlleles[[1]][[j]]=="W"){
          fAllele[[1]][[j]] <- c("W", "H", "E", "R", "B")
          fProbs[[1]][[j]] <- locus1Probs$homingMixed$W
        } else if(momAlleles[[1]][[j]]=="H"){
          fAllele[[1]][[j]] <- c("H", "E", "R", "B")
          fProbs[[1]][[j]] <- locus1Probs$homingMixed$H
        } else if(momAlleles[[1]][[j]]=="E"){
          fAllele[[1]][[j]] <- c("E", "R")
          fProbs[[1]][[j]] <- locus1Probs$homingMixed$E
        } else if(momAlleles[[1]][[j]]=="R"){
          fAllele[[1]][[j]] <- "R"
          fProbs[[1]][[j]] <- locus1Probs$homingMixed$R
        }else if(momAlleles[[1]][[j]]=="B"){
          fAllele[[1]][[j]] <- "B"
          fProbs[[1]][[j]] <- locus1Probs$homingMixed$B
        }#end locus 1 if statements

        fAllele[[2]][[j]] <- c("E", "R")
        fProbs[[2]][[j]] <- locus2Probs$mendelian$E
      }#end loop over alleles
    } else if(all(fScore[[1]]) && all(!fScore[[2]])){
      #TT and FF
      #loop over 2 alleles at each locus
      for(j in 1:2){
        #at locus one, what is there, fill it

        fAllele[[1]][[j]] <- c("H", "R")
        fProbs[[1]][[j]] <- locus1Probs$mendelian$H

        if(momAlleles[[2]][[j]]=="W"){
          fAllele[[2]][[j]] <- c("W", "R")
          fProbs[[2]][[j]] <- locus2Probs$mendelian$W
        } else if(momAlleles[[2]][[j]]=="R"){
          fAllele[[2]][[j]] <- "R"
          fProbs[[2]][[j]] <- locus2Probs$mendelian$R
        }else if(momAlleles[[2]][[j]]=="B"){
          fAllele[[2]][[j]] <- "B"
          fProbs[[2]][[j]] <- locus2Probs$mendelian$B
        }#end locus 2 if statements
      }#end loop over alleles
    } else if(all(fScore[[1]]) && xor(fScore[[2]][1],fScore[[2]][2])){
      #TT and TF||FT
      #loop over 2 alleles at each locus
      for(j in 1:2){
        #at locus one, what is there, fill it

        fAllele[[1]][[j]] <- c("H", "E", "R", "B")
        fProbs[[1]][[j]] <- locus1Probs$homingTwo$H

        if(momAlleles[[2]][[j]]=="W"){
          fAllele[[2]][[j]] <- c("W", "E", "R", "B")
          fProbs[[2]][[j]] <- locus2Probs$homing$W
        } else if(momAlleles[[2]][[j]]=="E"){
          fAllele[[2]][[j]] <- c("E", "R")
          fProbs[[2]][[j]] <- locus2Probs$homing$E
        } else if(momAlleles[[2]][[j]]=="R"){
          fAllele[[2]][[j]] <- "R"
          fProbs[[2]][[j]] <- locus2Probs$homing$R
        }else if(momAlleles[[2]][[j]]=="B"){
          fAllele[[2]][[j]] <- "B"
          fProbs[[2]][[j]] <- locus2Probs$homing$B
        }#end locus 2 if statements
      }#end loop over alleles
    } else if(all(fScore[[1]]) && all(fScore[[2]])){
      #TT and TT
      #loop over 2 alleles at each locus
      for(j in 1:2){
        #at locus one, what is there, fill it

        fAllele[[1]][[j]] <- c("H", "E", "R", "B")
        fProbs[[1]][[j]] <- locus1Probs$homingTwo$H

        fAllele[[2]][[j]] <- c("E", "R")
        fProbs[[2]][[j]] <- locus2Probs$mendelian$E
      }#end loop over alleles
    }#end female checks


    ###########################################################################
    ## loop over male mate. This is the inner loop
    ###########################################################################
    for (mi in 1:fi){ #make this mi in 1 to fi
      #do male stuff here
      #split male genotype
      #This splits all characters.
      mSplit <- strsplit(x = genotypes[mi], split = "")[[1]]

      #make a list of each allele at every locus. This list is length nmPlex, and each
      # sublist has length 2
      dadAlleles <- list(mSplit[1:2], mSplit[3:4])

      #Score them
      mScore[[1]] <- grepl(pattern = "H", x = mSplit[1:2], fixed = TRUE)
      mScore[[2]] <- grepl(pattern = "E", x = mSplit[3:4], fixed = TRUE)

      #setup offspring allele lists
      mAllele <- rep(x = list(vector(mode = "list", 2)), numAlleles)
      mProbs <- rep(x = list(vector(mode = "list", 2)), numAlleles)



      #set father alleles
      if(all(!mScore[[1]])){
        #FF and FF||FT||TF||TT
        #loop over 2 alleles at each locus
        for(j in 1:2){
          #at locus one, what is there, fill it

          if(dadAlleles[[1]][[j]]=="W"){
            mAllele[[1]][[j]] <- c("W", "R")
            mProbs[[1]][[j]] <- locus1Probs$mendelian$W
          } else if(dadAlleles[[1]][[j]]=="E"){
            mAllele[[1]][[j]] <- c("E", "R")
            mProbs[[1]][[j]] <- locus1Probs$mendelian$E
          } else if(dadAlleles[[1]][[j]]=="R"){
            mAllele[[1]][[j]] <- "R"
            mProbs[[1]][[j]] <- locus1Probs$mendelian$R
          }else if(dadAlleles[[1]][[j]]=="B"){
            mAllele[[1]][[j]] <- "B"
            mProbs[[1]][[j]] <- locus1Probs$mendelian$B
          }#end locus 1 if statements

          if(dadAlleles[[2]][[j]]=="W"){
            mAllele[[2]][[j]] <- c("W", "R")
            mProbs[[2]][[j]] <- locus2Probs$mendelian$W
          } else if(dadAlleles[[2]][[j]]=="E"){
            mAllele[[2]][[j]] <- c("E", "R")
            mProbs[[2]][[j]] <- locus2Probs$mendelian$E
          } else if(dadAlleles[[2]][[j]]=="R"){
            mAllele[[2]][[j]] <- "R"
            mProbs[[2]][[j]] <- locus2Probs$mendelian$R
          }else if(dadAlleles[[2]][[j]]=="B"){
            mAllele[[2]][[j]] <- "B"
            mProbs[[2]][[j]] <- locus2Probs$mendelian$B
          }#end locus 2 if statements
        }#end loop over alleles
      } else if(xor(mScore[[1]][1],mScore[[1]][2]) && all(!mScore[[2]])){
        #TF||FT and FF
        #loop over 2 alleles at each locus
        for(j in 1:2){
          #at locus one, what is there, fill it

          if(dadAlleles[[1]][[j]]=="W"){
            mAllele[[1]][[j]] <- c("W", "H", "R", "B")
            mProbs[[1]][[j]] <- locus1Probs$homingOne$W
          } else if(dadAlleles[[1]][[j]]=="H"){
            mAllele[[1]][[j]] <- c("H","R")
            mProbs[[1]][[j]] <- locus1Probs$homingOne$H
          } else if(dadAlleles[[1]][[j]]=="E"){
            mAllele[[1]][[j]] <- c("E", "R")
            mProbs[[1]][[j]] <- locus1Probs$homingOne$E
          } else if(dadAlleles[[1]][[j]]=="R"){
            mAllele[[1]][[j]] <- "R"
            mProbs[[1]][[j]] <- locus1Probs$homingOne$R
          }else if(dadAlleles[[1]][[j]]=="B"){
            mAllele[[1]][[j]] <- "B"
            mProbs[[1]][[j]] <- locus1Probs$homingOne$B
          }#end locus 1 if statements

          if(dadAlleles[[2]][[j]]=="W"){
            mAllele[[2]][[j]] <- c("W", "R")
            mProbs[[2]][[j]] <- locus2Probs$mendelian$W
          } else if(dadAlleles[[2]][[j]]=="R"){
            mAllele[[2]][[j]] <- "R"
            mProbs[[2]][[j]] <- locus2Probs$mendelian$R
          }else if(dadAlleles[[2]][[j]]=="B"){
            mAllele[[2]][[j]] <- "B"
            mProbs[[2]][[j]] <- locus2Probs$mendelian$B
          }#end locus 2 if statements
        }#end loop over alleles
      } else if(xor(mScore[[1]][1], mScore[[1]][2]) && xor(mScore[[2]][1],mScore[[2]][2])){
        #TF||FT and TF||FT
        #loop over 2 alleles at each locus
        for(j in 1:2){
          #at locus one, what is there, fill it

          if(dadAlleles[[1]][[j]]=="W"){
            mAllele[[1]][[j]] <- c("W", "H", "E", "R", "B")
            mProbs[[1]][[j]] <- locus1Probs$homingMixed$W
          } else if(dadAlleles[[1]][[j]]=="H"){
            mAllele[[1]][[j]] <- c("H", "E", "R", "B")
            mProbs[[1]][[j]] <- locus1Probs$homingMixed$H
          } else if(dadAlleles[[1]][[j]]=="E"){
            mAllele[[1]][[j]] <- c("E", "R")
            mProbs[[1]][[j]] <- locus1Probs$homingMixed$E
          } else if(dadAlleles[[1]][[j]]=="R"){
            mAllele[[1]][[j]] <- "R"
            mProbs[[1]][[j]] <- locus1Probs$homingMixed$R
          } else if(dadAlleles[[1]][[j]]=="B"){
            mAllele[[1]][[j]] <- "B"
            mProbs[[1]][[j]] <- locus1Probs$homingMixed$B
          }#end locus 1 if statements

          if(dadAlleles[[2]][[j]]=="W"){
            mAllele[[2]][[j]] <- c("W", "E", "R", "B")
            mProbs[[2]][[j]] <- locus2Probs$homing$W
          } else if(dadAlleles[[2]][[j]]=="E"){
            mAllele[[2]][[j]] <- c("E", "R")
            mProbs[[2]][[j]] <- locus2Probs$homing$E
          } else if(dadAlleles[[2]][[j]]=="R"){
            mAllele[[2]][[j]] <- "R"
            mProbs[[2]][[j]] <- locus2Probs$homing$R
          }else if(dadAlleles[[2]][[j]]=="B"){
            mAllele[[2]][[j]] <- "B"
            mProbs[[2]][[j]] <- locus2Probs$homing$B
          }#end locus 2 if statements
        }#end loop over alleles
      } else if(xor(mScore[[1]][1], mScore[[1]][2]) && all(mScore[[2]])){
        #TF||FT and TT
        #loop over 2 alleles at each locus
        for(j in 1:2){
          #at locus one, what is there, fill it

          if(dadAlleles[[1]][[j]]=="W"){
            mAllele[[1]][[j]] <- c("W", "H", "E", "R", "B")
            mProbs[[1]][[j]] <- locus1Probs$homingMixed$W
          } else if(dadAlleles[[1]][[j]]=="H"){
            mAllele[[1]][[j]] <- c("H", "E", "R", "B")
            mProbs[[1]][[j]] <- locus1Probs$homingMixed$H
          } else if(dadAlleles[[1]][[j]]=="E"){
            mAllele[[1]][[j]] <- c("E", "R")
            mProbs[[1]][[j]] <- locus1Probs$homingMixed$E
          } else if(dadAlleles[[1]][[j]]=="R"){
            mAllele[[1]][[j]] <- "R"
            mProbs[[1]][[j]] <- locus1Probs$homingMixed$R
          }else if(dadAlleles[[1]][[j]]=="B"){
            mAllele[[1]][[j]] <- "B"
            mProbs[[1]][[j]] <- locus1Probs$homingMixed$B
          }#end locus 1 if statements

          mAllele[[2]][[j]] <- c("E", "R")
          mProbs[[2]][[j]] <- locus2Probs$mendelian$E
        }#end loop over alleles
      } else if(all(mScore[[1]]) && all(!mScore[[2]])){
        #TT and FF
        #loop over 2 alleles at each locus
        for(j in 1:2){
          #at locus one, what is there, fill it

          mAllele[[1]][[j]] <- c("H", "R")
          mProbs[[1]][[j]] <- locus1Probs$mendelian$H

          if(dadAlleles[[2]][[j]]=="W"){
            mAllele[[2]][[j]] <- c("W", "R")
            mProbs[[2]][[j]] <- locus2Probs$mendelian$W
          } else if(dadAlleles[[2]][[j]]=="R"){
            mAllele[[2]][[j]] <- "R"
            mProbs[[2]][[j]] <- locus2Probs$mendelian$R
          }else if(dadAlleles[[2]][[j]]=="B"){
            mAllele[[2]][[j]] <- "B"
            mProbs[[2]][[j]] <- locus2Probs$mendelian$B
          }#end locus 2 if statements
        }#end loop over alleles
      } else if(all(mScore[[1]]) && xor(mScore[[2]][1], mScore[[2]][2])){
        #TT and TF||FT
        #loop over 2 alleles at each locus
        for(j in 1:2){
          #at locus one, what is there, fill it

          mAllele[[1]][[j]] <- c("H", "E", "R", "B")
          mProbs[[1]][[j]] <- locus1Probs$homingTwo$H

          if(dadAlleles[[2]][[j]]=="W"){
            mAllele[[2]][[j]] <- c("W", "E", "R", "B")
            mProbs[[2]][[j]] <- locus2Probs$homing$W
          } else if(dadAlleles[[2]][[j]]=="E"){
            mAllele[[2]][[j]] <- c("E", "R")
            mProbs[[2]][[j]] <- locus2Probs$homing$E
          } else if(dadAlleles[[2]][[j]]=="R"){
            mAllele[[2]][[j]] <- "R"
            mProbs[[2]][[j]] <- locus2Probs$homing$R
          }else if(dadAlleles[[2]][[j]]=="B"){
            mAllele[[2]][[j]] <- "B"
            mProbs[[2]][[j]] <- locus2Probs$homing$B
          }#end locus 2 if statements
        }#end loop over alleles
      } else if(all(mScore[[1]]) && all(mScore[[2]])){
        #TT and TT
        #loop over 2 alleles at each locus
        for(j in 1:2){
          #at locus one, what is there, fill it

          mAllele[[1]][[j]] <- c("H", "E", "R", "B")
          mProbs[[1]][[j]] <- locus1Probs$homingTwo$H

          mAllele[[2]][[j]] <- c("E", "R")
          mProbs[[2]][[j]] <- locus2Probs$mendelian$E
        }#end loop over alleles
      }#end male checks


      #########################################################################
      ## Get combinations and put them in the tMatrix. This must be done
      ##  inside the inner loop
      #########################################################################

      #combine each locus into single lists, so that alleles within a locus can't
      # be combined with each other, but do get combined with allelels for the
      # other loci.
      # ie, unlist the sublists in each allele/probs list. This will give single-
      # depth lists the same length as the number of Daisy components.
      fAllLoci <- lapply(X = fAllele, FUN = unlist, recursive=TRUE)
      fProbsLoci <- lapply(X = fProbs, FUN = unlist, recursive=TRUE)
      mAllLoci <- lapply(X = mAllele, FUN = unlist, recursive=TRUE)
      mProbsLoci <- lapply(X = mProbs, FUN = unlist, recursive=TRUE)

      #combine male and female alleles at each locus.
      # This requires looping through each locus, getting all combinations of
      lociAList <- vector(mode = "list", length = numAlleles)
      lociPList <- vector(mode = "list", length = numAlleles)

      for( i in 1:numAlleles){

        #get all combinationes of male/female for each allele
        holdAllOne <- expand.grid(fAllLoci[[i]], mAllLoci[[i]], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        holdProbOne <- expand.grid(fProbsLoci[[i]], mProbsLoci[[i]], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

        #sort each combination so they are the same.
        holdAllOne <- apply(X = holdAllOne, MARGIN = 1, FUN = sort, method = 'radix')

        #paste alleles togheter
        holdAllTwo <- do.call(what = "paste0", list(holdAllOne[1, ], holdAllOne[2, ]))
        holdProbTwo <- holdProbOne[ ,1]*holdProbOne[ ,2]

        #aggregate and return
        aggregateHold <- vapply(X = unique(holdAllTwo), FUN = function(x){
          sum(holdProbTwo[holdAllTwo==x])},
          FUN.VALUE = numeric(length = 1L))

        #fill lists
        lociAList[[i]] <- names(aggregateHold)
        lociPList[[i]] <- aggregateHold

      }

      #get all combinations of each loci. This gives the total genotype
      outAList <- expand.grid(lociAList, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      outPList <- expand.grid(lociPList, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

      #combine allele names and probabilities
      outAList <- apply(X = outAList, MARGIN = 1, FUN = paste0, collapse="")
      outPList <- apply(X = outPList, MARGIN = 1, FUN = prod)
      #can use matrixStats::rowProds(x = as.matrix(outPList))

      #normalize
      outPList <- outPList/sum(outPList)

      #set values in tMatrix
      tMatrix[fi,mi, outAList ] <- outPList

    }# end male loop
  }# end female loop


  ## set the other half of the matrix
  # Boolean matrix for subsetting, used several times
  boolMat <- upper.tri(x = tMatrix[ , ,1], diag = FALSE)
  # loop over depth, set upper triangle
  for(z in 1:numGen){tMatrix[ , ,z][boolMat] <- t(tMatrix[ , ,z])[boolMat]}

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
    releaseType = "HHWW"
  ))

}
