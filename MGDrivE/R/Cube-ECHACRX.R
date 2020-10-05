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
#    Why isn't there female homing into H?
#   January 2020
#    Update implementation for a linked construct, ie, both target loci on on the X-chromosome
#    Added crossover ability
#    "E" targets W allele at first locus - H gets first shot, then E attacks, is
#     this right? Should it not target the first locus if alone?
#    Deposition - no secondary deposition in females. ie, no female deposition into her own alleles.
#     John is questioning this in other instances, I think we can't estimate it.
#    Male deposition - assume female with H and E only targets H on male first locus
#     and W on male second locus. This assume CAS9 is predominantly loaded with H gRNAs
#    Added back locus 1 homing parameter, Ethan's construct doesn't target the locus 1 wild-type
#   20200228
#    Ethan says that homing efficiency of H allele, in the presence of the E allele,
#    is vastly reduced. So, for the "homingMixed" case in females, we need new parameters
#    for H targeting W.
#
###############################################################################

#' Inheritance Cube: ECHACRX
#'
#' This function creates an X-linked ECHACR construct, it has 5 alleles at the first locus
#' and 4 alleles at the second.
#'  * W: Wild-type
#'  * H: Homing allele
#'  * E: Eraser allele
#'  * R: No-cost resistance allele
#'  * B: Detrimental resistance allele
#'  * cHW: Rate of homing from H, W -> H transition
#'  * cEH: Rate of homing from E, H -> E transition
#'  * cEW2: Rate of homing from E, W -> E transition
#'
#' This inheritance pattern corresponds to the [Active Genetic Neutralizing Elements for Halting or Deleting Gene Drives](https://doi.org/10.1016/j.molcel.2020.09.003) publication.
#'
#' @param cHW Cutting efficiency of drive allele at locus 1
#' @param cEHW Cutting efficiency of drive allele, in the presence of ECHACR element, at locus 1
#' @param cEW1 Cutting efficiency of ECHACR element into W at locus 1
#' @param cEW2 Cutting efficiency of ECHACR element into W at locus 2
#' @param cEH Cutting efficiency of ECHACR element into H
#' @param chHW Homing efficiency of drive allele at locus 1
#' @param crHW Resistance allele efficiency of drive allele at locus 1
#' @param chEHW Homing efficiency of drive allele, in the presence of ECHACR element, at locus 1
#' @param crEHW Resistance allele efficiency of drive allele, in the presence of ECHACR element, at locus 1
#' @param ceEW1 Homing efficiency of ECHACR element into W at locus 1
#' @param crEW1 Resistance allele efficiency of ECHACR element into W at locus 1
#' @param ceEW2 Homing efficiency of ECHACR element into W at locus 2
#' @param crEW2 Resistance allele efficiency of ECHACR element into W at locus 2
#' @param ceEH Homing efficiency of ECHACR element into H
#' @param crEH Resistance allele efficiency of ECHACR element into H
#' @param d1 Background mutation rate from W into R allele
#' @param d2 Background mutation rate from H into R allele
#' @param d3 Background mutation rate from E into R allele
#' @param dHW Female H deposition rate against W
#' @param dEH Female E deposition rate against H
#' @param dEW Female E deposition rate against W
#' @param drHW Female resistance generation rate, from H allele
#' @param drEH Female resistance generation rate, from E allele
#' @param drEW Female resistance generation rate, from E allele
#' @param crossF Female crossover rate. 0 is fully linked, 0.5 is unlinked, 1 is negatively linked
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
cubeECHACRX <- function(cHW=1.0, cEHW=1.0, cEW1=1.0, cEW2=1.0, cEH=1.0,
                        chHW=0, crHW=0, chEHW=0, crEHW=0, ceEW1=0, crEW1=0,
                        ceEW2=0, crEW2=0, ceEH=0, crEH=0,
                        d1=0, d2=0, d3=0,
                        dHW=0, dEH=0, dEW=0, drHW=0, drEH=0, drEW=0, crossF=0,
                        eta=NULL, phi=NULL,omega=NULL, xiF=NULL, xiM=NULL, s=NULL){


  # cHW=1.0; cEW2=1.0; cEH=1.0; chHW=0; crHW=0; ceEW2=0; crEW2=0;
  # ceEH=0; crEH=0; d1=0; d2=0; d3=0; dHW=0; dEH=0; dEW=0;
  # drHW=0; drEH=0; drEW=0; crossF=0;



  ## safety checks
  if(any(c(cHW,cEW1,cEW2,cEH,chHW,crHW,ceEW1,crEW1,ceEW2,crEW2,ceEH,crEH,d1,d2,d3,dHW,dEH,dEW,drHW,drEH,drEW)>1)
     || any(c(cHW,cEW1,cEW2,cEH,chHW,crHW,ceEW1,crEW1,ceEW2,crEW2,ceEH,crEH,d1,d2,d3,dHW,dEH,dEW,drHW,drEH,drEW)<0)){
    stop("Parameters are rates must be between 0 and 1.")
  }


  #############################################################################
  ## generate all genotypes, set up vectors and matrices
  #############################################################################
  # I need more clarification on how this is setup, so writing stuff here
  #  Per Ethan, they actually created an X-linked drive, with locus 1 and locus 2
  #  about 1cM apart. This means they are effectively linked, but will crossover
  #  at ~1% per generation. So, this drive was rebuilt (read, totally reorganized in genotype space)
  #  to accomodate this.
  # Each genotype is 4 characters long. The first 2 character are one chromosome, the
  #  second 2 characters are the other chromosome. Consider this as an "XX"/"XY" system,
  #  where there are 2 loci of interest on the X chromosome. Thus, 2 chromosomes (X and Y, or X and X)
  #  and 2 loci per X (the Y will just be YY, because we don't care about stuff on it),
  #  thus the genotypes are length 4 (2 times 2)
  # In Ethan's current drive, the Homing construct has a Cas and gRNAs, while the
  #  eracing construct has 2 gRNAs. The homing construct exists at locus 1 and targets locus 1 on the
  #  X-chromosome, and the eracing construct exists at locus 2 and targets locus 1 and locus 2. Once
  #  the homing construct has been removed/damaged, there is no homing, even if
  #  the E construct is at locus 1.
  # IF the E construct is at locus 2, the gRNAs are present, and there is the possibility
  #  of homing, if the H construct is at locus 1.
  # IF the E construct is at locus 1, in this case, it actually implies a damanged H
  #  at that locus, so there are no E gRNAs and no CAS9, but it keeps the fluorescent
  #  marker of the H allele


  # run this stuff once, then set the variables for later.

  #list of possible alleles at each locus
  #the first locus has the drive, and can be erased
  # the second locus is the "CHACR" element. No drive, but eracing piece
  gTypes <- list(c("W", "H", "E", "R", "B"), c("W", "E", "R", "B"))

  # # Since locus 1 and locus 2 are potentially linked, generate each combination of
  # #  them at one place
  # alleles <- expand.grid(gTypes,KEEP.OUT.ATTRS = FALSE,
  #                        stringsAsFactors = FALSE)
  # # paste them together. Add male allele. Only occurs in YY format
  # alleles <- do.call(what = paste0, args = list(alleles[,1], alleles[,2]))
  # alleles <- c(alleles, "YY")
  #
  # # expand all combinations of alleles with locus 1 and locus 2
  # #  This provides our 2 copies, since diploids have are XX or XY, and each
  # #  X has locus 1 and locus 2
  # hold <- expand.grid(alleles,alleles,KEEP.OUT.ATTRS = FALSE,
  #                     stringsAsFactors = FALSE)
  # # sort, because order of alleles doesn't matter, paste, and keep unique
  # #  order of loci does matter, we do not sort that!
  # genotypes <- unique(vapply(X = 1:dim(hold)[1],
  #                            FUN = function(x){
  #                              paste0(sort(x = hold[x, ]), collapse = "")},
  #                            FUN.VALUE = character(1)
  # ))
  # # remove YYYY, not viable
  # genotypes <- genotypes[!genotypes=="YYYY"]
  #
  # #separate male/female genotypes, since this is sex specific now
  # index <- grepl(pattern = "YY", x = genotypes)
  # maleGen <- genotypes[index]
  # femaleGen <- genotypes[!index]

  femaleGen <- c("WWWW","HWWW","EWWW","RWWW","BWWW","WEWW","HEWW","EEWW","REWW",
                 "BEWW","WRWW","HRWW","ERWW","RRWW","BRWW","WBWW","HBWW","EBWW",
                 "RBWW","BBWW","HWHW","EWHW","HWRW","BWHW","HWWE","HEHW","EEHW",
                 "HWRE","BEHW","HWWR","HRHW","ERHW","HWRR","BRHW","HWWB","HBHW",
                 "EBHW","HWRB","BBHW","EWEW","EWRW","BWEW","EWWE","EWHE","EEEW",
                 "EWRE","BEEW","EWWR","EWHR","EREW","EWRR","BREW","EWWB","EWHB",
                 "EBEW","EWRB","BBEW","RWRW","BWRW","RWWE","HERW","EERW","RERW",
                 "BERW","RWWR","HRRW","ERRW","RRRW","BRRW","RWWB","HBRW","EBRW",
                 "RBRW","BBRW","BWBW","BWWE","BWHE","BWEE","BWRE","BEBW","BWWR",
                 "BWHR","BWER","BWRR","BRBW","BWWB","BWHB","BWEB","BWRB","BBBW",
                 "WEWE","HEWE","EEWE","REWE","BEWE","WEWR","HRWE","ERWE","RRWE",
                 "BRWE","WBWE","HBWE","EBWE","RBWE","BBWE","HEHE","EEHE","HERE",
                 "BEHE","HEWR","HEHR","ERHE","HERR","BRHE","HEWB","HBHE","EBHE",
                 "HERB","BBHE","EEEE","EERE","BEEE","EEWR","EEHR","EEER","EERR",
                 "BREE","EEWB","EEHB","EBEE","EERB","BBEE","RERE","BERE","REWR",
                 "HRRE","ERRE","RERR","BRRE","REWB","HBRE","EBRE","RBRE","BBRE",
                 "BEBE","BEWR","BEHR","BEER","BERR","BEBR","BEWB","BEHB","BEEB",
                 "BERB","BBBE","WRWR","HRWR","ERWR","RRWR","BRWR","WBWR","HBWR",
                 "EBWR","RBWR","BBWR","HRHR","ERHR","HRRR","BRHR","HRWB","HBHR",
                 "EBHR","HRRB","BBHR","ERER","ERRR","BRER","ERWB","ERHB","EBER",
                 "ERRB","BBER","RRRR","BRRR","RRWB","HBRR","EBRR","RBRR","BBRR",
                 "BRBR","BRWB","BRHB","BREB","BRRB","BBBR","WBWB","HBWB","EBWB",
                 "RBWB","BBWB","HBHB","EBHB","HBRB","BBHB","EBEB","EBRB","BBEB",
                 "RBRB","BBRB","BBBB")

  maleGen <- c("WWYY","HWYY","EWYY","RWYY","BWYY","WEYY","HEYY","EEYY","REYY","BEYY",
               "WRYY","HRYY","ERYY","RRYY","BRYY","WBYY","HBYY","EBYY","RBYY","BBYY")

  #############################################################################
  ## setup all probability lists
  #############################################################################

  ## female
  femaleLocus1 = femaleLocus2 = list()

  # locus 1, homing and eracing
  femaleLocus1$mendelian <-list("W"=c("W"=1-d1,"R"=d1),
                                "H"=c("H"=1-d2,"R"=d2),
                                "E"=c("E"=1-d3,"R"=d3),
                                "R"=c("R"=1),
                                "B"=c("B"=1))

  femaleLocus1$homing1 <- list("W"=c("W"=(1-d1)*(1-cHW),
                                     "H"=(1-d1)*cHW*chHW,
                                     "R"=d1 + (1-d1)*cHW*(1-chHW)*crHW,
                                     "B"=(1-d1)*cHW*(1-chHW)*(1-crHW)),
                               "H"=c("H"=1-d2,"R"=d2),
                               "E"=c("E"=1-d3,"R"=d3),
                               "R"=c("R"=1),
                               "B"=c("B"=1))

  femaleLocus1$homingMixed <- list("W"=c("W"=(1-d1)*(1-cEHW)*(1-cEW1),
                                         "H"=(1-d1)*cEHW*chEHW*(1-cEH),
                                         "E"=(1-d1)*(1-cEHW)*cEW1*ceEW1 + (1-d1)*cEHW*chEHW*cEH*ceEH,
                                         "R"=d1 + (1-d1)*cEHW*(1-chEHW)*crEHW + (1-d1)*(1-cEHW)*cEW1*(1-ceEW1)*crEW1 + (1-d1)*cEHW*chEHW*cEH*(1-ceEH)*crEH,
                                         "B"=(1-d1)*cEHW*(1-chEHW)*(1-crEHW) + (1-d1)*(1-cEHW)*cEW1*(1-ceEW1)*(1-crEW1) + (1-d1)*cEHW*chEHW*cEH*(1-ceEH)*(1-crEH)),
                                   "H"=c("H"=(1-d2)*(1-cEH),
                                         "E"=(1-d2)*cEH*ceEH,
                                         "R"=d2 + (1-d2)*cEH*(1-ceEH)*crEH,
                                         "B"=(1-d2)*cEH*(1-ceEH)*(1-crEH)),
                                   "E"=c("E"=1-d3,"R"=d3),
                                   "R"=c("R"=1),
                                   "B"=c("B"=1))

  femaleLocus1$homing2 <- list("W"=c("W"=1-d1,"R"=d1),
                               "H"=c("H"=(1-d2)*(1-cEH),
                                     "E"=(1-d2)*cEH*ceEH,
                                     "R"=d2 + (1-d2)*cEH*(1-ceEH)*crEH,
                                     "B"=(1-d2)*cEH*(1-ceEH)*(1-crEH)),
                               "E"=c("E"=1-d3,"R"=d3),
                               "R"=c("R"=1),
                               "B"=c("B"=1))

  # locus 2, CHACR element
  femaleLocus2$mendelian <- list("W"=c("W"=1-d1,"R"=d1),
                                 "E"=c("E"=1-d3,"R"=d3),
                                 "R"=c("R"=1),
                                 "B"=c("B"=1))

  femaleLocus2$homing <- list("W"=c("W"=(1-d1)*(1-cEW2),
                                    "E"=(1-d1)*cEW2*ceEW2,
                                    "R"=d1 + (1-d1)*cEW2*(1-ceEW2)*crEW2,
                                    "B"=(1-d1)*cEW2*(1-ceEW2)*(1-crEW2)),
                              "E"=c("E"=1-d3,"R"=d3),
                              "R"=c("R"=1),
                              "B"=c("B"=1))

  ## male
  maleLocus1 = maleLocus2 = list()

  # locus 1
  maleLocus1$mendelian <- list("W"=c("W"=1-d1,"R"=d1),
                               "H"=c("H"=1-d2,"R"=d2),
                               "E"=c("E"=1-d3,"R"=d3),
                               "R"=c("R"=1),
                               "B"=c("B"=1))

  maleLocus1$homing1 <- list("W"=c("W"=(1-d1)*(1-dHW),
                                   "R"=d1 + (1-d1)*dHW*drHW,
                                   "B"=(1-d1)*dHW*(1-drHW)),
                             "H"=c("H"=1-d2,"R"=d2),
                             "E"=c("E"=1-d3,"R"=d3),
                             "R"=c("R"=1),
                             "B"=c("B"=1))

  maleLocus1$homing2 <- list("W"=c("W"=1-d1,"R"=d1),
                             "H"=c("H"=(1-d2)*(1-dEH),
                                   "R"=d2 + (1-d2)*dEH*drEH,
                                   "B"=(1-d2)*dEH*(1-drEH)),
                             "E"=c("E"=1-d3,"R"=d3),
                             "R"=c("R"=1),
                             "B"=c("B"=1))

  # locus 2
  maleLocus2$mendelian <- list("W"=c("W"=1-d1,"R"=d1),
                               "E"=c("E"=1-d3,"R"=d3),
                               "R"=c("R"=1),
                               "B"=c("B"=1))

  maleLocus2$homing <- list("W"=c("W"=(1-d1)*(1-dEW),
                                  "R"=d1 + (1-d1)*dEW*drEW,
                                  "B"=(1-d1)*dEW*(1-drEW)),
                            "E"=c("E"=1-d3,"R"=d3),
                            "R"=c("R"=1),
                            "B"=c("B"=1))


  #############################################################################
  ## fill transition matrix
  #############################################################################
  #use this many times down below
  numGen <- length(femaleGen) + length(maleGen)

  #create transition matrix to fill
  tMatrix <- array(data = 0,
                   dim = c(numGen,numGen,numGen),
                   dimnames = list(c(femaleGen,maleGen),c(femaleGen,maleGen),c(femaleGen,maleGen)))

  #number of alleles, set score lists
  numAlleles <- 2
  fScore <- mScore <- vector(mode = "list", length = numAlleles)


  #############################################################################
  ## loop over all matings, female outer loop
  #############################################################################
  for (fi in 1:length(femaleGen)){
    #do female stuff here
    #This splits all characters.
    fSplit <- strsplit(x = femaleGen[fi], split = "")[[1]]
    #make a list of each allele at every locus. This list is length nmPlex, and each
    # sublist has length 2
    momAlleles <- list(fSplit[1:2], fSplit[3:4])
    #Score them
    # H allele implies a CAS9 and gRNAs at locus 1
    # E allele implies 2 gRNAs at locus 2
    fScore[[1]] <- grepl(pattern = "H", x = fSplit[c(1,3)], fixed = TRUE)
    fScore[[2]] <- grepl(pattern = "E", x = fSplit[c(2,4)], fixed = TRUE)
    #setup offspring allele lists
    fAllele <- rep(x = list(vector(mode = "list", 2)), numAlleles)
    fProbs <- rep(x = list(vector(mode = "list", 2)), numAlleles)


    ####################
    ## female alleles
    ####################
    if(all(!fScore[[1]])){
      #FF and FF||FT||TF||TT
      #loop over 2 alleles at each locus
      for(j in 1:numAlleles){

        # locus 1
        fProbs[[j]][[1]] <- femaleLocus1$mendelian[[ momAlleles[[j]][[1]] ]]
        fAllele[[j]][[1]] <- names(fProbs[[j]][[1]])

        # locus 2
        fProbs[[j]][[2]] <- femaleLocus2$mendelian[[ momAlleles[[j]][[2]] ]]
        fAllele[[j]][[2]] <- names(fProbs[[j]][[2]])

      }#end loop over alleles
    } else if(xor(fScore[[1]][1],fScore[[1]][2]) && all(!fScore[[2]])){
      #TF||FT and FF
      #loop over 2 alleles at each locus
      for(j in 1:numAlleles){

        # locus 1
        fProbs[[j]][[1]] <- femaleLocus1$homing1[[ momAlleles[[j]][[1]] ]]
        fAllele[[j]][[1]] <- names(fProbs[[j]][[1]])

        # locus 2
        fProbs[[j]][[2]] <- femaleLocus2$mendelian[[ momAlleles[[j]][[2]] ]]
        fAllele[[j]][[2]] <- names(fProbs[[j]][[2]])

      }#end loop over alleles
    } else if(xor(fScore[[1]][1],fScore[[1]][2]) && any(fScore[[2]])){
      #TF||FT and TF||FT||TT
      #loop over 2 alleles at each locus
      for(j in 1:numAlleles){

        # locus 1
        fProbs[[j]][[1]] <- femaleLocus1$homingMixed[[ momAlleles[[j]][[1]] ]]
        fAllele[[j]][[1]] <- names(fProbs[[j]][[1]])

        # locus 2
        fProbs[[j]][[2]] <- femaleLocus2$homing[[ momAlleles[[j]][[2]] ]]
        fAllele[[j]][[2]] <- names(fProbs[[j]][[2]])

      }#end loop over alleles
    } else if(all(fScore[[1]]) && all(!fScore[[2]])){
      #TT and FF
      #loop over 2 alleles at each locus
      for(j in 1:numAlleles){

        # locus 1
        fProbs[[j]][[1]] <- femaleLocus1$mendelian[[ momAlleles[[j]][[1]] ]]
        fAllele[[j]][[1]] <- names(fProbs[[j]][[1]])

        # locus 2
        fProbs[[j]][[2]] <- femaleLocus2$mendelian[[ momAlleles[[j]][[2]] ]]
        fAllele[[j]][[2]] <- names(fProbs[[j]][[2]])

      }#end loop over alleles
    } else if(all(fScore[[1]]) && any(fScore[[2]])){
      #TT and TF||FT||TT
      #loop over 2 alleles at each locus
      for(j in 1:numAlleles){

        # locus 1
        fProbs[[j]][[1]] <- femaleLocus1$homing2[[ momAlleles[[j]][[1]] ]]
        fAllele[[j]][[1]] <- names(fProbs[[j]][[1]])

        # locus 2
        fProbs[[j]][[2]] <- femaleLocus2$homing[[ momAlleles[[j]][[2]] ]]
        fAllele[[j]][[2]] <- names(fProbs[[j]][[2]])

      }#end loop over alleles
    }#end female checks


    ####################
    ## female crossover
    ####################
    #because this is done recursively, I need to hold the value out
    holdA <- fAllele[[1]][[2]]
    holdP <- fProbs[[1]][[2]]

    #swap allele 2's in both
    fAllele[[1]][[2]] <- c(holdA, fAllele[[2]][[2]])
    fAllele[[2]][[2]] <- c(fAllele[[2]][[2]], holdA)

    #set the probs for the new allele 2's in both
    fProbs[[1]][[2]] <- c(holdP*(1-crossF), fProbs[[2]][[2]]*crossF)
    fProbs[[2]][[2]] <- c(fProbs[[2]][[2]]*(1-crossF), holdP*crossF)


    ####################
    ## combine female allele
    ####################
    # This requires looping through each locus, getting all combinations of
    fLociA <- fLociP <- vector(mode = "list", length = numAlleles)

    for( j in 1:numAlleles){

      #get all combinations of loci for each allele, keep male/female separate still
      holdAllF <- expand.grid(fAllele[[j]], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      holdProbF <- expand.grid(fProbs[[j]], KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

      #paste alleles together and store in list
      fLociA[[j]] <- do.call(what = "paste0", list(holdAllF[ ,1], holdAllF[ ,2]))
      fLociP[[j]] <- holdProbF[ ,1]*holdProbF[ ,2]
    }

    # unlist
    #  All female alleles are now combined with their proper loci, with whatever
    #  probability of crossover there is. Now, each of these must pair with a male
    #  chromosome for offspring, thus unlist so we don't pair with themselves later
    fAllLoci <- unlist(fLociA)
    fProbsLoci <- unlist(fLociP)


    ###########################################################################
    ## loop over male mate. This is the inner loop
    ###########################################################################
    for (mi in 1:length(maleGen)){ #make this mi in 1 to fi
      #do male stuff here
      #split male genotype
      #This splits all characters.
      mSplit <- strsplit(x = maleGen[mi], split = "")[[1]]

      #setup offspring allele lists
      mAllele <- vector(mode = "list", 2)
      mProbs <- vector(mode = "list", 2)


      #set father alleles
      if(all(!fScore[[1]])){
        # female doesn't carry H, no deposition, just mendelian in father alleles

        # locus 1
        mProbs[[1]] <- c(maleLocus1$mendelian[[ mSplit[1] ]])
        mAllele[[1]] <- names(mProbs[[1]])

        # locus 2
        mProbs[[2]] <- c(maleLocus2$mendelian[[ mSplit[2] ]])
        mAllele[[2]] <- names(mProbs[[2]])

      } else if(any(fScore[[1]]) && all(!fScore[[2]]) ){
        # female has at least 1 H allele, but no E at second locus

        # locus 1
        mProbs[[1]] <- c(maleLocus1$homing1[[ mSplit[1] ]])
        mAllele[[1]] <- names(mProbs[[1]])

        # locus 2
        mProbs[[2]] <- c(maleLocus2$mendelian[[ mSplit[2] ]])
        mAllele[[2]] <- names(mProbs[[2]])

      } else if(any(fScore[[1]]) && any(fScore[[2]]) ){
        # female has at least 1 H at first locus and at least 1 E at second

        # locus 1
        mProbs[[1]] <- c(maleLocus1$homing2[[ mSplit[1] ]])
        mAllele[[1]] <- names(mProbs[[1]])

        # locus 2
        mProbs[[2]] <- c(maleLocus2$homing[[ mSplit[2] ]])
        mAllele[[2]] <- names(mProbs[[2]])

      }#end male checks


      ####################
      ## male crossover
      ####################
      # There is none


      ####################
      ## combine male allele
      ####################
      mAllLoci <- expand.grid(mAllele, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      mAllLoci <- c(do.call(what = "paste0", list(mAllLoci[ ,1], mAllLoci[ ,2])), "YY")

      mProbsLoci <- expand.grid(mProbs, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      mProbsLoci <- c(mProbsLoci[ ,1] * mProbsLoci[ ,2], 1)


      #########################################################################
      ## Get combinations and put them in the tMatrix. This must be done
      ##  inside the inner loop
      #########################################################################
      # male and female alleles/probs are already unlisted, so we have a vector of
      #  male and a vector of female alleles, with the loci already properly mixed.
      # Here, combine male with female, reduce by same genotype, then put into matrix

      ####################
      ## Combinations
      ####################
      # 1 - all combinations
      # 2 - sort, since order of chromosome doesn't matter
      # 3 - combine into genotypes
      allGens <- expand.grid(fAllLoci, mAllLoci, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      allGens <- apply(X = allGens, MARGIN = 1, FUN = sort.int, method = 'radix')
      allGens <- do.call(what = file.path, list(allGens[1, ], allGens[2, ], fsep=""))

      # 1 - all probability combinations
      # 2 - mutliply together
      allProbs <- expand.grid(fProbsLoci, mProbsLoci, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
      allProbs <- allProbs[ ,1]*allProbs[ ,2]


      ####################
      ## Aggregate/Normalize
      ####################
      # aggregate over duplicated genotypes
      aggregateHold <- vapply(X = unique(allGens), FUN = function(x){
        sum(allProbs[allGens==x])},
        FUN.VALUE = numeric(length = 1L))

      # normalize to a probability
      aggregateHold <- aggregateHold/sum(aggregateHold)


      ####################
      ## cube
      ####################
      #set values in tMatrix
      tMatrix[femaleGen[fi],maleGen[mi], names(aggregateHold) ] <- aggregateHold


    }# end male loop
  }# end female loop


  ## set the other half of the matrix
  tMatrix[tMatrix < .Machine$double.eps] <- 0 #protection from underflow errors

  ## initialize viability mask. No mother-specific death.
  viabilityMask <- array(data = 1, dim = c(numGen,numGen,numGen),
                         dimnames = list(c(femaleGen,maleGen), c(femaleGen,maleGen), c(femaleGen,maleGen)))

  ## genotype-specific modifiers
  phi = setNames(object = c(rep.int(x = 1, times = length(femaleGen)), rep.int(x = 0, times = length(maleGen))), nm = c(femaleGen, maleGen))
  modifiers = cubeModifiers(c(femaleGen,maleGen), eta = eta, phi = phi, omega = omega, xiF = xiF, xiM = xiM, s = s)


  ## put everything into a labeled list to return
  return(list(
    ih = tMatrix,
    tau = viabilityMask,
    genotypesID = c(femaleGen,maleGen),
    genotypesN = numGen,
    wildType = c("WWWW","WWYY"),
    eta = modifiers$eta,
    phi = modifiers$phi,
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = "HHYY"
  ))

}
