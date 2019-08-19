###############################################################################
#        ____        __       __
#       / __ \____ _/ /______/ /_
#      / /_/ / __ `/ __/ ___/ __ \
#     / ____/ /_/ / /_/ /__/ / / /
#    /_/    \__,_/\__/\___/_/ /_/
#
#   Patch Class Mosquito Population Simulation
#   Marshall Lab
#   November 2017
#
###############################################################################

###############################################################################
# Daily Simulation
###############################################################################

#' Daily Population Dynamics for a Patch
#'
#' Run population dynamics (not including migration) for this patch.
#'
oneDay_PopDynamics_Patch <- function(){

  ################
  #aquatic stages#
  ################

  # eggs
  self$oneDay_ovipositG1()

  # larvae
  self$oneDay_calcLarvalDensityDependentFactor() # calculate larDDMortal
  self$oneDay_larSurvival() # survival
  self$oneDay_hatchingFract() # hatching
  self$oneDay_eggReleases2eggs() # egg release to eggs
  self$oneDay_larHatching()

  # pupating
  self$oneDay_eggsFract2()
  self$oneDay_calcCumulativeLarvalDensityDependentFactor() # calculate f
  self$oneDay_larPupating()

  # balance
  private$LARnew = private$larSurvival + private$larHatching - private$larPupating

  # bastard fix
  private$LARnew[private$LARnew < 0]  = 0L

  # egg releases into adult pop
  self$oneDay_eggReleases2adults()

  #######################
  #adult male population#
  #######################

  self$oneDay_numMaleFemale()

  # survival
  self$oneDay_admSurvival()
  self$oneDay_calcCumulativePupaDensityDependentFactor() # calculate f
  self$oneDay_admPupating()

  self$oneDay_maleReleases()

  #########################
  #adult female population#
  #########################

  # first gonotrophic cycle
  self$oneDay_af1Survival()
  self$oneDay_af1Pupation()

  self$oneDay_femaleReleases()

  self$oneDay_af1Mating()

}

Patch$set(which = "public",name = "oneDay_PopDynamics",
          value = oneDay_PopDynamics_Patch, overwrite = TRUE
)




###############################################################################
# Aquatic Stage Dynamics
###############################################################################

###############################################################################
# ovipositG1
###############################################################################

#' Stochastic Oviposition
#'
#' Calculate the number of eggs oviposited by female mosquitoes following:
#' \deqn{\overline{O(T_x)} = \sum_{j=1}^{n} \Bigg( \bigg( (\beta*\overline{s} * \overline{ \overline{Af_{[t-T_x]}}}) * \overline{\overline{\overline{Ih}}} \bigg) * \Lambda  \Bigg)^{\top}_{ij}}
#' The deterministic result for number of eggs is used as the mean of a Poisson-distributed number of actual eggs oviposited.
#'
oneDay_ovipositG1_stochastic_Patch <- function(){

  #used multiple times/get more complicated in other uses
  genotypesN = private$NetworkPointer$get_genotypesN()
  timeIndex = private$NetworkPointer$get_timeAq()

  probsHold <- private$AF1dly[[timeIndex]] *
    private$NetworkPointer$get_beta()*
    private$NetworkPointer$get_s()

  for(depthInd in 1:genotypesN){
    probs <- probsHold * private$NetworkPointer$get_drivecubeindex(NULL,NULL,depthInd) *
      private$NetworkPointer$get_tau(NULL,NULL,depthInd)

    private$ovipositG1[depthInd] <- sum(rpois(n = genotypesN*genotypesN, lambda = probs))
  }

}

#' Deterministc Oviposition
#'
#' Calculate the number of eggs oviposited by female mosquitoes following:
#' \deqn{\overline{O(T_x)} = \sum_{j=1}^{n} \Bigg( \bigg( (\beta*\overline{s} * \overline{ \overline{Af_{[t-T_x]}}}) * \overline{\overline{\overline{Ih}}} \bigg) * \Lambda  \Bigg)^{\top}_{ij}}
#'
oneDay_ovipositG1_deterministic_Patch <- function(){

  genotypesN = private$NetworkPointer$get_genotypesN()
  timeIndex = private$NetworkPointer$get_timeAq()

  #fill offspring cube with parents
  for(slice in 1:genotypesN){

    private$ovipositG1[slice] = sum(private$AF1dly[[timeIndex]] *
                                      private$NetworkPointer$get_beta() *
                                      private$NetworkPointer$get_s() *
                                      private$NetworkPointer$get_drivecubeindex(NULL,NULL,slice) *
                                      private$NetworkPointer$get_tau(NULL,NULL,slice)
                                    )

  }

}

###############################################################################
# larSurvival
###############################################################################

#' Stochastic Larval Survival
#'
#' Calculate the number of larvae surviving from day to day, given by:
#' \deqn{\overline{L_{[t-1]}} * (1-\mu_{l}) * F(\overline{L_{[t-1]})}}
#' Stochasticity is introduced by assuming \eqn{(1-\mu_{l}) * F(\overline{L_{[t-1]})}} defines binomial likelihood of survival for each
#' genotype of larvae.
#'
oneDay_larSurvival_stochastic_Patch <- function(){

  size = length(private$LARdly[[1]])
  multipliedProbs = (1-private$NetworkPointer$get_muAq()) * private$larDDMortal
  genotypes = private$NetworkPointer$get_genotypesID()

  counter = which(private$LARdly[[1]] != 0)

  private$larSurvival = matrix(data=0L,nrow=2,ncol=size)
  if (length(multipliedProbs)!= size) { multipliedProbs = rep.int( x = multipliedProbs, times = size ) }

  fProbs = cbind(multipliedProbs, 1-multipliedProbs, deparse.level = 0)

  holder = vapply( counter , function(x){rmultinom(n = 1, size = private$LARdly[[1]][x], prob = fProbs[x,])}, FUN.VALUE = integer(2) )

  private$larSurvival[,counter] = holder

  private$larSurvival = private$larSurvival[1,]
}

#' Deterministc Larval Survival
#'
#' Calculate the number of larvae surviving from day to day, given by:
#' \deqn{\overline{L_{[t-1]}} * (1-\mu_{l}) * F(\overline{L_{[t-1]})}}
#'
oneDay_larSurvival_deterministic_Patch <- function(){
  private$larSurvival = private$LARdly[[1]] * (1-private$NetworkPointer$get_muAq()) * private$larDDMortal
}

###############################################################################
# hatchingFract
###############################################################################

#' Stochastic Fraction of Eggs Maturing to Hatch
#'
#' Calculate the fraction of hatching eggs accounting for delay given by: \eqn{\overline{O(T_e)}}
#' Stochasticity is introduced by assuming \eqn{\overline{O(T_e)}} defines the mean of a Poisson distributed number of eggs that can hatch.
#'
oneDay_hatchingFract_stochastic_Patch <- function(){

  #used multiple times/get more complicated in other uses
  genotypesN = private$NetworkPointer$get_genotypesN()
  timeIndex = private$NetworkPointer$get_timeAq("E")

  probsHold <- private$AF1dly[[timeIndex]] *
    private$NetworkPointer$get_beta()*
    private$NetworkPointer$get_s()

  for(depthInd in 1:genotypesN){
    probs <- probsHold * private$NetworkPointer$get_drivecubeindex(NULL,NULL,depthInd) *
      private$NetworkPointer$get_tau(NULL,NULL,depthInd)

    private$hatchingFract[depthInd] <- sum(rpois(n = genotypesN*genotypesN, lambda = probs))

  }

}

#' Deterministc Fraction of Eggs Maturing to Hatch
#'
#' Calculate the fraction of hatching eggs accounting for delay given by: \eqn{\overline{O(T_e)}}
#'
oneDay_hatchingFract_deterministic_Patch <- function(){

  genotypesN = private$NetworkPointer$get_genotypesN()
  timeIndex = private$NetworkPointer$get_timeAq("E")


  #fill offspring cube with parents
  for(slice in 1:genotypesN){

    private$hatchingFract[slice] = sum(private$AF1dly[[timeIndex]] *
                                      private$NetworkPointer$get_beta() *
                                      private$NetworkPointer$get_s() *
                                      private$NetworkPointer$get_drivecubeindex(NULL,NULL,slice) *
                                      private$NetworkPointer$get_tau(NULL,NULL,slice)
    )

  }

}

###############################################################################
# larHatching
###############################################################################

#' Stochastic Egg Hatching to Larval Stage
#'
#' Calculate the eggs that have survived and hatched during a day given by: \eqn{\overline{O(T_e)}* \theta_{e}}
#' The number of eggs that survive to hatch follows a binomial distribution.
#'
oneDay_larHatching_stochastic_Patch <- function(){

  size = length(private$hatchingFract)
  multipliedProbs = private$NetworkPointer$get_thetaAq("E")
  genotypes = private$NetworkPointer$get_genotypesID()

  counter = which(private$hatchingFract != 0)

  private$larHatching = matrix(data=0L,nrow=2,ncol=size)
  if (length(multipliedProbs)!= size) { multipliedProbs = rep.int( x = multipliedProbs, times = size ) }

  fProbs = cbind(multipliedProbs, 1-multipliedProbs, deparse.level = 0)

  holder = vapply( counter , function(x){rmultinom(n = 1, size = private$hatchingFract[x], prob = fProbs[x,])}, FUN.VALUE = integer(2) )

  private$larHatching[,counter] = holder

  private$larHatching = private$larHatching[1,]
}

#' Stochastic Egg Hatching to Larval Stage
#'
#' Calculate the eggs that have survived and hatched during a day given by: \eqn{\overline{O(T_e)}* \theta_{e}}
#'
oneDay_larHatching_deterministic_Patch <- function(){
  private$larHatching = private$NetworkPointer$get_thetaAq("E") * private$hatchingFract
}

###############################################################################
# eggsFract2
###############################################################################

#' Stochastic Larval Pupation
#'
#' Calculate the number of larvae that will pupate prior to calculating density-dependent effects in \code{\link{oneDay_calcCumulativeLarvalDensityDependentFactor_Patch}}.
#'
oneDay_eggsFract2_stochastic_Patch <- function(){

  #used multiple times/get more complicated in other uses
  genotypesN = private$NetworkPointer$get_genotypesN()
  timeIndex = private$NetworkPointer$get_timeAq("E") + private$NetworkPointer$get_timeAq("L")


  probsHold <- private$AF1dly[[timeIndex]] *
    private$NetworkPointer$get_beta()*
    private$NetworkPointer$get_s()

  for(depthInd in 1:genotypesN){
    probs <- probsHold * private$NetworkPointer$get_drivecubeindex(NULL,NULL,depthInd) *
      private$NetworkPointer$get_tau(NULL,NULL,depthInd)

    #private$eggsCube[ , ,depthInd] <- rpois(n = genotypesN*genotypesN, lambda = probs)
    private$eggsFract2[depthInd] <- sum(rpois(n = genotypesN*genotypesN, lambda = probs))
  }

}

#' Deterministc Larval Pupation
#'
#' Calculate the number of larvae that will pupate prior to calculating density-dependent effects in \code{\link{oneDay_calcCumulativeLarvalDensityDependentFactor_Patch}}.
#'
oneDay_eggsFract2_deterministic_Patch <- function(){

  genotypesN = private$NetworkPointer$get_genotypesN()
  timeIndex = private$NetworkPointer$get_timeAq("E") + private$NetworkPointer$get_timeAq("L")

  #fill offspring cube with parents
  for(slice in 1:genotypesN){

    private$eggsFract2[slice] = sum(private$AF1dly[[timeIndex]] *
                                         private$NetworkPointer$get_beta() *
                                         private$NetworkPointer$get_s() *
                                         private$NetworkPointer$get_drivecubeindex(NULL,NULL,slice) *
                                         private$NetworkPointer$get_tau(NULL,NULL,slice)
    )

  }

}

###############################################################################
# larPupating
###############################################################################

#' Stochastic Larval Pupation
#'
#' Calculate the number of larvae that have transformed into pupae given by \eqn{\overline{O(T_e+T_l)} * \theta_{e} * D(\theta_l,0)}
#' This number follows a binomial distribution.
#'
oneDay_larPupating_stochastic_Patch <- function(){

  size = length(private$eggsFract2)
  multipliedProbs = private$NetworkPointer$get_thetaAq("E") * private$f

  genotypes = private$NetworkPointer$get_genotypesID()

  counter = which(private$eggsFract2 != 0)

  private$larPupating = matrix(data=0L,nrow=2,ncol=size)
  if (length(multipliedProbs)!= size) { multipliedProbs = rep.int( x = multipliedProbs, times = size ) }

  fProbs = cbind(multipliedProbs, 1-multipliedProbs, deparse.level = 0)

  holder = vapply( counter , function(x){rmultinom(n = 1, size = private$eggsFract2[x], prob = fProbs[x,])}, FUN.VALUE = integer(2) )

  private$larPupating[,counter] = holder

  private$larPupating = private$larPupating[1,]

}

#' Deterministc Larval Pupation
#'
#' Calculate the number of larvae that have transformed into pupae given by \eqn{\overline{O(T_e+T_l)} * \theta_{e} * D(\theta_l,0)}
#'
oneDay_larPupating_deterministic_Patch <- function(){
  # change to match code in. not sure if right; ask HMSC
  # https://github.com/Chipdelmal/MGDrivE/commit/c5561a6a94d6f198b20ffb2c31d69726b371c4c3#diff-4c7332795848689b349b9716b1814472L240
  # private$larPupating = private$NetworkPointer$get_thetaAq("E") * private$NetworkPointer$get_thetaAq("L") * private$f * private$eggsFract2
  private$larPupating = private$NetworkPointer$get_thetaAq("E") * private$f * private$eggsFract2

}


###############################################################################
# Adult Male Population Dynamics
###############################################################################

###############################################################################
# numMaleFemale
###############################################################################

#' Stochastic Sex Ratio
#'
#' Calculate the number of males, \eqn{(1-\overline{\phi}) *  \overline{E^{'}}} and females, \eqn{\overline{\phi} *  \overline{E^{'}}}
#' These counts are follow a binomial distribution.
#'
oneDay_numMaleFemale_stochastic_Patch <- function(){

  phi = private$NetworkPointer$get_phi()
  sex_ratio = cbind(1-phi, phi)

  private$numMaleFemale[] = vapply(X=1:length(phi),
                                   FUN=function(x){rmultinom(1,private$ovipositG1[x],sex_ratio[x,])},
                                   FUN.VALUE=integer(2L))
}

#' Deterministc Sex Ratio
#'
#' Calculate the number of males, \eqn{(1-\overline{\phi}) *  \overline{E^{'}}} and females, \eqn{\overline{\phi} *  \overline{E^{'}}}
#'
oneDay_numMaleFemale_deterministic_Patch <- function(){

  phi = private$NetworkPointer$get_phi()

  private$numMaleFemale["M", ] = private$ovipositG1 * (1-phi)
  private$numMaleFemale["F", ] = private$ovipositG1 * phi

}

###############################################################################
# admSurvival
###############################################################################

#' Stochastic Adult Male Survival
#'
#' Daily adult male survival is sampled from a binomial distribution where survival probability is given by \deqn{(1-\mu_{ad}) * \overline{\omega_m}}, where \eqn{\mu_{ad}} is adult mortality rate and \eqn{\overline{\omega_m}} corresponds to
#' genotype-specific mortality effects.
#'
oneDay_admSurvival_stochastic_Patch <- function(){

  size = length(private$ADMdly[[1]])
  multipliedProbs = (1-private$NetworkPointer$get_muAd()) * private$NetworkPointer$getOmega()
  genotypes = private$NetworkPointer$get_genotypesID()

  counter = which(private$ADMdly[[1]] != 0)

  private$admSurvival = matrix(data=0L,nrow=2,ncol=size)
  if (length(multipliedProbs)!= size) { multipliedProbs = rep.int( x = multipliedProbs, times = size ) }

  fProbs = cbind(multipliedProbs, 1-multipliedProbs, deparse.level = 0)

  holder = vapply( counter , function(x){rmultinom(n = 1, size = private$ADMdly[[1]][x], prob = fProbs[x,])}, FUN.VALUE = integer(2) )

  private$admSurvival[,counter] = holder

  private$admSurvival = private$admSurvival[1,]

}

#' Deterministc Adult Male Survival
#'
#' Daily adult male survival is calculated according to \deqn{\overline{\overline{Am_{[t-1]}}} * (1-\mu_{ad}) * \overline{\omega_m}}, where \eqn{\mu_{ad}} is adult mortality rate and \eqn{\overline{\omega_m}} corresponds to
#' genotype-specific mortality effects.
#'
oneDay_admSurvival_deterministic_Patch <- function(){
  private$admSurvival = private$ADMdly[[1]] * (1-private$NetworkPointer$get_muAd()) * private$NetworkPointer$getOmega()
}

###############################################################################
# admPupating
###############################################################################

#' Stochastic Adult Male Pupation
#'
#' Pupation from male pupae to adults is sampled from a binomial distribution for each genotype where the probability of pupation is given by \deqn{\bigg(\overline{\xi_m} * (\theta_{e} * \theta_{p}) * (1-\mu_{ad}) * D(\theta_l,T_p) \bigg)},
#' where \eqn{\overline{\xi_{m}}} is the genotype-specific pupation success probability.
#'
oneDay_admPupating_stochastic_Patch <- function(){

  size = length(private$numMaleFemale["M",])
  multipliedProbs = private$NetworkPointer$get_thetaAq("E") * private$NetworkPointer$get_thetaAq("P") * (1-private$NetworkPointer$get_muAd()) * private$f * private$NetworkPointer$get_xiM()
  genotypes = private$NetworkPointer$get_genotypesID()

  counter = which(private$numMaleFemale["M",] != 0)

  private$admPupating = matrix(data=0L,nrow=2,ncol=size)
  if (length(multipliedProbs)!= size) { multipliedProbs = rep.int( x = multipliedProbs, times = size ) }

  fProbs = cbind(multipliedProbs, 1-multipliedProbs, deparse.level = 0)

  holder = vapply( counter , function(x){rmultinom(n = 1, size = private$numMaleFemale["M",][x], prob = fProbs[x,])}, FUN.VALUE = integer(2) )

  private$admPupating[,counter] = holder

  private$admPupating = private$admPupating[1,]

}

#' Deterministc Adult Male Pupation
#'
#' Adult male emergence is calculated based on the number of male pupae multiplied by \deqn{\bigg(\overline{\xi_m} * (\theta_{e} * \theta_{p}) * (1-\mu_{ad}) * D(\theta_l,T_p) \bigg)},
#' where \eqn{\overline{\xi_{m}}} is the genotype-specific pupation success probability.
#'
oneDay_admPupating_deterministic_Patch <- function(){
  private$admPupating = private$numMaleFemale["M",] * private$NetworkPointer$get_thetaAq("E") * private$NetworkPointer$get_thetaAq("P") * (1-private$NetworkPointer$get_muAd()) * private$f * private$NetworkPointer$get_xiM()
}


###############################################################################
# Adult Female Population Dynamics
###############################################################################

###############################################################################
# af1Survival
###############################################################################

#' Stochastic Adult Female Survival
#'
#' Daily adult female survival is sampled from a binomial distribution where survival probability is given by \deqn{(1-\mu_{ad}) * \overline{\omega_f}}, where \eqn{\mu_{ad}} is adult mortality rate and \eqn{\overline{\omega_f}} corresponds to
#' genotype-specific mortality effects.
#'
oneDay_af1Survival_stochastic_Patch <- function(){

  size = length(private$AF1dly[[1]])
  multipliedProbs = (1-private$NetworkPointer$get_muAd()) * private$NetworkPointer$getOmega()
  genotypes = private$NetworkPointer$get_genotypesID()

  counter = which(private$AF1dly[[1]] != 0)

  private$af1Survival = matrix(data=0L,nrow=2,ncol=size)
  if (length(multipliedProbs)!= size) { multipliedProbs = rep.int( x = multipliedProbs, times = size ) }

  fProbs = cbind(multipliedProbs, 1-multipliedProbs, deparse.level = 0)

  holder = vapply( counter , function(x){rmultinom(n = 1, size = private$AF1dly[[1]][x], prob = fProbs[x,])}, FUN.VALUE = integer(2) )

  private$af1Survival[,counter] = holder

  private$af1Survival = private$af1Survival[1,]

}

#' Deterministc Adult Female Survival
#'
#' Daily adult female survival is calculated according to \deqn{\overline{\overline{Af_{[t-1]}}} * (1-\mu_{ad}) * \overline{\omega_f}}, where \eqn{\mu_{ad}} is adult mortality rate and \eqn{\overline{\omega_f}} corresponds to
#' genotype-specific mortality effects.
#'
oneDay_af1Survival_deterministic_Patch <- function(){
  private$af1Survival = private$AF1dly[[1]] * (1-private$NetworkPointer$get_muAd()) * private$NetworkPointer$getOmega()
}

########################################################################
# af1Pupation
########################################################################

#' Stochastic Adult Female Pupation
#'
#' Pupation from female pupae to adults is sampled from a binomial distribution for each genotype where the probability of pupation is given by \deqn{\bigg(\overline{\xi_f} * (\theta_{e} * \theta_{p}) * (1-\mu_{ad}) * D(\theta_l,T_p) \bigg)},
#' where \eqn{\overline{\xi_{f}}} is the genotype-specific pupation success probability.
#'
oneDay_af1Pupation_stochastic_Patch <- function(){

  size = length(private$numMaleFemale["F",])
  multipliedProbs = private$NetworkPointer$get_thetaAq("E") * private$NetworkPointer$get_thetaAq("P") * (1-private$NetworkPointer$get_muAd()) * private$f * private$NetworkPointer$get_xiF()
  genotypes = private$NetworkPointer$get_genotypesID()

  counter = which(private$numMaleFemale["F",] != 0)

  private$af1Pupation = matrix(data=0L,nrow=2,ncol=size)
  if (length(multipliedProbs)!= size) { multipliedProbs = rep.int( x = multipliedProbs, times = size ) }

  fProbs = cbind(multipliedProbs, 1-multipliedProbs, deparse.level = 0)

  holder = vapply( counter , function(x){rmultinom(n = 1, size = private$numMaleFemale["F",][x], prob = fProbs[x,])}, FUN.VALUE = integer(2) )

  private$af1Pupation[,counter] = holder

  private$af1Pupation = private$af1Pupation[1,]

}

#' Deterministc Adult Female Pupation
#'
#' Adult female emergence is calculated based on the number of female pupae multiplied by \deqn{\bigg(\overline{\xi_f} * (\theta_{e} * \theta_{p}) * (1-\mu_{ad}) * D(\theta_l,T_p) \bigg)},
#' where \eqn{\overline{\xi_{f}}} is the genotype-specific pupation success probability.
#'
oneDay_af1Pupation_deterministic_Patch <- function(){
  private$af1Pupation = private$numMaleFemale["F",] * private$NetworkPointer$get_xiF() * private$NetworkPointer$get_thetaAq("E") * private$NetworkPointer$get_thetaAq("P") * (1-private$NetworkPointer$get_muAd()) * private$f

}

########################################################################
# af1Mating
########################################################################

#' Stochastic Mating
#'
#' Mating for each newly emerging adult female genotype is sampled from a multinomial distribution with probabilities equal to the adult male population vector multiplied by \eqn{\overline{\eta}}, genotype-specific male mating fitness.
#'
oneDay_af1Mating_stochastic_Patch <- function(){

  eta = private$NetworkPointer$get_eta()
  genotypesN = private$NetworkPointer$get_genotypesN()
  genotypesID = private$NetworkPointer$get_genotypesID()

  counter = which(private$af1Pupation!=0)
  private$af1Mating = matrix(data=0L,nrow=genotypesN,ncol=genotypesN,dimnames=list(genotypesID,genotypesID))

  # if no males to mate with stop calculation
  if(sum(private$ADMnew)==0){
    # update population
    private$AF1new = private$af1Survival + private$af1Mating
  # else proceed to calculate mating
  } else {
    mate_prob = private$ADMnew * eta
    holder = vapply(X=counter,FUN=function(x){rmultinom(n=1,size=private$af1Pupation[x],prob=mate_prob)},FUN.VALUE=integer(genotypesN))
    private$af1Mating[counter,] = t(holder)
    # update population
    private$AF1new = private$af1Survival + private$af1Mating
  }
}

#' Deterministc Mating
#'
#' Mating is calculated as the outer product of newly emerging adult females and adult males, modulated by \eqn{\overline{\eta}}, genotype-specific male mating fitness.
#'
oneDay_af1Mating_deterministic_Patch <- function(){

  private$af1Mating = private$af1Pupation %o% normalise(private$ADMnew * private$NetworkPointer$get_eta())
  private$AF1new = private$af1Survival + private$af1Mating
}


###############################################################################
# Helper Functions
###############################################################################

#' Shift Population
#'
#' Update larval and adult populations daily at end of time step, calls \code{\link{shiftAndUpdatePopVector}}
#'
oneDay_updatePopulation_Patch <- function(){

  # delayed population vectors
  shiftAndUpdatePopVector(popVector = private$LARdly,newPop = private$LARnew)
  shiftAndUpdatePopVector(popVector = private$ADMdly,newPop = private$ADMnew)
  shiftAndUpdatePopVector(popVector = private$AF1dly,newPop = private$AF1new)

  # population tracking arrays
  private$LAR[[private$NetworkPointer$get_tNow()]] = private$LARnew
  private$ADM[[private$NetworkPointer$get_tNow()]] = private$ADMnew
  private$AF1[[private$NetworkPointer$get_tNow()]] = private$AF1new

}

Patch$set(which = "public",name = "oneDay_updatePopulation",
          value = oneDay_updatePopulation_Patch, overwrite = TRUE
)

#' Calculate Larval density-dependent Mortality
#'
#' Calculate \deqn{D(\theta_l,T_x) = \left\{ \begin{array}{ll} \theta_{l[0]}^{'}=\theta_l	& \quad i = 0 \\ \theta_{l[i+1]}^{'} = \theta_{l[i]}^{'} *F(\overline{L_{[t-i-T_x]}})	& \quad i \leq T_l \end{array}    \right.},
#' the effect of density-dependence on larval mortality for one day.
#'
oneDay_calcLarvalDensityDependentFactor_Patch <- function(){
  # browser()
  summedPop = sum(private$LARdly[[2]])
  alpha = private$NetworkPointer$get_alpha(private$patchID)
  L = private$NetworkPointer$get_timeAq("L")
  private$larDDMortal = (alpha/(alpha+summedPop))^(1/L)
}

Patch$set(which = "public",name = "oneDay_calcLarvalDensityDependentFactor",
  value = oneDay_calcLarvalDensityDependentFactor_Patch,overwrite = TRUE
)

#' Calculate Cumulative Larval density-dependent Mortality
#'
#' Calculate \deqn{F(L[t])=\Bigg(\frac{\alpha}{\alpha+\sum{\overline{L[t]}}}\Bigg)^{1/T_l}\\}, the cumulative effect of density-dependence on larval mortality over the entire duration of larval stage.
#'
oneDay_calcCumulativeLarvalDensityDependentFactor_Patch <- function(){
  # browser()
  private$f = private$NetworkPointer$get_thetaAq("L")

  alpha = private$NetworkPointer$get_alpha(private$patchID)
  larvalDuration = private$NetworkPointer$get_timeAq("L")
  for(i in 1:larvalDuration){
    summedPop = sum(private$LARdly[[1+i]])
    private$f = private$f * (alpha/(alpha+summedPop))^(1/larvalDuration)
  }
  # print("hi")
}

Patch$set(which = "public",name = "oneDay_calcCumulativeLarvalDensityDependentFactor",
  value = oneDay_calcCumulativeLarvalDensityDependentFactor_Patch,overwrite = TRUE
)

#' Calculate Cumulative Pupae density-dependent Mortality
#'
#' Calculate \deqn{D(\theta_l,T_{P}) = \left\{ \begin{array}{ll} \theta_{l[0]}^{'}=\theta_l 								& \quad i = 0 \\ \theta_{l[i+1]}^{'} = \theta_{l[i]}^{'} *F(\overline{L_{[t-i-T_P]}})	& \quad i \leq T_{P} \end{array} \right.}, the cumulative effect of density-dependence on larval mortality over the entire duration of larval and pupal stages.
#'
oneDay_calcCumulativePupaDensityDependentFactor_Patch <- function(){
  # browser()
  private$f = private$NetworkPointer$get_thetaAq("L")

  alpha = private$NetworkPointer$get_alpha(private$patchID)
  larvalDuration = private$NetworkPointer$get_timeAq("L")
  pupalDuration = private$NetworkPointer$get_timeAq("P")

  for(i in 1:larvalDuration){
    summedPop = sum(private$LARdly[[1+i+pupalDuration]])
    private$f = private$f * (alpha/(alpha+summedPop))^(1/larvalDuration)
  }
}

Patch$set(which = "public",name = "oneDay_calcCumulativePupaDensityDependentFactor",
  value = oneDay_calcCumulativePupaDensityDependentFactor_Patch,overwrite = TRUE
)
