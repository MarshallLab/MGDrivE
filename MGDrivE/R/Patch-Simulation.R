###############################################################################
#        ____        __       __
#       / __ \____ _/ /______/ /_
#      / /_/ / __ `/ __/ ___/ __ \
#     / ____/ /_/ / /_/ /__/ / / /
#    /_/    \__,_/\__/\___/_/ /_/
#
#   Patch Class Mosquito Population Simulation
#   Marshall Lab
#   jared_bennett@berkeley.edu
#   December 2019
#
###############################################################################
###############################################################################
# Daily Simulation
###############################################################################

#' Daily Population Dynamics for a Patch
#'
#' Run population dynamics (including migration) for this patch. \cr
#' Performed in this order, see the following for each function: \cr
#' Adult Death: \code{\link{oneDay_adultDeath_deterministic_Patch}} or \code{\link{oneDay_adultDeath_stochastic_Patch}} \cr
#' Pupa Death/Maturation: \code{\link{oneDay_pupaDM_deterministic_Patch}} or \code{\link{oneDay_pupaDM_stochastic_Patch}} \cr
#' Larva Death/Maturation: \code{\link{oneDay_larvaDM_deterministic_Patch}} or \code{\link{oneDay_larvaDM_stochastic_Patch}} \cr
#' Egg Death/Maturation: \code{\link{oneDay_eggDM_deterministic_Patch}} or \code{\link{oneDay_eggDM_stochastic_Patch}} \cr
#' Pupation: \code{\link{oneDay_pupation_deterministic_Patch}} or \code{\link{oneDay_pupation_stochastic_Patch}} \cr
#' Releases: \code{\link{oneDay_releases_Patch}} \cr
#' Mating: \code{\link{oneDay_mating_deterministic_Patch}} or \code{\link{oneDay_mating_stochastic_Patch}} \cr
#' Lay Eggs: \code{\link{oneDay_oviposit_deterministic_Patch}} or \code{\link{oneDay_oviposit_stochastic_Patch}} \cr
#' Release Eggs: \code{\link{oneDay_eggReleases_Patch}} \cr
#'
#' @keywords internal
#'
oneDay_PopDynamics_Patch <- function(){

  ##################
  #Death/Maturation#
  ##################
  self$oneDay_adultD()
  self$oneDay_pupaDM()
  self$oneDay_larvaDM()
  self$oneDay_eggDM()

  ##########
  #Pupation#
  ##########
  self$oneDay_pupation()

  ##########
  #Releases#
  ##########
  self$oneDay_releases()

  ########
  #Mating#
  ########
  self$oneDay_mating()

  ##########
  #Lay Eggs#
  ##########
  self$oneDay_layEggs()
  self$oneDay_releaseEggs()

}


###############################################################################
# Death/Maturation
###############################################################################
##########
# Adult
##########
#' Deterministic Adult Survival
#'
#' Daily adult survival is calculated according to \deqn{\overline{\overline{Af_{[t-1]}}} * (1-\mu_{ad}) * \overline{\omega_{m/f}}},
#' where \eqn{\mu_{ad}} corresponds to adult mortality rate and \eqn{\overline{\omega_{m/f}}}
#' corresponds to genotype-specific male/female mortality effects.
#'
#' @keywords internal
#'
oneDay_adultDeath_deterministic_Patch <- function(){

  # probability of survival
  probHolder = (1-private$NetworkPointer$get_muAd()) * private$NetworkPointer$get_omega()

  # males that live through the day
  private$popMale[] = private$popMale * probHolder

  # females that live through the day
  # this works becase female genotypes are rows, and R applies column-wise, thereby
  #  properly applying the vector to each female genotype in order
  private$popFemale[] = private$popFemale * probHolder

}

#' Stochastic Adult Survival
#'
#' Daily adult survival is sampled from a binomial distribution where survival
#' probability is given by \deqn{(1-\mu_{ad}) * \overline{\omega_m/f}}.
#' \eqn{\mu_{ad}} corresponds to adult mortality rate and \eqn{\overline{\omega_m/f}}
#' corresponds to genotype-specific mortality effects.
#'
#' @keywords internal
#'
oneDay_adultDeath_stochastic_Patch <- function(){

  # probability of survival
  probHolder = (1-private$NetworkPointer$get_muAd()) * private$NetworkPointer$get_omega()

  # males that live through the day
  private$popMale[] <- rbinom(n = private$NetworkPointer$get_genotypesN(),
                              size = private$popMale,
                              prob = probHolder)

  # females that live through the day
  # this also works because of how R applies things and fills matrices
  private$popFemale[] = rbinom(n = (private$NetworkPointer$get_genotypesN())^2,
                               size = private$popFemale,
                               prob = probHolder)

}


##########
# Pupa
##########
#' Deterministic Pupa Death and Pupation
#'
#' Daily pupa survival is calculated according to \deqn{\overline{P_{[t-1]}} * (1-\mu_{aq})},
#' where \eqn{\mu_{aq}} corresponds to daily non-density-dependent aquatic mortality. \cr
#' See \code{\link{parameterizeMGDrivE}} for how these parameters are derived.
#'
#' @keywords internal
#'
oneDay_pupaDM_deterministic_Patch <- function(){

  # things to reuse
  life <- 1-private$NetworkPointer$get_muAq()
  pupaStart <- private$NetworkPointer$get_timeAq(stage = 'E') + private$NetworkPointer$get_timeAq(stage = 'L') + 1
  pupaEnd <- private$NetworkPointer$get_timeAq()

  # Treat last day differently because the pupa start to hatch
  #  This does not handle the continuous to discrete time conversion artifact, see
  #  the pupation function for that.
  private$popHolder[] <- private$popAquatic[ ,pupaEnd] * life

  # check if there are other days, then
  # run loop backwards to move populations as we go
  if((pupaEnd - pupaStart) > 0){
    for(i in (pupaEnd-1):pupaStart){
      private$popAquatic[ ,i+1] = private$popAquatic[ ,i] * life
    } # end loop
  }

}

#' Stochastic Pupa Death and Pupation
#'
#' Daily pupa survival is sampled from a binomial distribution, where survival
#' probability is given by \deqn{1-\mu_{aq}}. \eqn{\mu_{aq}} corresponds
#' to daily non-density-dependent aquatic mortality. \cr
#' See \code{\link{parameterizeMGDrivE}} for how these parameters are derived.
#'
#' @keywords internal
#'
oneDay_pupaDM_stochastic_Patch <- function(){

  # things to reuse
  nGeno <- private$NetworkPointer$get_genotypesN()
  life <- 1-private$NetworkPointer$get_muAq()
  pupaStart <- private$NetworkPointer$get_timeAq(stage = 'E') + private$NetworkPointer$get_timeAq(stage = 'L') + 1
  pupaEnd <- private$NetworkPointer$get_timeAq()

  # Treat last day differently because the pupa start to hatch
  #  This does not handle the continuous to discrete time conversion artifact, see
  #  the pupation function for that.
  private$popHolder[] <- rbinom(n = nGeno,
                                size = private$popAquatic[ ,pupaEnd],
                                prob = life)

  # check if there are other days, then
  # run loop backwards to move populations as we go
  if((pupaEnd - pupaStart) > 0){
    for(i in (pupaEnd-1):pupaStart){
      private$popAquatic[ ,i+1] <- rbinom(n = nGeno,
                                          size = private$popAquatic[ ,i],
                                          prob = life)
    } # end loop
  }

}


##########
# Larva
##########
#' Deterministic Larva Death and Pupation
#'
#' Calculate the number of larvae surviving from day to day, given by:
#' \deqn{\overline{L_{[t-1]}} * (1-\mu_{aq}) * F(L)}. F(L), the density dependence
#' is calculated as \deqn{F(L[t])=\Bigg(\frac{\alpha}{\alpha+\sum{\overline{L[t]}}}\Bigg)^{1/T_l}}.
#' See \code{\link{parameterizeMGDrivE}} for how these parameters are derived.
#' Pupation has no parameters, so the final day of larvae naturally enter the pupal state.
#'
#' @keywords internal
#'
oneDay_larvaDM_deterministic_Patch <- function(){

  # things to reuse
  larvaStart <- private$NetworkPointer$get_timeAq(stage = 'E') + 1
  larvaEnd <- private$NetworkPointer$get_timeAq(stage = 'E') + private$NetworkPointer$get_timeAq(stage = 'L')

  # calculate density dependence
  # (alpha/(alpha + L))^(1/t_l)
  density <- (private$NetworkPointer$get_alpha(ix = private$patchID)/
                (private$NetworkPointer$get_alpha(ix = private$patchID) + sum(x = private$popAquatic[,larvaStart:larvaEnd])) )^
             (1/private$NetworkPointer$get_timeAq(stage = 'L'))

  # total life rate
  life <- density*(1-private$NetworkPointer$get_muAq())

  # run loop backwards to move populations as we go
  for(i in larvaEnd:larvaStart){
    private$popAquatic[ ,i+1] = private$popAquatic[ ,i] * life
  } # end loop

}

#' Stochastic Larva Death and Pupation
#'
#' The daily number of larvae surviving is drawn from a binomial distribution, where
#' survival probability is given by \deqn{(1-\mu_{aq}) * F(L)}. F(L), the density dependence
#' is calculated as \deqn{F(L[t])=\Bigg(\frac{\alpha}{\alpha+\sum{\overline{L[t]}}}\Bigg)^{1/T_l}}.
#' See \code{\link{parameterizeMGDrivE}} for how these parameters are derived.
#' Pupation has no parameters, so the final day of larvae naturally enter the pupal state.
#'
#' @keywords internal
#'
oneDay_larvaDM_stochastic_Patch <- function(){

  # things to reuse
  nGeno <- private$NetworkPointer$get_genotypesN()
  larvaStart <- private$NetworkPointer$get_timeAq(stage = 'E') + 1
  larvaEnd <- private$NetworkPointer$get_timeAq(stage = 'E') + private$NetworkPointer$get_timeAq(stage = 'L')

  # calculate density dependence
  # (alpha/(alpha + L))^(1/t_l)
  density <- (private$NetworkPointer$get_alpha(ix = private$patchID)/
                (private$NetworkPointer$get_alpha(ix = private$patchID) + sum(x = private$popAquatic[,larvaStart:larvaEnd])) )^
    (1/private$NetworkPointer$get_timeAq(stage = 'L'))

  # total life rate
  life <- density*(1-private$NetworkPointer$get_muAq())

  # run loop backwards to move populations as we go
  for(i in larvaEnd:larvaStart){
    private$popAquatic[ ,i+1] = rbinom(n = nGeno,
                                       size = private$popAquatic[ ,i],
                                       prob = life)
  } # end loop

}


##########
# Egg
##########
#' Deterministic Egg Death and Pupation
#'
#' Daily egg survival is calculated according to \deqn{\overline{E_{[t-1]}} * (1-\mu_{aq})},
#' where \eqn{\mu_{aq}} corresponds to daily non-density-dependent aquatic mortality.
#' Eggs transition into larvae at the end of \eqn{T_e}. \cr
#' See \code{\link{parameterizeMGDrivE}} for how these parameters are derived.
#'
#' @keywords internal
#'
oneDay_eggDM_deterministic_Patch <- function(){

  # things to reuse
  life <- 1-private$NetworkPointer$get_muAq()
  eggEnd <- private$NetworkPointer$get_timeAq(stage = 'E')

  # run loop backwards to move populations as we go
  for(i in eggEnd:1){
    private$popAquatic[ ,i+1] = private$popAquatic[ ,i] * life
  } # end loop

}

#' Stochastic Egg Death and Pupation
#'
#' Daily egg survival is sampled from a binomial distribution, where survival
#' probability is given by \eqn{1-\mu_{aq}}. \eqn{\mu_{aq}} corresponds
#' to daily non-density-dependent aquatic mortality. \cr
#' Eggs transition into larvae at the end of \eqn{T_e}. \cr
#' See \code{\link{parameterizeMGDrivE}} for how these parameters are derived.
#'
#' @keywords internal
#'
oneDay_eggDM_stochastic_Patch <- function(){

  # things to reuse
  nGeno <- private$NetworkPointer$get_genotypesN()
  life <- 1-private$NetworkPointer$get_muAq()
  eggEnd <- private$NetworkPointer$get_timeAq(stage = 'E')

  # run loop backwards to move populations as we go
  for(i in eggEnd:1){
    private$popAquatic[ ,i+1] = rbinom(n = nGeno,
                                       size = private$popAquatic[ ,i],
                                       prob = life)
  } # end loop

}


###############################################################################
# Pupation
###############################################################################

#' Deterministic Pupation
#'
#' Pupa first undergo one extra day of survival, calculated as \deqn{\overline{P_{[t-1]}} * (1-\mu_{ad})}.
#' This is an artifact of the conversion from continuous to discrete time (as mentioned
#' in the original Hancock paper this model is derived from). \cr
#' Then, pupation into adult males is calculated as \deqn{(1-\overline{\phi}) * \overline{P_{[t]}}}
#' and into adult females as \deqn{\overline{\phi} * \overline{P_{[t]}}}
#'
#' @keywords internal
#'
oneDay_pupation_deterministic_Patch <- function(){

  # one extra death to match continuous time math
  #  This is an artifact of being discrete time
  private$popHolder[] = private$popHolder * (1-private$NetworkPointer$get_muAd())

  # get sex ratio for emergence
  phi = private$NetworkPointer$get_phi()

  # perform genotype-specific sex ratio and sex-dependent emergence
  private$popMale[] = private$popMale + private$popHolder * (1-phi) * private$NetworkPointer$get_xiM()
  private$popUnmated[] = private$popUnmated + private$popHolder * phi * private$NetworkPointer$get_xiF()

}

#' Stochastic Pupation
#'
#' Pupa first undergo one extra day of survival, calculated as a binomial over
#' \deqn{\overline{P_{[t-1]}} * (1-\mu_{ad})}.
#' This is an artifact of the conversion from continuous to discrete time (as mentioned
#' in the original Hancock paper this model is derived from). \cr
#' Then, pupation is sampled from a binomial, where \eqn{(1-\overline{\phi})} is
#' the genotype-specific probability of becoming male, and \eqn{\overline{\phi}}
#' is the genotype-specific of becoming female.
#'
#' @keywords internal
#'
oneDay_pupation_stochastic_Patch <- function(){

  # reuse
  nGeno = private$NetworkPointer$get_genotypesN()

  # there needs to be 1 extra death for discrete time
  private$popHolder[] = rbinom(n = nGeno,
                               size = private$popHolder,
                               prob =  1-private$NetworkPointer$get_muAd())

  # pull male/female distinction
  # only need to pull female, male is just popHolder - female
  private$popPupSex[] = rbinom(n = nGeno,
                             size = private$popHolder,
                             prob = private$NetworkPointer$get_phi() )

  # genotype-specific and sex-dependent emergence
  private$popMale[] = private$popMale + rbinom(n = nGeno,
                                               size = private$popHolder - private$popPupSex,
                                               prob = private$NetworkPointer$get_xiM() )
  private$popUnmated[] = private$popUnmated + rbinom(n = nGeno,
                                                     size = private$popPupSex,
                                                     prob = private$NetworkPointer$get_xiF() )

}


###############################################################################
# Releases
###############################################################################

#' Release Male/Female/Mated-Female Mosquitoes in a Patch
#'
#' Based on this patch's release schedule, \code{\link{generateReleaseVector}},
#' this function handles daily releases.
#'
#' @keywords internal
#'
oneDay_releases_Patch <- function(){

  # things to reuse
  tNow <- private$NetworkPointer$get_tNow()

  ##########
  # male
  ##########
  if((length(private$maleReleases) > 0) && (private$maleReleases[[1]]$tRelease <= tNow) ){
    # combine release
    # 1st column is index, 2nd column is amount
    idxRel <- private$maleReleases[[1]]$nRelease[ ,1]
    private$popMale[idxRel] = private$popMale[idxRel] + private$maleReleases[[1]]$nRelease[ ,2]
    # remove finished release
    private$maleReleases[[1]] = NULL
  }

  ##########
  # female
  ##########
  if((length(private$femaleReleases) > 0) && (private$femaleReleases[[1]]$tRelease <= tNow) ){
    # combine unmated females
    # 1st column is index, 2nd column is amount
    idxRel <- private$femaleReleases[[1]]$nRelease[ ,1]
    private$popUnmated[idxRel] = private$popUnmated[idxRel] + private$femaleReleases[[1]]$nRelease[ ,2]
    # remove finished release
    private$femaleReleases[[1]] = NULL
  }

  ##########
  # mated female
  ##########
  if((length(private$matedFemaleReleases) > 0) && (private$matedFemaleReleases[[1]]$tRelease <= tNow) ){
    # combine mated females
    # 1st column is female index, 2nd column is male index, 3rd column is amount
    idxRel <- private$matedFemaleReleases[[1]]$nRelease[ ,1:2,drop=FALSE]
    private$popFemale[idxRel] = private$popFemale[idxRel] + private$matedFemaleReleases[[1]]$nRelease[ ,3]
    # remove finished release
    private$matedFemaleReleases[[1]] = NULL
  }

}

#' Release Eggs in a Patch
#'
#' Based on this patch's release schedule, \code{\link{generateReleaseVector}},
#' this function handles daily egg releases.
#'
#' @keywords internal
#'
oneDay_eggReleases_Patch <- function(){

  # things to reuse
  tNow <- private$NetworkPointer$get_tNow()

  ##########
  # eggs
  ##########
  if((length(private$eggReleases) > 0) && (private$eggReleases[[1]]$tRelease <= tNow) ){
    # combine eggs
    # 1st column is index, 2nd column is amount
    idxRel <- private$eggReleases[[1]]$nRelease[ ,1]
    private$popAquatic[idxRel,1] = private$popAquatic[idxRel,1] + private$eggReleases[[1]]$nRelease[ ,2]
    # remove finished release
    private$eggReleases[[1]] = NULL
  }

}


###############################################################################
# Mating
###############################################################################

#' Deterministic Mating
#'
#' Mating is calculated as the outer product of newly emerging adult females and
#' all-current adult males, modulated by \eqn{\overline{\overline{\eta}}}, the genotype-specific
#' male mating fitness. \eqn{\overline{\overline{\eta}}} corresponds to all female (rows)
#' and male (columns) genotypes, to perform any type of assortative mating. \cr
#' If there are no adult males, the unmated females experience one day of death,
#' calculated as \deqn{\overline{Af_t} * (1-\mu_{ad}) * \overline{\omega_f}}, and
#' remain unmated until tomorrow.
#'
#' @keywords internal
#'
oneDay_mating_deterministic_Patch <- function(){

  # things to reuse
  nGeno = private$NetworkPointer$get_genotypesN()
  mProb = numeric(length = nGeno)
  mort = (1-private$NetworkPointer$get_muAd()) * private$NetworkPointer$get_omega()

  # loop over each female genotype, mate with available males
  for(i in 1:nGeno){

    # check if there are females to mate, skip if not
    if(private$popUnmated[i] > 0){
      # get mating prob for males for each female genotype
      mProb[] <- normalise(private$popMale * private$NetworkPointer$get_eta(i))

      # check if males have a chance of mating these females, if not then skip
      if(sum(mProb) > 0){
        # mate
        private$popFemale[i, ] = private$popFemale[i, ] + private$popUnmated[i] * mProb

        # clear unmated vector
        private$popUnmated[i] = 0
      } else {
        # male probs are 0, females experience death
        # females don't mate, make them die, don't clear vector
        private$popUnmated[i] = private$popUnmated[i] * mort[i]

      } # end male check
    } # end female check
  } # end mating loop

}

#' Stochastic Mating
#'
#' Mating for each newly emerging adult female genotype is sampled from a multinomial
#' distribution with probabilities equal to the adult male population vector
#' multiplied by \eqn{\overline{\overline{\eta}}}, the genotype-specific
#' male mating fitness. \eqn{\overline{\overline{\eta}}} corresponds to all female (rows)
#' and male (columns) genotypes, to perform any type of assortative mating. \cr
#' If there are no adult males, the unmated females experience one day of death,
#' sampled from a binomial distribution parameterized by \deqn{(1-\mu_{ad}) * \overline{\omega_f}},
#' and remain unmated until tomorrow.
#'
#' @keywords internal
#'
oneDay_mating_stochastic_Patch <- function(){

  # things to reuse
  nGeno = private$NetworkPointer$get_genotypesN()
  mProb = numeric(length = nGeno)
  mort = (1-private$NetworkPointer$get_muAd()) * private$NetworkPointer$get_omega()

  # loop over each female genotype, mate with available males
  for(i in 1:nGeno){

    # check if there are females to mate, skip if not
    if(private$popUnmated[i] > 0){
      # get mating prob for males for each female genotype
      mProb[] <- normalise(private$popMale * private$NetworkPointer$get_eta(i))

      # check if males have a chance of mating these females, if not then skip
      if(sum(mProb) > 0){
        # mate
        private$popFemale[i, ] = private$popFemale[i, ] + rmultinom(n = 1,
                                                                    size = private$popUnmated[i],
                                                                    prob = mProb)

        # clear unmated vector
        private$popUnmated[i] = 0
      } else {
        # male probs are 0, females experience death
        # females don't mate, make them die, don't clear vector
        private$popUnmated[i] = rbinom(n = 1, size = private$popUnmated[i], prob = mort[i])

      } # end male check
    } # end female check
  } # end mating loop

}


###############################################################################
# Lay Eggs
###############################################################################

#' Deterministic Oviposition
#'
#' Calculate the number of eggs oviposited by female mosquitoes following:
#' \deqn{\overline{O(T_x)} = \sum_{j=1}^{n} \Bigg( \bigg( (\beta*\overline{s} * \overline{ \overline{Af_{[t]}}}) * \overline{\overline{\overline{Ih}}} \bigg) * \Lambda \Bigg)^{\top}_{ij}}
#'
#' @keywords internal
#'
oneDay_oviposit_deterministic_Patch <- function(){

  # things to reuse
  nGeno = private$NetworkPointer$get_genotypesN()
  femBetaS = private$popFemale * (private$NetworkPointer$get_beta() * private$NetworkPointer$get_s())

  #fill offspring cube with parents
  for(slice in 1:nGeno){

    private$popAquatic[slice, 1] = sum(femBetaS *
                                       private$NetworkPointer$get_drivecubeindex(NULL,NULL,slice) *
                                       private$NetworkPointer$get_tau(NULL,NULL,slice) )

  } # end loop over genotypes

}

#' Stochastic Oviposition
#'
#' Calculate the number of eggs oviposited by female mosquitoes following:
#' \deqn{\overline{O(T_x)} = \sum_{j=1}^{n} \Bigg( \bigg( (\beta*\overline{s} * \overline{ \overline{Af_{[t]}}}) * \overline{\overline{\overline{Ih}}} \bigg) * \Lambda  \Bigg)^{\top}_{ij}}
#' The deterministic result for number of eggs is used as the mean of a Poisson-distributed
#' number of actual eggs oviposited.
#'
#' @keywords internal
#'
oneDay_oviposit_stochastic_Patch <- function(){

  # things to reuse
  nGeno = private$NetworkPointer$get_genotypesN()
  femBetaS = private$popFemale * (private$NetworkPointer$get_beta() * private$NetworkPointer$get_s())

  #fill offspring cube with parents
  for(slice in 1:nGeno){

    private$popAquatic[slice, 1] = sum(rpois(n = nGeno*nGeno,
                                             lambda = femBetaS *
                                                      private$NetworkPointer$get_drivecubeindex(NULL,NULL,slice) *
                                                      private$NetworkPointer$get_tau(NULL,NULL,slice) )
                                       )
  } # end loop over genotypes

}
