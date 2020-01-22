###############################################################################
#        ____        __       __
#       / __ \____ _/ /______/ /_
#      / /_/ / __ `/ __/ ___/ __ \
#     / ____/ /_/ / /_/ /__/ / / /
#    /_/    \__,_/\__/\___/_/ /_/#
#
#   Patch Class Implementation
#   Marshall Lab
#   jared_bennett@berkeley.edu
#   December 2019
#
###############################################################################

#' Calculate Distribution of Larval Population
#'
#' This hidden function calculates the distribution of larvae through time by
#' treating the larval-stage as a discrete-time markov chain, and solving for the
#' stationary distribution. As the only aquatic population known for initializing
#' MGDrivE is the equilibrium larval population, this acts as an anchor from which
#' to calculate the egg and pupae distributions from (see \code{\link{set_initialPopulation_Patch}}).
#'
#' @param mu Double, death rate
#' @param t Integer, stage time
#'
calcLarvalDist <- function(mu, t){
  # treat the world as a markov chain, solve for the stable distribution
  # set diagonals to 1
  M = diag(x = t)

  # set upper off-diagonal to negative of 1 minus deathrate
  M[cbind(1:(t-1),2:t)] = mu - 1

  # find inverse and take only first row
  #  not most efficient way to invert, but this is small and direct
  num = solve(M)[1, ]

  # normalize
  ans = num / sum(num)

  # return
  return(ans)
}

# # original function
# # mu: daily probability of death
# # t: number of days in larval state
# qsd <- function(mu,t){
#   M <- matrix(data = 0,nrow = t+1,ncol = t+1,dimnames = list(c(paste0(1:t),"D"),c(paste0(1:t),"D")))
#   i <- 1:(t-1); j <- 2:t
#   for(k in seq_along(i)){
#     M[i[k],j[k]] <- 1-mu
#   }
#   M[1:t,t+1] <- mu
#   M[t+1,t+1] <- 1
#   M[t,t+1] <- 1
#   MM <- M[1:t,1:t]
#   e <- rep(1,t)
#   pi <- c(1,rep(0,t-1))
#   num <- (pi %*% solve(diag(nrow(MM)) - MM))
#   ans <- num / as.numeric((num %*% e))
#   return(ans)
# }

# # this one sets entire aquatic population in one step
# popDR <- function(muAq, alpha, larvalEQ, timeAq){
#
#   epLife <- 1 - muAq
#   lLife <- (alpha/(alpha + larvalEQ))^(1/timeAq[['L']]) * epLife
#
#   M = diag(x = timeAq[['L']])
#
#   # set upper off-diagonal to negative of 1 minus deathrate
#   M[cbind(seq_len(timeAq[['L']]-1),seq_len(timeAq[['L']]-1)+1)] = (1-lLife) - 1
#
#   # find inverse and take only first row
#   num = solve(M)[1, ]
#
#   # normalize
#   ans = num / sum(num)
#
#   # extend to eggs and first pupae
#   ans <- c(ans[1]/epLife^(timeAq[['E']]:1),ans, ans[timeAq[['L']]]*lLife)
#   # extend to rest of pupae
#   ans <- c(ans, tail(x = ans, n = 1)*epLife^(seq_len(timeAq[['P']]-1)) )
#
#   return(ans)
# }

#' Set Initial Population
#'
#' This hidden function distributes the population at time 0 in the steady-state
#' conformation. This involves finding the number of mosquitoes in each day of the
#' aquatic stages, and then splitting adults into male and female. Each stage is
#' appropriately split amongst the initial population genotypes (see \code{\link{parameterizeMGDrivE}}).
#' It internally calls \code{\link{calcLarvalDist}} to determine the distribution
#' of larvae before setting the eggs and pupa from that.
#'
#' @param adultEQ Equilibrium number of adults
#' @param larvalEQ Equilibrium number of larvae
#' @param adultRatioF Genotype specific ratio for adult females
#' @param adultRatioM Genotype specific ratio for adult males
#' @param larvalRatio Genotype specific ratio for larvae
#' @param timeAq Time for each aquatic stage
#' @param muAq Aquatic death rate
#' @param alpha Density-dependent centering parameter
#'
set_initialPopulation_Patch <- function(adultEQ = adultEQ, larvalEQ = larvalEQ,
                                        adultRatioF = adultRatioF, adultRatioM = adultRatioM,
                                        larvalRatio = larvalRatio,
                                        timeAq = timeAq, muAq = muAq, alpha = alpha){

  # reuse these
  # epLife is the chance of not dying for eggs and pupa
  #  There is no density dependence in this one
  # lLife includes the density dependence for larvae
  epLife <- 1 - muAq
  lLife <- (alpha/(alpha + larvalEQ))^(1/timeAq[['L']]) * epLife


  ##########
  # set aquatic population breakdown
  ##########
  # initial larval pop based on equilibrium calcs
  larvalDist <- calcLarvalDist(mu = 1-lLife, t = timeAq[['L']])

  # set larvae
  for(i in 1:timeAq[['L']]){
    private$popAquatic[names(larvalRatio),timeAq[['E']] + i] = (larvalEQ * larvalDist[i]) * larvalRatio
  } # end larva loop

  # set eggs
  for(i in timeAq[['E']]:1){
    private$popAquatic[ ,i] = private$popAquatic[ ,i+1] / epLife
  } # end egg loop

  # set pupae
  # need one more density-dependent death to get first pupae time
  tPupa <- timeAq[['E']] + timeAq[['L']]
  private$popAquatic[ ,tPupa + 1] <- private$popAquatic[ ,tPupa] * lLife
  # if pupa is more than 1 day, loop over rest
  if(timeAq[['P']] > 1){
    for(i in (tPupa + 1):(sum(timeAq)-1)){
      private$popAquatic[ ,i+1] = private$popAquatic[ ,i] * epLife
    } # end pupa loop
  }

  ##########
  # set male population breakdown
  ##########
  private$popMale[names(adultRatioM)] = adultRatioM * adultEQ/2

  ##########
  # set mated female population breakdown
  ##########
  # this isn't exactly correct, it mates all males to all females, ignoring
  #  the genotype-specific male mating abilities
  private$popUnmated[names(adultRatioF)] = adultRatioF * adultEQ/2
  private$popFemale = private$popUnmated %o% normalise(private$popMale)
  private$popUnmated[] = 0

}

#' Set Initial Population Deterministic
#'
#' Calls \code{\link{set_initialPopulation_Patch}} to initialize a steady-state
#' population distribution.
#'
#' @param adultEQ Equilibrium number of adults
#' @param larvalEQ Equilibrium number of larvae
#' @param adultRatioF Genotype specific ratio for adult females
#' @param adultRatioM Genotype specific ratio for adult males
#' @param larvalRatio Genotype specific ratio for larvae
#' @param timeAq Time for each aquatic stage
#' @param muAq Aquatic death rate
#' @param alpha Density-dependent centering parameter
#'
set_population_deterministic_Patch <- function(adultEQ = adultEQ, larvalEQ = larvalEQ,
                                         adultRatioF = adultRatioF, adultRatioM = adultRatioM,
                                         larvalRatio = larvalRatio,
                                         timeAq = timeAq, muAq = muAq, alpha = alpha){

  self$initialPopulation(adultEQ = adultEQ, larvalEQ = larvalEQ,
                         adultRatioF = adultRatioF, adultRatioM = adultRatioM,
                         larvalRatio = larvalRatio,
                         timeAq = timeAq, muAq = muAq, alpha = alpha)

}

#' Set Initial Population Stochastic
#'
#' Calls \code{\link{set_initialPopulation_Patch}} to initialize a steady-state
#' population distribution. Populations are then rounded to integer values.
#'
#' @param adultEQ Equilibrium number of adults
#' @param larvalEQ Equilibrium number of larvae
#' @param adultRatioF Genotype specific ratio for adult females
#' @param adultRatioM Genotype specific ratio for adult males
#' @param larvalRatio Genotype specific ratio for larvae
#' @param timeAq Time for each aquatic stage
#' @param muAq Aquatic death rate
#' @param alpha Density-dependent centering parameter
#'
set_population_stochastic_Patch <- function(adultEQ = adultEQ, larvalEQ = larvalEQ,
                                           adultRatioF = adultRatioF, adultRatioM = adultRatioM,
                                           larvalRatio = larvalRatio,
                                           timeAq = timeAq, muAq = muAq, alpha = alpha){

  # set initial population
  self$initialPopulation(adultEQ = adultEQ, larvalEQ = larvalEQ,
                         adultRatioF = adultRatioF, adultRatioM = adultRatioM,
                         larvalRatio = larvalRatio,
                         timeAq = timeAq, muAq = muAq, alpha = alpha)

  ##########
  # make everything an integer
  ##########
  private$popAquatic[] <- round(private$popAquatic)
  private$popMale[] <- round(private$popMale)
  private$popFemale[] <- round(private$popFemale)

}

#' Reset Patch to Initial Conditions
#'
#' Resets a patch to its initial configuration so that a new one does not have
#' to be created and allocated in the network (for Monte Carlo simulation).
#'
#' @param verbose Chatty? Default is TRUE
#'
reset_Patch <- function(verbose = TRUE){

  if(verbose){cat("reset patch ",private$patchID,"\n",sep="")}

  # reset population
  private$popAquatic[] = private$popAquatict0
  private$popMale[] = private$popMalet0
  private$popFemale[] = private$popFemalet0

  # Reset Mosquito Releases
  private$eggReleases = private$NetworkPointer$get_patchReleases(private$patchID,"Egg")
  private$maleReleases = private$NetworkPointer$get_patchReleases(private$patchID,"M")
  private$femaleReleases = private$NetworkPointer$get_patchReleases(private$patchID,"F")
  private$matedFemaleReleases = private$NetworkPointer$get_patchReleases(private$patchID,"mF")

}

#' Initialize Output from Focal Patch
#'
#' Writes output to the text connections specified in the enclosing \code{\link{Network}}.
#'
oneDay_initOutput_Patch <- function(){

  ##########
  # headers
  ##########
  if(private$patchID == 1){
    # males
    writeLines(text = paste0(c("Time","Patch",private$NetworkPointer$get_genotypesID()), collapse = ","),
           con = private$NetworkPointer$get_conADM(), sep = "\n")
    # females
    femaleCrosses = c(t(outer(private$NetworkPointer$get_genotypesID(),private$NetworkPointer$get_genotypesID(),FUN = paste0)))
    writeLines(text = paste0(c("Time","Patch",femaleCrosses), collapse = ","),
               con = private$NetworkPointer$get_conADF(),sep = "\n")
  }

  ##########
  # males
  ##########
  writeLines(text = paste0(c(0,private$patchID,private$popMale),collapse = ","),
             con = private$NetworkPointer$get_conADM(), sep = "\n")

  ##########
  # females
  ##########
  writeLines(text = paste0(c(0,private$patchID,c(t(private$popFemale))),collapse = ","),
             con = private$NetworkPointer$get_conADF(), sep = "\n")

}

#' Write Output from Focal Patch
#'
#' Writes output to the text connections specified in the enclosing \code{\link{Network}}.
#'
oneDay_writeOutput_Patch <- function(){

  tNow = private$NetworkPointer$get_tNow()

  # write males
  ADMout = paste0(c(tNow,private$patchID,private$popMale),collapse = ",")
  writeLines(text = ADMout,con = private$NetworkPointer$get_conADM(),sep = "\n")

  # write females
  ADFout = paste0(c(tNow,private$patchID,c(t(private$popFemale))),collapse = ",")
  writeLines(text = ADFout,con = private$NetworkPointer$get_conADF(),sep = "\n")

}
