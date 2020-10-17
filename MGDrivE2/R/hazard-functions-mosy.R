################################################################################
#
#   MGDrivE2: hazard functions
#   Marshall Lab
#   Sean L. Wu (slwu89@berkeley.edu)
#   October 2019
#   Jared Bennett (jared_bennett@berkeley.edu)
#   April 2020
#
################################################################################

################################################################################
# OVIPOSITION
################################################################################

# make an oviposition hazard function
# the functions these make need to be stored in the same order as v
# so we put the transition at the place in v[t], where t comes from oviposit[i,j,k]
# i: row index into oviposit transition array
# j: col index into oviposit transition array
# k: slice index into oviposit transition array
# exact: check enabling (if not, always calc hazard and truncate to zero at small values)
make_oviposit_haz <- function(trans,u,cube,params,exact = TRUE,tol = 1e-8){

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # get the genotype-specific bits for this transition
  input <- strsplit(x = u[s], split = "_", fixed = TRUE)[[1]]
  output <- strsplit(x = u[trans$o[2]], split = "_", fixed = TRUE)[[1]]


  f_gen <- input[2]
  m_gen <- input[3]
  o_gen <- output[2]

  # offspring genotype probability
  o_prob <- cube$ih[f_gen,m_gen,o_gen] * cube$tau[f_gen,m_gen,o_gen]

  # egg laying rate
  beta <- params$beta * cube$s[f_gen]

  # return the hazard function
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      function(t,M){
        # check the transition is enabled to fire
        if(w <= M[s]){
          return(o_prob * beta * M[s])
        } else {
          return(0)
        }
      }
    )

  } else {

    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        haz <- o_prob * beta * M[s]
        if(haz < tol){
          return(0)
        } else {
          return(haz)
        }
      }
    )

  }
  # end of function
}


################################################################################
# EGG TRANSITIONS
################################################################################

make_egg_mort_haz <- function(trans,u,cube,params,exact = TRUE,tol = 1e-8){

  # if these are time-varying, chuck them into the returned function
  muE <- params$muE

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # return the hazard function
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      function(t,M){
        if(w <= M[s]){
          return(muE * M[s])
        } else {
          return(0)
        }
      }
    )

  } else {

    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        haz <- muE * M[s]
        if(haz < tol){
          return(0)
        } else {
          return(haz)
        }
      }
    )

  }
  # end of function
}

# make egg advancement function
make_egg_adv_haz <- function(trans,u,cube,params,exact = TRUE,tol = 1e-8){

  # if these are time-varying, chuck them into the returned function
  # nE is not allowed to vary with time
  qE <- params$qE
  nE <- params$nE

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # return the hazard function
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      function(t,M){
        if(w <= M[s]){
          return(qE*nE*M[s])
        } else {
          return(0)
        }
      }
    )

  } else {

    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        haz <- qE*nE*M[s]
        if(haz < tol){
          return(0)
        } else {
          return(haz)
        }
      }
    )

  }
  # end of function
}


################################################################################
# LARVAE TRANSITIONS
################################################################################

# this one is ~special~ it needs the places of the larvae it uses to compute the hazard.
# logistic style hazard
make_larvae_mort_haz_log <- function(trans,u,l_ix,node,cube,params,exact = TRUE,tol = 1e-8){

  # rate constants
  muL <- params$muL
  K <- params$K[node]

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # assign here so that each newly generated closure has the right indices
  l_ix <- l_ix

  # return the hazard function
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      function(t,M){
        if(w <= M[s]){
          L <- sum(M[l_ix])
          return(muL*(1 + (L/K))*M[s])
        } else {
          return(0)
        }
      }
    )

  } else {

    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        # get total larvae
        L <- sum(M[l_ix])
        haz <- muL*(1 + (L/K))*M[s]
        # check and return
        if(haz < tol){
          return(0)
        } else {
          return(haz)
        }
      }
    )

  }
  # end of function
}

# lotka-volterra style hazard
make_larvae_mort_haz_lk <- function(trans,u,l_ix,node,cube,params,exact = TRUE,tol = 1e-8){

  # rate constants
  muL <- params$muL
  gamma <- params$gamma[node]

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # assign here so that each newly generated closure has the right indices
  l_ix <- l_ix

  # return the hazard function
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      function(t,M){
        if(w <= M[s]){
          L <- sum(M[l_ix])
          return((muL + (gamma*L))*M[s])
        } else {
          return(0)
        }
      }
    )

  } else {

    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        L <- sum(M[l_ix])
        haz <- (muL + (gamma*L))*M[s]
        if(haz < tol){
          return(0)
        } else {
          return(haz)
        }
      }
    )

  }
  # end of function
}

# larval advancement
make_larvae_adv_haz <- function(trans,u,cube,params,exact = TRUE,tol = 1e-8){

  # if these are time-varying, chuck them into the returned function
  # nE is not allowed to vary with time
  qL <- params$qL
  nL <- params$nL

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # return the hazard function
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      function(t,M){
        if(w <= M[s]){
          return(qL*nL*M[s])
        } else {
          return(0)
        }
      }
    )

  } else {

    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        haz <- qL*nL*M[s]
        if(haz < tol){
          return(0)
        } else {
          return(haz)
        }
      }
    )

  }
  # end of function
}


################################################################################
# PUPAE TRANSITIONS
################################################################################

make_pupae_mort_haz <- function(trans,u,cube,params,exact = TRUE,tol = 1e-8){

  # rate constants
  muP <- params$muP

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # return the hazard function
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      function(t,M){
        if(w <= M[s]){
          return(muP*M[s])
        } else {
          return(0)
        }
      }
    )

  } else {

    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        haz <- muP*M[s]
        if(haz < tol){
          return(0)
        } else {
          return(haz)
        }
      }
    )

  }
  # end of function
}

# pupae advancement
make_pupae_adv_haz <- function(trans,u,cube,params,exact = TRUE,tol = 1e-8){

  # if these are time-varying, chuck them into the returned function
  # nE is not allowed to vary with time
  qP <- params$qP
  nP <- params$nP

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # return the hazard function
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      function(t,M){
        if(w <= M[s]){
          return(qP*nP*M[s])
        } else {
          return(0)
        }
      }
    )

  } else {

    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        haz <- qP*nP*M[s]
        if(haz < tol){
          return(0)
        } else {
          return(haz)
        }
      }
    )

  }
  # end of function
}

# pupae emerge to male
make_pupae_2male_haz <- function(trans,u,cube,params,exact = TRUE,tol = 1e-8){

  # if these are time-varying, chuck them into the returned function
  # nE is not allowed to vary with time
  qP <- params$qP
  nP <- params$nP

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # phi is dependent on genotype
  p_gen <- strsplit(x = u[s], split = "_", fixed = TRUE)[[1]][2]
  phi <- cube$phi[p_gen]

  # xiM is also dependent on genotype
  xi <- cube$xiM[p_gen]

  # return the hazard function
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      function(t,M){
        if(w <= M[s]){
          return(qP*nP*(1 - phi)*xi*M[s])
        } else {
          return(0)
        }
      }
    )

  } else {

    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        haz <- qP*nP*(1 - phi)*xi*M[s]
        if(haz < tol){
          return(0)
        } else {
          return(haz)
        }
      }
    )

  }
  # end of function
}

# pupae emerge to female
make_pupae_2female_haz <- function(trans,u,m_ix,cube,params,exact = TRUE,tol = 1e-8){

  # assign here so that each newly generated closure has the right indices
  m_ix <- m_ix

  # if these are time-varying, chuck them into the returned function
  # nE is not allowed to vary with time
  qP <- params$qP
  nP <- params$nP

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # phi is dependent on genotype
  p_gen <- strsplit(x = u[s[1]], split = "_", fixed = TRUE)[[1]][2]
  phi <- cube$phi[p_gen]

  # xiF is also dependent on genotype
  xi <- cube$xiF[p_gen]

  # mating "weights"
  eta <- cube$eta[p_gen,]

  # need to know the index of the male genotype
  m_gen <- strsplit(x = u[s[2]], split = "_", fixed = TRUE)[[1]][2]
  j <- match(x = m_gen, table = cube$genotypesID)

  # return the hazard function
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      function(t,M){
        if(all(w <= M[s])){
          # mating propensity
          mate_p <- M[m_ix] * eta
          mate_p <- mate_p / sum(mate_p)
          return(qP*nP*phi*xi*mate_p[j]*M[s[1]])
        } else {
          return(0)
        }
      }
    )

  } else {

    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        if(sum(M[m_ix]) < tol){
          return(0)
        }
        # mating propensity
        mate_p <- M[m_ix] * eta
        mate_p <- normalize(mate_p)
        haz <- qP*nP*phi*xi*mate_p[j]*M[s[1]]
        if(haz < tol){
          return(0)
        } else {
          return(haz)
        }
      }
    )

  }
  # end of function
}

# pupae emerge to unmated female
make_pupae_2unmated_haz <- function(trans,u,m_ix,cube,params,exact = TRUE,tol = 1e-8){

  # assign here so that each newly generated closure has the right indices
  m_ix <- m_ix

  # if these are time-varying, chuck them into the returned function
  # nE is not allowed to vary with time
  qP <- params$qP
  nP <- params$nP

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # phi is dependent on genotype
  p_gen <- strsplit(x = u[s], split = "_", fixed = TRUE)[[1]][2]
  phi <- cube$phi[p_gen]

  # xiF is also dependent on genotype
  xi <- cube$xiF[p_gen]

  # return the hazard function
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      function(t,M){
        if((sum(M[m_ix]) == 0) & (w <= M[s])){
          return(qP*nP*phi*xi*M[s])
        } else {
          return(0)
        }
      }
    )

  } else {

    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        if(sum(M[m_ix]) > tol){
          return(0)
        }
        haz <- qP*nP*phi*xi*M[s]
        if(haz < tol){
          return(0)
        } else {
          return(haz)
        }

      }
    )

  }
  # end of function
}


################################################################################
# MALE TRANSITIONS
################################################################################

make_male_mort_haz <- function(trans,u,cube,params,exact = TRUE,tol = 1e-8){

  # rate constants
  muM <- params$muM

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # omega is dependent on genotype
  m_gen <- strsplit(x = u[s], split = "_", fixed = TRUE)[[1]][2]
  omega <- cube$omega[m_gen]
  omega <- ifelse(omega == 0,1e3,1/omega)

  # return the hazard function
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      function(t,M){
        if(w <= M[s]){
          return(muM * omega * M[s])
        } else {
          return(0)
        }
      }
    )

  } else {

    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        haz <- muM * omega * M[s]
        if(haz < tol){
          return(0)
        } else {
          return(haz)
        }
      }
    )

  }
  # end of function
}


################################################################################
# FEMALE TRANSITIONS
################################################################################

make_female_mort_haz <- function(trans,u,cube,params,exact = TRUE,tol = 1e-8){

  # rate constants
  muF <- params$muF

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # omega is dependent on genotype
  f_gen <- strsplit(x = u[s], split = "_", fixed = TRUE)[[1]][2]
  omega <- cube$omega[f_gen]
  omega <- ifelse(omega == 0,1e3,1/omega)

  # return the hazard function
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      function(t,M){
        if(w <= M[s]){
          return(muF * omega * M[s])
        } else {
          return(0)
        }
      }
    )

  } else {

    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        haz <- muF * omega * M[s]
        if(haz < tol){
          return(0)
        } else {
          return(haz)
        }
      }
    )

  }
  # end of function
}

# unmated females mate
make_unmated_2female_haz <- function(trans,u,m_ix,cube,params,exact = TRUE,tol = 1e-8){

  # assign here so that each newly generated closure has the right indices
  m_ix <- m_ix

  nu <- params$nu

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # mating "weights"
  f_gen <- strsplit(x = u[s[1]], split = "_", fixed = TRUE)[[1]][2]
  eta <- cube$eta[f_gen,]

  # need to know the index of the male genotype
  m_gen <- strsplit(x = u[s[2]], split = "_", fixed = TRUE)[[1]][2]
  j <- match(x = m_gen, table = cube$genotypesID)

  # return the hazard function
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      function(t,M){
        if(all(w <= M[s])){
          # mating propensity
          mate_p <- M[m_ix] * eta
          mate_p <- mate_p / sum(mate_p)
          return(nu*mate_p[j]*M[s[1]])
        } else {
          return(0)
        }
      }
    )

  } else {

    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        if(sum(M[m_ix]) < tol){
          return(0)
        }
        # mating propensity
        mate_p <- M[m_ix] * eta
        mate_p <- normalize(mate_p)
        haz <- nu*mate_p[j]*M[s[1]]
        if(haz < tol){
          return(0)
        } else {
          return(haz)
        }
      }
    )

  }
  # end of function
}


################################################################################
# FEMALE MOSQUITO INFECTION HAZARD
################################################################################

make_female_inf_epi_haz <- function(trans,u,h_ix,cube,params,exact = TRUE,tol = 1e-8){

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # c is dependent on genotype
  f_gen <- strsplit(x = u[s[1]],split = "_",fixed = TRUE)[[1]][2]

  # rate constants
  a <- params$a
  c <- cube$c[f_gen]

  # assign here so that each newly generated closure has the right indices
  h_ix <- h_ix

  # safety check
  if(check_double(c)){
    stop("c is missing from cube list; called from 'make_female_inf_epi_haz'")
  }

  # hazard fn
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      function(t,M){
        if(all(w <= M[s])){
          NH <- sum(M[h_ix])
          return(a*c*(M[s[2]]/NH)*M[s[1]])
        } else {
          return(0)
        }
      }
    )

  } else {

    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        # set small numbers to 0
        # magic number is sqrt(.Machine$double.eps)
        M[h_ix[(M[h_ix] <= 1.490116e-08)]] <- 0
        # get total males
        NH <- sum(M[h_ix])
        # safety check
        if(NH <= 1.490116e-08){
          haz <- 0
        } else {
          haz <- a*c*(M[s[2]]/NH)*M[s[1]]
        }
        # return
        if(haz < tol){
          return(0)
        } else {
          return(haz)
        }
      }
    )

  }
  # end of function
}


################################################################################
# FEMALE MOSQUITO INCUBATION (EIP) HAZARD
# can reuse for EIP -> infectious
################################################################################

make_female_eip_epi_haz <- function(trans,u,params,exact = TRUE,tol = 1e-8){

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # rate constants
  qn <- params$qEIP * params$nEIP


  # hazard fn
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      function(t,M){
        if(w <= M[s]){
          return(qn * M[s])
        } else {
          return(0)
        }
      }
    )

  } else {

    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        haz <- qn * M[s]
        if(haz < tol){
          return(0)
        } else {
          return(haz)
        }
      }
    )

  }
  # end of function
}


################################################################################
# MOVEMENT HAZARD FUNCTIONS
################################################################################

# female movement
make_move_female_haz <- function(trans,u,params,exact=TRUE,tol=1e-8){

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # get the origin and destination
  # do we know the length of this always? - NO, changes with infection vs mosquito only
  #  ie, F_geno_geno_num, means we want the 4th split every time.
  od_pair <- as.integer(lapply(X = strsplit(x = trans$label,split = "->",fixed = TRUE)[[1]],
                               function(x){tail(x = strsplit(x = x,split = "_",fixed = TRUE)[[1]],
                                                n = 1) }
                               )
                        )

  # movement parameters
  move_r <- params$mosquito_move_rates[od_pair[[1]] ]
  move_p <- params$mosquito_move_probs[od_pair[[1]],od_pair[[2]] ]

  # safety check
  if(check_double(move_r) || check_double(move_p)){
    stop(c("mosquito_move_rates or mosquito_move_probs missing from parameters",
    "list; called from 'make_move_female_haz'"))
  }

  # return the hazard function
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      # rate of move events firing is: rate to jump out of the source *
      #  prob to go from source to dest * number of males who can make the jump
      function(t,M){
        # check the transition is enabled to fire
        if(w <= M[s]){
          return(move_r * move_p * M[s])
        } else {
          return(0)
        }
      }
    )

  } else {

    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        haz <- move_r * move_p * M[s]
        if(haz < tol){
          return(0)
        } else {
          return(haz)
        }
      }
    )

  }
  # end of function
}

# female movement
make_move_male_haz <- function(trans,u,params,exact=TRUE,tol=1e-8){

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # get the origin and destination
  # do we know the length of this always?
  #  ie, M_geno_num, means we want the 3rd split every time.
  od_pair <- as.integer(lapply(X = strsplit(x = trans$label,split = "->",fixed = TRUE)[[1]],
                               function(x){tail(x = strsplit(x = x,split = "_",fixed = TRUE)[[1]],
                                                n = 1) }
                               )
                        )

  # movement parameters
  move_r <- params$mosquito_move_rates[od_pair[[1]] ]
  move_p <- params$mosquito_move_probs[od_pair[[1]],od_pair[[2]] ]

  # safety check
  if(check_double(move_r) || check_double(move_p)){
    stop(c("mosquito_move_rates or mosquito_move_probs missing from parameters",
    "list; called from 'make_move_male_haz'"))
  }

  # return the hazard function
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      # rate of move events firing is: rate to jump out of the source *
      #  prob to go from source to dest * number of males who can make the jump
      function(t,M){
        # check the transition is enabled to fire
        if(w <= M[s]){
          return(move_r * move_p * M[s])
        } else {
          return(0)
        }
      }
    )

  } else {

    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        haz <- move_r * move_p * M[s]
        if(haz < tol){
          return(0)
        } else {
          return(haz)
        }
      }
    )

  }
  # end of function
}
