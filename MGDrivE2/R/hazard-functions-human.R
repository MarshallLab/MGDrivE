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
# HUMAN DEMOGRAPHY
################################################################################

# human birth
make_human_birth_sis_haz <- function(trans,u,params,exact=TRUE,tol=1e-8){

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # rate constants
  muH <- params$muH


  # hazard fn
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      function(t,M){
        if(w <= M[s]){
          return(muH * M[s])
        } else {
          return(0)
        }
      }
    )

  } else {

    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        haz <- muH * M[s]
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

# human death
make_human_death_sis_haz <- function(trans,u,params,exact=TRUE,tol=1e-8){

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # rate constants
  muH <- params$muH


  # hazard fn
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      function(t,M){
        if(w <= M[s]){
          return(muH * M[s])
        } else {
          return(0)
        }
      }
    )

  } else {

    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        haz <- muH * M[s]
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
# HUMAN INFECTION
################################################################################

# infection
make_human_inf_sis_haz <- function(trans,u,h_ix,cube,params,exact=TRUE,tol=1e-8){

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # c is dependent on genotype
  f_gen <- strsplit(x = u[s[1]],split = "_",fixed = TRUE)[[1]][2]

  # rate constants
  a <- params$a
  b <- cube$b[f_gen]

  # assign here so that each newly generated closure has the right indices
  h_ix <- h_ix

  # safety check
  if(check_double(b)){
    stop("b is missing from cube list; called from 'make_human_inf_sis_haz'")
  }

  # hazard fn
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      function(t,M){
        if(all(w <= M[s])){
          NH <- sum(M[h_ix])
          return(a * b * M[s[2]] * (1/NH) * M[s[1]])
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
          haz = 0
        } else {
          haz = a * b * M[s[2]] * (1/NH) * M[s[1]]
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

# recovery
make_human_rec_sis_haz <- function(trans,u,params,exact=TRUE,tol=1e-8){

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # rate constants
  r <- params$r


  # hazard fn
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      function(t,M){
        if(w <= M[s]){
          return(r * M[s])
        } else {
          return(0)
        }
      }
    )

  } else {

    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        haz <- r * M[s]
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
# HUMAN LATENT PERIOD
################################################################################

# latent, SEIR
make_human_latent_seir_haz <- function(trans,u,params,exact=TRUE,tol=1e-8){

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # rate constants
  delta <- params$delta

  # hazard fn
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      function(t,M){
        if(w <= M[s]){
          return(delta * M[s])
        } else {
          return(0)
        }
      }
    )

  } else {

    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        haz <- delta * M[s]
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
# HUMAN MOVEMENT HAZARD FUNCTIONS
################################################################################

# human movement
make_human_move_sis_haz <- function(trans,u,params,exact=TRUE,tol=1e-8){

  # which places have input arcs to this transition
  s <- trans$s

  # weights of those arcs
  w <- trans$s_w

  # get the origin and destination
  #  do we need the "tail"?
  #  should follow "infection_H_label", so the label is always teh 3rd element
  od_pair <- as.integer(vapply(X = strsplit(x = trans$label, split = "->", fixed = TRUE)[[1]],
                    FUN = function(x){tail(x = strsplit(x = x, split = "_", fixed = TRUE)[[1]], 1)},
                    FUN.VALUE = character(length = 1)))

  # movement parameters
  move_r <- params$human_move_rates[od_pair[1]]
  move_p <- params$human_move_probs[od_pair[1],od_pair[2]]

  # safety check
  if(check_double(move_r) || check_double(move_p)){
    stop("human_move_rates or human_move_probs missing from parameters list; called from 'make_human_move_haz'")
  }

  # return the hazard function
  if(exact){

    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      # rate of move events firing is: rate to jump out of the source * prob to go from source to dest * number of males who can make the jump
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
