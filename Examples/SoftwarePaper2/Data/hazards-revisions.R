
# --------------------------------------------------------------------------------
#   time varying hazards
# --------------------------------------------------------------------------------

make_female_mort_haz_inhom <- function(trans, u, cube, params,
                                       exact = TRUE, tol = 1e-8){
  
  # mortality is a time-dependent hazard
  muF <- params$muF
  if(typeof(muF) != "closure"){
    stop("Inhomogeneous hazard 'make_female_mort_haz_inhom', ",
         "'muF' in 'params' list needs to be a function")
  }
  
  # which places have input arcs to this transition
  s <- trans$s
  
  # weights of those arcs
  w <- trans$s_w
  
  # omega is dependent on genotype
  f_gen <- strsplit(x = u[s], split = "_", fixed = TRUE)[[1]][2]
  omega <- cube$omega[f_gen]
  
  # return the hazard function
  if(exact){
    
    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      function(t,M){
        if(w <= M[s]){
          # return(mu_ad_mean * omega * M[s])
          return(muF(t) * omega * M[s])
        } else {
          return(0)
        }
      }
    )
    
  } else {
    
    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        # haz <- mu_ad_mean * omega * M[s]
        haz <- muF(t) * omega * M[s]
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

make_male_mort_haz_inhom <- function(trans, u, cube, params,
                                     exact = TRUE, tol = 1e-8){
  
  # mortality is a time-dependent hazard
  muM <- params$muM
  if(typeof(muM) != "closure"){
    stop("Inhomogeneous hazard 'make_male_mort_haz_inhom', ",
         "value 'muM' in 'params' list needs to be a function")
  }
  
  # which places have input arcs to this transition
  s <- trans$s
  
  # weights of those arcs
  w <- trans$s_w
  
  # omega is dependent on genotype
  m_gen <- strsplit(x = u[s], split = "_", fixed = TRUE)[[1]][2]
  omega <- cube$omega[m_gen]
  
  # return the hazard function
  if(exact){
    
    # EXACT hazards (check enabling degree: for discrete simulation only)
    return(
      function(t,M){
        if(w <= M[s]){
          # return(mu_ad_mean * omega * M[s])
          return(muM(t) * omega * M[s])
        } else {
          return(0)
        }
      }
    )
    
  } else {
    
    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        # haz <- mu_ad_mean * omega * M[s]
        haz <- muM(t) * omega * M[s]
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

make_larvae_mort_haz_log_inhom <- function(trans,u,l_ix,node,cube,params,exact = TRUE,tol = 1e-8){
  
  # rate constants
  muL <- params$muL
  K <- params$K[[node]]
  if(typeof(K) != "closure"){
    stop("Inhomogeneous hazard 'make_larvae_mort_haz_log', ",
         "value 'K' in 'params' list needs to be a function")
  }
  
  
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
          return(muL*(1 + (L/K(t)))*M[s])
          # return(muL*(1 + (L/K_mean))*M[s])
        } else {
          return(0)
        }
      }
    )
    
  } else {
    
    # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
    return(
      function(t,M){
        # get total males
        L <- sum(M[l_ix])
        haz <- muL*(1 + (L/K(t)))*M[s]
        # haz <- muL*(1 + (L/K_mean))*M[s]
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

# make hazards by hand
make_hazards <- function(spn_P,spn_T,cube,par,log_dd=TRUE,exact=TRUE,tol=1e-12,verbose=TRUE){
  
  if(tol > 1e-6 & !exact){
    cat("warning: hazard function tolerance ",tol," is large; consider tolerance < 1e-6 for sufficient accuracy\n")
  }
  
  if(log_dd){
    if(!("K" %in% names(par))){
      stop("if using logistic (carrying capacity) based density-dependent larval mortality, please specify parameter 'K' in par")
    }
  } else {
    if(!("gamma" %in% names(par))){
      stop("if using Lotka-Volterra based density-dependent larval mortality, please specify parameter 'gamma' in par")
    }
  }
  
  # transitions and places
  v <- spn_T$v
  u <- spn_P$u
  
  n <- length(v)
  if(verbose){
    pb <- txtProgressBar(min = 1,max = n,style = 3)
    pp <- 1
  }
  
  # the hazard functions
  h <- vector("list",n)
  h <- setNames(h,v)
  
  # get male and larvae indices
  l_ix <- as.vector(spn_P$ix[[1]]$larvae)
  m_ix <- spn_P$ix[[1]]$males
  
  # human indices
  h_ix <- spn_P$ix[[1]]$humans
  
  cat(" --- generating hazard functions for SPN --- \n")
  
  # make the hazards
  for(t in 1:n){
    
    type <- spn_T$T[[t]]$class
    
    # make the correct type of hazard
    
    # MOSQUITO HAZARDS
    if(type == "oviposit"){
      h[[t]] <- MGDrivE2:::make_oviposit_haz(t = spn_T$T[[t]],u = u,cube = cube,par = par,exact = exact,tol = tol)
    } else if(type == "egg_adv"){
      h[[t]] <- MGDrivE2:::make_egg_adv_haz(t = spn_T$T[[t]],u = u,cube = cube,par = par,exact = exact,tol = tol)
    } else if(type == "egg_mort"){
      h[[t]] <- MGDrivE2:::make_egg_mort_haz(t = spn_T$T[[t]],u = u,cube = cube,par = par,exact = exact,tol = tol)
    } else if(type == "larvae_adv"){
      h[[t]] <- MGDrivE2:::make_larvae_adv_haz(t = spn_T$T[[t]],u = u,cube = cube,par = par,exact = exact,tol = tol)
      # INHOMOGENEOUS
    } else if(type == "larvae_mort"){
      h[[t]] <- make_larvae_mort_haz_log_inhom(t = spn_T$T[[t]],u = u,l_ix = l_ix,node=1,cube = cube,par = par,exact = exact,tol = tol)
    } else if(type == "pupae_adv"){
      h[[t]] <- MGDrivE2:::make_pupae_adv_haz(t = spn_T$T[[t]],u = u,cube = cube,par = par,exact = exact,tol = tol)
    } else if(type == "pupae_mort"){
      h[[t]] <- MGDrivE2:::make_pupae_mort_haz(t = spn_T$T[[t]],u = u,cube = cube,par = par,exact = exact,tol = tol)
    } else if(type == "pupae_2m"){
      h[[t]] <- MGDrivE2:::make_pupae_2male_haz(t = spn_T$T[[t]],u = u,cube = cube,par = par,exact = exact,tol = tol)
    } else if(type == "pupae_2f"){
      h[[t]] <- MGDrivE2:::make_pupae_2female_haz(t = spn_T$T[[t]],u = u,m_ix = m_ix,cube = cube,par = par,exact = exact,tol = tol)
    } else if(type == "pupae_2unmated"){
      h[[t]] <- MGDrivE2:::make_pupae_2unmated_haz(t = spn_T$T[[t]],u = u,m_ix = m_ix,cube = cube,par = par,exact = exact,tol = tol)
    } else if(type == "female_unmated_mate"){
      h[[t]] <- MGDrivE2:::make_unmated_2female_haz(t = spn_T$T[[t]],u = u,m_ix = m_ix,cube = cube,par = par,exact = exact,tol = tol)
      # INHOMOGENEOUS
    } else if(type == "male_mort"){
      h[[t]] <- make_male_mort_haz_inhom(t = spn_T$T[[t]],u = u,cube = cube,par = par,exact = exact,tol = tol)
      # INHOMOGENEOUS
    } else if(type %in% c("female_mort","female_unmated_mort")){
      h[[t]] <- make_female_mort_haz_inhom(t = spn_T$T[[t]],u = u,cube = cube,par = par,exact = exact,tol = tol)
    } else if(type == "female_inf"){
      h[[t]] <- MGDrivE2:::make_female_inf_epi_haz(t = spn_T$T[[t]],u = u,h_ix = h_ix,cube = cube,par = par,exact = exact,tol = tol)
    } else if(type == "female_eip"){
      h[[t]] <- MGDrivE2:::make_female_eip_epi_haz(t = spn_T$T[[t]],u = u,par = par,exact = exact,tol = tol)
    } else if(type == "female_inc"){
      # can reuse above hazard because transition hazard is the same
      h[[t]] <- MGDrivE2:::make_female_eip_epi_haz(t = spn_T$T[[t]],u = u,par = par,exact = exact,tol = tol)
      # HUMAN HAZARDS
    } else if(type == "H_birth"){
      h[[t]] <- MGDrivE2:::make_human_birth_sis_haz(t = spn_T$T[[t]],u = u,par = par,exact = exact,tol = tol)
    } else if(type == "H_mort"){
      h[[t]] <- MGDrivE2:::make_human_death_sis_haz(t = spn_T$T[[t]],u = u,par = par,exact = exact,tol = tol)
    } else if(type == "H_infection"){
      h[[t]] <- MGDrivE2:::make_human_inf_sis_haz(t = spn_T$T[[t]],u = u,h_ix = h_ix,cube = cube,par = par,exact = exact,tol = tol)
    } else if(type == "H_recovery"){
      h[[t]] <- MGDrivE2:::make_human_rec_sis_haz(t = spn_T$T[[t]],u = u,par = par,exact = exact,tol = tol)
    } else {
      stop(paste0("error in making hazard function for unknown class type: ",type))
    }
    
    if(verbose){setTxtProgressBar(pb,t)}
  }
  
  if(verbose){close(pb)}
  
  cat(" --- done generating hazard functions for SPN --- \n")
  
  return(list("hazards"=h,"flag"=exact))
}
