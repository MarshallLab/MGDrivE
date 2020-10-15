# --------------------------------------------------------------------------------
#
#   Simulations for MGDrivE2 paper
#   Comoros
#   September 2020
#   Sean L. Wu (slwu89@berkeley.edu)
#
# --------------------------------------------------------------------------------

library(here)
library(ggplot2)
library(MGDrivE)
library(MGDrivE2)
rm(list=ls());gc()

# average PfPR from MAP
COM_pfpr <- 0.3641667

# fitted cube
source(here::here("Data/cube.R"))


# --------------------------------------------------------------------------------
#   mean/variance for Erlang aquatic stages
# --------------------------------------------------------------------------------

# calibrating parameters
CV <- c(E=0.2,L=0.3,P=0.2)

get_shape <- function(cv,q){
  mu <- 1/q
  n <- 1 / ((q^2) * ((mu*cv)^2))
  return(round(n))
}

tEgg <- 3
tLarva <- 7
tPupa <- 2

qE <- 1/tEgg
nE <- get_shape(cv = CV[["E"]],q = qE)

qL <- 1/tLarva
nL <- get_shape(cv = CV[["L"]],q = qL)

qP <- 1/tPupa
nP <- get_shape(cv = CV[["P"]],q = qP)

parameters <- list(
  # lifecycle
  qE = qE,
  nE = nE,
  qL = qL,
  nL = nL,
  qP = qP,
  nP = nP,
  beta = 32,
  nu = 1/(4/24),
  phi = 0.5,
  # epi
  f = 1/3,
  Q = 0.9,
  b = 0.55,
  c = 0.15,
  r = 1/200,
  muH = 1/(62*365),
  qEIP = 1/10,
  nEIP = 6
)

parameters$a <- parameters$f * parameters$Q
parameters$X <- COM_pfpr
parameters$NH <- 350998


# --------------------------------------------------------------------------------
#   time varying parameters
# --------------------------------------------------------------------------------

# carrying capacity
carry <- read.csv2(file = here::here("Data/K_km_v3.csv"),sep = ",",stringsAsFactors = FALSE)
carry$Grande_Comore <- as.numeric(carry$Grande_Comore)
carry$Moheli <- as.numeric(carry$Moheli)
carry$Anjouan <- as.numeric(carry$Anjouan)

# adult death
ad <- read.csv2(file = here::here("Data/km_muEM.csv"),sep = ",",stringsAsFactors = FALSE)
ad$Grande_Comore <- as.numeric(ad$Grande_Comore)
ad$Moheli <- as.numeric(ad$Moheli)
ad$Anjouan <- as.numeric(ad$Anjouan)
ad$All <- as.numeric(ad$All)
ad$Day <- seq_along(ad$Grande_Comore)/24

# mean adult death for equilibrium
mu_ad_mean <- mean(ad$Grande_Comore)
parameters$muF <- mu_ad_mean
parameters$muM <- mu_ad_mean

# calculate equilibrium (with mean values)
ad_F_eq <- MGDrivE2:::base_female_SIS(params = parameters)
NF <- sum(ad_F_eq)

mu_aqua <- MGDrivE2::solve_muAqua(params = parameters,rm = 1.096)
parameters$muE <- parameters$muL <- parameters$muP <- mu_aqua

mosy_eq <- MGDrivE2:::basic_eq_life(params = parameters,NF = NF,phi = 0.5,log_dd = TRUE)

# set up time varying carrying capacity
K_mean <- mean(carry$Grande_Comore)
K_eq <- mosy_eq$params$K
adjustment <- K_eq / K_mean

K_ts <- carry$Grande_Comore
K_ts <- K_ts * adjustment

step_K <- stats::stepfun(x = carry$Day,y = c(K_ts[1],K_ts),f = 0,right = FALSE)
parameters$K <- c(step_K)

# set up time varying mortality
step_mort <- stats::stepfun(x = ad$Day,y = c(ad$Grande_Comore[1],ad$Grande_Comore),f = 0,right = FALSE)
parameters$muF <- parameters$muM <- step_mort


# --------------------------------------------------------------------------------
#   set up PN and simulation objects
# --------------------------------------------------------------------------------

# transmission paramters (b: mosy -> human, c: human -> mosy)
cube$c <- setNames(object = rep(x = parameters$c, times = cube$genotypesN), nm = cube$genotypesID)
cube$c[grep(pattern = "H",x = names(cube$c))] <- 0
cube$b <- setNames(object = rep(x = parameters$b, times = cube$genotypesN), nm = cube$genotypesID)
cube$b[grep(pattern = "H",x = names(cube$b))] <- 0

# PN
SPN_P <- spn_P_epiSIS_node(params = parameters, cube = cube)
SPN_T <- spn_T_epiSIS_node(spn_P = SPN_P, params = parameters, cube = cube)

S <- spn_S(spn_P = SPN_P, spn_T = SPN_T)

# initial condition
M0 <- setNames(object = numeric(length = length(SPN_P$u)), nm = SPN_P$u)

e_ix <- SPN_P$ix[[1]]$egg[,which(colnames(SPN_P$ix[[1]]$egg) == cube$wildType)]
l_ix <- SPN_P$ix[[1]]$larvae[,which(colnames(SPN_P$ix[[1]]$larvae) == cube$wildType)]
p_ix <- SPN_P$ix[[1]]$pupae[,which(colnames(SPN_P$ix[[1]]$pupae) == cube$wildType)]

wt_idx <- sapply(dimnames(SPN_P$ix[[1]]$females[,,1]),function(x){which(x == cube$wildType)})
f_ix  <- SPN_P$ix[[1]]$females[wt_idx[1],wt_idx[2],]

m_ix <- SPN_P$ix[[1]]$males[names(SPN_P$ix[[1]]$males) == cube$wildType]

M0[e_ix] <- mosy_eq$init[1, grep("E",names(mosy_eq$init[1,])) ]
M0[l_ix] <- mosy_eq$init[1, grep("L",names(mosy_eq$init[1,])) ]
M0[p_ix] <- mosy_eq$init[1, grep("P",names(mosy_eq$init[1,])) ]
M0[f_ix] <- ad_F_eq[1,]
M0[m_ix] <- mosy_eq$init[1, grep("M",names(mosy_eq$init[1,])) ]
M0[["H_S"]] <- parameters$NH * (1 - COM_pfpr)
M0[["H_I"]] <- parameters$NH * COM_pfpr


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

# hazard vector
hazards <- make_hazards(
  spn_P = SPN_P,spn_T = SPN_T,cube = cube,
  par = parameters,log_dd = TRUE,exact = TRUE,verbose = TRUE
)


# --------------------------------------------------------------------------------
#   release strategy
# --------------------------------------------------------------------------------

# releases
r_times <- seq(from = 365*3, length.out = 8, by = 7)
# r_size <- sum(mosy_eq$init[1,c("NM","NF")]) * 0.1
r_size <- 1e4

events <- data.frame(
  "var" = paste0("M_", cube$releaseType),
  "time" = r_times,
  "value" = r_size,
  "method" = "add",
  stringsAsFactors = FALSE
)

# tracking incidence
Sout <- track_hinf(spn_T = SPN_T, S = S)


# --------------------------------------------------------------------------------
#   simulation
# --------------------------------------------------------------------------------

# times
tmax <- min(tail(carry$Day,1),tail(ad$Day,1))
dt <- 1

set.seed(3213L)

system.time(
  sim_out <- sim_trajectory_R(
    x0 = M0, t0 = 0, tt = tmax, dt = dt, S = S, Sout = Sout,
    hazards = hazards, sampler = "tau", dt_stoch = 0.05,
    events = events, verbose = T,maxhaz = 1e12)
)

sim_female <- summarize_females_epi(out = sim_out$state, spn_P = SPN_P)
sim_male <- summarize_males(out = sim_out$state)
sim_humans <- summarize_humans_epiSIS(out = sim_out$state)

library(gridExtra)

p1 <- ggplot(data = sim_female) +
  geom_line(aes(x=time,y=value,color=genotype)) +
  facet_wrap(. ~ inf, scales = "free_y") +
  theme_bw()

p2 <- ggplot(data = sim_male) +
  geom_line(aes(x=time,y=value,color=genotype)) +
  theme_bw()

p3 <- ggplot(data = sim_humans) +
  geom_line(aes(x=time,y=value,color=inf)) +
  theme_bw()

grid.arrange(p1,p3,nrow=2)
