# --------------------------------------------------------------------------------
#
#   Simulations for MGDrivE2 paper
#   Paper Revisions (January 2021)
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
source(here::here("Data/cube-revisions.R"))

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

# exaggerated K
ex_factor <- 2
K_diff <- diff(K_ts)
K_ts_ex <- cumsum(c(K_ts[1],K_diff*ex_factor))
K_ts_ex_adj <- K_ts_ex * (mean(K_ts) / mean(K_ts_ex))

step_K <- stats::stepfun(x = carry$Day,y = c(K_ts_ex_adj[1],K_ts_ex_adj),f = 0,right = FALSE)
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

source(here::here("Data/hazards-revisions.R"))

# hazard vector
hazards <- make_hazards(
  spn_P = SPN_P,spn_T = SPN_T,cube = cube,
  par = parameters,log_dd = TRUE,exact = TRUE,verbose = TRUE
)


# --------------------------------------------------------------------------------
#   release strategy
# --------------------------------------------------------------------------------

# release scheme from https://elifesciences.org/articles/51701#s2 
r_times <- seq(from = 365*3, length.out = 5, by = 7)
r_size <- NF*2 # corresponds to a 1:1 release scheme
r_size <- 50000

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

# set.seed(3213L)

system.time(
  sim_out <- sim_trajectory_R(
    x0 = M0, tmax = tmax, dt = dt, S = S, Sout = Sout,
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


# --------------------------------------------------------------------------------
#   sweep release ratios
# --------------------------------------------------------------------------------

library(foreach)
library(iterators)
library(doParallel)

release_fracs <- seq(from=1,to=2-0.05,by=0.1)

cl <- parallel::makeCluster(10,outfile = "/home/slwu89/Desktop/outfile.txt")
doParallel::registerDoParallel(cl)

parallel::clusterSetRNGStream(cl = cl,iseed = 5837192L)

system.time(sweep_out <- foreach(i = 1:length(release_fracs), .export = c("cube","r_times","r_size","M0","tmax","dt","S","hazards","SPN_P"), .packages = "MGDrivE2") %dopar% {
  
  r_size <- NF*release_fracs[i] # corresponds to a 1:1 release scheme
  
  events <- data.frame(
    "var" = paste0("M_", cube$releaseType),
    "time" = r_times,
    "value" = r_size,
    "method" = "add",
    stringsAsFactors = FALSE
  )
  
  sim_out <- sim_trajectory_R(
    x0 = M0, tmax = tmax, dt = dt, S = S, Sout = NULL,
    hazards = hazards, sampler = "tau", dt_stoch = 0.05,
    events = events, verbose = T,maxhaz = 1e12
  )
  
  ret <- list()
  ret$frac <- release_fracs[i]
  ret$sim_female <- summarize_females_epi(out = sim_out$state, spn_P = SPN_P)
  ret$sim_male <- summarize_males(out = sim_out$state)
  ret$sim_humans <- summarize_humans_epiSIS(out = sim_out$state)
  
  return(ret)
  
})

parallel::stopCluster(cl)

library(gridExtra)

lapply(X = sweep_out,FUN = function(xx){
  
  p1 <- ggplot(data = xx$sim_female) +
    geom_line(aes(x=time,y=value,color=genotype)) +
    facet_wrap(. ~ inf, scales = "free_y") +
    theme_bw()

  p3 <- ggplot(data = xx$sim_humans) +
    geom_line(aes(x=time,y=value,color=inf)) +
    theme_bw()

  pp <- grid.arrange(p1,p3,nrow=2)
  
  ggsave(filename = paste0(here::here(),"/Figs/sweep",xx$frac,".pdf"),plot = pp,device = "pdf",width = 16,height = 18)
})
