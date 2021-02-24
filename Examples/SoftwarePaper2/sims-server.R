# --------------------------------------------------------------------------------
#
#   Simulations for MGDrivE2 paper
#   Comoros
#   September 2020
#   Sean L. Wu (slwu89@berkeley.edu)
#
# --------------------------------------------------------------------------------

#library(here)
library(MGDrivE)
library(MGDrivE2)
rm(list=ls());gc()

# average PfPR from MAP
COM_pfpr <- 0.3641667

# load external data
dataDir <- "~/MGDrivE/Main/SoftwarePaper2"
source(file.path(dataDir,"Data/cube-revisions.R"))
source(file.path(dataDir,"Data/hazards-revisions.R"))
carry <- read.csv2(file =file.path(dataDir,"Data/K_km_v3.csv"),sep = ",",stringsAsFactors = FALSE)
ad <- read.csv2(file = file.path(dataDir,"Data/km_muEM.csv"),sep = ",",stringsAsFactors = FALSE)


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
  qEIP = 1/11,
  nEIP = 6
)

parameters$a <- parameters$f * parameters$Q
parameters$X <- COM_pfpr
parameters$NH <- 350998


# --------------------------------------------------------------------------------
#   time varying parameters
# --------------------------------------------------------------------------------

# carrying capacity
carry$Grande_Comore <- as.numeric(carry$Grande_Comore)
carry$Moheli <- as.numeric(carry$Moheli)
carry$Anjouan <- as.numeric(carry$Anjouan)

# adult death
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

# hazard vector
hazards <- make_hazards(
  spn_P = SPN_P,spn_T = SPN_T,cube = cube,
  par = parameters,log_dd = TRUE,exact = TRUE,verbose = FALSE
)


# --------------------------------------------------------------------------------
#   release strategy
# --------------------------------------------------------------------------------

# releases
r_times <- seq(from = 365*3, length.out = 8, by = 7)
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
dt_stoch <- 0.01

# parallel options
# parameters for stochastic simulations
num_core <- 25

# total runs
n <- 100

# main output folders
# main_out <- "mgdrive2_paper"
# analysis_out <- c("raw", "traces", "analyzed")
# analysis_folders <- file.path(main_out, analysis_out)
main_out <- "/RAID0/mgdrive2_paper"
analysis_out <- c("raw", "traces", "analyzed")
analysis_folders <- file.path(main_out, analysis_out)
ltsDir <- "/RAID5/marshallShare/mgdrive2_paper"

# setup the cluster
#computer = "windows"
computer = "notWindows"
if(computer == "windows"){
  # windows can't run fork clusters
  #  it has to use sockets
  #
  # be very careful, this causing copying of some sort, and it eats up a ton of ram
  cl <- parallel::makePSOCKcluster(names = num_core)

  # load MGDrivE2 on each socket
  parallel::clusterEvalQ(cl = cl, expr = {library(MGDrivE2)})

  # export required objects to each socket
  parallel::clusterExport(
    cl = cl,
    varlist = c("M0", "tmax", "dt", "dt_stoch", "S", "Sout","hazards", "events", "SPN_P", "main_out","analysis_out","analysis_folders")
  )

} else {
  # *nix systems
  # no copying, no memory issues!
  cl <- parallel::makeForkCluster(nnodes = num_core)

}

# set parallel seed
parallel::clusterSetRNGStream(cl = cl, iseed = 682394L)

# run funcs
parallel::clusterApplyLB(cl = cl, x = 1:n, fun = function(x){

  # build analysis folders
  rep_out <- formatC(x = x, width=3, format='d', flag='0')
  for(i in analysis_folders){ dir.create(path = i, recursive = TRUE) }

  # build repetition folders
  rep_folders <- file.path(analysis_folders[1], rep_out)
  for(i in rep_folders){ dir.create(i) }

  # run sims
  sim_trajectory_CSV(
    x0 = M0, tmax = tmax, dt = dt,
    dt_stoch = dt_stoch, folders = rep_folders, sampler = "tau",
    stage = c("M", "F", "H"), S = S, Sout = Sout,
    hazards = hazards, events = events, verbose = FALSE, maxhaz = 1e12
  )

  # split everything by patch, aggregate by genotype
  split_aggregate_CSV(
    read_dir = analysis_folders[1], write_dir = analysis_folders[2],
    spn_P = SPN_P, tmax = tmax, dt = dt, verbose = FALSE, sum_fem = TRUE
  )

})

# stop cluster
parallel::stopCluster(cl)

# mean and 95% quantiles
summarize_stats_CSV(
  read_dir = analysis_folders[2], write_dir = analysis_folders[3],
  spn_P = SPN_P, tmax = tmax, dt = dt, mean = TRUE,
  quantiles = c(0.025, 0.975), verbose = FALSE
)

# move data to longTermStorage
if(!dir.exists(ltsDir)) dir.create(path = ltsDir, recursive = TRUE)

# copy analyzed and trace dirs
file.copy(from = analysis_folders[c(2,3)], to = ltsDir, overwrite = TRUE,
          recursive = TRUE, copy.mode = FALSE, copy.date = TRUE)

# put the incidence into the right folder
inc_folder <- file.path(paste0(ltsDir,"/epitraces"))
if(!dir.exists(inc_folder)){
  dir.create(inc_folder)
}

rep_folders <- list.files(analysis_folders[1])
for(r in seq_along(rep_folders)){
  path <- file.path(analysis_folders[1],rep_folders[r])
  file.copy(
    from = file.path(path,"Tracking.csv"),
    to = file.path(paste0(inc_folder,"/incidence_",rep_folders[r],".csv")),
    overwrite = TRUE,recursive = FALSE,copy.mode = FALSE,copy.date = TRUE
  )
}

# compress data for storage
for(i in 1:length(analysis_out)){

  # build terminal command
  # -C change to the following directory
  # -c output to this place (sent to stdout here)
  # -f read from this file
  # | lbzip2 pipe stdin to lbzip2
  # -9 lbzip2 level 9
  # -n number of cores to use
  # > output to here
  cmd <- file.path('-C ', main_out, ' --remove-files -cf - ', analysis_out[i],
                   ' | lbzip2 -9 -n ', num_core, ' > ', analysis_folders[i], '.tar.bz2', fsep = '')

  # store
  system2(command = 'tar', args = cmd )
} # end storage loop
