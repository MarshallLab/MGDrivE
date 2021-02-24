
CUTTINGF <- 0.788
CUTTINGM <- 0.4805
HGENF <- HOMING_RATE_FEMALE <- 0.7978
HGENM <- HOMING_RATE_MALE <- 0.7688

RGEN <- 1/3*1/2

CCOST <- 0.078
#FIT_COST=c(0.05,0.10)/11.11111 # 5% and 10% lifetime cost on gRNA and B alleles
FC <- FIT_COST <- c(0.11)/11.11111 #seq.int(from = 0.01, to = 0.15, by = 0.01)/11.11111 # fitness sweep.

# genotype keys:
# C: Cas9 allele
# H: gRNA/payload allele
# W: wild type
# B: broken allele
# R: functional but resistant

cube <- MGDrivE::cubeSplitDrive(
  cM = CUTTINGM, cF = CUTTINGF,
  chM = HGENM, crM = RGEN,
  chF = HGENF, crF = RGEN,
  s = c(
    "WCWW"=1-CCOST,"WCWH"=1-CCOST,"WCWR"=1-CCOST,"WCWB"=1-CCOST,"WCHH"=1-CCOST,
    "WCHR"=1-CCOST,"WCHB"=1-CCOST,"WCRR"=1-CCOST,"WCRB"=1-CCOST,"WCBB"=1-CCOST,
    "CCWW"=1-2*CCOST,"CCWH"=1-2*CCOST,"CCWR"=1-2*CCOST,"CCWB"=1-2*CCOST,"CCHH"=1-2*CCOST,
    "CCHR"=1-2*CCOST,"CCHB"=1-2*CCOST,"CCRR"=1-2*CCOST,"CCRB"=1-2*CCOST,"CCBB"=1-2*CCOST
    ),
  omega = c(
    'WWHH'=1-FC, 'WWHB'=1-FC,
    'WWBB'=1-FC, 'WCHH'=1-FC, 'WCHB'=1-FC,
    'WCBB'=1-FC, 'CCHH'=1-FC, 'CCHB'=1-FC,
    'CCBB'=1-FC
    )
)