# john fitted params, Adriana Adolfi cube
cutting=1.0
pHDR_F=0.98
pHDR_M=1.0
pRES=0.170
sHH=0
sHB=0
sBB=0.998
fMC=0.937
prMR=pRES
# fitted parameters implemented in the basic CRISPR cube
cube <- MGDrivE::cubeHomingDrive(cM = cutting, chM = pHDR_M, crM = pRES,
                                 cF = cutting, chF = pHDR_F, crF = pRES,
                                 dF = fMC, dhF = 0, drF = prMR,
                                 s = c('HH'=1-sHH, 'HB'=1-sHB, 'BB'=1-sBB) )
