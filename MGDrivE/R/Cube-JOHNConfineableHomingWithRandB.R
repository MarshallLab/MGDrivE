#################################################################################
## HOMING CONFINEMENT SWITCH MODEL IN R:                                       ##
## John Marshall, 3/14/2018, john.marshall@berkeley.edu                        ##
## Population genetics model for homing drive confinement strategies:          ##
#################################################################################

#' Inheritance Cube: Confinable Homing Drive, John
#'
#' write me
#'
#' @param eM Male homing rate
#' @param eF Female homing rate
#' @param prRF Female no-cost resistant allele generation rate
#' @param prRM Male no-cost resistant allele generation rate
#' @param r Crossover probability
#' @param eta genotype-specific mating fitness
#' @param phi genotype-specific sex ratio at emergence
#' @param omega genotype-specific multiplicative modifier of adult mortality
#' @param xiF genotype-specific female pupatory success
#' @param xiM genotype-specific male pupatory success
#' @param s genotype-specific fractional reduction(increase) in fertility
#'
#' @export
cubeConfinableHomingJOHN <- function(eM=1, eF=1, prRF=0, prRM=0, r=0,eta = NULL, phi = NULL,
                             omega = NULL, xiF = NULL, xiM = NULL, s = NULL){


  #Jared Notes
  # eM <- homing in males?
  # eF <- homing in females?
  # prRF <- R generation rate in females?
  # prRM <- R generatino rate in males?
  # r <- crossover rate?


  ## Resistant allele derivations:
  rhoRF <- (1-eF)*prRF
  rhoBF <- (1-eF)*(1-prRF)
  rhoRM <- (1-eM)*prRM
  rhoBM <- (1-eM)*(1-prRM)

  ## Enumerating the crosses:
  numCrosses <- 36*36
  cross <- vector(mode="list", numCrosses)

  ## Store parents for each cross in the list:
  cross[[ 1]]$father <- "HA_HA"; cross[[ 1]]$mother <- "HA_HA"
  cross[[ 2]]$father <- "HA_HA"; cross[[ 2]]$mother <- "HA_Ha"
  cross[[ 3]]$father <- "HA_HA"; cross[[ 3]]$mother <- "HA_hA"
  cross[[ 4]]$father <- "HA_HA"; cross[[ 4]]$mother <- "HA_ha"
  cross[[ 5]]$father <- "HA_HA"; cross[[ 5]]$mother <- "HA_RA"
  cross[[ 6]]$father <- "HA_HA"; cross[[ 6]]$mother <- "HA_Ra"
  cross[[ 7]]$father <- "HA_HA"; cross[[ 7]]$mother <- "HA_BA"
  cross[[ 8]]$father <- "HA_HA"; cross[[ 8]]$mother <- "HA_Ba"
  cross[[ 9]]$father <- "HA_HA"; cross[[ 9]]$mother <- "Ha_Ha"
  cross[[10]]$father <- "HA_HA"; cross[[10]]$mother <- "Ha_hA"
  cross[[11]]$father <- "HA_HA"; cross[[11]]$mother <- "Ha_ha"
  cross[[12]]$father <- "HA_HA"; cross[[12]]$mother <- "Ha_RA"
  cross[[13]]$father <- "HA_HA"; cross[[13]]$mother <- "Ha_Ra"
  cross[[14]]$father <- "HA_HA"; cross[[14]]$mother <- "Ha_BA"
  cross[[15]]$father <- "HA_HA"; cross[[15]]$mother <- "Ha_Ba"
  cross[[16]]$father <- "HA_HA"; cross[[16]]$mother <- "hA_hA"
  cross[[17]]$father <- "HA_HA"; cross[[17]]$mother <- "hA_ha"
  cross[[18]]$father <- "HA_HA"; cross[[18]]$mother <- "hA_RA"
  cross[[19]]$father <- "HA_HA"; cross[[19]]$mother <- "hA_Ra"
  cross[[20]]$father <- "HA_HA"; cross[[20]]$mother <- "hA_BA"
  cross[[21]]$father <- "HA_HA"; cross[[21]]$mother <- "hA_Ba"
  cross[[22]]$father <- "HA_HA"; cross[[22]]$mother <- "ha_ha"
  cross[[23]]$father <- "HA_HA"; cross[[23]]$mother <- "ha_RA"
  cross[[24]]$father <- "HA_HA"; cross[[24]]$mother <- "ha_Ra"
  cross[[25]]$father <- "HA_HA"; cross[[25]]$mother <- "ha_BA"
  cross[[26]]$father <- "HA_HA"; cross[[26]]$mother <- "ha_Ba"
  cross[[27]]$father <- "HA_HA"; cross[[27]]$mother <- "RA_RA"
  cross[[28]]$father <- "HA_HA"; cross[[28]]$mother <- "RA_Ra"
  cross[[29]]$father <- "HA_HA"; cross[[29]]$mother <- "RA_BA"
  cross[[30]]$father <- "HA_HA"; cross[[30]]$mother <- "RA_Ba"
  cross[[31]]$father <- "HA_HA"; cross[[31]]$mother <- "Ra_Ra"
  cross[[32]]$father <- "HA_HA"; cross[[32]]$mother <- "Ra_BA"
  cross[[33]]$father <- "HA_HA"; cross[[33]]$mother <- "Ra_Ba"
  cross[[34]]$father <- "HA_HA"; cross[[34]]$mother <- "BA_BA"
  cross[[35]]$father <- "HA_HA"; cross[[35]]$mother <- "BA_Ba"
  cross[[36]]$father <- "HA_HA"; cross[[36]]$mother <- "Ba_Ba"

  cross[[37]]$father <- "HA_Ha"; cross[[37]]$mother <- "HA_HA"
  cross[[38]]$father <- "HA_Ha"; cross[[38]]$mother <- "HA_Ha"
  cross[[39]]$father <- "HA_Ha"; cross[[39]]$mother <- "HA_hA"
  cross[[40]]$father <- "HA_Ha"; cross[[40]]$mother <- "HA_ha"
  cross[[41]]$father <- "HA_Ha"; cross[[41]]$mother <- "HA_RA"
  cross[[42]]$father <- "HA_Ha"; cross[[42]]$mother <- "HA_Ra"
  cross[[43]]$father <- "HA_Ha"; cross[[43]]$mother <- "HA_BA"
  cross[[44]]$father <- "HA_Ha"; cross[[44]]$mother <- "HA_Ba"
  cross[[45]]$father <- "HA_Ha"; cross[[45]]$mother <- "Ha_Ha"
  cross[[46]]$father <- "HA_Ha"; cross[[46]]$mother <- "Ha_hA"
  cross[[47]]$father <- "HA_Ha"; cross[[47]]$mother <- "Ha_ha"
  cross[[48]]$father <- "HA_Ha"; cross[[48]]$mother <- "Ha_RA"
  cross[[49]]$father <- "HA_Ha"; cross[[49]]$mother <- "Ha_Ra"
  cross[[50]]$father <- "HA_Ha"; cross[[50]]$mother <- "Ha_BA"
  cross[[51]]$father <- "HA_Ha"; cross[[51]]$mother <- "Ha_Ba"
  cross[[52]]$father <- "HA_Ha"; cross[[52]]$mother <- "hA_hA"
  cross[[53]]$father <- "HA_Ha"; cross[[53]]$mother <- "hA_ha"
  cross[[54]]$father <- "HA_Ha"; cross[[54]]$mother <- "hA_RA"
  cross[[55]]$father <- "HA_Ha"; cross[[55]]$mother <- "hA_Ra"
  cross[[56]]$father <- "HA_Ha"; cross[[56]]$mother <- "hA_BA"
  cross[[57]]$father <- "HA_Ha"; cross[[57]]$mother <- "hA_Ba"
  cross[[58]]$father <- "HA_Ha"; cross[[58]]$mother <- "ha_ha"
  cross[[59]]$father <- "HA_Ha"; cross[[59]]$mother <- "ha_RA"
  cross[[60]]$father <- "HA_Ha"; cross[[60]]$mother <- "ha_Ra"
  cross[[61]]$father <- "HA_Ha"; cross[[61]]$mother <- "ha_BA"
  cross[[62]]$father <- "HA_Ha"; cross[[62]]$mother <- "ha_Ba"
  cross[[63]]$father <- "HA_Ha"; cross[[63]]$mother <- "RA_RA"
  cross[[64]]$father <- "HA_Ha"; cross[[64]]$mother <- "RA_Ra"
  cross[[65]]$father <- "HA_Ha"; cross[[65]]$mother <- "RA_BA"
  cross[[66]]$father <- "HA_Ha"; cross[[66]]$mother <- "RA_Ba"
  cross[[67]]$father <- "HA_Ha"; cross[[67]]$mother <- "Ra_Ra"
  cross[[68]]$father <- "HA_Ha"; cross[[68]]$mother <- "Ra_BA"
  cross[[69]]$father <- "HA_Ha"; cross[[69]]$mother <- "Ra_Ba"
  cross[[70]]$father <- "HA_Ha"; cross[[70]]$mother <- "BA_BA"
  cross[[71]]$father <- "HA_Ha"; cross[[71]]$mother <- "BA_Ba"
  cross[[72]]$father <- "HA_Ha"; cross[[72]]$mother <- "Ba_Ba"

  cross[[73]]$father <- "HA_hA"; cross[[73]]$mother <- "HA_HA"
  cross[[74]]$father <- "HA_hA"; cross[[74]]$mother <- "HA_Ha"
  cross[[75]]$father <- "HA_hA"; cross[[75]]$mother <- "HA_hA"
  cross[[76]]$father <- "HA_hA"; cross[[76]]$mother <- "HA_ha"
  cross[[77]]$father <- "HA_hA"; cross[[77]]$mother <- "HA_RA"
  cross[[78]]$father <- "HA_hA"; cross[[78]]$mother <- "HA_Ra"
  cross[[79]]$father <- "HA_hA"; cross[[79]]$mother <- "HA_BA"
  cross[[80]]$father <- "HA_hA"; cross[[80]]$mother <- "HA_Ba"
  cross[[81]]$father <- "HA_hA"; cross[[81]]$mother <- "Ha_Ha"
  cross[[82]]$father <- "HA_hA"; cross[[82]]$mother <- "Ha_hA"
  cross[[83]]$father <- "HA_hA"; cross[[83]]$mother <- "Ha_ha"
  cross[[84]]$father <- "HA_hA"; cross[[84]]$mother <- "Ha_RA"
  cross[[85]]$father <- "HA_hA"; cross[[85]]$mother <- "Ha_Ra"
  cross[[86]]$father <- "HA_hA"; cross[[86]]$mother <- "Ha_BA"
  cross[[87]]$father <- "HA_hA"; cross[[87]]$mother <- "Ha_Ba"
  cross[[88]]$father <- "HA_hA"; cross[[88]]$mother <- "hA_hA"
  cross[[89]]$father <- "HA_hA"; cross[[89]]$mother <- "hA_ha"
  cross[[90]]$father <- "HA_hA"; cross[[90]]$mother <- "hA_RA"
  cross[[91]]$father <- "HA_hA"; cross[[91]]$mother <- "hA_Ra"
  cross[[92]]$father <- "HA_hA"; cross[[92]]$mother <- "hA_BA"
  cross[[93]]$father <- "HA_hA"; cross[[93]]$mother <- "hA_Ba"
  cross[[94]]$father <- "HA_hA"; cross[[94]]$mother <- "ha_ha"
  cross[[95]]$father <- "HA_hA"; cross[[95]]$mother <- "ha_RA"
  cross[[96]]$father <- "HA_hA"; cross[[96]]$mother <- "ha_Ra"
  cross[[97]]$father <- "HA_hA"; cross[[97]]$mother <- "ha_BA"
  cross[[98]]$father <- "HA_hA"; cross[[98]]$mother <- "ha_Ba"
  cross[[99]]$father <- "HA_hA"; cross[[99]]$mother <- "RA_RA"
  cross[[100]]$father <- "HA_hA"; cross[[100]]$mother <- "RA_Ra"
  cross[[101]]$father <- "HA_hA"; cross[[101]]$mother <- "RA_BA"
  cross[[102]]$father <- "HA_hA"; cross[[102]]$mother <- "RA_Ba"
  cross[[103]]$father <- "HA_hA"; cross[[103]]$mother <- "Ra_Ra"
  cross[[104]]$father <- "HA_hA"; cross[[104]]$mother <- "Ra_BA"
  cross[[105]]$father <- "HA_hA"; cross[[105]]$mother <- "Ra_Ba"
  cross[[106]]$father <- "HA_hA"; cross[[106]]$mother <- "BA_BA"
  cross[[107]]$father <- "HA_hA"; cross[[107]]$mother <- "BA_Ba"
  cross[[108]]$father <- "HA_hA"; cross[[108]]$mother <- "Ba_Ba"

  cross[[109]]$father <- "HA_ha"; cross[[109]]$mother <- "HA_HA"
  cross[[110]]$father <- "HA_ha"; cross[[110]]$mother <- "HA_Ha"
  cross[[111]]$father <- "HA_ha"; cross[[111]]$mother <- "HA_hA"
  cross[[112]]$father <- "HA_ha"; cross[[112]]$mother <- "HA_ha"
  cross[[113]]$father <- "HA_ha"; cross[[113]]$mother <- "HA_RA"
  cross[[114]]$father <- "HA_ha"; cross[[114]]$mother <- "HA_Ra"
  cross[[115]]$father <- "HA_ha"; cross[[115]]$mother <- "HA_BA"
  cross[[116]]$father <- "HA_ha"; cross[[116]]$mother <- "HA_Ba"
  cross[[117]]$father <- "HA_ha"; cross[[117]]$mother <- "Ha_Ha"
  cross[[118]]$father <- "HA_ha"; cross[[118]]$mother <- "Ha_hA"
  cross[[119]]$father <- "HA_ha"; cross[[119]]$mother <- "Ha_ha"
  cross[[120]]$father <- "HA_ha"; cross[[120]]$mother <- "Ha_RA"
  cross[[121]]$father <- "HA_ha"; cross[[121]]$mother <- "Ha_Ra"
  cross[[122]]$father <- "HA_ha"; cross[[122]]$mother <- "Ha_BA"
  cross[[123]]$father <- "HA_ha"; cross[[123]]$mother <- "Ha_Ba"
  cross[[124]]$father <- "HA_ha"; cross[[124]]$mother <- "hA_hA"
  cross[[125]]$father <- "HA_ha"; cross[[125]]$mother <- "hA_ha"
  cross[[126]]$father <- "HA_ha"; cross[[126]]$mother <- "hA_RA"
  cross[[127]]$father <- "HA_ha"; cross[[127]]$mother <- "hA_Ra"
  cross[[128]]$father <- "HA_ha"; cross[[128]]$mother <- "hA_BA"
  cross[[129]]$father <- "HA_ha"; cross[[129]]$mother <- "hA_Ba"
  cross[[130]]$father <- "HA_ha"; cross[[130]]$mother <- "ha_ha"
  cross[[131]]$father <- "HA_ha"; cross[[131]]$mother <- "ha_RA"
  cross[[132]]$father <- "HA_ha"; cross[[132]]$mother <- "ha_Ra"
  cross[[133]]$father <- "HA_ha"; cross[[133]]$mother <- "ha_BA"
  cross[[134]]$father <- "HA_ha"; cross[[134]]$mother <- "ha_Ba"
  cross[[135]]$father <- "HA_ha"; cross[[135]]$mother <- "RA_RA"
  cross[[136]]$father <- "HA_ha"; cross[[136]]$mother <- "RA_Ra"
  cross[[137]]$father <- "HA_ha"; cross[[137]]$mother <- "RA_BA"
  cross[[138]]$father <- "HA_ha"; cross[[138]]$mother <- "RA_Ba"
  cross[[139]]$father <- "HA_ha"; cross[[139]]$mother <- "Ra_Ra"
  cross[[140]]$father <- "HA_ha"; cross[[140]]$mother <- "Ra_BA"
  cross[[141]]$father <- "HA_ha"; cross[[141]]$mother <- "Ra_Ba"
  cross[[142]]$father <- "HA_ha"; cross[[142]]$mother <- "BA_BA"
  cross[[143]]$father <- "HA_ha"; cross[[143]]$mother <- "BA_Ba"
  cross[[144]]$father <- "HA_ha"; cross[[144]]$mother <- "Ba_Ba"

  cross[[145]]$father <- "HA_RA"; cross[[145]]$mother <- "HA_HA"
  cross[[146]]$father <- "HA_RA"; cross[[146]]$mother <- "HA_Ha"
  cross[[147]]$father <- "HA_RA"; cross[[147]]$mother <- "HA_hA"
  cross[[148]]$father <- "HA_RA"; cross[[148]]$mother <- "HA_ha"
  cross[[149]]$father <- "HA_RA"; cross[[149]]$mother <- "HA_RA"
  cross[[150]]$father <- "HA_RA"; cross[[150]]$mother <- "HA_Ra"
  cross[[151]]$father <- "HA_RA"; cross[[151]]$mother <- "HA_BA"
  cross[[152]]$father <- "HA_RA"; cross[[152]]$mother <- "HA_Ba"
  cross[[153]]$father <- "HA_RA"; cross[[153]]$mother <- "Ha_Ha"
  cross[[154]]$father <- "HA_RA"; cross[[154]]$mother <- "Ha_hA"
  cross[[155]]$father <- "HA_RA"; cross[[155]]$mother <- "Ha_ha"
  cross[[156]]$father <- "HA_RA"; cross[[156]]$mother <- "Ha_RA"
  cross[[157]]$father <- "HA_RA"; cross[[157]]$mother <- "Ha_Ra"
  cross[[158]]$father <- "HA_RA"; cross[[158]]$mother <- "Ha_BA"
  cross[[159]]$father <- "HA_RA"; cross[[159]]$mother <- "Ha_Ba"
  cross[[160]]$father <- "HA_RA"; cross[[160]]$mother <- "hA_hA"
  cross[[161]]$father <- "HA_RA"; cross[[161]]$mother <- "hA_ha"
  cross[[162]]$father <- "HA_RA"; cross[[162]]$mother <- "hA_RA"
  cross[[163]]$father <- "HA_RA"; cross[[163]]$mother <- "hA_Ra"
  cross[[164]]$father <- "HA_RA"; cross[[164]]$mother <- "hA_BA"
  cross[[165]]$father <- "HA_RA"; cross[[165]]$mother <- "hA_Ba"
  cross[[166]]$father <- "HA_RA"; cross[[166]]$mother <- "ha_ha"
  cross[[167]]$father <- "HA_RA"; cross[[167]]$mother <- "ha_RA"
  cross[[168]]$father <- "HA_RA"; cross[[168]]$mother <- "ha_Ra"
  cross[[169]]$father <- "HA_RA"; cross[[169]]$mother <- "ha_BA"
  cross[[170]]$father <- "HA_RA"; cross[[170]]$mother <- "ha_Ba"
  cross[[171]]$father <- "HA_RA"; cross[[171]]$mother <- "RA_RA"
  cross[[172]]$father <- "HA_RA"; cross[[172]]$mother <- "RA_Ra"
  cross[[173]]$father <- "HA_RA"; cross[[173]]$mother <- "RA_BA"
  cross[[174]]$father <- "HA_RA"; cross[[174]]$mother <- "RA_Ba"
  cross[[175]]$father <- "HA_RA"; cross[[175]]$mother <- "Ra_Ra"
  cross[[176]]$father <- "HA_RA"; cross[[176]]$mother <- "Ra_BA"
  cross[[177]]$father <- "HA_RA"; cross[[177]]$mother <- "Ra_Ba"
  cross[[178]]$father <- "HA_RA"; cross[[178]]$mother <- "BA_BA"
  cross[[179]]$father <- "HA_RA"; cross[[179]]$mother <- "BA_Ba"
  cross[[180]]$father <- "HA_RA"; cross[[180]]$mother <- "Ba_Ba"

  cross[[181]]$father <- "HA_Ra"; cross[[181]]$mother <- "HA_HA"
  cross[[182]]$father <- "HA_Ra"; cross[[182]]$mother <- "HA_Ha"
  cross[[183]]$father <- "HA_Ra"; cross[[183]]$mother <- "HA_hA"
  cross[[184]]$father <- "HA_Ra"; cross[[184]]$mother <- "HA_ha"
  cross[[185]]$father <- "HA_Ra"; cross[[185]]$mother <- "HA_RA"
  cross[[186]]$father <- "HA_Ra"; cross[[186]]$mother <- "HA_Ra"
  cross[[187]]$father <- "HA_Ra"; cross[[187]]$mother <- "HA_BA"
  cross[[188]]$father <- "HA_Ra"; cross[[188]]$mother <- "HA_Ba"
  cross[[189]]$father <- "HA_Ra"; cross[[189]]$mother <- "Ha_Ha"
  cross[[190]]$father <- "HA_Ra"; cross[[190]]$mother <- "Ha_hA"
  cross[[191]]$father <- "HA_Ra"; cross[[191]]$mother <- "Ha_ha"
  cross[[192]]$father <- "HA_Ra"; cross[[192]]$mother <- "Ha_RA"
  cross[[193]]$father <- "HA_Ra"; cross[[193]]$mother <- "Ha_Ra"
  cross[[194]]$father <- "HA_Ra"; cross[[194]]$mother <- "Ha_BA"
  cross[[195]]$father <- "HA_Ra"; cross[[195]]$mother <- "Ha_Ba"
  cross[[196]]$father <- "HA_Ra"; cross[[196]]$mother <- "hA_hA"
  cross[[197]]$father <- "HA_Ra"; cross[[197]]$mother <- "hA_ha"
  cross[[198]]$father <- "HA_Ra"; cross[[198]]$mother <- "hA_RA"
  cross[[199]]$father <- "HA_Ra"; cross[[199]]$mother <- "hA_Ra"
  cross[[200]]$father <- "HA_Ra"; cross[[200]]$mother <- "hA_BA"
  cross[[201]]$father <- "HA_Ra"; cross[[201]]$mother <- "hA_Ba"
  cross[[202]]$father <- "HA_Ra"; cross[[202]]$mother <- "ha_ha"
  cross[[203]]$father <- "HA_Ra"; cross[[203]]$mother <- "ha_RA"
  cross[[204]]$father <- "HA_Ra"; cross[[204]]$mother <- "ha_Ra"
  cross[[205]]$father <- "HA_Ra"; cross[[205]]$mother <- "ha_BA"
  cross[[206]]$father <- "HA_Ra"; cross[[206]]$mother <- "ha_Ba"
  cross[[207]]$father <- "HA_Ra"; cross[[207]]$mother <- "RA_RA"
  cross[[208]]$father <- "HA_Ra"; cross[[208]]$mother <- "RA_Ra"
  cross[[209]]$father <- "HA_Ra"; cross[[209]]$mother <- "RA_BA"
  cross[[210]]$father <- "HA_Ra"; cross[[210]]$mother <- "RA_Ba"
  cross[[211]]$father <- "HA_Ra"; cross[[211]]$mother <- "Ra_Ra"
  cross[[212]]$father <- "HA_Ra"; cross[[212]]$mother <- "Ra_BA"
  cross[[213]]$father <- "HA_Ra"; cross[[213]]$mother <- "Ra_Ba"
  cross[[214]]$father <- "HA_Ra"; cross[[214]]$mother <- "BA_BA"
  cross[[215]]$father <- "HA_Ra"; cross[[215]]$mother <- "BA_Ba"
  cross[[216]]$father <- "HA_Ra"; cross[[216]]$mother <- "Ba_Ba"

  cross[[217]]$father <- "HA_BA"; cross[[217]]$mother <- "HA_HA"
  cross[[218]]$father <- "HA_BA"; cross[[218]]$mother <- "HA_Ha"
  cross[[219]]$father <- "HA_BA"; cross[[219]]$mother <- "HA_hA"
  cross[[220]]$father <- "HA_BA"; cross[[220]]$mother <- "HA_ha"
  cross[[221]]$father <- "HA_BA"; cross[[221]]$mother <- "HA_RA"
  cross[[222]]$father <- "HA_BA"; cross[[222]]$mother <- "HA_Ra"
  cross[[223]]$father <- "HA_BA"; cross[[223]]$mother <- "HA_BA"
  cross[[224]]$father <- "HA_BA"; cross[[224]]$mother <- "HA_Ba"
  cross[[225]]$father <- "HA_BA"; cross[[225]]$mother <- "Ha_Ha"
  cross[[226]]$father <- "HA_BA"; cross[[226]]$mother <- "Ha_hA"
  cross[[227]]$father <- "HA_BA"; cross[[227]]$mother <- "Ha_ha"
  cross[[228]]$father <- "HA_BA"; cross[[228]]$mother <- "Ha_RA"
  cross[[229]]$father <- "HA_BA"; cross[[229]]$mother <- "Ha_Ra"
  cross[[230]]$father <- "HA_BA"; cross[[230]]$mother <- "Ha_BA"
  cross[[231]]$father <- "HA_BA"; cross[[231]]$mother <- "Ha_Ba"
  cross[[232]]$father <- "HA_BA"; cross[[232]]$mother <- "hA_hA"
  cross[[233]]$father <- "HA_BA"; cross[[233]]$mother <- "hA_ha"
  cross[[234]]$father <- "HA_BA"; cross[[234]]$mother <- "hA_RA"
  cross[[235]]$father <- "HA_BA"; cross[[235]]$mother <- "hA_Ra"
  cross[[236]]$father <- "HA_BA"; cross[[236]]$mother <- "hA_BA"
  cross[[237]]$father <- "HA_BA"; cross[[237]]$mother <- "hA_Ba"
  cross[[238]]$father <- "HA_BA"; cross[[238]]$mother <- "ha_ha"
  cross[[239]]$father <- "HA_BA"; cross[[239]]$mother <- "ha_RA"
  cross[[240]]$father <- "HA_BA"; cross[[240]]$mother <- "ha_Ra"
  cross[[241]]$father <- "HA_BA"; cross[[241]]$mother <- "ha_BA"
  cross[[242]]$father <- "HA_BA"; cross[[242]]$mother <- "ha_Ba"
  cross[[243]]$father <- "HA_BA"; cross[[243]]$mother <- "RA_RA"
  cross[[244]]$father <- "HA_BA"; cross[[244]]$mother <- "RA_Ra"
  cross[[245]]$father <- "HA_BA"; cross[[245]]$mother <- "RA_BA"
  cross[[246]]$father <- "HA_BA"; cross[[246]]$mother <- "RA_Ba"
  cross[[247]]$father <- "HA_BA"; cross[[247]]$mother <- "Ra_Ra"
  cross[[248]]$father <- "HA_BA"; cross[[248]]$mother <- "Ra_BA"
  cross[[249]]$father <- "HA_BA"; cross[[249]]$mother <- "Ra_Ba"
  cross[[250]]$father <- "HA_BA"; cross[[250]]$mother <- "BA_BA"
  cross[[251]]$father <- "HA_BA"; cross[[251]]$mother <- "BA_Ba"
  cross[[252]]$father <- "HA_BA"; cross[[252]]$mother <- "Ba_Ba"

  cross[[253]]$father <- "HA_Ba"; cross[[253]]$mother <- "HA_HA"
  cross[[254]]$father <- "HA_Ba"; cross[[254]]$mother <- "HA_Ha"
  cross[[255]]$father <- "HA_Ba"; cross[[255]]$mother <- "HA_hA"
  cross[[256]]$father <- "HA_Ba"; cross[[256]]$mother <- "HA_ha"
  cross[[257]]$father <- "HA_Ba"; cross[[257]]$mother <- "HA_RA"
  cross[[258]]$father <- "HA_Ba"; cross[[258]]$mother <- "HA_Ra"
  cross[[259]]$father <- "HA_Ba"; cross[[259]]$mother <- "HA_BA"
  cross[[260]]$father <- "HA_Ba"; cross[[260]]$mother <- "HA_Ba"
  cross[[261]]$father <- "HA_Ba"; cross[[261]]$mother <- "Ha_Ha"
  cross[[262]]$father <- "HA_Ba"; cross[[262]]$mother <- "Ha_hA"
  cross[[263]]$father <- "HA_Ba"; cross[[263]]$mother <- "Ha_ha"
  cross[[264]]$father <- "HA_Ba"; cross[[264]]$mother <- "Ha_RA"
  cross[[265]]$father <- "HA_Ba"; cross[[265]]$mother <- "Ha_Ra"
  cross[[266]]$father <- "HA_Ba"; cross[[266]]$mother <- "Ha_BA"
  cross[[267]]$father <- "HA_Ba"; cross[[267]]$mother <- "Ha_Ba"
  cross[[268]]$father <- "HA_Ba"; cross[[268]]$mother <- "hA_hA"
  cross[[269]]$father <- "HA_Ba"; cross[[269]]$mother <- "hA_ha"
  cross[[270]]$father <- "HA_Ba"; cross[[270]]$mother <- "hA_RA"
  cross[[271]]$father <- "HA_Ba"; cross[[271]]$mother <- "hA_Ra"
  cross[[272]]$father <- "HA_Ba"; cross[[272]]$mother <- "hA_BA"
  cross[[273]]$father <- "HA_Ba"; cross[[273]]$mother <- "hA_Ba"
  cross[[274]]$father <- "HA_Ba"; cross[[274]]$mother <- "ha_ha"
  cross[[275]]$father <- "HA_Ba"; cross[[275]]$mother <- "ha_RA"
  cross[[276]]$father <- "HA_Ba"; cross[[276]]$mother <- "ha_Ra"
  cross[[277]]$father <- "HA_Ba"; cross[[277]]$mother <- "ha_BA"
  cross[[278]]$father <- "HA_Ba"; cross[[278]]$mother <- "ha_Ba"
  cross[[279]]$father <- "HA_Ba"; cross[[279]]$mother <- "RA_RA"
  cross[[280]]$father <- "HA_Ba"; cross[[280]]$mother <- "RA_Ra"
  cross[[281]]$father <- "HA_Ba"; cross[[281]]$mother <- "RA_BA"
  cross[[282]]$father <- "HA_Ba"; cross[[282]]$mother <- "RA_Ba"
  cross[[283]]$father <- "HA_Ba"; cross[[283]]$mother <- "Ra_Ra"
  cross[[284]]$father <- "HA_Ba"; cross[[284]]$mother <- "Ra_BA"
  cross[[285]]$father <- "HA_Ba"; cross[[285]]$mother <- "Ra_Ba"
  cross[[286]]$father <- "HA_Ba"; cross[[286]]$mother <- "BA_BA"
  cross[[287]]$father <- "HA_Ba"; cross[[287]]$mother <- "BA_Ba"
  cross[[288]]$father <- "HA_Ba"; cross[[288]]$mother <- "Ba_Ba"

  cross[[289]]$father <- "Ha_Ha"; cross[[289]]$mother <- "HA_HA"
  cross[[290]]$father <- "Ha_Ha"; cross[[290]]$mother <- "HA_Ha"
  cross[[291]]$father <- "Ha_Ha"; cross[[291]]$mother <- "HA_hA"
  cross[[292]]$father <- "Ha_Ha"; cross[[292]]$mother <- "HA_ha"
  cross[[293]]$father <- "Ha_Ha"; cross[[293]]$mother <- "HA_RA"
  cross[[294]]$father <- "Ha_Ha"; cross[[294]]$mother <- "HA_Ra"
  cross[[295]]$father <- "Ha_Ha"; cross[[295]]$mother <- "HA_BA"
  cross[[296]]$father <- "Ha_Ha"; cross[[296]]$mother <- "HA_Ba"
  cross[[297]]$father <- "Ha_Ha"; cross[[297]]$mother <- "Ha_Ha"
  cross[[298]]$father <- "Ha_Ha"; cross[[298]]$mother <- "Ha_hA"
  cross[[299]]$father <- "Ha_Ha"; cross[[299]]$mother <- "Ha_ha"
  cross[[300]]$father <- "Ha_Ha"; cross[[300]]$mother <- "Ha_RA"
  cross[[301]]$father <- "Ha_Ha"; cross[[301]]$mother <- "Ha_Ra"
  cross[[302]]$father <- "Ha_Ha"; cross[[302]]$mother <- "Ha_BA"
  cross[[303]]$father <- "Ha_Ha"; cross[[303]]$mother <- "Ha_Ba"
  cross[[304]]$father <- "Ha_Ha"; cross[[304]]$mother <- "hA_hA"
  cross[[305]]$father <- "Ha_Ha"; cross[[305]]$mother <- "hA_ha"
  cross[[306]]$father <- "Ha_Ha"; cross[[306]]$mother <- "hA_RA"
  cross[[307]]$father <- "Ha_Ha"; cross[[307]]$mother <- "hA_Ra"
  cross[[308]]$father <- "Ha_Ha"; cross[[308]]$mother <- "hA_BA"
  cross[[309]]$father <- "Ha_Ha"; cross[[309]]$mother <- "hA_Ba"
  cross[[310]]$father <- "Ha_Ha"; cross[[310]]$mother <- "ha_ha"
  cross[[311]]$father <- "Ha_Ha"; cross[[311]]$mother <- "ha_RA"
  cross[[312]]$father <- "Ha_Ha"; cross[[312]]$mother <- "ha_Ra"
  cross[[313]]$father <- "Ha_Ha"; cross[[313]]$mother <- "ha_BA"
  cross[[314]]$father <- "Ha_Ha"; cross[[314]]$mother <- "ha_Ba"
  cross[[315]]$father <- "Ha_Ha"; cross[[315]]$mother <- "RA_RA"
  cross[[316]]$father <- "Ha_Ha"; cross[[316]]$mother <- "RA_Ra"
  cross[[317]]$father <- "Ha_Ha"; cross[[317]]$mother <- "RA_BA"
  cross[[318]]$father <- "Ha_Ha"; cross[[318]]$mother <- "RA_Ba"
  cross[[319]]$father <- "Ha_Ha"; cross[[319]]$mother <- "Ra_Ra"
  cross[[320]]$father <- "Ha_Ha"; cross[[320]]$mother <- "Ra_BA"
  cross[[321]]$father <- "Ha_Ha"; cross[[321]]$mother <- "Ra_Ba"
  cross[[322]]$father <- "Ha_Ha"; cross[[322]]$mother <- "BA_BA"
  cross[[323]]$father <- "Ha_Ha"; cross[[323]]$mother <- "BA_Ba"
  cross[[324]]$father <- "Ha_Ha"; cross[[324]]$mother <- "Ba_Ba"

  cross[[325]]$father <- "Ha_hA"; cross[[325]]$mother <- "HA_HA"
  cross[[326]]$father <- "Ha_hA"; cross[[326]]$mother <- "HA_Ha"
  cross[[327]]$father <- "Ha_hA"; cross[[327]]$mother <- "HA_hA"
  cross[[328]]$father <- "Ha_hA"; cross[[328]]$mother <- "HA_ha"
  cross[[329]]$father <- "Ha_hA"; cross[[329]]$mother <- "HA_RA"
  cross[[330]]$father <- "Ha_hA"; cross[[330]]$mother <- "HA_Ra"
  cross[[331]]$father <- "Ha_hA"; cross[[331]]$mother <- "HA_BA"
  cross[[332]]$father <- "Ha_hA"; cross[[332]]$mother <- "HA_Ba"
  cross[[333]]$father <- "Ha_hA"; cross[[333]]$mother <- "Ha_Ha"
  cross[[334]]$father <- "Ha_hA"; cross[[334]]$mother <- "Ha_hA"
  cross[[335]]$father <- "Ha_hA"; cross[[335]]$mother <- "Ha_ha"
  cross[[336]]$father <- "Ha_hA"; cross[[336]]$mother <- "Ha_RA"
  cross[[337]]$father <- "Ha_hA"; cross[[337]]$mother <- "Ha_Ra"
  cross[[338]]$father <- "Ha_hA"; cross[[338]]$mother <- "Ha_BA"
  cross[[339]]$father <- "Ha_hA"; cross[[339]]$mother <- "Ha_Ba"
  cross[[340]]$father <- "Ha_hA"; cross[[340]]$mother <- "hA_hA"
  cross[[341]]$father <- "Ha_hA"; cross[[341]]$mother <- "hA_ha"
  cross[[342]]$father <- "Ha_hA"; cross[[342]]$mother <- "hA_RA"
  cross[[343]]$father <- "Ha_hA"; cross[[343]]$mother <- "hA_Ra"
  cross[[344]]$father <- "Ha_hA"; cross[[344]]$mother <- "hA_BA"
  cross[[345]]$father <- "Ha_hA"; cross[[345]]$mother <- "hA_Ba"
  cross[[346]]$father <- "Ha_hA"; cross[[346]]$mother <- "ha_ha"
  cross[[347]]$father <- "Ha_hA"; cross[[347]]$mother <- "ha_RA"
  cross[[348]]$father <- "Ha_hA"; cross[[348]]$mother <- "ha_Ra"
  cross[[349]]$father <- "Ha_hA"; cross[[349]]$mother <- "ha_BA"
  cross[[350]]$father <- "Ha_hA"; cross[[350]]$mother <- "ha_Ba"
  cross[[351]]$father <- "Ha_hA"; cross[[351]]$mother <- "RA_RA"
  cross[[352]]$father <- "Ha_hA"; cross[[352]]$mother <- "RA_Ra"
  cross[[353]]$father <- "Ha_hA"; cross[[353]]$mother <- "RA_BA"
  cross[[354]]$father <- "Ha_hA"; cross[[354]]$mother <- "RA_Ba"
  cross[[355]]$father <- "Ha_hA"; cross[[355]]$mother <- "Ra_Ra"
  cross[[356]]$father <- "Ha_hA"; cross[[356]]$mother <- "Ra_BA"
  cross[[357]]$father <- "Ha_hA"; cross[[357]]$mother <- "Ra_Ba"
  cross[[358]]$father <- "Ha_hA"; cross[[358]]$mother <- "BA_BA"
  cross[[359]]$father <- "Ha_hA"; cross[[359]]$mother <- "BA_Ba"
  cross[[360]]$father <- "Ha_hA"; cross[[360]]$mother <- "Ba_Ba"

  cross[[361]]$father <- "Ha_ha"; cross[[361]]$mother <- "HA_HA"
  cross[[362]]$father <- "Ha_ha"; cross[[362]]$mother <- "HA_Ha"
  cross[[363]]$father <- "Ha_ha"; cross[[363]]$mother <- "HA_hA"
  cross[[364]]$father <- "Ha_ha"; cross[[364]]$mother <- "HA_ha"
  cross[[365]]$father <- "Ha_ha"; cross[[365]]$mother <- "HA_RA"
  cross[[366]]$father <- "Ha_ha"; cross[[366]]$mother <- "HA_Ra"
  cross[[367]]$father <- "Ha_ha"; cross[[367]]$mother <- "HA_BA"
  cross[[368]]$father <- "Ha_ha"; cross[[368]]$mother <- "HA_Ba"
  cross[[369]]$father <- "Ha_ha"; cross[[369]]$mother <- "Ha_Ha"
  cross[[370]]$father <- "Ha_ha"; cross[[370]]$mother <- "Ha_hA"
  cross[[371]]$father <- "Ha_ha"; cross[[371]]$mother <- "Ha_ha"
  cross[[372]]$father <- "Ha_ha"; cross[[372]]$mother <- "Ha_RA"
  cross[[373]]$father <- "Ha_ha"; cross[[373]]$mother <- "Ha_Ra"
  cross[[374]]$father <- "Ha_ha"; cross[[374]]$mother <- "Ha_BA"
  cross[[375]]$father <- "Ha_ha"; cross[[375]]$mother <- "Ha_Ba"
  cross[[376]]$father <- "Ha_ha"; cross[[376]]$mother <- "hA_hA"
  cross[[377]]$father <- "Ha_ha"; cross[[377]]$mother <- "hA_ha"
  cross[[378]]$father <- "Ha_ha"; cross[[378]]$mother <- "hA_RA"
  cross[[379]]$father <- "Ha_ha"; cross[[379]]$mother <- "hA_Ra"
  cross[[380]]$father <- "Ha_ha"; cross[[380]]$mother <- "hA_BA"
  cross[[381]]$father <- "Ha_ha"; cross[[381]]$mother <- "hA_Ba"
  cross[[382]]$father <- "Ha_ha"; cross[[382]]$mother <- "ha_ha"
  cross[[383]]$father <- "Ha_ha"; cross[[383]]$mother <- "ha_RA"
  cross[[384]]$father <- "Ha_ha"; cross[[384]]$mother <- "ha_Ra"
  cross[[385]]$father <- "Ha_ha"; cross[[385]]$mother <- "ha_BA"
  cross[[386]]$father <- "Ha_ha"; cross[[386]]$mother <- "ha_Ba"
  cross[[387]]$father <- "Ha_ha"; cross[[387]]$mother <- "RA_RA"
  cross[[388]]$father <- "Ha_ha"; cross[[388]]$mother <- "RA_Ra"
  cross[[389]]$father <- "Ha_ha"; cross[[389]]$mother <- "RA_BA"
  cross[[390]]$father <- "Ha_ha"; cross[[390]]$mother <- "RA_Ba"
  cross[[391]]$father <- "Ha_ha"; cross[[391]]$mother <- "Ra_Ra"
  cross[[392]]$father <- "Ha_ha"; cross[[392]]$mother <- "Ra_BA"
  cross[[393]]$father <- "Ha_ha"; cross[[393]]$mother <- "Ra_Ba"
  cross[[394]]$father <- "Ha_ha"; cross[[394]]$mother <- "BA_BA"
  cross[[395]]$father <- "Ha_ha"; cross[[395]]$mother <- "BA_Ba"
  cross[[396]]$father <- "Ha_ha"; cross[[396]]$mother <- "Ba_Ba"

  cross[[397]]$father <- "Ha_RA"; cross[[397]]$mother <- "HA_HA"
  cross[[398]]$father <- "Ha_RA"; cross[[398]]$mother <- "HA_Ha"
  cross[[399]]$father <- "Ha_RA"; cross[[399]]$mother <- "HA_hA"
  cross[[400]]$father <- "Ha_RA"; cross[[400]]$mother <- "HA_ha"
  cross[[401]]$father <- "Ha_RA"; cross[[401]]$mother <- "HA_RA"
  cross[[402]]$father <- "Ha_RA"; cross[[402]]$mother <- "HA_Ra"
  cross[[403]]$father <- "Ha_RA"; cross[[403]]$mother <- "HA_BA"
  cross[[404]]$father <- "Ha_RA"; cross[[404]]$mother <- "HA_Ba"
  cross[[405]]$father <- "Ha_RA"; cross[[405]]$mother <- "Ha_Ha"
  cross[[406]]$father <- "Ha_RA"; cross[[406]]$mother <- "Ha_hA"
  cross[[407]]$father <- "Ha_RA"; cross[[407]]$mother <- "Ha_ha"
  cross[[408]]$father <- "Ha_RA"; cross[[408]]$mother <- "Ha_RA"
  cross[[409]]$father <- "Ha_RA"; cross[[409]]$mother <- "Ha_Ra"
  cross[[410]]$father <- "Ha_RA"; cross[[410]]$mother <- "Ha_BA"
  cross[[411]]$father <- "Ha_RA"; cross[[411]]$mother <- "Ha_Ba"
  cross[[412]]$father <- "Ha_RA"; cross[[412]]$mother <- "hA_hA"
  cross[[413]]$father <- "Ha_RA"; cross[[413]]$mother <- "hA_ha"
  cross[[414]]$father <- "Ha_RA"; cross[[414]]$mother <- "hA_RA"
  cross[[415]]$father <- "Ha_RA"; cross[[415]]$mother <- "hA_Ra"
  cross[[416]]$father <- "Ha_RA"; cross[[416]]$mother <- "hA_BA"
  cross[[417]]$father <- "Ha_RA"; cross[[417]]$mother <- "hA_Ba"
  cross[[418]]$father <- "Ha_RA"; cross[[418]]$mother <- "ha_ha"
  cross[[419]]$father <- "Ha_RA"; cross[[419]]$mother <- "ha_RA"
  cross[[420]]$father <- "Ha_RA"; cross[[420]]$mother <- "ha_Ra"
  cross[[421]]$father <- "Ha_RA"; cross[[421]]$mother <- "ha_BA"
  cross[[422]]$father <- "Ha_RA"; cross[[422]]$mother <- "ha_Ba"
  cross[[423]]$father <- "Ha_RA"; cross[[423]]$mother <- "RA_RA"
  cross[[424]]$father <- "Ha_RA"; cross[[424]]$mother <- "RA_Ra"
  cross[[425]]$father <- "Ha_RA"; cross[[425]]$mother <- "RA_BA"
  cross[[426]]$father <- "Ha_RA"; cross[[426]]$mother <- "RA_Ba"
  cross[[427]]$father <- "Ha_RA"; cross[[427]]$mother <- "Ra_Ra"
  cross[[428]]$father <- "Ha_RA"; cross[[428]]$mother <- "Ra_BA"
  cross[[429]]$father <- "Ha_RA"; cross[[429]]$mother <- "Ra_Ba"
  cross[[430]]$father <- "Ha_RA"; cross[[430]]$mother <- "BA_BA"
  cross[[431]]$father <- "Ha_RA"; cross[[431]]$mother <- "BA_Ba"
  cross[[432]]$father <- "Ha_RA"; cross[[432]]$mother <- "Ba_Ba"

  cross[[433]]$father <- "Ha_Ra"; cross[[433]]$mother <- "HA_HA"
  cross[[434]]$father <- "Ha_Ra"; cross[[434]]$mother <- "HA_Ha"
  cross[[435]]$father <- "Ha_Ra"; cross[[435]]$mother <- "HA_hA"
  cross[[436]]$father <- "Ha_Ra"; cross[[436]]$mother <- "HA_ha"
  cross[[437]]$father <- "Ha_Ra"; cross[[437]]$mother <- "HA_RA"
  cross[[438]]$father <- "Ha_Ra"; cross[[438]]$mother <- "HA_Ra"
  cross[[439]]$father <- "Ha_Ra"; cross[[439]]$mother <- "HA_BA"
  cross[[440]]$father <- "Ha_Ra"; cross[[440]]$mother <- "HA_Ba"
  cross[[441]]$father <- "Ha_Ra"; cross[[441]]$mother <- "Ha_Ha"
  cross[[442]]$father <- "Ha_Ra"; cross[[442]]$mother <- "Ha_hA"
  cross[[443]]$father <- "Ha_Ra"; cross[[443]]$mother <- "Ha_ha"
  cross[[444]]$father <- "Ha_Ra"; cross[[444]]$mother <- "Ha_RA"
  cross[[445]]$father <- "Ha_Ra"; cross[[445]]$mother <- "Ha_Ra"
  cross[[446]]$father <- "Ha_Ra"; cross[[446]]$mother <- "Ha_BA"
  cross[[447]]$father <- "Ha_Ra"; cross[[447]]$mother <- "Ha_Ba"
  cross[[448]]$father <- "Ha_Ra"; cross[[448]]$mother <- "hA_hA"
  cross[[449]]$father <- "Ha_Ra"; cross[[449]]$mother <- "hA_ha"
  cross[[450]]$father <- "Ha_Ra"; cross[[450]]$mother <- "hA_RA"
  cross[[451]]$father <- "Ha_Ra"; cross[[451]]$mother <- "hA_Ra"
  cross[[452]]$father <- "Ha_Ra"; cross[[452]]$mother <- "hA_BA"
  cross[[453]]$father <- "Ha_Ra"; cross[[453]]$mother <- "hA_Ba"
  cross[[454]]$father <- "Ha_Ra"; cross[[454]]$mother <- "ha_ha"
  cross[[455]]$father <- "Ha_Ra"; cross[[455]]$mother <- "ha_RA"
  cross[[456]]$father <- "Ha_Ra"; cross[[456]]$mother <- "ha_Ra"
  cross[[457]]$father <- "Ha_Ra"; cross[[457]]$mother <- "ha_BA"
  cross[[458]]$father <- "Ha_Ra"; cross[[458]]$mother <- "ha_Ba"
  cross[[459]]$father <- "Ha_Ra"; cross[[459]]$mother <- "RA_RA"
  cross[[460]]$father <- "Ha_Ra"; cross[[460]]$mother <- "RA_Ra"
  cross[[461]]$father <- "Ha_Ra"; cross[[461]]$mother <- "RA_BA"
  cross[[462]]$father <- "Ha_Ra"; cross[[462]]$mother <- "RA_Ba"
  cross[[463]]$father <- "Ha_Ra"; cross[[463]]$mother <- "Ra_Ra"
  cross[[464]]$father <- "Ha_Ra"; cross[[464]]$mother <- "Ra_BA"
  cross[[465]]$father <- "Ha_Ra"; cross[[465]]$mother <- "Ra_Ba"
  cross[[466]]$father <- "Ha_Ra"; cross[[466]]$mother <- "BA_BA"
  cross[[467]]$father <- "Ha_Ra"; cross[[467]]$mother <- "BA_Ba"
  cross[[468]]$father <- "Ha_Ra"; cross[[468]]$mother <- "Ba_Ba"

  cross[[469]]$father <- "Ha_BA"; cross[[469]]$mother <- "HA_HA"
  cross[[470]]$father <- "Ha_BA"; cross[[470]]$mother <- "HA_Ha"
  cross[[471]]$father <- "Ha_BA"; cross[[471]]$mother <- "HA_hA"
  cross[[472]]$father <- "Ha_BA"; cross[[472]]$mother <- "HA_ha"
  cross[[473]]$father <- "Ha_BA"; cross[[473]]$mother <- "HA_RA"
  cross[[474]]$father <- "Ha_BA"; cross[[474]]$mother <- "HA_Ra"
  cross[[475]]$father <- "Ha_BA"; cross[[475]]$mother <- "HA_BA"
  cross[[476]]$father <- "Ha_BA"; cross[[476]]$mother <- "HA_Ba"
  cross[[477]]$father <- "Ha_BA"; cross[[477]]$mother <- "Ha_Ha"
  cross[[478]]$father <- "Ha_BA"; cross[[478]]$mother <- "Ha_hA"
  cross[[479]]$father <- "Ha_BA"; cross[[479]]$mother <- "Ha_ha"
  cross[[480]]$father <- "Ha_BA"; cross[[480]]$mother <- "Ha_RA"
  cross[[481]]$father <- "Ha_BA"; cross[[481]]$mother <- "Ha_Ra"
  cross[[482]]$father <- "Ha_BA"; cross[[482]]$mother <- "Ha_BA"
  cross[[483]]$father <- "Ha_BA"; cross[[483]]$mother <- "Ha_Ba"
  cross[[484]]$father <- "Ha_BA"; cross[[484]]$mother <- "hA_hA"
  cross[[485]]$father <- "Ha_BA"; cross[[485]]$mother <- "hA_ha"
  cross[[486]]$father <- "Ha_BA"; cross[[486]]$mother <- "hA_RA"
  cross[[487]]$father <- "Ha_BA"; cross[[487]]$mother <- "hA_Ra"
  cross[[488]]$father <- "Ha_BA"; cross[[488]]$mother <- "hA_BA"
  cross[[489]]$father <- "Ha_BA"; cross[[489]]$mother <- "hA_Ba"
  cross[[490]]$father <- "Ha_BA"; cross[[490]]$mother <- "ha_ha"
  cross[[491]]$father <- "Ha_BA"; cross[[491]]$mother <- "ha_RA"
  cross[[492]]$father <- "Ha_BA"; cross[[492]]$mother <- "ha_Ra"
  cross[[493]]$father <- "Ha_BA"; cross[[493]]$mother <- "ha_BA"
  cross[[494]]$father <- "Ha_BA"; cross[[494]]$mother <- "ha_Ba"
  cross[[495]]$father <- "Ha_BA"; cross[[495]]$mother <- "RA_RA"
  cross[[496]]$father <- "Ha_BA"; cross[[496]]$mother <- "RA_Ra"
  cross[[497]]$father <- "Ha_BA"; cross[[497]]$mother <- "RA_BA"
  cross[[498]]$father <- "Ha_BA"; cross[[498]]$mother <- "RA_Ba"
  cross[[499]]$father <- "Ha_BA"; cross[[499]]$mother <- "Ra_Ra"
  cross[[500]]$father <- "Ha_BA"; cross[[500]]$mother <- "Ra_BA"
  cross[[501]]$father <- "Ha_BA"; cross[[501]]$mother <- "Ra_Ba"
  cross[[502]]$father <- "Ha_BA"; cross[[502]]$mother <- "BA_BA"
  cross[[503]]$father <- "Ha_BA"; cross[[503]]$mother <- "BA_Ba"
  cross[[504]]$father <- "Ha_BA"; cross[[504]]$mother <- "Ba_Ba"

  cross[[505]]$father <- "Ha_Ba"; cross[[505]]$mother <- "HA_HA"
  cross[[506]]$father <- "Ha_Ba"; cross[[506]]$mother <- "HA_Ha"
  cross[[507]]$father <- "Ha_Ba"; cross[[507]]$mother <- "HA_hA"
  cross[[508]]$father <- "Ha_Ba"; cross[[508]]$mother <- "HA_ha"
  cross[[509]]$father <- "Ha_Ba"; cross[[509]]$mother <- "HA_RA"
  cross[[510]]$father <- "Ha_Ba"; cross[[510]]$mother <- "HA_Ra"
  cross[[511]]$father <- "Ha_Ba"; cross[[511]]$mother <- "HA_BA"
  cross[[512]]$father <- "Ha_Ba"; cross[[512]]$mother <- "HA_Ba"
  cross[[513]]$father <- "Ha_Ba"; cross[[513]]$mother <- "Ha_Ha"
  cross[[514]]$father <- "Ha_Ba"; cross[[514]]$mother <- "Ha_hA"
  cross[[515]]$father <- "Ha_Ba"; cross[[515]]$mother <- "Ha_ha"
  cross[[516]]$father <- "Ha_Ba"; cross[[516]]$mother <- "Ha_RA"
  cross[[517]]$father <- "Ha_Ba"; cross[[517]]$mother <- "Ha_Ra"
  cross[[518]]$father <- "Ha_Ba"; cross[[518]]$mother <- "Ha_BA"
  cross[[519]]$father <- "Ha_Ba"; cross[[519]]$mother <- "Ha_Ba"
  cross[[520]]$father <- "Ha_Ba"; cross[[520]]$mother <- "hA_hA"
  cross[[521]]$father <- "Ha_Ba"; cross[[521]]$mother <- "hA_ha"
  cross[[522]]$father <- "Ha_Ba"; cross[[522]]$mother <- "hA_RA"
  cross[[523]]$father <- "Ha_Ba"; cross[[523]]$mother <- "hA_Ra"
  cross[[524]]$father <- "Ha_Ba"; cross[[524]]$mother <- "hA_BA"
  cross[[525]]$father <- "Ha_Ba"; cross[[525]]$mother <- "hA_Ba"
  cross[[526]]$father <- "Ha_Ba"; cross[[526]]$mother <- "ha_ha"
  cross[[527]]$father <- "Ha_Ba"; cross[[527]]$mother <- "ha_RA"
  cross[[528]]$father <- "Ha_Ba"; cross[[528]]$mother <- "ha_Ra"
  cross[[529]]$father <- "Ha_Ba"; cross[[529]]$mother <- "ha_BA"
  cross[[530]]$father <- "Ha_Ba"; cross[[530]]$mother <- "ha_Ba"
  cross[[531]]$father <- "Ha_Ba"; cross[[531]]$mother <- "RA_RA"
  cross[[532]]$father <- "Ha_Ba"; cross[[532]]$mother <- "RA_Ra"
  cross[[533]]$father <- "Ha_Ba"; cross[[533]]$mother <- "RA_BA"
  cross[[534]]$father <- "Ha_Ba"; cross[[534]]$mother <- "RA_Ba"
  cross[[535]]$father <- "Ha_Ba"; cross[[535]]$mother <- "Ra_Ra"
  cross[[536]]$father <- "Ha_Ba"; cross[[536]]$mother <- "Ra_BA"
  cross[[537]]$father <- "Ha_Ba"; cross[[537]]$mother <- "Ra_Ba"
  cross[[538]]$father <- "Ha_Ba"; cross[[538]]$mother <- "BA_BA"
  cross[[539]]$father <- "Ha_Ba"; cross[[539]]$mother <- "BA_Ba"
  cross[[540]]$father <- "Ha_Ba"; cross[[540]]$mother <- "Ba_Ba"

  cross[[541]]$father <- "hA_hA"; cross[[541]]$mother <- "HA_HA"
  cross[[542]]$father <- "hA_hA"; cross[[542]]$mother <- "HA_Ha"
  cross[[543]]$father <- "hA_hA"; cross[[543]]$mother <- "HA_hA"
  cross[[544]]$father <- "hA_hA"; cross[[544]]$mother <- "HA_ha"
  cross[[545]]$father <- "hA_hA"; cross[[545]]$mother <- "HA_RA"
  cross[[546]]$father <- "hA_hA"; cross[[546]]$mother <- "HA_Ra"
  cross[[547]]$father <- "hA_hA"; cross[[547]]$mother <- "HA_BA"
  cross[[548]]$father <- "hA_hA"; cross[[548]]$mother <- "HA_Ba"
  cross[[549]]$father <- "hA_hA"; cross[[549]]$mother <- "Ha_Ha"
  cross[[550]]$father <- "hA_hA"; cross[[550]]$mother <- "Ha_hA"
  cross[[551]]$father <- "hA_hA"; cross[[551]]$mother <- "Ha_ha"
  cross[[552]]$father <- "hA_hA"; cross[[552]]$mother <- "Ha_RA"
  cross[[553]]$father <- "hA_hA"; cross[[553]]$mother <- "Ha_Ra"
  cross[[554]]$father <- "hA_hA"; cross[[554]]$mother <- "Ha_BA"
  cross[[555]]$father <- "hA_hA"; cross[[555]]$mother <- "Ha_Ba"
  cross[[556]]$father <- "hA_hA"; cross[[556]]$mother <- "hA_hA"
  cross[[557]]$father <- "hA_hA"; cross[[557]]$mother <- "hA_ha"
  cross[[558]]$father <- "hA_hA"; cross[[558]]$mother <- "hA_RA"
  cross[[559]]$father <- "hA_hA"; cross[[559]]$mother <- "hA_Ra"
  cross[[560]]$father <- "hA_hA"; cross[[560]]$mother <- "hA_BA"
  cross[[561]]$father <- "hA_hA"; cross[[561]]$mother <- "hA_Ba"
  cross[[562]]$father <- "hA_hA"; cross[[562]]$mother <- "ha_ha"
  cross[[563]]$father <- "hA_hA"; cross[[563]]$mother <- "ha_RA"
  cross[[564]]$father <- "hA_hA"; cross[[564]]$mother <- "ha_Ra"
  cross[[565]]$father <- "hA_hA"; cross[[565]]$mother <- "ha_BA"
  cross[[566]]$father <- "hA_hA"; cross[[566]]$mother <- "ha_Ba"
  cross[[567]]$father <- "hA_hA"; cross[[567]]$mother <- "RA_RA"
  cross[[568]]$father <- "hA_hA"; cross[[568]]$mother <- "RA_Ra"
  cross[[569]]$father <- "hA_hA"; cross[[569]]$mother <- "RA_BA"
  cross[[570]]$father <- "hA_hA"; cross[[570]]$mother <- "RA_Ba"
  cross[[571]]$father <- "hA_hA"; cross[[571]]$mother <- "Ra_Ra"
  cross[[572]]$father <- "hA_hA"; cross[[572]]$mother <- "Ra_BA"
  cross[[573]]$father <- "hA_hA"; cross[[573]]$mother <- "Ra_Ba"
  cross[[574]]$father <- "hA_hA"; cross[[574]]$mother <- "BA_BA"
  cross[[575]]$father <- "hA_hA"; cross[[575]]$mother <- "BA_Ba"
  cross[[576]]$father <- "hA_hA"; cross[[576]]$mother <- "Ba_Ba"

  cross[[577]]$father <- "hA_ha"; cross[[577]]$mother <- "HA_HA"
  cross[[578]]$father <- "hA_ha"; cross[[578]]$mother <- "HA_Ha"
  cross[[579]]$father <- "hA_ha"; cross[[579]]$mother <- "HA_hA"
  cross[[580]]$father <- "hA_ha"; cross[[580]]$mother <- "HA_ha"
  cross[[581]]$father <- "hA_ha"; cross[[581]]$mother <- "HA_RA"
  cross[[582]]$father <- "hA_ha"; cross[[582]]$mother <- "HA_Ra"
  cross[[583]]$father <- "hA_ha"; cross[[583]]$mother <- "HA_BA"
  cross[[584]]$father <- "hA_ha"; cross[[584]]$mother <- "HA_Ba"
  cross[[585]]$father <- "hA_ha"; cross[[585]]$mother <- "Ha_Ha"
  cross[[586]]$father <- "hA_ha"; cross[[586]]$mother <- "Ha_hA"
  cross[[587]]$father <- "hA_ha"; cross[[587]]$mother <- "Ha_ha"
  cross[[588]]$father <- "hA_ha"; cross[[588]]$mother <- "Ha_RA"
  cross[[589]]$father <- "hA_ha"; cross[[589]]$mother <- "Ha_Ra"
  cross[[590]]$father <- "hA_ha"; cross[[590]]$mother <- "Ha_BA"
  cross[[591]]$father <- "hA_ha"; cross[[591]]$mother <- "Ha_Ba"
  cross[[592]]$father <- "hA_ha"; cross[[592]]$mother <- "hA_hA"
  cross[[593]]$father <- "hA_ha"; cross[[593]]$mother <- "hA_ha"
  cross[[594]]$father <- "hA_ha"; cross[[594]]$mother <- "hA_RA"
  cross[[595]]$father <- "hA_ha"; cross[[595]]$mother <- "hA_Ra"
  cross[[596]]$father <- "hA_ha"; cross[[596]]$mother <- "hA_BA"
  cross[[597]]$father <- "hA_ha"; cross[[597]]$mother <- "hA_Ba"
  cross[[598]]$father <- "hA_ha"; cross[[598]]$mother <- "ha_ha"
  cross[[599]]$father <- "hA_ha"; cross[[599]]$mother <- "ha_RA"
  cross[[600]]$father <- "hA_ha"; cross[[600]]$mother <- "ha_Ra"
  cross[[601]]$father <- "hA_ha"; cross[[601]]$mother <- "ha_BA"
  cross[[602]]$father <- "hA_ha"; cross[[602]]$mother <- "ha_Ba"
  cross[[603]]$father <- "hA_ha"; cross[[603]]$mother <- "RA_RA"
  cross[[604]]$father <- "hA_ha"; cross[[604]]$mother <- "RA_Ra"
  cross[[605]]$father <- "hA_ha"; cross[[605]]$mother <- "RA_BA"
  cross[[606]]$father <- "hA_ha"; cross[[606]]$mother <- "RA_Ba"
  cross[[607]]$father <- "hA_ha"; cross[[607]]$mother <- "Ra_Ra"
  cross[[608]]$father <- "hA_ha"; cross[[608]]$mother <- "Ra_BA"
  cross[[609]]$father <- "hA_ha"; cross[[609]]$mother <- "Ra_Ba"
  cross[[610]]$father <- "hA_ha"; cross[[610]]$mother <- "BA_BA"
  cross[[611]]$father <- "hA_ha"; cross[[611]]$mother <- "BA_Ba"
  cross[[612]]$father <- "hA_ha"; cross[[612]]$mother <- "Ba_Ba"

  cross[[613]]$father <- "hA_RA"; cross[[613]]$mother <- "HA_HA"
  cross[[614]]$father <- "hA_RA"; cross[[614]]$mother <- "HA_Ha"
  cross[[615]]$father <- "hA_RA"; cross[[615]]$mother <- "HA_hA"
  cross[[616]]$father <- "hA_RA"; cross[[616]]$mother <- "HA_ha"
  cross[[617]]$father <- "hA_RA"; cross[[617]]$mother <- "HA_RA"
  cross[[618]]$father <- "hA_RA"; cross[[618]]$mother <- "HA_Ra"
  cross[[619]]$father <- "hA_RA"; cross[[619]]$mother <- "HA_BA"
  cross[[620]]$father <- "hA_RA"; cross[[620]]$mother <- "HA_Ba"
  cross[[621]]$father <- "hA_RA"; cross[[621]]$mother <- "Ha_Ha"
  cross[[622]]$father <- "hA_RA"; cross[[622]]$mother <- "Ha_hA"
  cross[[623]]$father <- "hA_RA"; cross[[623]]$mother <- "Ha_ha"
  cross[[624]]$father <- "hA_RA"; cross[[624]]$mother <- "Ha_RA"
  cross[[625]]$father <- "hA_RA"; cross[[625]]$mother <- "Ha_Ra"
  cross[[626]]$father <- "hA_RA"; cross[[626]]$mother <- "Ha_BA"
  cross[[627]]$father <- "hA_RA"; cross[[627]]$mother <- "Ha_Ba"
  cross[[628]]$father <- "hA_RA"; cross[[628]]$mother <- "hA_hA"
  cross[[629]]$father <- "hA_RA"; cross[[629]]$mother <- "hA_ha"
  cross[[630]]$father <- "hA_RA"; cross[[630]]$mother <- "hA_RA"
  cross[[631]]$father <- "hA_RA"; cross[[631]]$mother <- "hA_Ra"
  cross[[632]]$father <- "hA_RA"; cross[[632]]$mother <- "hA_BA"
  cross[[633]]$father <- "hA_RA"; cross[[633]]$mother <- "hA_Ba"
  cross[[634]]$father <- "hA_RA"; cross[[634]]$mother <- "ha_ha"
  cross[[635]]$father <- "hA_RA"; cross[[635]]$mother <- "ha_RA"
  cross[[636]]$father <- "hA_RA"; cross[[636]]$mother <- "ha_Ra"
  cross[[637]]$father <- "hA_RA"; cross[[637]]$mother <- "ha_BA"
  cross[[638]]$father <- "hA_RA"; cross[[638]]$mother <- "ha_Ba"
  cross[[639]]$father <- "hA_RA"; cross[[639]]$mother <- "RA_RA"
  cross[[640]]$father <- "hA_RA"; cross[[640]]$mother <- "RA_Ra"
  cross[[641]]$father <- "hA_RA"; cross[[641]]$mother <- "RA_BA"
  cross[[642]]$father <- "hA_RA"; cross[[642]]$mother <- "RA_Ba"
  cross[[643]]$father <- "hA_RA"; cross[[643]]$mother <- "Ra_Ra"
  cross[[644]]$father <- "hA_RA"; cross[[644]]$mother <- "Ra_BA"
  cross[[645]]$father <- "hA_RA"; cross[[645]]$mother <- "Ra_Ba"
  cross[[646]]$father <- "hA_RA"; cross[[646]]$mother <- "BA_BA"
  cross[[647]]$father <- "hA_RA"; cross[[647]]$mother <- "BA_Ba"
  cross[[648]]$father <- "hA_RA"; cross[[648]]$mother <- "Ba_Ba"

  cross[[649]]$father <- "hA_Ra"; cross[[649]]$mother <- "HA_HA"
  cross[[650]]$father <- "hA_Ra"; cross[[650]]$mother <- "HA_Ha"
  cross[[651]]$father <- "hA_Ra"; cross[[651]]$mother <- "HA_hA"
  cross[[652]]$father <- "hA_Ra"; cross[[652]]$mother <- "HA_ha"
  cross[[653]]$father <- "hA_Ra"; cross[[653]]$mother <- "HA_RA"
  cross[[654]]$father <- "hA_Ra"; cross[[654]]$mother <- "HA_Ra"
  cross[[655]]$father <- "hA_Ra"; cross[[655]]$mother <- "HA_BA"
  cross[[656]]$father <- "hA_Ra"; cross[[656]]$mother <- "HA_Ba"
  cross[[657]]$father <- "hA_Ra"; cross[[657]]$mother <- "Ha_Ha"
  cross[[658]]$father <- "hA_Ra"; cross[[658]]$mother <- "Ha_hA"
  cross[[659]]$father <- "hA_Ra"; cross[[659]]$mother <- "Ha_ha"
  cross[[660]]$father <- "hA_Ra"; cross[[660]]$mother <- "Ha_RA"
  cross[[661]]$father <- "hA_Ra"; cross[[661]]$mother <- "Ha_Ra"
  cross[[662]]$father <- "hA_Ra"; cross[[662]]$mother <- "Ha_BA"
  cross[[663]]$father <- "hA_Ra"; cross[[663]]$mother <- "Ha_Ba"
  cross[[664]]$father <- "hA_Ra"; cross[[664]]$mother <- "hA_hA"
  cross[[665]]$father <- "hA_Ra"; cross[[665]]$mother <- "hA_ha"
  cross[[666]]$father <- "hA_Ra"; cross[[666]]$mother <- "hA_RA"
  cross[[667]]$father <- "hA_Ra"; cross[[667]]$mother <- "hA_Ra"
  cross[[668]]$father <- "hA_Ra"; cross[[668]]$mother <- "hA_BA"
  cross[[669]]$father <- "hA_Ra"; cross[[669]]$mother <- "hA_Ba"
  cross[[670]]$father <- "hA_Ra"; cross[[670]]$mother <- "ha_ha"
  cross[[671]]$father <- "hA_Ra"; cross[[671]]$mother <- "ha_RA"
  cross[[672]]$father <- "hA_Ra"; cross[[672]]$mother <- "ha_Ra"
  cross[[673]]$father <- "hA_Ra"; cross[[673]]$mother <- "ha_BA"
  cross[[674]]$father <- "hA_Ra"; cross[[674]]$mother <- "ha_Ba"
  cross[[675]]$father <- "hA_Ra"; cross[[675]]$mother <- "RA_RA"
  cross[[676]]$father <- "hA_Ra"; cross[[676]]$mother <- "RA_Ra"
  cross[[677]]$father <- "hA_Ra"; cross[[677]]$mother <- "RA_BA"
  cross[[678]]$father <- "hA_Ra"; cross[[678]]$mother <- "RA_Ba"
  cross[[679]]$father <- "hA_Ra"; cross[[679]]$mother <- "Ra_Ra"
  cross[[680]]$father <- "hA_Ra"; cross[[680]]$mother <- "Ra_BA"
  cross[[681]]$father <- "hA_Ra"; cross[[681]]$mother <- "Ra_Ba"
  cross[[682]]$father <- "hA_Ra"; cross[[682]]$mother <- "BA_BA"
  cross[[683]]$father <- "hA_Ra"; cross[[683]]$mother <- "BA_Ba"
  cross[[684]]$father <- "hA_Ra"; cross[[684]]$mother <- "Ba_Ba"

  cross[[685]]$father <- "hA_BA"; cross[[685]]$mother <- "HA_HA"
  cross[[686]]$father <- "hA_BA"; cross[[686]]$mother <- "HA_Ha"
  cross[[687]]$father <- "hA_BA"; cross[[687]]$mother <- "HA_hA"
  cross[[688]]$father <- "hA_BA"; cross[[688]]$mother <- "HA_ha"
  cross[[689]]$father <- "hA_BA"; cross[[689]]$mother <- "HA_RA"
  cross[[690]]$father <- "hA_BA"; cross[[690]]$mother <- "HA_Ra"
  cross[[691]]$father <- "hA_BA"; cross[[691]]$mother <- "HA_BA"
  cross[[692]]$father <- "hA_BA"; cross[[692]]$mother <- "HA_Ba"
  cross[[693]]$father <- "hA_BA"; cross[[693]]$mother <- "Ha_Ha"
  cross[[694]]$father <- "hA_BA"; cross[[694]]$mother <- "Ha_hA"
  cross[[695]]$father <- "hA_BA"; cross[[695]]$mother <- "Ha_ha"
  cross[[696]]$father <- "hA_BA"; cross[[696]]$mother <- "Ha_RA"
  cross[[697]]$father <- "hA_BA"; cross[[697]]$mother <- "Ha_Ra"
  cross[[698]]$father <- "hA_BA"; cross[[698]]$mother <- "Ha_BA"
  cross[[699]]$father <- "hA_BA"; cross[[699]]$mother <- "Ha_Ba"
  cross[[700]]$father <- "hA_BA"; cross[[700]]$mother <- "hA_hA"
  cross[[701]]$father <- "hA_BA"; cross[[701]]$mother <- "hA_ha"
  cross[[702]]$father <- "hA_BA"; cross[[702]]$mother <- "hA_RA"
  cross[[703]]$father <- "hA_BA"; cross[[703]]$mother <- "hA_Ra"
  cross[[704]]$father <- "hA_BA"; cross[[704]]$mother <- "hA_BA"
  cross[[705]]$father <- "hA_BA"; cross[[705]]$mother <- "hA_Ba"
  cross[[706]]$father <- "hA_BA"; cross[[706]]$mother <- "ha_ha"
  cross[[707]]$father <- "hA_BA"; cross[[707]]$mother <- "ha_RA"
  cross[[708]]$father <- "hA_BA"; cross[[708]]$mother <- "ha_Ra"
  cross[[709]]$father <- "hA_BA"; cross[[709]]$mother <- "ha_BA"
  cross[[710]]$father <- "hA_BA"; cross[[710]]$mother <- "ha_Ba"
  cross[[711]]$father <- "hA_BA"; cross[[711]]$mother <- "RA_RA"
  cross[[712]]$father <- "hA_BA"; cross[[712]]$mother <- "RA_Ra"
  cross[[713]]$father <- "hA_BA"; cross[[713]]$mother <- "RA_BA"
  cross[[714]]$father <- "hA_BA"; cross[[714]]$mother <- "RA_Ba"
  cross[[715]]$father <- "hA_BA"; cross[[715]]$mother <- "Ra_Ra"
  cross[[716]]$father <- "hA_BA"; cross[[716]]$mother <- "Ra_BA"
  cross[[717]]$father <- "hA_BA"; cross[[717]]$mother <- "Ra_Ba"
  cross[[718]]$father <- "hA_BA"; cross[[718]]$mother <- "BA_BA"
  cross[[719]]$father <- "hA_BA"; cross[[719]]$mother <- "BA_Ba"
  cross[[720]]$father <- "hA_BA"; cross[[720]]$mother <- "Ba_Ba"

  cross[[721]]$father <- "hA_Ba"; cross[[721]]$mother <- "HA_HA"
  cross[[722]]$father <- "hA_Ba"; cross[[722]]$mother <- "HA_Ha"
  cross[[723]]$father <- "hA_Ba"; cross[[723]]$mother <- "HA_hA"
  cross[[724]]$father <- "hA_Ba"; cross[[724]]$mother <- "HA_ha"
  cross[[725]]$father <- "hA_Ba"; cross[[725]]$mother <- "HA_RA"
  cross[[726]]$father <- "hA_Ba"; cross[[726]]$mother <- "HA_Ra"
  cross[[727]]$father <- "hA_Ba"; cross[[727]]$mother <- "HA_BA"
  cross[[728]]$father <- "hA_Ba"; cross[[728]]$mother <- "HA_Ba"
  cross[[729]]$father <- "hA_Ba"; cross[[729]]$mother <- "Ha_Ha"
  cross[[730]]$father <- "hA_Ba"; cross[[730]]$mother <- "Ha_hA"
  cross[[731]]$father <- "hA_Ba"; cross[[731]]$mother <- "Ha_ha"
  cross[[732]]$father <- "hA_Ba"; cross[[732]]$mother <- "Ha_RA"
  cross[[733]]$father <- "hA_Ba"; cross[[733]]$mother <- "Ha_Ra"
  cross[[734]]$father <- "hA_Ba"; cross[[734]]$mother <- "Ha_BA"
  cross[[735]]$father <- "hA_Ba"; cross[[735]]$mother <- "Ha_Ba"
  cross[[736]]$father <- "hA_Ba"; cross[[736]]$mother <- "hA_hA"
  cross[[737]]$father <- "hA_Ba"; cross[[737]]$mother <- "hA_ha"
  cross[[738]]$father <- "hA_Ba"; cross[[738]]$mother <- "hA_RA"
  cross[[739]]$father <- "hA_Ba"; cross[[739]]$mother <- "hA_Ra"
  cross[[740]]$father <- "hA_Ba"; cross[[740]]$mother <- "hA_BA"
  cross[[741]]$father <- "hA_Ba"; cross[[741]]$mother <- "hA_Ba"
  cross[[742]]$father <- "hA_Ba"; cross[[742]]$mother <- "ha_ha"
  cross[[743]]$father <- "hA_Ba"; cross[[743]]$mother <- "ha_RA"
  cross[[744]]$father <- "hA_Ba"; cross[[744]]$mother <- "ha_Ra"
  cross[[745]]$father <- "hA_Ba"; cross[[745]]$mother <- "ha_BA"
  cross[[746]]$father <- "hA_Ba"; cross[[746]]$mother <- "ha_Ba"
  cross[[747]]$father <- "hA_Ba"; cross[[747]]$mother <- "RA_RA"
  cross[[748]]$father <- "hA_Ba"; cross[[748]]$mother <- "RA_Ra"
  cross[[749]]$father <- "hA_Ba"; cross[[749]]$mother <- "RA_BA"
  cross[[750]]$father <- "hA_Ba"; cross[[750]]$mother <- "RA_Ba"
  cross[[751]]$father <- "hA_Ba"; cross[[751]]$mother <- "Ra_Ra"
  cross[[752]]$father <- "hA_Ba"; cross[[752]]$mother <- "Ra_BA"
  cross[[753]]$father <- "hA_Ba"; cross[[753]]$mother <- "Ra_Ba"
  cross[[754]]$father <- "hA_Ba"; cross[[754]]$mother <- "BA_BA"
  cross[[755]]$father <- "hA_Ba"; cross[[755]]$mother <- "BA_Ba"
  cross[[756]]$father <- "hA_Ba"; cross[[756]]$mother <- "Ba_Ba"

  cross[[757]]$father <- "ha_ha"; cross[[757]]$mother <- "HA_HA"
  cross[[758]]$father <- "ha_ha"; cross[[758]]$mother <- "HA_Ha"
  cross[[759]]$father <- "ha_ha"; cross[[759]]$mother <- "HA_hA"
  cross[[760]]$father <- "ha_ha"; cross[[760]]$mother <- "HA_ha"
  cross[[761]]$father <- "ha_ha"; cross[[761]]$mother <- "HA_RA"
  cross[[762]]$father <- "ha_ha"; cross[[762]]$mother <- "HA_Ra"
  cross[[763]]$father <- "ha_ha"; cross[[763]]$mother <- "HA_BA"
  cross[[764]]$father <- "ha_ha"; cross[[764]]$mother <- "HA_Ba"
  cross[[765]]$father <- "ha_ha"; cross[[765]]$mother <- "Ha_Ha"
  cross[[766]]$father <- "ha_ha"; cross[[766]]$mother <- "Ha_hA"
  cross[[767]]$father <- "ha_ha"; cross[[767]]$mother <- "Ha_ha"
  cross[[768]]$father <- "ha_ha"; cross[[768]]$mother <- "Ha_RA"
  cross[[769]]$father <- "ha_ha"; cross[[769]]$mother <- "Ha_Ra"
  cross[[770]]$father <- "ha_ha"; cross[[770]]$mother <- "Ha_BA"
  cross[[771]]$father <- "ha_ha"; cross[[771]]$mother <- "Ha_Ba"
  cross[[772]]$father <- "ha_ha"; cross[[772]]$mother <- "hA_hA"
  cross[[773]]$father <- "ha_ha"; cross[[773]]$mother <- "hA_ha"
  cross[[774]]$father <- "ha_ha"; cross[[774]]$mother <- "hA_RA"
  cross[[775]]$father <- "ha_ha"; cross[[775]]$mother <- "hA_Ra"
  cross[[776]]$father <- "ha_ha"; cross[[776]]$mother <- "hA_BA"
  cross[[777]]$father <- "ha_ha"; cross[[777]]$mother <- "hA_Ba"
  cross[[778]]$father <- "ha_ha"; cross[[778]]$mother <- "ha_ha"
  cross[[779]]$father <- "ha_ha"; cross[[779]]$mother <- "ha_RA"
  cross[[780]]$father <- "ha_ha"; cross[[780]]$mother <- "ha_Ra"
  cross[[781]]$father <- "ha_ha"; cross[[781]]$mother <- "ha_BA"
  cross[[782]]$father <- "ha_ha"; cross[[782]]$mother <- "ha_Ba"
  cross[[783]]$father <- "ha_ha"; cross[[783]]$mother <- "RA_RA"
  cross[[784]]$father <- "ha_ha"; cross[[784]]$mother <- "RA_Ra"
  cross[[785]]$father <- "ha_ha"; cross[[785]]$mother <- "RA_BA"
  cross[[786]]$father <- "ha_ha"; cross[[786]]$mother <- "RA_Ba"
  cross[[787]]$father <- "ha_ha"; cross[[787]]$mother <- "Ra_Ra"
  cross[[788]]$father <- "ha_ha"; cross[[788]]$mother <- "Ra_BA"
  cross[[789]]$father <- "ha_ha"; cross[[789]]$mother <- "Ra_Ba"
  cross[[790]]$father <- "ha_ha"; cross[[790]]$mother <- "BA_BA"
  cross[[791]]$father <- "ha_ha"; cross[[791]]$mother <- "BA_Ba"
  cross[[792]]$father <- "ha_ha"; cross[[792]]$mother <- "Ba_Ba"

  cross[[793]]$father <- "ha_RA"; cross[[793]]$mother <- "HA_HA"
  cross[[794]]$father <- "ha_RA"; cross[[794]]$mother <- "HA_Ha"
  cross[[795]]$father <- "ha_RA"; cross[[795]]$mother <- "HA_hA"
  cross[[796]]$father <- "ha_RA"; cross[[796]]$mother <- "HA_ha"
  cross[[797]]$father <- "ha_RA"; cross[[797]]$mother <- "HA_RA"
  cross[[798]]$father <- "ha_RA"; cross[[798]]$mother <- "HA_Ra"
  cross[[799]]$father <- "ha_RA"; cross[[799]]$mother <- "HA_BA"
  cross[[800]]$father <- "ha_RA"; cross[[800]]$mother <- "HA_Ba"
  cross[[801]]$father <- "ha_RA"; cross[[801]]$mother <- "Ha_Ha"
  cross[[802]]$father <- "ha_RA"; cross[[802]]$mother <- "Ha_hA"
  cross[[803]]$father <- "ha_RA"; cross[[803]]$mother <- "Ha_ha"
  cross[[804]]$father <- "ha_RA"; cross[[804]]$mother <- "Ha_RA"
  cross[[805]]$father <- "ha_RA"; cross[[805]]$mother <- "Ha_Ra"
  cross[[806]]$father <- "ha_RA"; cross[[806]]$mother <- "Ha_BA"
  cross[[807]]$father <- "ha_RA"; cross[[807]]$mother <- "Ha_Ba"
  cross[[808]]$father <- "ha_RA"; cross[[808]]$mother <- "hA_hA"
  cross[[809]]$father <- "ha_RA"; cross[[809]]$mother <- "hA_ha"
  cross[[810]]$father <- "ha_RA"; cross[[810]]$mother <- "hA_RA"
  cross[[811]]$father <- "ha_RA"; cross[[811]]$mother <- "hA_Ra"
  cross[[812]]$father <- "ha_RA"; cross[[812]]$mother <- "hA_BA"
  cross[[813]]$father <- "ha_RA"; cross[[813]]$mother <- "hA_Ba"
  cross[[814]]$father <- "ha_RA"; cross[[814]]$mother <- "ha_ha"
  cross[[815]]$father <- "ha_RA"; cross[[815]]$mother <- "ha_RA"
  cross[[816]]$father <- "ha_RA"; cross[[816]]$mother <- "ha_Ra"
  cross[[817]]$father <- "ha_RA"; cross[[817]]$mother <- "ha_BA"
  cross[[818]]$father <- "ha_RA"; cross[[818]]$mother <- "ha_Ba"
  cross[[819]]$father <- "ha_RA"; cross[[819]]$mother <- "RA_RA"
  cross[[820]]$father <- "ha_RA"; cross[[820]]$mother <- "RA_Ra"
  cross[[821]]$father <- "ha_RA"; cross[[821]]$mother <- "RA_BA"
  cross[[822]]$father <- "ha_RA"; cross[[822]]$mother <- "RA_Ba"
  cross[[823]]$father <- "ha_RA"; cross[[823]]$mother <- "Ra_Ra"
  cross[[824]]$father <- "ha_RA"; cross[[824]]$mother <- "Ra_BA"
  cross[[825]]$father <- "ha_RA"; cross[[825]]$mother <- "Ra_Ba"
  cross[[826]]$father <- "ha_RA"; cross[[826]]$mother <- "BA_BA"
  cross[[827]]$father <- "ha_RA"; cross[[827]]$mother <- "BA_Ba"
  cross[[828]]$father <- "ha_RA"; cross[[828]]$mother <- "Ba_Ba"

  cross[[829]]$father <- "ha_Ra"; cross[[829]]$mother <- "HA_HA"
  cross[[830]]$father <- "ha_Ra"; cross[[830]]$mother <- "HA_Ha"
  cross[[831]]$father <- "ha_Ra"; cross[[831]]$mother <- "HA_hA"
  cross[[832]]$father <- "ha_Ra"; cross[[832]]$mother <- "HA_ha"
  cross[[833]]$father <- "ha_Ra"; cross[[833]]$mother <- "HA_RA"
  cross[[834]]$father <- "ha_Ra"; cross[[834]]$mother <- "HA_Ra"
  cross[[835]]$father <- "ha_Ra"; cross[[835]]$mother <- "HA_BA"
  cross[[836]]$father <- "ha_Ra"; cross[[836]]$mother <- "HA_Ba"
  cross[[837]]$father <- "ha_Ra"; cross[[837]]$mother <- "Ha_Ha"
  cross[[838]]$father <- "ha_Ra"; cross[[838]]$mother <- "Ha_hA"
  cross[[839]]$father <- "ha_Ra"; cross[[839]]$mother <- "Ha_ha"
  cross[[840]]$father <- "ha_Ra"; cross[[840]]$mother <- "Ha_RA"
  cross[[841]]$father <- "ha_Ra"; cross[[841]]$mother <- "Ha_Ra"
  cross[[842]]$father <- "ha_Ra"; cross[[842]]$mother <- "Ha_BA"
  cross[[843]]$father <- "ha_Ra"; cross[[843]]$mother <- "Ha_Ba"
  cross[[844]]$father <- "ha_Ra"; cross[[844]]$mother <- "hA_hA"
  cross[[845]]$father <- "ha_Ra"; cross[[845]]$mother <- "hA_ha"
  cross[[846]]$father <- "ha_Ra"; cross[[846]]$mother <- "hA_RA"
  cross[[847]]$father <- "ha_Ra"; cross[[847]]$mother <- "hA_Ra"
  cross[[848]]$father <- "ha_Ra"; cross[[848]]$mother <- "hA_BA"
  cross[[849]]$father <- "ha_Ra"; cross[[849]]$mother <- "hA_Ba"
  cross[[850]]$father <- "ha_Ra"; cross[[850]]$mother <- "ha_ha"
  cross[[851]]$father <- "ha_Ra"; cross[[851]]$mother <- "ha_RA"
  cross[[852]]$father <- "ha_Ra"; cross[[852]]$mother <- "ha_Ra"
  cross[[853]]$father <- "ha_Ra"; cross[[853]]$mother <- "ha_BA"
  cross[[854]]$father <- "ha_Ra"; cross[[854]]$mother <- "ha_Ba"
  cross[[855]]$father <- "ha_Ra"; cross[[855]]$mother <- "RA_RA"
  cross[[856]]$father <- "ha_Ra"; cross[[856]]$mother <- "RA_Ra"
  cross[[857]]$father <- "ha_Ra"; cross[[857]]$mother <- "RA_BA"
  cross[[858]]$father <- "ha_Ra"; cross[[858]]$mother <- "RA_Ba"
  cross[[859]]$father <- "ha_Ra"; cross[[859]]$mother <- "Ra_Ra"
  cross[[860]]$father <- "ha_Ra"; cross[[860]]$mother <- "Ra_BA"
  cross[[861]]$father <- "ha_Ra"; cross[[861]]$mother <- "Ra_Ba"
  cross[[862]]$father <- "ha_Ra"; cross[[862]]$mother <- "BA_BA"
  cross[[863]]$father <- "ha_Ra"; cross[[863]]$mother <- "BA_Ba"
  cross[[864]]$father <- "ha_Ra"; cross[[864]]$mother <- "Ba_Ba"

  cross[[865]]$father <- "ha_BA"; cross[[865]]$mother <- "HA_HA"
  cross[[866]]$father <- "ha_BA"; cross[[866]]$mother <- "HA_Ha"
  cross[[867]]$father <- "ha_BA"; cross[[867]]$mother <- "HA_hA"
  cross[[868]]$father <- "ha_BA"; cross[[868]]$mother <- "HA_ha"
  cross[[869]]$father <- "ha_BA"; cross[[869]]$mother <- "HA_RA"
  cross[[870]]$father <- "ha_BA"; cross[[870]]$mother <- "HA_Ra"
  cross[[871]]$father <- "ha_BA"; cross[[871]]$mother <- "HA_BA"
  cross[[872]]$father <- "ha_BA"; cross[[872]]$mother <- "HA_Ba"
  cross[[873]]$father <- "ha_BA"; cross[[873]]$mother <- "Ha_Ha"
  cross[[874]]$father <- "ha_BA"; cross[[874]]$mother <- "Ha_hA"
  cross[[875]]$father <- "ha_BA"; cross[[875]]$mother <- "Ha_ha"
  cross[[876]]$father <- "ha_BA"; cross[[876]]$mother <- "Ha_RA"
  cross[[877]]$father <- "ha_BA"; cross[[877]]$mother <- "Ha_Ra"
  cross[[878]]$father <- "ha_BA"; cross[[878]]$mother <- "Ha_BA"
  cross[[879]]$father <- "ha_BA"; cross[[879]]$mother <- "Ha_Ba"
  cross[[880]]$father <- "ha_BA"; cross[[880]]$mother <- "hA_hA"
  cross[[881]]$father <- "ha_BA"; cross[[881]]$mother <- "hA_ha"
  cross[[882]]$father <- "ha_BA"; cross[[882]]$mother <- "hA_RA"
  cross[[883]]$father <- "ha_BA"; cross[[883]]$mother <- "hA_Ra"
  cross[[884]]$father <- "ha_BA"; cross[[884]]$mother <- "hA_BA"
  cross[[885]]$father <- "ha_BA"; cross[[885]]$mother <- "hA_Ba"
  cross[[886]]$father <- "ha_BA"; cross[[886]]$mother <- "ha_ha"
  cross[[887]]$father <- "ha_BA"; cross[[887]]$mother <- "ha_RA"
  cross[[888]]$father <- "ha_BA"; cross[[888]]$mother <- "ha_Ra"
  cross[[889]]$father <- "ha_BA"; cross[[889]]$mother <- "ha_BA"
  cross[[890]]$father <- "ha_BA"; cross[[890]]$mother <- "ha_Ba"
  cross[[891]]$father <- "ha_BA"; cross[[891]]$mother <- "RA_RA"
  cross[[892]]$father <- "ha_BA"; cross[[892]]$mother <- "RA_Ra"
  cross[[893]]$father <- "ha_BA"; cross[[893]]$mother <- "RA_BA"
  cross[[894]]$father <- "ha_BA"; cross[[894]]$mother <- "RA_Ba"
  cross[[895]]$father <- "ha_BA"; cross[[895]]$mother <- "Ra_Ra"
  cross[[896]]$father <- "ha_BA"; cross[[896]]$mother <- "Ra_BA"
  cross[[897]]$father <- "ha_BA"; cross[[897]]$mother <- "Ra_Ba"
  cross[[898]]$father <- "ha_BA"; cross[[898]]$mother <- "BA_BA"
  cross[[899]]$father <- "ha_BA"; cross[[899]]$mother <- "BA_Ba"
  cross[[900]]$father <- "ha_BA"; cross[[900]]$mother <- "Ba_Ba"

  cross[[901]]$father <- "ha_Ba"; cross[[901]]$mother <- "HA_HA"
  cross[[902]]$father <- "ha_Ba"; cross[[902]]$mother <- "HA_Ha"
  cross[[903]]$father <- "ha_Ba"; cross[[903]]$mother <- "HA_hA"
  cross[[904]]$father <- "ha_Ba"; cross[[904]]$mother <- "HA_ha"
  cross[[905]]$father <- "ha_Ba"; cross[[905]]$mother <- "HA_RA"
  cross[[906]]$father <- "ha_Ba"; cross[[906]]$mother <- "HA_Ra"
  cross[[907]]$father <- "ha_Ba"; cross[[907]]$mother <- "HA_BA"
  cross[[908]]$father <- "ha_Ba"; cross[[908]]$mother <- "HA_Ba"
  cross[[909]]$father <- "ha_Ba"; cross[[909]]$mother <- "Ha_Ha"
  cross[[910]]$father <- "ha_Ba"; cross[[910]]$mother <- "Ha_hA"
  cross[[911]]$father <- "ha_Ba"; cross[[911]]$mother <- "Ha_ha"
  cross[[912]]$father <- "ha_Ba"; cross[[912]]$mother <- "Ha_RA"
  cross[[913]]$father <- "ha_Ba"; cross[[913]]$mother <- "Ha_Ra"
  cross[[914]]$father <- "ha_Ba"; cross[[914]]$mother <- "Ha_BA"
  cross[[915]]$father <- "ha_Ba"; cross[[915]]$mother <- "Ha_Ba"
  cross[[916]]$father <- "ha_Ba"; cross[[916]]$mother <- "hA_hA"
  cross[[917]]$father <- "ha_Ba"; cross[[917]]$mother <- "hA_ha"
  cross[[918]]$father <- "ha_Ba"; cross[[918]]$mother <- "hA_RA"
  cross[[919]]$father <- "ha_Ba"; cross[[919]]$mother <- "hA_Ra"
  cross[[920]]$father <- "ha_Ba"; cross[[920]]$mother <- "hA_BA"
  cross[[921]]$father <- "ha_Ba"; cross[[921]]$mother <- "hA_Ba"
  cross[[922]]$father <- "ha_Ba"; cross[[922]]$mother <- "ha_ha"
  cross[[923]]$father <- "ha_Ba"; cross[[923]]$mother <- "ha_RA"
  cross[[924]]$father <- "ha_Ba"; cross[[924]]$mother <- "ha_Ra"
  cross[[925]]$father <- "ha_Ba"; cross[[925]]$mother <- "ha_BA"
  cross[[926]]$father <- "ha_Ba"; cross[[926]]$mother <- "ha_Ba"
  cross[[927]]$father <- "ha_Ba"; cross[[927]]$mother <- "RA_RA"
  cross[[928]]$father <- "ha_Ba"; cross[[928]]$mother <- "RA_Ra"
  cross[[929]]$father <- "ha_Ba"; cross[[929]]$mother <- "RA_BA"
  cross[[930]]$father <- "ha_Ba"; cross[[930]]$mother <- "RA_Ba"
  cross[[931]]$father <- "ha_Ba"; cross[[931]]$mother <- "Ra_Ra"
  cross[[932]]$father <- "ha_Ba"; cross[[932]]$mother <- "Ra_BA"
  cross[[933]]$father <- "ha_Ba"; cross[[933]]$mother <- "Ra_Ba"
  cross[[934]]$father <- "ha_Ba"; cross[[934]]$mother <- "BA_BA"
  cross[[935]]$father <- "ha_Ba"; cross[[935]]$mother <- "BA_Ba"
  cross[[936]]$father <- "ha_Ba"; cross[[936]]$mother <- "Ba_Ba"

  cross[[937]]$father <- "RA_RA"; cross[[937]]$mother <- "HA_HA"
  cross[[938]]$father <- "RA_RA"; cross[[938]]$mother <- "HA_Ha"
  cross[[939]]$father <- "RA_RA"; cross[[939]]$mother <- "HA_hA"
  cross[[940]]$father <- "RA_RA"; cross[[940]]$mother <- "HA_ha"
  cross[[941]]$father <- "RA_RA"; cross[[941]]$mother <- "HA_RA"
  cross[[942]]$father <- "RA_RA"; cross[[942]]$mother <- "HA_Ra"
  cross[[943]]$father <- "RA_RA"; cross[[943]]$mother <- "HA_BA"
  cross[[944]]$father <- "RA_RA"; cross[[944]]$mother <- "HA_Ba"
  cross[[945]]$father <- "RA_RA"; cross[[945]]$mother <- "Ha_Ha"
  cross[[946]]$father <- "RA_RA"; cross[[946]]$mother <- "Ha_hA"
  cross[[947]]$father <- "RA_RA"; cross[[947]]$mother <- "Ha_ha"
  cross[[948]]$father <- "RA_RA"; cross[[948]]$mother <- "Ha_RA"
  cross[[949]]$father <- "RA_RA"; cross[[949]]$mother <- "Ha_Ra"
  cross[[950]]$father <- "RA_RA"; cross[[950]]$mother <- "Ha_BA"
  cross[[951]]$father <- "RA_RA"; cross[[951]]$mother <- "Ha_Ba"
  cross[[952]]$father <- "RA_RA"; cross[[952]]$mother <- "hA_hA"
  cross[[953]]$father <- "RA_RA"; cross[[953]]$mother <- "hA_ha"
  cross[[954]]$father <- "RA_RA"; cross[[954]]$mother <- "hA_RA"
  cross[[955]]$father <- "RA_RA"; cross[[955]]$mother <- "hA_Ra"
  cross[[956]]$father <- "RA_RA"; cross[[956]]$mother <- "hA_BA"
  cross[[957]]$father <- "RA_RA"; cross[[957]]$mother <- "hA_Ba"
  cross[[958]]$father <- "RA_RA"; cross[[958]]$mother <- "ha_ha"
  cross[[959]]$father <- "RA_RA"; cross[[959]]$mother <- "ha_RA"
  cross[[960]]$father <- "RA_RA"; cross[[960]]$mother <- "ha_Ra"
  cross[[961]]$father <- "RA_RA"; cross[[961]]$mother <- "ha_BA"
  cross[[962]]$father <- "RA_RA"; cross[[962]]$mother <- "ha_Ba"
  cross[[963]]$father <- "RA_RA"; cross[[963]]$mother <- "RA_RA"
  cross[[964]]$father <- "RA_RA"; cross[[964]]$mother <- "RA_Ra"
  cross[[965]]$father <- "RA_RA"; cross[[965]]$mother <- "RA_BA"
  cross[[966]]$father <- "RA_RA"; cross[[966]]$mother <- "RA_Ba"
  cross[[967]]$father <- "RA_RA"; cross[[967]]$mother <- "Ra_Ra"
  cross[[968]]$father <- "RA_RA"; cross[[968]]$mother <- "Ra_BA"
  cross[[969]]$father <- "RA_RA"; cross[[969]]$mother <- "Ra_Ba"
  cross[[970]]$father <- "RA_RA"; cross[[970]]$mother <- "BA_BA"
  cross[[971]]$father <- "RA_RA"; cross[[971]]$mother <- "BA_Ba"
  cross[[972]]$father <- "RA_RA"; cross[[972]]$mother <- "Ba_Ba"

  cross[[973]]$father <- "RA_Ra"; cross[[973]]$mother <- "HA_HA"
  cross[[974]]$father <- "RA_Ra"; cross[[974]]$mother <- "HA_Ha"
  cross[[975]]$father <- "RA_Ra"; cross[[975]]$mother <- "HA_hA"
  cross[[976]]$father <- "RA_Ra"; cross[[976]]$mother <- "HA_ha"
  cross[[977]]$father <- "RA_Ra"; cross[[977]]$mother <- "HA_RA"
  cross[[978]]$father <- "RA_Ra"; cross[[978]]$mother <- "HA_Ra"
  cross[[979]]$father <- "RA_Ra"; cross[[979]]$mother <- "HA_BA"
  cross[[980]]$father <- "RA_Ra"; cross[[980]]$mother <- "HA_Ba"
  cross[[981]]$father <- "RA_Ra"; cross[[981]]$mother <- "Ha_Ha"
  cross[[982]]$father <- "RA_Ra"; cross[[982]]$mother <- "Ha_hA"
  cross[[983]]$father <- "RA_Ra"; cross[[983]]$mother <- "Ha_ha"
  cross[[984]]$father <- "RA_Ra"; cross[[984]]$mother <- "Ha_RA"
  cross[[985]]$father <- "RA_Ra"; cross[[985]]$mother <- "Ha_Ra"
  cross[[986]]$father <- "RA_Ra"; cross[[986]]$mother <- "Ha_BA"
  cross[[987]]$father <- "RA_Ra"; cross[[987]]$mother <- "Ha_Ba"
  cross[[988]]$father <- "RA_Ra"; cross[[988]]$mother <- "hA_hA"
  cross[[989]]$father <- "RA_Ra"; cross[[989]]$mother <- "hA_ha"
  cross[[990]]$father <- "RA_Ra"; cross[[990]]$mother <- "hA_RA"
  cross[[991]]$father <- "RA_Ra"; cross[[991]]$mother <- "hA_Ra"
  cross[[992]]$father <- "RA_Ra"; cross[[992]]$mother <- "hA_BA"
  cross[[993]]$father <- "RA_Ra"; cross[[993]]$mother <- "hA_Ba"
  cross[[994]]$father <- "RA_Ra"; cross[[994]]$mother <- "ha_ha"
  cross[[995]]$father <- "RA_Ra"; cross[[995]]$mother <- "ha_RA"
  cross[[996]]$father <- "RA_Ra"; cross[[996]]$mother <- "ha_Ra"
  cross[[997]]$father <- "RA_Ra"; cross[[997]]$mother <- "ha_BA"
  cross[[998]]$father <- "RA_Ra"; cross[[998]]$mother <- "ha_Ba"
  cross[[999]]$father <- "RA_Ra"; cross[[999]]$mother <- "RA_RA"
  cross[[1000]]$father <- "RA_Ra"; cross[[1000]]$mother <- "RA_Ra"
  cross[[1001]]$father <- "RA_Ra"; cross[[1001]]$mother <- "RA_BA"
  cross[[1002]]$father <- "RA_Ra"; cross[[1002]]$mother <- "RA_Ba"
  cross[[1003]]$father <- "RA_Ra"; cross[[1003]]$mother <- "Ra_Ra"
  cross[[1004]]$father <- "RA_Ra"; cross[[1004]]$mother <- "Ra_BA"
  cross[[1005]]$father <- "RA_Ra"; cross[[1005]]$mother <- "Ra_Ba"
  cross[[1006]]$father <- "RA_Ra"; cross[[1006]]$mother <- "BA_BA"
  cross[[1007]]$father <- "RA_Ra"; cross[[1007]]$mother <- "BA_Ba"
  cross[[1008]]$father <- "RA_Ra"; cross[[1008]]$mother <- "Ba_Ba"

  cross[[1009]]$father <- "RA_BA"; cross[[1009]]$mother <- "HA_HA"
  cross[[1010]]$father <- "RA_BA"; cross[[1010]]$mother <- "HA_Ha"
  cross[[1011]]$father <- "RA_BA"; cross[[1011]]$mother <- "HA_hA"
  cross[[1012]]$father <- "RA_BA"; cross[[1012]]$mother <- "HA_ha"
  cross[[1013]]$father <- "RA_BA"; cross[[1013]]$mother <- "HA_RA"
  cross[[1014]]$father <- "RA_BA"; cross[[1014]]$mother <- "HA_Ra"
  cross[[1015]]$father <- "RA_BA"; cross[[1015]]$mother <- "HA_BA"
  cross[[1016]]$father <- "RA_BA"; cross[[1016]]$mother <- "HA_Ba"
  cross[[1017]]$father <- "RA_BA"; cross[[1017]]$mother <- "Ha_Ha"
  cross[[1018]]$father <- "RA_BA"; cross[[1018]]$mother <- "Ha_hA"
  cross[[1019]]$father <- "RA_BA"; cross[[1019]]$mother <- "Ha_ha"
  cross[[1020]]$father <- "RA_BA"; cross[[1020]]$mother <- "Ha_RA"
  cross[[1021]]$father <- "RA_BA"; cross[[1021]]$mother <- "Ha_Ra"
  cross[[1022]]$father <- "RA_BA"; cross[[1022]]$mother <- "Ha_BA"
  cross[[1023]]$father <- "RA_BA"; cross[[1023]]$mother <- "Ha_Ba"
  cross[[1024]]$father <- "RA_BA"; cross[[1024]]$mother <- "hA_hA"
  cross[[1025]]$father <- "RA_BA"; cross[[1025]]$mother <- "hA_ha"
  cross[[1026]]$father <- "RA_BA"; cross[[1026]]$mother <- "hA_RA"
  cross[[1027]]$father <- "RA_BA"; cross[[1027]]$mother <- "hA_Ra"
  cross[[1028]]$father <- "RA_BA"; cross[[1028]]$mother <- "hA_BA"
  cross[[1029]]$father <- "RA_BA"; cross[[1029]]$mother <- "hA_Ba"
  cross[[1030]]$father <- "RA_BA"; cross[[1030]]$mother <- "ha_ha"
  cross[[1031]]$father <- "RA_BA"; cross[[1031]]$mother <- "ha_RA"
  cross[[1032]]$father <- "RA_BA"; cross[[1032]]$mother <- "ha_Ra"
  cross[[1033]]$father <- "RA_BA"; cross[[1033]]$mother <- "ha_BA"
  cross[[1034]]$father <- "RA_BA"; cross[[1034]]$mother <- "ha_Ba"
  cross[[1035]]$father <- "RA_BA"; cross[[1035]]$mother <- "RA_RA"
  cross[[1036]]$father <- "RA_BA"; cross[[1036]]$mother <- "RA_Ra"
  cross[[1037]]$father <- "RA_BA"; cross[[1037]]$mother <- "RA_BA"
  cross[[1038]]$father <- "RA_BA"; cross[[1038]]$mother <- "RA_Ba"
  cross[[1039]]$father <- "RA_BA"; cross[[1039]]$mother <- "Ra_Ra"
  cross[[1040]]$father <- "RA_BA"; cross[[1040]]$mother <- "Ra_BA"
  cross[[1041]]$father <- "RA_BA"; cross[[1041]]$mother <- "Ra_Ba"
  cross[[1042]]$father <- "RA_BA"; cross[[1042]]$mother <- "BA_BA"
  cross[[1043]]$father <- "RA_BA"; cross[[1043]]$mother <- "BA_Ba"
  cross[[1044]]$father <- "RA_BA"; cross[[1044]]$mother <- "Ba_Ba"

  cross[[1045]]$father <- "RA_Ba"; cross[[1045]]$mother <- "HA_HA"
  cross[[1046]]$father <- "RA_Ba"; cross[[1046]]$mother <- "HA_Ha"
  cross[[1047]]$father <- "RA_Ba"; cross[[1047]]$mother <- "HA_hA"
  cross[[1048]]$father <- "RA_Ba"; cross[[1048]]$mother <- "HA_ha"
  cross[[1049]]$father <- "RA_Ba"; cross[[1049]]$mother <- "HA_RA"
  cross[[1050]]$father <- "RA_Ba"; cross[[1050]]$mother <- "HA_Ra"
  cross[[1051]]$father <- "RA_Ba"; cross[[1051]]$mother <- "HA_BA"
  cross[[1052]]$father <- "RA_Ba"; cross[[1052]]$mother <- "HA_Ba"
  cross[[1053]]$father <- "RA_Ba"; cross[[1053]]$mother <- "Ha_Ha"
  cross[[1054]]$father <- "RA_Ba"; cross[[1054]]$mother <- "Ha_hA"
  cross[[1055]]$father <- "RA_Ba"; cross[[1055]]$mother <- "Ha_ha"
  cross[[1056]]$father <- "RA_Ba"; cross[[1056]]$mother <- "Ha_RA"
  cross[[1057]]$father <- "RA_Ba"; cross[[1057]]$mother <- "Ha_Ra"
  cross[[1058]]$father <- "RA_Ba"; cross[[1058]]$mother <- "Ha_BA"
  cross[[1059]]$father <- "RA_Ba"; cross[[1059]]$mother <- "Ha_Ba"
  cross[[1060]]$father <- "RA_Ba"; cross[[1060]]$mother <- "hA_hA"
  cross[[1061]]$father <- "RA_Ba"; cross[[1061]]$mother <- "hA_ha"
  cross[[1062]]$father <- "RA_Ba"; cross[[1062]]$mother <- "hA_RA"
  cross[[1063]]$father <- "RA_Ba"; cross[[1063]]$mother <- "hA_Ra"
  cross[[1064]]$father <- "RA_Ba"; cross[[1064]]$mother <- "hA_BA"
  cross[[1065]]$father <- "RA_Ba"; cross[[1065]]$mother <- "hA_Ba"
  cross[[1066]]$father <- "RA_Ba"; cross[[1066]]$mother <- "ha_ha"
  cross[[1067]]$father <- "RA_Ba"; cross[[1067]]$mother <- "ha_RA"
  cross[[1068]]$father <- "RA_Ba"; cross[[1068]]$mother <- "ha_Ra"
  cross[[1069]]$father <- "RA_Ba"; cross[[1069]]$mother <- "ha_BA"
  cross[[1070]]$father <- "RA_Ba"; cross[[1070]]$mother <- "ha_Ba"
  cross[[1071]]$father <- "RA_Ba"; cross[[1071]]$mother <- "RA_RA"
  cross[[1072]]$father <- "RA_Ba"; cross[[1072]]$mother <- "RA_Ra"
  cross[[1073]]$father <- "RA_Ba"; cross[[1073]]$mother <- "RA_BA"
  cross[[1074]]$father <- "RA_Ba"; cross[[1074]]$mother <- "RA_Ba"
  cross[[1075]]$father <- "RA_Ba"; cross[[1075]]$mother <- "Ra_Ra"
  cross[[1076]]$father <- "RA_Ba"; cross[[1076]]$mother <- "Ra_BA"
  cross[[1077]]$father <- "RA_Ba"; cross[[1077]]$mother <- "Ra_Ba"
  cross[[1078]]$father <- "RA_Ba"; cross[[1078]]$mother <- "BA_BA"
  cross[[1079]]$father <- "RA_Ba"; cross[[1079]]$mother <- "BA_Ba"
  cross[[1080]]$father <- "RA_Ba"; cross[[1080]]$mother <- "Ba_Ba"

  cross[[1081]]$father <- "Ra_Ra"; cross[[1081]]$mother <- "HA_HA"
  cross[[1082]]$father <- "Ra_Ra"; cross[[1082]]$mother <- "HA_Ha"
  cross[[1083]]$father <- "Ra_Ra"; cross[[1083]]$mother <- "HA_hA"
  cross[[1084]]$father <- "Ra_Ra"; cross[[1084]]$mother <- "HA_ha"
  cross[[1085]]$father <- "Ra_Ra"; cross[[1085]]$mother <- "HA_RA"
  cross[[1086]]$father <- "Ra_Ra"; cross[[1086]]$mother <- "HA_Ra"
  cross[[1087]]$father <- "Ra_Ra"; cross[[1087]]$mother <- "HA_BA"
  cross[[1088]]$father <- "Ra_Ra"; cross[[1088]]$mother <- "HA_Ba"
  cross[[1089]]$father <- "Ra_Ra"; cross[[1089]]$mother <- "Ha_Ha"
  cross[[1090]]$father <- "Ra_Ra"; cross[[1090]]$mother <- "Ha_hA"
  cross[[1091]]$father <- "Ra_Ra"; cross[[1091]]$mother <- "Ha_ha"
  cross[[1092]]$father <- "Ra_Ra"; cross[[1092]]$mother <- "Ha_RA"
  cross[[1093]]$father <- "Ra_Ra"; cross[[1093]]$mother <- "Ha_Ra"
  cross[[1094]]$father <- "Ra_Ra"; cross[[1094]]$mother <- "Ha_BA"
  cross[[1095]]$father <- "Ra_Ra"; cross[[1095]]$mother <- "Ha_Ba"
  cross[[1096]]$father <- "Ra_Ra"; cross[[1096]]$mother <- "hA_hA"
  cross[[1097]]$father <- "Ra_Ra"; cross[[1097]]$mother <- "hA_ha"
  cross[[1098]]$father <- "Ra_Ra"; cross[[1098]]$mother <- "hA_RA"
  cross[[1099]]$father <- "Ra_Ra"; cross[[1099]]$mother <- "hA_Ra"
  cross[[1100]]$father <- "Ra_Ra"; cross[[1100]]$mother <- "hA_BA"
  cross[[1101]]$father <- "Ra_Ra"; cross[[1101]]$mother <- "hA_Ba"
  cross[[1102]]$father <- "Ra_Ra"; cross[[1102]]$mother <- "ha_ha"
  cross[[1103]]$father <- "Ra_Ra"; cross[[1103]]$mother <- "ha_RA"
  cross[[1104]]$father <- "Ra_Ra"; cross[[1104]]$mother <- "ha_Ra"
  cross[[1105]]$father <- "Ra_Ra"; cross[[1105]]$mother <- "ha_BA"
  cross[[1106]]$father <- "Ra_Ra"; cross[[1106]]$mother <- "ha_Ba"
  cross[[1107]]$father <- "Ra_Ra"; cross[[1107]]$mother <- "RA_RA"
  cross[[1108]]$father <- "Ra_Ra"; cross[[1108]]$mother <- "RA_Ra"
  cross[[1109]]$father <- "Ra_Ra"; cross[[1109]]$mother <- "RA_BA"
  cross[[1110]]$father <- "Ra_Ra"; cross[[1110]]$mother <- "RA_Ba"
  cross[[1111]]$father <- "Ra_Ra"; cross[[1111]]$mother <- "Ra_Ra"
  cross[[1112]]$father <- "Ra_Ra"; cross[[1112]]$mother <- "Ra_BA"
  cross[[1113]]$father <- "Ra_Ra"; cross[[1113]]$mother <- "Ra_Ba"
  cross[[1114]]$father <- "Ra_Ra"; cross[[1114]]$mother <- "BA_BA"
  cross[[1115]]$father <- "Ra_Ra"; cross[[1115]]$mother <- "BA_Ba"
  cross[[1116]]$father <- "Ra_Ra"; cross[[1116]]$mother <- "Ba_Ba"

  cross[[1117]]$father <- "Ra_BA"; cross[[1117]]$mother <- "HA_HA"
  cross[[1118]]$father <- "Ra_BA"; cross[[1118]]$mother <- "HA_Ha"
  cross[[1119]]$father <- "Ra_BA"; cross[[1119]]$mother <- "HA_hA"
  cross[[1120]]$father <- "Ra_BA"; cross[[1120]]$mother <- "HA_ha"
  cross[[1121]]$father <- "Ra_BA"; cross[[1121]]$mother <- "HA_RA"
  cross[[1122]]$father <- "Ra_BA"; cross[[1122]]$mother <- "HA_Ra"
  cross[[1123]]$father <- "Ra_BA"; cross[[1123]]$mother <- "HA_BA"
  cross[[1124]]$father <- "Ra_BA"; cross[[1124]]$mother <- "HA_Ba"
  cross[[1125]]$father <- "Ra_BA"; cross[[1125]]$mother <- "Ha_Ha"
  cross[[1126]]$father <- "Ra_BA"; cross[[1126]]$mother <- "Ha_hA"
  cross[[1127]]$father <- "Ra_BA"; cross[[1127]]$mother <- "Ha_ha"
  cross[[1128]]$father <- "Ra_BA"; cross[[1128]]$mother <- "Ha_RA"
  cross[[1129]]$father <- "Ra_BA"; cross[[1129]]$mother <- "Ha_Ra"
  cross[[1130]]$father <- "Ra_BA"; cross[[1130]]$mother <- "Ha_BA"
  cross[[1131]]$father <- "Ra_BA"; cross[[1131]]$mother <- "Ha_Ba"
  cross[[1132]]$father <- "Ra_BA"; cross[[1132]]$mother <- "hA_hA"
  cross[[1133]]$father <- "Ra_BA"; cross[[1133]]$mother <- "hA_ha"
  cross[[1134]]$father <- "Ra_BA"; cross[[1134]]$mother <- "hA_RA"
  cross[[1135]]$father <- "Ra_BA"; cross[[1135]]$mother <- "hA_Ra"
  cross[[1136]]$father <- "Ra_BA"; cross[[1136]]$mother <- "hA_BA"
  cross[[1137]]$father <- "Ra_BA"; cross[[1137]]$mother <- "hA_Ba"
  cross[[1138]]$father <- "Ra_BA"; cross[[1138]]$mother <- "ha_ha"
  cross[[1139]]$father <- "Ra_BA"; cross[[1139]]$mother <- "ha_RA"
  cross[[1140]]$father <- "Ra_BA"; cross[[1140]]$mother <- "ha_Ra"
  cross[[1141]]$father <- "Ra_BA"; cross[[1141]]$mother <- "ha_BA"
  cross[[1142]]$father <- "Ra_BA"; cross[[1142]]$mother <- "ha_Ba"
  cross[[1143]]$father <- "Ra_BA"; cross[[1143]]$mother <- "RA_RA"
  cross[[1144]]$father <- "Ra_BA"; cross[[1144]]$mother <- "RA_Ra"
  cross[[1145]]$father <- "Ra_BA"; cross[[1145]]$mother <- "RA_BA"
  cross[[1146]]$father <- "Ra_BA"; cross[[1146]]$mother <- "RA_Ba"
  cross[[1147]]$father <- "Ra_BA"; cross[[1147]]$mother <- "Ra_Ra"
  cross[[1148]]$father <- "Ra_BA"; cross[[1148]]$mother <- "Ra_BA"
  cross[[1149]]$father <- "Ra_BA"; cross[[1149]]$mother <- "Ra_Ba"
  cross[[1150]]$father <- "Ra_BA"; cross[[1150]]$mother <- "BA_BA"
  cross[[1151]]$father <- "Ra_BA"; cross[[1151]]$mother <- "BA_Ba"
  cross[[1152]]$father <- "Ra_BA"; cross[[1152]]$mother <- "Ba_Ba"

  cross[[1153]]$father <- "Ra_Ba"; cross[[1153]]$mother <- "HA_HA"
  cross[[1154]]$father <- "Ra_Ba"; cross[[1154]]$mother <- "HA_Ha"
  cross[[1155]]$father <- "Ra_Ba"; cross[[1155]]$mother <- "HA_hA"
  cross[[1156]]$father <- "Ra_Ba"; cross[[1156]]$mother <- "HA_ha"
  cross[[1157]]$father <- "Ra_Ba"; cross[[1157]]$mother <- "HA_RA"
  cross[[1158]]$father <- "Ra_Ba"; cross[[1158]]$mother <- "HA_Ra"
  cross[[1159]]$father <- "Ra_Ba"; cross[[1159]]$mother <- "HA_BA"
  cross[[1160]]$father <- "Ra_Ba"; cross[[1160]]$mother <- "HA_Ba"
  cross[[1161]]$father <- "Ra_Ba"; cross[[1161]]$mother <- "Ha_Ha"
  cross[[1162]]$father <- "Ra_Ba"; cross[[1162]]$mother <- "Ha_hA"
  cross[[1163]]$father <- "Ra_Ba"; cross[[1163]]$mother <- "Ha_ha"
  cross[[1164]]$father <- "Ra_Ba"; cross[[1164]]$mother <- "Ha_RA"
  cross[[1165]]$father <- "Ra_Ba"; cross[[1165]]$mother <- "Ha_Ra"
  cross[[1166]]$father <- "Ra_Ba"; cross[[1166]]$mother <- "Ha_BA"
  cross[[1167]]$father <- "Ra_Ba"; cross[[1167]]$mother <- "Ha_Ba"
  cross[[1168]]$father <- "Ra_Ba"; cross[[1168]]$mother <- "hA_hA"
  cross[[1169]]$father <- "Ra_Ba"; cross[[1169]]$mother <- "hA_ha"
  cross[[1170]]$father <- "Ra_Ba"; cross[[1170]]$mother <- "hA_RA"
  cross[[1171]]$father <- "Ra_Ba"; cross[[1171]]$mother <- "hA_Ra"
  cross[[1172]]$father <- "Ra_Ba"; cross[[1172]]$mother <- "hA_BA"
  cross[[1173]]$father <- "Ra_Ba"; cross[[1173]]$mother <- "hA_Ba"
  cross[[1174]]$father <- "Ra_Ba"; cross[[1174]]$mother <- "ha_ha"
  cross[[1175]]$father <- "Ra_Ba"; cross[[1175]]$mother <- "ha_RA"
  cross[[1176]]$father <- "Ra_Ba"; cross[[1176]]$mother <- "ha_Ra"
  cross[[1177]]$father <- "Ra_Ba"; cross[[1177]]$mother <- "ha_BA"
  cross[[1178]]$father <- "Ra_Ba"; cross[[1178]]$mother <- "ha_Ba"
  cross[[1179]]$father <- "Ra_Ba"; cross[[1179]]$mother <- "RA_RA"
  cross[[1180]]$father <- "Ra_Ba"; cross[[1180]]$mother <- "RA_Ra"
  cross[[1181]]$father <- "Ra_Ba"; cross[[1181]]$mother <- "RA_BA"
  cross[[1182]]$father <- "Ra_Ba"; cross[[1182]]$mother <- "RA_Ba"
  cross[[1183]]$father <- "Ra_Ba"; cross[[1183]]$mother <- "Ra_Ra"
  cross[[1184]]$father <- "Ra_Ba"; cross[[1184]]$mother <- "Ra_BA"
  cross[[1185]]$father <- "Ra_Ba"; cross[[1185]]$mother <- "Ra_Ba"
  cross[[1186]]$father <- "Ra_Ba"; cross[[1186]]$mother <- "BA_BA"
  cross[[1187]]$father <- "Ra_Ba"; cross[[1187]]$mother <- "BA_Ba"
  cross[[1188]]$father <- "Ra_Ba"; cross[[1188]]$mother <- "Ba_Ba"

  cross[[1189]]$father <- "BA_BA"; cross[[1189]]$mother <- "HA_HA"
  cross[[1190]]$father <- "BA_BA"; cross[[1190]]$mother <- "HA_Ha"
  cross[[1191]]$father <- "BA_BA"; cross[[1191]]$mother <- "HA_hA"
  cross[[1192]]$father <- "BA_BA"; cross[[1192]]$mother <- "HA_ha"
  cross[[1193]]$father <- "BA_BA"; cross[[1193]]$mother <- "HA_RA"
  cross[[1194]]$father <- "BA_BA"; cross[[1194]]$mother <- "HA_Ra"
  cross[[1195]]$father <- "BA_BA"; cross[[1195]]$mother <- "HA_BA"
  cross[[1196]]$father <- "BA_BA"; cross[[1196]]$mother <- "HA_Ba"
  cross[[1197]]$father <- "BA_BA"; cross[[1197]]$mother <- "Ha_Ha"
  cross[[1198]]$father <- "BA_BA"; cross[[1198]]$mother <- "Ha_hA"
  cross[[1199]]$father <- "BA_BA"; cross[[1199]]$mother <- "Ha_ha"
  cross[[1200]]$father <- "BA_BA"; cross[[1200]]$mother <- "Ha_RA"
  cross[[1201]]$father <- "BA_BA"; cross[[1201]]$mother <- "Ha_Ra"
  cross[[1202]]$father <- "BA_BA"; cross[[1202]]$mother <- "Ha_BA"
  cross[[1203]]$father <- "BA_BA"; cross[[1203]]$mother <- "Ha_Ba"
  cross[[1204]]$father <- "BA_BA"; cross[[1204]]$mother <- "hA_hA"
  cross[[1205]]$father <- "BA_BA"; cross[[1205]]$mother <- "hA_ha"
  cross[[1206]]$father <- "BA_BA"; cross[[1206]]$mother <- "hA_RA"
  cross[[1207]]$father <- "BA_BA"; cross[[1207]]$mother <- "hA_Ra"
  cross[[1208]]$father <- "BA_BA"; cross[[1208]]$mother <- "hA_BA"
  cross[[1209]]$father <- "BA_BA"; cross[[1209]]$mother <- "hA_Ba"
  cross[[1210]]$father <- "BA_BA"; cross[[1210]]$mother <- "ha_ha"
  cross[[1211]]$father <- "BA_BA"; cross[[1211]]$mother <- "ha_RA"
  cross[[1212]]$father <- "BA_BA"; cross[[1212]]$mother <- "ha_Ra"
  cross[[1213]]$father <- "BA_BA"; cross[[1213]]$mother <- "ha_BA"
  cross[[1214]]$father <- "BA_BA"; cross[[1214]]$mother <- "ha_Ba"
  cross[[1215]]$father <- "BA_BA"; cross[[1215]]$mother <- "RA_RA"
  cross[[1216]]$father <- "BA_BA"; cross[[1216]]$mother <- "RA_Ra"
  cross[[1217]]$father <- "BA_BA"; cross[[1217]]$mother <- "RA_BA"
  cross[[1218]]$father <- "BA_BA"; cross[[1218]]$mother <- "RA_Ba"
  cross[[1219]]$father <- "BA_BA"; cross[[1219]]$mother <- "Ra_Ra"
  cross[[1220]]$father <- "BA_BA"; cross[[1220]]$mother <- "Ra_BA"
  cross[[1221]]$father <- "BA_BA"; cross[[1221]]$mother <- "Ra_Ba"
  cross[[1222]]$father <- "BA_BA"; cross[[1222]]$mother <- "BA_BA"
  cross[[1223]]$father <- "BA_BA"; cross[[1223]]$mother <- "BA_Ba"
  cross[[1224]]$father <- "BA_BA"; cross[[1224]]$mother <- "Ba_Ba"

  cross[[1225]]$father <- "BA_Ba"; cross[[1225]]$mother <- "HA_HA"
  cross[[1226]]$father <- "BA_Ba"; cross[[1226]]$mother <- "HA_Ha"
  cross[[1227]]$father <- "BA_Ba"; cross[[1227]]$mother <- "HA_hA"
  cross[[1228]]$father <- "BA_Ba"; cross[[1228]]$mother <- "HA_ha"
  cross[[1229]]$father <- "BA_Ba"; cross[[1229]]$mother <- "HA_RA"
  cross[[1230]]$father <- "BA_Ba"; cross[[1230]]$mother <- "HA_Ra"
  cross[[1231]]$father <- "BA_Ba"; cross[[1231]]$mother <- "HA_BA"
  cross[[1232]]$father <- "BA_Ba"; cross[[1232]]$mother <- "HA_Ba"
  cross[[1233]]$father <- "BA_Ba"; cross[[1233]]$mother <- "Ha_Ha"
  cross[[1234]]$father <- "BA_Ba"; cross[[1234]]$mother <- "Ha_hA"
  cross[[1235]]$father <- "BA_Ba"; cross[[1235]]$mother <- "Ha_ha"
  cross[[1236]]$father <- "BA_Ba"; cross[[1236]]$mother <- "Ha_RA"
  cross[[1237]]$father <- "BA_Ba"; cross[[1237]]$mother <- "Ha_Ra"
  cross[[1238]]$father <- "BA_Ba"; cross[[1238]]$mother <- "Ha_BA"
  cross[[1239]]$father <- "BA_Ba"; cross[[1239]]$mother <- "Ha_Ba"
  cross[[1240]]$father <- "BA_Ba"; cross[[1240]]$mother <- "hA_hA"
  cross[[1241]]$father <- "BA_Ba"; cross[[1241]]$mother <- "hA_ha"
  cross[[1242]]$father <- "BA_Ba"; cross[[1242]]$mother <- "hA_RA"
  cross[[1243]]$father <- "BA_Ba"; cross[[1243]]$mother <- "hA_Ra"
  cross[[1244]]$father <- "BA_Ba"; cross[[1244]]$mother <- "hA_BA"
  cross[[1245]]$father <- "BA_Ba"; cross[[1245]]$mother <- "hA_Ba"
  cross[[1246]]$father <- "BA_Ba"; cross[[1246]]$mother <- "ha_ha"
  cross[[1247]]$father <- "BA_Ba"; cross[[1247]]$mother <- "ha_RA"
  cross[[1248]]$father <- "BA_Ba"; cross[[1248]]$mother <- "ha_Ra"
  cross[[1249]]$father <- "BA_Ba"; cross[[1249]]$mother <- "ha_BA"
  cross[[1250]]$father <- "BA_Ba"; cross[[1250]]$mother <- "ha_Ba"
  cross[[1251]]$father <- "BA_Ba"; cross[[1251]]$mother <- "RA_RA"
  cross[[1252]]$father <- "BA_Ba"; cross[[1252]]$mother <- "RA_Ra"
  cross[[1253]]$father <- "BA_Ba"; cross[[1253]]$mother <- "RA_BA"
  cross[[1254]]$father <- "BA_Ba"; cross[[1254]]$mother <- "RA_Ba"
  cross[[1255]]$father <- "BA_Ba"; cross[[1255]]$mother <- "Ra_Ra"
  cross[[1256]]$father <- "BA_Ba"; cross[[1256]]$mother <- "Ra_BA"
  cross[[1257]]$father <- "BA_Ba"; cross[[1257]]$mother <- "Ra_Ba"
  cross[[1258]]$father <- "BA_Ba"; cross[[1258]]$mother <- "BA_BA"
  cross[[1259]]$father <- "BA_Ba"; cross[[1259]]$mother <- "BA_Ba"
  cross[[1260]]$father <- "BA_Ba"; cross[[1260]]$mother <- "Ba_Ba"

  cross[[1261]]$father <- "Ba_Ba"; cross[[1261]]$mother <- "HA_HA"
  cross[[1262]]$father <- "Ba_Ba"; cross[[1262]]$mother <- "HA_Ha"
  cross[[1263]]$father <- "Ba_Ba"; cross[[1263]]$mother <- "HA_hA"
  cross[[1264]]$father <- "Ba_Ba"; cross[[1264]]$mother <- "HA_ha"
  cross[[1265]]$father <- "Ba_Ba"; cross[[1265]]$mother <- "HA_RA"
  cross[[1266]]$father <- "Ba_Ba"; cross[[1266]]$mother <- "HA_Ra"
  cross[[1267]]$father <- "Ba_Ba"; cross[[1267]]$mother <- "HA_BA"
  cross[[1268]]$father <- "Ba_Ba"; cross[[1268]]$mother <- "HA_Ba"
  cross[[1269]]$father <- "Ba_Ba"; cross[[1269]]$mother <- "Ha_Ha"
  cross[[1270]]$father <- "Ba_Ba"; cross[[1270]]$mother <- "Ha_hA"
  cross[[1271]]$father <- "Ba_Ba"; cross[[1271]]$mother <- "Ha_ha"
  cross[[1272]]$father <- "Ba_Ba"; cross[[1272]]$mother <- "Ha_RA"
  cross[[1273]]$father <- "Ba_Ba"; cross[[1273]]$mother <- "Ha_Ra"
  cross[[1274]]$father <- "Ba_Ba"; cross[[1274]]$mother <- "Ha_BA"
  cross[[1275]]$father <- "Ba_Ba"; cross[[1275]]$mother <- "Ha_Ba"
  cross[[1276]]$father <- "Ba_Ba"; cross[[1276]]$mother <- "hA_hA"
  cross[[1277]]$father <- "Ba_Ba"; cross[[1277]]$mother <- "hA_ha"
  cross[[1278]]$father <- "Ba_Ba"; cross[[1278]]$mother <- "hA_RA"
  cross[[1279]]$father <- "Ba_Ba"; cross[[1279]]$mother <- "hA_Ra"
  cross[[1280]]$father <- "Ba_Ba"; cross[[1280]]$mother <- "hA_BA"
  cross[[1281]]$father <- "Ba_Ba"; cross[[1281]]$mother <- "hA_Ba"
  cross[[1282]]$father <- "Ba_Ba"; cross[[1282]]$mother <- "ha_ha"
  cross[[1283]]$father <- "Ba_Ba"; cross[[1283]]$mother <- "ha_RA"
  cross[[1284]]$father <- "Ba_Ba"; cross[[1284]]$mother <- "ha_Ra"
  cross[[1285]]$father <- "Ba_Ba"; cross[[1285]]$mother <- "ha_BA"
  cross[[1286]]$father <- "Ba_Ba"; cross[[1286]]$mother <- "ha_Ba"
  cross[[1287]]$father <- "Ba_Ba"; cross[[1287]]$mother <- "RA_RA"
  cross[[1288]]$father <- "Ba_Ba"; cross[[1288]]$mother <- "RA_Ra"
  cross[[1289]]$father <- "Ba_Ba"; cross[[1289]]$mother <- "RA_BA"
  cross[[1290]]$father <- "Ba_Ba"; cross[[1290]]$mother <- "RA_Ba"
  cross[[1291]]$father <- "Ba_Ba"; cross[[1291]]$mother <- "Ra_Ra"
  cross[[1292]]$father <- "Ba_Ba"; cross[[1292]]$mother <- "Ra_BA"
  cross[[1293]]$father <- "Ba_Ba"; cross[[1293]]$mother <- "Ra_Ba"
  cross[[1294]]$father <- "Ba_Ba"; cross[[1294]]$mother <- "BA_BA"
  cross[[1295]]$father <- "Ba_Ba"; cross[[1295]]$mother <- "BA_Ba"
  cross[[1296]]$father <- "Ba_Ba"; cross[[1296]]$mother <- "Ba_Ba"

  ## Initialize list entries:
  for (i in 1:numCrosses) {

    ## Initialize parental gamete frequencies:
    cross[[i]]$maternal_HA <- 0
    cross[[i]]$maternal_Ha <- 0
    cross[[i]]$maternal_hA <- 0
    cross[[i]]$maternal_ha <- 1
    cross[[i]]$maternal_RA <- 0
    cross[[i]]$maternal_Ra <- 0
    cross[[i]]$maternal_BA <- 0
    cross[[i]]$maternal_Ba <- 0

    cross[[i]]$paternal_HA <- 0
    cross[[i]]$paternal_Ha <- 0
    cross[[i]]$paternal_hA <- 0
    cross[[i]]$paternal_ha <- 1
    cross[[i]]$paternal_RA <- 0
    cross[[i]]$paternal_Ra <- 0
    cross[[i]]$paternal_BA <- 0
    cross[[i]]$paternal_Ba <- 0

    ## Initialize offspring genotype frequencies:
    cross[[i]]$HA_HA_m <- 0
    cross[[i]]$HA_Ha_m <- 0
    cross[[i]]$HA_hA_m <- 0
    cross[[i]]$HA_ha_m <- 0
    cross[[i]]$HA_RA_m <- 0
    cross[[i]]$HA_Ra_m <- 0
    cross[[i]]$HA_BA_m <- 0
    cross[[i]]$HA_Ba_m <- 0
    cross[[i]]$Ha_Ha_m <- 0
    cross[[i]]$Ha_hA_m <- 0
    cross[[i]]$Ha_ha_m <- 0
    cross[[i]]$Ha_RA_m <- 0
    cross[[i]]$Ha_Ra_m <- 0
    cross[[i]]$Ha_BA_m <- 0
    cross[[i]]$Ha_Ba_m <- 0
    cross[[i]]$hA_hA_m <- 0
    cross[[i]]$hA_ha_m <- 0
    cross[[i]]$hA_RA_m <- 0
    cross[[i]]$hA_Ra_m <- 0
    cross[[i]]$hA_BA_m <- 0
    cross[[i]]$hA_Ba_m <- 0
    cross[[i]]$ha_ha_m <- 1
    cross[[i]]$ha_RA_m <- 0
    cross[[i]]$ha_Ra_m <- 0
    cross[[i]]$ha_BA_m <- 0
    cross[[i]]$ha_Ba_m <- 0
    cross[[i]]$RA_RA_m <- 0
    cross[[i]]$RA_Ra_m <- 0
    cross[[i]]$RA_BA_m <- 0
    cross[[i]]$RA_Ba_m <- 0
    cross[[i]]$Ra_Ra_m <- 0
    cross[[i]]$Ra_BA_m <- 0
    cross[[i]]$Ra_Ba_m <- 0
    cross[[i]]$BA_BA_m <- 0
    cross[[i]]$BA_Ba_m <- 0
    cross[[i]]$Ba_Ba_m <- 0

    cross[[i]]$HA_HA_f <- 0
    cross[[i]]$HA_Ha_f <- 0
    cross[[i]]$HA_hA_f <- 0
    cross[[i]]$HA_ha_f <- 0
    cross[[i]]$HA_RA_f <- 0
    cross[[i]]$HA_Ra_f <- 0
    cross[[i]]$HA_BA_f <- 0
    cross[[i]]$HA_Ba_f <- 0
    cross[[i]]$Ha_Ha_f <- 0
    cross[[i]]$Ha_hA_f <- 0
    cross[[i]]$Ha_ha_f <- 0
    cross[[i]]$Ha_RA_f <- 0
    cross[[i]]$Ha_Ra_f <- 0
    cross[[i]]$Ha_BA_f <- 0
    cross[[i]]$Ha_Ba_f <- 0
    cross[[i]]$hA_hA_f <- 0
    cross[[i]]$hA_ha_f <- 0
    cross[[i]]$hA_RA_f <- 0
    cross[[i]]$hA_Ra_f <- 0
    cross[[i]]$hA_BA_f <- 0
    cross[[i]]$hA_Ba_f <- 0
    cross[[i]]$ha_ha_f <- 1
    cross[[i]]$ha_RA_f <- 0
    cross[[i]]$ha_Ra_f <- 0
    cross[[i]]$ha_BA_f <- 0
    cross[[i]]$ha_Ba_f <- 0
    cross[[i]]$RA_RA_f <- 0
    cross[[i]]$RA_Ra_f <- 0
    cross[[i]]$RA_BA_f <- 0
    cross[[i]]$RA_Ba_f <- 0
    cross[[i]]$Ra_Ra_f <- 0
    cross[[i]]$Ra_BA_f <- 0
    cross[[i]]$Ra_Ba_f <- 0
    cross[[i]]$BA_BA_f <- 0
    cross[[i]]$BA_Ba_f <- 0
    cross[[i]]$Ba_Ba_f <- 0
  }

  ## Calculate parental gamete frequencies according to homing rules:
  for (i in 1:numCrosses) {

    ## Maternal gamete frequencies:
    if (cross[[i]]$mother == "HA_HA") {
      cross[[i]]$maternal_HA <- 1
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- 0
      cross[[i]]$maternal_Ra <- 0
      cross[[i]]$maternal_BA <- 0
      cross[[i]]$maternal_Ba <- 0
    }
    if (cross[[i]]$mother == "HA_Ha") {
      cross[[i]]$maternal_HA <- 1/2
      cross[[i]]$maternal_Ha <- 1/2
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- 0
      cross[[i]]$maternal_Ra <- 0
      cross[[i]]$maternal_BA <- 0
      cross[[i]]$maternal_Ba <- 0
    }
    if (cross[[i]]$mother == "HA_hA") {
      cross[[i]]$maternal_HA <- (1+eF)/2
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- (1-eF-rhoRF-rhoBF)/2
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- rhoRF/2
      cross[[i]]$maternal_Ra <- 0
      cross[[i]]$maternal_BA <- rhoBF/2
      cross[[i]]$maternal_Ba <- 0
    }
    if (cross[[i]]$mother == "HA_ha") {
      cross[[i]]$maternal_HA <- r*(eF/2) + (1-r)*(1/2)
      cross[[i]]$maternal_Ha <- r*(1/2) + (1-r)*(eF/2)
      cross[[i]]$maternal_hA <- r*((1-eF-rhoRF-rhoBF)/2)
      cross[[i]]$maternal_ha <- (1-r)*((1-eF-rhoRF-rhoBF)/2)
      cross[[i]]$maternal_RA <- r*(rhoRF/2)
      cross[[i]]$maternal_Ra <- (1-r)*(rhoRF/2)
      cross[[i]]$maternal_BA <- r*(rhoBF/2)
      cross[[i]]$maternal_Ba <- (1-r)*(rhoBF/2)
    }
    if (cross[[i]]$mother == "HA_RA") {
      cross[[i]]$maternal_HA <- 1/2
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- 1/2
      cross[[i]]$maternal_Ra <- 0
      cross[[i]]$maternal_BA <- 0
      cross[[i]]$maternal_Ba <- 0
    }
    if (cross[[i]]$mother == "HA_Ra") {
      cross[[i]]$maternal_HA <- (1-r)*(1/2)
      cross[[i]]$maternal_Ha <- r*(1/2)
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- r*(1/2)
      cross[[i]]$maternal_Ra <- (1-r)*(1/2)
      cross[[i]]$maternal_BA <- 0
      cross[[i]]$maternal_Ba <- 0
    }
    if (cross[[i]]$mother == "HA_BA") {
      cross[[i]]$maternal_HA <- 1/2
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- 0
      cross[[i]]$maternal_Ra <- 0
      cross[[i]]$maternal_BA <- 1/2
      cross[[i]]$maternal_Ba <- 0
    }
    if (cross[[i]]$mother == "HA_Ba") {
      cross[[i]]$maternal_HA <- (1-r)*(1/2)
      cross[[i]]$maternal_Ha <- r*(1/2)
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- 0
      cross[[i]]$maternal_Ra <- 0
      cross[[i]]$maternal_BA <- r*(1/2)
      cross[[i]]$maternal_Ba <- (1-r)*(1/2)
    }
    if (cross[[i]]$mother == "Ha_Ha") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 1
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- 0
      cross[[i]]$maternal_Ra <- 0
      cross[[i]]$maternal_BA <- 0
      cross[[i]]$maternal_Ba <- 0
    }
    if (cross[[i]]$mother == "Ha_hA") {
      cross[[i]]$maternal_HA <- (1-r)*(eF/2) + r*(1/2)
      cross[[i]]$maternal_Ha <- (1-r)*(1/2) + r*(eF/2)
      cross[[i]]$maternal_hA <- (1-r)*((1-eF-rhoRF-rhoBF)/2)
      cross[[i]]$maternal_ha <- r*((1-eF-rhoRF-rhoBF)/2)
      cross[[i]]$maternal_RA <- (1-r)*(rhoRF/2)
      cross[[i]]$maternal_Ra <- r*(rhoRF/2)
      cross[[i]]$maternal_BA <- (1-r)*(rhoBF/2)
      cross[[i]]$maternal_Ba <- r*(rhoBF/2)
    }
    if (cross[[i]]$mother == "Ha_ha") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- (1+eF)/2
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- (1-eF-rhoRF-rhoBF)/2
      cross[[i]]$maternal_RA <- 0
      cross[[i]]$maternal_Ra <- rhoRF/2
      cross[[i]]$maternal_BA <- 0
      cross[[i]]$maternal_Ba <- rhoBF/2
    }
    if (cross[[i]]$mother == "Ha_RA") {
      cross[[i]]$maternal_HA <- r*(1/2)
      cross[[i]]$maternal_Ha <- (1-r)*(1/2)
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- (1-r)*(1/2)
      cross[[i]]$maternal_Ra <- r*(1/2)
      cross[[i]]$maternal_BA <- 0
      cross[[i]]$maternal_Ba <- 0
    }
    if (cross[[i]]$mother == "Ha_Ra") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 1/2
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- 0
      cross[[i]]$maternal_Ra <- 1/2
      cross[[i]]$maternal_BA <- 0
      cross[[i]]$maternal_Ba <- 0
    }
    if (cross[[i]]$mother == "Ha_BA") {
      cross[[i]]$maternal_HA <- r*(1/2)
      cross[[i]]$maternal_Ha <- (1-r)*(1/2)
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- 0
      cross[[i]]$maternal_Ra <- 0
      cross[[i]]$maternal_BA <- (1-r)*(1/2)
      cross[[i]]$maternal_Ba <- r*(1/2)
    }
    if (cross[[i]]$mother == "Ha_Ba") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 1/2
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- 0
      cross[[i]]$maternal_Ra <- 0
      cross[[i]]$maternal_BA <- 0
      cross[[i]]$maternal_Ba <- 1/2
    }
    if (cross[[i]]$mother == "hA_hA") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- 1
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- 0
      cross[[i]]$maternal_Ra <- 0
      cross[[i]]$maternal_BA <- 0
      cross[[i]]$maternal_Ba <- 0
    }
    if (cross[[i]]$mother == "hA_ha") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- 1/2
      cross[[i]]$maternal_ha <- 1/2
      cross[[i]]$maternal_RA <- 0
      cross[[i]]$maternal_Ra <- 0
      cross[[i]]$maternal_BA <- 0
      cross[[i]]$maternal_Ba <- 0
    }
    if (cross[[i]]$mother == "hA_RA") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- 1/2
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- 1/2
      cross[[i]]$maternal_Ra <- 0
      cross[[i]]$maternal_BA <- 0
      cross[[i]]$maternal_Ba <- 0
    }
    if (cross[[i]]$mother == "hA_Ra") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- (1-r)*(1/2)
      cross[[i]]$maternal_ha <- r*(1/2)
      cross[[i]]$maternal_RA <- r*(1/2)
      cross[[i]]$maternal_Ra <- (1-r)*(1/2)
      cross[[i]]$maternal_BA <- 0
      cross[[i]]$maternal_Ba <- 0
    }
    if (cross[[i]]$mother == "hA_BA") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- 1/2
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- 0
      cross[[i]]$maternal_Ra <- 0
      cross[[i]]$maternal_BA <- 1/2
      cross[[i]]$maternal_Ba <- 0
    }
    if (cross[[i]]$mother == "hA_Ba") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- (1-r)*(1/2)
      cross[[i]]$maternal_ha <- r*(1/2)
      cross[[i]]$maternal_RA <- 0
      cross[[i]]$maternal_Ra <- 0
      cross[[i]]$maternal_BA <- r*(1/2)
      cross[[i]]$maternal_Ba <- (1-r)*(1/2)
    }
    if (cross[[i]]$mother == "ha_ha") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 1
      cross[[i]]$maternal_RA <- 0
      cross[[i]]$maternal_Ra <- 0
      cross[[i]]$maternal_BA <- 0
      cross[[i]]$maternal_Ba <- 0
    }
    if (cross[[i]]$mother == "ha_RA") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- r*(1/2)
      cross[[i]]$maternal_ha <- (1-r)*(1/2)
      cross[[i]]$maternal_RA <- (1-r)*(1/2)
      cross[[i]]$maternal_Ra <- r*(1/2)
      cross[[i]]$maternal_BA <- 0
      cross[[i]]$maternal_Ba <- 0
    }
    if (cross[[i]]$mother == "ha_Ra") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 1/2
      cross[[i]]$maternal_RA <- 0
      cross[[i]]$maternal_Ra <- 1/2
      cross[[i]]$maternal_BA <- 0
      cross[[i]]$maternal_Ba <- 0
    }
    if (cross[[i]]$mother == "ha_BA") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- r*(1/2)
      cross[[i]]$maternal_ha <- (1-r)*(1/2)
      cross[[i]]$maternal_RA <- 0
      cross[[i]]$maternal_Ra <- 0
      cross[[i]]$maternal_BA <- (1-r)*(1/2)
      cross[[i]]$maternal_Ba <- r*(1/2)
    }
    if (cross[[i]]$mother == "ha_Ba") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 1/2
      cross[[i]]$maternal_RA <- 0
      cross[[i]]$maternal_Ra <- 0
      cross[[i]]$maternal_BA <- 0
      cross[[i]]$maternal_Ba <- 1/2
    }
    if (cross[[i]]$mother == "RA_RA") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- 1
      cross[[i]]$maternal_Ra <- 0
      cross[[i]]$maternal_BA <- 0
      cross[[i]]$maternal_Ba <- 0
    }
    if (cross[[i]]$mother == "RA_Ra") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- 1/2
      cross[[i]]$maternal_Ra <- 1/2
      cross[[i]]$maternal_BA <- 0
      cross[[i]]$maternal_Ba <- 0
    }
    if (cross[[i]]$mother == "RA_BA") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- 1/2
      cross[[i]]$maternal_Ra <- 0
      cross[[i]]$maternal_BA <- 1/2
      cross[[i]]$maternal_Ba <- 0
    }
    if (cross[[i]]$mother == "RA_Ba") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- (1-r)*(1/2)
      cross[[i]]$maternal_Ra <- r*(1/2)
      cross[[i]]$maternal_BA <- r*(1/2)
      cross[[i]]$maternal_Ba <- (1-r)*(1/2)
    }
    if (cross[[i]]$mother == "Ra_Ra") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- 0
      cross[[i]]$maternal_Ra <- 1
      cross[[i]]$maternal_BA <- 0
      cross[[i]]$maternal_Ba <- 0
    }
    if (cross[[i]]$mother == "Ra_BA") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- r*(1/2)
      cross[[i]]$maternal_Ra <- (1-r)*(1/2)
      cross[[i]]$maternal_BA <- (1-r)*(1/2)
      cross[[i]]$maternal_Ba <- r*(1/2)
    }
    if (cross[[i]]$mother == "Ra_Ba") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- 0
      cross[[i]]$maternal_Ra <- 1/2
      cross[[i]]$maternal_BA <- 0
      cross[[i]]$maternal_Ba <- 1/2
    }
    if (cross[[i]]$mother == "BA_BA") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- 0
      cross[[i]]$maternal_Ra <- 0
      cross[[i]]$maternal_BA <- 1
      cross[[i]]$maternal_Ba <- 0
    }
    if (cross[[i]]$mother == "BA_Ba") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- 0
      cross[[i]]$maternal_Ra <- 0
      cross[[i]]$maternal_BA <- 1/2
      cross[[i]]$maternal_Ba <- 1/2
    }
    if (cross[[i]]$mother == "Ba_Ba") {
      cross[[i]]$maternal_HA <- 0
      cross[[i]]$maternal_Ha <- 0
      cross[[i]]$maternal_hA <- 0
      cross[[i]]$maternal_ha <- 0
      cross[[i]]$maternal_RA <- 0
      cross[[i]]$maternal_Ra <- 0
      cross[[i]]$maternal_BA <- 0
      cross[[i]]$maternal_Ba <- 1
    }

    ## Paternal gamete frequencies:
    if (cross[[i]]$father == "HA_HA") {
      cross[[i]]$paternal_HA <- 1
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- 0
      cross[[i]]$paternal_Ra <- 0
      cross[[i]]$paternal_BA <- 0
      cross[[i]]$paternal_Ba <- 0
    }
    if (cross[[i]]$father == "HA_Ha") {
      cross[[i]]$paternal_HA <- 1/2
      cross[[i]]$paternal_Ha <- 1/2
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- 0
      cross[[i]]$paternal_Ra <- 0
      cross[[i]]$paternal_BA <- 0
      cross[[i]]$paternal_Ba <- 0
    }
    if (cross[[i]]$father == "HA_hA") {
      cross[[i]]$paternal_HA <- (1+eM)/2
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- (1-eM-rhoRM-rhoBM)/2
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- rhoRM/2
      cross[[i]]$paternal_Ra <- 0
      cross[[i]]$paternal_BA <- rhoBM/2
      cross[[i]]$paternal_Ba <- 0
    }
    if (cross[[i]]$father == "HA_ha") {
      cross[[i]]$paternal_HA <- r*(eM/2) + (1-r)*(1/2)
      cross[[i]]$paternal_Ha <- r*(1/2) + (1-r)*(eM/2)
      cross[[i]]$paternal_hA <- r*((1-eM-rhoRM-rhoBM)/2)
      cross[[i]]$paternal_ha <- (1-r)*((1-eM-rhoRM-rhoBM)/2)
      cross[[i]]$paternal_RA <- r*(rhoRM/2)
      cross[[i]]$paternal_Ra <- (1-r)*(rhoRM/2)
      cross[[i]]$paternal_BA <- r*(rhoBM/2)
      cross[[i]]$paternal_Ba <- (1-r)*(rhoBM/2)
    }
    if (cross[[i]]$father == "HA_RA") {
      cross[[i]]$paternal_HA <- 1/2
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- 1/2
      cross[[i]]$paternal_Ra <- 0
      cross[[i]]$paternal_BA <- 0
      cross[[i]]$paternal_Ba <- 0
    }
    if (cross[[i]]$father == "HA_Ra") {
      cross[[i]]$paternal_HA <- (1-r)*(1/2)
      cross[[i]]$paternal_Ha <- r*(1/2)
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- r*(1/2)
      cross[[i]]$paternal_Ra <- (1-r)*(1/2)
      cross[[i]]$paternal_BA <- 0
      cross[[i]]$paternal_Ba <- 0
    }
    if (cross[[i]]$father == "HA_BA") {
      cross[[i]]$paternal_HA <- 1/2
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- 0
      cross[[i]]$paternal_Ra <- 0
      cross[[i]]$paternal_BA <- 1/2
      cross[[i]]$paternal_Ba <- 0
    }
    if (cross[[i]]$father == "HA_Ba") {
      cross[[i]]$paternal_HA <- (1-r)*(1/2)
      cross[[i]]$paternal_Ha <- r*(1/2)
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- 0
      cross[[i]]$paternal_Ra <- 0
      cross[[i]]$paternal_BA <- r*(1/2)
      cross[[i]]$paternal_Ba <- (1-r)*(1/2)
    }
    if (cross[[i]]$father == "Ha_Ha") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 1
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- 0
      cross[[i]]$paternal_Ra <- 0
      cross[[i]]$paternal_BA <- 0
      cross[[i]]$paternal_Ba <- 0
    }
    if (cross[[i]]$father == "Ha_hA") {
      cross[[i]]$paternal_HA <- (1-r)*(eM/2) + r*(1/2)
      cross[[i]]$paternal_Ha <- (1-r)*(1/2) + r*(eM/2)
      cross[[i]]$paternal_hA <- (1-r)*((1-eM-rhoRM-rhoBM)/2)
      cross[[i]]$paternal_ha <- r*((1-eM-rhoRM-rhoBM)/2)
      cross[[i]]$paternal_RA <- (1-r)*(rhoRM/2)
      cross[[i]]$paternal_Ra <- r*(rhoRM/2)
      cross[[i]]$paternal_BA <- (1-r)*(rhoBM/2)
      cross[[i]]$paternal_Ba <- r*(rhoBM/2)
    }
    if (cross[[i]]$father == "Ha_ha") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- (1+eM)/2
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- (1-eM-rhoRM-rhoBM)/2
      cross[[i]]$paternal_RA <- 0
      cross[[i]]$paternal_Ra <- rhoRM/2
      cross[[i]]$paternal_BA <- 0
      cross[[i]]$paternal_Ba <- rhoBM/2
    }
    if (cross[[i]]$father == "Ha_RA") {
      cross[[i]]$paternal_HA <- r*(1/2)
      cross[[i]]$paternal_Ha <- (1-r)*(1/2)
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- (1-r)*(1/2)
      cross[[i]]$paternal_Ra <- r*(1/2)
      cross[[i]]$paternal_BA <- 0
      cross[[i]]$paternal_Ba <- 0
    }
    if (cross[[i]]$father == "Ha_Ra") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 1/2
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- 0
      cross[[i]]$paternal_Ra <- 1/2
      cross[[i]]$paternal_BA <- 0
      cross[[i]]$paternal_Ba <- 0
    }
    if (cross[[i]]$father == "Ha_BA") {
      cross[[i]]$paternal_HA <- r*(1/2)
      cross[[i]]$paternal_Ha <- (1-r)*(1/2)
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- 0
      cross[[i]]$paternal_Ra <- 0
      cross[[i]]$paternal_BA <- (1-r)*(1/2)
      cross[[i]]$paternal_Ba <- r*(1/2)
    }
    if (cross[[i]]$father == "Ha_Ba") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 1/2
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- 0
      cross[[i]]$paternal_Ra <- 0
      cross[[i]]$paternal_BA <- 0
      cross[[i]]$paternal_Ba <- 1/2
    }
    if (cross[[i]]$father == "hA_hA") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- 1
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- 0
      cross[[i]]$paternal_Ra <- 0
      cross[[i]]$paternal_BA <- 0
      cross[[i]]$paternal_Ba <- 0
    }
    if (cross[[i]]$father == "hA_ha") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- 1/2
      cross[[i]]$paternal_ha <- 1/2
      cross[[i]]$paternal_RA <- 0
      cross[[i]]$paternal_Ra <- 0
      cross[[i]]$paternal_BA <- 0
      cross[[i]]$paternal_Ba <- 0
    }
    if (cross[[i]]$father == "hA_RA") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- 1/2
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- 1/2
      cross[[i]]$paternal_Ra <- 0
      cross[[i]]$paternal_BA <- 0
      cross[[i]]$paternal_Ba <- 0
    }
    if (cross[[i]]$father == "hA_Ra") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- (1-r)*(1/2)
      cross[[i]]$paternal_ha <- r*(1/2)
      cross[[i]]$paternal_RA <- r*(1/2)
      cross[[i]]$paternal_Ra <- (1-r)*(1/2)
      cross[[i]]$paternal_BA <- 0
      cross[[i]]$paternal_Ba <- 0
    }
    if (cross[[i]]$father == "hA_BA") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- 1/2
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- 0
      cross[[i]]$paternal_Ra <- 0
      cross[[i]]$paternal_BA <- 1/2
      cross[[i]]$paternal_Ba <- 0
    }
    if (cross[[i]]$father == "hA_Ba") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- (1-r)*(1/2)
      cross[[i]]$paternal_ha <- r*(1/2)
      cross[[i]]$paternal_RA <- 0
      cross[[i]]$paternal_Ra <- 0
      cross[[i]]$paternal_BA <- r*(1/2)
      cross[[i]]$paternal_Ba <- (1-r)*(1/2)
    }
    if (cross[[i]]$father == "ha_ha") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 1
      cross[[i]]$paternal_RA <- 0
      cross[[i]]$paternal_Ra <- 0
      cross[[i]]$paternal_BA <- 0
      cross[[i]]$paternal_Ba <- 0
    }
    if (cross[[i]]$father == "ha_RA") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- r*(1/2)
      cross[[i]]$paternal_ha <- (1-r)*(1/2)
      cross[[i]]$paternal_RA <- (1-r)*(1/2)
      cross[[i]]$paternal_Ra <- r*(1/2)
      cross[[i]]$paternal_BA <- 0
      cross[[i]]$paternal_Ba <- 0
    }
    if (cross[[i]]$father == "ha_Ra") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 1/2
      cross[[i]]$paternal_RA <- 0
      cross[[i]]$paternal_Ra <- 1/2
      cross[[i]]$paternal_BA <- 0
      cross[[i]]$paternal_Ba <- 0
    }
    if (cross[[i]]$father == "ha_BA") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- r*(1/2)
      cross[[i]]$paternal_ha <- (1-r)*(1/2)
      cross[[i]]$paternal_RA <- 0
      cross[[i]]$paternal_Ra <- 0
      cross[[i]]$paternal_BA <- (1-r)*(1/2)
      cross[[i]]$paternal_Ba <- r*(1/2)
    }
    if (cross[[i]]$father == "ha_Ba") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 1/2
      cross[[i]]$paternal_RA <- 0
      cross[[i]]$paternal_Ra <- 0
      cross[[i]]$paternal_BA <- 0
      cross[[i]]$paternal_Ba <- 1/2
    }
    if (cross[[i]]$father == "RA_RA") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- 1
      cross[[i]]$paternal_Ra <- 0
      cross[[i]]$paternal_BA <- 0
      cross[[i]]$paternal_Ba <- 0
    }
    if (cross[[i]]$father == "RA_Ra") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- 1/2
      cross[[i]]$paternal_Ra <- 1/2
      cross[[i]]$paternal_BA <- 0
      cross[[i]]$paternal_Ba <- 0
    }
    if (cross[[i]]$father == "RA_BA") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- 1/2
      cross[[i]]$paternal_Ra <- 0
      cross[[i]]$paternal_BA <- 1/2
      cross[[i]]$paternal_Ba <- 0
    }
    if (cross[[i]]$father == "RA_Ba") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- (1-r)*(1/2)
      cross[[i]]$paternal_Ra <- r*(1/2)
      cross[[i]]$paternal_BA <- r*(1/2)
      cross[[i]]$paternal_Ba <- (1-r)*(1/2)
    }
    if (cross[[i]]$father == "Ra_Ra") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- 0
      cross[[i]]$paternal_Ra <- 1
      cross[[i]]$paternal_BA <- 0
      cross[[i]]$paternal_Ba <- 0
    }
    if (cross[[i]]$father == "Ra_BA") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- r*(1/2)
      cross[[i]]$paternal_Ra <- (1-r)*(1/2)
      cross[[i]]$paternal_BA <- (1-r)*(1/2)
      cross[[i]]$paternal_Ba <- r*(1/2)
    }
    if (cross[[i]]$father == "Ra_Ba") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- 0
      cross[[i]]$paternal_Ra <- 1/2
      cross[[i]]$paternal_BA <- 0
      cross[[i]]$paternal_Ba <- 1/2
    }
    if (cross[[i]]$father == "BA_BA") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- 0
      cross[[i]]$paternal_Ra <- 0
      cross[[i]]$paternal_BA <- 1
      cross[[i]]$paternal_Ba <- 0
    }
    if (cross[[i]]$father == "BA_Ba") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- 0
      cross[[i]]$paternal_Ra <- 0
      cross[[i]]$paternal_BA <- 1/2
      cross[[i]]$paternal_Ba <- 1/2
    }
    if (cross[[i]]$father == "Ba_Ba") {
      cross[[i]]$paternal_HA <- 0
      cross[[i]]$paternal_Ha <- 0
      cross[[i]]$paternal_hA <- 0
      cross[[i]]$paternal_ha <- 0
      cross[[i]]$paternal_RA <- 0
      cross[[i]]$paternal_Ra <- 0
      cross[[i]]$paternal_BA <- 0
      cross[[i]]$paternal_Ba <- 1
    }
  }

  ## Calculate offspring genotype frequencies:
  for (i in 1:numCrosses) {

    cross[[i]]$HA_HA_f <- 0.5 * cross[[i]]$maternal_HA * cross[[i]]$paternal_HA
    cross[[i]]$HA_HA_m <- 0.5 * cross[[i]]$maternal_HA * cross[[i]]$paternal_HA
    cross[[i]]$HA_Ha_f <- 0.5 * (cross[[i]]$maternal_Ha * cross[[i]]$paternal_HA
                                 + cross[[i]]$maternal_HA * cross[[i]]$paternal_Ha)
    cross[[i]]$HA_Ha_m <- 0.5 * (cross[[i]]$maternal_Ha * cross[[i]]$paternal_HA
                                 + cross[[i]]$maternal_HA * cross[[i]]$paternal_Ha)
    cross[[i]]$HA_hA_f <- 0.5 * (cross[[i]]$maternal_HA * cross[[i]]$paternal_hA
                                 + cross[[i]]$maternal_hA * cross[[i]]$paternal_HA)
    cross[[i]]$HA_hA_m <- 0.5 * (cross[[i]]$maternal_HA * cross[[i]]$paternal_hA
                                 + cross[[i]]$maternal_hA * cross[[i]]$paternal_HA)
    cross[[i]]$HA_ha_f <- 0.5 * (cross[[i]]$maternal_HA * cross[[i]]$paternal_ha
                                 + cross[[i]]$maternal_ha * cross[[i]]$paternal_HA)
    cross[[i]]$HA_ha_m <- 0.5 * (cross[[i]]$maternal_HA * cross[[i]]$paternal_ha
                                 + cross[[i]]$maternal_ha * cross[[i]]$paternal_HA)
    cross[[i]]$HA_RA_f <- 0.5 * (cross[[i]]$maternal_HA * cross[[i]]$paternal_RA
                                 + cross[[i]]$maternal_RA * cross[[i]]$paternal_HA)
    cross[[i]]$HA_RA_m <- 0.5 * (cross[[i]]$maternal_HA * cross[[i]]$paternal_RA
                                 + cross[[i]]$maternal_RA * cross[[i]]$paternal_HA)
    cross[[i]]$HA_Ra_f <- 0.5 * (cross[[i]]$maternal_HA * cross[[i]]$paternal_Ra
                                 + cross[[i]]$maternal_Ra * cross[[i]]$paternal_HA)
    cross[[i]]$HA_Ra_m <- 0.5 * (cross[[i]]$maternal_HA * cross[[i]]$paternal_Ra
                                 + cross[[i]]$maternal_Ra * cross[[i]]$paternal_HA)
    cross[[i]]$HA_BA_f <- 0.5 * (cross[[i]]$maternal_HA * cross[[i]]$paternal_BA
                                 + cross[[i]]$maternal_BA * cross[[i]]$paternal_HA)
    cross[[i]]$HA_BA_m <- 0.5 * (cross[[i]]$maternal_HA * cross[[i]]$paternal_BA
                                 + cross[[i]]$maternal_BA * cross[[i]]$paternal_HA)
    cross[[i]]$HA_Ba_f <- 0.5 * (cross[[i]]$maternal_HA * cross[[i]]$paternal_Ba
                                 + cross[[i]]$maternal_Ba * cross[[i]]$paternal_HA)
    cross[[i]]$HA_Ba_m <- 0.5 * (cross[[i]]$maternal_HA * cross[[i]]$paternal_Ba
                                 + cross[[i]]$maternal_Ba * cross[[i]]$paternal_HA)

    cross[[i]]$Ha_Ha_f <- 0.5 * cross[[i]]$maternal_Ha * cross[[i]]$paternal_Ha
    cross[[i]]$Ha_Ha_m <- 0.5 * cross[[i]]$maternal_Ha * cross[[i]]$paternal_Ha
    cross[[i]]$Ha_hA_f <- 0.5 * (cross[[i]]$maternal_Ha * cross[[i]]$paternal_hA
                                 + cross[[i]]$maternal_hA * cross[[i]]$paternal_Ha)
    cross[[i]]$Ha_hA_m <- 0.5 * (cross[[i]]$maternal_Ha * cross[[i]]$paternal_hA
                                 + cross[[i]]$maternal_hA * cross[[i]]$paternal_Ha)
    cross[[i]]$Ha_ha_f <- 0.5 * (cross[[i]]$maternal_Ha * cross[[i]]$paternal_ha
                               + cross[[i]]$maternal_ha * cross[[i]]$paternal_Ha)
    cross[[i]]$Ha_ha_m <- 0.5 * (cross[[i]]$maternal_Ha * cross[[i]]$paternal_ha
                               + cross[[i]]$maternal_ha * cross[[i]]$paternal_Ha)
    cross[[i]]$Ha_RA_f <- 0.5 * (cross[[i]]$maternal_Ha * cross[[i]]$paternal_RA
                                 + cross[[i]]$maternal_RA * cross[[i]]$paternal_Ha)
    cross[[i]]$Ha_RA_m <- 0.5 * (cross[[i]]$maternal_Ha * cross[[i]]$paternal_RA
                                 + cross[[i]]$maternal_RA * cross[[i]]$paternal_Ha)
    cross[[i]]$Ha_Ra_f <- 0.5 * (cross[[i]]$maternal_Ha * cross[[i]]$paternal_Ra
                                 + cross[[i]]$maternal_Ra * cross[[i]]$paternal_Ha)
    cross[[i]]$Ha_Ra_m <- 0.5 * (cross[[i]]$maternal_Ha * cross[[i]]$paternal_Ra
                                 + cross[[i]]$maternal_Ra * cross[[i]]$paternal_Ha)
    cross[[i]]$Ha_BA_f <- 0.5 * (cross[[i]]$maternal_Ha * cross[[i]]$paternal_BA
                                 + cross[[i]]$maternal_BA * cross[[i]]$paternal_Ha)
    cross[[i]]$Ha_BA_m <- 0.5 * (cross[[i]]$maternal_Ha * cross[[i]]$paternal_BA
                                 + cross[[i]]$maternal_BA * cross[[i]]$paternal_Ha)
    cross[[i]]$Ha_Ba_f <- 0.5 * (cross[[i]]$maternal_Ha * cross[[i]]$paternal_Ba
                                 + cross[[i]]$maternal_Ba * cross[[i]]$paternal_Ha)
    cross[[i]]$Ha_Ba_m <- 0.5 * (cross[[i]]$maternal_Ha * cross[[i]]$paternal_Ba
                                 + cross[[i]]$maternal_Ba * cross[[i]]$paternal_Ha)

    cross[[i]]$hA_hA_f <- 0.5 * cross[[i]]$maternal_hA * cross[[i]]$paternal_hA
    cross[[i]]$hA_hA_m <- 0.5 * cross[[i]]$maternal_hA * cross[[i]]$paternal_hA
    cross[[i]]$hA_ha_f <- 0.5 * (cross[[i]]$maternal_ha * cross[[i]]$paternal_hA
                                 + cross[[i]]$maternal_hA * cross[[i]]$paternal_ha)
    cross[[i]]$hA_ha_m <- 0.5 * (cross[[i]]$maternal_ha * cross[[i]]$paternal_hA
                                 + cross[[i]]$maternal_hA * cross[[i]]$paternal_ha)
    cross[[i]]$hA_RA_f <- 0.5 * (cross[[i]]$maternal_RA * cross[[i]]$paternal_hA
                                 + cross[[i]]$maternal_hA * cross[[i]]$paternal_RA)
    cross[[i]]$hA_RA_m <- 0.5 * (cross[[i]]$maternal_RA * cross[[i]]$paternal_hA
                                 + cross[[i]]$maternal_hA * cross[[i]]$paternal_RA)
    cross[[i]]$hA_Ra_f <- 0.5 * (cross[[i]]$maternal_Ra * cross[[i]]$paternal_hA
                                 + cross[[i]]$maternal_hA * cross[[i]]$paternal_Ra)
    cross[[i]]$hA_Ra_m <- 0.5 * (cross[[i]]$maternal_Ra * cross[[i]]$paternal_hA
                                 + cross[[i]]$maternal_hA * cross[[i]]$paternal_Ra)
    cross[[i]]$hA_BA_f <- 0.5 * (cross[[i]]$maternal_BA * cross[[i]]$paternal_hA
                                 + cross[[i]]$maternal_hA * cross[[i]]$paternal_BA)
    cross[[i]]$hA_BA_m <- 0.5 * (cross[[i]]$maternal_BA * cross[[i]]$paternal_hA
                                 + cross[[i]]$maternal_hA * cross[[i]]$paternal_BA)
    cross[[i]]$hA_Ba_f <- 0.5 * (cross[[i]]$maternal_Ba * cross[[i]]$paternal_hA
                                 + cross[[i]]$maternal_hA * cross[[i]]$paternal_Ba)
    cross[[i]]$hA_Ba_m <- 0.5 * (cross[[i]]$maternal_Ba * cross[[i]]$paternal_hA
                                 + cross[[i]]$maternal_hA * cross[[i]]$paternal_Ba)

    cross[[i]]$ha_ha_f <- 0.5 * cross[[i]]$maternal_ha * cross[[i]]$paternal_ha
    cross[[i]]$ha_ha_m <- 0.5 * cross[[i]]$maternal_ha * cross[[i]]$paternal_ha
    cross[[i]]$ha_RA_f <- 0.5 * (cross[[i]]$maternal_RA * cross[[i]]$paternal_ha
                                 + cross[[i]]$maternal_ha * cross[[i]]$paternal_RA)
    cross[[i]]$ha_RA_m <- 0.5 * (cross[[i]]$maternal_RA * cross[[i]]$paternal_ha
                                 + cross[[i]]$maternal_ha * cross[[i]]$paternal_RA)
    cross[[i]]$ha_Ra_f <- 0.5 * (cross[[i]]$maternal_Ra * cross[[i]]$paternal_ha
                                 + cross[[i]]$maternal_ha * cross[[i]]$paternal_Ra)
    cross[[i]]$ha_Ra_m <- 0.5 * (cross[[i]]$maternal_Ra * cross[[i]]$paternal_ha
                                 + cross[[i]]$maternal_ha * cross[[i]]$paternal_Ra)
    cross[[i]]$ha_BA_f <- 0.5 * (cross[[i]]$maternal_BA * cross[[i]]$paternal_ha
                                 + cross[[i]]$maternal_ha * cross[[i]]$paternal_BA)
    cross[[i]]$ha_BA_m <- 0.5 * (cross[[i]]$maternal_BA * cross[[i]]$paternal_ha
                                 + cross[[i]]$maternal_ha * cross[[i]]$paternal_BA)
    cross[[i]]$ha_Ba_f <- 0.5 * (cross[[i]]$maternal_Ba * cross[[i]]$paternal_ha
                                 + cross[[i]]$maternal_ha * cross[[i]]$paternal_Ba)
    cross[[i]]$ha_Ba_m <- 0.5 * (cross[[i]]$maternal_Ba * cross[[i]]$paternal_ha
                                 + cross[[i]]$maternal_ha * cross[[i]]$paternal_Ba)

    cross[[i]]$RA_RA_f <- 0.5 * cross[[i]]$maternal_RA * cross[[i]]$paternal_RA
    cross[[i]]$RA_RA_m <- 0.5 * cross[[i]]$maternal_RA * cross[[i]]$paternal_RA
    cross[[i]]$RA_Ra_f <- 0.5 * (cross[[i]]$maternal_Ra * cross[[i]]$paternal_RA
                                 + cross[[i]]$maternal_RA * cross[[i]]$paternal_Ra)
    cross[[i]]$RA_Ra_m <- 0.5 * (cross[[i]]$maternal_Ra * cross[[i]]$paternal_RA
                                 + cross[[i]]$maternal_RA * cross[[i]]$paternal_Ra)
    cross[[i]]$RA_BA_f <- 0.5 * (cross[[i]]$maternal_BA * cross[[i]]$paternal_RA
                                 + cross[[i]]$maternal_RA * cross[[i]]$paternal_BA)
    cross[[i]]$RA_BA_m <- 0.5 * (cross[[i]]$maternal_BA * cross[[i]]$paternal_RA
                                 + cross[[i]]$maternal_RA * cross[[i]]$paternal_BA)
    cross[[i]]$RA_Ba_f <- 0.5 * (cross[[i]]$maternal_Ba * cross[[i]]$paternal_RA
                                 + cross[[i]]$maternal_RA * cross[[i]]$paternal_Ba)
    cross[[i]]$RA_Ba_m <- 0.5 * (cross[[i]]$maternal_Ba * cross[[i]]$paternal_RA
                                 + cross[[i]]$maternal_RA * cross[[i]]$paternal_Ba)

    cross[[i]]$Ra_Ra_f <- 0.5 * cross[[i]]$maternal_Ra * cross[[i]]$paternal_Ra
    cross[[i]]$Ra_Ra_m <- 0.5 * cross[[i]]$maternal_Ra * cross[[i]]$paternal_Ra
    cross[[i]]$Ra_BA_f <- 0.5 * (cross[[i]]$maternal_BA * cross[[i]]$paternal_Ra
                                 + cross[[i]]$maternal_Ra * cross[[i]]$paternal_BA)
    cross[[i]]$Ra_BA_m <- 0.5 * (cross[[i]]$maternal_BA * cross[[i]]$paternal_Ra
                                 + cross[[i]]$maternal_Ra * cross[[i]]$paternal_BA)
    cross[[i]]$Ra_Ba_f <- 0.5 * (cross[[i]]$maternal_Ba * cross[[i]]$paternal_Ra
                                 + cross[[i]]$maternal_Ra * cross[[i]]$paternal_Ba)
    cross[[i]]$Ra_Ba_m <- 0.5 * (cross[[i]]$maternal_Ba * cross[[i]]$paternal_Ra
                                 + cross[[i]]$maternal_Ra * cross[[i]]$paternal_Ba)

    cross[[i]]$BA_BA_f <- 0.5 * cross[[i]]$maternal_BA * cross[[i]]$paternal_BA
    cross[[i]]$BA_BA_m <- 0.5 * cross[[i]]$maternal_BA * cross[[i]]$paternal_BA
    cross[[i]]$BA_Ba_f <- 0.5 * (cross[[i]]$maternal_Ba * cross[[i]]$paternal_BA
                                 + cross[[i]]$maternal_BA * cross[[i]]$paternal_Ba)
    cross[[i]]$BA_Ba_m <- 0.5 * (cross[[i]]$maternal_Ba * cross[[i]]$paternal_BA
                                 + cross[[i]]$maternal_BA * cross[[i]]$paternal_Ba)

    cross[[i]]$Ba_Ba_f <- 0.5 * cross[[i]]$maternal_Ba * cross[[i]]$paternal_Ba
    cross[[i]]$Ba_Ba_m <- 0.5 * cross[[i]]$maternal_Ba * cross[[i]]$paternal_Ba

  }


  #############################################################################
  ## Jared Addition
  #############################################################################

  #these are here as notes. First thing was created by hand from the above code.
  #Second thing was written to cut, sort, and paste them so the alleles are in
  #alphabetical order
  # #copied and built by hand from above naming
  # genotypes  <- c("HAHA","HAHW","HAWA","HAWW","HARA","HARW","HABA","HABW","HWHW","HWWA","HWWW",
  #                 "HWRA","HWRW","HWBA","HWBW","WAWA","WAWW","WARA","WARW","WABA","WABW","WWWW",
  #                 "WWRA","WWRW","WWBA","WWBW","RARA","RARW","RABA","RABW","RWRW","RWBA","RWBW",
  #                 "BABA","BABW","BWBW")
  #
  # genotypes <- vapply(X = genotypes, FUN = function(x){
  #   paste0(
  #     sort(
  #       c(substr(x = x, start = 1, stop = 2),
  #         substr(x = x, start = 3, stop = 4))
  #     ),
  #     collapse = "")
  # },
  # FUN.VALUE = character(length = 1L))

  genotypes <- c("HAHA","HAHW","HAWA","HAWW","HARA","HARW","BAHA","BWHA","HWHW",
                 "HWWA","HWWW","HWRA","HWRW","BAHW","BWHW","WAWA","WAWW","RAWA",
                 "RWWA","BAWA","BWWA","WWWW","RAWW","RWWW","BAWW","BWWW","RARA",
                 "RARW","BARA","BWRA","RWRW","BARW","BWRW","BABA","BABW","BWBW")
  numGenotypes <- length(genotypes)
  offset <- 18 # empirically determined from his list



  # offset is 2+num female alleles + num male alleles
  # or you can count it from his cross object
  tMatrix <- array(data = 0,
                   dim = list(numGenotypes,numGenotypes,numGenotypes),
                   dimnames = list(genotypes,genotypes,genotypes))

  #fill tMatrix from list above
  for(column in 1:numGenotypes){
    for(row in 1:numGenotypes){
      for(depth in 1:numGenotypes){

        tMatrix[row, column, depth] <- cross[[row+((column-1)%%numGenotypes)*numGenotypes]][[offset+depth]]+
          cross[[row+((column-1)%%numGenotypes)*numGenotypes]][[offset+numGenotypes+depth]]

      }
    }
  }#end fill loop

  #protection from underflow errors
  tMatrix[tMatrix < .Machine$double.eps/2] <- 0

  ## initialize viability mask. No mother-specific death.
  viabilityMask <- array(data = 1, dim = c(numGenotypes, numGenotypes, numGenotypes),
                         dimnames = list(genotypes, genotypes, genotypes))

  ## genotype-specific modifiers
  if(is.null(eta)){eta = setNames(object = rep.int(x = 1, times = numGenotypes), nm = genotypes)}        # genotype-specific mating fitness
  if(is.null(phi)){phi = setNames(object = rep.int(x = 0.5, times = numGenotypes), nm = genotypes)}      # genotype-specific sex ratio at emergence
  if(is.null(omega)){omega = setNames(object = rep.int(x = 1, times = numGenotypes), nm = genotypes)}    # genotype-specific multiplicative modifier of adult mortality
  if(is.null(xiF)){xiF = setNames(object = rep.int(x = 1, times = numGenotypes), nm = genotypes)}        # genotype-specific female pupatory success
  if(is.null(xiM)){xiM = setNames(object = rep.int(x = 1, times = numGenotypes), nm = genotypes)}        # genotype-specific male pupatory success
  if(is.null(s)){s = setNames(object = rep.int(x = 1, times = numGenotypes), nm = genotypes)}            # genotype-specific fractional reduction(increase) in fertility


  ## put everything into a labeled list to return
  return(list(
    ih = tMatrix,
    tau = viabilityMask,
    genotypesID = genotypes,
    genotypesN = numGenotypes,
    wildType = "WWWW",
    eta = eta,
    phi = phi,
    omega = omega,
    xiF = xiF,
    xiM = xiM,
    s = s,
    releaseType = "HAHA"
  ))

}
