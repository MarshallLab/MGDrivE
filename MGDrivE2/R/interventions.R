################################################################################
#
#   MGDrivE2: interventions
#   Marshall Lab
#   Agastya Mondal (agastya_mondal@berkeley.edu)
#   March 2022
#
################################################################################


#' This set of functions modifies mosquito life history parameters in the presence of
#' adult interventions - indoor residual spraying (IRS) and insecticide treated nets (ITN)
#' This is based on the work of Le Menach et al (2007) and Griffin et al (2010).
#' We vary three parameters in the presence of interventions:
#'  Egg laying rate (beta)
#'  Adult mortality (muF)
#'  Mosquito biting rate (av0)
#'
#' @name add_interventions
#'
#' @param params a named list of parameters
#' @param IRS_cov proportion of humans in the node receiving IRS
#' @param LLIN_cov proportion of humans in the node receiving LLIN
#'
#' @return a vector of the equilibrium number of females in each SEI stage
#'
#' @export
add_interventions <- function(params, IRS_cov, LLIN_cov) {
  if((LLIN_cov == 0) && (IRS_cov == 0)){
    return(params)
  }
  # TODO check that LLIN_cov and ITS_cov are proportions
  params_vc <-
    params # create new vector to hold parameters under vector control
  coverage <- generate_coverage_parameters(IRS_cov, LLIN_cov)

  # feeding rate modifications
  feeding_params <-
    modify_feeding_rate(IRS_cov, LLIN_cov, params_vc, coverage)
  params_vc$av0 <- feeding_params$av0
  params_vc$fv0 <- feeding_params$fv0
  params_vc$Q0 <- feeding_params$Q0
  params_vc$w <- feeding_params$w
  params_vc$z <- feeding_params$z

  # adult mortality modifications
  mortality_params <-
    modify_adult_mortality(IRS_cov, LLIN_cov, params_vc, coverage)
  params_vc$muF <- mortality_params$muF
  params_vc$muM <- mortality_params$muM

  # egg laying rate modifications
  egg_params <- modify_egg_laying(params_vc)
  params_vc$beta <- egg_params$beta

  # modify survival probability
  params_vc$Surv0 <- (params$nEIP * params$qEIP/ (params$nEIP * params$qEIP + mortality_params$muF))^params$nEIP

  return(params_vc)
}

# function which generates the correct coverage values for IRS and ITN alone and together
# using conditional probabilities
generate_coverage_parameters <- function(IRS_cov, LLIN_cov) {
  coverage <- list()
  coverage$c_llin <- LLIN_cov - LLIN_cov * IRS_cov
  coverage$c_irs <- IRS_cov - LLIN_cov * IRS_cov
  coverage$c_irs_llin <- LLIN_cov * IRS_cov
  coverage$c0 <- 1 - LLIN_cov - IRS_cov + LLIN_cov * IRS_cov
  return(coverage)
}

# function that modifies the mosquito feeding rate in the presence of interventions
modify_feeding_rate <-
  function(IRS_cov, LLIN_cov, params, coverage) {
    # first modify the gonotrophic cycle length
    z <- with(c(params, coverage), {
      return(
        Q0 * c_llin * thetaB * r_llin + Q0 * c_irs * thetaI * r_irs + Q0 * c_irs_llin *
          (thetaI - thetaB) * r_irs
        + Q0 * c_irs_llin * thetaB * (r_irs + (1 - r_irs) * r_llin)
      )
    })
    # modified cycle length
    fv0 <- 1 / ((params$tau1 / (1 - z)) + params$tau2)

    # modified biting probability
    w <- with(c(params, coverage), {
      return(
        1 - Q0 + Q0 * c0 + Q0 * c_llin * (1 - thetaB + thetaB * s_llin)
        + Q0 * c_irs * (1 - thetaI + thetaI * s_irs)
        + Q0 * c_irs_llin*((thetaI - thetaB) * s_irs + 1 - thetaI + thetaB *
                            ((1 - r_irs) * s_llin * s_irs))
      )
    })
    Q0 <- 1 - ((1 - params$Q0) / w)

    # finally, return modified biting rate
    # and other parameters for further calculations
    return(list(
      fv0 = fv0,
      Q0 = Q0,
      av0 = fv0 * Q0,
      w = w,
      z = z
    ))
  }

# function that modifies the adult mortality rate in the presence of interventions
modify_adult_mortality <-
  function(IRS_cov, LLIN_cov, params, coverage) {
    p1 <- exp(-params$muF * params$tau1)
    p2 <- exp(-params$muF * params$tau2)
    p1_vc <- (p1 * params$w) / (1 - params$z * p1)
    p_vc <- (p1_vc * p2) ^ params$fv0
    mu <- -log(p_vc)
    return(list(muF = mu,
                muM = mu))
  }

# function to modify egg-laying rate of female mosquitoes
# in presence of interventions
modify_egg_laying <- function(params) {
  beta <- with(params, {
    a <- eps * muF
    b <- exp(muF / fv0) - 1
    return(a / b)
  })
  return(list(beta = beta))
}
