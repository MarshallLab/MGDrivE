#' Kernels Parameters
#'
#' A named list containing maximum likelihood fitted parameter values from
#' mosquito dispersal estimates.
#'
#' @name kernels
#' @docType data
#'
#' @usage data(kernels)
#'
#' @format named list with 5 elements:
#' \describe{
#'   \item{lnorm_mean}{log mean of log-normal density}
#'   \item{lnorm_sd}{log standard deviation of log-normal density}
#'   \item{gamma_shape}{shape parameter of gamma density}
#'   \item{gamma_sd}{rate parameter of gamma density}
#'   \item{exp_rate}{rate parameter of exponential density}
#' }
"kernels"

#' Movement Matrix: All 2
#'
#' A movement matrix for simulation with 3 patches.
#'
#' @name moveMatAll2
#' @docType data
#'
#' @usage data(moveMatAll2)
#'
#'
#' @format A matrix with 3 rows and 3 columns:
#' \describe{
#'   Patches 1 and 3 are sources for patch 2 which is a sink.
#' }
"moveMatAll2"

#' Movement Matrix: Diagonal
#'
#' A movement matrix for simulation with 3 patches.
#'
#' @name moveMatDiag
#' @docType data
#'
#' @usage data(moveMatDiag)
#'
#'
#' @format A matrix with 3 rows and 3 columns:
#' \describe{
#'   3 independent patches.
#' }
"moveMatDiag"

#' Movement Matrix: Cascade 3
#'
#' A movement matrix for simulation with 3 patches.
#'
#'
#' @name moveMatCascade3
#' @docType data
#'
#' @usage data(moveMatCascade3)
#'
#' @format A matrix with 3 rows and 3 columns:
#' \describe{
#'   Mosquitoes in patch 1 have equal probability to stay or move to 2; mosquitoes
#'   in patch 2 have equal probability to stay or move to 3; mosquitoes in patch 3 stay there.
#' }
"moveMatCascade3"

#' Movement Matrix: Independent 3
#'
#' A movement matrix for simulation with 3 patches.
#'
#' @name moveMatIndependent3
#' @docType data
#'
#' @usage data(moveMatIndependent3)
#'
#'
#' @format A matrix with 3 rows and 3 columns:
#' \describe{
#'   Mosquitoes in patch 1 stay with probability 0.975, move to patch 2 with probability 0.025,
#'   mosquitoes in patch 2 and 3 stay in their patches.
#' }
"moveMatIndependent3"

#' Movement Matrix: Die
#'
#' A movement matrix for simulation with 3 patches.
#'
#' @name moveMatDie
#' @docType data
#'
#' @usage data(moveMatDie)
#'
#' @format A matrix with 3 rows and 3 columns:
#' \describe{
#'   All entries of matrix are 0 for testing that all mosquitoes will be killed.
#' }
"moveMatDie"

#' Movement Matrix: Triple
#'
#' A movement matrix for simulation with 3 patches.
#'
#' @name moveMatTriple
#' @docType data
#'
#' @usage data(moveMatTriple)
#'
#' @format A matrix with 3 rows and 3 columns:
#' \describe{
#'   All entries of matrix are 1 for testing that mosquitoes will be produced.
#' }
"moveMatTriple"

#' Movement Matrix: Mixed Spill
#'
#' A movement matrix for simulation with 3 patches.
#'
#' @name moveMatMixedSpil
#' @docType data
#'
#' @usage data(moveMatMixedSpil)
#'
#' @format A matrix with 3 rows and 3 columns:
#' \describe{
#'   Mosquitoes in patch 1 stay with probability 0.999, move to patch 2 with probability 0.001,
#'   mosquitoes in patch 2 and 3 stay in their patches.

#' }
"moveMatMixedSpil"

#' Movement Matrix: Diagonal One City
#'
#' A movement matrix for simulation with 1 patch.
#'
#' @name moveMatDiagOneCity
#' @docType data
#'
#' @usage data(moveMatDiagOneCity)
#'
#' @format A matrix with 1 rows and 1 columns:
#' \describe{
#'   A 1 by 1 matrix with entry 1.
#' }
"moveMatDiagOneCity"

#' Movement Matrix: Tri-diagonal
#'
#' A movement matrix for simulation with 12 patches.
#'
#' @name moveMatTriDiagonal
#' @docType data
#'
#' @usage data(moveMatTriDiagonal)
#'
#' @format A matrix with 12 rows and 12 columns:
#' \describe{
#'   Tri-diagonal matrix with approximately 0.985 probability on diagonal and rest of probability mass on k-1 and k+1 off-diagonal elements.
#' }
"moveMatTriDiagonal"

#' Movement Matrix: Tale of Two Cities
#'
#' A movement matrix for simulation with 2 patches.
#'
#' @name moveMatTaleOfTwoCities
#' @docType data
#'
#' @usage data(moveMatTaleOfTwoCities)
#'
#' @format A matrix with 2 rows and 2 columns:
#' \describe{
#'   Mosquitoes do not move between the two patches.
#' }
"moveMatTaleOfTwoCities"
