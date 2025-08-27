###############################################################################
#        ____        __       __
#       / __ \____ _/ /______/ /_
#      / /_/ / __ `/ __/ ___/ __ \
#     / ____/ /_/ / /_/ /__/ / / /
#    /_/    \__,_/\__/\___/_/ /_/
#
#   Patch Migration
#   Marshall Lab
#   jared_bennett@berkeley.edu
#   December 2019
#
###############################################################################
###############################################################################
# Outbound Migration
###############################################################################

# This has been moved to the network level.
# see Network-Migration.R


###############################################################################
# Inbound Migration
###############################################################################

#' Inbound Migration
#'
#' Accumulate all inbound migration to this patch.
#'
#' @param maleIn Vector of inbound migration
#' @param femaleIn Matrix of inbound migration
#'
#' @keywords internal
#'
oneDay_migrationIn_Patch <- function(maleIn, femaleIn){
  private$popMale[] = maleIn
  private$popFemale[] = femaleIn
}
