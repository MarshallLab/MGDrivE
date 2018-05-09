########################################################################
#       __  _____________       _       ______
#      /  |/  / ____/ __ \_____(_)   __/ ____/
#     / /|_/ / / __/ / / / ___/ / | / / __/
#    / /  / / /_/ / /_/ / /  / /| |/ / /___
#   /_/  /_/\____/_____/_/  /_/ |___/_____/
#
#   MGDrivE Graphics
#   Marshall Lab
#   November 2017
#
########################################################################

#' Utility to Imitate ggplot2 Colors
#'
#' Sample at equally spaced intervals along the color wheel
#'
#' @param n number of colors
#' @param alpha transparency
#'
#' @export
ggCol_utility <- function(n, alpha = 1) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n]
}
